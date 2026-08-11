/* app.js -- UI layer: canvases, sliders, and the live readout.
 * All mathematics lives in model.js (kept DOM-free); this file only
 * renders state and routes input.
 *
 * Controls:
 *   - each camera has an explicit 6-DOF pose: position x/y/z sliders and
 *     pan/tilt/roll sliders; dragging a camera view turns that camera
 *     first-person style; "aim at molecule" re-points it.
 *   - the molecule has yaw/pitch/roll + x/y/z sliders, and can be dragged
 *     directly in the world view (drag = translate in the view plane,
 *     shift-drag = rotate, alt-drag = translate in depth y).
 *
 * Panels: world overview; camera 1 display; camera 2 display; and the
 * transfer panel -- camera 1's display warped into camera 2's frame by
 * the plane-induced homography (doc section 5.1), compared against where
 * camera 2 actually sees the atoms (yellow rings). The mismatch is
 * parallax: what no single 2D matrix can transfer, and what depth is
 * made of.
 */
'use strict';

(function () {
  var DISP_W = 64, DISP_H = 48, ZOOM = 5;   /* coarse displays, magnified */
  var FOCAL = 52;                            /* pixels, on the 64x48 sensor */
  var BASE_SIGMA_D = Math.sqrt(2) * (1 / Math.sqrt(12)); /* quantization */
  var DEG = Math.PI / 180;
  var OSCALE = 36;                           /* world-view px per metre */

  var mol = MV.makeMolecule();
  var state = {
    cam1: { x: 3.72, y: -1.80, z: 0.75, pan: 154 * DEG, tilt: -10 * DEG, roll: 0 },
    cam2: { x: 3.76, y: 1.82, z: 0.42, pan: -154 * DEG, tilt: -6 * DEG, roll: 0 },
    pose: { yaw: 0.5, pitch: 0.3, roll: 0, x: 0, y: 0, z: 0 },
    autoAim: { 1: true, 2: true },            /* cameras track the molecule */
    selected: 6                               /* the O atom */
  };
  var bufs = {};                              /* low-res buffers per view */

  function camera(c) {
    return MV.makeCameraPose({ pos: [c.x, c.y, c.z],
                               pan: c.pan, tilt: c.tilt, roll: c.roll,
                               f: FOCAL, W: DISP_W, H: DISP_H });
  }
  function molCenter() {
    return [state.pose.x, state.pose.y, state.pose.z];
  }

  /* ---------------- coarse pixel displays ------------------------------- */

  function pixelGrid(cx, cv) {
    cx.strokeStyle = 'rgba(255,255,255,0.07)';
    cx.lineWidth = 1;
    var i;
    for (i = 0; i <= DISP_W; i++) {
      cx.beginPath(); cx.moveTo(i * ZOOM, 0); cx.lineTo(i * ZOOM, cv.height); cx.stroke();
    }
    for (i = 0; i <= DISP_H; i++) {
      cx.beginPath(); cx.moveTo(0, i * ZOOM); cx.lineTo(cv.width, i * ZOOM); cx.stroke();
    }
  }

  function drawDisplay(canvasId, cam, atoms, F, otherCam, otherAtoms) {
    var cv = document.getElementById(canvasId), cx = cv.getContext('2d');
    var buf = document.createElement('canvas');
    buf.width = DISP_W; buf.height = DISP_H;
    var bx = buf.getContext('2d');
    bx.fillStyle = '#10141a'; bx.fillRect(0, 0, DISP_W, DISP_H);

    /* bonds first, then atoms, at sensor resolution: what the pixels see */
    bx.strokeStyle = '#9aa3ad'; bx.lineWidth = 1;
    mol.bonds.forEach(function (b) {
      var p = MV.project(cam, atoms[b[0]]), q = MV.project(cam, atoms[b[1]]);
      if (!p || !q) return;
      bx.beginPath(); bx.moveTo(p.u, p.v); bx.lineTo(q.u, q.v); bx.stroke();
    });
    atoms.forEach(function (X, i) {
      var p = MV.project(cam, X);
      if (!p) return;
      var r = Math.max(1.2, 5.5 / p.z);
      bx.fillStyle = mol.atoms[i].color;
      bx.beginPath(); bx.arc(p.u, p.v, r, 0, 2 * Math.PI); bx.fill();
    });
    bufs[canvasId] = buf;

    /* magnify with hard pixel edges: the display, at the pixel level */
    cx.imageSmoothingEnabled = false;
    cx.clearRect(0, 0, cv.width, cv.height);
    cx.drawImage(buf, 0, 0, cv.width, cv.height);
    pixelGrid(cx, cv);

    /* overlay: selected atom's own pixel + epipolar line from other view */
    var sel = state.selected;
    var p = MV.project(cam, atoms[sel]);
    if (p) {
      var q = MV.quantize(p);
      cx.strokeStyle = '#ffd75e'; cx.lineWidth = 2;
      cx.strokeRect((q.u - 0.5) * ZOOM, (q.v - 0.5) * ZOOM, ZOOM, ZOOM);
    }
    var po = MV.project(otherCam, otherAtoms[sel]);
    if (po) {
      var l = MV.epipolarLine(F, MV.quantize(po)); /* line here from other pixel */
      drawEpiLine(cx, l);
    }
  }

  /* draw a u + b v + c = 0 clipped to the display, in display pixels */
  function drawEpiLine(cx, l) {
    var pts = [], a = l[0], b = l[1], c = l[2];
    function tryPt(u, v) {
      if (u >= -1e-6 && u <= DISP_W + 1e-6 && v >= -1e-6 && v <= DISP_H + 1e-6)
        pts.push([u, v]);
    }
    if (Math.abs(b) > 1e-12) { tryPt(0, -c / b); tryPt(DISP_W, -(c + a * DISP_W) / b); }
    if (Math.abs(a) > 1e-12) { tryPt(-c / a, 0); tryPt(-(c + b * DISP_H) / a, DISP_H); }
    if (pts.length < 2) return;
    cx.strokeStyle = '#5ec8ff'; cx.lineWidth = 1.5;
    cx.setLineDash([6, 4]);
    cx.beginPath();
    cx.moveTo(pts[0][0] * ZOOM, pts[0][1] * ZOOM);
    cx.lineTo(pts[1][0] * ZOOM, pts[1][1] * ZOOM);
    cx.stroke();
    cx.setLineDash([]);
  }

  /* ---------------- transfer panel: view 1 warped into camera 2 --------- */

  function drawTransfer(cam1, cam2, atoms) {
    var cv = document.getElementById('xfer'), cx = cv.getContext('2d');
    var pc = MV.project(cam1, molCenter());
    var stats = { d: NaN, planeMean: 0, planeN: 0, depthSel: NaN };
    cx.clearRect(0, 0, cv.width, cv.height);
    if (!pc || !bufs.view1) {
      cx.fillStyle = '#10141a'; cx.fillRect(0, 0, cv.width, cv.height);
      cx.fillStyle = '#9aa3ad'; cx.font = '13px system-ui, sans-serif';
      cx.fillText('molecule center is behind camera 1', 14, 20);
      return stats;
    }
    var d = pc.z;                             /* assumed scene-plane depth */
    stats.d = d;
    var H = MV.planeHomography(cam1, cam2, d);
    var Hinv = MV.inv3(H);

    /* inverse-warp view 1's low-res buffer into camera 2's frame */
    var src = bufs.view1.getContext('2d').getImageData(0, 0, DISP_W, DISP_H);
    var buf = document.createElement('canvas');
    buf.width = DISP_W; buf.height = DISP_H;
    var bx = buf.getContext('2d');
    var dst = bx.createImageData(DISP_W, DISP_H);
    var x, y, p1, sx, sy, si, di;
    for (y = 0; y < DISP_H; y++) {
      for (x = 0; x < DISP_W; x++) {
        di = 4 * (y * DISP_W + x);
        p1 = MV.applyH(Hinv, { u: x + 0.5, v: y + 0.5 });
        if (p1 && p1.u >= 0 && p1.u < DISP_W && p1.v >= 0 && p1.v < DISP_H) {
          sx = Math.floor(p1.u); sy = Math.floor(p1.v);
          si = 4 * (sy * DISP_W + sx);
          dst.data[di] = src.data[si];
          dst.data[di + 1] = src.data[si + 1];
          dst.data[di + 2] = src.data[si + 2];
        } else {                              /* outside view 1: darker */
          dst.data[di] = 8; dst.data[di + 1] = 9; dst.data[di + 2] = 11;
        }
        dst.data[di + 3] = 255;
      }
    }
    bx.putImageData(dst, 0, 0);
    bufs.xfer = buf;                          /* raw warp, for the diff panel */
    cx.imageSmoothingEnabled = false;
    cx.drawImage(buf, 0, 0, cv.width, cv.height);
    pixelGrid(cx, cv);

    /* where camera 2 ACTUALLY sees each atom: yellow rings to compare */
    atoms.forEach(function (X, i) {
      var p2 = MV.project(cam2, X);
      if (!p2) return;
      var r = Math.max(1.2, 5.5 / p2.z) * ZOOM;
      cx.strokeStyle = i === state.selected ? '#ffd75e' : 'rgba(255,215,94,0.45)';
      cx.lineWidth = i === state.selected ? 2 : 1;
      cx.beginPath(); cx.arc(p2.u * ZOOM, p2.v * ZOOM, r, 0, 2 * Math.PI); cx.stroke();
    });

    /* depth-true transfer of the selected atom's quantized pixel:
     * back-project at its camera-1 depth, reproject into camera 2 */
    var sel = state.selected;
    var p1sel = MV.project(cam1, atoms[sel]);
    var p2sel = MV.project(cam2, atoms[sel]);
    if (p1sel && p2sel) {
      var Xd = MV.backproject(cam1, MV.quantize(p1sel), p1sel.z);
      var pd = MV.project(cam2, Xd);
      if (pd) {
        cx.strokeStyle = '#5ec8ff'; cx.lineWidth = 2;
        cx.beginPath();
        cx.moveTo((pd.u - 1.2) * ZOOM, pd.v * ZOOM);
        cx.lineTo((pd.u + 1.2) * ZOOM, pd.v * ZOOM);
        cx.moveTo(pd.u * ZOOM, (pd.v - 1.2) * ZOOM);
        cx.lineTo(pd.u * ZOOM, (pd.v + 1.2) * ZOOM);
        cx.stroke();
        stats.depthSel = Math.hypot(pd.u - p2sel.u, pd.v - p2sel.v);
      }
    }

    /* plane-transfer error over the atoms seen by both cameras */
    atoms.forEach(function (X) {
      var p1 = MV.project(cam1, X), p2 = MV.project(cam2, X);
      if (!p1 || !p2) return;
      var ph = MV.applyH(H, p1);
      if (!ph) return;
      stats.planeMean += Math.hypot(ph.u - p2.u, ph.v - p2.v);
      stats.planeN++;
    });
    if (stats.planeN) stats.planeMean /= stats.planeN;
    return stats;
  }

  /* ---------------- diff panel: |actual - predicted|, amplified ---------- */

  function drawDiff() {
    var cv = document.getElementById('diff'), cx = cv.getContext('2d');
    var stats = { mean: NaN, max: 0 };
    cx.clearRect(0, 0, cv.width, cv.height);
    if (!bufs.view2 || !bufs.xfer) {
      cx.fillStyle = '#10141a'; cx.fillRect(0, 0, cv.width, cv.height);
      return stats;
    }
    var a = bufs.view2.getContext('2d').getImageData(0, 0, DISP_W, DISP_H);
    var b = bufs.xfer.getContext('2d').getImageData(0, 0, DISP_W, DISP_H);
    var buf = document.createElement('canvas');
    buf.width = DISP_W; buf.height = DISP_H;
    var bx = buf.getContext('2d');
    var dst = bx.createImageData(DISP_W, DISP_H);
    var i, d, amp, sum = 0, n = DISP_W * DISP_H;
    for (i = 0; i < n; i++) {
      d = (Math.abs(a.data[4 * i] - b.data[4 * i]) +
           Math.abs(a.data[4 * i + 1] - b.data[4 * i + 1]) +
           Math.abs(a.data[4 * i + 2] - b.data[4 * i + 2])) / 3;
      sum += d;
      if (d > stats.max) stats.max = d;
      amp = Math.min(255, d * 2);             /* x2 gain for visibility */
      dst.data[4 * i] = amp;
      dst.data[4 * i + 1] = amp * 0.55;
      dst.data[4 * i + 2] = amp * 0.12;
      dst.data[4 * i + 3] = 255;
    }
    stats.mean = sum / n;
    bx.putImageData(dst, 0, 0);
    cx.imageSmoothingEnabled = false;
    cx.drawImage(buf, 0, 0, cv.width, cv.height);
    pixelGrid(cx, cv);
    return stats;
  }

  /* ---------------- error graph: mean diff over interaction time -------- */

  var diffHist = [];                          /* rolling mean-|diff| samples */
  var HIST_MAX = 240;

  function pushDiff(v) {
    if (isNaN(v)) return;
    diffHist.push(v);
    if (diffHist.length > HIST_MAX) diffHist.shift();
  }

  function drawErrGraph(hoverIdx) {
    var cv = document.getElementById('errgraph'), cx = cv.getContext('2d');
    var W = cv.width, H = cv.height;
    var padL = 6, padR = 46, padT = 12, padB = 6;
    cx.clearRect(0, 0, W, H);
    if (diffHist.length < 2) {
      cx.fillStyle = '#4a5563'; cx.font = '10px system-ui, sans-serif';
      cx.fillText('interact to accumulate history', padL + 2, H / 2);
      return;
    }
    var top = Math.max(4, Math.max.apply(null, diffHist) * 1.1);
    var n = diffHist.length;
    function px(i) { return padL + (W - padL - padR) * i / (HIST_MAX - 1); }
    function py(v) { return padT + (H - padT - padB) * (1 - v / top); }
    /* recessive grid: top-of-scale and midpoint */
    cx.strokeStyle = '#d7dce2'; cx.lineWidth = 1;
    cx.fillStyle = '#4a5563'; cx.font = '9px system-ui, sans-serif';
    [top, top / 2].forEach(function (v) {
      cx.beginPath(); cx.moveTo(padL, py(v)); cx.lineTo(W - padR, py(v)); cx.stroke();
      cx.fillText(v.toFixed(v < 10 ? 1 : 0), padL + 1, py(v) - 2);
    });
    /* the series: one thin accent line */
    cx.strokeStyle = '#2d6cdf'; cx.lineWidth = 2;
    cx.lineJoin = 'round';
    cx.beginPath();
    var i;
    for (i = 0; i < n; i++) {
      if (i === 0) cx.moveTo(px(i), py(diffHist[i]));
      else cx.lineTo(px(i), py(diffHist[i]));
    }
    cx.stroke();
    /* current value: accent dot + ink label at the line's end */
    var lastX = px(n - 1), lastY = py(diffHist[n - 1]);
    cx.fillStyle = '#2d6cdf';
    cx.beginPath(); cx.arc(lastX, lastY, 3, 0, 2 * Math.PI); cx.fill();
    cx.fillStyle = '#1c2733'; cx.font = '11px system-ui, sans-serif';
    cx.fillText('mean |Δ| ' + diffHist[n - 1].toFixed(1),
                Math.min(lastX + 6, W - padR + 2),
                Math.max(padT + 8, Math.min(H - 4, lastY + 4)));
    /* hover crosshair + value */
    if (hoverIdx != null && hoverIdx >= 0 && hoverIdx < n) {
      var hx = px(hoverIdx), hy = py(diffHist[hoverIdx]);
      cx.strokeStyle = '#9aa3ad'; cx.lineWidth = 1;
      cx.setLineDash([3, 3]);
      cx.beginPath(); cx.moveTo(hx, padT); cx.lineTo(hx, H - padB); cx.stroke();
      cx.setLineDash([]);
      cx.fillStyle = '#2d6cdf';
      cx.beginPath(); cx.arc(hx, hy, 3, 0, 2 * Math.PI); cx.fill();
      cx.fillStyle = '#1c2733'; cx.font = '10px system-ui, sans-serif';
      cx.fillText(diffHist[hoverIdx].toFixed(1),
                  Math.min(hx + 5, W - 30), Math.max(padT + 8, hy - 5));
    }
  }

  (function () {
    var cv = document.getElementById('errgraph');
    cv.addEventListener('mousemove', function (e) {
      var rect = cv.getBoundingClientRect();
      var frac = (e.clientX - rect.left - 6) / (cv.width - 6 - 46);
      var idx = Math.round(frac * (HIST_MAX - 1));
      drawErrGraph(Math.max(0, Math.min(diffHist.length - 1, idx)));
    });
    cv.addEventListener('mouseleave', function () { drawErrGraph(null); });
  })();

  /* ---------------- world overview -------------------------------------- */

  function oproj(X) {
    var cv = document.getElementById('overview');
    return [cv.width / 2 + OSCALE * (X[0] - 0.45 * X[1]),
            cv.height / 2 - OSCALE * (X[2] * 0.9 + 0.28 * X[1])];
  }

  function drawOverview(atoms, cam1, cam2) {
    var cv = document.getElementById('overview'), cx = cv.getContext('2d');
    cx.clearRect(0, 0, cv.width, cv.height);
    /* working volume: 2 m cube wireframe */
    cx.strokeStyle = '#c3cad2'; cx.lineWidth = 1;
    var s = 1.1, k, corners = [], edges = [
      [0, 1], [1, 3], [3, 2], [2, 0], [4, 5], [5, 7], [7, 6], [6, 4],
      [0, 4], [1, 5], [2, 6], [3, 7]];
    for (k = 0; k < 8; k++)
      corners.push([(k & 1 ? s : -s), (k & 2 ? s : -s), (k & 4 ? s : -s)]);
    edges.forEach(function (e) {
      var p = oproj(corners[e[0]]), q = oproj(corners[e[1]]);
      cx.beginPath(); cx.moveTo(p[0], p[1]); cx.lineTo(q[0], q[1]); cx.stroke();
    });
    /* molecule */
    cx.strokeStyle = '#6b7280'; cx.lineWidth = 2;
    mol.bonds.forEach(function (b) {
      var p = oproj(atoms[b[0]]), q = oproj(atoms[b[1]]);
      cx.beginPath(); cx.moveTo(p[0], p[1]); cx.lineTo(q[0], q[1]); cx.stroke();
    });
    atoms.forEach(function (X, i) {
      var p = oproj(X);
      cx.fillStyle = mol.atoms[i].color;
      cx.beginPath();
      cx.arc(p[0], p[1], i === state.selected ? 6 : 4, 0, 2 * Math.PI);
      cx.fill();
      if (i === state.selected) {
        cx.strokeStyle = '#ffd75e'; cx.lineWidth = 2;
        cx.stroke();
      }
    });
    /* cameras: center dot, viewing frustum out to 1.3 m, label */
    [[cam1, 'cam 1'], [cam2, 'cam 2']].forEach(function (cc) {
      var cam = cc[0], p = oproj(cam.C);
      /* frustum: rays through the four sensor corners (K^-1 corners),
       * so the funnel is the camera's true field of view */
      var Ki = MV.inv3(cam.K), depth = 1.3;
      var far = [[0, 0], [DISP_W, 0], [DISP_W, DISP_H], [0, DISP_H]]
        .map(function (uv) {
          var dc = MV.matVec(Ki, [uv[0], uv[1], 1]);
          return oproj(MV.add(cam.C,
            MV.matVec(MV.transpose(cam.R), MV.scale(dc, depth))));
        });
      var i, j;
      for (i = 0; i < 4; i++) {              /* side faces, faintly filled */
        j = (i + 1) % 4;
        cx.fillStyle = 'rgba(45,108,223,0.07)';
        cx.beginPath();
        cx.moveTo(p[0], p[1]);
        cx.lineTo(far[i][0], far[i][1]);
        cx.lineTo(far[j][0], far[j][1]);
        cx.closePath(); cx.fill();
      }
      cx.strokeStyle = 'rgba(45,108,223,0.55)'; cx.lineWidth = 1;
      for (i = 0; i < 4; i++) {              /* funnel edges + far rim */
        j = (i + 1) % 4;
        cx.beginPath(); cx.moveTo(p[0], p[1]);
        cx.lineTo(far[i][0], far[i][1]); cx.stroke();
        cx.beginPath(); cx.moveTo(far[i][0], far[i][1]);
        cx.lineTo(far[j][0], far[j][1]); cx.stroke();
      }
      cx.fillStyle = '#1c2733';
      cx.beginPath(); cx.arc(p[0], p[1], 5, 0, 2 * Math.PI); cx.fill();
      cx.font = '12px system-ui, sans-serif';
      cx.fillText(cc[1], p[0] + 8, p[1] + 14);
    });
    drawWorldAxes(cx, cv);
  }

  /* labeled world-frame axis triad, anchored in the lower-left corner */
  function drawWorldAxes(cx, cv) {
    var ax = 34, ay = cv.height - 30, len = 0.55;
    /* screen direction of a world unit vector, from the oproj mapping */
    function dir(X) {
      return [OSCALE * (X[0] - 0.45 * X[1]) * len,
              -OSCALE * (X[2] * 0.9 + 0.28 * X[1]) * len];
    }
    [[[1, 0, 0], 'x'], [[0, 1, 0], 'y'], [[0, 0, 1], 'z']].forEach(function (a) {
      var d = dir(a[0]);
      cx.strokeStyle = '#4a5563'; cx.lineWidth = 1.5;
      cx.beginPath(); cx.moveTo(ax, ay); cx.lineTo(ax + d[0], ay + d[1]); cx.stroke();
      cx.fillStyle = '#4a5563'; cx.font = 'italic 12px Georgia, serif';
      cx.fillText(a[1], ax + d[0] * 1.25 - 3, ay + d[1] * 1.25 + 4);
    });
  }

  /* ---------------- readout ---------------------------------------------- */

  function fmt(x, d) { return (x >= 0 ? ' ' : '') + x.toFixed(d); }
  function fmtMat(M, d) {
    return M.map(function (r) {
      return r.map(function (x) { return fmt(x, d); }).join('  ');
    }).join('\n');
  }

  function updateReadout(cam1, cam2, atoms, F, xstats) {
    var sel = state.selected;
    var name = mol.atoms[sel].el + '(' + sel + ')';
    var p1 = MV.project(cam1, atoms[sel]);
    var p2 = MV.project(cam2, atoms[sel]);
    var out = document.getElementById('readout');
    if (!p1 || !p2 ||
        p1.u < 0 || p1.u >= DISP_W || p1.v < 0 || p1.v >= DISP_H ||
        p2.u < 0 || p2.u >= DISP_W || p2.v < 0 || p2.v >= DISP_H) {
      out.textContent = 'selected atom ' + name +
        ' is not visible in both cameras — move or re-aim a camera\n' +
        '(the "aim at molecule" buttons re-point a camera at the molecule).';
      return;
    }
    var q1 = MV.quantize(p1), q2 = MV.quantize(p2);
    var res = MV.epipolarResidual(F, q1, q2);
    var Xhat = MV.triangulate(cam1, q1, cam2, q2);
    var Xtrue = atoms[sel];
    var err = Xhat ? MV.norm(MV.sub(Xhat, Xtrue)) : NaN;
    var rms = Xhat ? MV.reprojError(
      [{ cam: cam1, p: q1 }, { cam: cam2, p: q2 }], Xhat) : NaN;
    var B = MV.norm(MV.sub(cam2.C, cam1.C));
    var Z = (p1.z + p2.z) / 2;
    var sZ = MV.depthSigma(Z, FOCAL, B, BASE_SIGMA_D);

    out.textContent =
      'selected atom ' + name + '   true X = (' +
        Xtrue.map(function (x) { return x.toFixed(3); }).join(', ') + ') m\n' +
      '\n' +
      'projection (doc §4, eq. pinhole):\n' +
      '  camera 1 pixel: continuous (' + p1.u.toFixed(2) + ', ' + p1.v.toFixed(2) +
        ')  → display pixel center (' + q1.u + ', ' + q1.v + ')\n' +
      '  camera 2 pixel: continuous (' + p2.u.toFixed(2) + ', ' + p2.v.toFixed(2) +
        ')  → display pixel center (' + q2.u + ', ' + q2.v + ')\n' +
      '\n' +
      'fundamental matrix F (doc §6, from E = [t]×R):\n' +
      fmtMat(F.map(function (r) {
        return r.map(function (x) { return x * 100; });
      }), 4) + '   (×10⁻²)\n' +
      '  epipolar residual x₂ᵀF x₁ = ' + res.toExponential(2) +
        '   (0 exactly for the continuous pixels; nonzero = quantization)\n' +
      '\n' +
      'view transfer 1 → 2 (predicted panel; doc §5.1, eq. planehom):\n' +
      '  plane homography at camera-1 depth d = ' + xstats.d.toFixed(2) + ' m:  ' +
        'mean atom error ' + xstats.planeMean.toFixed(2) + ' px over ' +
        xstats.planeN + ' atoms — parallax no 2D matrix can transfer\n' +
      '  depth-true transfer of ' + name + ' (back-project at z₁, reproject): ' +
        (isNaN(xstats.depthSel) ? 'n/a' :
          xstats.depthSel.toFixed(2) + ' px — only quantization remains') + '\n' +
      '  display diff |actual − predicted| (diff panel): mean ' +
        (xstats.diff && !isNaN(xstats.diff.mean) ?
          xstats.diff.mean.toFixed(1) + ' / max ' +
          xstats.diff.max.toFixed(0) + '  (of 255)' : 'n/a') +
        ' — bright pixels are the homography’s failures\n' +
      '\n' +
      'DLT triangulation from the two display pixels (doc §7):\n' +
      (Xhat ?
        '  X̂ = (' + Xhat.map(function (x) { return x.toFixed(3); }).join(', ') +
          ') m   error ‖X̂−X‖ = ' + (err * 1000).toFixed(1) + ' mm' +
          '   RMS reprojection = ' + rms.toFixed(3) + ' px\n' :
        '  (degenerate)\n') +
      '\n' +
      'rig geometry (doc §8): baseline B = ' + B.toFixed(2) + ' m,  depth Z ≈ ' +
        Z.toFixed(2) + ' m,  f = ' + FOCAL + ' px\n' +
      '  predicted σ_Z = Z²σ_d/(fB) = ' + (sZ * 1000).toFixed(1) +
        ' mm  for σ_d = ' + BASE_SIGMA_D.toFixed(2) + ' px (pixel quantization)';
  }

  /* ---------------- main render ------------------------------------------ */

  function render() {
    /* auto-aim: a camera with tracking on always points at the molecule,
     * so dragging its position never loses the subject */
    ['1', '2'].forEach(function (n) {
      if (!state.autoAim[n]) return;
      var c = state['cam' + n];
      var a = MV.aimAngles([c.x, c.y, c.z], molCenter());
      c.pan = a.pan; c.tilt = a.tilt;
    });
    var cam1 = camera(state.cam1), cam2 = camera(state.cam2);
    var atoms = MV.poseMolecule(mol, state.pose);
    var F = MV.fundamental(cam1, cam2);
    var F21 = MV.fundamental(cam2, cam1);
    drawOverview(atoms, cam1, cam2);
    drawDisplay('view1', cam1, atoms, F21, cam2, atoms);
    drawDisplay('view2', cam2, atoms, F, cam1, atoms);
    var xstats = drawTransfer(cam1, cam2, atoms);
    xstats.diff = drawDiff();
    pushDiff(xstats.diff.mean);
    drawErrGraph(null);
    updateReadout(cam1, cam2, atoms, F, xstats);
    syncSliders();                            /* keep fine controls honest */
  }

  /* ---------------- controls --------------------------------------------- */

  var clamp = function (x, lo, hi) { return Math.max(lo, Math.min(hi, x)); };

  function bindSlider(id, obj, key, scale) {
    var el = document.getElementById(id);
    el.addEventListener('input', function () {
      obj[key] = parseFloat(el.value) * (scale || 1);
      render();
    });
  }

  ['1', '2'].forEach(function (n) {
    var c = state['cam' + n];
    bindSlider('cam' + n + 'x', c, 'x');
    bindSlider('cam' + n + 'y', c, 'y');
    bindSlider('cam' + n + 'z', c, 'z');
    bindSlider('cam' + n + 'pan', c, 'pan', DEG);
    bindSlider('cam' + n + 'tilt', c, 'tilt', DEG);
    bindSlider('cam' + n + 'roll', c, 'roll', DEG);
    document.getElementById('aim' + n).addEventListener('click', function () {
      var a = MV.aimAngles([c.x, c.y, c.z], molCenter());
      c.pan = a.pan; c.tilt = a.tilt;
      render();
    });
    document.getElementById('aimauto' + n).addEventListener('change', function (e) {
      state.autoAim[n] = e.target.checked;
      render();
    });
  });
  bindSlider('molyaw', state.pose, 'yaw', DEG);
  bindSlider('molpitch', state.pose, 'pitch', DEG);
  bindSlider('molroll', state.pose, 'roll', DEG);
  bindSlider('molx', state.pose, 'x');
  bindSlider('moly', state.pose, 'y');
  bindSlider('molz', state.pose, 'z');

  /* drag on a camera view turns that camera (first-person pan/tilt);
   * the scroll wheel rolls it about its optical axis */
  ['view1', 'view2'].forEach(function (id, idx) {
    var cv = document.getElementById(id), drag = null;
    var n = idx === 0 ? '1' : '2';
    cv.addEventListener('mousedown', function (e) {
      drag = { x: e.clientX, y: e.clientY };
    });
    window.addEventListener('mousemove', function (e) {
      if (!drag) return;
      if (e.buttons === 0) { drag = null; return; }  /* lost mouseup */
      var c = state['cam' + n];
      state.autoAim[n] = false;               /* manual look = tracking off */
      document.getElementById('aimauto' + n).checked = false;
      c.pan -= (e.clientX - drag.x) * 0.004;   /* drag right = look right */
      c.tilt = clamp(c.tilt - (e.clientY - drag.y) * 0.004, -1.5, 1.5);
      drag = { x: e.clientX, y: e.clientY };
      render();
    });
    window.addEventListener('mouseup', function () { drag = null; });
    cv.addEventListener('wheel', function (e) {
      e.preventDefault();
      state['cam' + n].roll += e.deltaY * 0.002;
      render();
    }, { passive: false });
  });

  /* drag things in the world view. Grabbing a camera dot moves that
   * camera (view-plane; alt = depth y); grabbing anywhere else moves the
   * molecule (shift = rotate, alt = depth y). */
  (function () {
    var cv = document.getElementById('overview'), drag = null;
    cv.addEventListener('mousedown', function (e) {
      var rect = cv.getBoundingClientRect();
      var mx = e.clientX - rect.left, my = e.clientY - rect.top;
      var target = 'mol';
      ['1', '2'].forEach(function (n) {
        var c = state['cam' + n];
        var p = oproj([c.x, c.y, c.z]);
        if (Math.hypot(p[0] - mx, p[1] - my) < 12) target = 'cam' + n;
      });
      drag = { x: e.clientX, y: e.clientY, target: target };
      e.preventDefault();
    });
    window.addEventListener('mousemove', function (e) {
      if (!drag) return;
      if (e.buttons === 0) { drag = null; return; }  /* lost mouseup */
      var dx = e.clientX - drag.x, dy = e.clientY - drag.y;
      if (drag.target === 'mol') {
        var p = state.pose;
        if (e.shiftKey) {
          p.yaw += dx * 0.01;
          p.pitch = clamp(p.pitch + dy * 0.01, -Math.PI / 2, Math.PI / 2);
        } else if (e.altKey) {
          p.y = clamp(p.y + dx / OSCALE / -0.45, -0.8, 0.8);
        } else {
          p.x = clamp(p.x + dx / OSCALE, -0.8, 0.8);
          p.z = clamp(p.z - dy / OSCALE / 0.9, -0.8, 0.8);
        }
      } else {
        var c = state[drag.target];
        if (e.altKey) {
          c.y = clamp(c.y + dx / OSCALE / -0.45, -6, 6);
        } else {
          c.x = clamp(c.x + dx / OSCALE, -6, 6);
          c.z = clamp(c.z - dy / OSCALE / 0.9, -3, 3);
        }
      }
      drag = { x: e.clientX, y: e.clientY, target: drag.target };
      render();
    });
    window.addEventListener('mouseup', function () { drag = null; });
  })();

  /* click an atom (in either view) to select it */
  ['view1', 'view2'].forEach(function (id, idx) {
    var cv = document.getElementById(id), down = null;
    cv.addEventListener('mousedown', function (e) { down = { x: e.clientX, y: e.clientY }; });
    cv.addEventListener('click', function (e) {
      /* ignore clicks that were drags */
      if (down && (Math.abs(e.clientX - down.x) > 3 || Math.abs(e.clientY - down.y) > 3))
        return;
      var rect = cv.getBoundingClientRect();
      var u = (e.clientX - rect.left) / ZOOM, v = (e.clientY - rect.top) / ZOOM;
      var cam = camera(state[idx === 0 ? 'cam1' : 'cam2']);
      var atoms = MV.poseMolecule(mol, state.pose);
      var best = -1, bd = 25;
      atoms.forEach(function (X, i) {
        var p = MV.project(cam, X);
        if (!p) return;
        var d = (p.u - u) * (p.u - u) + (p.v - v) * (p.v - v);
        if (d < bd) { bd = d; best = i; }
      });
      if (best >= 0) { state.selected = best; render(); }
    });
  });

  function setVal(id, v) { document.getElementById(id).value = v; }
  /* --- URL-addressable scenes (doc section 1.1) ------------------------ */
  /* ?scene=<preset> applies a named configuration; any individual
   * parameter (slider ids, angles in degrees; sel=<atom>; aim1/aim2=0|1)
   * overrides on top. The document references these by name. */
  var PRESETS = {
    'default': {},
    /* cameras nearly together: plane transfer becomes exact, DLT
     * triangulation degrades (doc eq. parallax vs eq. deptherr) */
    narrow: { cam1x: 3.95, cam1y: -0.25, cam1z: 0.60,
              cam2x: 3.95, cam2y: 0.25, cam2z: 0.55, aim1: 1, aim2: 1 },
    /* camera 2 rolled: its epipolar lines slant */
    roll: { cam2roll: 35 }
  };

  function applyParams() {
    var q = new URLSearchParams(window.location.search);
    var vals = {};
    var preset = PRESETS[q.get('scene') || ''];
    if (preset) Object.keys(preset).forEach(function (k) { vals[k] = preset[k]; });
    q.forEach(function (v, k) { if (k !== 'scene') vals[k] = parseFloat(v); });
    Object.keys(vals).forEach(function (k) {
      var v = vals[k];
      if (isNaN(v)) return;
      var m = k.match(/^cam([12])(x|y|z|pan|tilt|roll)$/);
      if (m) {
        var c = state['cam' + m[1]];
        if (m[2] === 'pan' || m[2] === 'tilt') {
          c[m[2]] = v * DEG;
          state.autoAim[m[1]] = false;        /* explicit look direction */
        } else if (m[2] === 'roll') {
          c.roll = v * DEG;
        } else {
          c[m[2]] = v;
        }
        return;
      }
      var mm = k.match(/^mol(yaw|pitch|roll|x|y|z)$/);
      if (mm) {
        var key = mm[1];
        state.pose[key] =
          (key === 'yaw' || key === 'pitch' || key === 'roll') ? v * DEG : v;
        return;
      }
      if (k === 'sel')
        state.selected = Math.max(0, Math.min(mol.atoms.length - 1, Math.round(v)));
      if (k === 'aim1') state.autoAim[1] = !!v;
      if (k === 'aim2') state.autoAim[2] = !!v;
    });
  }

  function sceneLink() {
    var q = [];
    ['1', '2'].forEach(function (n) {
      var c = state['cam' + n];
      q.push('cam' + n + 'x=' + c.x.toFixed(2),
             'cam' + n + 'y=' + c.y.toFixed(2),
             'cam' + n + 'z=' + c.z.toFixed(2),
             'cam' + n + 'roll=' + (c.roll / DEG).toFixed(1),
             'aim' + n + '=' + (state.autoAim[n] ? 1 : 0));
      if (!state.autoAim[n])
        q.push('cam' + n + 'pan=' + (c.pan / DEG).toFixed(1),
               'cam' + n + 'tilt=' + (c.tilt / DEG).toFixed(1));
    });
    var p = state.pose;
    q.push('molyaw=' + (p.yaw / DEG).toFixed(1),
           'molpitch=' + (p.pitch / DEG).toFixed(1),
           'molroll=' + (p.roll / DEG).toFixed(1),
           'molx=' + p.x.toFixed(2), 'moly=' + p.y.toFixed(2),
           'molz=' + p.z.toFixed(2), 'sel=' + state.selected);
    return window.location.href.split('?')[0] + '?' + q.join('&');
  }

  document.getElementById('copylink').addEventListener('click', function () {
    var btn = this, url = sceneLink();
    function done(ok) {
      btn.textContent = ok ? 'copied ✓' : 'copy failed';
      setTimeout(function () { btn.textContent = 'copy scene link'; }, 1600);
    }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(url).then(function () { done(true); },
                                              function () { window.prompt('Scene link:', url); done(true); });
    } else {
      window.prompt('Scene link:', url);
      done(true);
    }
  });

  function syncSliders() {
    ['1', '2'].forEach(function (n) {
      var c = state['cam' + n];
      setVal('cam' + n + 'x', c.x);
      setVal('cam' + n + 'y', c.y);
      setVal('cam' + n + 'z', c.z);
      setVal('cam' + n + 'pan', wrapDeg(c.pan / DEG));
      setVal('cam' + n + 'tilt', c.tilt / DEG);
      setVal('cam' + n + 'roll', wrapDeg(c.roll / DEG));
    });
    setVal('molyaw', wrapDeg(state.pose.yaw / DEG));
    setVal('molpitch', state.pose.pitch / DEG);
    setVal('molroll', wrapDeg(state.pose.roll / DEG));
    setVal('molx', state.pose.x);
    setVal('moly', state.pose.y);
    setVal('molz', state.pose.z);
  }
  function wrapDeg(d) {
    while (d > 180) d -= 360;
    while (d < -180) d += 360;
    return d;
  }

  applyParams();
  document.getElementById('aimauto1').checked = state.autoAim[1];
  document.getElementById('aimauto2').checked = state.autoAim[2];
  syncSliders();
  render();
})();
