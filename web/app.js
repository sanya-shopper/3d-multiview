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
 */
'use strict';

(function () {
  var DISP_W = 64, DISP_H = 48, ZOOM = 6;   /* coarse displays, magnified */
  var FOCAL = 52;                            /* pixels, on the 64x48 sensor */
  var BASE_SIGMA_D = Math.sqrt(2) * (1 / Math.sqrt(12)); /* quantization */
  var DEG = Math.PI / 180;

  var mol = MV.makeMolecule();
  var state = {
    /* poses chosen to match the old defaults: both cameras ~4.2 m out,
     * aimed at the origin */
    cam1: { x: 3.72, y: -1.80, z: 0.75, pan: 154 * DEG, tilt: -10 * DEG, roll: 0 },
    cam2: { x: 3.76, y: 1.82, z: 0.42, pan: -154 * DEG, tilt: -6 * DEG, roll: 0 },
    pose: { yaw: 0.5, pitch: 0.3, roll: 0, x: 0, y: 0, z: 0 },
    selected: 6                               /* the O atom */
  };

  function camera(c) {
    return MV.makeCameraPose({ pos: [c.x, c.y, c.z],
                               pan: c.pan, tilt: c.tilt, roll: c.roll,
                               f: FOCAL, W: DISP_W, H: DISP_H });
  }
  function molCenter() {
    return [state.pose.x, state.pose.y, state.pose.z];
  }

  /* ---------------- coarse pixel displays ------------------------------- */

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

    /* magnify with hard pixel edges: the display, at the pixel level */
    cx.imageSmoothingEnabled = false;
    cx.clearRect(0, 0, cv.width, cv.height);
    cx.drawImage(buf, 0, 0, cv.width, cv.height);

    /* faint pixel grid so the quantization is visible */
    cx.strokeStyle = 'rgba(255,255,255,0.07)';
    cx.lineWidth = 1;
    var i;
    for (i = 0; i <= DISP_W; i++) {
      cx.beginPath(); cx.moveTo(i * ZOOM, 0); cx.lineTo(i * ZOOM, cv.height); cx.stroke();
    }
    for (i = 0; i <= DISP_H; i++) {
      cx.beginPath(); cx.moveTo(0, i * ZOOM); cx.lineTo(cv.width, i * ZOOM); cx.stroke();
    }

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
      drawLine(cx, l, cv.width, cv.height);
    }
  }

  /* draw a u + b v + c = 0 clipped to the display, in display pixels */
  function drawLine(cx, l, w, h) {
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

  /* ---------------- world overview -------------------------------------- */

  var OSCALE = 52;
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
      cx.arc(p[0], p[1], i === state.selected ? 7 : 5, 0, 2 * Math.PI);
      cx.fill();
      if (i === state.selected) {
        cx.strokeStyle = '#ffd75e'; cx.lineWidth = 2;
        cx.stroke();
      }
    });
    /* cameras: center dot, true optical axis (row 3 of R), label */
    [[cam1, 'camera 1'], [cam2, 'camera 2']].forEach(function (cc) {
      var cam = cc[0], p = oproj(cam.C);
      var tip = oproj(MV.add(cam.C, MV.scale(cam.R[2], 0.9)));
      cx.strokeStyle = '#2d6cdf'; cx.lineWidth = 1.5;
      cx.beginPath(); cx.moveTo(p[0], p[1]); cx.lineTo(tip[0], tip[1]); cx.stroke();
      cx.fillStyle = '#1c2733';
      cx.beginPath(); cx.arc(p[0], p[1], 5, 0, 2 * Math.PI); cx.fill();
      cx.font = '12px system-ui, sans-serif';
      cx.fillText(cc[1], p[0] + 8, p[1] + 4);
    });
  }

  /* ---------------- readout ---------------------------------------------- */

  function fmt(x, d) { return (x >= 0 ? ' ' : '') + x.toFixed(d); }
  function fmtMat(M, d) {
    return M.map(function (r) {
      return r.map(function (x) { return fmt(x, d); }).join('  ');
    }).join('\n');
  }

  function updateReadout(cam1, cam2, atoms, F) {
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
    var cam1 = camera(state.cam1), cam2 = camera(state.cam2);
    var atoms = MV.poseMolecule(mol, state.pose);
    var F = MV.fundamental(cam1, cam2);
    var F21 = MV.fundamental(cam2, cam1);
    drawOverview(atoms, cam1, cam2);
    drawDisplay('view1', cam1, atoms, F21, cam2, atoms);
    drawDisplay('view2', cam2, atoms, F, cam1, atoms);
    updateReadout(cam1, cam2, atoms, F);
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
      syncSliders(); render();
    });
  });
  bindSlider('molyaw', state.pose, 'yaw', DEG);
  bindSlider('molpitch', state.pose, 'pitch', DEG);
  bindSlider('molroll', state.pose, 'roll', DEG);
  bindSlider('molx', state.pose, 'x');
  bindSlider('moly', state.pose, 'y');
  bindSlider('molz', state.pose, 'z');

  /* drag on a camera view turns that camera (first-person pan/tilt) */
  ['view1', 'view2'].forEach(function (id, idx) {
    var cv = document.getElementById(id), drag = null;
    cv.addEventListener('mousedown', function (e) {
      drag = { x: e.clientX, y: e.clientY };
    });
    window.addEventListener('mousemove', function (e) {
      if (!drag) return;
      var c = state[idx === 0 ? 'cam1' : 'cam2'];
      c.pan -= (e.clientX - drag.x) * 0.004;   /* drag right = look right */
      c.tilt = clamp(c.tilt - (e.clientY - drag.y) * 0.004, -1.5, 1.5);
      drag = { x: e.clientX, y: e.clientY };
      syncSliders(); render();
    });
    window.addEventListener('mouseup', function () { drag = null; });
  });

  /* drag the molecule in the world view: translate in the view plane;
   * shift-drag rotates (yaw/pitch); alt-drag translates in y (depth) */
  (function () {
    var cv = document.getElementById('overview'), drag = null;
    cv.addEventListener('mousedown', function (e) {
      drag = { x: e.clientX, y: e.clientY };
      e.preventDefault();
    });
    window.addEventListener('mousemove', function (e) {
      if (!drag) return;
      var dx = e.clientX - drag.x, dy = e.clientY - drag.y;
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
      drag = { x: e.clientX, y: e.clientY };
      syncSliders(); render();
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

  syncSliders();
  render();
})();
