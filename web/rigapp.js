/* rigapp.js -- UI for "earn the rig": solve two static cameras from a
 * parallelepiped of known dimensions, then measure with the solved rig.
 * All mathematics lives in model.js; this file draws and routes input.
 *
 * Phase 1 (solve): the box is posed in the volume; each banked pose adds
 * its 20 mark correspondences (8 corners + 12 edge midpoints); the rig
 * is solved from ALL banked, deduplicated correspondences (calibrated
 * 8-point on normalized coordinates -> essential matrix -> chirality ->
 * R,t at unit baseline) and the known 1.20 m edge anchors the scale. The
 * readout compares solved baseline/rotation and the OTHER box dimensions
 * against ground truth.
 *
 * Phase 2 (measure): the box moves freely; its corners are triangulated
 * with the SOLVED rig (never the true one), accumulating a point cloud
 * and a scatter-updated TSDF (doc sections 3.5, 9.2).
 */
'use strict';

(function () {
  var DISP_W = 160, DISP_H = 120, ZOOM = 2;
  var FOCAL = 130;
  var DEG = Math.PI / 180;
  var OSCALE = 36;
  var I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

  /* ---- the true rig (static; the page's "reality") -------------------- */
  function fixedCam(pos) {
    var a = MV.aimAngles(pos, [0, 0, 0]);
    return MV.makeCameraPose({ pos: pos, pan: a.pan, tilt: a.tilt, roll: 0,
                               f: FOCAL, W: DISP_W, H: DISP_H });
  }
  var cam1 = fixedCam([3.72, -1.80, 0.75]);
  var cam2 = fixedCam([3.76, 1.82, 0.42]);
  var Btrue = MV.norm(MV.sub(cam2.C, cam1.C));
  var RtrueRel = MV.matMul(cam2.R, MV.transpose(cam1.R));

  /* ---- the box of known dimensions ------------------------------------ */
  var box = MV.makeBox(1.20, 0.80, 0.50);

  /* six faces, each its own color so every view shows which side is
   * which; corner indices wind around each face quad */
  var FACES = [
    { name: 'orange', color: '230,85,13',  axis: 0, sign: 1,  bit: 1 },
    { name: 'purple', color: '117,107,177', axis: 0, sign: -1, bit: 1 },
    { name: 'green',  color: '49,163,84',  axis: 1, sign: 1,  bit: 2 },
    { name: 'red',    color: '214,39,40',  axis: 1, sign: -1, bit: 2 },
    { name: 'blue',   color: '31,120,180', axis: 2, sign: 1,  bit: 4 },
    { name: 'gold',   color: '181,137,0',  axis: 2, sign: -1, bit: 4 }
  ];
  FACES.forEach(function (f) {
    var others = [1, 2, 4].filter(function (b) { return b !== f.bit; });
    var base = f.sign > 0 ? f.bit : 0;
    f.corners = [base, base | others[0], base | others[0] | others[1],
                 base | others[1]];
  });
  function faceGeom(f, atoms) {
    var c = [0, 0, 0];
    f.corners.forEach(function (i) { c = MV.add(c, atoms[i]); });
    c = MV.scale(c, 0.25);
    var n = [0, 0, 0]; n[f.axis] = f.sign;
    return { center: c, normal: MV.matVec(poseR(state.pose), n) };
  }
  var state = {
    pose: { yaw: 25 * DEG, pitch: 15 * DEG, roll: 0, x: 0, y: 0, z: 0 },
    phase: 'solve',                           /* 'solve' | 'measure' */
    subpixel: false,
    autoBank: true,                           /* bank on release, if novel */
    bankMsg: 'drag the box, release to auto-bank',
    banks: [],                                /* each: {obs, pose, R} */
    solved: null,                             /* {F, R, t, s, camS1, camS2} */
    baselineHist: [],                         /* % error per bank */
    faceCov: [0, 0, 0, 0, 0, 0],              /* measure-phase coverage */
    prevS: null, stableCount: 0,
    solveStatus: 'none',                      /* none | rough | solved */
    cloud: [],                                /* world-frame points */
    tsdf: null,
    sliceK: 16,
    lastPoseKey: ''
  };

  /* solved-frame (camera-1 frame) -> true world, for display only */
  function toWorld(X) {
    return MV.add(MV.matVec(MV.transpose(cam1.R), X), cam1.C);
  }

  /* the box's detectable marks: 8 corners + 12 printed edge-midpoint
   * dots (the midpoints break the near-degenerate configuration that
   * bare parallelepiped corners present to the two-view solve) */
  function markers(atoms) {
    var pts = atoms.slice();
    box.bonds.forEach(function (b) {
      pts.push(MV.scale(MV.add(atoms[b[0]], atoms[b[1]]), 0.5));
    });
    return pts;
  }

  /* detect the box marks on both displays under the current pose */
  function observe() {
    var atoms = MV.poseMolecule(box, state.pose);
    var out = [];
    markers(atoms).forEach(function (X, i) {
      var p1 = MV.project(cam1, X), p2 = MV.project(cam2, X);
      if (!p1 || !p2) return;
      if (p1.u < 0 || p1.u >= DISP_W || p1.v < 0 || p1.v >= DISP_H) return;
      if (p2.u < 0 || p2.u >= DISP_W || p2.v < 0 || p2.v >= DISP_H) return;
      if (!state.subpixel) { p1 = MV.quantize(p1); p2 = MV.quantize(p2); }
      out.push({ i: i, p1: { u: p1.u, v: p1.v }, p2: { u: p2.u, v: p2.v } });
    });
    return out;
  }

  /* rotation matrix of a box pose (same composition as poseMolecule) */
  function poseR(pose) {
    var cy = Math.cos(pose.yaw), sy = Math.sin(pose.yaw);
    var cp = Math.cos(pose.pitch), sp = Math.sin(pose.pitch);
    var cr = Math.cos(pose.roll), sr = Math.sin(pose.roll);
    return MV.matMul([[cr, -sr, 0], [sr, cr, 0], [0, 0, 1]],
      MV.matMul([[1, 0, 0], [0, cp, -sp], [0, sp, cp]],
                [[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]]));
  }

  /* how different the current pose is from the closest banked one:
   * >= 1 means "novel enough to add information" (15 deg or 15 cm) */
  function novelty(pose) {
    if (!state.banks.length) return Infinity;
    var R = poseR(pose), c = [pose.x, pose.y, pose.z], best = Infinity;
    state.banks.forEach(function (b) {
      var d = MV.rotationAngle(R, b.R) / (15 * DEG) +
              MV.norm(MV.sub(c, [b.pose.x, b.pose.y, b.pose.z])) / 0.15;
      if (d < best) best = d;
    });
    return best;
  }

  /* bank the current pose (returns true if banked) */
  function bankPose(auto) {
    var obs = observe();
    if (obs.length < 8) {
      state.bankMsg = 'box partly out of view — not banked';
      return false;
    }
    if (auto && novelty(state.pose) < 1) {
      state.bankMsg = 'pose too similar to a banked one — move or rotate more';
      return false;
    }
    state.banks.push({ obs: obs,
                       pose: JSON.parse(JSON.stringify(state.pose)),
                       R: poseR(state.pose) });
    runSolve();
    var msg = (auto ? 'auto-' : '') + 'banked pose ' + state.banks.length;
    var h = state.baselineHist;
    if (h.length >= 2 && h[h.length - 1] > h[h.length - 2]) {
      msg += ' — error ticked UP: that pose\'s pixel noise pulled the fit; ' +
        'normal, it falls on average (best so far ' +
        Math.min.apply(null, h).toFixed(2) + '%)';
    }
    state.bankMsg = msg;
    return true;
  }

  /* ---- the solve ------------------------------------------------------- */

  function runSolve() {
    var all = [], seen = {};
    state.banks.forEach(function (b) { b.obs.forEach(function (c) {
      var key = c.p1.u + ',' + c.p1.v + '|' + c.p2.u + ',' + c.p2.v;
      if (seen[key]) return;                  /* quantization duplicates */
      seen[key] = 1;
      all.push(c);
    }); });
    if (all.length < 8) { state.solved = null; return; }
    var F = MV.essentialFromPairs(all, cam1.K, cam2.K);
    if (!F) { state.solved = null; return; }
    var rel = MV.relativePose(F, cam1.K, cam2.K, all, DISP_W, DISP_H);
    if (!rel) { state.solved = null; return; }
    var camU1 = MV.camFromRt(cam1.K, I3, [0, 0, 0], DISP_W, DISP_H);
    var camU2 = MV.camFromRt(cam2.K, rel.R, rel.t, DISP_W, DISP_H);
    /* metric scale from the KNOWN 1.20 m edges (bit-1 corner pairs) of
     * every banked pose, in the unit-baseline reconstruction */
    var ratios = [];
    state.banks.forEach(function (b) {
      var X = {};
      b.obs.forEach(function (c) {
        X[c.i] = MV.triangulate(camU1, c.p1, camU2, c.p2);
      });
      b.obs.forEach(function (c) {
        var j = c.i ^ 1;                      /* corner pairs only */
        if (c.i < 8 && j > c.i && X[c.i] && X[j])
          ratios.push(box.dims[0] / MV.norm(MV.sub(X[j], X[c.i])));
      });
    });
    if (!ratios.length) { state.solved = null; return; }
    var s = ratios.reduce(function (a, b) { return a + b; }, 0) / ratios.length;
    state.solved = {
      F: F, R: rel.R, t: rel.t, s: s,
      camS1: camU1,
      camS2: MV.camFromRt(cam2.K, rel.R, MV.scale(rel.t, s), DISP_W, DISP_H),
      rms: MV.epipolarRMS(F, all), npts: all.length
    };
    state.baselineHist.push(Math.abs(s - Btrue) / Btrue * 100);
    /* IN-SYSTEM convergence: the solver cannot see the truth, so
     * "solved" is declared from solve-to-solve stability */
    if (state.prevS !== null && Math.abs(s - state.prevS) / state.prevS < 0.01)
      state.stableCount++;
    else
      state.stableCount = 0;
    state.prevS = s;
    state.solveStatus =
      (state.banks.length >= 4 && state.stableCount >= 2) ? 'solved' : 'rough';
  }

  /* triangulate the currently visible corners with the SOLVED metric rig */
  function measureCorners(obs) {
    if (!state.solved) return {};
    var X = {};
    obs.forEach(function (c) {
      var p = MV.triangulate(state.solved.camS1, c.p1, state.solved.camS2, c.p2);
      if (p) X[c.i] = p;                      /* solved frame = cam1 frame */
    });
    return X;
  }

  /* mean measured edge length per box axis, from solved-frame corners */
  function measureDims(X) {
    var dims = [NaN, NaN, NaN];
    [1, 2, 4].forEach(function (bit, axis) {
      var acc = 0, n = 0, i;
      for (i = 0; i < 8; i++) {
        var j = i ^ bit;
        if (j > i && X[i] && X[j]) { acc += MV.norm(MV.sub(X[j], X[i])); n++; }
      }
      if (n) dims[axis] = acc / n;
    });
    return dims;
  }

  /* ---- drawing --------------------------------------------------------- */

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

  function drawView(canvasId, cam, atoms, obs, obsKey, reproj) {
    var cv = document.getElementById(canvasId), cx = cv.getContext('2d');
    var buf = document.createElement('canvas');
    buf.width = DISP_W; buf.height = DISP_H;
    var bx = buf.getContext('2d');
    bx.fillStyle = '#10141a'; bx.fillRect(0, 0, DISP_W, DISP_H);
    FACES.forEach(function (f) {
      var g = faceGeom(f, atoms);
      if (MV.dot(g.normal, MV.sub(g.center, cam.C)) >= 0) return; /* back */
      var pts = f.corners.map(function (i) { return MV.project(cam, atoms[i]); });
      if (pts.some(function (p) { return !p; })) return;
      bx.fillStyle = 'rgba(' + f.color + ',0.50)';
      bx.beginPath();
      pts.forEach(function (p, k) {
        if (k === 0) bx.moveTo(p.u, p.v); else bx.lineTo(p.u, p.v);
      });
      bx.closePath(); bx.fill();
    });
    bx.strokeStyle = '#9aa3ad'; bx.lineWidth = 1;
    box.bonds.forEach(function (b) {
      var p = MV.project(cam, atoms[b[0]]), q = MV.project(cam, atoms[b[1]]);
      if (!p || !q) return;
      bx.beginPath(); bx.moveTo(p.u, p.v); bx.lineTo(q.u, q.v); bx.stroke();
    });
    markers(atoms).forEach(function (X, i) {
      var p = MV.project(cam, X);
      if (!p) return;
      bx.fillStyle = i < 8 ? '#2d6cdf' : '#7aa2e8';
      bx.beginPath(); bx.arc(p.u, p.v, i < 8 ? 2.4 : 1.4, 0, 2 * Math.PI);
      bx.fill();
    });
    cx.imageSmoothingEnabled = false;
    cx.clearRect(0, 0, cv.width, cv.height);
    cx.drawImage(buf, 0, 0, cv.width, cv.height);
    /* detected mark pixels */
    obs.forEach(function (c) {
      var p = c[obsKey];
      cx.strokeStyle = '#ffd75e'; cx.lineWidth = 1;
      cx.strokeRect((p.u - 1) * ZOOM, (p.v - 1) * ZOOM, ZOOM * 2, ZOOM * 2);
    });
    /* measure phase: reprojections of the solved-rig triangulations */
    if (reproj) reproj.forEach(function (p) {
      cx.strokeStyle = '#5ec8ff'; cx.lineWidth = 1.5;
      cx.beginPath();
      cx.moveTo((p.u - 1) * ZOOM, p.v * ZOOM); cx.lineTo((p.u + 1) * ZOOM, p.v * ZOOM);
      cx.moveTo(p.u * ZOOM, (p.v - 1) * ZOOM); cx.lineTo(p.u * ZOOM, (p.v + 1) * ZOOM);
      cx.stroke();
    });
  }

  function oproj(X) {
    var cv = document.getElementById('overview');
    return [cv.width / 2 + OSCALE * (X[0] - 0.45 * X[1]),
            cv.height / 2 - OSCALE * (X[2] * 0.9 + 0.28 * X[1])];
  }

  function drawOverview(atoms) {
    var cv = document.getElementById('overview'), cx = cv.getContext('2d');
    cx.clearRect(0, 0, cv.width, cv.height);
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
    /* slice-height marker: the TSDF panel shows THIS dashed level */
    if (state.tsdf && state.phase === 'measure') {
      var zk = sliceZ(), sq = 1.1;
      cx.strokeStyle = '#b45309'; cx.lineWidth = 1;
      cx.setLineDash([4, 3]);
      cx.beginPath();
      [[-sq, -sq], [sq, -sq], [sq, sq], [-sq, sq]].forEach(function (c, k) {
        var p = oproj([c[0], c[1], zk]);
        if (k === 0) cx.moveTo(p[0], p[1]); else cx.lineTo(p[0], p[1]);
      });
      cx.closePath(); cx.stroke();
      cx.setLineDash([]);
    }
    /* the box, true pose: front faces filled in their colors */
    var VDIR = MV.unit([0.45, 1, -0.311]);    /* oblique-view kernel */
    FACES.forEach(function (f) {
      var g = faceGeom(f, atoms);
      if (MV.dot(g.normal, VDIR) >= 0) return;
      var pts = f.corners.map(function (i) { return oproj(atoms[i]); });
      cx.fillStyle = 'rgba(' + f.color + ',0.30)';
      cx.beginPath();
      pts.forEach(function (p, k) {
        if (k === 0) cx.moveTo(p[0], p[1]); else cx.lineTo(p[0], p[1]);
      });
      cx.closePath(); cx.fill();
    });
    cx.strokeStyle = '#6b7280'; cx.lineWidth = 2;
    box.bonds.forEach(function (b) {
      var p = oproj(atoms[b[0]]), q = oproj(atoms[b[1]]);
      cx.beginPath(); cx.moveTo(p[0], p[1]); cx.lineTo(q[0], q[1]); cx.stroke();
    });
    atoms.forEach(function (X) {
      var p = oproj(X);
      cx.fillStyle = '#2d6cdf';
      cx.beginPath(); cx.arc(p[0], p[1], 3.5, 0, 2 * Math.PI); cx.fill();
    });
    /* cameras with frusta */
    [[cam1, 'cam 1'], [cam2, 'cam 2']].forEach(function (cc) {
      var cam = cc[0], p = oproj(cam.C);
      var Ki = MV.inv3(cam.K), depth = 1.3;
      var far = [[0, 0], [DISP_W, 0], [DISP_W, DISP_H], [0, DISP_H]]
        .map(function (uv) {
          var dc = MV.matVec(Ki, [uv[0], uv[1], 1]);
          return oproj(MV.add(cam.C,
            MV.matVec(MV.transpose(cam.R), MV.scale(dc, depth))));
        });
      var i, j;
      for (i = 0; i < 4; i++) {
        j = (i + 1) % 4;
        cx.fillStyle = 'rgba(45,108,223,0.07)';
        cx.beginPath(); cx.moveTo(p[0], p[1]);
        cx.lineTo(far[i][0], far[i][1]); cx.lineTo(far[j][0], far[j][1]);
        cx.closePath(); cx.fill();
      }
      cx.strokeStyle = 'rgba(45,108,223,0.55)'; cx.lineWidth = 1;
      for (i = 0; i < 4; i++) {
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
    /* the SOLVED camera 2, ghosted, with its offset from truth */
    if (state.solved) {
      var C2s = toWorld(MV.scale(MV.matVec(MV.transpose(state.solved.R),
                                           MV.scale(state.solved.t, state.solved.s)), -1));
      var pg = oproj(C2s), pt = oproj(cam2.C);
      var ghostCol = state.solveStatus === 'solved' ? '#2e7d32' : '#b45309';
      cx.strokeStyle = ghostCol; cx.lineWidth = 1.5;
      cx.beginPath(); cx.arc(pg[0], pg[1], 6, 0, 2 * Math.PI); cx.stroke();
      cx.setLineDash([3, 3]);
      cx.beginPath(); cx.moveTo(pg[0], pg[1]); cx.lineTo(pt[0], pt[1]); cx.stroke();
      cx.setLineDash([]);
      cx.fillStyle = ghostCol; cx.font = '10px system-ui, sans-serif';
      cx.fillText('solved cam 2', pg[0] + 8, pg[1] - 6);
    }
    drawWorldAxes(cx, cv);
  }

  function drawWorldAxes(cx, cv) {
    var ax = 34, ay = cv.height - 30, len = 0.55;
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

  /* TSDF slice: diverging map, dark at the zero crossing (the surface).
   * Overlaid: the box's TRUE cross-section at this height (dashed) and
   * the directions the camera rays come from, so the streaks read. */
  function sliceZ() {
    return state.tsdf ?
      state.tsdf.min[2] + (state.sliceK + 0.5) * state.tsdf.cell : 0;
  }
  function boxCrossSection(zk, atoms) {
    var pts = [];
    box.bonds.forEach(function (b) {
      var A = atoms[b[0]], B = atoms[b[1]], dz = B[2] - A[2];
      if (Math.abs(dz) < 1e-9) return;
      var t = (zk - A[2]) / dz;
      if (t < 0 || t > 1) return;
      pts.push([A[0] + t * (B[0] - A[0]), A[1] + t * (B[1] - A[1])]);
    });
    if (pts.length < 3) return null;
    var cx0 = 0, cy0 = 0;
    pts.forEach(function (p) { cx0 += p[0]; cy0 += p[1]; });
    cx0 /= pts.length; cy0 /= pts.length;
    pts.sort(function (a, b) {
      return Math.atan2(a[1] - cy0, a[0] - cx0) -
             Math.atan2(b[1] - cy0, b[0] - cx0);
    });
    return pts;
  }
  function drawSlice() {
    var cv = document.getElementById('tsdfslice'), cx = cv.getContext('2d');
    cx.clearRect(0, 0, cv.width, cv.height);
    if (!state.tsdf) {
      cx.fillStyle = '#4a5563'; cx.font = '11px system-ui, sans-serif';
      cx.fillText('solve the rig, then move the box to fuse', 12, cv.height / 2);
      return;
    }
    var t = state.tsdf, n = t.n, px = cv.width / n;
    var s = t.slice(state.sliceK), i, j;
    function lerp(a, b, f) { return Math.round(a + (b - a) * f); }
    for (j = 0; j < n; j++) for (i = 0; i < n; i++) {
      var v = s[j][i], col;
      if (isNaN(v)) col = '#eceff3';
      else {
        var f = Math.max(-1, Math.min(1, v / t.tau));
        col = f >= 0
          ? 'rgb(' + lerp(30, 219, f) + ',' + lerp(70, 231, f) + ',' + lerp(160, 248, f) + ')'
          : 'rgb(' + lerp(150, 248, -f) + ',' + lerp(85, 227, -f) + ',' + lerp(20, 207, -f) + ')';
      }
      cx.fillStyle = col;
      /* screen y up = world +y: flip rows */
      cx.fillRect(i * px, (n - 1 - j) * px, Math.ceil(px), Math.ceil(px));
    }
    /* world (x,y) -> slice-panel pixels, matching the voxel fill above */
    var size = t.cell * n;
    function sp(x, y) {
      return [(x - t.min[0]) / size * cv.width,
              cv.height - (y - t.min[1]) / size * cv.height];
    }
    /* the box's true cross-section at this height, for orientation */
    var zk = sliceZ();
    var atoms = MV.poseMolecule(box, state.pose);
    var poly = boxCrossSection(zk, atoms);
    if (poly) {
      cx.strokeStyle = '#1c2733'; cx.lineWidth = 1.5;
      cx.setLineDash([5, 3]);
      cx.beginPath();
      poly.forEach(function (p, k) {
        var q = sp(p[0], p[1]);
        if (k === 0) cx.moveTo(q[0], q[1]); else cx.lineTo(q[0], q[1]);
      });
      cx.closePath(); cx.stroke();
      cx.setLineDash([]);
    } else {
      cx.fillStyle = '#4a5563'; cx.font = '10px system-ui, sans-serif';
      cx.fillText('box does not reach this height', 6, cv.height - 6);
    }
    /* where the camera rays come from */
    [[cam1, 'cam 1'], [cam2, 'cam 2']].forEach(function (cc) {
      var C = cc[0].C;
      var d = MV.unit([C[0], C[1], 0]);
      var bx0 = cv.width / 2 + d[0] * (cv.width / 2 - 16);
      var by0 = cv.height / 2 - d[1] * (cv.height / 2 - 16);
      cx.strokeStyle = '#2d6cdf'; cx.lineWidth = 1.5;
      cx.beginPath(); cx.moveTo(bx0, by0);
      cx.lineTo(bx0 - d[0] * 10, by0 + d[1] * 10); cx.stroke();
      cx.fillStyle = '#2d6cdf'; cx.font = '9px system-ui, sans-serif';
      cx.fillText(cc[1], bx0 - 12, by0 + (d[1] > 0 ? 12 : -5));
    });
    document.getElementById('slicelabel').textContent =
      'z = ' + zk.toFixed(2) + ' m';
  }

  /* dedicated point-cloud panel: the world view stays clean */
  function drawCloud(atoms) {
    var cv = document.getElementById('cloudview'), cx = cv.getContext('2d');
    var S = 58;
    function cp(X) {
      return [cv.width / 2 + S * (X[0] - 0.45 * X[1]),
              cv.height / 2 - S * (X[2] * 0.9 + 0.28 * X[1])];
    }
    cx.clearRect(0, 0, cv.width, cv.height);
    /* faint volume outline for orientation */
    cx.strokeStyle = '#e2e6ea'; cx.lineWidth = 1;
    var s = 1.1, k, corners = [], edges = [
      [0, 1], [1, 3], [3, 2], [2, 0], [4, 5], [5, 7], [7, 6], [6, 4],
      [0, 4], [1, 5], [2, 6], [3, 7]];
    for (k = 0; k < 8; k++)
      corners.push([(k & 1 ? s : -s), (k & 2 ? s : -s), (k & 4 ? s : -s)]);
    edges.forEach(function (e) {
      var p = cp(corners[e[0]]), q = cp(corners[e[1]]);
      cx.beginPath(); cx.moveTo(p[0], p[1]); cx.lineTo(q[0], q[1]); cx.stroke();
    });
    if (!state.cloud.length) {
      cx.fillStyle = '#4a5563'; cx.font = '10px system-ui, sans-serif';
      cx.fillText('cloud accumulates in phase 2', 10, cv.height / 2);
      return;
    }
    /* ghost of the box's current true pose, for reference */
    cx.strokeStyle = '#c9ced4'; cx.lineWidth = 1;
    box.bonds.forEach(function (b) {
      var p = cp(atoms[b[0]]), q = cp(atoms[b[1]]);
      cx.beginPath(); cx.moveTo(p[0], p[1]); cx.lineTo(q[0], q[1]); cx.stroke();
    });
    cx.fillStyle = 'rgba(15,118,110,0.6)';
    state.cloud.forEach(function (X) {
      var p = cp(X);
      cx.fillRect(p[0] - 1, p[1] - 1, 2, 2);
    });
    cx.fillStyle = '#4a5563'; cx.font = '10px system-ui, sans-serif';
    cx.fillText(state.cloud.length + ' pts', 6, 12);
  }

  /* face-coverage meter: how clearly each colored side has been seen */
  var COV_TARGET = 8;                         /* ~8 frontal observations */
  function drawFaceCov() {
    var cv = document.getElementById('facecov'), cx = cv.getContext('2d');
    cx.clearRect(0, 0, cv.width, cv.height);
    if (state.phase !== 'measure') {
      cx.fillStyle = '#4a5563'; cx.font = '10px system-ui, sans-serif';
      cx.fillText('face coverage appears in phase 2', 6, cv.height / 2);
      return;
    }
    var rowH = cv.height / 6;
    FACES.forEach(function (f, i) {
      var y = i * rowH + 2;
      cx.fillStyle = 'rgb(' + f.color + ')';
      cx.fillRect(2, y, 10, rowH - 4);
      var frac = Math.min(1, state.faceCov[i] / COV_TARGET);
      cx.fillStyle = '#eceff3';
      cx.fillRect(16, y, cv.width - 20, rowH - 4);
      cx.fillStyle = 'rgba(' + f.color + ',0.75)';
      cx.fillRect(16, y, (cv.width - 20) * frac, rowH - 4);
      if (frac >= 1) {
        cx.fillStyle = '#1c2733'; cx.font = '9px system-ui, sans-serif';
        cx.fillText('✓', cv.width - 12, y + rowH - 7);
      }
    });
  }

  /* convergence graph: baseline error % after each banked pose */
  function drawGraph() {
    var cv = document.getElementById('convgraph'), cx = cv.getContext('2d');
    var W = cv.width, H = cv.height;
    var padL = 6, padR = 40, padT = 12, padB = 14;
    cx.clearRect(0, 0, W, H);
    var h = state.baselineHist;
    if (!h.length) {
      cx.fillStyle = '#4a5563'; cx.font = '10px system-ui, sans-serif';
      cx.fillText('bank poses to accumulate solves', padL + 2, H / 2);
      return;
    }
    var top = Math.max(2, Math.max.apply(null, h) * 1.15);
    function px(i) { return padL + (W - padL - padR) * i / Math.max(9, h.length - 1); }
    function py(v) { return padT + (H - padT - padB) * (1 - v / top); }
    cx.strokeStyle = '#d7dce2'; cx.lineWidth = 1;
    cx.fillStyle = '#4a5563'; cx.font = '9px system-ui, sans-serif';
    [top, top / 2].forEach(function (v) {
      cx.beginPath(); cx.moveTo(padL, py(v)); cx.lineTo(W - padR, py(v)); cx.stroke();
      cx.fillText(v.toFixed(1) + '%', padL + 1, py(v) - 2);
    });
    cx.strokeStyle = '#2d6cdf'; cx.lineWidth = 2; cx.lineJoin = 'round';
    cx.beginPath();
    h.forEach(function (v, i) {
      if (i === 0) cx.moveTo(px(i), py(v)); else cx.lineTo(px(i), py(v));
    });
    cx.stroke();
    cx.fillStyle = '#2d6cdf';
    h.forEach(function (v, i) {
      cx.beginPath(); cx.arc(px(i), py(v), 3, 0, 2 * Math.PI); cx.fill();
    });
    var best = Math.min.apply(null, h);
    cx.strokeStyle = '#16a085'; cx.lineWidth = 1;
    cx.setLineDash([4, 3]);
    cx.beginPath(); cx.moveTo(padL, py(best)); cx.lineTo(W - padR, py(best)); cx.stroke();
    cx.setLineDash([]);
    cx.fillStyle = '#16a085'; cx.font = '9px system-ui, sans-serif';
    cx.fillText('best ' + best.toFixed(2) + '%', W - padR + 2, py(best) + 3);
    var last = h[h.length - 1];
    cx.fillStyle = '#1c2733'; cx.font = '11px system-ui, sans-serif';
    cx.fillText('baseline err ' + last.toFixed(2) + '%',
                Math.min(px(h.length - 1) + 6, W - padR - 60),
                Math.max(padT + 8, py(last) - 6));
    cx.fillStyle = '#4a5563'; cx.font = '9px system-ui, sans-serif';
    cx.fillText('solves (one per banked pose)', padL, H - 3);
  }

  /* ---- readout ---------------------------------------------------------- */

  function fmtDims(d) {
    return d.map(function (x) { return isNaN(x) ? '—' : (x * 1000).toFixed(0); })
      .join(' × ') + ' mm';
  }

  function updateReadout(obs) {
    var out = document.getElementById('readout');
    var lines = [];
    var trueDims = 'true box: ' + fmtDims(box.dims) +
      '   (the 1200 mm edge is the known anchor)';
    if (!state.solved) {
      lines.push('phase 1 — SOLVE THE RIG (doc §6, §1: a known object anchors the unit)');
      lines.push('');
      lines.push(trueDims);
      lines.push('');
      lines.push('banked poses: ' + state.banks.length +
        '   correspondences: ' + state.banks.reduce(function (a, b) { return a + b.obs.length; }, 0));
      lines.push('');
      lines.push('bank a pose to run the first solve: the 8-point machinery on');
      lines.push('the 20 mark correspondences (§6) → essential matrix → (R, t) at');
      lines.push('unit baseline →');
      lines.push('the known edge sets the metre. Vary the box pose between banks;');
      lines.push('more and more-varied poses average the pixel-quantization noise.');
    } else {
      var sol = state.solved;
      var rotErr = MV.rotationAngle(sol.R, RtrueRel) / DEG;
      var X = measureCorners(obs);
      var dims = measureDims(X);
      var rmsPts = [], rms = NaN;
      var marksTrue = markers(MV.poseMolecule(box, state.pose));
      Object.keys(X).forEach(function (i) {
        var Xt = marksTrue[i];               /* all 20 marks, not 8 atoms */
        if (!Xt) return;
        rmsPts.push(MV.norm(MV.sub(toWorld(X[i]), Xt)));
      });
      if (rmsPts.length)
        rms = Math.sqrt(rmsPts.reduce(function (a, b) { return a + b * b; }, 0) / rmsPts.length);
      lines.push('phase ' + (state.phase === 'solve' ? '1 — SOLVE THE RIG' : '2 — MEASURE WITH THE SOLVED RIG'));
      lines.push('');
      lines.push(trueDims);
      lines.push('');
      lines.push('rig solve from ' + sol.npts + ' correspondences over ' +
        state.banks.length + ' pose(s)  (epipolar RMS ' + sol.rms.toFixed(2) + ' px):');
      /* computed vs actual poses and distances, world frame; camera 1
       * anchors the solved frame so its pose is identity by construction */
      var C2s = toWorld(MV.scale(MV.matVec(MV.transpose(sol.R),
        MV.scale(sol.t, sol.s)), -1));
      var posErr = MV.norm(MV.sub(C2s, cam2.C));
      var R2sW = MV.matMul(sol.R, cam1.R);
      function ptr(R) {
        var zc = R[2], xc = R[0];
        var pan = Math.atan2(zc[1], zc[0]) / DEG;
        var tilt = Math.asin(Math.max(-1, Math.min(1, zc[2]))) / DEG;
        var xc0 = MV.unit(MV.cross(zc, [0, 0, 1]));
        var yc0 = MV.cross(zc, xc0);
        var roll = Math.atan2(MV.dot(xc, yc0), MV.dot(xc, xc0)) / DEG;
        return pan.toFixed(1) + '/' + tilt.toFixed(1) + '/' + roll.toFixed(1) + '°';
      }
      function fmtVec(v) {
        return '(' + v.map(function (x) { return x.toFixed(2); }).join(', ') + ')';
      }
      lines.push('computed vs actual (camera 1 anchors the solved frame):');
      lines.push('  baseline |C₂−C₁|:  solved ' + (sol.s * 1000).toFixed(0) +
        ' mm   true ' + (Btrue * 1000).toFixed(0) + ' mm   (' +
        (Math.abs(sol.s - Btrue) / Btrue * 100).toFixed(2) + '% off)');
      lines.push('  camera 2 position: solved ' + fmtVec(C2s) + '   true ' +
        fmtVec(cam2.C) + ' m   off ' + (posErr * 1000).toFixed(0) + ' mm');
      lines.push('  camera 2 pan/tilt/roll: solved ' + ptr(R2sW) +
        '   true ' + ptr(cam2.R) + '   (3D error ' + rotErr.toFixed(2) + '°)');
      var cKeys = Object.keys(X).filter(function (i) { return +i < 8; });
      if (cKeys.length >= 4) {
        var cS = [0, 0, 0];
        cKeys.forEach(function (i) { cS = MV.add(cS, X[i]); });
        cS = MV.scale(cS, 1 / cKeys.length);
        var d1s = MV.norm(cS);                 /* cam 1 at solved origin */
        var d1t = MV.norm(MV.sub(
          [state.pose.x, state.pose.y, state.pose.z], cam1.C));
        lines.push('  camera 1 → box centre: solved ' + d1s.toFixed(3) +
          ' m   true ' + d1t.toFixed(3) + ' m');
      }
      lines.push('');
      lines.push('box measured through the solved rig (triangulated corners, §7):');
      lines.push('  estimated ' + fmtDims(dims) + '   vs true ' + fmtDims(box.dims));
      lines.push('  (the 1200 mm edge is the anchor — its agreement is by');
      lines.push('   construction; the other two dimensions are honest tests)');
      if (!isNaN(rms))
        lines.push('  corner position RMS vs truth: ' + (rms * 1000).toFixed(1) + ' mm');
      if (state.phase === 'solve') {
        lines.push('');
        coach(obs).forEach(function (l) { lines.push(l); });
      }
      if (state.phase === 'measure') {
        lines.push('');
        lines.push('model accumulation (§9): cloud ' + state.cloud.length +
          ' points; TSDF ' + (state.tsdf ? state.tsdf.observed() : 0) +
          ' observed voxels (scatter/ray updates, §3.5)');
        var unseen = [];
        FACES.forEach(function (f, i) {
          if (state.faceCov[i] < COV_TARGET * 0.5) unseen.push(f.name);
        });
        if (unseen.length)
          lines.push('faces not yet seen clearly: ' + unseen.join(', ') +
            ' — rotate the box to front them toward a camera');
        else
          lines.push('all six faces covered — the model has seen the whole box');
        lines.push('hold the box still and quantization noise averages down (§9.2);');
        lines.push('a moving box smears the TSDF — the reason the doc separates');
        lines.push('background from movers (§9.3–9.4).');
      }
    }
    out.textContent = lines.join('\n');
  }

  /* actionable advice for driving the solve error down */
  function coach(obs) {
    var lines = ['driving the error down:'];
    var n = state.banks.length, i, j;
    /* the estimate wanders: each bank adds a fresh draw of quantization
     * noise, so expect improvement on AVERAGE, not every bank */
    var h = state.baselineHist;
    if (h.length >= 2 && h[h.length - 1] > h[h.length - 2])
      lines.push('  · error went UP on that bank — normal: each pose adds a fresh');
    else
      lines.push('  · error falls only on average — each pose adds a fresh');
    lines.push('    random draw of quantization noise, so the estimate wanders');
    lines.push('    (~1/√N in expectation; doc §13.2 mindset)');
    /* rotation diversity */
    if (n >= 2) {
      var meanRot = 0, cnt = 0;
      for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) {
        meanRot += MV.rotationAngle(state.banks[i].R, state.banks[j].R);
        cnt++;
      }
      meanRot = meanRot / cnt / DEG;
      if (meanRot < 35)
        lines.push('  · banked orientations are similar (mean ' + meanRot.toFixed(0) +
          '°) — shift-drag to ROTATE the box between banks');
      /* translation diversity */
      var cs = state.banks.map(function (b) { return [b.pose.x, b.pose.y, b.pose.z]; });
      var mean = [0, 0, 0];
      cs.forEach(function (c) { mean = MV.add(mean, c); });
      mean = MV.scale(mean, 1 / n);
      var spread = 0;
      cs.forEach(function (c) { spread += MV.dot(MV.sub(c, mean), MV.sub(c, mean)); });
      spread = Math.sqrt(spread / n);
      if (spread < 0.22)
        lines.push('  · banked positions cluster (spread ' + (spread * 100).toFixed(0) +
          ' cm) — drag the box around the volume, incl. alt-drag for depth');
    }
    /* apparent size on the displays */
    if (obs.length) {
      var lo = Infinity, hi = -Infinity;
      obs.forEach(function (c) { lo = Math.min(lo, c.p1.u); hi = Math.max(hi, c.p1.u); });
      if ((hi - lo) < 0.30 * DISP_W)
        lines.push('  · the box looks small in camera 1 (' + Math.round(hi - lo) +
          ' of ' + DISP_W + ' px) — bring it closer so each pixel error matters less');
    }
    if (n >= 6 && !state.subpixel)
      lines.push('  · at this point you are near the quantization floor — the');
    if (n >= 6 && !state.subpixel)
      lines.push('    subpixel-corners toggle shows the noise-free limit');
    if (lines.length === 4) lines.push('  · keep banking varied poses — you are doing it right');
    return lines;
  }

  /* ---- main render ------------------------------------------------------ */

  function render() {
    var atoms = MV.poseMolecule(box, state.pose);
    var obs = observe();
    var reproj1 = null, reproj2 = null;
    if (state.solved && state.phase === 'measure') {
      var X = measureCorners(obs);
      reproj1 = []; reproj2 = [];
      Object.keys(X).forEach(function (i) {
        var p1 = MV.project(state.solved.camS1, X[i]);
        var p2 = MV.project(state.solved.camS2, X[i]);
        if (p1) reproj1.push(p1);
        if (p2) reproj2.push(p2);
      });
      /* accumulate the model only when the pose actually changed */
      var key = JSON.stringify(state.pose);
      if (key !== state.lastPoseKey) {
        state.lastPoseKey = key;
        if (!state.tsdf)
          state.tsdf = MV.makeTSDF({ center: [0, 0, 0], size: 2.4, n: 32, tau: 0.09 });
        var C1w = toWorld([0, 0, 0]);
        var C2w = toWorld(MV.scale(MV.matVec(MV.transpose(state.solved.R),
          MV.scale(state.solved.t, state.solved.s)), -1));
        Object.keys(X).forEach(function (i) {
          var Xw = toWorld(X[i]);
          state.cloud.push(Xw);
          state.tsdf.integrate(Xw, C1w);
          state.tsdf.integrate(Xw, C2w);
        });
        /* per-face scan coverage: a face is being captured when it
         * fronts a camera; frontal views count more than grazing ones */
        FACES.forEach(function (f, fi) {
          var g = faceGeom(f, atoms);
          [cam1.C, cam2.C].forEach(function (C) {
            var cos = MV.dot(g.normal, MV.unit(MV.sub(C, g.center)));
            if (cos > 0.25) state.faceCov[fi] += cos;
          });
        });
        while (state.cloud.length > 4000) state.cloud.shift();
      }
    }
    drawOverview(atoms);
    drawView('view1', cam1, atoms, obs, 'p1', reproj1);
    drawView('view2', cam2, atoms, obs, 'p2', reproj2);
    drawCloud(atoms);
    drawSlice();
    drawGraph();
    drawFaceCov();
    updateReadout(obs);
    syncButtons();
  }

  /* ---- controls --------------------------------------------------------- */

  var clamp = function (x, lo, hi) { return Math.max(lo, Math.min(hi, x)); };

  document.getElementById('bank').addEventListener('click', function () {
    bankPose(false);
    render();
  });
  document.getElementById('autobank').addEventListener('change', function (e) {
    state.autoBank = e.target.checked;
  });
  document.getElementById('tomeasure').addEventListener('click', function () {
    if (state.solved) { state.phase = 'measure'; render(); }
  });
  document.getElementById('tosolve').addEventListener('click', function () {
    state.phase = 'solve'; render();
  });
  document.getElementById('resetsolve').addEventListener('click', function () {
    state.banks = []; state.solved = null; state.baselineHist = [];
    state.cloud = []; state.tsdf = null; state.phase = 'solve';
    state.faceCov = [0, 0, 0, 0, 0, 0];
    state.prevS = null; state.stableCount = 0; state.solveStatus = 'none';
    state.bankMsg = 'drag the box, release to auto-bank';
    render();
  });
  /* click the cloud panel: wipe the model and re-estimate immediately
   * from the box's current pose */
  document.getElementById('cloudview').addEventListener('click', function () {
    if (state.phase !== 'measure') return;
    state.cloud = []; state.tsdf = null;
    state.faceCov = [0, 0, 0, 0, 0, 0];
    state.lastPoseKey = '';                   /* forces fresh accumulation */
    render();
  });
  document.getElementById('clearmodel').addEventListener('click', function () {
    state.cloud = []; state.tsdf = null; state.lastPoseKey = '';
    state.faceCov = [0, 0, 0, 0, 0, 0];
    render();
  });
  document.getElementById('subpixel').addEventListener('change', function (e) {
    state.subpixel = e.target.checked;
    render();
  });
  document.getElementById('slicek').addEventListener('input', function (e) {
    state.sliceK = parseInt(e.target.value, 10);
    drawSlice();
  });

  function syncButtons() {
    document.getElementById('bank').disabled = observe().length < 8;
    document.getElementById('tomeasure').disabled = !state.solved || state.phase === 'measure';
    document.getElementById('tosolve').disabled = state.phase === 'solve';
    document.getElementById('bankcount').textContent =
      state.banks.length + ' pose' + (state.banks.length === 1 ? '' : 's') + ' banked';
    var badge = document.getElementById('phasebadge');
    var toMeasureBtn = document.getElementById('tomeasure');
    if (state.phase === 'measure') {
      badge.textContent = 'phase 2: measure';
      badge.className = 'badge badge-solved';
      toMeasureBtn.classList.remove('btn-next');
      document.getElementById('bankmsg').textContent =
        'rotate the box so every colored face fills its coverage bar';
    } else if (state.solveStatus === 'solved') {
      badge.textContent = '✓ RIG SOLVED (stable to <1%)';
      badge.className = 'badge badge-solved';
      toMeasureBtn.classList.add('btn-next');
      document.getElementById('bankmsg').textContent =
        'next: accept solve → measure (more banks still refine it)';
    } else {
      badge.textContent = state.solved
        ? 'phase 1: solving — not yet stable' : 'phase 1: solve';
      badge.className = 'badge' + (state.solved ? ' badge-rough' : '');
      toMeasureBtn.classList.remove('btn-next');
      document.getElementById('bankmsg').textContent = state.bankMsg;
    }
  }

  /* drag the box in the world view: translate / shift-rotate / alt-depth */
  (function () {
    var cv = document.getElementById('overview'), drag = null;
    cv.addEventListener('mousedown', function (e) {
      drag = { x: e.clientX, y: e.clientY };
      e.preventDefault();
    });
    function endDrag() {
      if (drag && state.autoBank && state.phase === 'solve') {
        bankPose(true);
        render();
      }
      drag = null;
    }
    window.addEventListener('mousemove', function (e) {
      if (!drag) return;
      /* releases outside the window never deliver mouseup: treat a
       * button-up move as the release, or the box sticks to the cursor */
      if (e.buttons === 0) { endDrag(); return; }
      var dx = e.clientX - drag.x, dy = e.clientY - drag.y;
      var p = state.pose;
      if (e.shiftKey) {
        p.yaw += dx * 0.01;
        p.pitch = clamp(p.pitch + dy * 0.01, -Math.PI / 2, Math.PI / 2);
      } else if (e.altKey) {
        p.y = clamp(p.y + dx / OSCALE / -0.45, -0.7, 0.7);
      } else {
        p.x = clamp(p.x + dx / OSCALE, -0.7, 0.7);
        p.z = clamp(p.z - dy / OSCALE / 0.9, -0.7, 0.7);
      }
      drag = { x: e.clientX, y: e.clientY };
      render();
    });
    window.addEventListener('mouseup', endDrag);
  })();

  render();
})();
