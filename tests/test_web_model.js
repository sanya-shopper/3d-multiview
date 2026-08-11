/* test_web_model.js -- headless unit test for web/model.js (node).
 * The web page's math layer must agree with the document's equations:
 * these checks mirror the worked examples of doc/multiview.pdf. */
'use strict';

var MV = require('../web/model.js');
var failures = 0;

function ok(cond, name) {
  console.log((cond ? 'ok:   ' : 'FAIL: ') + name);
  if (!cond) failures++;
}
function approx(a, b, tol) { return Math.abs(a - b) <= tol; }

/* --- rotation sanity: R rows orthonormal, R^-1 = R' (doc 3.2) ---------- */
var cam1 = MV.makeCamera({ az: -0.4, el: 0.2, dist: 4, f: 52, W: 64, H: 48 });
var cam2 = MV.makeCamera({ az: 0.5, el: 0.1, dist: 4.5, f: 52, W: 64, H: 48 });
(function () {
  var RtR = MV.matMul(MV.transpose(cam1.R), cam1.R), i, j, worst = 0;
  for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
    worst = Math.max(worst, Math.abs(RtR[i][j] - (i === j ? 1 : 0)));
  ok(worst < 1e-12, 'camera rotation is orthonormal');
})();

/* --- camera center: projecting C is degenerate, C + z axis maps ahead -- */
(function () {
  var Xc = MV.add(MV.matVec(cam1.R, cam1.C), cam1.t);
  ok(MV.norm(Xc) < 1e-12, 'camera center maps to camera origin (doc ex. center)');
})();

/* --- skew matrix: [v]x w = v x w (doc eq. skew) ------------------------ */
(function () {
  var v = [1, 2, 3], w = [-2, 0.5, 4];
  var d = MV.sub(MV.matVec(MV.skew(v), w), MV.cross(v, w));
  ok(MV.norm(d) < 1e-12, '[v]x w equals v x w');
})();

/* --- epipolar constraint holds exactly for continuous projections ------ */
var F = MV.fundamental(cam1, cam2);
(function () {
  var X = [0.3, -0.2, 0.4], p1 = MV.project(cam1, X), p2 = MV.project(cam2, X);
  var r = MV.epipolarResidual(F, p1, p2);
  ok(approx(r, 0, 1e-9), 'x2\' F x1 = 0 for continuous pixels (doc eq. epiconstraint)');
})();

/* --- F has rank 2: F e1 = 0 at the epipole ----------------------------- */
(function () {
  var e1 = MV.project(cam1, cam2.C);       /* image of the other center */
  var l = MV.matVec(F, [e1.u, e1.v, 1]);
  ok(MV.norm(l) < 1e-6 * MV.norm([F[0][0], F[1][1], F[2][2]]) + 1e-9,
     'epipole is the null vector of F (doc 3.3/6)');
})();

/* --- DLT triangulation recovers the point from continuous pixels ------- */
(function () {
  var X = [-0.25, 0.4, -0.3];
  var Xhat = MV.triangulate(cam1, MV.project(cam1, X),
                            cam2, MV.project(cam2, X));
  ok(Xhat && MV.norm(MV.sub(Xhat, X)) < 1e-9,
     'DLT recovers a point exactly from continuous pixels (doc eq. dlt)');
})();

/* --- quantized pixels: error small and reprojection consistent --------- */
(function () {
  var X = [0.1, 0.2, 0.05];
  var q1 = MV.quantize(MV.project(cam1, X));
  var q2 = MV.quantize(MV.project(cam2, X));
  var Xhat = MV.triangulate(cam1, q1, cam2, q2);
  var err = MV.norm(MV.sub(Xhat, X));
  ok(err > 0 && err < 0.2,
     'quantized triangulation error is nonzero but bounded (' +
     (err * 1000).toFixed(1) + ' mm)');
  var rms = MV.reprojError([{ cam: cam1, p: q1 }, { cam: cam2, p: q2 }], Xhat);
  ok(rms < Math.sqrt(0.5) + 1e-6,
     'RMS reprojection of the DLT point is sub-pixel');
})();

/* --- depth error law: sigma_Z = Z^2 sigma_d / (f B) (doc eq. 14) ------- */
(function () {
  var s = MV.depthSigma(4.7, 800, 0.5, 0.42);
  ok(approx(s, 0.023, 0.0005),
     'depth error law reproduces doc ex. deptherr (23 mm)');
})();

/* --- back-projection round trip ---------------------------------------- */
(function () {
  var X = [0.4, -0.15, 0.22];
  var p = MV.project(cam1, X);
  var Xb = MV.backproject(cam1, p, p.z);
  ok(MV.norm(MV.sub(Xb, X)) < 1e-12,
     'backproject inverts project at the observed depth (doc section 4)');
})();

/* --- plane-induced homography (doc eq. planehom) ----------------------- */
(function () {
  var d = 3.8;                       /* plane depth in camera-1 frame */
  var H = MV.planeHomography(cam1, cam2, d);
  /* a point ON the plane transfers exactly */
  var Xon = MV.backproject(cam1, { u: 21.3, v: 30.1 }, d);
  var ph = MV.applyH(H, MV.project(cam1, Xon));
  var p2 = MV.project(cam2, Xon);
  ok(Math.hypot(ph.u - p2.u, ph.v - p2.v) < 1e-9,
     'plane homography transfers on-plane points exactly');
  /* a point OFF the plane misses by parallax */
  var Xoff = MV.backproject(cam1, { u: 21.3, v: 30.1 }, d - 0.8);
  var ph2 = MV.applyH(H, MV.project(cam1, Xoff));
  var p22 = MV.project(cam2, Xoff);
  ok(Math.hypot(ph2.u - p22.u, ph2.v - p22.v) > 0.5,
     'off-plane points miss by parallax (what depth is made of)');
})();

/* --- explicit 6-DOF camera pose ---------------------------------------- */
(function () {
  var pos = [3.1, -1.4, 0.8], target = [0.2, 0.3, -0.1];
  var a = MV.aimAngles(pos, target);
  var cam = MV.makeCameraPose({ pos: pos, pan: a.pan, tilt: a.tilt, roll: 0,
                                f: 52, W: 64, H: 48 });
  var p = MV.project(cam, target);
  ok(p && approx(p.u, 32, 1e-9) && approx(p.v, 24, 1e-9),
     'aimAngles + makeCameraPose put the target on the principal point');

  /* orbit wrapper and pose form agree */
  var orb = MV.makeCamera({ az: 0.7, el: -0.25, dist: 3.7, f: 52, W: 64, H: 48 });
  var a2 = MV.aimAngles(orb.C, [0, 0, 0]);
  var pose = MV.makeCameraPose({ pos: orb.C, pan: a2.pan, tilt: a2.tilt,
                                 roll: 0, f: 52, W: 64, H: 48 });
  var i, j, worst = 0;
  for (i = 0; i < 3; i++) for (j = 0; j < 4; j++)
    worst = Math.max(worst, Math.abs(orb.P[i][j] - pose.P[i][j]));
  ok(worst < 1e-9, 'orbit camera equals pose camera aimed the same way');

  /* a 90-degree roll turns image u-offsets into v-offsets */
  var rolled = MV.makeCameraPose({ pos: pos, pan: a.pan, tilt: a.tilt,
                                   roll: Math.PI / 2, f: 52, W: 64, H: 48 });
  var off = [target[0], target[1], target[2] + 0.3];   /* world offset */
  var p0 = MV.project(cam, off), p90 = MV.project(rolled, off);
  /* right-hand-rule roll about the optical axis: u' = v, v' = -u */
  ok(approx(p90.u - 32, (p0.v - 24), 1e-6) &&
     approx(p90.v - 24, -(p0.u - 32), 1e-6),
     'roll of 90 degrees rotates the image plane as expected');
})();

/* --- svd3 reconstructs and is orthonormal ------------------------------ */
(function () {
  var M = [[2, -1, 0.5], [0.3, 1.7, -2.2], [1.1, 0.4, 0.9]];
  var s = MV.svd3(M), i, j, k, worst = 0;
  for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
    var m = 0;
    for (k = 0; k < 3; k++) m += s.S[k] * s.U[k][i] * s.V[k][j];
    worst = Math.max(worst, Math.abs(m - M[i][j]));
  }
  ok(worst < 1e-9, 'svd3: U S V\' reconstructs the matrix');
  ok(s.S[0] >= s.S[1] && s.S[1] >= s.S[2], 'svd3: singular values descending');
})();

/* --- rig solving: 8-point + E decomposition + chirality ---------------- */
(function () {
  /* synthesize exact correspondences of a box seen by the two cameras */
  var box = MV.makeBox(0.8, 0.5, 0.3);
  var pairs = [];
  [{ yaw: 0.4, pitch: 0.2, roll: 0.1, x: 0, y: 0, z: 0 },
   { yaw: -0.7, pitch: -0.3, roll: 0.9, x: 0.3, y: -0.2, z: 0.2 },
   { yaw: 1.4, pitch: 0.5, roll: -0.5, x: -0.3, y: 0.3, z: -0.2 }]
    .forEach(function (pose) {
      MV.poseMolecule(box, pose).forEach(function (X) {
        var p1 = MV.project(cam1, X), p2 = MV.project(cam2, X);
        if (p1 && p2) pairs.push({ p1: p1, p2: p2 });
      });
    });
  var F = MV.eightPoint(pairs);
  ok(F && MV.epipolarRMS(F, pairs) < 1e-6,
     '8-point: epipolar RMS ~ 0 on exact correspondences');
  var Fe = MV.essentialFromPairs(pairs, cam1.K, cam2.K);
  ok(Fe && MV.epipolarRMS(Fe, pairs) < 1e-6,
     'calibrated E solve: epipolar RMS ~ 0 on exact correspondences');
  F = Fe;                                    /* downstream uses the E path */
  var rel = MV.relativePose(F, cam1.K, cam2.K, pairs, 64, 48);
  var Rtrue = MV.matMul(cam2.R, MV.transpose(cam1.R));
  var ttrue = MV.sub(cam2.t, MV.matVec(Rtrue, cam1.t));
  ok(MV.rotationAngle(rel.R, Rtrue) < 1e-6,
     'E decomposition recovers the relative rotation');
  ok(Math.abs(MV.dot(rel.t, MV.unit(ttrue))) > 1 - 1e-9,
     'E decomposition recovers the baseline direction (unit scale)');
  /* scale anchoring: triangulate in the unit-baseline rig, measure one
   * known edge, and the recovered scale is the true baseline length */
  var camS1 = MV.camFromRt(cam1.K, [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0], 64, 48);
  var camS2 = MV.camFromRt(cam2.K, rel.R, rel.t, 64, 48);
  var X0 = MV.triangulate(camS1, pairs[0].p1, camS2, pairs[0].p2);
  var X1 = MV.triangulate(camS1, pairs[1].p1, camS2, pairs[1].p2);
  var estEdge = MV.norm(MV.sub(X1, X0));      /* corners 0-1: the a edge */
  var scale = 0.8 / estEdge;
  ok(Math.abs(scale - MV.norm(ttrue)) < 1e-6,
     'known edge anchors the metric scale to the true baseline');
})();

/* --- TSDF scatter integration ------------------------------------------ */
(function () {
  var t = MV.makeTSDF({ center: [0, 0, 0], size: 2, n: 20, tau: 0.15 });
  var X = [0.4, 0.1, 0], C = [-0.9, 0.1, 0];
  t.integrate(X, C);
  var k = Math.floor((0 - t.min[2]) / t.cell);
  var s = t.slice(k);
  function at(x, y) {
    return s[Math.floor((y - t.min[1]) / t.cell)][Math.floor((x - t.min[0]) / t.cell)];
  }
  ok(Math.abs(at(0.4, 0.1)) < t.tau, 'TSDF: band voxel holds a small signed distance');
  ok(at(-0.4, 0.1) === t.tau, 'TSDF: free-space voxel carved at +tau');
  ok(isNaN(at(0.8, 0.1)), 'TSDF: voxel behind the surface untouched');
  ok(t.observed() > 0, 'TSDF: observed-voxel count positive');
})();

/* --- molecule/pose plumbing -------------------------------------------- */
(function () {
  var mol = MV.makeMolecule();
  var atoms = MV.poseMolecule(mol,
    { yaw: 0.3, pitch: -0.2, roll: 0.7, x: 0.1, y: 0, z: -0.1 });
  ok(atoms.length === mol.atoms.length, 'pose preserves atom count');
  var d0 = MV.norm(MV.sub(mol.atoms[1].p, mol.atoms[0].p));
  var d1 = MV.norm(MV.sub(atoms[1], atoms[0]));
  ok(approx(d0, d1, 1e-12), 'pose is rigid (bond lengths preserved)');
})();

if (failures) {
  console.log(failures + ' web-model test(s) FAILED');
  process.exit(1);
}
console.log('all web-model tests passed');
