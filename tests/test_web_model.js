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

/* --- molecule/pose plumbing -------------------------------------------- */
(function () {
  var mol = MV.makeMolecule();
  var atoms = MV.poseMolecule(mol, { yaw: 0.3, pitch: -0.2, x: 0.1, y: 0, z: -0.1 });
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
