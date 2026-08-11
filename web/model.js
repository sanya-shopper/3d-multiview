/* model.js -- the data model and the mathematics, DOM-free.
 *
 * Everything here mirrors the learning document doc/multiview.pdf:
 *   - projection x ~ K [R|t] X            (doc eq. "pinhole", section 4)
 *   - fundamental matrix F = K2^-T [t]x R K1^-1   (section 6)
 *   - epipolar constraint x2' F x1 = 0            (section 6)
 *   - DLT triangulation from N=2 views           (section 7)
 *   - depth uncertainty sigma_Z = Z^2/(f B) sigma_d (section 8)
 *
 * No DOM access: this file is loaded by the page (app.js drives it) and
 * by the headless unit test tests/test_web_model.js.
 */
'use strict';

var MV = (function () {

  /* ---------------- small linear algebra (row-major arrays) ------------- */

  function v3(x, y, z) { return [x, y, z]; }
  function add(a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
  function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
  function scale(a, s) { return [a[0] * s, a[1] * s, a[2] * s]; }
  function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function cross(a, b) {
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]];
  }
  function norm(a) { return Math.sqrt(dot(a, a)); }
  function unit(a) { return scale(a, 1 / norm(a)); }

  /* 3x3 matrices as arrays of 3 rows */
  function matVec(M, v) {
    return [dot(M[0], v), dot(M[1], v), dot(M[2], v)];
  }
  function matMul(A, B) {
    var C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], i, j, k;
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        for (k = 0; k < 3; k++) C[i][j] += A[i][k] * B[k][j];
    return C;
  }
  function transpose(M) {
    return [[M[0][0], M[1][0], M[2][0]],
            [M[0][1], M[1][1], M[2][1]],
            [M[0][2], M[1][2], M[2][2]]];
  }
  function inv3(M) {
    var a = M[0][0], b = M[0][1], c = M[0][2],
        d = M[1][0], e = M[1][1], f = M[1][2],
        g = M[2][0], h = M[2][1], i = M[2][2];
    var A = e * i - f * h, B = c * h - b * i, C = b * f - c * e;
    var det = a * A + d * B + g * C;
    return [[A / det, B / det, C / det],
            [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
            [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det]];
  }
  /* cross-product matrix [v]x  (doc eq. "skew", section 3.3) */
  function skew(v) {
    return [[0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0]];
  }

  /* full eigendecomposition of a symmetric NxN matrix by cyclic Jacobi:
   * returns {d: eigenvalues, V: matrix whose COLUMN j is eigenvector j} */
  function eigSym(S) {
    var n = S.length, V = [], A = [], i, j, p, q, sweep;
    for (i = 0; i < n; i++) {
      V.push([]); A.push([]);
      for (j = 0; j < n; j++) { V[i].push(i === j ? 1 : 0); A[i].push(S[i][j]); }
    }
    for (sweep = 0; sweep < 30; sweep++) {
      var off = 0;
      for (p = 0; p < n; p++) for (q = p + 1; q < n; q++) off += A[p][q] * A[p][q];
      if (off < 1e-24) break;
      for (p = 0; p < n; p++) for (q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-30) continue;
        var theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        var t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        var cth = 1 / Math.sqrt(t * t + 1), sth = t * cth, k;
        for (k = 0; k < n; k++) {
          var akp = A[k][p], akq = A[k][q];
          A[k][p] = cth * akp - sth * akq;
          A[k][q] = sth * akp + cth * akq;
        }
        for (k = 0; k < n; k++) {
          var apk = A[p][k], aqk = A[q][k];
          A[p][k] = cth * apk - sth * aqk;
          A[q][k] = sth * apk + cth * aqk;
        }
        for (k = 0; k < n; k++) {
          var vkp = V[k][p], vkq = V[k][q];
          V[k][p] = cth * vkp - sth * vkq;
          V[k][q] = sth * vkp + cth * vkq;
        }
      }
    }
    var d = [];
    for (i = 0; i < n; i++) d.push(A[i][i]);
    return { d: d, V: V };
  }

  /* smallest eigenvector of a symmetric NxN matrix -- the homogeneous
   * least-squares solver (doc section 3.4) */
  function smallestEigvec(S) {
    var e = eigSym(S), n = S.length, i;
    var best = 0;
    for (i = 1; i < n; i++) if (e.d[i] < e.d[best]) best = i;
    var v = [];
    for (i = 0; i < n; i++) v.push(e.V[i][best]);
    return v;
  }

  /* SVD of a 3x3 matrix via the eigendecomposition of M'M:
   * returns {U, S, V} with U, V arrays of COLUMN vectors, S descending. */
  function svd3(M) {
    var e = eigSym(matMul(transpose(M), M));
    var order = [0, 1, 2].sort(function (a, b) { return e.d[b] - e.d[a]; });
    var S = [], V = [], U = [null, null, null], k;
    for (k = 0; k < 3; k++) {
      var oi = order[k];
      V.push([e.V[0][oi], e.V[1][oi], e.V[2][oi]]);
      S.push(Math.sqrt(Math.max(0, e.d[oi])));
    }
    /* left vectors: u = M v / s where s is significant RELATIVE to s1;
     * a tiny s makes M v pure noise, so complete the basis instead */
    var floor = 1e-6 * (S[0] || 1);
    for (k = 0; k < 3; k++)
      if (S[k] > floor) U[k] = unit(matVec(M, V[k]));
    if (!U[0]) return { U: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], S: S, V: V };
    if (!U[1]) {
      var any = Math.abs(U[0][0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
      U[1] = unit(cross(U[0], any));
    }
    if (!U[2]) U[2] = unit(cross(U[0], U[1]));
    return { U: U, S: S, V: V };
  }

  function colsToMat(c) {
    return [[c[0][0], c[1][0], c[2][0]],
            [c[0][1], c[1][1], c[2][1]],
            [c[0][2], c[1][2], c[2][2]]];
  }
  function det3(M) {
    return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
         - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
         + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  }

  /* ---------------- the molecule: a stylized 3D graph ------------------- */
  /* A cyclohexane-like chair ring with two substituents -- chosen because
   * the chair is genuinely three-dimensional, so moving the cameras
   * visibly changes depth ordering.  Geometry is stylized, not chemistry;
   * units are metres because the point is the rig, not the molecule. */
  function makeMolecule() {
    var atoms = [], bonds = [], i;
    var R = 0.55, dz = 0.22;
    for (i = 0; i < 6; i++) {
      var a = Math.PI / 3 * i;
      atoms.push({ el: 'C', color: '#555b63',
                   p: v3(R * Math.cos(a), R * Math.sin(a), (i % 2 ? dz : -dz)) });
      bonds.push([i, (i + 1) % 6]);
    }
    atoms.push({ el: 'O', color: '#c0392b', p: v3(0.95, 0.0, 0.45) });  /* 6 */
    bonds.push([0, 6]);
    atoms.push({ el: 'N', color: '#2d6cdf', p: v3(-0.55, -0.75, -0.5) }); /* 7 */
    bonds.push([4, 7]);
    /* more colored substituents: distinct landmarks to watch while the
     * cameras and the molecule move */
    atoms.push({ el: 'S', color: '#b8860b', p: v3(-0.50, 1.00, 0.15) });  /* 8 */
    bonds.push([2, 8]);
    atoms.push({ el: 'F', color: '#27ae60', p: v3(0.55, 0.95, 0.62) });   /* 9 */
    bonds.push([1, 9]);
    atoms.push({ el: 'P', color: '#8e44ad', p: v3(0.50, -0.95, 0.62) });  /* 10 */
    bonds.push([5, 10]);
    atoms.push({ el: 'Mg', color: '#d81b60', p: v3(-1.05, 0.15, 0.62) }); /* 11 */
    bonds.push([3, 11]);
    return { atoms: atoms, bonds: bonds };
  }

  /* pose the molecule: yaw/pitch/roll (radians) then translate */
  function poseMolecule(mol, pose) {
    var cy = Math.cos(pose.yaw), sy = Math.sin(pose.yaw);
    var cp = Math.cos(pose.pitch), sp = Math.sin(pose.pitch);
    var cr = Math.cos(pose.roll || 0), sr = Math.sin(pose.roll || 0);
    var Ry = [[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]];
    var Rx = [[1, 0, 0], [0, cp, -sp], [0, sp, cp]];
    var Rz = [[cr, -sr, 0], [sr, cr, 0], [0, 0, 1]];
    var R = matMul(Rz, matMul(Rx, Ry));
    return mol.atoms.map(function (a) {
      return add(matVec(R, a.p), [pose.x, pose.y, pose.z]);
    });
  }

  /* ---------------- cameras --------------------------------------------- */

  /* Build a camera from an explicit pose: center C, optical axis zc,
   * roll about that axis. Returns {K, R, t, C, P(3x4), f, W, H} --
   * doc section 4: rows of R are the camera axes in world coordinates. */
  function buildCamera(C, zc, roll, par) {
    var up = v3(0, 0, 1);
    var xc = cross(zc, up);
    if (norm(xc) < 1e-9) xc = v3(1, 0, 0);    /* looking straight down z */
    xc = unit(xc);
    var yc = cross(zc, xc);
    if (roll) {                                /* rotate image axes about zc */
      var cr = Math.cos(roll), sr = Math.sin(roll);
      var xr = add(scale(xc, cr), scale(yc, sr));
      yc = cross(zc, xr);
      xc = xr;
    }
    var R = [xc, yc, zc];                     /* rows = camera axes (doc 3.2) */
    var t = scale(matVec(R, C), -1);          /* t = -R C */
    var K = [[par.f, 0, par.W / 2],
             [0, par.f, par.H / 2],
             [0, 0, 1]];
    var P = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], i, j;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        P[i][j] = K[i][0] * R[0][j] + K[i][1] * R[1][j] + K[i][2] * R[2][j];
      P[i][3] = K[i][0] * t[0] + K[i][1] * t[1] + K[i][2] * t[2];
    }
    return { K: K, R: R, t: t, C: C, P: P, f: par.f, W: par.W, H: par.H };
  }

  /* Explicit 6-DOF pose: position pos = [x,y,z]; pan (heading of the
   * optical axis in the ground plane), tilt (its elevation), roll about
   * the axis. pan/tilt/roll in radians. */
  function makeCameraPose(par) {
    var cp = Math.cos(par.pan), sp = Math.sin(par.pan);
    var ct = Math.cos(par.tilt), st = Math.sin(par.tilt);
    var zc = v3(ct * cp, ct * sp, st);
    return buildCamera(par.pos.slice(), zc, par.roll || 0, par);
  }

  /* pan/tilt that point a camera at C toward target (roll unchanged) */
  function aimAngles(C, target) {
    var d = unit(sub(target, C));
    return { pan: Math.atan2(d[1], d[0]), tilt: Math.asin(d[2]) };
  }

  /* Orbit parameterization around a target: azimuth az, elevation el,
   * distance dist; the camera looks at the target with zero roll.
   * Kept as a thin wrapper over the pose form. */
  function makeCamera(par) {
    var target = par.target || v3(0, 0, 0);
    var caz = Math.cos(par.az), saz = Math.sin(par.az);
    var cel = Math.cos(par.el), sel = Math.sin(par.el);
    var C = add(target, scale(v3(caz * cel, saz * cel, sel), par.dist));
    return buildCamera(C, unit(sub(target, C)), 0, par);
  }

  /* project world point -> {u,v,z} (continuous pixels; z = camera depth) */
  function project(cam, X) {
    var Xc = add(matVec(cam.R, X), cam.t);
    if (Xc[2] <= 1e-6) return null;           /* behind the camera */
    var xn = Xc[0] / Xc[2], yn = Xc[1] / Xc[2];
    return { u: cam.K[0][0] * xn + cam.K[0][2],
             v: cam.K[1][1] * yn + cam.K[1][2],
             z: Xc[2] };
  }

  /* quantize to the pixel the sample lands in; center of that pixel */
  function quantize(p) {
    return { u: Math.floor(p.u) + 0.5, v: Math.floor(p.v) + 0.5 };
  }

  /* back-project pixel {u,v} at camera depth z to the world point on that
   * ray (doc section 4: the ray C + s d, here with the depth known) */
  function backproject(cam, p, z) {
    var Ki = inv3(cam.K);
    var dc = matVec(Ki, [p.u, p.v, 1]);       /* camera-frame dir, dc[2]=1 */
    var Xc = scale(dc, z);
    return matVec(transpose(cam.R), sub(Xc, cam.t));
  }

  /* ---------------- two-view relations ---------------------------------- */

  /* F = K2^-T [t21]x R21 K1^-1 with R21 = R2 R1', t21 = t2 - R21 t1
   * (essential matrix E = [t]x R, doc section 6) */
  function fundamental(cam1, cam2) {
    var R21 = matMul(cam2.R, transpose(cam1.R));
    var t21 = sub(cam2.t, matVec(R21, cam1.t));
    var E = matMul(skew(t21), R21);
    return matMul(transpose(inv3(cam2.K)), matMul(E, inv3(cam1.K)));
  }

  /* epipolar residual x2' F x1 for pixel points {u,v} */
  function epipolarResidual(F, p1, p2) {
    var l = matVec(F, [p1.u, p1.v, 1]);
    return (l[0] * p2.u + l[1] * p2.v + l[2]);
  }

  /* the epipolar line l2 = F x1 as [a,b,c], a u + b v + c = 0 */
  function epipolarLine(F, p1) { return matVec(F, [p1.u, p1.v, 1]); }

  /* Plane-induced homography from camera 1 pixels to camera 2 pixels
   * for the fronto-parallel (to camera 1) plane at camera-1 depth d:
   * H = K2 (R21 - t21 n' / d) K1^-1 with n = camera-1 optical axis
   * (doc eq. planehom, section 5.1 -- there induced by the target plane,
   * here by an assumed scene plane). Exact only for points ON the plane;
   * everything off it misses by parallax, which is the information
   * stereo depth is made of. */
  function planeHomography(cam1, cam2, d) {
    var R21 = matMul(cam2.R, transpose(cam1.R));
    var t21 = sub(cam2.t, matVec(R21, cam1.t));
    /* on the plane n.Xc1 = d (n = (0,0,1)):
     * Xc2 = R21 Xc1 + t21 (n.Xc1 / d) = (R21 + t21 n'/d) Xc1,
     * so the n'/d term adds t21/d to column 3 of R21 */
    var M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], i, j;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) M[i][j] = R21[i][j];
      M[i][2] += t21[i] / d;
    }
    return matMul(cam2.K, matMul(M, inv3(cam1.K)));
  }

  /* apply a homography to a pixel point */
  function applyH(H, p) {
    var x = matVec(H, [p.u, p.v, 1]);
    if (Math.abs(x[2]) < 1e-12) return null;
    return { u: x[0] / x[2], v: x[1] / x[2] };
  }

  /* DLT triangulation from two views (doc eq. "dlt", section 7):
   * rows u P3 - P1, v P3 - P2 per view; smallest right singular vector. */
  function triangulate(cam1, p1, cam2, p2) {
    var rows = [], A;
    [[cam1, p1], [cam2, p2]].forEach(function (vp) {
      var P = vp[0].P, u = vp[1].u, v = vp[1].v, j;
      var r1 = [], r2 = [];
      for (j = 0; j < 4; j++) {
        r1.push(u * P[2][j] - P[0][j]);
        r2.push(v * P[2][j] - P[1][j]);
      }
      rows.push(r1); rows.push(r2);
    });
    /* normal matrix A'A (4x4 symmetric) */
    A = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
    rows.forEach(function (r) {
      var i, j;
      for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) A[i][j] += r[i] * r[j];
    });
    var h = smallestEigvec(A);
    if (Math.abs(h[3]) < 1e-12) return null;  /* point at infinity */
    return v3(h[0] / h[3], h[1] / h[3], h[2] / h[3]);
  }

  /* RMS reprojection error (px) of world point X against observations */
  function reprojError(views, X) {
    var s = 0, n = 0;
    views.forEach(function (vp) {
      var p = project(vp.cam, X);
      if (!p) return;
      var du = p.u - vp.p.u, dv = p.v - vp.p.v;
      s += du * du + dv * dv; n += 2;
    });
    return n ? Math.sqrt(s / n) : NaN;
  }

  /* predicted depth uncertainty sigma_Z = Z^2/(f B) sigma_d (doc eq. 14) */
  function depthSigma(Z, f, B, sigmaD) { return Z * Z * sigmaD / (f * B); }

  /* ---------------- solving a rig from correspondences ------------------ */

  /* camera object from explicit K, R, t (same shape as makeCamera) */
  function camFromRt(K, R, t, W, H) {
    var C = scale(matVec(transpose(R), t), -1);
    var P = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], i, j;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        P[i][j] = K[i][0] * R[0][j] + K[i][1] * R[1][j] + K[i][2] * R[2][j];
      P[i][3] = K[i][0] * t[0] + K[i][1] * t[1] + K[i][2] * t[2];
    }
    return { K: K, R: R, t: t, C: C, P: P, f: K[0][0], W: W, H: H };
  }

  /* Hartley normalization: translate to centroid, scale mean dist to
   * sqrt(2); returns {T, pts} (doc section 6, the "normalized" in the
   * normalized 8-point algorithm) */
  function hartley(pts) {
    var n = pts.length, cx = 0, cy = 0, i;
    for (i = 0; i < n; i++) { cx += pts[i].u; cy += pts[i].v; }
    cx /= n; cy /= n;
    var md = 0;
    for (i = 0; i < n; i++) md += Math.hypot(pts[i].u - cx, pts[i].v - cy);
    md /= n;
    var s = Math.SQRT2 / (md || 1);
    var out = pts.map(function (p) {
      return { u: s * (p.u - cx), v: s * (p.v - cy) };
    });
    return { T: [[s, 0, -s * cx], [0, s, -s * cy], [0, 0, 1]], pts: out };
  }

  /* Normalized 8-point algorithm (doc section 6): pairs of pixel points
   * [{p1:{u,v}, p2:{u,v}}, ...], n >= 8 -> fundamental matrix with rank-2
   * enforcement, unit Frobenius norm. */
  function eightPoint(pairs) {
    if (pairs.length < 8) return null;
    var n1 = hartley(pairs.map(function (c) { return c.p1; }));
    var n2 = hartley(pairs.map(function (c) { return c.p2; }));
    var A = [], i, j, k;
    for (i = 0; i < 9; i++) { A.push([]); for (j = 0; j < 9; j++) A[i].push(0); }
    for (k = 0; k < pairs.length; k++) {
      var a = n1.pts[k], b = n2.pts[k];
      var r = [b.u * a.u, b.u * a.v, b.u, b.v * a.u, b.v * a.v, b.v, a.u, a.v, 1];
      for (i = 0; i < 9; i++) for (j = 0; j < 9; j++) A[i][j] += r[i] * r[j];
    }
    var f = smallestEigvec(A);
    var Fn = [[f[0], f[1], f[2]], [f[3], f[4], f[5]], [f[6], f[7], f[8]]];
    /* rank-2 enforcement: drop the smallest singular value */
    var sv = svd3(Fn), F2 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (k = 0; k < 2; k++)
      for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
        F2[i][j] += sv.S[k] * sv.U[k][i] * sv.V[k][j];
    var F = matMul(transpose(n2.T), matMul(F2, n1.T));
    var fr = 0;
    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) fr += F[i][j] * F[i][j];
    fr = Math.sqrt(fr);
    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) F[i][j] /= fr;
    return F;
  }

  /* mean point-to-epipolar-line distance in pixels (solve quality) */
  function epipolarRMS(F, pairs) {
    var s = 0, n = 0;
    pairs.forEach(function (c) {
      var l = matVec(F, [c.p1.u, c.p1.v, 1]);
      var d = (l[0] * c.p2.u + l[1] * c.p2.v + l[2]) / Math.hypot(l[0], l[1]);
      s += d * d; n++;
    });
    return n ? Math.sqrt(s / n) : NaN;
  }

  /* Decompose E = K2' F K1 into the relative pose (R, t) of camera 2 in
   * camera 1's frame, unit baseline; the 4-fold ambiguity is resolved by
   * the chirality test -- triangulated points must be in front of both
   * cameras (doc section 6). */
  function relativePose(F, K1, K2, pairs, W, H) {
    var E = matMul(transpose(K2), matMul(F, K1));
    var sv = svd3(E);
    var U = sv.U.slice(), V = sv.V.slice(), i;
    if (det3(colsToMat(U)) < 0) U = U.map(function (c) { return scale(c, -1); });
    if (det3(colsToMat(V)) < 0) V = V.map(function (c) { return scale(c, -1); });
    var Um = colsToMat(U), Vm = colsToMat(V);
    var Wm = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
    var Ra = matMul(Um, matMul(Wm, transpose(Vm)));
    var Rb = matMul(Um, matMul(transpose(Wm), transpose(Vm)));
    var tv = U[2];
    var best = null;
    [[Ra, tv], [Ra, scale(tv, -1)], [Rb, tv], [Rb, scale(tv, -1)]]
      .forEach(function (cand) {
        var R = cand[0], t = cand[1];
        var cam1 = camFromRt(K1, [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0], W, H);
        var cam2 = camFromRt(K2, R, t, W, H);
        var good = 0, tried = 0;
        for (i = 0; i < pairs.length && tried < 24; i += Math.max(1, Math.floor(pairs.length / 24))) {
          tried++;
          var X = triangulate(cam1, pairs[i].p1, cam2, pairs[i].p2);
          if (!X) continue;
          if (X[2] > 0 && (add(matVec(R, X), t))[2] > 0) good++;
        }
        if (!best || good > best.good) best = { R: R, t: t, good: good };
      });
    return best;
  }

  /* Calibrated two-view solve: run the 8-point machinery on K-normalized
   * coordinates and project the result onto the essential manifold
   * (sigma1 = sigma2, sigma3 = 0) -- much better conditioned than
   * estimating F in raw pixels when the correspondences cluster in a
   * small part of the image. Returns F (pixels) for downstream reuse. */
  function essentialFromPairs(pairs, K1, K2) {
    var K1i = inv3(K1), K2i = inv3(K2);
    var norm = pairs.map(function (c) {
      var a = matVec(K1i, [c.p1.u, c.p1.v, 1]);
      var b = matVec(K2i, [c.p2.u, c.p2.v, 1]);
      return { p1: { u: a[0] / a[2], v: a[1] / a[2] },
               p2: { u: b[0] / b[2], v: b[1] / b[2] } };
    });
    var Fn = eightPoint(norm);
    if (!Fn) return null;
    var sv = svd3(Fn);
    var s = (sv.S[0] + sv.S[1]) / 2;
    var E = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], i, j, k;
    for (k = 0; k < 2; k++)
      for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)
        E[i][j] += s * sv.U[k][i] * sv.V[k][j];
    return matMul(transpose(K2i), matMul(E, K1i));
  }

  /* geodesic angle between two rotations, radians */
  function rotationAngle(Ra, Rb) {
    var M = matMul(transpose(Ra), Rb);
    var c = (M[0][0] + M[1][1] + M[2][2] - 1) / 2;
    return Math.acos(Math.max(-1, Math.min(1, c)));
  }

  /* ---------------- a parallelepiped of known dimensions ---------------- */
  /* Same shape as makeMolecule so poseMolecule works on it: corner i has
   * bit0/1/2 of i selecting the +x/+y/+z face; edges connect indices
   * differing in one bit, and an edge's axis is that bit. */
  function makeBox(a, b, c) {
    var atoms = [], bonds = [], i;
    for (i = 0; i < 8; i++) {
      atoms.push({ el: 'corner', color: '#2d6cdf',
                   p: v3((i & 1 ? a : -a) / 2, (i & 2 ? b : -b) / 2,
                         (i & 4 ? c : -c) / 2) });
    }
    for (i = 0; i < 8; i++) [1, 2, 4].forEach(function (bit) {
      if ((i ^ bit) > i) bonds.push([i, i ^ bit]);
    });
    return { atoms: atoms, bonds: bonds, dims: [a, b, c] };
  }

  /* ---------------- a small TSDF (doc section 9.2), scatter form -------- */
  /* Cubic grid over [center - size/2, center + size/2]^3, n^3 voxels.
   * integrate(X, C) marches the ray C -> X (the scatter update of doc
   * section 3.5): voxels along the free segment get the clamped +tau,
   * voxels in the band get the signed offset, nothing behind X + tau is
   * touched. Weighted running average per voxel (doc eq. tsdf). */
  function makeTSDF(opts) {
    var n = opts.n, size = opts.size, tau = opts.tau;
    var min = [opts.center[0] - size / 2, opts.center[1] - size / 2,
               opts.center[2] - size / 2];
    var cell = size / n;
    var D = new Float64Array(n * n * n);
    var Wt = new Float64Array(n * n * n);
    function update(p, d) {
      var i = Math.floor((p[0] - min[0]) / cell);
      var j = Math.floor((p[1] - min[1]) / cell);
      var k = Math.floor((p[2] - min[2]) / cell);
      if (i < 0 || j < 0 || k < 0 || i >= n || j >= n || k >= n) return;
      var x = (k * n + j) * n + i;
      D[x] = (Wt[x] * D[x] + d) / (Wt[x] + 1);
      Wt[x] += 1;
    }
    function integrate(X, C) {
      var dir = sub(X, C), L = norm(dir);
      if (L < 1e-9) return;
      dir = scale(dir, 1 / L);
      var s;
      for (s = 0; s <= L + tau; s += cell / 2) {
        var eta = L - s;                      /* + in front of X, - behind */
        if (eta < -tau) break;
        update(add(C, scale(dir, s)), Math.min(tau, eta));
      }
    }
    function slice(k) {
      var out = [], i, j;
      for (j = 0; j < n; j++) {
        var row = [];
        for (i = 0; i < n; i++) {
          var x = (k * n + j) * n + i;
          row.push(Wt[x] > 0 ? D[x] : NaN);
        }
        out.push(row);
      }
      return out;
    }
    function observed() {
      var c = 0, i;
      for (i = 0; i < Wt.length; i++) if (Wt[i] > 0) c++;
      return c;
    }
    return { integrate: integrate, slice: slice, observed: observed,
             n: n, cell: cell, min: min, tau: tau };
  }

  return {
    v3: v3, add: add, sub: sub, scale: scale, dot: dot, cross: cross,
    norm: norm, unit: unit, matVec: matVec, matMul: matMul,
    transpose: transpose, inv3: inv3, skew: skew,
    smallestEigvec: smallestEigvec,
    makeMolecule: makeMolecule, poseMolecule: poseMolecule,
    makeCamera: makeCamera, makeCameraPose: makeCameraPose,
    aimAngles: aimAngles, project: project, quantize: quantize,
    backproject: backproject,
    fundamental: fundamental, epipolarResidual: epipolarResidual,
    epipolarLine: epipolarLine, triangulate: triangulate,
    planeHomography: planeHomography, applyH: applyH,
    reprojError: reprojError, depthSigma: depthSigma,
    svd3: svd3, camFromRt: camFromRt, eightPoint: eightPoint,
    essentialFromPairs: essentialFromPairs,
    epipolarRMS: epipolarRMS, relativePose: relativePose,
    rotationAngle: rotationAngle, makeBox: makeBox, makeTSDF: makeTSDF
  };
})();

/* headless test hook (node); harmless in the browser */
if (typeof module !== 'undefined' && module.exports) module.exports = MV;
