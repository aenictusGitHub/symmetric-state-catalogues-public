/* Pure symmetric six-qubit witnesses for clicked surface coordinates.
 * Sources: Julia/N6_C3_fixed_S1S2_core.jl, N6_C2_fixed_S1S2_core.jl,
 * N6_gray_search_attainment.jl, and the complete Q=0 side-face proof.
 * A failed constructor is not an exclusion. Every returned state is checked
 * with physical reduced purities by SectorMath, independently of the inverse.
 */
(function (global) {
  'use strict';
  const api = global.N6SurfaceConstructors = global.N6SurfaceConstructors || {};
  const sqrt = x => Math.sqrt(Math.max(0, x));
  const clamp = (x, lo, hi) => Math.max(lo, Math.min(hi, x));
  const zero = () => Array.from({length: 7}, () => [0, 0]);
  const cis = (a, r = 1) => [r * Math.cos(a), r * Math.sin(a)];
  const abs = z => Math.hypot(z[0], z[1]);
  const multiply = (z, a) => [z[0] * Math.cos(a) - z[1] * Math.sin(a),
    z[0] * Math.sin(a) + z[1] * Math.cos(a)];
  const number = x => Number.isFinite(x) ? Number(x.toPrecision(13)).toString() : '';
  const evenFormula = String.raw`A e^{i\alpha}D_0+B D_2+C D_4+D e^{i\delta}D_6`;

  function checked(coefficients, point, tolerance, family, parameters, latex, parameterLatex) {
    if (!coefficients || coefficients.length !== 7 ||
        coefficients.some(z => !Array.isArray(z) || z.length !== 2 || !z.every(Number.isFinite))) return null;
    const norm = coefficients.reduce((sum, z) => sum + z[0] ** 2 + z[1] ** 2, 0);
    if (!(norm > 0) || Math.abs(norm - 1) > 2e-8) return null;
    const c = coefficients.map(z => [z[0] / Math.sqrt(norm), z[1] / Math.sqrt(norm)]);
    const sectors = global.SectorMath.sectorLengths(c);
    const residual = Math.max(...point.map((v, i) => Math.abs(v - sectors[i])));
    if (!Number.isFinite(residual) || residual > tolerance) return null;
    return {coefficients: c, sectors, family, parameters: {...parameters, residual},
      analyticalExpression: {latex, parameterLatex}};
  }

  function c3(point, tol) {
    const [x, y, z] = point, v = sqrt(2.5 * (y - 5));
    if (v > 5 + tol) return null;
    const s = (4 + Math.min(v, 5)) / 9;
    if (x <= 6 * s * s + tol) {
      const d = Math.min(s, sqrt(x / 6)), h = sqrt(s * s - d * d), b = 1 - s;
      const low = 8 + 2 * x - 8 * (h + b) ** 2;
      const high = 8 + 2 * x - 8 * (h - b) ** 2;
      if (z >= low - tol && z <= high + tol) {
        const q = (8 + 2 * x - z) / 8;
        const phase = h * b < 1e-15 ? 0 : .5 * Math.acos(clamp((h * h + b * b - q) / (2 * h * b), -1, 1));
        const c = zero(); c[0][0] = sqrt((s + d) / 2); c[3] = cis(phase, sqrt(b)); c[6][0] = sqrt((s - d) / 2);
        const out = checked(c, point, tol, 'C₃ three-level family', {s, d, phase},
          String.raw`\sqrt{\frac{s+d}{2}}D_0+e^{i\phi}\sqrt{1-s}D_3+\sqrt{\frac{s-d}{2}}D_6`,
          String.raw`s\simeq ${number(s)},\quad d\simeq ${number(d)},\quad\phi\simeq ${number(phase)}`);
        if (out) return out;
      }
    }
    if (v <= 3 + tol && Math.abs(x - 2 * (2 - v) ** 2 / 3) <= tol && Math.abs(z - 8 - 2 * x) <= tol) {
      const p = clamp(1 - v / 3, 0, 1), c = zero(); c[1][0] = sqrt(p); c[4][0] = sqrt(1 - p);
      return checked(c, point, tol, 'C₃ two-level character', {p},
        String.raw`\sqrt p\,D_1+\sqrt{1-p}\,D_4`, String.raw`p\simeq ${number(p)}`);
    }
    return null;
  }

  function oddC2(point, tol) {
    const [x, y, z] = point, f = sqrt(1.5 * x), delta = f / 2, gap = 15 + 5 * x - y - z;
    if (delta > 1 + tol || gap < -tol) return null;
    const radius = sqrt(5 * Math.max(gap, 0) / 128), w = (y - 5) / 16;
    for (const beta0 of [5 / 8 - radius, 5 / 8 + radius]) {
      if (beta0 < -tol || beta0 > 1 - delta + tol) continue;
      const beta = clamp(beta0, 0, Math.max(1 - delta, 0)), r2 = (1 - beta) ** 2 - delta * delta;
      if (r2 < -tol) continue;
      const denominator = beta * sqrt(r2), numerator = w - beta + 3 * beta * beta / 5;
      if (Math.abs(numerator) > denominator + tol) continue;
      const p1 = (1 - beta + delta) / 2, p5 = (1 - beta - delta) / 2;
      const phase = denominator < 1e-15 ? 0 : .5 * Math.acos(clamp(numerator / denominator, -1, 1));
      const c = zero(); c[1][0] = sqrt(p1); c[3] = cis(phase, sqrt(beta)); c[5][0] = sqrt(p5);
      const out = checked(c, point, tol, 'Odd C₂ three-level family', {f, beta, phase},
        String.raw`\sqrt{\frac{1-\beta+f/2}{2}}D_1+e^{i\phi}\sqrt\beta\,D_3+\sqrt{\frac{1-\beta-f/2}{2}}D_5`,
        String.raw`f\simeq ${number(f)},\quad\beta\simeq ${number(beta)},\quad\phi\simeq ${number(phase)}`);
      if (out) return out;
    }
    return null;
  }

  // Necessary and sufficient two-triangle reconstruction at fixed populations.
  function evenState(x, y, s, a, q, tol) {
    const f = sqrt(1.5 * x), H = 7.5 * (y - 5), b = f - 3 * a;
    const populations = [(s + a) / 2, (1 - s + b) / 2, (1 - s - b) / 2, (s - a) / 2];
    if (Math.min(...populations) < -tol) return null;
    const [A, B, C, D] = populations.map(sqrt), T = H - 3 * (8 * s - 3) ** 2;
    if (T < -tol * Math.max(1, H)) return null;
    const r = 2 * A * D, w = 2 * B * C, h = sqrt(T), K0 = (f - 2 * a) ** 2 - (2 * s - 1) ** 2;
    const low = Math.max(0, (r - w) ** 2, (h - 6 * w) ** 2 / 15 - K0);
    const high = Math.min(1, (r + w) ** 2, (h + 6 * w) ** 2 / 15 - K0);
    if (q < low - tol || q > high + tol || low > high + tol) return null;
    const product = 2 * A * B * C * D;
    const lambda = product < 1e-16 ? 0 : Math.acos(clamp((q / 4 - A * A * D * D - B * B * C * C) / product, -1, 1));
    const U = [A * B + C * D * Math.cos(lambda), C * D * Math.sin(lambda)];
    const L = 2 * Math.sqrt(15) * abs(U), K = 12 * B * C;
    let alpha = 0;
    if (L * K > 1e-16) {
      const cosine = (T - L * L - K * K) / (2 * L * K);
      if (Math.abs(cosine) > 1 + 10 * tol) return null;
      alpha = Math.atan2(U[1], U[0]) - Math.acos(clamp(cosine, -1, 1));
    }
    const d = lambda - alpha, c = zero(); c[0] = cis(alpha, A); c[2][0] = B; c[4][0] = C; c[6] = cis(d, D);
    return {c, A, B, C, D, alpha, delta: d, s, a, q};
  }

  function evenRecord(data, point, tol, family) {
    if (!data) return null;
    const {A, B, C, D, alpha, delta} = data;
    return checked(data.c, point, tol, family, {...data, c: undefined}, evenFormula,
      String.raw`(A,B,C,D)\simeq(${[A, B, C, D].map(number).join(',')}),\quad\alpha\simeq ${number(alpha)},\quad\delta\simeq ${number(delta)}`);
  }

  function qzeroSide(point, tol) {
    const [x, y, z] = point;
    if (Math.abs(y - 5) > tol || x > 8 / 3 + tol) return null;
    const q = clamp((8 + 2 * x - z) / 8, 0, 1), f = sqrt(1.5 * x), r = Math.sqrt(q);
    // Explicit real three-mode span. It supplies the complete side for f>=1/2.
    if (q >= 1 - 4 * f * f - tol && q <= 1 - f * f / 4 + tol) {
      if (q >= 1 - 1e-15 && x < tol) {
        const c = zero(); c[0][0] = c[6][0] = Math.SQRT2 / 3; c[3][1] = Math.sqrt(5) / 3;
        const out = checked(c, point, tol, 'Quadrupole-free three-mode span', {alpha: 1, beta: 0, gamma: 0},
          String.raw`\frac{\sqrt2(D_0+D_6)+i\sqrt5D_3}{3}`, '');
        if (out) return out;
      }
      if (q < 1) {
        const A = (1 + r) / 2, lambda0 = (6 * x / (1 - q) - 1) / 15;
        if (lambda0 >= -tol && lambda0 <= 1 + tol) {
          const lambda = clamp(lambda0, 0, 1), alpha = sqrt(A);
          const beta = sqrt(1 - A) * (sqrt(1 - lambda) + 2 * sqrt(2 * lambda)) / 3;
          const gamma = sqrt(1 - A) * (2 * sqrt(2 * (1 - lambda)) - sqrt(lambda)) / 3;
          const c = zero();
          c[0][0] = alpha * Math.SQRT2 / 3 + beta * Math.sqrt(3) / 4;
          c[6][0] = alpha * Math.SQRT2 / 3 - beta * Math.sqrt(3) / 4;
          c[2][0] = beta * Math.sqrt(5) / 4; c[4][0] = -beta * Math.sqrt(5) / 4;
          c[3][1] = alpha * Math.sqrt(5) / 3; c[1][1] = -gamma / Math.SQRT2; c[5][1] = gamma / Math.SQRT2;
          const out = checked(c, point, tol, 'Quadrupole-free three-mode span', {alpha, beta, gamma},
            String.raw`\alpha e_0+\beta e_1+\gamma e_2`,
            String.raw`e_0=\frac{\sqrt2(D_0+D_6)+i\sqrt5D_3}{3},\quad e_1=\frac{\sqrt3(D_0-D_6)+\sqrt5(D_2-D_4)}4,\quad e_2=-\frac{i(D_1-D_5)}{\sqrt2},\quad(\alpha,\beta,\gamma)\simeq(${[alpha, beta, gamma].map(number).join(',')})`);
          if (out) return out;
        }
      }
    }
    // The missing small-dipole region has a quadratic population inverse.
    const rad = (5 / 128) * (1 - 3 * f * f / 32 - q);
    if (rad >= -tol) for (const sign of [-1, 1]) {
      const a = 23 * f / 64 + sign * sqrt(rad);
      const out = evenRecord(evenState(x, 5, 3 / 8, a, q, tol), point, tol, 'Quadrupole-free even C₂ family');
      if (out) return out;
    }
    return null;
  }

  function realC2(f, r, e) {
    if (r >= 1 || r < 0) return null;
    const g = sqrt(1 - r * r), nu = f / g - 2 * e;
    if (Math.abs(e) > 1 + 1e-12 || Math.abs(nu) > 1 + 1e-12) return null;
    const A = sqrt((1 + g) / 2), B = r / sqrt(2 * (1 + g));
    const ca = sqrt((1 + e) / 2), sa = sqrt((1 - e) / 2), cb = sqrt((1 + nu) / 2), sb = sqrt((1 - nu) / 2);
    return [[A * ca * cb + B * sa * sb, 0], [0, 0], [A * ca * sb - B * sa * cb, 0], [0, 0],
      [A * sa * cb - B * ca * sb, 0], [0, 0], [-A * sa * sb - B * ca * cb, 0]];
  }
  const rc = Math.sqrt((648 * Math.sqrt(5) - 1445) / 931);
  function lowStationary(r) {
    function derivative(e) {
      const ss = sqrt(1 - e * e), tt = sqrt(1 - 4 * e * e), st = ss * tt;
      const sp = (-10 * e + 16 * e ** 3) / (2 * st);
      const z = 1 - 8 * e * e + 4 * r * st, ze = -16 * e + 4 * r * sp;
      const h = Math.sqrt(15) * (e * tt + 2 * r * e * ss) + 3 * (st - r * (1 + 2 * e * e));
      const he = Math.sqrt(15) * (tt - 4 * e * e / tt + 2 * r * (ss - e * e / ss)) + 3 * (sp - 4 * r * e);
      return 3 * z * ze + h * he;
    }
    let lo = .116, hi = .130;
    for (let k = 0; k < 58; k++) { const e = (lo + hi) / 2; if (derivative(e) > 0) lo = e; else hi = e; }
    const e = (lo + hi) / 2, c = realC2(0, r, e), z = 1 - 8 * e * e + 4 * r * sqrt((1 - e * e) * (1 - 4 * e * e));
    const p = 2 * Math.sqrt(15) * (c[0][0] * c[2][0] + c[4][0] * c[6][0]), t = 12 * c[2][0] * c[4][0];
    return {c, e, z, p, t, tau: 1.5 * z * z + (p + t) ** 2 / 2};
  }

  function zeroDipoleSide(point, tol) {
    const [x, y, z] = point;
    if (Math.abs(x) > tol || z < -tol || z > 8 + tol) return null;
    const r = sqrt(clamp((8 - z) / 8, 0, 1)), low = 5 + Math.max(1 - 9 * r, 0) ** 2 / 10;
    const data = r < rc ? lowStationary(r) : null;
    const high = data ? 5 + 4 * data.tau / 15 : r <= 1 / 7 ? 239 / 35 + 64 * r * r / 5 :
      r <= 1 / 4 ? 33 / 5 + 16 * r / 5 + 8 * r * r / 5 : 7 + 8 * r * r;
    if (y < low - tol || y > high + tol) return null;
    const tau = 15 * (y - 5) / 4;
    if (data && tau >= 1.5 * data.z ** 2 + (data.p - data.t) ** 2 / 2 - 2e-11) {
      const theta = Math.acos(clamp((tau - 1.5 * data.z ** 2 - (data.p ** 2 + data.t ** 2) / 2) / (data.p * data.t), -1, 1));
      const c = data.c.map(v => v.slice()); c[4] = multiply(c[4], theta); c[6] = multiply(c[6], theta);
      const out = checked(c, point, tol, 'Zero-dipole low C₂ phase family', {r, e: data.e, theta},
        String.raw`aD_0+bD_2+e^{i\theta}(cD_4+dD_6)`,
        String.raw`(a,b,c,d)\simeq(${[0, 2, 4, 6].map(i => number(data.c[i][0])).join(',')}),\quad\theta\simeq ${number(theta)}`);
      if (out) return out;
    }
    if (r >= 1 / 7 - 1e-14 && tau >= 6 * (1 - r) ** 2 - 2e-11) {
      const beta = clamp((15 - sqrt(45 - 6 * tau + 180 * r * r)) / 24, 0, 5 / 8);
      const phase = beta < 1e-15 ? 0 : .5 * Math.acos(clamp((r * r - (1 - beta) ** 2 - beta * beta) / (2 * beta * (1 - beta)), -1, 1));
      const c = zero(); c[1][0] = c[5][0] = sqrt((1 - beta) / 2); c[3] = cis(phase, sqrt(beta));
      const out = checked(c, point, tol, 'Zero-dipole odd C₂ family', {r, beta, phase},
        String.raw`\sqrt{\frac{1-\beta}{2}}(D_1+D_5)+e^{i\phi}\sqrt\beta D_3`,
        String.raw`\beta\simeq ${number(beta)},\quad\phi\simeq ${number(phase)}`);
      if (out) return out;
    }
    if (r <= 1 / 7 + 1e-14 && tau >= 3351 / 784 + 363 * r * r / 16 - 2e-11) {
      const w2 = (243 / 128) * (5 / 49 + r * r), a = 10 * sqrt(w2) / 9, b = 6 * sqrt(w2);
      const theta = Math.acos(clamp((tau - 96 / 49 - (a * a + b * b) / 2) / (a * b), -1, 1));
      const s = 29 / 56, d = sqrt((27 / 128) * (1 / 49 - r * r));
      const A = (s + d) / 2, D = (s - d) / 2, B = (1 - s - 3 * d) / 2, C = (1 - s + 3 * d) / 2;
      const v = 2 * sqrt(A * D), w = 2 * sqrt(B * C), delta = Math.acos(clamp((r * r - v * v - w * w) / (2 * v * w), -1, 1));
      const alpha = Math.atan2(sqrt(C * D) * Math.sin(delta), sqrt(A * B) + sqrt(C * D) * Math.cos(delta)) + theta;
      const c = zero(); c[0][0] = sqrt(A); c[2][0] = sqrt(B); c[4] = cis(alpha, sqrt(C)); c[6] = cis(alpha + delta, sqrt(D));
      const out = checked(c, point, tol, 'Zero-dipole middle C₂ phase family', {r, theta, alpha, delta},
        String.raw`\sqrt A D_0+\sqrt B D_2+e^{i\alpha}\sqrt C D_4+e^{i(\alpha+\delta)}\sqrt D D_6`,
        String.raw`(A,B,C,D)\simeq(${[A, B, C, D].map(number).join(',')}),\quad\alpha\simeq ${number(alpha)},\quad\delta\simeq ${number(delta)}`);
      if (out) return out;
    }
    if (tau >= 3 * (1 - 9 * r) ** 2 / 8 - 2e-11) {
      const top = 6 - 3 * r + 57 * r * r / 2, a = 15 * (1 + r) ** 2 / 8, b = -7.5 * (1 + r), cc = top - tau;
      const t = clamp(2 * cc / (-b + sqrt(b * b - 4 * a * cc)), 0, 1), c = zero();
      c[0][0] = c[6][0] = sqrt((1 - r) / 4); c[2][1] = c[4][1] = sqrt((1 + r) * (1 - t) / 4); c[3][0] = -sqrt((1 + r) * t / 2);
      const out = checked(c, point, tol, 'Zero-dipole five-level half-turn family', {r, t},
        String.raw`\frac{\sqrt{1-r}}2(D_0+D_6)+\frac{i\sqrt{(1+r)(1-t)}}2(D_2+D_4)-\sqrt{\frac{(1+r)t}{2}}D_3`,
        String.raw`r\simeq ${number(r)},\quad t\simeq ${number(t)}`);
      if (out) return out;
    }
    return null; // Remaining lower C3 branch was already tested by c3().
  }

  function stationaryP(u, f) {
    const d = 1 + 15 * u * u;
    return (3 * u + 1) * (5 * u - 1) * d * d * (15 * u * u - 30 * u - 1) * f * f +
      4 * d * (3375 * u ** 6 - 8325 * u ** 5 - 2475 * u ** 4 + 1410 * u ** 3 + 165 * u * u - 37 * u - 1) * f +
      (15 * u * u + 6 * u - 1) * (45 * u * u - 10 * u - 3) * (225 * u ** 4 - 900 * u ** 3 - 90 * u * u + 60 * u + 1);
  }
  function stationaryRoot(f, lo, hi) {
    let fl = stationaryP(lo, f);
    if (fl * stationaryP(hi, f) > 0) return null;
    for (let k = 0; k < 65; k++) { const mid = (lo + hi) / 2, fm = stationaryP(mid, f); if (fm * fl > 0) { lo = mid; fl = fm; } else hi = mid; }
    return (lo + hi) / 2;
  }
  function nullInvariants(u, f) {
    const d = 1 + 15 * u * u, e = (1 - 15 * u * u) / d, nu = f - 2 * e, H = (1 + 6 * u - 15 * u * u) / d;
    return {e, nu, tau: 1.5 * (1 + 4 * e * nu) ** 2 + 7.5 * (1 - nu * nu) * H * H};
  }
  function realNullBand(point, tol) {
    const [x, y, z] = point, f = sqrt(1.5 * x);
    if (Math.abs(z - 8 - 2 * x) > tol || f > .366379159399607379 + 1e-14) return null;
    const ul = stationaryRoot(f, .35, .4), uu = stationaryRoot(f, .19, .23);
    if (ul === null || uu === null) return null;
    const target = 15 * (y - 5) / 4, low = nullInvariants(ul, f).tau, high = nullInvariants(uu, f).tau;
    if (target < low - 4 * tol || target > high + 4 * tol) return null;
    let lo = uu, hi = ul;
    for (let k = 0; k < 65; k++) { const mid = (lo + hi) / 2; if (nullInvariants(mid, f).tau > target) lo = mid; else hi = mid; }
    const u = (lo + hi) / 2, nu = nullInvariants(u, f).nu;
    if (Math.abs(nu) > 1 + tol) return null;
    const ca = sqrt((1 + nu) / 2), sa = sqrt((1 - nu) / 2), N = sqrt(1 + 15 * u * u), c = zero();
    c[0][0] = ca / N; c[2][0] = sa / N; c[4][0] = Math.sqrt(15) * u * ca / N; c[6][0] = -Math.sqrt(15) * u * sa / N;
    return checked(c, point, tol, 'Real C₂ null-roof band', {f, u, nu},
      String.raw`\frac{cD_0+sD_2+\sqrt{15}u(cD_4-sD_6)}{\sqrt{1+15u^2}}`,
      String.raw`u\simeq ${number(u)},\quad\nu\simeq ${number(nu)},\quad c=\sqrt{\frac{1+\nu}{2}},\quad s=\sqrt{\frac{1-\nu}{2}}`);
  }

  const addPoly = (a, b) => Array.from({length: Math.max(a.length, b.length)}, (_, i) => (a[i] || 0) + (b[i] || 0));
  const scalePoly = (a, k) => a.map(v => v * k);
  function mulPoly(a, b) { const c = Array(a.length + b.length - 1).fill(0); a.forEach((x, i) => b.forEach((y, j) => c[i + j] += x * y)); return c; }
  const evalPoly = (a, x) => a.reduceRight((v, c) => v * x + c, 0);
  // Derivative isolation partitions a bounded real interval into monotone pieces.
  // Endpoint/even-multiplicity roots are retained; final physical verification
  // prevents a small polynomial residual from certifying a false state.
  function realRoots(coefficients, lo, hi) {
    const magnitude = Math.max(...coefficients.map(Math.abs));
    if (!(magnitude > 0)) return [];
    const a = coefficients.map(v => v / magnitude);
    while (a.length > 1 && Math.abs(a[a.length - 1]) < 1e-15) a.pop();
    if (a.length === 1) return [];
    if (a.length === 2) { const r = -a[0] / a[1]; return r >= lo - 1e-12 && r <= hi + 1e-12 ? [clamp(r, lo, hi)] : []; }
    const critical = realRoots(a.slice(1).map((v, i) => v * (i + 1)), lo, hi), cuts = [lo, ...critical, hi].sort((x, y) => x - y), roots = [];
    for (const x of cuts) if (Math.abs(evalPoly(a, x)) < 2e-11) roots.push(x);
    for (let i = 1; i < cuts.length; i++) {
      let l = cuts[i - 1], r = cuts[i], fl = evalPoly(a, l), fr = evalPoly(a, r);
      if (fl * fr >= 0) continue;
      for (let k = 0; k < 65; k++) { const m = (l + r) / 2, fm = evalPoly(a, m); if (fl * fm > 0) { l = m; fl = fm; } else r = m; }
      roots.push((l + r) / 2);
    }
    return roots.sort((x, y) => x - y).filter((v, i, arr) => !i || v - arr[i - 1] > 1e-10);
  }

  function complexNull(point, tol) {
    const [x, y, z] = point;
    if (Math.abs(z - 8 - 2 * x) > tol) return null;
    const f = sqrt(1.5 * x), H = 7.5 * (y - 5), lo = Math.max(-1, (f - 1) / 2), hi = Math.min(1, (f + 1) / 2);
    if (lo > hi) return null;
    const v = [f, -2], u = [0, 1], v2 = addPoly([1], scalePoly(mulPoly(v, v), -1));
    const uv = addPoly([1], scalePoly(mulPoly(u, v), 4));
    const Z = addPoly([H], addPoly(scalePoly(mulPoly(uv, uv), -3), scalePoly(mulPoly(v2, [9, 0, 6]), -1)));
    const rhs = scalePoly(mulPoly([0, 0, 1, 0, -1], mulPoly(v2, v2)), 540);
    const roots = realRoots(addPoly(rhs, scalePoly(mulPoly(Z, Z), -1)), lo, hi), cuts = [lo, ...roots, hi].sort((a, b) => a - b);
    const candidates = cuts.slice(); for (let i = 1; i < cuts.length; i++) candidates.push((cuts[i - 1] + cuts[i]) / 2);
    for (const e of candidates) {
      const nu = clamp(f - 2 * e, -1, 1), value = H - 3 * (1 + 4 * e * nu) ** 2 - (1 - nu * nu) * (9 + 6 * e * e);
      const denominator = 6 * Math.sqrt(15) * e * sqrt(1 - e * e) * (1 - nu * nu);
      if (Math.abs(denominator) < 1e-15 ? Math.abs(value) > 1e-9 : Math.abs(value / denominator) > 1 + 1e-9) continue;
      const chi = Math.abs(denominator) < 1e-15 ? 0 : Math.acos(clamp(value / denominator, -1, 1));
      const ca = sqrt((1 + e) / 2), sa = sqrt((1 - e) / 2), cb = sqrt((1 + nu) / 2), sb = sqrt((1 - nu) / 2), c = zero();
      c[0][0] = ca * cb; c[2][0] = ca * sb; c[4] = cis(chi, sa * cb); c[6] = cis(chi + Math.PI, sa * sb);
      const out = checked(c, point, tol, 'Complex C₂ null-roof family', {f, e, nu, chi},
        String.raw`c_a c_bD_0+c_a s_bD_2+e^{i\chi}s_a(c_bD_4-s_bD_6)`,
        String.raw`e\simeq ${number(e)},\quad\nu\simeq ${number(nu)},\quad\chi\simeq ${number(chi)},\quad c_a=\sqrt{\frac{1+e}{2}},\ s_a=\sqrt{\frac{1-e}{2}},\ c_b=\sqrt{\frac{1+\nu}{2}},\ s_b=\sqrt{\frac{1-\nu}{2}}`);
      if (out) return out;
    }
    return null;
  }

  api.cyclic = function cyclic(point, tolerance = 1e-7) {
    if (!Array.isArray(point) || point.length !== 3 || !point.every(Number.isFinite) ||
        !Number.isFinite(tolerance) || tolerance <= 0 || !global.SectorMath) return null;
    const [x, y, z] = point;
    if (x < -tolerance || x > 6 + tolerance || y < 5 - tolerance || y > 15 + tolerance || z < -tolerance || z > 20 + tolerance) return null;
    for (const constructor of [c3, oddC2, qzeroSide, zeroDipoleSide, realNullBand, complexNull]) {
      const result = constructor(point, tolerance);
      if (result) return result;
    }
    return null;
  };
})(globalThis);
