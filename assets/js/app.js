(() => {
'use strict';
// ================================================================== physics
const C = { g: 9.81, Patm: 101325, rhoW: 1000, rhoA: 1.2, gam: 1.4, R: 287, T0: 293 };
const AIR_END = 1.0005;   // the air phase ends when the chamber is this close to ambient pressure
const FINS = {
  kucuk: { cp: 0.66, cna: 4.0, mass: 0.03, cd: 0.03, span: 0.05, chord: 0.12, name: 'Küçük' },
  orta:  { cp: 0.74, cna: 7.0, mass: 0.06, cd: 0.05, span: 0.075, chord: 0.14, name: 'Orta' },
  buyuk: { cp: 0.82, cna: 11.0, mass: 0.10, cd: 0.08, span: 0.105, chord: 0.17, name: 'Büyük' },
};
const WING_RHO = 0.6;   // kg/m², foam wing with a spar
const tankLen = (VL, D) => (VL / 1000) / (Math.PI * (D / 2) ** 2) + 0.06;
const wrapA = (a) => Math.atan2(Math.sin(a), Math.cos(a));

function geometry(p) {
  const fin = FINS[p.fins];
  const D = (p.Dmm || 100) / 1000, R = D / 2;
  const Ln = 0.15, Lp = tankLen(2 * p.payload, D), Lc = tankLen(p.Vc, D), Ltail = 0.04;
  const L = Ln + Lp + Lc + Ltail;
  const xPay = Ln + Lp / 2, xChEnd = Ln + Lp + Lc, xDry = 0.55 * L, xCp = fin.cp * L;
  const Vc = p.Vc / 1000, Vw0 = Vc * p.wf;
  const mPay = p.payload, mProp = Vw0 * C.rhoW;
  const wing = p.wing === 'sabit' || p.wing === 'acilir' ? p.wing : null;
  const Sw = wing ? p.wspan * p.wchord : 0, AR = wing ? p.wspan * p.wspan / Sw : 0;
  const mStruct = p.dry + fin.mass;
  const xW = (mStruct * xDry + mPay * xPay) / (mStruct + mPay);   // wing aerodynamic centre placed at the empty CG
  // a winged design needs a horizontal tail; size it for a tail volume of 0.45 on the arm to the fins
  const lt = Math.max(0.15, (L - 0.06) - xW), St = wing ? 0.45 * Sw * p.wchord / lt : 0;
  const mWing = wing ? WING_RHO * (Sw + St) + (wing === 'acilir' ? 0.04 : 0) : 0;
  const mDry = mStruct + mWing;
  const cg = (Vw) => {
    const mw = Vw * C.rhoW, lw = Lc * (Vw / Vc), xw = xChEnd - lw / 2;
    return (mStruct * xDry + mWing * xW + mPay * xPay + mw * xw) / (mDry + mPay + mw);
  };
  const xcg0 = cg(Vw0), xcg1 = cg(0);
  const An = Math.PI * (p.nozzle / 2000) ** 2, Ab = Math.PI * R * R, cd = p.cd || 0.96;
  // neutral point: body + fins at xCp, horizontal tail at xT, wing at xW (only when open)
  const xT = L - 0.06, CLa = wing ? 2 * Math.PI * AR / (AR + 2) : 0, CLaT = 2 * Math.PI * 4 / 6;
  const kB = Ab * (2.0 + fin.cna), kT = St * CLaT, kW = Sw * CLa;
  const xNP = (ew) => (kB * xCp + kT * xT + ew * kW * xW) / (kB + kT + ew * kW);
  const ew0 = wing === 'sabit' ? 1 : 0, np0 = xNP(ew0);
  const Tw0 = 2 * cd * cd * An * p.P * 1e5;     // F = ṁ·v = ρ·A_n·(C_d·√(2ΔP/ρ))² = 2·C_d²·A_n·ΔP
  const Atube = Math.PI * ((p.nozzle - 2) / 2000) ** 2;
  const T0 = p.tube ? p.P * 1e5 * Atube : Tw0;   // with a launch tube the first push is ΔP·A_tube
  return { fin, D, R, Ln, Lp, Lc, Ltail, L, xPay, xChEnd, xDry, xCp, Vc, Vw0, mPay, mDry, mStruct, mWing, mProp, TOW: mDry + mPay + mProp, cg, xcg0, xcg1,
    sm0: (np0 - xcg0) / D, sm1: (np0 - xcg1) / D, smW: wing ? (xNP(1) - xcg1) / D : null, xNP0: np0, xNP1: wing ? xNP(1) : np0, An, Ab, cd, CD0: (p.CD == null ? 0.42 : p.CD) + fin.cd, T0, Tw0, Atube, tw: T0 / ((mDry + mPay + mProp) * C.g),
    wing, Sw, AR, CLa, xW, winc: (p.winc || 0) * Math.PI / 180, chord: p.wchord || 0,
    St, xT, CLaT };
}

// small 3-vector helpers (arrays)
const v3n = (a) => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const v3dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const v3cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
function v3slerp(a, b, f) {
  // rotate unit vector a toward unit vector b by fraction f of the angle between them
  const d = Math.max(-1, Math.min(1, v3dot(a, b))), ang = Math.acos(d);
  if (ang < 1e-7 || f <= 0) return a.slice();
  let ax = v3cross(a, b);
  if (Math.hypot(ax[0], ax[1], ax[2]) < 1e-9) ax = Math.abs(a[2]) < 0.9 ? v3cross(a, [0, 0, 1]) : v3cross(a, [1, 0, 0]);
  ax = v3n(ax);
  const th = ang * Math.min(1, f), c = Math.cos(th), s = Math.sin(th), k = v3cross(ax, a), kd = v3dot(ax, a);
  return v3n([a[0] * c + k[0] * s + ax[0] * kd * (1 - c), a[1] * c + k[1] * s + ax[1] * kd * (1 - c), a[2] * c + k[2] * s + ax[2] * kd * (1 - c)]);
}

// Flight model: x downrange, y up, z to the right of the launch direction.
// Rigid body with 5 degrees of freedom (position + pitch/yaw, no roll); the sloshing water is a
// spring-mass pendulum that can swing in both lateral directions.
function simulate(p, record) {
  const G = geometry(p);
  const fin = G.fin;
  const S = G.Ab;
  const Vc = G.Vc, Vw0 = G.Vw0, Va0 = Vc - Vw0;
  const P0 = C.Patm + p.P * 1e5;
  const An = G.An;
  const mPay = G.mPay, mDry = G.mDry, L = G.L, xPay = G.xPay, xCp = G.xCp;
  const Cd0 = G.CD0, CdSide = 1.2 * (L * G.D) / S, cnaTot = 2.0 + fin.cna, cdN = G.cd;
  const Rt = G.R, xi = 1.841;
  const lerpS = (a, c, f) => a + (c - a) * f;
  let tWing = null;
  const wingFrac = (t) => (!G.wing ? 0 : G.wing === 'sabit' ? 1 : tWing == null ? 0 : Math.min(1, (t - tWing) / 0.35));
  const h = (p.payload / 1000) / (Math.PI * Rt * Rt);
  const kappa = (xi / Rt) * Math.tanh(xi * h / Rt);
  const mu = Math.min(0.6, (2 * Rt / (h * xi * (xi * xi - 1))) * Math.tanh(xi * h / Rt));
  const b = p.baffle;
  const mS = mu * mPay * (1 - 0.6 * Math.min(b, 1));
  const zeta = 0.008 + 0.25 * Math.pow(b, 1.2);
  const sMax = 0.4 * Rt;
  const tubeLen = p.tube ? 0.30 : 0;
  const Atube = Math.PI * ((p.nozzle - 2) / 2000) ** 2;
  const railLen = 1.0;
  const th0 = Math.min(p.angle, 89.95) * Math.PI / 180;   // exactly vertical never turns over in still air
  const windZ = p.windZ || 0;

  let phase = tubeLen > 0 ? 0 : 1;
  let Pstart = P0, VaStart = Va0, P1 = 0, Tair1 = 0, ma1 = 0;
  let chute = false, chuteA = 0, onRail = true;
  const s0 = -(Math.cos(th0)) / (kappa * Math.sin(th0));
  let st = { x: 0, y: 0, z: 0, vx: 0, vy: 0, vz: 0, bx: Math.cos(th0), by: Math.sin(th0), bz: 0, wx: 0, wy: 0, wz: 0,
    sx: -Math.sin(th0) * s0, sy: Math.cos(th0) * s0, sz: 0, ux: 0, uy: 0, uz: 0, Vw: Vw0, ma: 0, xt: 0 };
  // smoothed direction of apparent gravity (world frame) for the payload tank and the chamber
  let gS = [0, -1, 0], gC = [0, -1, 0], thp = th0;

  const maAir0 = P0 * Va0 / (C.R * C.T0);     // compressed air in the chamber (constant until the water is gone)
  const airState = (ma) => { const r = ma / ma1; return { P: P1 * Math.pow(r, C.gam), T: Tair1 * Math.pow(r, C.gam - 1) }; };
  function thrustOf(s_, ph) {
    if (ph === 0) {
      const Va = Va0 + Atube * s_.xt, P = P0 * Math.pow(Va0 / Va, C.gam);
      return { T: Math.max(0, (P - C.Patm) * Atube), P, dVw: 0, dma: 0 };
    }
    if (ph === 1) {
      const Va = Vc - s_.Vw, P = Pstart * Math.pow(VaStart / Va, C.gam);
      // v_su = C_d·√(2ΔP/ρ),  ṁ = ρ·A_n·v_su,  F = ṁ·v_su
      const ve = cdN * Math.sqrt(2 * Math.max(0, P - C.Patm) / C.rhoW);
      return { T: C.rhoW * An * ve * ve, P, dVw: -An * ve, dma: 0 };
    }
    if (ph === 2) {
      const { P, T } = airState(Math.max(s_.ma, 1e-9));
      if (P <= C.Patm * AIR_END) return { T: 0, P, dVw: 0, dma: 0 };
      const g = C.gam; let mdot, Th;
      if (P / C.Patm > 1.893) {
        const Pe = P * 0.5283, Te = T * 0.8333, ve = Math.sqrt(g * C.R * Te), rhoE = Pe / (C.R * Te);
        mdot = cdN * rhoE * An * ve; Th = mdot * ve + (Pe - C.Patm) * An;
      } else {
        const ratio = Math.pow(C.Patm / P, (g - 1) / g);
        const ve = Math.sqrt(Math.max(0, 2 * g / (g - 1) * C.R * T * (1 - ratio)));
        const rhoE = C.Patm / (C.R * T * ratio);
        mdot = cdN * rhoE * An * ve; Th = mdot * ve;
      }
      return { T: Th, P, dVw: 0, dma: -mdot };
    }
    return { T: 0, P: C.Patm, dVw: 0, dma: 0 };
  }
  const windAt = (t, y) => (y < 0.3 ? [0, 0] : [
    p.wind * (1 + 0.25 * Math.sin(1.3 * t) + 0.15 * Math.sin(3.7 * t + 1)),
    windZ * (1 + 0.25 * Math.sin(1.1 * t + 2) + 0.15 * Math.sin(3.1 * t + 0.5)),
  ]);

  function deriv(s_, t, ph) {
    const mw = s_.Vw * C.rhoW, ma = ph === 2 ? s_.ma : ph <= 1 ? maAir0 : 0;
    const m = mDry + mPay + mw + ma;
    const xcg = G.cg(s_.Vw);
    const Lt = xcg - xPay, lCp = xCp - xcg, lFin = Math.max(0.05, L - 0.06 - xcg);
    const I = m * L * L / 12;
    const tr = thrustOf(s_, ph);
    const bl = Math.hypot(s_.bx, s_.by, s_.bz), bx = s_.bx / bl, by = s_.by / bl, bz = s_.bz / bl;
    const w = windAt(t, s_.y);
    const rvx = s_.vx - w[0], rvy = s_.vy, rvz = s_.vz - w[1], V = Math.hypot(rvx, rvy, rvz);
    let Fx = tr.T * bx, Fy = tr.T * by, Fz = tr.T * bz, Mx = 0, My = 0, Mz = 0, alpha = 0, lift = 0;
    if (V > 0.3) {
      const va = rvx * bx + rvy * by + rvz * bz;
      const px = rvx - va * bx, py = rvy - va * by, pz = rvz - va * bz, vp = Math.hypot(px, py, pz);
      alpha = Math.atan2(vp, va);
      const qd = 0.5 * C.rhoA * V * V, ca = va / V, sa = vp / V;
      let Dm = qd * S * (Cd0 * ca * ca + CdSide * sa * sa);
      if (chute) Dm += qd * 1.5 * chuteA;
      Fx -= Dm * rvx / V; Fy -= Dm * rvy / V; Fz -= Dm * rvz / V;
      if (!chute) {
        const kn = qd * S * cnaTot / V;            // normal force = -kn * v_perp
        Fx -= kn * px; Fy -= kn * py; Fz -= kn * pz;
        const c = lCp * kn, cd = 0.5 * C.rhoA * V * S * fin.cna * lFin * lFin;
        Mx += c * (by * pz - bz * py) - cd * s_.wx;
        My += c * (bz * px - bx * pz) - cd * s_.wy;
        Mz += c * (bx * py - by * px) - cd * s_.wz;
      }
      const ew = chute ? 0 : wingFrac(t);
      if (G.wing && !chute) {
        let ux = -by * bx, uy = 1 - by * by, uz = -by * bz, ul = Math.hypot(ux, uy, uz);
        if (ul < 1e-4) { ux = -Math.sin(th0); uy = Math.cos(th0); uz = 0; ul = 1; }
        ux /= ul; uy /= ul; uz /= ul;
        const at = Math.atan2(-(rvx * ux + rvy * uy + rvz * uz), va);
        const e = Math.asin(Math.min(1, Math.abs(Math.sin(at)))), st_ = Math.min(1, Math.max(0, (e - 0.26) / 0.35));
        const sct = Math.sin(at) * Math.abs(Math.cos(at));   // lift is measured toward the wing's upper side, so reversed flow keeps the sign of sin
        const CLt = lerpS(G.CLaT * sct, 1.2 * sct, st_), CDt = 0.015 + 1.2 * Math.pow(Math.sin(at), 2) * st_;
        const qt = 0.5 * C.rhoA * V * V * G.St;
        const vhx = rvx / V, vhy = rvy / V, vhz = rvz / V, ud = ux * vhx + uy * vhy + uz * vhz;
        let lx = ux - ud * vhx, ly = uy - ud * vhy, lz = uz - ud * vhz; const ll = Math.hypot(lx, ly, lz) || 1; lx /= ll; ly /= ll; lz /= ll;
        const ftx = qt * (CLt * lx - CDt * vhx), fty = qt * (CLt * ly - CDt * vhy), ftz = qt * (CLt * lz - CDt * vhz);
        Fx += ftx; Fy += fty; Fz += ftz;
        const levT = xcg - G.xT, cdt = 0.5 * C.rhoA * V * G.St * G.CLaT * levT * levT;
        Mx += levT * (by * ftz - bz * fty) - cdt * s_.wx;
        My += levT * (bz * ftx - bx * ftz) - cdt * s_.wy;
        Mz += levT * (bx * fty - by * ftx) - cdt * s_.wz;
      }
      if (ew > 0) {
        // wing kept level by dihedral: its "up" is world-up made perpendicular to the body axis
        let ux = -by * bx, uy = 1 - by * by, uz = -by * bz, ul = Math.hypot(ux, uy, uz);
        if (ul < 1e-4) { ux = -Math.sin(th0); uy = Math.cos(th0); uz = 0; ul = 1; }
        ux /= ul; uy /= ul; uz /= ul;
        const aw = Math.atan2(-(rvx * ux + rvy * uy + rvz * uz), va) + G.winc;
        // attached flow up to ~12°, then a flat plate (also covers reversed flow in a tail slide)
        const e = Math.asin(Math.min(1, Math.abs(Math.sin(aw)))), st_ = Math.min(1, Math.max(0, (e - 0.21) / 0.35));
        const scw = Math.sin(aw) * Math.abs(Math.cos(aw)), CLatt = G.CLa * scw;
        const CL = lerpS(CLatt, 1.2 * scw, st_);
        const CDw = 0.02 + (1 - st_) * CLatt * CLatt / (Math.PI * 0.8 * G.AR) + 1.2 * Math.pow(Math.sin(aw), 2) * st_;
        const qw = 0.5 * C.rhoA * V * V * G.Sw * ew;
        const vhx = rvx / V, vhy = rvy / V, vhz = rvz / V, ud = ux * vhx + uy * vhy + uz * vhz;
        let lx = ux - ud * vhx, ly = uy - ud * vhy, lz = uz - ud * vhz; const ll = Math.hypot(lx, ly, lz) || 1; lx /= ll; ly /= ll; lz /= ll;
        lift = qw * CL; const Dw = qw * CDw;
        const fwx = lift * lx - Dw * vhx, fwy = lift * ly - Dw * vhy, fwz = lift * lz - Dw * vhz;
        Fx += fwx; Fy += fwy; Fz += fwz;
        const lev = xcg - (G.xW + 0.25 * G.chord * st_), cdw = 0.5 * C.rhoA * V * G.Sw * ew * G.CLa * lev * lev;
        Mx += lev * (by * fwz - bz * fwy) - cdw * s_.wx;
        My += lev * (bz * fwx - bx * fwz) - cdw * s_.wy;
        Mz += lev * (bx * fwy - by * fwx) - cdw * s_.wz;
      }
    }
    const fx = Fx / m, fy = Fy / m, fz = Fz / m;
    let aAx = fx * bx + fy * by + fz * bz, aAl = 0, fpx, fpy, fpz;
    if (onRail) {
      const vAl = s_.vx * bx + s_.vy * by + s_.vz * bz, sAl = s_.x * bx + s_.y * by + s_.z * bz;
      aAl = aAx - C.g * by;
      if (sAl <= 1e-6 && vAl <= 0 && aAl < 0) aAl = 0;      // resting on the launcher; elsewhere it slides back
      fpx = -C.g * by * bx; fpy = C.g - C.g * by * by; fpz = -C.g * by * bz;   // rail cancels the normal part of gravity
      aAx = aAl + C.g * by;
    } else {
      fpx = fx - aAx * bx; fpy = fy - aAx * by; fpz = fz - aAx * bz;
    }
    const w2 = kappa * Math.max(Math.abs(aAx), 0.3), ws = Math.sqrt(w2);
    const Fsx = mS * (w2 * s_.sx + 2 * zeta * ws * s_.ux), Fsy = mS * (w2 * s_.sy + 2 * zeta * ws * s_.uy), Fsz = mS * (w2 * s_.sz + 2 * zeta * ws * s_.uz);
    const kT = mS / m * tr.T;
    Mx += Lt * (by * Fsz - bz * Fsy) + kT * (by * s_.sz - bz * s_.sy);
    My += Lt * (bz * Fsx - bx * Fsz) + kT * (bz * s_.sx - bx * s_.sz);
    Mz += Lt * (bx * Fsy - by * Fsx) + kT * (bx * s_.sy - by * s_.sx);
    let alx = 0, aly = 0, alz = 0;
    const d = { x: s_.vx, y: s_.vy, z: s_.vz, vx: 0, vy: 0, vz: 0, bx: 0, by: 0, bz: 0, wx: 0, wy: 0, wz: 0,
      sx: s_.ux, sy: s_.uy, sz: s_.uz, ux: 0, uy: 0, uz: 0, Vw: tr.dVw, ma: tr.dma, xt: 0 };
    if (onRail) {
      d.vx = aAl * bx; d.vy = aAl * by; d.vz = aAl * bz;
    } else {
      d.vx = fx; d.vy = fy - C.g; d.vz = fz;
      if (!chute) {
        const mb = Mx * bx + My * by + Mz * bz;   // no roll
        alx = (Mx - mb * bx) / I; aly = (My - mb * by) / I; alz = (Mz - mb * bz) / I;
        d.bx = s_.wy * bz - s_.wz * by; d.by = s_.wz * bx - s_.wx * bz; d.bz = s_.wx * by - s_.wy * bx;
        d.wx = alx; d.wy = aly; d.wz = alz;
      }
    }
    d.xt = ph === 0 ? (s_.vx * bx + s_.vy * by + s_.vz * bz) : 0;
    const cx = aly * bz - alz * by, cy = alz * bx - alx * bz, cz = alx * by - aly * bx;
    d.ux = -w2 * s_.sx - 2 * zeta * ws * s_.ux - fpx - Lt * cx;
    d.uy = -w2 * s_.sy - 2 * zeta * ws * s_.uy - fpy - Lt * cy;
    d.uz = -w2 * s_.sz - 2 * zeta * ws * s_.uz - fpz - Lt * cz;
    return { d, T: tr.T, P: tr.P, alpha, aAx, fp: [fpx, fpy, fpz], lift };
  }
  const add = (a, d, k) => { const o = {}; for (const key in a) o[key] = a[key] + d[key] * k; return o; };
  function tidy(s_) {
    // keep the body axis a unit vector, and spin and slosh perpendicular to it
    const bl = Math.hypot(s_.bx, s_.by, s_.bz); s_.bx /= bl; s_.by /= bl; s_.bz /= bl;
    const pr = (kx, ky, kz) => { const q = s_[kx] * s_.bx + s_[ky] * s_.by + s_[kz] * s_.bz; s_[kx] -= q * s_.bx; s_[ky] -= q * s_.by; s_[kz] -= q * s_.bz; };
    pr('wx', 'wy', 'wz'); pr('sx', 'sy', 'sz'); pr('ux', 'uy', 'uz');
    const sl = Math.hypot(s_.sx, s_.sy, s_.sz);
    if (sl > sMax) {
      const nx = s_.sx / sl, ny = s_.sy / sl, nz = s_.sz / sl;
      s_.sx = nx * sMax; s_.sy = ny * sMax; s_.sz = nz * sMax;
      const vr = s_.ux * nx + s_.uy * ny + s_.uz * nz;
      if (vr > 0) { s_.ux -= 1.3 * vr * nx; s_.uy -= 1.3 * vr * ny; s_.uz -= 1.3 * vr * nz; }
    }
  }

  const rec = [], ev = [];
  ev.push({ t: 0, k: 'start', txt: `Fırlatma: ${nf(p.P, 1)} bar, ${nf(Vw0 * 1e6, 0)} ml itki suyu` });
  let t = 0, maxY = 0, tApo = 0, maxV = 0, tMaxV = 0, maxA = 0, maxAlphaUp = 0, impulse = 0, peakT = 0;
  let tWaterEnd = null, tAirEnd = null, vRail = null, tRail = null, tChute = null, lastRec = -1, landed = false, xLand = 0, zLand = 0, tLand = 0, guard = 0;
  let waterLeft = 0, sloshSum = 0, sloshN = 0, xApo = 0, zApo = 0, vBurn = 0, maxDecel = 0, vLand = 0, maxLift = 0;
  const pushRec = (k1, last) => {
    thp += wrapA(Math.atan2(st.by, st.bx) - thp);
    const nx = -Math.sin(thp), ny = Math.cos(thp);
    rec.push({ t, x: st.x, y: last ? Math.max(0, st.y) : st.y, z: st.z, th: thp, bx: st.bx, by: st.by, bz: st.bz,
      s: (st.sx * nx + st.sy * ny) / sMax, sa: Math.hypot(st.sx, st.sy, st.sz) / sMax, sx: st.sx / sMax, sy: st.sy / sMax, sz: st.sz / sMax,
      ph: last ? 3 : phase, T: last ? 0 : k1.T, v: Math.hypot(st.vx, st.vy, st.vz), vx: st.vx, vy: st.vy, vz: st.vz, chute,
      P: last ? 0 : Math.max(0, (k1.P - C.Patm) / 1e5), Vw: st.Vw * 1e6, m: mDry + mPay + st.Vw * C.rhoW + (phase === 2 ? st.ma : phase <= 1 ? maAir0 : 0), ax: last ? 0 : k1.aAx / C.g, al: last ? 0 : k1.alpha * 180 / Math.PI,
      gSx: gS[0], gSy: gS[1], gSz: gS[2], gCx: gC[0], gCy: gC[1], gCz: gC[2],
      gm: last ? 1 : Math.hypot(k1.aAx, k1.fp[0], k1.fp[1], k1.fp[2]) / C.g, wd: chute ? 0 : wingFrac(t), Lw: last ? 0 : k1.lift });
  };
  while (!landed && t < 120 && guard < 600000) {
    guard++;
    const stalled = onRail && t > 0.01 && Math.hypot(st.vx, st.vy, st.vz) < 1e-6 && Math.hypot(st.x, st.y, st.z) < 1e-6;
    const dt = stalled ? 0.002 : phase <= 2 ? 0.0005 : 0.01;
    const k1 = deriv(st, t, phase);
    const k2 = deriv(add(st, k1.d, dt / 2), t + dt / 2, phase);
    const k3 = deriv(add(st, k2.d, dt / 2), t + dt / 2, phase);
    const k4 = deriv(add(st, k3.d, dt), t + dt, phase);
    const nst = {};
    for (const key in st) nst[key] = st[key] + dt / 6 * (k1.d[key] + 2 * k2.d[key] + 2 * k3.d[key] + k4.d[key]);
    impulse += k1.T * dt; peakT = Math.max(peakT, k1.T);
    const prevPhase = phase;
    st = nst; t += dt;
    tidy(st);
    if (onRail && st.x * st.bx + st.y * st.by + st.z * st.bz < 0) { st.x = 0; st.y = 0; st.z = 0; st.vx = 0; st.vy = 0; st.vz = 0; }
    {
      const bv = [st.bx, st.by, st.bz], f = [k1.aAx * bv[0] + k1.fp[0], k1.aAx * bv[1] + k1.fp[1], k1.aAx * bv[2] + k1.fp[2]];
      const gm = Math.hypot(f[0], f[1], f[2]) / C.g, rate = 1 - Math.exp(-dt * (2 + 6 * Math.sqrt(gm)));
      if (gm > 1e-4) gC = v3slerp(gC, v3n([-f[0], -f[1], -f[2]]), rate);
      const sg = k1.aAx > 0 ? -1 : 1;
      gS = v3slerp(gS, v3n([sg * bv[0] + kappa * st.sx, sg * bv[1] + kappa * st.sy, sg * bv[2] + kappa * st.sz]), rate);
    }
    if (phase === 0 && st.xt >= tubeLen) {
      phase = 1;
      const Va = Va0 + Atube * tubeLen;
      Pstart = P0 * Math.pow(Va0 / Va, C.gam); VaStart = Vc - st.Vw;
      ev.push({ t, k: 'tube', txt: `Fırlatma tüpünden çıktı · ${nf(Math.hypot(st.vx, st.vy, st.vz), 1)} m/s` });
    }
    if (phase === 1) {
      const Va = Vc - st.Vw, P = Pstart * Math.pow(VaStart / Va, C.gam);
      if (st.Vw <= 0 || P <= C.Patm) {
        const left_ = st.Vw > 1e-6 ? st.Vw : 0;       // water still inside when the air has no pressure left
        st.Vw = left_; phase = 2; tWaterEnd = t; waterLeft = left_;
        P1 = Math.max(P, C.Patm);
        Tair1 = C.T0 * Math.pow(P1 / P0, (C.gam - 1) / C.gam);
        ma1 = P1 * (Vc - left_) / (C.R * Tair1); st.ma = ma1;
        if (left_ > 0) ev.push({ t, k: 'water', txt: `Basınç bitti · ${nf(left_ * 1e6, 0)} ml su haznede kaldı` });
        else ev.push({ t, k: 'water', txt: `Su bitti, hava evresi · hazne ${nf((P1 - C.Patm) / 1e5, 1)} bar` });
        if (P1 <= C.Patm * AIR_END) { phase = 3; tAirEnd = t; }
      }
    } else if (phase === 2) {
      const { P } = airState(Math.max(st.ma, 1e-9));
      if (P <= C.Patm * AIR_END || st.ma <= 0) { phase = 3; tAirEnd = t; st.ma = 0; }
    }
    if (prevPhase <= 2 && phase === 3) { vBurn = Math.hypot(st.vx, st.vy, st.vz); ev.push({ t, k: 'burnout', txt: `İtki bitti · toplam ${nf(impulse, 1)} N·s` }); }
    if (onRail && Math.hypot(st.x, st.y, st.z) >= railLen) { onRail = false; vRail = Math.hypot(st.vx, st.vy, st.vz); tRail = t; ev.push({ t, k: 'rail', txt: `Rampadan ayrıldı · ${nf(vRail, 1)} m/s` }); }
    const V = Math.hypot(st.vx, st.vy, st.vz);
    if (V > maxV) { maxV = V; tMaxV = t; }
    maxA = Math.max(maxA, Math.abs(k1.aAx)); if (phase === 3 && !onRail) maxDecel = Math.max(maxDecel, -k1.aAx); maxLift = Math.max(maxLift, k1.lift);
    if (!chute && !onRail && st.vy > 0) maxAlphaUp = Math.max(maxAlphaUp, Math.abs(k1.alpha));
    if (!chute) { sloshSum += Math.hypot(st.sx, st.sy, st.sz) / sMax * dt; sloshN += dt; }
    if (st.y > maxY) { maxY = st.y; tApo = t; xApo = st.x; zApo = st.z; }
    if (G.wing === 'acilir' && tWing == null && !chute && !onRail && phase === 3 && st.vy < 0 && t > tApo + 0.15) {
      tWing = t; ev.push({ t, k: 'wing', txt: 'Kanatlar açıldı, süzülme başladı' });
    }
    if (p.chute && !chute && !onRail && phase === 3 && st.vy < 0 && t > tApo + 0.4) {
      chute = true; chuteA = 2 * (mDry + mPay) * C.g / (C.rhoA * 1.5 * 25); tChute = t;
      ev.push({ t, k: 'chute', txt: 'Paraşüt açıldı' });
    }
    if (record && (lastRec < 0 || t - lastRec >= (phase <= 2 ? 0.002 : 0.02))) { lastRec = t; pushRec(k1, false); }
    if (t > 0.05 && st.y < 0 && !onRail) { landed = true; xLand = st.x; zLand = st.z; tLand = t; vLand = Math.hypot(st.vx, st.vy, st.vz); }
    if (onRail && t > 3) break;
  }
  const left = !onRail;
  let truncated = false;
  if (left && !landed) {   // still flying when the 120 s limit was reached
    truncated = true; xLand = st.x; zLand = st.z; tLand = t; vLand = Math.hypot(st.vx, st.vy, st.vz);
    ev.push({ t, k: 'cut', txt: `Hesap ${nf(t, 0)} s'de durduruldu · roket hâlâ ${nf(Math.max(0, st.y), 1)} m yüksekte` });
  }
  if (phase <= 2) {       // landed (or stopped) while the bottle was still pushing
    if (tWaterEnd == null) tWaterEnd = t;
    tAirEnd = t; vBurn = Math.hypot(st.vx, st.vy, st.vz);
  }
  if (record && rec.length) pushRec(null, true);
  const dist = Math.hypot(xLand, zLand);
  if (left) {
    ev.push({ t: tMaxV, k: 'maxv', txt: `Azami hız · ${nf(maxV, 1)} m/s` });
    ev.push({ t: tApo, k: 'apo', txt: `Tepe noktası · ${nf(maxY, 1)} m` });
    if (!truncated) ev.push({ t: tLand, k: 'land', txt: `Yere indi · ${nf(dist, 1)} m uzakta` + (Math.abs(zLand) >= 0.3 ? ` (yana ${nf(Math.abs(zLand), 1)} m)` : '') });
  } else {
    ev.push({ t, k: 'fail', txt: 'Roket rampadan ayrılamadı' });
  }
  ev.sort((a, c) => a.t - c.t);
  const ratio = mPay / G.TOW;
  const root = Math.sqrt(dist * dist + maxY * maxY);
  const score = left ? (root + tLand) * ratio : 0;
  return { p, G, rec, ev, score, root, dist, apo: maxY, tApo, xApo, zApo, xLand, zLand, time: left ? tLand : t, maxV, maxA, impulse, peakT,
    tWaterEnd, tAirEnd, vRail, tRail, tChute, maxAlphaUp, slosh: sloshN ? sloshSum / sloshN : 0, ratio, left, zeta, mS,
    kappa, sMax, th0, P1g: Math.max(0.01, (P1 - C.Patm) / 1e5), vBurn, maxDecel, vLand, tWing, maxLift, maAir0, waterLeft, Pstart, truncated };
}

// ================================================================== helpers
const $ = (s) => document.querySelector(s);
const $$ = (s) => Array.from(document.querySelectorAll(s));
function nf(v, d = 1) { let n = Number(v); if (Math.abs(n) < 0.5 * Math.pow(10, -d)) n = 0; return n.toLocaleString('tr-TR', { minimumFractionDigits: d, maximumFractionDigits: d }); }
const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
const store = {
  get(k, d) { try { const v = localStorage.getItem(k); return v === null ? d : JSON.parse(v); } catch (e) { return d; } },
  set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { /* ignore */ } },
};

// control values are kept in the units shown next to each slider (L, kg, %, bar, mm, cm, °)
const DEFAULTS = { Vc: 3, payload: 1, dry: 0.7, Dmm: 100, CD: 0.42, fins: 'orta', baffle: 40, P: 8, Vsu: 0.9, nozzle: 22, cd: 0.96, tube: false,
  wing: 'yok', wspan: 1.2, wchord: 18, winc: 6, angle: 85, wind: 2, windZ: 0, chute: false };
const PRESETS = {
  ornek: {},
  bolmesiz: { baffle: 0 },
  agir: { Vc: 5, payload: 2, dry: 0.9, Vsu: 2.0 },
  hafif: { dry: 0.4, fins: 'kucuk' },
  kuraldisi: { Vsu: 1.35, baffle: 60, P: 11 },
  yanruzgar: { windZ: 3.5 },
  suzulen: { dry: 0.4, payload: 0.5, Vc: 2, Vsu: 0.5, P: 10, tube: true, wing: 'acilir', wspan: 1.6, wchord: 20, winc: 6 },
};
// number boxes: stored value × scale = shown value; min/max are for what can be typed
const NUMS = {
  Vc: { scale: 1, dec: 2, min: 0.3, max: 10 }, payload: { scale: 1, dec: 2, min: 0.1, max: 5 }, dry: { scale: 1000, dec: 0, min: 0.05, max: 3 },
  Dmm: { scale: 1, dec: 0, min: 40, max: 200 }, CD: { scale: 1, dec: 2, min: 0.05, max: 1.5 }, baffle: { scale: 1, dec: 1, min: 0, max: 95 },
  P: { scale: 1, dec: 2, min: 0.5, max: 15 }, Vsu: { scale: 1, dec: 2, min: 0.01, max: 8 }, nozzle: { scale: 1, dec: 1, min: 3, max: 40 },
  cd: { scale: 1, dec: 2, min: 0.5, max: 1 }, wspan: { scale: 1, dec: 2, min: 0.2, max: 3 }, wchord: { scale: 1, dec: 0, min: 3, max: 60 },
  winc: { scale: 1, dec: 1, min: -5, max: 15 }, angle: { scale: 1, dec: 1, min: 60, max: 90 }, wind: { scale: 1, dec: 1, min: -15, max: 15 },
  windZ: { scale: 1, dec: 1, min: -15, max: 15 },
};
const VAL = Object.assign({}, DEFAULTS);
const radioVal = (name) => { const el = document.querySelector(`input[name=${name}]:checked`); return el ? el.value : DEFAULTS[name]; };
function paintRange(el) { el.style.setProperty('--p', clampN((el.value - el.min) / (el.max - el.min) * 100, 0, 100).toFixed(1) + '%'); }
function clampN(v, a, b) { return Math.max(a, Math.min(b, v)); }
const fmtNum = (k, v) => { const c = NUMS[k]; return String(+(v * c.scale).toFixed(c.dec)).replace('.', ','); };
function showNum(k, force) {
  const el = $('#n-' + k); if (!el || (!force && document.activeElement === el)) return;
  el.value = fmtNum(k, VAL[k]);
}
function setVal(k, v, from) {
  const c = NUMS[k]; VAL[k] = clampN(v, c.min, c.max);
  const r = $('#c-' + k); if (r && from !== 'range') { r.value = VAL[k]; paintRange(r); }
  if (from !== 'num') showNum(k);
}
// physics parameters (SI fractions) from the control values
function fromCtl(v) {
  const Vc = v.Vc, Vsu = Math.min(v.Vsu, Vc * 0.95);
  return { Vc, payload: v.payload, dry: v.dry, fins: v.fins, baffle: v.baffle / 100, P: v.P, wf: Vsu / Vc, nozzle: v.nozzle,
    tube: !!v.tube, angle: v.angle, wind: v.wind, windZ: v.windZ || 0, chute: !!v.chute, cd: v.cd, Dmm: v.Dmm, CD: v.CD,
    wing: v.wing || 'yok', wspan: v.wspan, wchord: v.wchord / 100, winc: v.winc };
}
function toCtl(p) {
  return { Vc: p.Vc, payload: p.payload, dry: p.dry, fins: p.fins, baffle: Math.round(p.baffle * 1000) / 10, P: p.P, Vsu: +(p.Vc * p.wf).toFixed(3), nozzle: p.nozzle,
    tube: p.tube, angle: p.angle, wind: p.wind, windZ: p.windZ, chute: p.chute, cd: p.cd, Dmm: p.Dmm, CD: p.CD,
    wing: p.wing, wspan: p.wspan, wchord: Math.round(p.wchord * 1000) / 10, winc: p.winc };
}
function readParams() {
  VAL.fins = radioVal('fins'); VAL.wing = radioVal('wing');
  VAL.tube = $('#c-tube').checked; VAL.chute = $('#c-chute').checked;
  return fromCtl(VAL);
}
function writeControls(v) {
  Object.keys(NUMS).forEach((k) => setVal(k, v[k] == null ? DEFAULTS[k] : v[k]));
  ['fins', 'wing'].forEach((n) => { const el = document.querySelector(`input[name=${n}][value="${v[n] || DEFAULTS[n]}"]`); if (el) el.checked = true; });
  $('#c-tube').checked = !!v.tube; $('#c-chute').checked = !!v.chute;
}

function checkRules(p, r) {
  const G = r.G;
  return [
    { label: 'Kalkış ağırlığı ≤ 5 kg', val: nf(G.TOW, 2) + ' kg', ok: G.TOW <= 5.0001, ref: 'TY 2.1.8' },
    { label: 'Roket boyu ≤ 1,5 m', val: nf(G.L, 2) + ' m', ok: G.L <= 1.5, ref: 'TY 2.1.1' },
    { label: 'Basınç ≤ 10 bar', val: nf(p.P, 2) + ' bar', ok: p.P <= 10 + 1e-9, ref: 'Tİ 3 · TY 2.2.1' },
    { label: 'İtki suyu ≤ çalkalanan su', val: nf(G.mProp, 2) + ' / ' + nf(G.mPay, 2) + ' kg', ok: G.mProp <= G.mPay + 1e-6, ref: 'TY 2.1.7' },
    { label: 'Çalkalanma tankında ≥ 500 ml su', val: nf(p.payload * 1000, 0) + ' ml', ok: p.payload >= 0.5, ref: 'TY 2.1.5' },
    { label: 'Bölme levhası ≤ kesitin %50\'si', val: '%' + nf(p.baffle * 100, 1), ok: p.baffle <= 0.5 + 1e-9, ref: 'TY 2.3.1' },
    { label: 'Fırlatma açısı 85° ± 3°', val: nf(p.angle, 1) + '°', ok: p.angle >= 82 && p.angle <= 88, ref: 'TY 2.1.3' },
    { label: 'Kanat açıklığı ≤ 2,25 m', val: G.wing ? nf(p.wspan, 2) + ' m' : 'kanat yok', ok: !G.wing || p.wspan <= 2.25 + 1e-9, ref: 'TY 2.1.2' },
  ];
}
const rulesOk = (p, r) => checkRules(p, r).every((x) => x.ok);
const allOk = (p, r) => rulesOk(p, r) && r.left;

// ================================================================== state
const S = {
  view: 'konu', topic: 'genel', tab: 'ucus',
  cur: null, ref: null, refName: 'İKARUS-1 (örnek)',
  t: 0, tEnd: 0, playing: false, lastFrame: 0, dim: '2',
};

// ================================================================== guide
const topics = $$('#topics .topic');
const topicIds = topics.map((a) => a.id.slice(2));
function buildToc() {
  const groups = {};
  topics.forEach((a) => { (groups[a.dataset.group] = groups[a.dataset.group] || []).push(a); });
  const toc = $('#toc'), sel = $('#toc-sel');
  toc.innerHTML = ''; sel.innerHTML = '';
  Object.entries(groups).forEach(([g, list]) => {
    const div = document.createElement('div'); div.className = 'toc-g';
    div.innerHTML = `<span class="eyebrow">${g}</span>`;
    const og = document.createElement('optgroup'); og.label = g;
    list.forEach((a) => {
      const id = a.id.slice(2);
      const b = document.createElement('button'); b.type = 'button'; b.dataset.topic = id; b.textContent = a.dataset.title;
      b.addEventListener('click', () => showTopic(id, true));
      div.appendChild(b);
      const o = document.createElement('option'); o.value = id; o.textContent = a.dataset.title; og.appendChild(o);
    });
    toc.appendChild(div); sel.appendChild(og);
  });
  sel.addEventListener('change', () => showTopic(sel.value, true));
  // prev/next
  topics.forEach((a, i) => {
    const nav = document.createElement('footer'); nav.className = 't-nav';
    const prev = topics[i - 1], next = topics[i + 1];
    if (prev) nav.innerHTML += `<button type="button" class="prev" data-goto="${prev.id.slice(2)}"><span>← Önceki</span><b>${prev.dataset.title}</b></button>`;
    if (next) nav.innerHTML += `<button type="button" class="next" data-goto="${next.id.slice(2)}"><span>Sonraki →</span><b>${next.dataset.title}</b></button>`;
    else nav.innerHTML += `<button type="button" class="next" data-gosim="ucus"><span>Sonraki →</span><b>Simülasyona geç</b></button>`;
    a.appendChild(nav);
  });
}
function showTopic(id, scroll) {
  if (!topicIds.includes(id)) id = 'genel';
  S.topic = id;
  topics.forEach((a) => { a.hidden = a.id !== 't-' + id; });
  $$('#toc button').forEach((b) => b.setAttribute('aria-current', b.dataset.topic === id ? 'true' : 'false'));
  $('#toc-sel').value = id;
  const tb = $('#toc button[aria-current="true"]'), tc = $('#toc');
  if (tb && tc && tc.scrollHeight > tc.clientHeight + 4) {
    const top = tb.offsetTop;
    if (top < tc.scrollTop + 10 || top > tc.scrollTop + tc.clientHeight - 40) tc.scrollTop = Math.max(0, top - tc.clientHeight / 2);
  }
  if (S.view !== 'konu') showView('konu');
  setHash(id);
  if (scroll) window.scrollTo({ top: 0, behavior: 'auto' });
}

// ================================================================== views
function setHash(h) { try { history.replaceState(null, '', '#' + h); } catch (e) { /* ignore */ } }
function showView(v) {
  S.view = v;
  $('#view-konu').hidden = v !== 'konu';
  $('#view-sim').hidden = v !== 'sim';
  $('#nav-konu').setAttribute('aria-selected', v === 'konu' ? 'true' : 'false');
  $('#nav-sim').setAttribute('aria-selected', v === 'sim' ? 'true' : 'false');
  if (v === 'sim') { setHash('sim-' + S.tab); requestAnimationFrame(() => renderTab(true)); }
  else { setHash(S.topic); pause(); }
}
function showTab(tab) {
  S.tab = tab;
  $$('.subtabs button').forEach((b) => b.setAttribute('aria-selected', b.dataset.tab === tab ? 'true' : 'false'));
  $$('[data-panel]').forEach((p) => { p.hidden = p.dataset.panel !== tab; });
  if (S.view !== 'sim') showView('sim'); else setHash('sim-' + tab);
  requestAnimationFrame(() => renderTab(true));
}

// ================================================================== compute & render
function outputs(p, G) {
  Object.keys(NUMS).forEach(showNum);
  const H = (id, txt, bad) => { const el = $('#' + id); if (!el) return; el.textContent = txt; el.className = 'hint' + (bad ? ' bad' : ''); };
  H('h-payload', `Tank hacmi ${nf(p.payload * 2, 1)} L (%50 su, %50 hava) · ${nf(p.payload, 2)} kg.`);
  H('h-dry', `Gövde, burun, altimetre, tesisat. Kanatçıklar ${nf(G.fin.mass * 1000, 0)} g${G.wing ? `, kanat ve kuyruk ${nf(G.mWing * 1000, 0)} g` : ''} ayrıca eklenir → toplam boş ${nf(G.mDry * 1000, 0)} g.`);
  H('h-D', `Kesit alanı A_b = ${nf(G.Ab * 1e4, 1)} cm² · roket boyu ${nf(G.L, 2)} m.`, G.L > 1.5);
  H('h-baffle', p.baffle > 0.5 ? 'Kural dışı: tek levha kesitin en fazla %50\'sini kapatabilir.' : 'Kural: tek bir levha tank kesitinin en fazla %50\'sini kapatabilir.', p.baffle > 0.5);
  H('h-P', `Mutlak P₀ = P + 1 atm = ${nf(p.P + 1.01325, 2)} bar (${nf((p.P + 1.01325) * 1e5, 0)} Pa) · ${nf(p.P * 14.5038, 0)} psi gösterge.` + (p.P > 10 ? ' Kural dışı: en fazla 10 bar.' : p.P > 9.5 ? ' Sınıra çok yakın; manometre hatası için pay bırakın.' : ''), p.P > 10);
  const okW = G.mProp <= G.mPay + 1e-6, Vsu = p.Vc * p.wf, clipped = VAL.Vsu > p.Vc * 0.95 + 1e-9;
  H('h-wf', clipped ? `Hazneden fazla su yazılamaz; ${nf(Vsu, 2)} L ile sınırlandı.` : `V_hava0 = V_tank − V_su = ${nf(p.Vc - Vsu, 2)} L · haznenin %${nf(p.wf * 100, 0)}'i · ${okW ? 'çalkalanan sudan fazla değil' : 'çalkalanan sudan fazla (kural dışı)'}.`, !okW || clipped);
  H('h-noz', `Nozul alanı A_n = π·(D/2)² = ${nf(G.An * 1e6, 0)} mm².`);
  const g = $('#g-wing'); if (g) { g.classList.toggle('off', !G.wing); g.querySelectorAll('.wing-params input').forEach((el) => { el.disabled = !G.wing; }); }
  H('h-wspan', G.wing && p.wspan > 2.25 ? 'Kural dışı: kanat açıklığı en fazla 2,25 m.' : 'Kural: kanat açıklığı en fazla 2,25 m.', G.wing && p.wspan > 2.25);
  H('h-wing', G.wing ? `Kanat alanı S = b·c = ${nf(G.Sw, 3)} m², en-boy oranı ${nf(G.AR, 1)}. Modelde kuyruğa ${nf(G.St * 1e4, 0)} cm²'lik bir yatay dengeleyici eklendiği ve kanadı dengelediği varsayılır.` + (p.chute ? ' Paraşüt açılınca modelde kanat devre dışı kalır.' : '') : 'Kanat tipini Sabit veya Açılır yapınca bu ayarlar açılır.');
  $$('#controls input[type=range]').forEach(paintRange);
}

let heavyTimer = null;
function recompute(opts = {}) {
  const p = readParams();
  S.cur = simulate(p, true);
  S.cur.emit = buildEmit(S.cur);
  S.tEnd = S.cur.left ? S.cur.time + 1.4 : S.cur.time;
  outputs(p, S.cur.G);
  renderToolbar();
  pause();
  S.t = opts.t !== undefined ? opts.t : S.tEnd;
  setScrubMax(); setPlayUI();
  dirty.charts = true; dirty.sweep = true; dirty.tips = true; dirty.roket = true;
  renderTab(false);
}

const dirty = { charts: true, sweep: true, tips: true, roket: true };
function renderTab(shown) {
  if (S.view !== 'sim' || !S.cur) return;
  if (S.tab === 'ucus') { renderEvents(); if (S.dim === '3') size3(); else sizeCanvases(); drawFrame(); }
  else if (S.tab === 'puan') { renderScore(); clearTimeout(heavyTimer); heavyTimer = setTimeout(renderTips, shown ? 0 : 120); }
  else if (S.tab === 'grafik') { if (dirty.charts || shown) renderCharts(); }
  else if (S.tab === 'rapor') renderReport();
  else if (S.tab === 'roket') { drawAnatomy(S.cur.p, S.cur.G); renderMass(); }
  else if (S.tab === 'tarama') { clearTimeout(heavyTimer); heavyTimer = setTimeout(renderSweep, shown ? 0 : 150); }
}

function deltaTo(d, a, b) {
  const diff = a - b;
  if (b <= 0) { d.className = 'delta ' + (a > 0 ? 'up' : 'same'); d.textContent = (a > 0 ? '▲ +' + nf(a, 2) : '≈ 0,00') + ' · ' + S.refName + ' (uçmuyor)'; return; }
  const pct = diff / b * 100;
  d.className = 'delta ' + (Math.abs(pct) < 0.5 ? 'same' : diff > 0 ? 'up' : 'down');
  d.textContent = (Math.abs(pct) < 0.5 ? '≈ ' : diff > 0 ? '▲ ' : '▼ ') + (diff >= 0 ? '+' : '') + nf(diff, 2) + ' · ' + S.refName;
}
function renderToolbar() {
  const r = S.cur, ok = allOk(r.p, r);
  $('#tb-score').textContent = nf(r.score, 2);
  const v = $('#tb-verdict');
  const bad = checkRules(r.p, r).filter((x) => !x.ok).length;
  v.className = 'chip ' + (ok ? 'ok' : 'no');
  v.textContent = ok ? '✓ kurallara uygun' : r.left ? `✕ ${bad} kural ihlali` : '✕ rampadan ayrılamadı';
  const d = $('#tb-delta');
  deltaTo(d, r.score, S.ref.score);
}

// ------------------------------------------------------------------ playback
function setScrubMax() {
  const sc = $('#scrub');
  sc.max = String(Math.max(1, Math.round(S.tEnd * 1000)));
  sc.value = String(Math.round(S.t * 1000));
  sc.style.setProperty('--p', (S.t / Math.max(1e-6, S.tEnd) * 100).toFixed(2) + '%');
}
function powerEnd() { return (S.cur.tAirEnd || 0) + 0.05; }
function setT(t, fromScrub) {
  S.t = Math.max(0, Math.min(t, S.tEnd));
  const sc = $('#scrub');
  if (!fromScrub) sc.value = String(Math.round(S.t * 1000));
  sc.style.setProperty('--p', (S.t / Math.max(1e-6, S.tEnd) * 100).toFixed(2) + '%');
  if (S.tab === 'ucus') drawFrame();
  else if (S.tab === 'grafik') redrawCursor();
}
function setPlayUI() {
  const atEnd = S.t >= S.tEnd - 1e-6;
  $('#play-t').textContent = S.playing ? 'Duraklat' : (S.t > 0 && !atEnd ? 'Devam et' : 'Fırlat');
  const d = S.playing ? 'M5 4h3.5v12H5zM11.5 4H15v12h-3.5z' : 'M6 4l10 6-10 6z';
  $('#play-ico').setAttribute('d', d);
  $('#pb-ico').setAttribute('d', d);
  $('#rec').classList.toggle('on', S.playing);
  $('#rec3').classList.toggle('on', S.playing);
}
function play() {
  if (!S.cur || !S.cur.rec.length) return;
  if (S.t >= S.tEnd - 1e-6) S.t = 0;
  if (S.tab !== 'ucus' && S.tab !== 'grafik') showTab('ucus');
  S.playing = true; S.lastFrame = performance.now(); setPlayUI();
  requestAnimationFrame(frame);
}
function pause() { S.playing = false; if (S.cur) setPlayUI(); }
function slowFactor(t) {
  if (!$('#slowmo').checked || !S.cur.left) return 1;
  const pe = powerEnd();
  if (t < pe) return 0.06;
  const u = (t - pe) / 0.4;
  return u >= 1 ? 1 : 0.06 + 0.94 * u * u;
}
function frame(now) {
  if (!S.playing) return;
  const dtReal = Math.min(0.05, (now - S.lastFrame) / 1000); S.lastFrame = now;
  const spd = +$('#speed').value;
  let nt = S.t + dtReal * spd * slowFactor(S.t);
  if (nt >= S.tEnd) { nt = S.tEnd; S.playing = false; }
  setT(nt);
  setPlayUI();
  if (S.playing) requestAnimationFrame(frame);
}

// ------------------------------------------------------------------ flight view: helpers
const PAD_H = 0.1, NOZ_L = 0.025, EX_CLOSE = 1.6;
const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
const lerp = (a, b, f) => a + (b - a) * f;
const smooth = (u) => { u = clamp(u, 0, 1); return u * u * (3 - 2 * u); };
const hash = (n) => { const s = Math.sin(n * 127.1 + 311.7) * 43758.5453123; return s - Math.floor(s); };
function mixHex(a, b, f) {
  const pa = parseInt(a.slice(1), 16), pb = parseInt(b.slice(1), 16);
  const c = (s) => Math.round(((pa >> s) & 255) * (1 - f) + ((pb >> s) & 255) * f);
  return `rgb(${c(16)},${c(8)},${c(0)})`;
}
function rrect(ctx, x, y, w, h, r) {
  if (w < 0) { x += w; w = -w; }
  r = Math.max(0, Math.min(r, w / 2, h / 2));
  ctx.beginPath(); ctx.moveTo(x + r, y);
  ctx.arcTo(x + w, y, x + w, y + h, r); ctx.arcTo(x + w, y + h, x, y + h, r);
  ctx.arcTo(x, y + h, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath();
}
let FONTS = null;
const F = () => FONTS || (FONTS = { d: css('--f-display'), m: css('--f-mono'), b: css('--f-body') });
function pill(ctx, x, y, text, o = {}) {
  // x,y: anchor (left/right/center at the vertical middle)
  ctx.font = o.font || `500 11px ${F().m}`;
  const padX = o.padX || 7, h = o.h || 20, w = ctx.measureText(text).width + padX * 2;
  const x0 = o.align === 'right' ? x - w : o.align === 'center' ? x - w / 2 : x;
  rrect(ctx, x0, y - h / 2, w, h, o.r || 6);
  ctx.fillStyle = o.bg || 'rgba(4, 12, 22, 0.62)'; ctx.fill();
  if (o.stroke) { ctx.strokeStyle = o.stroke; ctx.lineWidth = 1; ctx.stroke(); }
  ctx.fillStyle = o.color || '#eaf3fa'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
  ctx.fillText(text, x0 + padX, y + 0.5);
  return w;
}

const scene = $('#scene'), sctx = scene.getContext('2d');
const cu = $('#closeup'), cctx = cu.getContext('2d');
const CAM = { mode: store.get('ikarus-cam', 'takip') === 'genel' ? 'genel' : 'takip' };

function setCanvas(cv, ctx, w, h) {
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  const W = Math.round(w * dpr), H = Math.round(h * dpr);
  if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
  cv.style.height = h + 'px';
  cv._w = w; cv._h = h;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
}
function sizeCanvases() {
  const w = scene.parentElement.clientWidth;
  if (!w) return;
  const h = Math.round(clamp(w * 0.66, 300, 500));
  setCanvas(scene, sctx, w, h);
  const cw = cu.parentElement.clientWidth;
  if (!cw) return;
  const side = getComputedStyle($('.feed-grid')).gridTemplateColumns.trim().split(/\s+/).length > 1;
  const ch = side ? Math.max(260, h - $('#cap').offsetHeight) : Math.round(clamp(cw * 0.95, 340, 480));
  setCanvas(cu, cctx, cw, ch);
}

const SPR = {};
function makeSprites() {
  const mk = (w, h, f) => { const c = document.createElement('canvas'); c.width = w; c.height = h; f(c.getContext('2d')); return c; };
  SPR.puff = mk(64, 64, (g) => {
    const gr = g.createRadialGradient(32, 32, 0, 32, 32, 32);
    gr.addColorStop(0, 'rgba(255,255,255,0.95)'); gr.addColorStop(0.45, 'rgba(240,247,255,0.55)'); gr.addColorStop(1, 'rgba(240,247,255,0)');
    g.fillStyle = gr; g.fillRect(0, 0, 64, 64);
  });
  SPR.clouds = [0, 1, 2, 3].map((v) => mk(300, 130, (g) => {
    for (let i = 0; i < 16; i++) {
      const x = 40 + 220 * hash(v * 50 + i), y = 78 + 22 * (hash(v * 31 + i * 3) - 0.5) - 24 * Math.sin(Math.PI * (x - 40) / 220);
      const rr = 20 + 26 * hash(v * 17 + i * 7) * Math.sin(Math.PI * (x - 30) / 240);
      const gr = g.createRadialGradient(x, y - rr * 0.2, 0, x, y, rr);
      gr.addColorStop(0, 'rgba(255,255,255,0.95)'); gr.addColorStop(0.55, 'rgba(250,252,255,0.6)'); gr.addColorStop(1, 'rgba(235,243,252,0)');
      g.fillStyle = gr; g.beginPath(); g.arc(x, y, rr, 0, Math.PI * 2); g.fill();
    }
    const sh = g.createLinearGradient(0, 70, 0, 120);
    sh.addColorStop(0, 'rgba(160,180,205,0)'); sh.addColorStop(1, 'rgba(160,180,205,0.35)');
    g.globalCompositeOperation = 'source-atop'; g.fillStyle = sh; g.fillRect(0, 0, 300, 130);
  }));
}

// ------------------------------------------------------------------ flight state at time t
const REC_KEYS = ['x', 'y', 'z', 'th', 'bx', 'by', 'bz', 's', 'sa', 'sx', 'sy', 'sz', 'v', 'T', 'P', 'Vw', 'ax', 'al', 'vx', 'vy', 'vz',
  'gSx', 'gSy', 'gSz', 'gCx', 'gCy', 'gCz', 'gm', 'wd', 'Lw', 'm'];
function recAt(r, t) {
  const a = r.rec; if (!a.length) return null;
  if (t >= a[a.length - 1].t) return a[a.length - 1];
  if (t <= a[0].t) return Object.assign({}, a[0], { t });
  let lo = 0, hi = a.length - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (a[mid].t <= t) lo = mid; else hi = mid; }
  const A = a[lo], B = a[hi], f = (t - A.t) / Math.max(1e-9, B.t - A.t);
  const o = { t, ph: A.ph, chute: A.chute };
  REC_KEYS.forEach((k) => { o[k] = A[k] + (B[k] - A[k]) * f; });
  return o;
}
function padRec(r) {
  const th0 = r.th0, s0 = clamp(-(Math.cos(th0) / (r.kappa * Math.sin(th0))) / r.sMax, -1, 1);
  return { t: 0, x: 0, y: 0, z: 0, th: th0, bx: Math.cos(th0), by: Math.sin(th0), bz: 0, s: s0, sa: Math.abs(s0),
    sx: -Math.sin(th0) * s0, sy: Math.cos(th0) * s0, sz: 0, ph: r.p.tube ? 0 : 1, T: 0, v: 0, vx: 0, vy: 0, vz: 0,
    chute: false, P: r.p.P, Vw: r.G.Vw0 * 1e6, ax: 0, al: 0, gSx: 0, gSy: -1, gSz: 0, gCx: 0, gCy: -1, gCz: 0, gm: 1, wd: r.G.wing === 'sabit' ? 1 : 0, Lw: 0, m: r.G.TOW + r.maAir0 };
}
function waterFill(r, t, now, pad) {
  const G = r.G;
  const fc = pad ? G.Vw0 / G.Vc : clamp(now.Vw / (G.Vc * 1e6), 0, 1);
  let mist = 0;
  if (!pad && r.tWaterEnd != null && t >= r.tWaterEnd) {
    if (t < (r.tAirEnd || 0)) mist = 0.3 + 0.5 * clamp(now.P / r.P1g, 0, 1);
    else mist = 0.32 * (1 - clamp((t - (r.tAirEnd || 0)) / 0.6, 0, 1));
  }
  return { fc, mist };
}
// side view (x–y plane)
function poseAt(r, t) {
  const G = r.G, L = G.L;
  const pad = t <= 0 || !r.rec.length;
  const now = pad ? padRec(r) : recAt(r, Math.min(t, r.time));
  let th = now.th;
  if (now.chute && r.tChute != null) {
    // hanging nose-up under the canopy
    const thc = recAt(r, r.tChute).th, e = smooth((t - r.tChute) / 0.9);
    th = thc + wrapA(Math.PI / 2 - thc) * e + 0.05 * Math.sin(2.1 * t) * e;
  }
  let bx = Math.cos(th), by = Math.sin(th);
  // with a side wind the body also points out of the view plane: draw it shortened (orthographic view)
  let fs = pad || now.chute ? 1 : clamp(Math.hypot(now.bx, now.by), 0.2, 1), Le = L * fs;
  let tx = now.x, ty = now.y + PAD_H;
  const low = Math.min(ty - NOZ_L * by, ty + Le * by);
  if (low < 0) ty -= low;
  const landed = r.left && !r.truncated && t > r.time;
  let cX = 0, cY = 0;
  if (landed) {
    // tip over after touchdown, pivoting on the contact point
    const ee = Math.pow(clamp((t - r.time) / 0.55, 0, 1), 2);
    const noseLow = by < 0;
    cX = noseLow ? tx + Le * bx : tx - NOZ_L * bx; cY = noseLow ? ty + Le * by : ty - NOZ_L * by;
    fs += (1 - fs) * ee; Le = L * fs;
    const thF = noseLow ? (bx >= 0 ? -0.1 : Math.PI + 0.1) : (bx >= 0 ? 0.06 : Math.PI - 0.06);
    th = th + wrapA(thF - th) * ee;
    bx = Math.cos(th); by = Math.sin(th);
    if (noseLow) { tx = cX - Le * bx; ty = cY - Le * by; } else { tx = cX + NOZ_L * bx; ty = cY + NOZ_L * by; }
    ty += ee * (G.R + 0.02);
  }
  // water directions: world angle of apparent gravity (projected on the view plane) minus the drawn body angle
  let aS = Math.atan2(now.gSy, now.gSx), aC = Math.atan2(now.gCy, now.gCx);
  if (landed) {
    const f = smooth((t - r.time) / 0.5);
    aS += wrapA(-Math.PI / 2 - aS) * f; aC += wrapA(-Math.PI / 2 - aC) * f;
  }
  const wf = waterFill(r, t, now, pad);
  return { now, pad, landed, th, bx, by, fs, x: tx, y: ty, nx: tx + Le * bx, ny: ty + Le * by,
    cx: tx + 0.5 * Le * bx, cy: ty + 0.5 * Le * by, gs: aS - th, gc: aC - th, fc: wf.fc, mist: wf.mist, s: now.s, sa: now.sa, cX, cY };
}
// full 3D pose
function pose3At(r, t) {
  const G = r.G, L = G.L;
  const pad = t <= 0 || !r.rec.length;
  const now = pad ? padRec(r) : recAt(r, Math.min(t, r.time));
  let b = v3n([now.bx, now.by, now.bz]);
  if (now.chute && r.tChute != null) {
    const c = recAt(r, r.tChute), e = smooth((t - r.tChute) / 0.9);
    b = v3slerp(v3n([c.bx, c.by, c.bz]), v3n([0.05 * Math.sin(2.1 * t), 1, 0.05 * Math.cos(1.7 * t)]), e);
  }
  const tail = [now.x, now.y + PAD_H, now.z];
  const low = Math.min(tail[1] - NOZ_L * b[1], tail[1] + L * b[1]);
  if (low < 0) tail[1] -= low;
  const landed = r.left && !r.truncated && t > r.time;
  let contact = null;
  if (landed) {
    const ee = Math.pow(clamp((t - r.time) / 0.55, 0, 1), 2);
    const noseLow = b[1] < 0;
    contact = noseLow ? [tail[0] + L * b[0], tail[1] + L * b[1], tail[2] + L * b[2]] : [tail[0] - NOZ_L * b[0], tail[1] - NOZ_L * b[1], tail[2] - NOZ_L * b[2]];
    let hdir = [b[0], 0, b[2]]; if (Math.hypot(hdir[0], hdir[2]) < 1e-3) hdir = [1, 0, 0]; hdir = v3n(hdir);
    const el = noseLow ? -0.1 : 0.06;
    const bF = [hdir[0] * Math.cos(el), Math.sin(el), hdir[2] * Math.cos(el)];
    b = v3slerp(b, bF, ee);
    if (noseLow) { tail[0] = contact[0] - L * b[0]; tail[1] = contact[1] - L * b[1]; tail[2] = contact[2] - L * b[2]; }
    else { tail[0] = contact[0] + NOZ_L * b[0]; tail[1] = contact[1] + NOZ_L * b[1]; tail[2] = contact[2] + NOZ_L * b[2]; }
    tail[1] += ee * (G.R + 0.02);
  }
  let gS = v3n([now.gSx, now.gSy, now.gSz]), gC = v3n([now.gCx, now.gCy, now.gCz]);
  if (landed) { const f = smooth((t - r.time) / 0.5); gS = v3slerp(gS, [0, -1, 0], f); gC = v3slerp(gC, [0, -1, 0], f); }
  const wf = waterFill(r, t, now, pad);
  const at = (u) => [tail[0] + u * b[0], tail[1] + u * b[1], tail[2] + u * b[2]];
  return { now, pad, landed, b, tail, nose: at(L), center: at(L / 2), at, gS, gC, fc: wf.fc, mist: wf.mist, sa: now.sa, contact };
}
function buildEmit(r) {
  const out = [];
  const tEnd = r.tAirEnd || 0;
  if (!r.rec.length || tEnd <= 0) return out;
  let i = 0;
  for (let tau = 0.0008; tau <= tEnd; tau += 0.0012) {
    const e = recAt(r, tau); i++;
    if (e.ph !== 1 && e.ph !== 2) continue;
    const bx = Math.cos(e.th), by = Math.sin(e.th), b3 = v3n([e.bx, e.by, e.bz]);
    out.push({ tau, i, x: e.x - bx * NOZ_L, y: e.y + PAD_H - by * NOZ_L, vx: e.vx, vy: e.vy, bx, by, ph: e.ph, P: e.P,
      p3: [e.x - b3[0] * NOZ_L, e.y + PAD_H - b3[1] * NOZ_L, e.z - b3[2] * NOZ_L], v3: [e.vx, e.vy, e.vz], b3 });
  }
  return out;
}
function boundsOf(r) {
  if (r._b) return r._b;
  let x0 = 0, x1 = 0, y1 = 0, z0 = 0, z1 = 0;
  r.rec.forEach((e) => { x0 = Math.min(x0, e.x); x1 = Math.max(x1, e.x); y1 = Math.max(y1, e.y); z0 = Math.min(z0, e.z); z1 = Math.max(z1, e.z); });
  return (r._b = { x0, x1, y1, z0, z1 });
}

// ------------------------------------------------------------------ cameras
function camOverview(W, H, r) {
  let xmin = -7.5, xmax = 4, ymax = 6;
  [S.ref, r].forEach((q) => { if (!q) return; const b = boundsOf(q); xmin = Math.min(xmin, b.x0 - 2); xmax = Math.max(xmax, b.x1 + 3); ymax = Math.max(ymax, b.y1 * 1.15 + 1.5); });
  const k = Math.min((W - 30) / (xmax - xmin), 0.78 * H / ymax);
  return { k, ox: W / 2 - (xmin + xmax) / 2 * k, oy: 0.88 * H };
}
function camFollow(W, H, pose) {
  const viewH = 9, k = H / viewH;
  const cyw = Math.max(pose.cy, 0.36 * viewH);
  return { k, ox: W / 2 - pose.cx * k, oy: H / 2 + cyw * k };
}
function mainCam(W, H, r, t, pose) {
  const ov = camOverview(W, H, r);
  if (CAM.mode === 'genel') return Object.assign(ov, { blend: 1 });
  const fo = camFollow(W, H, pose);
  const f = r.left && t > r.time ? smooth((t - r.time - 0.3) / 0.9) : 0;
  if (f <= 0) return Object.assign(fo, { blend: 0 });
  // glide from the tracking shot to the full-flight view after touchdown
  const c = (q) => [(W / 2 - q.ox) / q.k, (q.oy - H / 2) / q.k];
  const a = c(fo), b = c(ov), k = Math.exp(lerp(Math.log(fo.k), Math.log(ov.k), f));
  const cx = lerp(a[0], b[0], f), cy = lerp(a[1], b[1], f);
  return { k, ox: W / 2 - cx * k, oy: H / 2 + cy * k, blend: f };
}
function camClose(W, H, r, pose) {
  const k = Math.min(0.84 * W, H - 62) / (r.G.L + 0.06);
  const cyw = Math.max(pose.cy, (H / 2 - 30) / k);
  return { k, ox: W / 2 - pose.cx * k, oy: H / 2 + cyw * k };
}

// ------------------------------------------------------------------ backdrop: sky, sun, clouds, hills, trees
function drawBackdrop(ctx, W, H, cam, fs, alt) {
  const k = cam.k, gy = cam.oy, restG = 0.86 * H, hs = H / 440;
  const cxw = (W / 2 - cam.ox) / k;
  const lay = (f, d) => restG + (gy - restG) * f * fs - d * H;
  const shiftX = (f) => cxw * k * f * fs;
  const mtnY = lay(0.06, 0.075), hillY = lay(0.16, 0.045), treeY = lay(0.38, 0.016);
  const a = clamp(alt / 70, 0, 1);
  const sky = ctx.createLinearGradient(0, 0, 0, Math.max(mtnY, 60));
  sky.addColorStop(0, mixHex('#2b6dc1', '#123a7e', a));
  sky.addColorStop(0.6, mixHex('#6aa9e0', '#3d80cc', a));
  sky.addColorStop(1, '#d9ecf7');
  ctx.fillStyle = sky; ctx.fillRect(0, 0, W, H);
  // sun
  const sx = W * 0.82, sy = Math.min(H * 0.16, mtnY - 40);
  const sg = ctx.createRadialGradient(sx, sy, 0, sx, sy, 120 * hs);
  sg.addColorStop(0, 'rgba(255,250,232,1)'); sg.addColorStop(0.08, 'rgba(255,246,220,0.95)'); sg.addColorStop(0.2, 'rgba(255,236,190,0.35)'); sg.addColorStop(1, 'rgba(255,236,190,0)');
  ctx.fillStyle = sg; ctx.fillRect(sx - 120 * hs, sy - 120 * hs, 240 * hs, 240 * hs);
  // clouds
  const tileW = W + 700;
  for (let i = 0; i < 9; i++) {
    const spr = SPR.clouds[i % 4], sc = (0.45 + 0.7 * hash(i * 9.1)) * hs, f = 0.05 + 0.06 * hash(i * 3.3);
    let x = (hash(i * 7.7) * tileW - shiftX(f)) % tileW; if (x < 0) x += tileW; x -= 350;
    const y = (-0.55 + 1.15 * hash(i * 5.1)) * H * 0.62 + (gy - restG) * f * fs;
    if (y > mtnY || y + 130 * sc < 0) continue;
    ctx.globalAlpha = 0.6 + 0.35 * hash(i * 2.2); ctx.drawImage(spr, x, y, 300 * sc, 130 * sc);
  }
  ctx.globalAlpha = 1;
  if (mtnY > H + 10) return;
  const ridge = (X, seed, amp) => amp * (0.55 + 0.25 * Math.sin(X * 0.0047 / hs + seed) + 0.14 * Math.sin(X * 0.013 / hs + seed * 2.1) + 0.06 * Math.sin(X * 0.041 / hs + seed * 3.7));
  const band = (baseY, f, seed, amp, c0, c1) => {
    const sxx = shiftX(f);
    const g = ctx.createLinearGradient(0, baseY - amp, 0, baseY);
    g.addColorStop(0, c0); g.addColorStop(1, c1);
    ctx.beginPath(); ctx.moveTo(0, baseY + 2);
    for (let x = 0; x <= W + 8; x += 8) ctx.lineTo(x, baseY - ridge(x + sxx, seed, amp));
    ctx.lineTo(W + 8, baseY + 2); ctx.closePath(); ctx.fillStyle = g; ctx.fill();
  };
  band(mtnY, 0.06, 1.3, 0.14 * H, '#8eb0cb', '#b7cfe0');
  band(hillY, 0.16, 4.2, 0.075 * H, '#5d8f86', '#88ad9a');
  // far field
  const fld = ctx.createLinearGradient(0, hillY, 0, Math.max(hillY + 1, gy));
  fld.addColorStop(0, '#9dbb83'); fld.addColorStop(1, '#5f9849');
  ctx.fillStyle = fld; ctx.fillRect(0, hillY, W, Math.max(0, Math.min(H, gy) - hillY) + 1);
  const hz = ctx.createLinearGradient(0, hillY - 4, 0, hillY + 0.05 * H);
  hz.addColorStop(0, 'rgba(220,236,247,0.55)'); hz.addColorStop(1, 'rgba(220,236,247,0)');
  ctx.fillStyle = hz; ctx.fillRect(0, hillY - 4, W, 0.05 * H + 4);
  // mowing stripes for depth
  ctx.fillStyle = 'rgba(255,255,255,0.05)';
  for (let j = 1; j < 4; j++) { const yy = lerp(treeY, gy, j / 4); ctx.fillRect(0, yy, W, Math.max(1, (gy - treeY) / 14)); }
  // tree line
  if (treeY < gy && treeY > -40) {
    const sxT = shiftX(0.38), sp = 22 * hs;
    const i0 = Math.floor(sxT / sp) - 1, i1 = Math.floor((sxT + W) / sp) + 1;
    ctx.fillStyle = '#2f5a3a';
    for (let i = i0; i <= i1; i++) {
      if (hash(i * 1.7) < 0.3) continue;
      const x = i * sp - sxT + (hash(i * 2.3) - 0.5) * sp * 0.8, h = (0.018 + 0.028 * hash(i * 4.1)) * H;
      ctx.beginPath(); ctx.ellipse(x, treeY - h * 0.58, h * 0.4, h * 0.55, 0, 0, Math.PI * 2); ctx.fill();
      ctx.fillRect(x - 1, treeY - h * 0.18, 2, h * 0.2);
    }
    ctx.fillStyle = 'rgba(40,80,50,0.35)'; ctx.fillRect(0, treeY - 1, W, 3);
  }
}

// ------------------------------------------------------------------ ground: grass, pad, rail, crew
function person(ctx, X, Y, s, h, col, armTo) {
  // X,Y: feet on ground (px); s: px per m
  const lw = Math.max(1.2, 0.09 * s);
  ctx.strokeStyle = col; ctx.fillStyle = col; ctx.lineCap = 'round'; ctx.lineWidth = lw;
  ctx.beginPath(); ctx.moveTo(X - 0.09 * s, Y); ctx.lineTo(X - 0.03 * s, Y - 0.47 * h * s); ctx.moveTo(X + 0.09 * s, Y); ctx.lineTo(X + 0.03 * s, Y - 0.47 * h * s); ctx.stroke();
  rrect(ctx, X - 0.15 * s, Y - 0.82 * h * s, 0.3 * s, 0.37 * h * s, 0.08 * s); ctx.fill();
  ctx.beginPath(); ctx.arc(X, Y - 0.9 * h * s, 0.1 * s, 0, Math.PI * 2); ctx.fill();
  const shY = Y - 0.78 * h * s;
  ctx.beginPath();
  if (armTo) { ctx.moveTo(X + 0.1 * s, shY); ctx.lineTo(X + 0.1 * s + armTo[0] * 0.6 * s, shY - armTo[1] * 0.6 * s); }
  else { ctx.moveTo(X + 0.12 * s, shY); ctx.lineTo(X + 0.16 * s, shY + 0.32 * h * s); }
  ctx.moveTo(X - 0.12 * s, shY); ctx.lineTo(X - 0.16 * s, shY + 0.32 * h * s);
  ctx.stroke(); ctx.lineCap = 'butt';
}
function drawGround(ctx, W, H, cam, t, r, close, pose) {
  const k = cam.k, gy = cam.oy, f = F();
  const SX = (x) => cam.ox + x * k, SY = (y) => gy - y * k;
  const xL = -cam.ox / k - 1, xR = (W - cam.ox) / k + 1;
  const vis = (x0, x1, yTop) => x1 >= xL && x0 <= xR && SY(yTop) < H + 2;
  if (gy < H) {
    const g = ctx.createLinearGradient(0, gy, 0, H);
    g.addColorStop(0, '#4f8c3c'); g.addColorStop(1, '#2b5524');
    ctx.fillStyle = g; ctx.fillRect(0, gy, W, H - gy);
    ctx.fillStyle = 'rgba(190,230,140,0.35)'; ctx.fillRect(0, gy, W, 1);
  }
  // wet patch from the water jet
  if (t > 0 && r.tWaterEnd && gy < H + 6) {
    const wa = 0.4 * clamp(t / Math.max(0.02, r.tWaterEnd), 0, 1);
    ctx.fillStyle = `rgba(30,60,40,${wa})`;
    ctx.beginPath(); ctx.ellipse(SX(-0.2), gy + Math.max(1.5, 0.02 * k), 1.1 * k, Math.max(2, 0.05 * k), 0, 0, Math.PI * 2); ctx.fill();
  }
  // distance stakes
  const step = 5 * k >= 48 ? 5 : 10 * k >= 48 ? 10 : 20;
  ctx.font = `500 ${close ? 12 : 10.5}px ${f.m}`; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (let n = Math.ceil(xL / step); n <= Math.floor(xR / step); n++) {
    if (n === 0) continue;
    const X = SX(n * step), hpx = Math.max(7, 0.45 * k);
    if (gy - hpx > H) continue;
    ctx.fillStyle = '#e9eef2'; ctx.fillRect(X - Math.max(1, 0.012 * k), gy - hpx, Math.max(2, 0.024 * k), hpx);
    ctx.fillStyle = '#ff6b4a'; ctx.fillRect(X - Math.max(1, 0.012 * k), gy - hpx, Math.max(2, 0.024 * k), hpx * 0.3);
    if (gy + 16 < H && (close || X > 175 || gy + 5 < H - 50)) { ctx.fillStyle = 'rgba(240,248,255,0.85)'; ctx.fillText(nf(n * step, 0) + ' m', X, gy + 5); }
  }
  // compressor
  if (vis(-4.1, -3.3, 0.4)) {
    ctx.fillStyle = '#e1a43a'; rrect(ctx, SX(-3.95), SY(0.38), 0.5 * k, 0.3 * k, 0.05 * k); ctx.fill();
    ctx.fillStyle = 'rgba(0,0,0,0.25)'; ctx.fillRect(SX(-3.95), SY(0.12), 0.5 * k, 0.04 * k);
    ctx.fillStyle = '#1d2630';
    [-3.85, -3.55].forEach((x) => { ctx.beginPath(); ctx.arc(SX(x), SY(0.06), Math.max(1.5, 0.06 * k), 0, Math.PI * 2); ctx.fill(); });
    const gx = SX(-3.62), gyy = SY(0.27), gr = Math.max(1.5, 0.055 * k);
    ctx.fillStyle = '#f4f7fa'; ctx.beginPath(); ctx.arc(gx, gyy, gr, 0, Math.PI * 2); ctx.fill();
    const pr = t <= 0 ? r.p.P / 12 : 0, ang = Math.PI * (0.8 + 1.4 * pr);
    ctx.strokeStyle = '#c0392b'; ctx.lineWidth = Math.max(1, 0.01 * k);
    ctx.beginPath(); ctx.moveTo(gx, gyy); ctx.lineTo(gx + Math.cos(ang) * gr * 0.8, gyy + Math.sin(ang) * gr * 0.8); ctx.stroke();
  }
  // hose
  if (vis(-3.5, 0.1, 0.2)) {
    ctx.strokeStyle = '#151b22'; ctx.lineWidth = Math.max(1.2, 0.016 * k);
    ctx.beginPath(); ctx.moveTo(SX(-3.45), SY(0.2)); ctx.bezierCurveTo(SX(-2.4), SY(-0.02), SX(-1), SY(0.03), SX(-0.12), SY(0.05)); ctx.stroke();
  }
  // wind sock
  if (vis(-3.4, -1.6, 2.4)) {
    const px = SX(-2.4), top = SY(2.2);
    ctx.strokeStyle = '#c9d3dc'; ctx.lineWidth = Math.max(1.2, 0.02 * k);
    ctx.beginPath(); ctx.moveTo(px, gy); ctx.lineTo(px, top); ctx.stroke();
    const w = r.p.wind, dir = w >= 0 ? 1 : -1, m = Math.min(1, Math.abs(w) / 7);
    const droop = (1 - m) * 1.25 + 0.05 * Math.sin(9 * t + 1) * m;
    const len = 0.7 * k, r0 = Math.max(2, 0.09 * k), r1 = Math.max(1, 0.04 * k);
    ctx.save(); ctx.translate(px, top); ctx.rotate(dir > 0 ? droop : Math.PI - droop); if (dir < 0) ctx.scale(1, -1);
    for (let i = 0; i < 5; i++) {
      const u0 = len * i / 5, u1 = len * (i + 1) / 5, a0 = lerp(r0, r1, i / 5), a1 = lerp(r0, r1, (i + 1) / 5);
      const wob = 0.05 * k * Math.sin(12 * t + i) * m;
      ctx.beginPath(); ctx.moveTo(u0, -a0); ctx.lineTo(u1, -a1 + wob); ctx.lineTo(u1, a1 + wob); ctx.lineTo(u0, a0); ctx.closePath();
      ctx.fillStyle = i % 2 ? '#f4f6f8' : '#ff7a2f'; ctx.fill();
    }
    ctx.restore();
  }
  // crew behind the 5 m line
  if (vis(-6.9, -4.8, 2)) {
    const post = SX(-5);
    ctx.fillStyle = '#f4f6f8'; ctx.fillRect(post - Math.max(1, 0.015 * k), SY(0.9), Math.max(2, 0.03 * k), 0.9 * k);
    for (let i = 0; i < 3; i++) { ctx.fillStyle = '#e63946'; ctx.fillRect(post - Math.max(1, 0.015 * k), SY(0.9 - i * 0.3), Math.max(2, 0.03 * k), 0.15 * k); }
    const look = t > 0 && !pose.landed;
    const dx = pose.cx + 5.5, dy = pose.cy - 1.4, dl = Math.hypot(dx, dy) || 1;
    person(ctx, SX(-5.5), gy, k, 1.76, '#18222e', look ? [dx / dl, dy / dl] : null);
    person(ctx, SX(-6.0), gy, k, 1.64, '#223246', null);
    person(ctx, SX(-6.55), gy, k, 1.82, '#1b2735', look ? [0.55, 0.83] : null);
    if (!close && k > 18) {
      ctx.font = `500 10.5px ${f.m}`; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillStyle = 'rgba(240,248,255,0.9)';
      ctx.fillText('ekip · 5 m', SX(-6), SY(2.05));
    }
  }
  // launcher: legs, coupler and the guide rail (on the +normal side of the rocket)
  if (vis(-0.6, 1.4, 1.3)) {
    const th0 = r.th0, lw = Math.max(1.4, 0.022 * k);
    ctx.strokeStyle = '#5b6a78'; ctx.lineWidth = lw;
    ctx.beginPath(); ctx.moveTo(SX(-0.24), gy); ctx.lineTo(SX(0), SY(PAD_H)); ctx.lineTo(SX(0.24), gy); ctx.stroke();
    ctx.fillStyle = '#3a4756'; ctx.fillRect(SX(-0.09), SY(PAD_H), 0.18 * k, Math.max(2, 0.03 * k));
    const nx = -Math.sin(th0), ny = Math.cos(th0), off = r.G.R * (close ? EX_CLOSE : 1) + 0.014;
    const rx0 = 0 + nx * off, ry0 = PAD_H + ny * off;
    ctx.strokeStyle = '#d6dee6'; ctx.lineWidth = Math.max(1.5, 0.012 * k);
    ctx.beginPath(); ctx.moveTo(SX(rx0), SY(ry0 - 0.05)); ctx.lineTo(SX(rx0 + Math.cos(th0)), SY(ry0 + Math.sin(th0))); ctx.stroke();
    ctx.strokeStyle = '#7d8a97'; ctx.lineWidth = Math.max(1, 0.006 * k);
    ctx.beginPath(); ctx.moveTo(SX(rx0), SY(ry0 - 0.05)); ctx.lineTo(SX(0.05), gy); ctx.stroke();
    if (r.p.tube) {
      // launch tube: 0,30 m pipe from the launcher into the nozzle (diameter nozzle − 2 mm)
      const y0 = PAD_H - NOZ_L, tl = 0.30, tw = Math.max(1.6, (r.p.nozzle - 2) / 1000 * k * (close ? EX_CLOSE : 1));
      ctx.strokeStyle = '#a7b6c4'; ctx.lineWidth = tw; ctx.lineCap = 'butt';
      ctx.beginPath(); ctx.moveTo(SX(0), SY(y0)); ctx.lineTo(SX(Math.cos(th0) * tl), SY(y0 + Math.sin(th0) * tl)); ctx.stroke();
      ctx.lineCap = 'round';
    }
  }
}

// ------------------------------------------------------------------ exhaust: water jet, droplets, mist
function drawExhaust(ctx, cam, r, t, pose, close, jetR) {
  const em = r.emit; if (!em || !em.length || t <= 0) return;
  const k = cam.k, w = r.p.wind, W = ctx.canvas._w, H = ctx.canvas._h;
  const SX = (x) => cam.ox + x * k, SY = (y) => cam.oy - y * k;
  const tW = r.tWaterEnd || 0, tA = r.tAirEnd || 0, P1 = r.P1g || 1;
  const n = pose.now;
  // coherent jet / plume at the nozzle
  if (!pose.pad && t < tA && n.ph >= 1) {
    const water = n.ph === 1, P = Math.max(0, n.P);
    const ve = r.G.cd * Math.sqrt(200 * P);
    const len = (water ? Math.min(ve * 0.045, 1.6) : 0.5 * clamp(P / P1, 0.1, 1)) * k;
    const bx = pose.bx, by = pose.by, X = SX(pose.x - bx * NOZ_L), Y = SY(pose.y - by * NOZ_L);
    const dx = -bx, dy = by, px = -dy, py = dx;
    const r0 = jetR, r1 = jetR * (water ? 2.6 : 5) + 0.03 * k;
    const g = ctx.createLinearGradient(X, Y, X + dx * len, Y + dy * len);
    if (water) { g.addColorStop(0, 'rgba(240,250,255,0.95)'); g.addColorStop(0.35, 'rgba(205,236,255,0.75)'); g.addColorStop(1, 'rgba(200,232,255,0)'); }
    else { const a = 0.6 * clamp(P / P1, 0.15, 1); g.addColorStop(0, `rgba(255,255,255,${a})`); g.addColorStop(1, 'rgba(255,255,255,0)'); }
    ctx.save(); ctx.beginPath(); ctx.rect(0, 0, W, cam.oy); ctx.clip();
    ctx.fillStyle = g; ctx.beginPath();
    ctx.moveTo(X + px * r0, Y + py * r0); ctx.lineTo(X + dx * len + px * r1, Y + dy * len + py * r1);
    ctx.lineTo(X + dx * len - px * r1, Y + dy * len - py * r1); ctx.lineTo(X - px * r0, Y - py * r0); ctx.closePath(); ctx.fill();
    ctx.restore();
  }
  // mist puffs from the air phase
  const lifeM = 2.0;
  for (const s of em) {
    if (s.ph !== 2) continue;
    const a = t - s.tau; if (a <= 0 || a > lifeM) continue;
    const burst = s.tau - tW < 0.005 ? 6 : 3;
    for (let j = 0; j < burst; j++) {
      const h1 = hash(s.i * 3.1 + j * 17.7), h2 = hash(s.i * 7.9 + j * 5.3 + 2), h3 = hash(s.i * 1.3 + j * 9.1 + 4);
      const sp = (h1 - 0.5) * 0.8, cs = Math.cos(sp), sn = Math.sin(sp);
      const dx = -(s.bx * cs - s.by * sn), dy = -(s.bx * sn + s.by * cs);
      const spd = (5 + 16 * h2) * (0.35 + 0.65 * clamp(s.P / P1, 0, 1));
      const v0x = s.vx * 0.15 + dx * spd, v0y = s.vy * 0.15 + dy * spd;
      const kd = 3.2, vix = w * 0.8, viy = 0.25, E = (1 - Math.exp(-kd * a)) / kd;
      const x = s.x + vix * a + (v0x - vix) * E, y = Math.max(0.05, s.y + viy * a + (v0y - viy) * E);
      const rad = Math.min(close ? 170 : 90, (0.05 + 0.5 * (1 - Math.exp(-1.8 * a)) * (0.6 + 0.7 * h3)) * k);
      const X = SX(x), Y = SY(y);
      if (X + rad < 0 || X - rad > W || Y + rad < 0 || Y - rad > H) continue;
      ctx.globalAlpha = 0.34 * Math.pow(1 - a / lifeM, 1.6);
      ctx.drawImage(SPR.puff, X - rad, Y - rad, rad * 2, rad * 2);
    }
  }
  ctx.globalAlpha = 1;
  // water droplets: ballistic with drag, splashing when they reach the ground
  const lifeD = 1.5, nd = close ? 5 : 4, rmin = close ? 1.3 : 0.9;
  ctx.fillStyle = '#e4f4ff';
  for (const s of em) {
    if (s.ph !== 1) continue;
    const a = t - s.tau; if (a <= 0 || a > lifeD) continue;
    const ve = r.G.cd * Math.sqrt(200 * Math.max(0, s.P));
    for (let j = 0; j < nd; j++) {
      const h1 = hash(s.i * 13.1 + j * 7.7), h2 = hash(s.i * 3.7 + j * 19.3 + 1), h3 = hash(s.i * 5.3 + j * 2.9 + 2), h4 = hash(s.i * 1.9 + j * 11.3 + 3);
      const sp = (h1 - 0.5) * 0.24, cs = Math.cos(sp), sn = Math.sin(sp);
      const dx = -(s.bx * cs - s.by * sn), dy = -(s.bx * sn + s.by * cs);
      const spd = ve * (0.6 + 0.45 * h2);
      const v0x = s.vx + dx * spd, v0y = s.vy + dy * spd;
      const kd = 1.8 + 3.5 * h3, vix = w, viy = -C.g / kd;
      const pos = (aa) => { const E = (1 - Math.exp(-kd * aa)) / kd; return [s.x + vix * aa + (v0x - vix) * E, s.y + viy * aa + (v0y - viy) * E]; };
      let [x, y] = pos(a);
      let alpha = 0.9 * (1 - a / lifeD);
      if (y < 0) {
        let lo = 0, hi = a;
        for (let it = 0; it < 14; it++) { const m = (lo + hi) / 2; if (pos(m)[1] > 0) lo = m; else hi = m; }
        const after = a - lo; if (after > 0.5) continue;
        const ex = Math.exp(-kd * lo), vxi = vix + (v0x - vix) * ex, vyi = viy + (v0y - viy) * ex;
        const spread = 1 - Math.exp(-9 * after);
        x = pos(lo)[0] + vxi * 0.05 * (0.4 + h4) * spread + (h1 - 0.5) * 1.3 * spread;
        y = Math.abs(vyi) * 0.012 * (0.3 + h2) * Math.sin(Math.PI * Math.min(1, after / 0.32));
        alpha = 0.75 * (1 - after / 0.5);
      }
      const X = SX(x), Y = SY(y);
      if (X < -6 || X > W + 6 || Y < -6 || Y > H + 6) continue;
      const rad = clamp((0.005 + 0.012 * h4) * k, rmin, 6);
      ctx.globalAlpha = alpha;
      if (rad < 1.6) ctx.fillRect(X - rad, Y - rad, rad * 2, rad * 2);
      else { ctx.beginPath(); ctx.arc(X, Y, rad, 0, Math.PI * 2); ctx.fill(); }
    }
  }
  ctx.globalAlpha = 1;
}
function drawImpact(ctx, cam, r, t, pose) {
  if (!pose.landed) return;
  const a = t - r.time; if (a > 1.3) return;
  const k = cam.k, last = r.rec[r.rec.length - 1], vx0 = last ? last.vx : 0;
  for (let j = 0; j < 22; j++) {
    const h1 = hash(j * 4.3 + 7), h2 = hash(j * 9.7 + 1), h3 = hash(j * 2.1 + 5);
    const v0x = (h1 - 0.5) * 3.5 + vx0 * 0.15, v0y = 0.8 + 2.8 * h2, kd = 2.4;
    const E = (1 - Math.exp(-kd * a)) / kd, viy = -C.g / kd;
    const x = pose.cX + v0x * E, y = Math.max(0, pose.cY + viy * a + (v0y - viy) * E);
    const rad = clamp((0.01 + 0.02 * h3) * k, 1, 5);
    ctx.globalAlpha = 0.85 * (1 - a / 1.3);
    ctx.fillStyle = j % 3 ? '#6e5a3e' : '#6aa84f';
    ctx.fillRect(cam.ox + x * k - rad / 2, cam.oy - y * k - rad / 2, rad, rad);
  }
  ctx.globalAlpha = 1;
}

// ------------------------------------------------------------------ the rocket (cutaway)
function clipPoly(poly, nx, ny, d) {
  const out = [];
  for (let i = 0; i < poly.length; i++) {
    const A = poly[i], B = poly[(i + 1) % poly.length];
    const da = nx * A[0] + ny * A[1] - d, db = nx * B[0] + ny * B[1] - d;
    if (da >= 0) out.push(A);
    if ((da >= 0) !== (db >= 0)) { const f = da / (da - db); out.push([A[0] + (B[0] - A[0]) * f, A[1] + (B[1] - A[1]) * f]); }
  }
  return out;
}
function polyArea(P) { let a = 0; for (let i = 0; i < P.length; i++) { const A = P[i], B = P[(i + 1) % P.length]; a += A[0] * B[1] - B[0] * A[1]; } return Math.abs(a) / 2; }
function levelFor(hl, hw, nx, ny, frac) {
  // offset d so that the part of the rectangle with n·p >= d holds `frac` of its area
  const P = [[-hl, -hw], [hl, -hw], [hl, hw], [-hl, hw]], A = 4 * hl * hw;
  let lo = -(Math.abs(nx) * hl + Math.abs(ny) * hw), hi = -lo;
  for (let i = 0; i < 26; i++) { const mid = (lo + hi) / 2; if (polyArea(clipPoly(P, nx, ny, mid)) / A > frac) lo = mid; else hi = mid; }
  return (lo + hi) / 2;
}
function drawWater(ctx, cx, hl, hw, nx, ny, frac, rip, phase, cr) {
  if (frac <= 0.002) return;
  ctx.save();
  rrect(ctx, cx - hl, -hw, 2 * hl, 2 * hw, cr); ctx.clip();
  const D = Math.hypot(hl, hw) + 6;
  if (frac >= 0.998) {
    ctx.fillStyle = 'rgba(38,140,215,0.82)'; ctx.fillRect(cx - hl, -hw, 2 * hl, 2 * hw); ctx.restore(); return;
  }
  const d = Math.abs(frac - 0.5) < 1e-6 ? 0 : levelFor(hl, hw, nx, ny, frac);
  const qx = cx + nx * d, qy = ny * d, tx = -ny, ty = nx;
  const kw = Math.PI / Math.max(4, hw) * 0.9;
  const pts = [];
  for (let i = 0; i <= 28; i++) {
    const s = -D + 2 * D * i / 28, wv = rip * Math.sin(s * kw + phase);
    pts.push([qx + tx * s - nx * wv, qy + ty * s - ny * wv]);
  }
  const g = ctx.createLinearGradient(qx, qy, qx + nx * 2 * hw, qy + ny * 2 * hw);
  g.addColorStop(0, 'rgba(110,205,250,0.9)'); g.addColorStop(0.25, 'rgba(40,150,225,0.86)'); g.addColorStop(1, 'rgba(14,82,160,0.92)');
  ctx.fillStyle = g; ctx.beginPath();
  pts.forEach((p, i) => (i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1])));
  ctx.lineTo(qx + tx * D + nx * 2 * D, qy + ty * D + ny * 2 * D); ctx.lineTo(qx - tx * D + nx * 2 * D, qy - ty * D + ny * 2 * D);
  ctx.closePath(); ctx.fill();
  ctx.strokeStyle = 'rgba(225,246,255,0.95)'; ctx.lineWidth = hw > 10 ? 1.6 : 1;
  ctx.beginPath(); pts.forEach((p, i) => (i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]))); ctx.stroke();
  ctx.restore();
  return { qx, qy, tx, ty };
}
// Wing and tail planforms in the rocket's local frame (u from the tail toward the nose, v across).
// The wing is level, so its span points toward the viewer; it is drawn as seen from slightly above
// and to the side (a fixed screen direction), which keeps it visible at every pitch angle.
const WING_VIEW = [-0.3, 0.42];
function wingQuads(o, L, R) {
  const G = o.G, k = o.k, lod = o.lod;
  const e = o.w.wd == null ? 1 : Math.max(0, Math.min(1, o.w.wd));
  const ca = Math.cos(o.ang), sa = Math.sin(o.ang);
  const P = [WING_VIEW[0] * ca + WING_VIEW[1] * sa, -WING_VIEW[0] * sa + WING_VIEW[1] * ca];
  const quad = (uAc, c, half, sgn, fill, stroke) => {
    const le = uAc + 0.25 * c, te = uAc - 0.75 * c, tle = le - 0.1 * c, tte = tle - 0.75 * c;
    const r0 = R * sgn, r1 = (R + half) * sgn;
    return { pts: [[le + P[0] * r0, P[1] * r0], [tle + P[0] * r1, P[1] * r1], [tte + P[0] * r1, P[1] * r1], [te + P[0] * r0, P[1] * r0]], fill, stroke, lw: lod ? 1.2 : 0.8 };
  };
  const aN = lod ? 0.42 : 0.95, aF = lod ? 0.25 : 0.7;
  const uW = L - G.xW * k, c = G.chord * k, half = o.p.wspan / 2 * k * e;
  const tc = Math.sqrt(G.St / 4) * k, th = Math.sqrt(G.St) * k, uT = Math.max(0.75 * tc, L - G.xT * k);
  const near = [quad(uT, tc, th, 1, `rgba(226,234,242,${aN})`, 'rgba(143,161,179,0.95)')];
  const far = [quad(uT, tc, th, -1, `rgba(170,184,198,${aF})`, 'rgba(120,138,156,0.8)')];
  if (half > 0.5) {
    near.push(quad(uW, c, half, 1, `rgba(255,176,74,${aN})`, 'rgba(201,120,42,0.95)'));
    far.push(quad(uW, c, half, -1, `rgba(196,128,52,${aF})`, 'rgba(160,100,40,0.8)'));
  }
  return { near, far, e };
}
function drawQuad(ctx, q) {
  ctx.beginPath(); q.pts.forEach((pt, i) => (i ? ctx.lineTo(pt[0], pt[1]) : ctx.moveTo(pt[0], pt[1]))); ctx.closePath();
  ctx.fillStyle = q.fill; ctx.fill(); ctx.strokeStyle = q.stroke; ctx.lineWidth = q.lw; ctx.stroke();
}
function drawRocket(ctx, o) {
  const G = o.G, p = o.p, k = o.k, ex = o.ex, lod = o.lod, w = o.w;
  const R = G.R * k * ex;
  const L = G.L * k, Lt = G.Ltail * k, Lc = G.Lc * k, Lp = G.Lp * k, Ln = G.Ln * k;
  const uc0 = Lt, uc1 = Lt + Lc, up0 = uc1, up1 = uc1 + Lp, un0 = up1, un1 = L;
  const lw = lod ? 1.2 : 0.8;
  ctx.save(); ctx.translate(o.X, o.Y); ctx.rotate(o.ang);
  if (o.fs && o.fs < 0.999) ctx.scale(o.fs, 1);
  const hi = Math.cos(o.ang) >= 0 ? -1 : 1;
  if (G.wing) wingQuads(o, L, R).far.forEach((q) => drawQuad(ctx, q));
  // fins (profile pair)
  const fs = G.fin.span * k * ex, fc = G.fin.chord * k, fb = 0.012 * k;
  [-1, 1].forEach((sd) => {
    const g = ctx.createLinearGradient(0, sd * R, 0, sd * (R + fs));
    g.addColorStop(0, '#ff9a4f'); g.addColorStop(1, '#e8691f');
    ctx.fillStyle = g; ctx.beginPath();
    ctx.moveTo(fc, sd * R * 0.98); ctx.lineTo(0, sd * R * 0.98); ctx.lineTo(-fb, sd * (R + fs)); ctx.lineTo(0.45 * fc, sd * (R + fs)); ctx.closePath(); ctx.fill();
    if (lod) { ctx.strokeStyle = '#b4521a'; ctx.lineWidth = 1; ctx.stroke(); }
  });
  // nozzle
  const rn1 = clamp(p.nozzle / 2000 * k * ex, 0.6, R * 0.95), rn0 = Math.max(rn1 + 0.002 * k * ex, R * 0.42), nl = NOZ_L * k;
  const ng = ctx.createLinearGradient(0, -rn0, 0, rn0);
  ng.addColorStop(0, '#e1e8ee'); ng.addColorStop(0.5, '#8d9cab'); ng.addColorStop(1, '#4b5865');
  ctx.fillStyle = ng; ctx.beginPath(); ctx.moveTo(1, -rn0); ctx.lineTo(-nl, -rn1); ctx.lineTo(-nl, rn1); ctx.lineTo(1, rn0); ctx.closePath(); ctx.fill();
  // tail section
  ctx.fillStyle = '#25303c'; ctx.fillRect(0, -R, Lt + 2, 2 * R);
  // pressure chamber (PET bottle)
  const cr = Math.min(R * 0.7, Lc * 0.22), ins = lod ? 1.6 : 0.6;
  rrect(ctx, uc0, -R, Lc, 2 * R, cr); ctx.fillStyle = 'rgba(185,222,255,0.16)'; ctx.fill();
  const nC = [Math.cos(w.gc), -Math.sin(w.gc) / ex], mC = Math.hypot(nC[0], nC[1]);
  drawWater(ctx, (uc0 + uc1) / 2, Lc / 2 - ins, R - ins, nC[0] / mC, nC[1] / mC, w.fc, lod ? (R - ins) * 0.05 : 0, w.ph, Math.max(0, cr - ins));
  if (w.mist > 0.01) {
    ctx.save(); rrect(ctx, uc0 + ins, -R + ins, Lc - 2 * ins, 2 * (R - ins), Math.max(0, cr - ins)); ctx.clip();
    ctx.fillStyle = `rgba(236,245,252,${Math.min(0.85, w.mist * 0.95)})`; ctx.fillRect(uc0, -R, Lc, 2 * R);
    if (lod) {
      for (let i = 0; i < 7; i++) {
        const u = uc0 + Lc * (0.1 + 0.8 * hash(i * 3.3)), v = (hash(i * 5.7) - 0.5) * R, rr = R * (0.6 + 0.6 * hash(i));
        ctx.globalAlpha = w.mist * 0.8; ctx.drawImage(SPR.puff, u - rr + Math.sin(w.ph * 0.3 + i) * R * 0.2, v - rr, rr * 2, rr * 2);
      }
      ctx.globalAlpha = 1;
    }
    ctx.restore();
  }
  if (lod && w.fc < 0.02 && w.dry) {
    ctx.fillStyle = 'rgba(210,238,255,0.45)';
    for (let i = 0; i < 18; i++) {
      const u = uc0 + Lc * (0.08 + 0.84 * hash(i * 1.7 + 3)), v = (hash(i * 2.9 + 1) * 2 - 1) * (R - 3);
      ctx.beginPath(); ctx.arc(u, v, 0.8 + 1.4 * hash(i * 6.1), 0, Math.PI * 2); ctx.fill();
    }
  }
  rrect(ctx, uc0, -R, Lc, 2 * R, cr); ctx.strokeStyle = 'rgba(226,242,255,0.8)'; ctx.lineWidth = lw; ctx.stroke();
  // payload (sloshing) tank
  const pr = Math.min(R * 0.5, Lp * 0.18);
  rrect(ctx, up0, -R, Lp, 2 * R, pr); ctx.fillStyle = 'rgba(185,222,255,0.16)'; ctx.fill();
  const nS = [Math.cos(w.gs), -Math.sin(w.gs) / ex], mS = Math.hypot(nS[0], nS[1]);
  const rip = lod ? (R - ins) * (0.04 + 0.2 * Math.min(1, Math.abs(w.s))) : 0;
  const surf = drawWater(ctx, (up0 + up1) / 2, Lp / 2 - ins, R - ins, nS[0] / mS, nS[1] / mS, 0.5, rip, w.ph, Math.max(0, pr - ins));
  if (lod && surf && Math.abs(w.s) > 0.15) {
    // foam where the wave breaks
    ctx.save(); rrect(ctx, up0 + ins, -R + ins, Lp - 2 * ins, 2 * (R - ins), pr); ctx.clip();
    ctx.fillStyle = `rgba(240,250,255,${Math.min(0.8, Math.abs(w.s) * 0.9)})`;
    for (let i = 0; i < 12; i++) {
      const sOff = (hash(i * 7.3) - 0.5) * 2 * Math.hypot(Lp / 2, R), dep = 1.5 + 6 * hash(i * 1.1);
      const fx = surf.qx + surf.tx * sOff + nS[0] / mS * dep, fy = surf.qy + surf.ty * sOff + nS[1] / mS * dep;
      ctx.beginPath(); ctx.arc(fx, fy, 0.8 + 1.6 * hash(i * 3.9), 0, Math.PI * 2); ctx.fill();
    }
    ctx.restore();
  }
  if (p.baffle > 0 && R >= 3) {
    const stub = (R - ins) * (1 - Math.sqrt(Math.max(0, 1 - p.baffle))), thk = Math.max(lod ? 2.2 : 1, 0.004 * k);
    ctx.fillStyle = 'rgba(238,244,250,0.95)';
    [1 / 3, 2 / 3].forEach((f) => {
      const bx = up0 + Lp * f - thk / 2;
      ctx.fillRect(bx, -R + ins, thk, stub); ctx.fillRect(bx, R - ins - stub, thk, stub);
    });
  }
  rrect(ctx, up0, -R, Lp, 2 * R, pr); ctx.strokeStyle = 'rgba(226,242,255,0.8)'; ctx.lineWidth = lw; ctx.stroke();
  // coupler between the two bottles
  const cw = Math.max(2, 0.026 * k);
  ctx.fillStyle = '#34424f'; ctx.fillRect(uc1 - cw / 2, -R * 1.05, cw, R * 2.1);
  if (lod) { ctx.fillStyle = '#9fb0bf'; [-0.6, 0, 0.6].forEach((v) => { ctx.beginPath(); ctx.arc(uc1, v * R, Math.max(0.8, cw * 0.16), 0, Math.PI * 2); ctx.fill(); }); }
  // nose cone
  const nosePath = () => {
    ctx.beginPath(); ctx.moveTo(un0, -R);
    ctx.bezierCurveTo(un0 + 0.5 * Ln, -R, un1 - 0.1 * Ln, -0.32 * R, un1, 0);
    ctx.bezierCurveTo(un1 - 0.1 * Ln, 0.32 * R, un0 + 0.5 * Ln, R, un0, R); ctx.closePath();
  };
  nosePath();
  const ngr = ctx.createLinearGradient(0, hi * R, 0, -hi * R);
  ngr.addColorStop(0, '#ffffff'); ngr.addColorStop(0.5, '#e3e9ef'); ngr.addColorStop(1, '#9eabb8');
  ctx.fillStyle = ngr; ctx.fill();
  ctx.fillStyle = '#34424f'; ctx.fillRect(un0 - 1, -R * 1.02, Math.max(1.5, 0.008 * k), R * 2.04);
  if (lod) {
    // altimeter window and vent holes
    rrect(ctx, un0 + 0.16 * Ln, -0.44 * R, 0.3 * Ln, 0.88 * R, 3); ctx.fillStyle = '#0d1b28'; ctx.fill();
    ctx.fillStyle = '#1f8a4c'; ctx.fillRect(un0 + 0.19 * Ln, -0.3 * R, 0.24 * Ln, 0.6 * R);
    ctx.fillStyle = '#d8b65a'; ctx.fillRect(un0 + 0.26 * Ln, -0.12 * R, 0.08 * Ln, 0.24 * R);
    ctx.fillStyle = '#1b2632'; [-0.35, 0, 0.35].forEach((v) => { ctx.beginPath(); ctx.arc(un0 + 0.6 * Ln, v * R * 0.8, Math.max(0.8, 0.0015 * k), 0, Math.PI * 2); ctx.fill(); });
  }
  // cylindrical shading over the whole body
  ctx.save();
  ctx.beginPath(); ctx.rect(0, -R, un0, 2 * R); ctx.clip();
  const sg = ctx.createLinearGradient(0, hi * R, 0, -hi * R);
  sg.addColorStop(0, 'rgba(255,255,255,0.32)'); sg.addColorStop(0.22, 'rgba(255,255,255,0.06)'); sg.addColorStop(0.62, 'rgba(0,0,0,0)'); sg.addColorStop(1, 'rgba(0,10,25,0.3)');
  ctx.fillStyle = sg; ctx.fillRect(0, -R, un0, 2 * R);
  if (lod) { ctx.fillStyle = 'rgba(255,255,255,0.45)'; ctx.fillRect(uc0 + cr, hi * R * 0.6 - 0.75, up1 - uc0 - cr - pr, 1.5); }
  ctx.restore();
  // near half of the wing and tail (the far half was drawn behind the body); folded panels lie along the body
  if (G.wing) {
    const wq = wingQuads(o, L, R);
    wq.near.forEach((q) => drawQuad(ctx, q));
    if (wq.e < 0.999) {
      const half = o.p.wspan / 2 * k, pack = Math.min(half, 0.4 * L), c = G.chord * k, uW = L - G.xW * k;
      const x1 = Math.min(L - 4, uW + 0.25 * c), h = Math.min(c, 2 * R) * (lod ? 0.9 : 0.7);
      ctx.globalAlpha = 1 - wq.e;
      rrect(ctx, x1 - pack, -h / 2, pack, h, Math.min(4, h / 2));
      ctx.fillStyle = `rgba(255,176,74,${lod ? 0.3 : 0.75})`; ctx.fill();
      ctx.strokeStyle = 'rgba(255,176,74,0.95)'; ctx.lineWidth = lod ? 1.4 : 0.8; ctx.stroke();
      ctx.globalAlpha = 1;
    }
  }
  ctx.restore();
}
function rocketParts(G, X, Y, ang, k) {
  const at = (u) => [X + Math.cos(ang) * u * k, Y + Math.sin(ang) * u * k];
  return {
    pay: at(G.Ltail + G.Lc + G.Lp / 2),
    ch: at(G.Ltail + G.Lc / 2),
  };
}

// ------------------------------------------------------------------ parachute
function drawChute(ctx, cam, r, t, pose, minR) {
  if (r.tChute == null || t < r.tChute) return;
  const k = cam.k, e = smooth((t - r.tChute) / 0.5);
  let ax = pose.nx + 0.12 * Math.sin(1.7 * t), ay = pose.ny + 1.6, flat = 1;
  if (pose.landed) {
    const f = smooth((t - r.time) / 0.9);
    ax = lerp(ax, pose.nx + 0.9 * (r.p.wind >= 0 ? 1 : -1), f); ay = lerp(ay, 0.08, f); flat = 1 - 0.75 * f;
  }
  const Rr = Math.max(minR, 0.55 * k) * (0.3 + 0.7 * e);
  const X = cam.ox + ax * k, Y = cam.oy - ay * k, NX = cam.ox + pose.nx * k, NY = cam.oy - pose.ny * k;
  ctx.strokeStyle = 'rgba(232,238,244,0.8)'; ctx.lineWidth = 1;
  ctx.beginPath(); [-1, -0.5, 0, 0.5, 1].forEach((f) => { ctx.moveTo(X + f * Rr, Y); ctx.lineTo(NX, NY); }); ctx.stroke();
  for (let i = 0; i < 6; i++) {
    const a0 = Math.PI + i * Math.PI / 6, a1 = a0 + Math.PI / 6;
    ctx.beginPath(); ctx.moveTo(X, Y); ctx.ellipse(X, Y, Rr, Rr * 0.62 * flat, 0, a0, a1); ctx.closePath();
    ctx.fillStyle = i % 2 ? '#f3f6f9' : '#ff7a2f'; ctx.fill();
  }
  ctx.strokeStyle = 'rgba(0,0,0,0.25)'; ctx.beginPath(); ctx.moveTo(X - Rr, Y); ctx.lineTo(X + Rr, Y); ctx.stroke();
}

// ------------------------------------------------------------------ main camera
const PH_COL = ['#c7d3de', '#62d2ff', '#ffc061', '#ffffff'];
function drawPaths(ctx, cam, r, t) {
  const k = cam.k, SX = (x) => cam.ox + x * k, SY = (y) => cam.oy - (y + PAD_H) * k;
  const ref = S.ref;
  ctx.lineCap = 'round'; ctx.lineJoin = 'round';
  if (ref && ref !== r && ref.rec.length && (Math.abs(ref.apo - r.apo) > 1e-6 || Math.abs(ref.dist - r.dist) > 1e-6)) {
    ctx.setLineDash([7, 6]); ctx.strokeStyle = 'rgba(255,255,255,0.55)'; ctx.lineWidth = 1.6; ctx.beginPath();
    ref.rec.forEach((e, i) => (i ? ctx.lineTo(SX(e.x), SY(e.y)) : ctx.moveTo(SX(e.x), SY(e.y)))); ctx.stroke();
  }
  ctx.setLineDash([2, 6]); ctx.strokeStyle = 'rgba(255,255,255,0.4)'; ctx.lineWidth = 1.4; ctx.beginPath();
  r.rec.forEach((e, i) => (i ? ctx.lineTo(SX(e.x), SY(e.y)) : ctx.moveTo(SX(e.x), SY(e.y)))); ctx.stroke();
  ctx.setLineDash([]);
  const tEnd = Math.min(t, r.time);
  ctx.strokeStyle = 'rgba(255,255,255,0.22)'; ctx.lineWidth = 6; ctx.beginPath();
  let started = false;
  for (const e of r.rec) { if (e.t > tEnd) break; if (started) ctx.lineTo(SX(e.x), SY(e.y)); else { ctx.moveTo(SX(e.x), SY(e.y)); started = true; } }
  ctx.stroke();
  ctx.lineWidth = 2.4;
  let prev = null;
  for (const e of r.rec) {
    if (e.t > tEnd) break;
    const P = [SX(e.x), SY(e.y)];
    if (prev) { ctx.strokeStyle = PH_COL[e.chute ? 3 : e.ph] || '#fff'; ctx.beginPath(); ctx.moveTo(prev[0], prev[1]); ctx.lineTo(P[0], P[1]); ctx.stroke(); }
    prev = P;
  }
  ctx.lineCap = 'butt';
}
function drawMarkers(ctx, W, H, cam, r, t) {
  const k = cam.k, f = F(), SX = (x) => cam.ox + x * k, SY = (y) => cam.oy - y * k;
  if (!r.left) return;
  if (t >= r.tApo) {
    const a = recAt(r, r.tApo), X = SX(a.x), Y = SY(a.y + PAD_H);
    if (cam.blend > 0.5) {
      ctx.setLineDash([3, 5]); ctx.strokeStyle = 'rgba(255,255,255,0.45)'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(X, Y); ctx.lineTo(X, cam.oy); ctx.stroke(); ctx.setLineDash([]);
    }
    ctx.fillStyle = '#ffffff'; ctx.beginPath(); ctx.moveTo(X, Y - 5); ctx.lineTo(X + 5, Y); ctx.lineTo(X, Y + 5); ctx.lineTo(X - 5, Y); ctx.closePath(); ctx.fill();
    pill(ctx, X + (X > W * 0.7 ? -10 : 10), Y - 14, 'tepe ' + nf(r.apo, 1) + ' m', { align: X > W * 0.7 ? 'right' : 'left', font: `500 11px ${f.m}` });
  }
  if (t >= r.time) {
    const last = r.rec[r.rec.length - 1], X = SX(last.x), X0 = SX(0), Y = cam.oy - 4;
    ctx.strokeStyle = '#ffd25e'; ctx.lineWidth = 2; ctx.setLineDash([8, 5]);
    ctx.beginPath(); ctx.moveTo(X0, Y); ctx.lineTo(X, Y); ctx.stroke(); ctx.setLineDash([]);
    ctx.beginPath(); ctx.moveTo(X0, Y - 6); ctx.lineTo(X0, Y + 6); ctx.moveTo(X, Y - 6); ctx.lineTo(X, Y + 6); ctx.stroke();
    const mid = (X0 + X) / 2;
    if (mid > -40 && mid < W + 40) pill(ctx, clamp(mid, 110, W - 110), Y - 16, Math.abs(r.zLand) >= 0.3 ? `ileri ${nf(Math.abs(r.xLand), 1)} · yana ${nf(Math.abs(r.zLand), 1)} → mesafe ${nf(r.dist, 1)} m` : 'mesafe ' + nf(r.dist, 1) + ' m', { align: 'center', color: '#ffe3a0', font: `500 11px ${f.m}` });
  }
}
function drawMainHUD(ctx, W, H, cam, r, t, pose) {
  const f = F(), n = pose.now, k = cam.k;
  // altitude ruler (right edge)
  const yTop = cam.oy / k, yBot = (cam.oy - H) / k;
  const st = k > 14 ? 1 : k > 5 ? 5 : 10, lab = st === 1 ? 5 : st * 2;
  ctx.strokeStyle = 'rgba(255,255,255,0.55)'; ctx.fillStyle = 'rgba(255,255,255,0.85)'; ctx.lineWidth = 1;
  ctx.font = `500 10px ${f.m}`; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  for (let y = Math.max(st, Math.ceil(yBot / st) * st); y <= yTop; y += st) {
    const Y = cam.oy - (y + PAD_H) * k; if (Y < 92) continue;
    const major = Math.abs(y / lab - Math.round(y / lab)) < 1e-6;
    ctx.beginPath(); ctx.moveTo(W - (major ? 12 : 6), Y); ctx.lineTo(W, Y); ctx.stroke();
    if (major) ctx.fillText(nf(y, 0) + ' m', W - 15, Y);
  }
  // time and phase
  const ph = phaseOf(n);
  const pc = ph.code === 1 ? '#7fdcff' : ph.code === 2 ? '#ffc874' : '#d7e4ef';
  ctx.font = `600 19px ${f.m}`;
  const tt = 'T+' + nf(Math.min(t, r.time), 2) + ' s';
  ctx.font = `500 10.5px ${f.m}`;
  const pw = ctx.measureText(ph.label.toLocaleUpperCase('tr-TR')).width;
  ctx.font = `600 19px ${f.m}`;
  const bw = Math.max(ctx.measureText(tt).width, pw) + 22;
  rrect(ctx, 10, 10, bw, 52, 9); ctx.fillStyle = 'rgba(4,12,22,0.62)'; ctx.fill();
  ctx.fillStyle = '#ffffff'; ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic'; ctx.fillText(tt, 21, 34);
  ctx.font = `500 10.5px ${f.m}`; ctx.fillStyle = pc; ctx.fillText(ph.label.toLocaleUpperCase('tr-TR'), 21, 52);
  // readouts
  const onPad = pose.pad, done = pose.landed;
  const rows = [['İRTİFA', nf(onPad || done ? 0 : Math.max(0, n.y), 1) + ' m'], ['HIZ', nf(onPad || done ? 0 : n.v, 1) + ' m/s'], ['İVME', nf(onPad || done ? 0 : n.ax, 1) + ' g']];
  if (W >= 420) {
    ctx.font = `500 12px ${f.m}`;
    const vw = Math.max(...rows.map((q) => ctx.measureText(q[1]).width)) + 66;
    rrect(ctx, W - vw - 10, 10, vw, 66, 9); ctx.fillStyle = 'rgba(4,12,22,0.62)'; ctx.fill();
    rows.forEach((q, i) => {
      const Y = 28 + i * 18;
      ctx.font = `500 9.5px ${f.m}`; ctx.fillStyle = 'rgba(200,220,235,0.75)'; ctx.textAlign = 'left'; ctx.fillText(q[0], W - vw, Y);
      ctx.font = `500 12px ${f.m}`; ctx.fillStyle = '#fff'; ctx.textAlign = 'right'; ctx.fillText(q[1], W - 20, Y);
    });
  }
  // wind and scale
  const w = r.p.wind, wz = r.p.windZ || 0;
  const sb = [0.5, 1, 2, 5, 10, 20, 50].find((m) => m * k >= 46) || 50, sbp = sb * k;
  ctx.font = `500 10px ${f.m}`;
  const wtxt = (w === 0 && wz === 0 ? 'RÜZGÂR YOK' : `RÜZGÂR ${nf(Math.abs(w), 1)} m/s ${w > 0 ? '→' : w < 0 ? '←' : ''}`) + (wz !== 0 ? ` · YAN ${nf(Math.abs(wz), 1)} m/s (3B'de görünür)` : '');
  rrect(ctx, 10, H - 44, Math.max(sbp, 120, ctx.measureText(wtxt).width) + 20, 34, 8); ctx.fillStyle = 'rgba(4,12,22,0.55)'; ctx.fill();
  ctx.fillStyle = '#e6f1fa'; ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
  ctx.fillText(wtxt, 20, H - 29);
  ctx.strokeStyle = '#e6f1fa'; ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(20, H - 21); ctx.lineTo(20, H - 16); ctx.lineTo(20 + sbp, H - 16); ctx.lineTo(20 + sbp, H - 21); ctx.stroke();
  ctx.textAlign = 'left'; ctx.fillText(nf(sb, sb < 1 ? 1 : 0) + ' m', 26 + sbp, H - 14);
}
function drawScene(r, t, pose) {
  const W = scene._w, H = scene._h; if (!W) return;
  const ctx = sctx, cam = mainCam(W, H, r, t, pose);
  ctx.clearRect(0, 0, W, H);
  drawBackdrop(ctx, W, H, cam, 1, Math.max(0, (cam.oy - H / 2) / cam.k));
  if (cam.blend > 0.4) {
    const k = cam.k, ymax = (cam.oy) / k, st = ymax > 60 ? 20 : ymax > 25 ? 10 : 5;
    ctx.strokeStyle = `rgba(255,255,255,${0.12 * cam.blend})`; ctx.lineWidth = 1;
    for (let y = st; y < ymax; y += st) { const Y = cam.oy - y * k; ctx.beginPath(); ctx.moveTo(0, Y); ctx.lineTo(W - 40, Y); ctx.stroke(); }
  }
  drawGround(ctx, W, H, cam, t, r, false, pose);
  drawPaths(ctx, cam, r, t);
  // shadow under the rocket
  if (pose.cy < 12) {
    const a = 0.3 * (1 - pose.cy / 12), X = cam.ox + pose.cx * cam.k;
    ctx.fillStyle = `rgba(10,30,15,${a})`; ctx.beginPath();
    ctx.ellipse(X, cam.oy + 1, Math.max(4, Math.abs(Math.cos(pose.th)) * r.G.L * cam.k * 0.5 + 3), 2.5, 0, 0, Math.PI * 2); ctx.fill();
  }
  const kAx = Math.max(cam.k, 30 / r.G.L), ex = Math.max(1, 3.2 / (r.G.R * kAx));
  drawExhaust(ctx, cam, r, t, pose, false, Math.max(0.6, r.p.nozzle / 2000 * kAx * ex));
  drawChute(ctx, cam, r, t, pose, 9);
  drawRocket(ctx, { X: cam.ox + pose.x * cam.k, Y: cam.oy - pose.y * cam.k, ang: -pose.th, k: kAx, ex, lod: 0, G: r.G, p: r.p, fs: pose.fs,
    w: { gs: pose.gs, gc: pose.gc, fc: pose.fc, mist: pose.mist, s: pose.s, ph: t * 14, dry: false, wd: pose.now.wd } });
  drawImpact(ctx, cam, r, t, pose);
  drawMarkers(ctx, W, H, cam, r, t);
  drawMainHUD(ctx, W, H, cam, r, t, pose);
  const vg = ctx.createRadialGradient(W / 2, H / 2, Math.min(W, H) * 0.4, W / 2, H / 2, Math.hypot(W, H) * 0.62);
  vg.addColorStop(0, 'rgba(2,8,16,0)'); vg.addColorStop(1, 'rgba(2,8,16,0.42)');
  ctx.fillStyle = vg; ctx.fillRect(0, 0, W, H);
}

// ------------------------------------------------------------------ close-up camera
function drawStreaks(ctx, W, H, cam, pose) {
  const n = pose.now, v = Math.hypot(n.vx, n.vy);
  if (pose.pad || pose.landed || v < 2) return;
  const k = cam.k, sh = 1 / 40, lx = n.vx * sh, ly = n.vy * sh, cell = 0.3;
  const xL = -cam.ox / k, xR = (W - cam.ox) / k, yT = cam.oy / k, yB = (cam.oy - H) / k;
  const x0 = Math.min(xL, xL - lx), x1 = Math.max(xR, xR - lx), y0 = Math.max(0, Math.min(yB, yB - ly)), y1 = Math.max(yT, yT - ly);
  const i0 = Math.floor(x0 / cell), i1 = Math.ceil(x1 / cell), j0 = Math.floor(y0 / cell), j1 = Math.ceil(y1 / cell);
  if ((i1 - i0) * (j1 - j0) > 900) return;
  ctx.strokeStyle = `rgba(255,255,255,${0.1 + 0.22 * clamp((v - 2) / 30, 0, 1)})`; ctx.lineWidth = 1; ctx.beginPath();
  for (let i = i0; i <= i1; i++) for (let j = j0; j <= j1; j++) {
    const h = hash(i * 31.7 + j * 17.3);
    if (h > 0.55) continue;
    const x = (i + hash(i * 3.1 + j * 7.9)) * cell, y = (j + hash(i * 5.3 + j * 2.7)) * cell;
    const X = cam.ox + x * k, Y = cam.oy - y * k;
    ctx.moveTo(X, Y); ctx.lineTo(X + lx * k, Y - ly * k);
  }
  ctx.stroke();
}
function drawCloseHUD(ctx, W, H, r, t, pose, X, Y, ang, k) {
  const f = F(), n = pose.now;
  const side = pose.landed || pose.pad ? 0 : Math.abs(n.sz || 0);
  pill(ctx, 10, 20, 'ÇALKALANMA %' + nf(pose.landed ? 0 : pose.sa * 100, 0) + (side > 0.15 ? ' · YANA %' + nf(side * 100, 0) : ''), { font: `500 10.5px ${f.m}`, color: '#bfe9ff', stroke: 'rgba(76,195,242,0.45)' });
  pill(ctx, 10, 45, 'SÖNÜM ζ ' + nf(r.zeta, 2) + ' · LEVHA %' + nf(r.p.baffle * 100, 0) + ' · ÇAP ×1,6', { font: `500 10.5px ${f.m}`, color: '#d7e4ef' });
  // labels for the two tanks (placed beside the body, nudged apart if they collide)
  const R = r.G.R * k * EX_CLOSE, parts = rocketParts(r.G, X, Y, ang, k * (pose.fs || 1));
  let px = -Math.sin(ang), py = Math.cos(ang);
  if (px < -0.15 || (Math.abs(px) <= 0.15 && py > 0)) { px = -px; py = -py; }
  const boxes = [];
  [['Çalkalanma tankı', parts.pay], ['Basınç haznesi', parts.ch]].forEach(([txt, P]) => {
    ctx.font = `600 11px ${f.b}`;
    const tw = ctx.measureText(txt).width + 14, align = px >= 0 ? 'left' : 'right';
    let off = R + 26, lx, ly, bx0;
    for (let tries = 0; tries < 4; tries++) {
      lx = P[0] + px * off; ly = P[1] + py * off;
      if (align === 'left') lx = clamp(lx, 6, W - 6 - tw); else lx = clamp(lx, 6 + tw, W - 6);
      bx0 = align === 'left' ? lx : lx - tw;
      if (!boxes.some((q) => bx0 < q[0] + q[2] + 4 && bx0 + tw + 4 > q[0] && Math.abs(ly - q[1]) < 24)) break;
      off += 26;
    }
    if (ly < 62 || ly > H - 14) return;
    boxes.push([bx0, ly, tw]);
    const ax = P[0] + px * (R + 3), ay = P[1] + py * (R + 3);
    ctx.strokeStyle = 'rgba(255,255,255,0.7)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(align === 'left' ? bx0 : bx0 + tw, ly); ctx.stroke();
    ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(ax, ay, 2, 0, Math.PI * 2); ctx.fill();
    pill(ctx, bx0, ly, txt, { font: `600 11px ${f.b}` });
  });
  // apparent-gravity gauge
  const gx = W - 38, gyy = 52, gr = 22;
  ctx.fillStyle = 'rgba(4,12,22,0.62)'; ctx.beginPath(); ctx.arc(gx, gyy, gr + 6, 0, Math.PI * 2); ctx.fill();
  ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(gx, gyy, gr, 0, Math.PI * 2); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(gx, gyy + gr - 4); ctx.lineTo(gx, gyy + gr + 2); ctx.stroke();
  const wa = pose.th + pose.gc, dx = Math.cos(wa), dy = -Math.sin(wa);
  ctx.strokeStyle = '#ffb547'; ctx.fillStyle = '#ffb547'; ctx.lineWidth = 2.2;
  ctx.beginPath(); ctx.moveTo(gx - dx * 8, gyy - dy * 8); ctx.lineTo(gx + dx * (gr - 6), gyy + dy * (gr - 6)); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(gx + dx * gr, gyy + dy * gr); ctx.lineTo(gx + dx * (gr - 9) - dy * 5, gyy + dy * (gr - 9) + dx * 5); ctx.lineTo(gx + dx * (gr - 9) + dy * 5, gyy + dy * (gr - 9) - dx * 5); ctx.closePath(); ctx.fill();
  const gm = pose.pad || pose.landed ? 1 : n.gm;
  pill(ctx, gx, gyy + gr + 17, nf(gm, gm < 10 ? 1 : 0) + ' g', { align: 'center', font: `500 10px ${f.m}`, h: 17, color: '#ffd59a' });
  pill(ctx, gx - gr - 10, gyy, 'GÖRÜNÜR g', { align: 'right', font: `500 9px ${f.m}`, h: 17, color: '#ffd59a' });
}
function drawCloseup(r, t, pose) {
  const W = cu._w, H = cu._h; if (!W) return;
  const ctx = cctx, cam = camClose(W, H, r, pose);
  ctx.clearRect(0, 0, W, H);
  drawBackdrop(ctx, W, H, cam, 0.15, pose.cy);
  drawStreaks(ctx, W, H, cam, pose);
  drawGround(ctx, W, H, cam, t, r, true, pose);
  const k = cam.k;
  drawExhaust(ctx, cam, r, t, pose, true, Math.max(0.8, r.p.nozzle / 2000 * k * EX_CLOSE));
  drawChute(ctx, cam, r, t, pose, 0);
  const X = cam.ox + pose.x * k, Y = cam.oy - pose.y * k, ang = -pose.th;
  const powered = !pose.pad && t < (r.tAirEnd || 0);
  drawRocket(ctx, { X, Y, ang, k, ex: EX_CLOSE, lod: 1, G: r.G, p: r.p, fs: pose.fs,
    w: { gs: pose.gs, gc: pose.gc, fc: pose.fc, mist: pose.mist, s: pose.s, ph: t * (powered ? 60 : 9), dry: !pose.pad && t > (r.tWaterEnd || 0), wd: pose.now.wd } });
  drawImpact(ctx, cam, r, t, pose);
  drawCloseHUD(ctx, W, H, r, t, pose, X, Y, ang, k);
  const vg = ctx.createRadialGradient(W / 2, H / 2, Math.min(W, H) * 0.42, W / 2, H / 2, Math.hypot(W, H) * 0.62);
  vg.addColorStop(0, 'rgba(2,8,16,0)'); vg.addColorStop(1, 'rgba(2,8,16,0.4)');
  ctx.fillStyle = vg; ctx.fillRect(0, 0, W, H);
}

// ================================================================== 3D view (three.js r128, loaded on demand)
const THREE_URLS = ['https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js', 'https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js'];
const CAM3 = {
  takip: { yaw: -18, pitch: 8, dist: 5.2, zmin: 0.25, zmax: 8 },
  yakin: { yaw: -24, pitch: 10, dist: 1.55, zmin: 0.4, zmax: 4 },
  genel: { yaw: -26, pitch: 15, dist: 0, zmin: 0.3, zmax: 3 },
};
const V3 = { state: 'idle', loading: null, cam: CAM3[store.get('ikarus-cam3', 'takip')] ? store.get('ikarus-cam3', 'takip') : 'takip', yaw: 0, pitch: 0, zoom: 1,
  built: null, builtRef: null, w: 0, h: 0 };
const gl3 = $('#gl'), hud3 = $('#hud3'), hctx = hud3.getContext('2d');
function msg3(text) { const m = $('#msg3'); m.hidden = !text; m.textContent = text || ''; }
function loadThree() {
  if (window.THREE) return Promise.resolve();
  if (V3.loading) return V3.loading;
  V3.state = 'loading'; msg3('3B sahne hazırlanıyor…');
  V3.loading = new Promise((res, rej) => {
    let i = 0;
    const next = () => {
      if (i >= THREE_URLS.length) { rej(new Error('three.js')); return; }
      const sc = document.createElement('script'); sc.src = THREE_URLS[i++]; sc.async = true;
      sc.onload = () => (window.THREE ? res() : next()); sc.onerror = next;
      document.head.appendChild(sc);
    };
    next();
  });
  return V3.loading;
}
function ensure3() {
  if (V3.state === 'ready') return Promise.resolve(true);
  if (V3.state === 'failed') return Promise.resolve(false);
  return loadThree().then(() => {
    if (V3.state !== 'ready') { if (!init3()) return false; }
    msg3(''); return true;
  }, () => {
    V3.state = 'failed';
    msg3('3B görünüm için gereken kütüphane yüklenemedi. İnternet bağlantınızı kontrol edip sayfayı yenileyin; bu sırada 2B görünüm çalışmaya devam eder.');
    return false;
  });
}

function canvasTex(w, h, draw) {
  const c = document.createElement('canvas'); c.width = w; c.height = h; draw(c.getContext('2d'), w, h);
  const t = new THREE.CanvasTexture(c); return t;
}
function init3() {
  const T = THREE;
  let renderer;
  try { renderer = new T.WebGLRenderer({ canvas: gl3, antialias: true, powerPreference: 'high-performance' }); }
  catch (e) { V3.state = 'failed'; msg3('Bu tarayıcıda WebGL çalışmıyor; 3B görünüm açılamadı. 2B görünüm kullanılabilir.'); return false; }
  renderer.shadowMap.enabled = true; renderer.shadowMap.type = T.PCFSoftShadowMap;
  renderer.localClippingEnabled = true;
  const scene = new T.Scene();
  const HOR = 0xd3e7f3;
  scene.fog = new T.Fog(HOR, 110, 700);
  const camera = new T.PerspectiveCamera(50, 1, 0.03, 2500);
  const E = { T, renderer, scene, camera };
  V3.E = E;
  // sky dome
  const sky = new T.Mesh(new T.SphereGeometry(1000, 32, 16), new T.ShaderMaterial({
    side: T.BackSide, depthWrite: false, fog: false,
    uniforms: { top: { value: new T.Color(0x1f5fb4) }, mid: { value: new T.Color(0x67a6df) }, hor: { value: new T.Color(HOR) }, low: { value: new T.Color(0x8fb38a) } },
    vertexShader: 'varying vec3 vP; void main(){ vP = normalize(position); gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); }',
    fragmentShader: 'uniform vec3 top; uniform vec3 mid; uniform vec3 hor; uniform vec3 low; varying vec3 vP;' +
      'void main(){ float h = vP.y; vec3 c = h > 0.0 ? mix(hor, mix(mid, top, smoothstep(0.12, 0.75, h)), smoothstep(0.0, 0.16, h)) : mix(hor, low, smoothstep(0.0, -0.06, h)); gl_FragColor = vec4(c, 1.0); }',
  }));
  sky.frustumCulled = false; sky.renderOrder = -10; scene.add(sky); E.sky = sky;
  E.sunDir = new T.Vector3(-0.42, 0.74, 0.52).normalize();
  const sunTex = canvasTex(128, 128, (g) => {
    const gr = g.createRadialGradient(64, 64, 0, 64, 64, 64);
    gr.addColorStop(0, 'rgba(255,252,240,1)'); gr.addColorStop(0.12, 'rgba(255,248,225,1)'); gr.addColorStop(0.3, 'rgba(255,236,190,0.35)'); gr.addColorStop(1, 'rgba(255,236,190,0)');
    g.fillStyle = gr; g.fillRect(0, 0, 128, 128);
  });
  E.sun = new T.Sprite(new T.SpriteMaterial({ map: sunTex, fog: false, depthWrite: false, transparent: true }));
  E.sun.scale.set(260, 260, 1); scene.add(E.sun);
  // lights
  scene.add(new T.HemisphereLight(0xd6eaff, 0x4d6b42, 0.68));
  const light = new T.DirectionalLight(0xfff1dc, 0.88);
  light.castShadow = true; light.shadow.mapSize.set(2048, 2048);
  Object.assign(light.shadow.camera, { left: -12, right: 12, top: 12, bottom: -12, near: 1, far: 160 });
  light.shadow.bias = -0.0004; if ('normalBias' in light.shadow) light.shadow.normalBias = 0.02;
  scene.add(light, light.target); E.light = light;
  // ground
  const grass = canvasTex(256, 256, (g) => {
    g.fillStyle = '#6fa456'; g.fillRect(0, 0, 256, 256);
    for (let i = 0; i < 2600; i++) {
      const x = hash(i * 1.3) * 256, y = hash(i * 2.7 + 9) * 256, v = hash(i * 5.1 + 3);
      g.fillStyle = v < 0.33 ? 'rgba(120,175,85,0.55)' : v < 0.66 ? 'rgba(60,110,45,0.5)' : 'rgba(150,190,100,0.35)';
      g.fillRect(x, y, 1 + 2 * hash(i), 1 + 3 * hash(i + 7));
    }
  });
  grass.wrapS = grass.wrapT = T.RepeatWrapping; grass.repeat.set(240, 240);
  if (renderer.capabilities.getMaxAnisotropy) grass.anisotropy = Math.min(8, renderer.capabilities.getMaxAnisotropy());
  const gg = new T.PlaneGeometry(1500, 1500, 150, 150), gpos = gg.attributes.position, gcol = new Float32Array(gpos.count * 3);
  for (let i = 0; i < gpos.count; i++) {
    const x = gpos.getX(i), y = gpos.getY(i);
    const nz = 0.5 + 0.22 * Math.sin(x * 0.045 + 1.3) * Math.cos(y * 0.038) + 0.14 * Math.sin(x * 0.11 + y * 0.07 + 2) + 0.08 * Math.sin(x * 0.27 - y * 0.19);
    const c = new T.Color().setHSL(0.25 + 0.03 * nz, 0.42 + 0.1 * nz, 0.3 + 0.12 * nz);
    gcol[i * 3] = c.r; gcol[i * 3 + 1] = c.g; gcol[i * 3 + 2] = c.b;
  }
  gg.setAttribute('color', new T.BufferAttribute(gcol, 3));
  const ground = new T.Mesh(gg, new T.MeshLambertMaterial({ map: grass, vertexColors: true }));
  ground.rotation.x = -Math.PI / 2; ground.receiveShadow = true; scene.add(ground);
  // distance rings and the downrange line
  [5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150].forEach((rr) => {
    const m = new T.Mesh(new T.RingGeometry(rr - 0.035, rr + 0.035, 180), new T.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: rr % 10 === 0 ? 0.5 : 0.28, depthWrite: false }));
    m.rotation.x = -Math.PI / 2; m.position.y = 0.012; scene.add(m);
  });
  const dl = new T.Line(new T.BufferGeometry().setFromPoints([new T.Vector3(0, 0.015, 0), new T.Vector3(60, 0.015, 0)]),
    new T.LineDashedMaterial({ color: 0xffffff, dashSize: 0.6, gapSize: 0.4, transparent: true, opacity: 0.55 }));
  dl.computeLineDistances(); scene.add(dl);
  // far hills and trees
  for (let i = 0; i < 18; i++) {
    const a = i / 18 * Math.PI * 2 + hash(i * 3.1) * 0.3, rr = 240 + 110 * hash(i * 7.3);
    const hm = new T.Mesh(new T.SphereGeometry(1, 12, 8), new T.MeshStandardMaterial({ color: new T.Color().setHSL(0.27 + 0.08 * hash(i), 0.22, 0.36 + 0.08 * hash(i * 2)), flatShading: true, roughness: 1 }));
    hm.scale.set(60 + 60 * hash(i * 1.7), 16 + 22 * hash(i * 4.4), 60 + 60 * hash(i * 9.2));
    hm.position.set(Math.cos(a) * rr, -4, Math.sin(a) * rr); scene.add(hm);
  }
  const NT = 220, cone = new T.ConeGeometry(1.1, 3, 7), trunk = new T.CylinderGeometry(0.12, 0.17, 1, 5);
  const crowns = new T.InstancedMesh(cone, new T.MeshStandardMaterial({ color: 0xffffff, flatShading: true, roughness: 0.9 }), NT);
  const trunks = new T.InstancedMesh(trunk, new T.MeshLambertMaterial({ color: 0x5b4632 }), NT);
  const mtx = new T.Matrix4(), q = new T.Quaternion(), sc = new T.Vector3(), ps = new T.Vector3(), col = new T.Color();
  for (let i = 0; i < NT; i++) {
    const a = hash(i * 12.9) * Math.PI * 2, rr = 55 + 170 * Math.pow(hash(i * 7.7), 1.2), s = 0.9 + 1.1 * hash(i * 3.3);
    ps.set(Math.cos(a) * rr, s * 2.4, Math.sin(a) * rr); sc.set(s, s, s); mtx.compose(ps, q, sc); crowns.setMatrixAt(i, mtx);
    if (crowns.setColorAt) crowns.setColorAt(i, col.setHSL(0.3 + 0.06 * hash(i), 0.38, 0.2 + 0.08 * hash(i * 2)));
    ps.set(Math.cos(a) * rr, s * 0.5, Math.sin(a) * rr); mtx.compose(ps, q, sc); trunks.setMatrixAt(i, mtx);
  }
  scene.add(crowns, trunks);
  // clouds
  SPR.clouds.forEach((cv, k) => { cv._tex = new T.CanvasTexture(cv); });
  for (let i = 0; i < 16; i++) {
    const spr = new T.Sprite(new T.SpriteMaterial({ map: SPR.clouds[i % 4]._tex, fog: false, transparent: true, depthWrite: false, opacity: 0.75 + 0.2 * hash(i) }));
    const a = hash(i * 4.4) * Math.PI * 2, rr = 120 + 380 * hash(i * 8.1), s = 90 + 120 * hash(i * 2.9);
    spr.position.set(Math.cos(a) * rr, 70 + 110 * hash(i * 6.6), Math.sin(a) * rr); spr.scale.set(s, s * 0.43, 1); scene.add(spr);
  }
  // launch site props
  const M = (c, o = {}) => new T.MeshStandardMaterial(Object.assign({ color: c, roughness: 0.6, metalness: 0.1 }, o));
  const add = (geo, mat, x, y, z, shadow = true) => { const m = new T.Mesh(geo, mat); m.position.set(x, y, z); m.castShadow = shadow; m.receiveShadow = true; scene.add(m); return m; };
  const leg = M(0x5b6a78, { metalness: 0.5, roughness: 0.4 });
  [0, 2.094, 4.189].forEach((a) => {
    const g = new T.CylinderGeometry(0.012, 0.012, 0.3, 6); const m = add(g, leg, Math.cos(a) * 0.12, PAD_H / 2, Math.sin(a) * 0.12);
    m.lookAt(0, PAD_H * 2, 0); m.rotateX(Math.PI / 2);
  });
  add(new T.CylinderGeometry(0.05, 0.05, 0.04, 16), M(0x3a4756), 0, PAD_H - 0.02, 0);
  add(new T.BoxGeometry(0.52, 0.3, 0.34), M(0xe1a43a, { roughness: 0.5 }), -3.7, 0.23, 1.1);
  const wheel = M(0x1d2630);
  [[-3.85, 1.29], [-3.55, 1.29], [-3.85, 0.91], [-3.55, 0.91]].forEach(([x, z]) => { const w = add(new T.CylinderGeometry(0.06, 0.06, 0.04, 14), wheel, x, 0.06, z); w.rotation.x = Math.PI / 2; });
  const hose = new T.CatmullRomCurve3([new T.Vector3(-3.45, 0.2, 1.1), new T.Vector3(-2.4, 0.02, 0.9), new T.Vector3(-1, 0.03, 0.3), new T.Vector3(-0.08, PAD_H - 0.03, 0)]);
  add(new T.TubeGeometry(hose, 50, 0.012, 6), M(0x151b22), 0, 0, 0);
  // safety line at 5 m
  const red = M(0xe63946), wht = M(0xf1f4f7);
  [-4, -2, 0, 2, 4].forEach((z) => { add(new T.CylinderGeometry(0.02, 0.02, 0.9, 8), wht, -5, 0.45, z); });
  add(new T.BoxGeometry(0.01, 0.06, 8), red, -5, 0.84, 0);
  // crew
  const person = (x, z, h, shirt) => {
    const body = add(new T.CylinderGeometry(0.17, 0.15, 0.62 * h / 1.75, 10), M(shirt), x, 1.1 * h / 1.75, z);
    add(new T.SphereGeometry(0.11, 12, 10), M(0xd9a77f), x, 1.6 * h / 1.75, z);
    [-0.08, 0.08].forEach((dz) => add(new T.CylinderGeometry(0.06, 0.06, 0.82 * h / 1.75, 8), M(0x1d2733), x, 0.41 * h / 1.75, z + dz));
    return body;
  };
  person(-5.6, -0.7, 1.76, 0xff7a2f); person(-6.1, 0.4, 1.66, 0x223246); person(-6.7, -0.2, 1.82, 0x1b2735);
  // wind sock
  add(new T.CylinderGeometry(0.02, 0.025, 2.2, 8), M(0xc9d3dc, { metalness: 0.5 }), -2.4, 1.1, -1.6);
  const stripes = canvasTex(64, 16, (g) => { for (let i = 0; i < 5; i++) { g.fillStyle = i % 2 ? '#f4f6f8' : '#ff7a2f'; g.fillRect(i * 64 / 5, 0, 64 / 5 + 1, 16); } });
  const sockGeo = new T.CylinderGeometry(0.09, 0.04, 0.7, 14, 1, true); sockGeo.translate(0, -0.35, 0);
  // remap so the stripes run along the sock
  const uv = sockGeo.attributes.uv; for (let i = 0; i < uv.count; i++) uv.setXY(i, 1 - uv.getY(i), uv.getX(i));
  E.sock = add(sockGeo, new T.MeshStandardMaterial({ map: stripes, side: T.DoubleSide, roughness: 0.8 }), -2.4, 2.18, -1.6);

  // rail (rebuilt with the launch angle)
  E.rail = add(new T.CylinderGeometry(0.006, 0.006, 1, 8), M(0xd6dee6, { metalness: 0.8, roughness: 0.25 }), 0, 0, 0);
  // particles
  const dot = canvasTex(32, 32, (g) => { const gr = g.createRadialGradient(16, 16, 0, 16, 16, 16); gr.addColorStop(0, 'rgba(255,255,255,1)'); gr.addColorStop(0.5, 'rgba(255,255,255,0.8)'); gr.addColorStop(1, 'rgba(255,255,255,0)'); g.fillStyle = gr; g.fillRect(0, 0, 32, 32); });
  const puffTex = new T.CanvasTexture(SPR.puff);
  const mkPts = (n, map, color, blend) => {
    const geo = new T.BufferGeometry();
    geo.setAttribute('position', new T.BufferAttribute(new Float32Array(n * 3), 3));
    geo.setAttribute('aSize', new T.BufferAttribute(new Float32Array(n), 1));
    geo.setAttribute('aAlpha', new T.BufferAttribute(new Float32Array(n), 1));
    const mat = new T.ShaderMaterial({
      transparent: true, depthWrite: false, blending: blend || T.NormalBlending,
      uniforms: { map: { value: map }, color: { value: new T.Color(color) }, uScale: { value: 400 } },
      vertexShader: 'attribute float aSize; attribute float aAlpha; varying float vA; uniform float uScale;' +
        'void main(){ vA = aAlpha; vec4 mv = modelViewMatrix * vec4(position, 1.0); gl_PointSize = clamp(aSize * uScale / -mv.z, 2.0, 600.0); gl_Position = projectionMatrix * mv; }',
      fragmentShader: 'uniform sampler2D map; uniform vec3 color; varying float vA;' +
        'void main(){ vec4 t = texture2D(map, gl_PointCoord); float a = t.a * vA; if (a < 0.01) discard; gl_FragColor = vec4(color * mix(0.85, 1.0, t.r), a); }',
    });
    const pts = new T.Points(geo, mat); pts.frustumCulled = false; scene.add(pts);
    return { pts, geo, mat, n };
  };
  E.puffs = mkPts(900, puffTex, 0xf4f8fc);
  E.drops = mkPts(1200, dot, 0xe2f3ff);
  E.dust = mkPts(40, dot, 0x6e5a3e);
  // paths and markers (geometry filled in rebuild3)
  E.pathGroup = new T.Group(); scene.add(E.pathGroup);
  // parachute
  const gores = canvasTex(128, 16, (g) => { for (let i = 0; i < 8; i++) { g.fillStyle = i % 2 ? '#f3f6f9' : '#ff7a2f'; g.fillRect(i * 16, 0, 16, 16); } });
  E.canopy = new T.Mesh(new T.SphereGeometry(1, 24, 10, 0, Math.PI * 2, 0, Math.PI / 2), new T.MeshStandardMaterial({ map: gores, side: T.DoubleSide, roughness: 0.7 }));
  E.canopy.castShadow = true; scene.add(E.canopy);
  const lg = new T.BufferGeometry(); lg.setAttribute('position', new T.BufferAttribute(new Float32Array(8 * 2 * 3), 3));
  E.lines = new T.LineSegments(lg, new T.LineBasicMaterial({ color: 0xe8eef4 })); E.lines.frustumCulled = false; scene.add(E.lines);
  // hud interaction
  wire3();
  V3.state = 'ready';
  return true;
}

// rocket model with see-through tanks and clipped water volumes
function lathe(T, pts, seg) { return new T.LatheGeometry(pts.map((p) => new T.Vector2(p[0], p[1])), seg || 36); }
function bottle(T, y0, y1, R) {
  const rc = Math.min(R * 0.45, (y1 - y0) * 0.22), pts = [[0.0005, y0]];
  for (let i = 0; i <= 6; i++) { const a = -Math.PI / 2 + i / 6 * Math.PI / 2; pts.push([R - rc + rc * Math.cos(a), y0 + rc + rc * Math.sin(a)]); }
  for (let i = 0; i <= 6; i++) { const a = i / 6 * Math.PI / 2; pts.push([R - rc + rc * Math.cos(a), y1 - rc + rc * Math.sin(a)]); }
  pts.push([0.0005, y1]);
  return lathe(T, pts);
}
function rebuild3(r) {
  const E = V3.E, T = E.T, G = r.G, p = r.p, scene = E.scene;
  const drop = (o) => { if (o.geometry) o.geometry.dispose(); if (o.material) { if (o.material.alphaMap) o.material.alphaMap.dispose(); o.material.dispose(); } };
  if (E.rocket) { scene.remove(E.rocket); E.rocket.traverse(drop); }
  [E.surfP, E.surfC, E.wingG, E.tailG].forEach((m) => { if (m) { scene.remove(m); m.traverse(drop); } });
  E.wingG = E.tailG = null;
  const R = G.R, Rin = R - 0.003, Lt = G.Ltail, Lc = G.Lc, Lp = G.Lp, Ln = G.Ln, L = G.L;
  const uc0 = Lt, uc1 = Lt + Lc, up0 = uc1, up1 = uc1 + Lp, un0 = up1;
  const g = new T.Group(); g.matrixAutoUpdate = true;
  const M = (c, o = {}) => new T.MeshStandardMaterial(Object.assign({ color: c, roughness: 0.45, metalness: 0.1 }, o));
  const put = (geo, mat, y, shadow = true) => { const m = new T.Mesh(geo, mat); m.position.y = y; m.castShadow = shadow; g.add(m); return m; };
  const rN = p.nozzle / 2000;
  put(new T.CylinderGeometry(Math.min(R * 0.9, Math.max(0.0135, rN + 0.002)), Math.min(R * 0.85, rN), NOZ_L, 20), M(0x9aa8b6, { metalness: 0.8, roughness: 0.3 }), -NOZ_L / 2);
  put(new T.CylinderGeometry(R, R, Lt, 32), M(0x25303c, { metalness: 0.3 }), Lt / 2);
  put(new T.CylinderGeometry(R * 1.05, R * 1.05, 0.026, 32), M(0x34424f, { metalness: 0.4 }), uc1);
  // nose cone (ogive) and altimeter window
  const np = [];
  for (let i = 0; i <= 18; i++) { const v = i / 18; np.push([Math.max(0.0005, R * Math.pow(Math.cos(v * Math.PI / 2), 0.7)), un0 + v * Ln]); }
  put(lathe(T, np), M(0xf2f5f8, { roughness: 0.35 }), 0);
  const win = put(new T.BoxGeometry(0.024, 0.03, 0.01), M(0x0d1b28, { roughness: 0.2 }), un0 + Ln * 0.28, false);
  win.position.z = R * 0.9;
  // fins
  const fin = G.fin, sh = new T.Shape();
  sh.moveTo(0, 0); sh.lineTo(fin.chord, 0); sh.lineTo(0.45 * fin.chord, fin.span); sh.lineTo(-0.012, fin.span); sh.lineTo(0, 0);
  const fg = new T.ExtrudeGeometry(sh, { depth: 0.003, bevelEnabled: false }); fg.translate(0, 0, -0.0015);
  const fm = M(0xff7a2f, { roughness: 0.5 });
  [0, 1, 2, 3].forEach((k) => {
    const a = k * Math.PI / 2, d = new T.Vector3(Math.cos(a), 0, Math.sin(a));
    const ax = new T.Vector3(0, 1, 0), zz = new T.Vector3().crossVectors(ax, d);
    const m = new T.Mesh(fg, fm); m.matrixAutoUpdate = false;
    m.matrix.makeBasis(ax, d, zz).setPosition(d.x * R * 0.98, 0, d.z * R * 0.98);
    m.castShadow = true; g.add(m);
  });
  // water: closed cylinders clipped by a plane that follows the apparent gravity
  E.clipC = new T.Plane(); E.clipP = new T.Plane();
  E.capC0 = new T.Plane(); E.capC1 = new T.Plane(); E.capP0 = new T.Plane(); E.capP1 = new T.Plane();
  const wMat = (pl) => new T.MeshStandardMaterial({ color: 0x1f86d6, roughness: 0.12, metalness: 0.05, transparent: true, opacity: 0.86, side: T.DoubleSide, clippingPlanes: [pl] });
  const sMat = (a, b) => new T.MeshStandardMaterial({ color: 0x86d8ff, roughness: 0.08, metalness: 0.05, transparent: true, opacity: 0.92, side: T.DoubleSide, clippingPlanes: [a, b] });
  E.waterC = put(new T.CylinderGeometry(Rin, Rin, Lc - 0.008, 36), wMat(E.clipC), (uc0 + uc1) / 2, false);
  E.waterP = put(new T.CylinderGeometry(Rin, Rin, Lp - 0.008, 36), wMat(E.clipP), (up0 + up1) / 2, false);
  E.waterC.renderOrder = 1; E.waterP.renderOrder = 1;
  E.surfC = new T.Mesh(new T.CircleGeometry(1, 56), sMat(E.capC0, E.capC1)); E.surfC.matrixAutoUpdate = false; E.surfC.renderOrder = 2; scene.add(E.surfC);
  E.surfP = new T.Mesh(new T.CircleGeometry(1, 56), sMat(E.capP0, E.capP1)); E.surfP.matrixAutoUpdate = false; E.surfP.renderOrder = 2; scene.add(E.surfP);
  E.mist = put(new T.CylinderGeometry(Rin * 0.98, Rin * 0.98, Lc - 0.01, 28), new T.MeshBasicMaterial({ color: 0xf2f7fb, transparent: true, opacity: 0, depthWrite: false }), (uc0 + uc1) / 2, false);
  E.mist.renderOrder = 3;
  // baffle rings
  if (p.baffle > 0) {
    const bm = new T.MeshStandardMaterial({ color: 0xeef4fa, roughness: 0.4, transparent: true, opacity: 0.9, side: T.DoubleSide });
    [1 / 3, 2 / 3].forEach((f) => { const rg = new T.RingGeometry(Rin * Math.sqrt(Math.max(0.01, 1 - p.baffle)), Rin, 40); rg.rotateX(-Math.PI / 2); const m = put(rg, bm, up0 + Lp * f, false); m.renderOrder = 2; });
  }
  // PET bottles (drawn last, see-through)
  const glass = new T.MeshPhysicalMaterial({ color: 0xdcefff, roughness: 0.05, metalness: 0, transparent: true, opacity: 0.24, depthWrite: false, clearcoat: 1, clearcoatRoughness: 0.1 });
  const gc = put(bottle(T, uc0, uc1, R), glass, 0); gc.renderOrder = 4;
  const gp = put(bottle(T, up0, up1, R), glass, 0); gp.renderOrder = 4;
  // water jet
  const jetAlpha = canvasTex(4, 64, (c) => { const gr = c.createLinearGradient(0, 0, 0, 64); gr.addColorStop(0, '#fff'); gr.addColorStop(0.45, '#999'); gr.addColorStop(1, '#000'); c.fillStyle = gr; c.fillRect(0, 0, 4, 64); });
  const jg = new T.CylinderGeometry(0.22, 1, 1, 20, 1, true); jg.translate(0, -0.5, 0);
  E.jet = put(jg, new T.MeshBasicMaterial({ color: 0xeaf7ff, alphaMap: jetAlpha, transparent: true, opacity: 0.85, depthWrite: false, side: T.DoubleSide }), -NOZ_L, false);
  E.jet.renderOrder = 5; E.jetR = rN;
  E.rocket = g; scene.add(g);
  if (E.ltube) { scene.remove(E.ltube); drop(E.ltube); E.ltube = null; }
  if (p.tube) {
    const tg = new T.CylinderGeometry((p.nozzle - 2) / 2000, (p.nozzle - 2) / 2000, 0.30, 14); tg.translate(0, 0.15, 0);
    E.ltube = new T.Mesh(tg, M(0xb8c4cf, { metalness: 0.6, roughness: 0.35 }));
    E.ltube.position.set(0, PAD_H - NOZ_L, 0);
    E.ltube.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), new T.Vector3(Math.cos(r.th0), Math.sin(r.th0), 0));
    scene.add(E.ltube);
  }
  E.dims = { R, Rin, uc0, uc1, up0, up1, Lc, Lp };
  // wing panels hinge at the body sides; the group frame is (span, body axis, wing-up)
  if (G.wing) {
    const half = p.wspan / 2, c = G.chord, wm = M(0xffb04a, { roughness: 0.55 });
    const wg = new T.Group(); wg.matrixAutoUpdate = false;
    E.wingPiv = [1, -1].map((sd) => {
      const piv = new T.Group(); piv.position.set(sd * R, 0, 0);
      const inc = new T.Group(); inc.rotation.x = G.winc;
      const geo = new T.BoxGeometry(half, c, 0.008); geo.translate(sd * half / 2, -0.25 * c, 0);
      const m = new T.Mesh(geo, wm); m.castShadow = true;
      inc.add(m); piv.add(inc); wg.add(piv); piv.userData.sd = sd; return piv;
    });
    E.wingG = wg; scene.add(wg);
    const bt = Math.sqrt(4 * G.St), ct = G.St / bt, tg = new T.Group(); tg.matrixAutoUpdate = false;
    const tm = new T.Mesh(new T.BoxGeometry(bt, ct, 0.005), M(0xe8eef4, { roughness: 0.5 })); tm.castShadow = true; tg.add(tm);
    E.tailG = tg; E.tailU = Math.max(ct / 2, L - G.xT); scene.add(tg);
  }
  // launch rail beside the rocket
  const th0 = r.th0, n0 = new T.Vector3(-Math.sin(th0), Math.cos(th0), 0), b0 = new T.Vector3(Math.cos(th0), Math.sin(th0), 0);
  E.rail.position.set(0, PAD_H, 0).addScaledVector(n0, R + 0.014).addScaledVector(b0, 0.5);
  E.rail.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), b0);
  buildPaths3(r);
  V3.built = r;
}
function buildPaths3(r) {
  const E = V3.E, T = E.T, grp = E.pathGroup;
  while (grp.children.length) { const o = grp.children.pop(); if (o.geometry) o.geometry.dispose(); if (o.material) o.material.dispose(); }
  const pts = r.rec.map((e) => new T.Vector3(e.x, e.y + PAD_H, e.z));
  if (pts.length < 2) return;
  // flown path: tube with one ring per record, coloured by phase; revealed with a draw range
  const n = pts.length, rs = 6, rad = 0.022;
  const pos = new Float32Array(n * rs * 3), colr = new Float32Array(n * rs * 3), idx = [];
  const PH = [new T.Color(0xc7d3de), new T.Color(0x62d2ff), new T.Color(0xffc061), new T.Color(0xffffff)];
  const up = new T.Vector3(0, 1, 0), tng = new T.Vector3(), a = new T.Vector3(), bb = new T.Vector3();
  for (let i = 0; i < n; i++) {
    tng.subVectors(pts[Math.min(n - 1, i + 1)], pts[Math.max(0, i - 1)]); if (tng.lengthSq() < 1e-12) tng.set(1, 0, 0); tng.normalize();
    a.crossVectors(tng, Math.abs(tng.y) > 0.95 ? new T.Vector3(1, 0, 0) : up).normalize(); bb.crossVectors(tng, a).normalize();
    const c = PH[r.rec[i].chute ? 3 : r.rec[i].ph] || PH[3];
    for (let j = 0; j < rs; j++) {
      const ang = j / rs * Math.PI * 2, k = (i * rs + j) * 3;
      pos[k] = pts[i].x + (a.x * Math.cos(ang) + bb.x * Math.sin(ang)) * rad;
      pos[k + 1] = pts[i].y + (a.y * Math.cos(ang) + bb.y * Math.sin(ang)) * rad;
      pos[k + 2] = pts[i].z + (a.z * Math.cos(ang) + bb.z * Math.sin(ang)) * rad;
      colr[k] = c.r; colr[k + 1] = c.g; colr[k + 2] = c.b;
    }
    if (i < n - 1) for (let j = 0; j < rs; j++) {
      const a0 = i * rs + j, a1 = i * rs + (j + 1) % rs, b0 = a0 + rs, b1 = a1 + rs;
      idx.push(a0, b0, a1, a1, b0, b1);
    }
  }
  const tg = new T.BufferGeometry();
  tg.setAttribute('position', new T.BufferAttribute(pos, 3)); tg.setAttribute('color', new T.BufferAttribute(colr, 3)); tg.setIndex(idx);
  E.tube = new T.Mesh(tg, new T.MeshBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.8 })); E.tube.frustumCulled = false; grp.add(E.tube);
  const lc = new Float32Array(n * 3); for (let i = 0; i < n; i++) { const c = PH[r.rec[i].chute ? 3 : r.rec[i].ph] || PH[3]; lc[i * 3] = c.r; lc[i * 3 + 1] = c.g; lc[i * 3 + 2] = c.b; }
  const lg = new T.BufferGeometry().setFromPoints(pts); lg.setAttribute('color', new T.BufferAttribute(lc, 3));
  E.trail = new T.Line(lg, new T.LineBasicMaterial({ vertexColors: true })); E.trail.frustumCulled = false; grp.add(E.trail);
  const dash = (points, color, op) => { const l = new T.Line(new T.BufferGeometry().setFromPoints(points), new T.LineDashedMaterial({ color, dashSize: 0.45, gapSize: 0.3, transparent: true, opacity: op })); l.computeLineDistances(); grp.add(l); return l; };
  dash(pts, 0xffffff, 0.32);
  const ref = S.ref;
  if (ref && ref !== r && ref.rec.length && (Math.abs(ref.apo - r.apo) > 1e-6 || Math.abs(ref.dist - r.dist) > 1e-6)) dash(ref.rec.map((e) => new T.Vector3(e.x, e.y + PAD_H, e.z)), 0x9fb0c2, 0.8);
  V3.builtRef = ref;
  // apogee and landing markers
  if (r.left) {
    const apex = new T.Vector3(r.xApo, r.apo + PAD_H, r.zApo);
    E.apoLine = dash([apex, new T.Vector3(r.xApo, 0.02, r.zApo)], 0xffffff, 0.6);
    E.apoDot = new T.Mesh(new T.OctahedronGeometry(0.12), new T.MeshBasicMaterial({ color: 0xffffff })); E.apoDot.position.copy(apex); grp.add(E.apoDot);
    const land = new T.Vector3(r.xLand, 0.02, r.zLand);
    E.landLine = new T.Line(new T.BufferGeometry().setFromPoints([new T.Vector3(0, 0.03, 0), land]), new T.LineDashedMaterial({ color: 0xffd25e, dashSize: 0.5, gapSize: 0.3 }));
    E.landLine.computeLineDistances(); grp.add(E.landLine);
    E.landRing = new T.Mesh(new T.RingGeometry(0.3, 0.42, 40), new T.MeshBasicMaterial({ color: 0xffd25e, transparent: true, opacity: 0.9, depthWrite: false }));
    E.landRing.rotation.x = -Math.PI / 2; E.landRing.position.copy(land); grp.add(E.landRing);
  } else { E.apoLine = E.apoDot = E.landLine = E.landRing = null; }
}

// ------------------------------------------------------------------ 3D frame
function size3() {
  if (V3.state !== 'ready') return;
  const st = $('#stage3'), w = st.clientWidth; if (!w) return;
  const h = Math.round(clamp(w * 0.56, 320, 620));
  if (w === V3.w && h === V3.h) return;
  V3.w = w; V3.h = h;
  const E = V3.E;
  E.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
  E.renderer.setSize(w, h, false); gl3.style.height = h + 'px';
  E.camera.aspect = w / h; E.camera.updateProjectionMatrix();
  setCanvas(hud3, hctx, w, h);
}
function setPlaneThrough(pl, n, p) { pl.normal.set(n[0], n[1], n[2]); pl.constant = -(n[0] * p[0] + n[1] * p[1] + n[2] * p[2]); }
function placeSurface(mesh, n, b, cen, Rin, frac, len) {
  // ellipse where the water plane cuts the tank: centred on the axis, semi-axes Rin and Rin/|n·b|
  const T = V3.E.T, nb = v3dot(n, b), cb = Math.max(0.04, Math.abs(nb));
  const h = Math.abs(nb) * len / 2 + Math.sqrt(Math.max(0, 1 - nb * nb)) * Rin, dlt = (1 - 2 * frac) * h;
  const u = clamp(dlt / (Math.abs(nb) < 0.04 ? (nb < 0 ? -0.04 : 0.04) : nb), -len, len);
  const c = [cen[0] + b[0] * u, cen[1] + b[1] * u, cen[2] + b[2] * u];
  let mj = [b[0] - nb * n[0], b[1] - nb * n[1], b[2] - nb * n[2]];
  if (Math.hypot(mj[0], mj[1], mj[2]) < 1e-6) mj = Math.abs(n[0]) < 0.9 ? v3cross(n, [1, 0, 0]) : v3cross(n, [0, 0, 1]);
  mj = v3n(mj); const mn = v3cross(n, mj), a = Rin / cb;
  mesh.matrix.set(mj[0] * a, mn[0] * Rin, n[0], c[0], mj[1] * a, mn[1] * Rin, n[1], c[1], mj[2] * a, mn[2] * Rin, n[2], c[2], 0, 0, 0, 1);
  mesh.matrixWorldNeedsUpdate = true;
  return dlt;
}
function particles3(E, r, t, pose) {
  const em = r.emit || [], wx = r.p.wind, wz = r.p.windZ || 0, P1 = r.P1g || 1, tW = r.tWaterEnd || 0;
  const fill = (sys, list) => {
    const pa = sys.geo.attributes.position.array, sa = sys.geo.attributes.aSize.array, aa = sys.geo.attributes.aAlpha.array;
    const n = Math.min(list.length, sys.n);
    for (let i = 0; i < n; i++) { const q = list[i]; pa[i * 3] = q[0]; pa[i * 3 + 1] = q[1]; pa[i * 3 + 2] = q[2]; sa[i] = q[3]; aa[i] = q[4]; }
    sys.geo.setDrawRange(0, n);
    sys.geo.attributes.position.needsUpdate = true; sys.geo.attributes.aSize.needsUpdate = true; sys.geo.attributes.aAlpha.needsUpdate = true;
  };
  const drops = [], puffs = [], dust = [];
  if (t > 0) for (const s of em) {
    const a = t - s.tau; if (a <= 0) continue;
    const b = s.b3, e1 = v3n(Math.abs(b[1]) < 0.9 ? v3cross(b, [0, 1, 0]) : v3cross(b, [1, 0, 0])), e2 = v3cross(b, e1);
    if (s.ph === 2) {
      if (a > 2) continue;
      const burst = s.tau - tW < 0.005 ? 6 : 3;
      for (let j = 0; j < burst; j++) {
        const h1 = hash(s.i * 3.1 + j * 17.7), h2 = hash(s.i * 7.9 + j * 5.3 + 2), h3 = hash(s.i * 1.3 + j * 9.1 + 4), h5 = hash(s.i * 2.3 + j * 3.3 + 6);
        const sp1 = (h1 - 0.5) * 0.8, sp2 = (h5 - 0.5) * 0.8;
        const d = v3n([-b[0] + e1[0] * sp1 + e2[0] * sp2, -b[1] + e1[1] * sp1 + e2[1] * sp2, -b[2] + e1[2] * sp1 + e2[2] * sp2]);
        const spd = (5 + 16 * h2) * (0.35 + 0.65 * clamp(s.P / P1, 0, 1)), kd = 3.2, E_ = (1 - Math.exp(-kd * a)) / kd;
        const v0 = [s.v3[0] * 0.15 + d[0] * spd, s.v3[1] * 0.15 + d[1] * spd, s.v3[2] * 0.15 + d[2] * spd], vi = [wx * 0.8, 0.25, wz * 0.8];
        puffs.push([s.p3[0] + vi[0] * a + (v0[0] - vi[0]) * E_, Math.max(0.05, s.p3[1] + vi[1] * a + (v0[1] - vi[1]) * E_), s.p3[2] + vi[2] * a + (v0[2] - vi[2]) * E_,
          2 * (0.05 + 0.5 * (1 - Math.exp(-1.8 * a)) * (0.6 + 0.7 * h3)), 0.34 * Math.pow(1 - a / 2, 1.6)]);
      }
    } else if (s.ph === 1) {
      if (a > 1.5) continue;
      const ve = r.G.cd * Math.sqrt(200 * Math.max(0, s.P));
      for (let j = 0; j < 5; j++) {
        const h1 = hash(s.i * 13.1 + j * 7.7), h2 = hash(s.i * 3.7 + j * 19.3 + 1), h3 = hash(s.i * 5.3 + j * 2.9 + 2), h4 = hash(s.i * 1.9 + j * 11.3 + 3), h5 = hash(s.i * 8.3 + j * 4.1 + 5);
        const sp1 = (h1 - 0.5) * 0.24, sp2 = (h5 - 0.5) * 0.24;
        const d = v3n([-b[0] + e1[0] * sp1 + e2[0] * sp2, -b[1] + e1[1] * sp1 + e2[1] * sp2, -b[2] + e1[2] * sp1 + e2[2] * sp2]);
        const spd = ve * (0.6 + 0.45 * h2), kd = 1.8 + 3.5 * h3;
        const v0 = [s.v3[0] + d[0] * spd, s.v3[1] + d[1] * spd, s.v3[2] + d[2] * spd], vi = [wx, -C.g / kd, wz];
        const at = (aa) => { const E_ = (1 - Math.exp(-kd * aa)) / kd; return [s.p3[0] + vi[0] * aa + (v0[0] - vi[0]) * E_, s.p3[1] + vi[1] * aa + (v0[1] - vi[1]) * E_, s.p3[2] + vi[2] * aa + (v0[2] - vi[2]) * E_]; };
        let q = at(a), al = 0.9 * (1 - a / 1.5);
        if (q[1] < 0) {
          let lo = 0, hi = a; for (let it = 0; it < 14; it++) { const m = (lo + hi) / 2; if (at(m)[1] > 0) lo = m; else hi = m; }
          const after = a - lo; if (after > 0.5) continue;
          const ex = Math.exp(-kd * lo), vy = Math.abs(vi[1] + (v0[1] - vi[1]) * ex), qi = at(lo), spread = 1 - Math.exp(-9 * after);
          const ang = h1 * Math.PI * 2, rr = (0.25 + 0.9 * h4) * spread;
          q = [qi[0] + Math.cos(ang) * rr, vy * 0.012 * (0.3 + h2) * Math.sin(Math.PI * Math.min(1, after / 0.32)), qi[2] + Math.sin(ang) * rr];
          al = 0.75 * (1 - after / 0.5);
          if (j < 2) puffs.push([qi[0] + Math.cos(ang) * rr * 0.6, 0.08 + 0.5 * after * (0.5 + h2), qi[2] + Math.sin(ang) * rr * 0.6, 0.25 + 0.9 * after, 0.3 * (1 - after / 0.5)]);
        }
        drops.push([q[0], q[1], q[2], 2 * (0.005 + 0.012 * h4), al]);
      }
    }
  }
  if (pose.landed && pose.contact) {
    const a = t - r.time;
    if (a < 1.3) for (let j = 0; j < 30; j++) {
      const h1 = hash(j * 4.3 + 7), h2 = hash(j * 9.7 + 1), h3 = hash(j * 2.1 + 5), ang = h1 * Math.PI * 2, sp = 1 + 2.5 * hash(j * 6.1);
      const kd = 2.4, E_ = (1 - Math.exp(-kd * a)) / kd, vy0 = 0.8 + 2.8 * h2, viy = -C.g / kd;
      dust.push([pose.contact[0] + Math.cos(ang) * sp * E_, Math.max(0, pose.contact[1] + viy * a + (vy0 - viy) * E_), pose.contact[2] + Math.sin(ang) * sp * E_, 0.03 + 0.03 * h3, 0.9 * (1 - a / 1.3)]);
    }
  }
  fill(E.drops, drops); fill(E.puffs, puffs); fill(E.dust, dust);
}
function draw3D(r, t) {
  if (V3.state !== 'ready' || !V3.w) return;
  const E = V3.E, T = E.T;
  if (V3.built !== r || V3.builtRef !== S.ref) rebuild3(r);
  const G = r.G, pose = pose3At(r, t), b = pose.b, D = E.dims;
  // rocket
  E.rocket.position.set(pose.tail[0], pose.tail[1], pose.tail[2]);
  E.rocket.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), new T.Vector3(b[0], b[1], b[2]));
  E.rocket.updateMatrixWorld(true);
  // payload tank water (always half full)
  const cP = pose.at((D.up0 + D.up1) / 2), cC = pose.at((D.uc0 + D.uc1) / 2);
  setPlaneThrough(E.clipP, pose.gS, cP);
  placeSurface(E.surfP, pose.gS, b, cP, D.Rin, 0.5, D.Lp - 0.008);
  const negB = [-b[0], -b[1], -b[2]];
  setPlaneThrough(E.capP0, b, pose.at(D.up0 + 0.004)); setPlaneThrough(E.capP1, negB, pose.at(D.up1 - 0.004));
  // chamber water
  const fc = pose.fc;
  E.waterC.visible = E.surfC.visible = fc > 0.003;
  if (fc > 0.003) {
    const dlt = placeSurface(E.surfC, pose.gC, b, cC, D.Rin, fc, D.Lc - 0.008);
    setPlaneThrough(E.clipC, pose.gC, [cC[0] + pose.gC[0] * dlt, cC[1] + pose.gC[1] * dlt, cC[2] + pose.gC[2] * dlt]);
    setPlaneThrough(E.capC0, b, pose.at(D.uc0 + 0.004)); setPlaneThrough(E.capC1, negB, pose.at(D.uc1 - 0.004));
    E.surfC.visible = fc < 0.997;
  }
  E.mist.material.opacity = pose.mist * 0.8; E.mist.visible = pose.mist > 0.01;
  if (E.wingG) {
    let nu = [-b[1] * b[0], 1 - b[1] * b[1], -b[1] * b[2]];
    if (Math.hypot(nu[0], nu[1], nu[2]) < 1e-3) nu = E.lastNu || [-Math.sin(r.th0), Math.cos(r.th0), 0];
    nu = v3n(nu); E.lastNu = nu;
    const lat = v3cross(b, nu), place = (grp, u) => {
      const q = pose.at(u);
      grp.matrix.set(lat[0], b[0], nu[0], q[0], lat[1], b[1], nu[1], q[1], lat[2], b[2], nu[2], q[2], 0, 0, 0, 1);
      grp.matrixWorldNeedsUpdate = true;
    };
    place(E.wingG, G.L - G.xW); place(E.tailG, E.tailU);
    const wd = pose.pad ? (G.wing === 'sabit' ? 1 : 0) : (pose.now.wd == null ? 1 : pose.now.wd);
    const pack = Math.min(1, 0.4 * G.L / (r.p.wspan / 2));
    E.wingPiv.forEach((pv) => { pv.rotation.z = -pv.userData.sd * (1 - wd) * Math.PI / 2; pv.scale.x = lerp(pack, 1, wd); });
  }
  // jet at the nozzle
  const n = pose.now, powered = !pose.pad && t < (r.tAirEnd || 0) && n.ph >= 1;
  E.jet.visible = powered;
  if (powered) {
    const water = n.ph === 1, P = Math.max(0, n.P);
    const rj = (E.jetR || 0.011) / 0.22, len = water ? Math.min(r.G.cd * Math.sqrt(200 * P) * 0.045, 1.6) : 0.5 * clamp(P / r.P1g, 0.1, 1), rr = water ? rj : rj * 2;
    E.jet.scale.set(rr, len, rr); E.jet.material.opacity = water ? 0.85 : 0.5 * clamp(P / r.P1g, 0.15, 1);
    E.jet.material.color.setHex(water ? 0xeaf7ff : 0xffffff);
  }
  // parachute
  const showChute = r.tChute != null && t >= r.tChute;
  E.canopy.visible = E.lines.visible = showChute;
  if (showChute) {
    const e = smooth((t - r.tChute) / 0.5);
    let cx = pose.nose[0] + 0.12 * Math.sin(1.7 * t), cy = pose.nose[1] + 1.6, cz = pose.nose[2] + 0.1 * Math.cos(1.3 * t), flat = 1;
    if (pose.landed) {
      const f = smooth((t - r.time) / 0.9), wx = r.p.wind || 0.01, wz = r.p.windZ || 0, wl = Math.hypot(wx, wz);
      cx = lerp(cx, pose.nose[0] + 0.9 * wx / wl, f); cy = lerp(cy, 0.06, f); cz = lerp(cz, pose.nose[2] + 0.9 * wz / wl, f); flat = 1 - 0.85 * f;
    }
    const Rc = 0.55 * (0.3 + 0.7 * e);
    E.canopy.position.set(cx, cy, cz); E.canopy.scale.set(Rc, Rc * 0.62 * flat, Rc);
    const la = E.lines.geometry.attributes.position.array;
    for (let i = 0; i < 8; i++) {
      const a = i / 8 * Math.PI * 2;
      la.set([cx + Math.cos(a) * Rc, cy, cz + Math.sin(a) * Rc, pose.nose[0], pose.nose[1], pose.nose[2]], i * 6);
    }
    E.lines.geometry.attributes.position.needsUpdate = true;
  }
  // paths
  let iCur = 0; const rc = r.rec, tt = Math.min(t, r.time);
  { let lo = 0, hi = rc.length - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (rc[m].t <= tt) lo = m; else hi = m; } iCur = t <= 0 ? 0 : (rc[hi] && rc[hi].t <= tt ? hi : lo); }
  if (E.tube) { E.tube.geometry.setDrawRange(0, Math.max(0, iCur) * 36); E.trail.geometry.setDrawRange(0, iCur + 1); E.trail.visible = iCur > 0; }
  if (E.apoDot) { const on = t >= r.tApo; E.apoDot.visible = on; E.apoLine.visible = on; const ld = t >= r.time; E.landLine.visible = ld; E.landRing.visible = ld; }
  // wind sock
  {
    const wx = r.p.wind, wz = r.p.windZ || 0, wl = Math.hypot(wx, wz), m = Math.min(1, wl / 7);
    const droop = (1 - m) * 1.25 + 0.05 * Math.sin(9 * t + 1) * m, hd = wl > 0.01 ? Math.atan2(wz, wx) : 0;
    const dir = new T.Vector3(Math.cos(hd) * Math.cos(droop), -Math.sin(droop), Math.sin(hd) * Math.cos(droop));
    E.sock.quaternion.setFromUnitVectors(new T.Vector3(0, -1, 0), dir);
  }
  particles3(E, r, t, pose);
  // camera
  const cfg = CAM3[V3.cam];
  let target, dist;
  if (V3.cam === 'genel') {
    const bs = [r, S.ref].filter(Boolean).map(boundsOf);
    const x0 = Math.min(...bs.map((q) => q.x0)), x1 = Math.max(...bs.map((q) => q.x1)), z0 = Math.min(...bs.map((q) => q.z0)), z1 = Math.max(...bs.map((q) => q.z1)), y1 = Math.max(...bs.map((q) => q.y1));
    const span = Math.max(8, y1, x1 - x0, z1 - z0);
    target = [(x0 + x1) / 2, y1 * 0.42, (z0 + z1) / 2]; dist = span * 1.25 + 10;
  } else if (V3.cam === 'yakin') { target = pose.at((D.uc0 + D.up1) / 2); dist = cfg.dist; }
  else { target = pose.center; dist = cfg.dist; }
  dist *= V3.zoom;
  const yaw = (cfg.yaw + V3.yaw) * Math.PI / 180, pit = clamp(cfg.pitch + V3.pitch, -4, 86) * Math.PI / 180;
  const cp = [target[0] - Math.sin(yaw) * Math.cos(pit) * dist, target[1] + Math.sin(pit) * dist, target[2] + Math.cos(yaw) * Math.cos(pit) * dist];
  cp[1] = Math.max(0.12, cp[1]);
  E.camera.position.set(cp[0], cp[1], cp[2]); E.camera.lookAt(target[0], target[1], target[2]);
  E.sky.position.copy(E.camera.position);
  E.sun.position.copy(E.camera.position).addScaledVector(E.sunDir, 800);
  // shadow box follows the action
  const sh = V3.cam === 'genel' ? clamp(dist * 0.5, 14, 60) : 12, L = E.light, sc = L.shadow.camera;
  if (sc.right !== sh) { sc.left = -sh; sc.right = sh; sc.top = sh; sc.bottom = -sh; sc.updateProjectionMatrix(); }
  const fx = V3.cam === 'genel' ? target[0] : pose.center[0], fz = V3.cam === 'genel' ? target[2] : pose.center[2];
  L.target.position.set(fx, 0, fz); L.position.set(fx + E.sunDir.x * 60, E.sunDir.y * 60, fz + E.sunDir.z * 60);
  const us = V3.h * E.renderer.getPixelRatio() / (2 * Math.tan(E.camera.fov * Math.PI / 360));
  E.drops.mat.uniforms.uScale.value = us; E.puffs.mat.uniforms.uScale.value = us; E.dust.mat.uniforms.uScale.value = us;
  E.renderer.render(E.scene, E.camera);
  hud3Draw(r, t, pose, target);
}
function proj3(v) {
  const E = V3.E, p = new E.T.Vector3(v[0], v[1], v[2]).project(E.camera);
  return { x: (p.x + 1) / 2 * V3.w, y: (1 - p.y) / 2 * V3.h, ok: p.z > -1 && p.z < 1 };
}
function hud3Draw(r, t, pose, target) {
  const ctx = hctx, W = V3.w, H = V3.h, f = F(), n = pose.now;
  ctx.clearRect(0, 0, W, H);
  // world labels
  const lab = (v, text, o) => { const q = proj3(v); if (q.ok && q.x > -60 && q.x < W + 60 && q.y > 64 && q.y < H - 50) pill(ctx, q.x, q.y, text, Object.assign({ align: 'center', font: `500 11px ${f.m}` }, o)); };
  [10, 20, 30, 40].forEach((d) => lab([d, 0.05, 0.9], d + ' m', { bg: 'rgba(4,12,22,0.45)', font: `500 10px ${f.m}` }));
  lab([-6.1, 2.15, 0], 'ekip · 5 m', { bg: 'rgba(4,12,22,0.45)', font: `500 10px ${f.m}` });
  if (r.left && t >= r.tApo) lab([r.xApo, r.apo + PAD_H + 0.6, r.zApo], 'tepe ' + nf(r.apo, 1) + ' m');
  if (r.left && t >= r.time) lab([r.xLand, 0.6, r.zLand], 'mesafe ' + nf(r.dist, 1) + ' m', { color: '#ffe3a0' });
  // tank labels and apparent gravity in the close view
  if (V3.cam === 'yakin') {
    const D = V3.E.dims;
    [['Çalkalanma tankı', (D.up0 + D.up1) / 2], ['Basınç haznesi', (D.uc0 + D.uc1) / 2]].forEach(([txt, u]) => {
      const q = proj3(pose.at(u)); if (!q.ok || q.y < 64 || q.y > H - 50) return;
      ctx.font = `600 11px ${f.b}`; const tw = ctx.measureText(txt).width + 14;
      const lx = clamp(q.x + 46, 6, W - 6 - tw);
      ctx.strokeStyle = 'rgba(255,255,255,0.7)'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(q.x + 4, q.y); ctx.lineTo(lx, q.y); ctx.stroke();
      ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(q.x + 4, q.y, 2, 0, Math.PI * 2); ctx.fill();
      pill(ctx, lx, q.y, txt, { font: `600 11px ${f.b}` });
    });
    const c = pose.at((D.up0 + D.up1) / 2), q0 = proj3(c), q1 = proj3([c[0] + pose.gS[0] * 0.2, c[1] + pose.gS[1] * 0.2, c[2] + pose.gS[2] * 0.2]);
    let dx = q1.x - q0.x, dy = q1.y - q0.y; const dl = Math.hypot(dx, dy) || 1; dx /= dl; dy /= dl;
    const gx = W - 38, gyy = 104, gr = 22;
    ctx.fillStyle = 'rgba(4,12,22,0.62)'; ctx.beginPath(); ctx.arc(gx, gyy, gr + 6, 0, Math.PI * 2); ctx.fill();
    ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(gx, gyy, gr, 0, Math.PI * 2); ctx.stroke();
    ctx.strokeStyle = '#ffb547'; ctx.fillStyle = '#ffb547'; ctx.lineWidth = 2.2;
    ctx.beginPath(); ctx.moveTo(gx - dx * 8, gyy - dy * 8); ctx.lineTo(gx + dx * (gr - 6), gyy + dy * (gr - 6)); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(gx + dx * gr, gyy + dy * gr); ctx.lineTo(gx + dx * (gr - 9) - dy * 5, gyy + dy * (gr - 9) + dx * 5); ctx.lineTo(gx + dx * (gr - 9) + dy * 5, gyy + dy * (gr - 9) - dx * 5); ctx.closePath(); ctx.fill();
    const gm = pose.pad || pose.landed ? 1 : n.gm;
    pill(ctx, gx, gyy + gr + 17, nf(gm, gm < 10 ? 1 : 0) + ' g', { align: 'center', font: `500 10px ${f.m}`, h: 17, color: '#ffd59a' });
    pill(ctx, gx - gr - 10, gyy, 'GÖRÜNÜR g', { align: 'right', font: `500 9px ${f.m}`, h: 17, color: '#ffd59a' });
    pill(ctx, 10, 84, 'ÇALKALANMA %' + nf(pose.sa * 100, 0), { font: `500 10.5px ${f.m}`, color: '#bfe9ff', stroke: 'rgba(76,195,242,0.45)' });
  }
  // time and phase
  const ph = phaseOf(n), pc = ph.code === 1 ? '#7fdcff' : ph.code === 2 ? '#ffc874' : '#d7e4ef';
  const tt = 'T+' + nf(Math.min(t, r.time), 2) + ' s', pl = ph.label.toLocaleUpperCase('tr-TR');
  ctx.font = `500 10.5px ${f.m}`; const pw = ctx.measureText(pl).width; ctx.font = `600 19px ${f.m}`;
  const bw = Math.max(ctx.measureText(tt).width, pw) + 22;
  rrect(ctx, 10, 10, bw, 52, 9); ctx.fillStyle = 'rgba(4,12,22,0.62)'; ctx.fill();
  ctx.fillStyle = '#fff'; ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic'; ctx.fillText(tt, 21, 34);
  ctx.font = `500 10.5px ${f.m}`; ctx.fillStyle = pc; ctx.fillText(pl, 21, 52);
  // readouts
  const still = pose.pad || pose.landed;
  const rows = [['İRTİFA', nf(still ? 0 : Math.max(0, n.y), 1) + ' m'], ['HIZ', nf(still ? 0 : n.v, 1) + ' m/s'], ['YANA', nf(n.z, 1) + ' m']];
  if (W >= 420) {
    ctx.font = `500 12px ${f.m}`;
    const vw = Math.max(...rows.map((q) => ctx.measureText(q[1]).width)) + 70;
    rrect(ctx, W - vw - 10, 10, vw, 66, 9); ctx.fillStyle = 'rgba(4,12,22,0.62)'; ctx.fill();
    rows.forEach((q, i) => {
      const Y = 28 + i * 18;
      ctx.font = `500 9.5px ${f.m}`; ctx.fillStyle = 'rgba(200,220,235,0.75)'; ctx.textAlign = 'left'; ctx.fillText(q[0], W - vw, Y);
      ctx.font = `500 12px ${f.m}`; ctx.fillStyle = '#fff'; ctx.textAlign = 'right'; ctx.fillText(q[1], W - 20, Y);
    });
  }
  // wind arrow as seen from the current camera
  const wx = r.p.wind, wz = r.p.windZ || 0, wl = Math.hypot(wx, wz);
  rrect(ctx, 10, H - 48, 176, 38, 8); ctx.fillStyle = 'rgba(4,12,22,0.55)'; ctx.fill();
  ctx.font = `500 10px ${f.m}`; ctx.fillStyle = '#e6f1fa'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
  ctx.fillText(wl < 0.05 ? 'RÜZGÂR YOK' : `RÜZGÂR ${nf(wl, 1)} m/s`, 46, H - 35);
  ctx.fillStyle = 'rgba(200,220,235,0.75)'; ctx.fillText(`ileri ${nf(wx, 1)} · yan ${nf(wz, 1)}`, 46, H - 21);
  ctx.strokeStyle = 'rgba(255,255,255,0.35)'; ctx.beginPath(); ctx.arc(28, H - 29, 12, 0, Math.PI * 2); ctx.stroke();
  if (wl > 0.05) {
    const g0 = proj3([target[0], 0, target[2]]), g1 = proj3([target[0] + wx / wl, 0, target[2] + wz / wl]);
    let dx = g1.x - g0.x, dy = g1.y - g0.y; const dl = Math.hypot(dx, dy) || 1; dx /= dl; dy /= dl;
    ctx.strokeStyle = '#ffd25e'; ctx.fillStyle = '#ffd25e'; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(28 - dx * 8, H - 29 - dy * 8); ctx.lineTo(28 + dx * 6, H - 29 + dy * 6); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(28 + dx * 11, H - 29 + dy * 11); ctx.lineTo(28 + dx * 4 - dy * 4, H - 29 + dy * 4 + dx * 4); ctx.lineTo(28 + dx * 4 + dy * 4, H - 29 + dy * 4 - dx * 4); ctx.closePath(); ctx.fill();
  }
  if (W >= 520) { ctx.font = `500 10px ${f.m}`; ctx.fillStyle = 'rgba(255,255,255,0.7)'; ctx.textAlign = 'right'; ctx.fillText('sürükle: döndür · tekerlek: yakınlaş · çift tık: sıfırla', W - 12, H - 14); }
  const vg = ctx.createRadialGradient(W / 2, H / 2, Math.min(W, H) * 0.45, W / 2, H / 2, Math.hypot(W, H) * 0.62);
  vg.addColorStop(0, 'rgba(2,8,16,0)'); vg.addColorStop(1, 'rgba(2,8,16,0.35)');
  ctx.fillStyle = vg; ctx.fillRect(0, 0, W, H);
}
function wire3() {
  const ptrs = new Map();
  let pinch = 0;
  const redraw = () => { if (!S.playing) drawFrame(); };
  hud3.addEventListener('pointerdown', (e) => { hud3.setPointerCapture(e.pointerId); ptrs.set(e.pointerId, [e.clientX, e.clientY]); hud3.classList.add('drag'); if (ptrs.size === 2) { const [a, b] = [...ptrs.values()]; pinch = Math.hypot(a[0] - b[0], a[1] - b[1]); } });
  hud3.addEventListener('pointermove', (e) => {
    if (!ptrs.has(e.pointerId)) return;
    const prev = ptrs.get(e.pointerId); ptrs.set(e.pointerId, [e.clientX, e.clientY]);
    if (ptrs.size === 1) { V3.yaw -= (e.clientX - prev[0]) * 0.35; V3.pitch = clamp(V3.pitch + (e.clientY - prev[1]) * 0.25, -60, 80); }
    else if (ptrs.size === 2) { const [a, b] = [...ptrs.values()]; const d = Math.hypot(a[0] - b[0], a[1] - b[1]); if (pinch > 0) zoom3(pinch / d); pinch = d; }
    redraw();
  });
  const up = (e) => { ptrs.delete(e.pointerId); if (!ptrs.size) hud3.classList.remove('drag'); pinch = 0; };
  hud3.addEventListener('pointerup', up); hud3.addEventListener('pointercancel', up);
  hud3.addEventListener('wheel', (e) => { e.preventDefault(); zoom3(Math.exp(e.deltaY * 0.0012)); redraw(); }, { passive: false });
  hud3.addEventListener('dblclick', () => { V3.yaw = 0; V3.pitch = 0; V3.zoom = 1; redraw(); });
}
function zoom3(f) { const c = CAM3[V3.cam]; V3.zoom = clamp(V3.zoom * f, c.zmin, c.zmax); }
function syncCam3() { $$('[data-cam3]').forEach((b) => b.setAttribute('aria-pressed', b.dataset.cam3 === V3.cam ? 'true' : 'false')); }
function setDim(d) {
  S.dim = d === '3' ? '3' : '2'; store.set('ikarus-dim', S.dim);
  $('#dim-' + S.dim).checked = true;
  $('#feed2').hidden = S.dim === '3'; $('#feed3').hidden = S.dim !== '3';
  $('#dim-hint').textContent = S.dim === '3'
    ? 'Sahneyi sürükleyerek çevirin, tekerlek veya iki parmakla yakınlaşın. "Yakın" kamera tanklardaki suyu gösterir.'
    : 'Yandan iki kamera: fırlatma alanı ve roketin içi (kesit).';
  if (S.dim === '3') ensure3().then((ok) => { if (ok && S.view === 'sim' && S.tab === 'ucus') { V3.w = 0; size3(); drawFrame(); } });
  else if (S.view === 'sim' && S.tab === 'ucus') { sizeCanvases(); drawFrame(); }
}

// ------------------------------------------------------------------ captions, telemetry
const PH_NAME = ['Fırlatma tüpü', 'Su evresi', 'Hava evresi', 'Süzülme'];
const noseDeg = (th) => { let d = th * 180 / Math.PI; while (d > 180) d -= 360; while (d <= -180) d += 360; return d; };
function headingTxt(bx, bz) {
  const a = Math.atan2(bz, bx) * 180 / Math.PI;
  return Math.abs(a) < 15 ? 'ileri' : Math.abs(a) > 165 ? 'geri' : `${a > 0 ? 'sağa' : 'sola'} ${nf(Math.abs(a), 0)}°`;
}
function phaseOf(now) {
  const r = S.cur;
  if (!r.left) return { label: 'Rampadan ayrılamadı', code: 9 };
  if (S.t <= 0) return { label: 'Rampada, basınçlandırıldı', code: 9 };
  if (S.t >= r.time - 1e-6) return { label: r.truncated ? 'Hesap sınırı' : 'Yere indi', code: 9 };
  if (now.chute) return { label: 'Paraşütle iniş', code: 3 };
  if (S.t < (r.tAirEnd || 0)) return { label: PH_NAME[now.ph] + (S.playing && $('#slowmo').checked ? ' · ağır çekim' : ''), code: now.ph };
  if (r.G.wing && now.wd > 0.5 && S.t > r.tApo) return { label: 'Süzülme', code: 3 };
  return { label: S.t < r.tApo ? 'İtkisiz tırmanış' : 'Düşüş', code: 3 };
}
function waterSide(gs) {
  const c = Math.cos(gs);
  return c > 0.35 ? 'burun tarafına' : c < -0.35 ? 'kıç (nozul) tarafına' : 'yan duvara';
}
function captionFor(r, t, pose) {
  const p = r.p, n = pose.now;
  if (!r.left) return '<b>Rampadan ayrılamadı.</b> İtki, roketi 1 m\'lik rampa boyunca yeterince hızlandıramadan sönüyor. Basıncı veya nozul çapını artırın ya da kütleyi azaltın. Hesap 3 s\'de durduruldu.';
  if (t <= 0) return `<b>Rampada, ${nf(p.P, 1)} bar.</b> İtki suyu nozul tarafında. Rampa ${nf(p.angle, 1)}° eğik; çalkalanma tankındaki suyun yüzeyi yere paralel.`;
  if (pose.landed && r.truncated) return `<b>Hesap ${nf(r.time, 0)} s'de durduruldu.</b> Roket hâlâ süzülüyordu; süre ve mesafe bu ana kadarki değerler. Puan ${nf(r.score, 2)} (en az).`;
  if (pose.landed) return `<b>Yere indi:</b> ${nf(r.dist, 1)} m uzakta, ${nf(r.time, 2)} s sonra. Tanktaki su artık yalnız yerçekimini hissediyor. Puan ${nf(r.score, 2)}.`;
  if (n.chute) return `<b>Paraşütle iniş.</b> Yavaş düşüş süreyi uzatıyor. Hava direnci ağırlığı dengeliyor; su tankın ${waterSide(pose.gs)} çöküyor.`;
  if (t < (r.tAirEnd || 0)) {
    if (n.ph === 0) return '<b>Fırlatma tüpü.</b> Su henüz çıkmıyor; basınçlı hava tüpün içinde roketi itiyor.';
    if (n.ph === 1) return `<b>Su evresi · ${nf(n.T, 0)} N.</b> Su nozuldan ~${nf(r.G.cd * Math.sqrt(200 * Math.max(0, n.P)), 0)} m/s ile fışkırıyor. Çalkalanan su ${nf(n.gm, 0)} g ile kıç tarafına bastırılıyor.`;
    return '<b>Hava evresi.</b> Kalan basınçlı hava, genleşip soğuyarak sis bulutu hâlinde çıkıyor; itki milisaniyeler içinde sönüyor.';
  }
  const pct = nf(pose.sa * 100, 0);
  if (r.G.wing && n.wd > 0.5 && t > r.tApo && !n.chute) {
    const vh = Math.hypot(n.vx, n.vz), sink = -n.vy;
    return `<b>Süzülme.</b> Kanatlar ${nf(Math.max(0, n.Lw), 0)} N kaldırma üretiyor; roket ${nf(vh, 1)} m/s ileri giderken ${sink >= 0 ? `${nf(sink, 1)} m/s alçalıyor` : `${nf(-sink, 1)} m/s yükseliyor`}${sink > 0.3 ? ` (süzülme oranı ${nf(vh / sink, 1)})` : ''}. Su tankın ${waterSide(pose.gs)} yaslanıyor.`;
  }
  if (t < r.tApo) return `<b>İtkisiz tırmanış.</b> İtki bitti; hava direnci roketi yavaşlatıyor ve su ${waterSide(pose.gs)} savruluyor (çalkalanma %${pct}). ${p.baffle > 0 ? 'Levhalar salınımı sönümlüyor.' : 'Bölme levhası olmadığı için salınım uzun sürüyor.'}`;
  return `<b>Düşüş.</b> Roket tepe noktasını geçti. Su ${waterSide(pose.gs)} yaslanıyor (çalkalanma %${pct}); havada geçen her saniye puana ekleniyor.`;
}
let lastCap = '';
function drawFrame() {
  const r = S.cur; if (!r) return;
  const t = S.t, pose = poseAt(r, t);
  if (S.dim === '3') draw3D(r, t);
  else { drawScene(r, t, pose); drawCloseup(r, t, pose); }
  renderTele(pose.now);
  const cap = captionFor(r, t, pose);
  if (cap !== lastCap) { $('#cap').innerHTML = cap; $('#cap3').innerHTML = cap; lastCap = cap; }
  $('#pb-time').textContent = `${nf(Math.min(t, r.time), 2)} / ${nf(r.time, 2)} s`;
  markEvents();
}

function renderTele(now) {
  const r = S.cur;
  const ph = phaseOf(now);
  const el = $('#t-ph'); el.textContent = ph.label; el.dataset.ph = String(ph.code);
  const onPad = S.t <= 0;
  const items = [
    ['Zaman', nf(Math.min(S.t, r.time), 3) + ' s'],
    ['İrtifa', nf(Math.max(0, now.y), 1) + ' m'],
    ['Konum (ileri / yan)', nf(now.x, 1) + ' / ' + nf(now.z, 1) + ' m'],
    ['Hız', nf(onPad || S.t >= r.time ? 0 : now.v, 1) + ' m/s'],
    ['Eksenel ivme', nf(onPad ? 0 : now.ax, 1) + ' g'],
    ['Hazne basıncı', nf(onPad ? r.p.P : now.P, 2) + ' bar'],
    ['Kalan itki suyu', nf(onPad ? r.G.mProp * 1000 : now.Vw, 0) + ' ml'],
    ['İtki', nf(onPad ? 0 : now.T, 0) + ' N'],
    ['Burun açısı (yataydan)', nf(Math.asin(clamp(now.by, -1, 1)) * 180 / Math.PI, 1) + '°'],
    ['Burun yönü', Math.abs(now.by) > 0.995 ? 'dikey' : headingTxt(now.bx, now.bz)],
    ['Hücum açısı', nf(onPad || now.chute ? 0 : Math.abs(now.al), 1) + '°'],
    ['Kütle', nf(onPad ? r.G.TOW + r.maAir0 : now.m, 2) + ' kg'],
    ['Çalkalanma', '%' + nf(S.t >= r.time ? 0 : now.sa * 100, 0)],
  ];
  if (r.G.wing) items.push(['Kanat kaldırması', nf(onPad ? 0 : Math.max(0, now.Lw), 1) + ' N'], ['Alçalma hızı', nf(onPad || S.t >= r.time ? 0 : Math.max(0, -now.vy), 1) + ' m/s']);
  $('#tele').innerHTML = items.map(([k, v]) => `<div><span class="k">${k}</span><span class="v">${v}</span></div>`).join('');
}

function renderEvents() {
  $('#events').innerHTML = S.cur.ev.map((e, i) => `<li data-i="${i}"><button type="button" data-t="${e.t}"><span class="et">${nf(e.t, e.t < 1 ? 3 : 2)} s</span><span>${e.txt}</span></button></li>`).join('');
  markEvents();
}
function markEvents() {
  const lis = $$('#events li');
  let lastPast = -1;
  lis.forEach((li, i) => { const e = S.cur.ev[i]; const past = !!e && e.t <= S.t + 1e-6; li.classList.toggle('past', past); if (past) lastPast = i; });
  lis.forEach((li, i) => li.classList.toggle('now', i === lastPast));
}

// ------------------------------------------------------------------ score tab
function renderScore() {
  const r = S.cur, p = r.p, G = r.G, ref = S.ref;
  $('#score').textContent = nf(r.score, 2);
  const d = $('#delta');
  deltaTo(d, r.score, ref.score);
  $('#calc').innerHTML = r.left
    ? `<span>( √(${nf(r.dist, 1)}² + ${nf(r.apo, 1)}²) + ${nf(r.time, 2)} ) × ${nf(G.mPay, 2)} / ${nf(G.TOW, 2)}</span><span>= ( ${nf(r.root, 2)} + ${nf(r.time, 2)} ) × ${nf(r.ratio, 3)} = ${nf(r.score, 2)}</span>`
    : '<span>Roket rampadan ayrılamadı: itki, ağırlığı yenmeye yetmiyor. Basıncı veya nozulu artırın ya da kütleyi azaltın.</span>';
  const parts = [
    ['Mesafe-irtifa terimi', r.root, r.root * r.ratio, css('--accent')],
    ['Süre terimi', r.time, r.time * r.ratio, css('--water')],
  ];
  const mx = Math.max(1, r.score);
  $('#parts').innerHTML = !r.left ? '<span class="small">Roket rampadan ayrılmadığı için puan katkısı yok.</span>' : `<span class="small">Puana katkı (çarpanla birlikte)</span>` + parts.map(([k, raw, sc, col]) => `<div class="part"><span>${k}</span><span class="part-bar"><i style="width:${Math.min(100, sc / mx * 100).toFixed(1)}%;background:${col}"></i></span><span class="pv">${nf(sc, 2)}</span></div>`).join('')
    + `<div class="part"><span>Çarpan (yük/ağırlık)</span><span class="part-bar"><i style="width:${(r.ratio * 100).toFixed(1)}%;background:var(--c-prop)"></i></span><span class="pv">${nf(r.ratio, 3)}</span></div>`;
  const m = [
    ['Tepe irtifa', nf(r.apo, 1) + ' m', nf(ref.apo, 1)],
    ['Yatay mesafe', nf(r.dist, 1) + ' m', nf(ref.dist, 1)],
    ['Uçuş süresi', r.left ? nf(r.time, 2) + ' s' : '–', ref.left ? nf(ref.time, 2) : '–'],
    ['Azami hız', nf(r.maxV, 1) + ' m/s', nf(ref.maxV, 1)],
    ['Azami ivme', nf(r.maxA / C.g, 0) + ' g', nf(ref.maxA / C.g, 0)],
    ['Toplam itki', nf(r.impulse, 1) + ' N·s', nf(ref.impulse, 1)],
    ['İtki süresi', r.left ? nf((r.tAirEnd || 0) * 1000, 0) + ' ms' : '–', ref.left ? nf((ref.tAirEnd || 0) * 1000, 0) : '–'],
    ['Rampa çıkış hızı', r.vRail ? nf(r.vRail, 1) + ' m/s' : '–', ref.vRail ? nf(ref.vRail, 1) : '–'],
    ['Yük / kalkış ağırlığı', nf(r.ratio, 3), nf(ref.ratio, 3)],
  ];
  $('#metrics').innerHTML = m.map(([k, v, rv]) => `<div class="metric"><div class="k">${k}</div><div class="v">${v}</div><div class="r">karş.: ${rv}</div></div>`).join('');
  const rules = checkRules(p, r);
  $('#rules').innerHTML = rules.map((x) => `<li><span class="chip ${x.ok ? 'ok' : 'no'}">${x.ok ? '✓ uygun' : '✕ ihlal'}</span><span>${x.label} <span class="ref">${x.ref}</span></span><span class="val">${x.val}</span></li>`).join('');
  const bad = rules.filter((x) => !x.ok).length;
  const v = $('#verdict');
  v.className = 'verdict ' + (bad || !r.left ? 'no' : 'ok');
  v.textContent = bad ? `${bad} kural ihlali var. Bu tasarım yarışmada kabul edilmez veya cezalandırılır; puan yalnızca karşılaştırma içindir.` : r.left ? 'Tasarım, modelin kontrol edebildiği bütün kurallara uyuyor.' : 'Roket uçmuyor.';
}

function tweaks(p) {
  const Vsu = p.Vc * p.wf;
  const L = [
    { lbl: 'Boş yapıyı 100 g hafiflet', ok: p.dry - 0.1 >= 0.15 - 1e-9, f: (q) => { q.dry = +(q.dry - 0.1).toFixed(3); } },
    { lbl: 'Basıncı 0,5 bar artır', ok: p.P + 0.5 <= 12, f: (q) => { q.P += 0.5; } },
    { lbl: 'Basıncı 0,5 bar azalt', ok: p.P - 0.5 >= 1, f: (q) => { q.P -= 0.5; } },
    { lbl: 'İtki suyunu 0,1 L artır', ok: Vsu + 0.1 <= p.Vc * 0.95, f: (q) => { q.wf = (q.Vc * q.wf + 0.1) / q.Vc; } },
    { lbl: 'İtki suyunu 0,1 L azalt', ok: Vsu - 0.1 >= 0.05, f: (q) => { q.wf = (q.Vc * q.wf - 0.1) / q.Vc; } },
    { lbl: 'Çalkalanan suyu 0,25 L artır', ok: p.payload + 0.25 <= 4, f: (q) => { q.payload += 0.25; } },
    { lbl: 'Hazneyi 0,5 L büyüt (su aynı)', ok: p.Vc + 0.5 <= 8, f: (q) => { const v = q.Vc * q.wf; q.Vc += 0.5; q.wf = v / q.Vc; } },
    { lbl: 'Hazneyi 0,5 L küçült (su aynı)', ok: p.Vc - 0.5 >= Vsu / 0.95 && p.Vc - 0.5 >= 0.5, f: (q) => { const v = q.Vc * q.wf; q.Vc -= 0.5; q.wf = v / q.Vc; } },
    { lbl: 'Nozulu 2 mm genişlet', ok: p.nozzle + 2 <= 30, f: (q) => { q.nozzle += 2; } },
    { lbl: 'Nozulu 2 mm daralt', ok: p.nozzle - 2 >= 4, f: (q) => { q.nozzle -= 2; } },
    { lbl: p.fins === 'kucuk' ? 'Kanatçıkları orta yap' : 'Kanatçıkları küçült', ok: true, f: (q) => { q.fins = q.fins === 'kucuk' ? 'orta' : (q.fins === 'buyuk' ? 'orta' : 'kucuk'); } },
    { lbl: p.tube ? 'Fırlatma tüpünü kaldır' : 'Fırlatma tüpü ekle', ok: true, f: (q) => { q.tube = !q.tube; } },
    { lbl: p.chute ? 'Paraşütü kaldır' : 'Paraşüt ekle', ok: true, f: (q) => { q.chute = !q.chute; } },
    p.wing === 'acilir' || p.wing === 'sabit'
      ? { lbl: 'Kanadı kaldır', ok: true, f: (q) => { q.wing = 'yok'; } }
      : { lbl: 'Açılır kanat ekle (1,6 m × 20 cm, 6°)', ok: true, f: (q) => { q.wing = 'acilir'; q.wspan = 1.6; q.wchord = 0.2; q.winc = 6; } },
  ];
  return L.filter((x) => x.ok);
}
// slow work (many simulations) is done in slices between frames; a newer request cancels an older one
let jobId = 0;
function sliced(items, each, done) {
  const id = ++jobId; let i = 0;
  const step = () => {
    if (id !== jobId) return;
    const t0 = performance.now();
    do { each(items[i], i + 1); i++; } while (i < items.length && performance.now() - t0 < 30);
    if (i < items.length) setTimeout(step, 0); else done();
  };
  if (items.length) setTimeout(step, 0); else setTimeout(done, 0);
}
function renderTips() {
  if (S.tab !== 'puan') return;
  const base = S.cur.p, b0 = S.cur.score, list = tweaks(base), rows = [];
  const ul = $('#tips');
  ul.innerHTML = `<li class="small">Öneriler hesaplanıyor…</li>`;
  sliced(list, (tw) => {
    const q = Object.assign({}, base); tw.f(q);
    const r = simulate(q, false);
    rows.push({ lbl: tw.lbl, d: r.score - b0, ok: allOk(q, r), fly: r.left, q });
  }, () => {
    if (S.tab !== 'puan') return;
    rows.sort((a, b) => (b.ok - a.ok) || (b.d - a.d));
    ul.innerHTML = rows.map((x, i) => `<li><span>${x.lbl}${x.ok ? '' : x.fly ? ' <span class="chip no">kural dışı</span>' : ' <span class="chip mid">uçmuyor</span>'}</span><span class="dv ${Math.abs(x.d) < 0.005 ? '' : x.d > 0 ? 'up' : 'down'}">${x.d >= 0 ? '+' : ''}${nf(x.d, 2)}</span><button type="button" data-tip="${i}">Uygula</button></li>`).join('');
    ul.querySelectorAll('button[data-tip]').forEach((btn) => btn.addEventListener('click', () => {
      writeControls(toCtl(rows[+btn.dataset.tip].q)); recompute();
    }));
  });
}

// ------------------------------------------------------------------ calculation report
function renderReport() {
  if (S.tab !== 'rapor' || !S.cur) return;
  const r = S.cur, p = r.p, G = r.G, g = C.g;
  const tB = r.tAirEnd || 0, tW = r.tWaterEnd || tB;
  const kp = [
    ['Uçuş puanı', nf(r.score, 2), '', allOk(p, r) ? 'kurallara uygun' : r.left ? 'kural ihlali var' : 'rampadan ayrılamadı', true],
    ['Tepe noktası', nf(r.apo, 1), 'm', `${nf(r.tApo, 2)} s'de`],
    ['Havada kalma süresi', r.left ? nf(r.time, 2) : '–', r.left ? 's' : '', !r.left ? 'roket rampada kaldı' : r.truncated ? 'hesap sınırında kesildi' : r.tChute ? 'paraşüt dahil' : r.tWing ? 'süzülme dahil' : ''],
    ['Yatay mesafe', nf(r.dist, 1), 'm', Math.abs(r.zLand) >= 0.3 ? `yana ${nf(Math.abs(r.zLand), 1)} m` : 'fırlatma noktasından'],
    ['Azami itki', nf(r.peakT, 0), 'N', `ortalama ${nf(tB > 0 ? r.impulse / tB : 0, 0)} N`],
    ['Toplam itki (impuls)', nf(r.impulse, 1), 'N·s', `özgül itki I_sp ${nf(G.mProp > 0 ? r.impulse / (G.mProp * g) : 0, 1)} s`],
    ['İtki süresi', r.left ? nf(tB * 1000, 0) : '–', r.left ? 'ms' : '', !r.left ? 'roket rampada kaldı' : r.waterLeft > 0 ? `basınç ${nf(tW * 1000, 0)} ms'de bitti · ${nf(r.waterLeft * 1e6, 0)} ml su kaldı` : `su ${nf(tW * 1000, 0)} ms · hava ${nf((tB - tW) * 1000, 0)} ms`],
    ['Azami ivme (G yükü)', nf(r.maxA / g, 1), 'g', `itki boyunca ortalama ${nf(tB > 0 ? r.vBurn / tB / g : 0, 1)} g`],
    ['Azami hız', nf(r.maxV, 1), 'm/s', `itki bitiminde ${nf(r.vBurn, 1)} m/s`],
    ['Rampa çıkış hızı', r.vRail ? nf(r.vRail, 1) : '–', 'm/s', `itki/ağırlık ${nf(G.tw, 1)} : 1`],
    ['Azami yavaşlama', nf(r.maxDecel / g, 1), 'g', 'itki bittikten sonra, sürtünmeyle'],
    ['Yere çarpma hızı', nf(r.vLand, 1), 'm/s', r.tChute ? 'paraşütle' : r.tWing ? 'süzülmeden sonra' : 'serbest düşüş'],
    ['Kalkış kütlesi', nf(G.TOW, 2), 'kg', `boş ${nf(G.mDry, 2)} kg`],
    ['Yük / kalkış ağırlığı', nf(r.ratio, 3), '', 'puan çarpanı'],
  ];
  if (G.wing) kp.push(['Azami kanat kaldırması', nf(r.maxLift, 0), 'N', r.tWing ? `kanatlar ${nf(r.tWing, 2)} s'de açıldı` : 'sabit kanat']);
  kp.push(['Ortalama çalkalanma', '%' + nf(r.slosh * 100, 0), '', `sönüm ζ ${nf(r.zeta, 2)}`], ['Kararlılık payı', nf(G.sm0, 1), 'kalibre', `su bitince ${nf(G.sm1, 1)}` + (G.wing === 'acilir' ? ` · kanat açık ${nf(G.smW, 1)}` : '')]);
  $('#rep-kpi').innerHTML = kp.map(([k, v, u, d, hi]) => `<div class="kpi${hi ? ' hi' : ''}"><span class="k">${k}</span><span class="v">${v}${u ? `<small>${u}</small>` : ''}</span><span class="d">${d}</span></div>`).join('');
  // formulas with the numbers at t = 0
  // with a launch tube the water phase starts after the tube push, at a lower pressure
  const tubeOn = !!p.tube, t0l = tubeOn ? 'Su evresi başında' : 't = 0';
  const Pabs = tubeOn ? r.Pstart : p.P * 1e5 + C.Patm, dP = Pabs - C.Patm;
  const vs = G.cd * Math.sqrt(2 * dP / C.rhoW), md = C.rhoW * G.An * vs, F0 = md * vs;
  const m0 = G.TOW + r.maAir0, a0 = (F0 - m0 * g * Math.sin(p.angle * Math.PI / 180)) / m0;
  const vmx = r.maxV, Fd = 0.5 * C.rhoA * G.CD0 * G.Ab * vmx * vmx, Va0 = (G.Vc - G.Vw0 + (tubeOn ? G.Atube * 0.30 : 0)) * 1000;
  const eq = [
    ['İç basınç (havanın adyabatik genleşmesi)', 'P_anlık = P₀ · (V_hava0 / V_anlık_hava)^1,4', tubeOn
      ? `Önce fırlatma tüpü evresi: F = ΔP · A_tüp = ${nf(p.P, 2)} bar · ${nf(G.Atube * 1e6, 0)} mm² = ${nf(G.T0, 0)} N. Roket tüpten çıkınca su evresi P = ${nf(Pabs / 1e5, 3)} bar (mutlak), V_hava = ${nf(Va0, 2)} L ile başlar.`
      : `t = 0: P₀ = ${nf(Pabs / 1e5, 3)} bar (mutlak), V_hava0 = ${nf(Va0, 2)} L`],
    ['Basınç farkı', 'ΔP = P_anlık − P_atmosfer', `= ${nf(Pabs / 1e5, 3)} − 1,013 = ${nf(dP / 1e5, 2)} bar = ${nf(dP, 0)} Pa`],
    ['Suyun çıkış hızı', 'v_su = C_d · √(2 · ΔP / ρ_su)', `= ${nf(G.cd, 2)} · √(2 · ${nf(dP, 0)} / 1000) = ${nf(vs, 1)} m/s`],
    ['Su kaybı (debi)', 'ṁ = ρ_su · A_n · v_su', `= 1000 · ${nf(G.An, 6)} · ${nf(vs, 1)} = ${nf(md, 2)} kg/s`],
    ['Üretilen itki', 'F_itki = ṁ · v_su', `= ${nf(md, 2)} · ${nf(vs, 1)} = ${nf(F0, 0)} N`],
    ['Kalan su ve yeni hava hacmi', 'V_su ← V_su − (ṁ / ρ_su) · Δt ;  V_hava = V_tank − V_su', `Δt = 0,5 ms'lik ilk adımda ${nf(md / C.rhoW * 0.0005 * 1e6, 2)} ml su çıkar.`],
    ['Su bitince', 'V_su ≤ 0  ⇒  su itkisi = 0', !r.left ? 'Roket rampadan ayrılamadığı için hesap 3 s\'de durduruldu.' : r.waterLeft > 0 ? `Su bitmeden basınç atmosfere indi (${nf(tW * 1000, 0)} ms): haznede ${nf(r.waterLeft * 1e6, 0)} ml su kaldı, hava evresi olmadı. Daha az su ya da daha yüksek basınç kullanın.` : `Su ${nf(tW * 1000, 0)} ms'de bitti. Ardından hazne ${nf(r.P1g, 1)} bar'dan boşalırken hava ${nf((tB - tW) * 1000, 0)} ms daha itti.`],
    ['Anlık toplam kütle', 'm_toplam = m_boş + m_yük + V_su · ρ_su + m_hava', `t = 0: ${nf(G.mDry, 3)} + ${nf(G.mPay, 3)} + ${nf(G.mProp, 3)} + ${nf(r.maAir0, 3)} = ${nf(m0, 3)} kg. Sıkıştırılmış havanın kütlesi m_hava = P₀ · V_hava0 / (R · T) uçuşta hesaba katılır; jürinin tarttığı kalkış ağırlığına (${nf(G.TOW, 3)} kg) dahil değildir.`],
    ['Hava sürtünmesi', 'F_sürtünme = ½ · ρ_hava · C_D · A_b · v²', `Azami hızda (${nf(vmx, 1)} m/s): ½ · 1,2 · ${nf(G.CD0, 2)} · ${nf(G.Ab, 5)} · ${nf(vmx, 1)}² = ${nf(Fd, 2)} N`],
    ['Net kuvvet ve ivme (roket ekseni boyunca)', 'F_net = F_itki − m_toplam · g · sinθ − F_sürtünme ;  a = F_net / m_toplam', `${t0l}: (${nf(F0, 0)} − ${nf(m0, 2)} · 9,81 · sin ${nf(p.angle, 1)}°) / ${nf(m0, 2)} = ${nf(a0, 0)} m/s² = ${nf(a0 / g, 1)} g. Roketin hissettiği G yükü yerçekimini içermez: F_itki / m = ${nf(F0 / m0 / g, 1)} g.`],
    ['Yeni hız ve yükseklik', 'v ← v + a · Δt ;  h ← h + v · Δt', 'Aynı denklemler, Euler yerine daha hassas olan RK4 yöntemiyle çözülür.'],
  ];
  if (G.wing) eq.push(['Kanat kaldırması ve sürtünmesi', 'L = ½ · ρ · v² · S · C_L ;  C_L = C_Lα · (α + i) ;  D = ½ · ρ · v² · S · (0,02 + C_L² / (π · 0,8 · AR))',
    `S = ${nf(G.Sw, 3)} m², AR = ${nf(G.AR, 1)}, C_Lα = ${nf(G.CLa, 2)} /rad, i = ${nf(p.winc, 1)}°. Uçuştaki en büyük kaldırma ${nf(r.maxLift, 0)} N.`]);
  $('#rep-eq').innerHTML = eq.map(([n, f, x]) => `<li><div><div class="n">${n}</div><div class="f">${f}</div><div class="x">${x}</div></div></li>`).join('');
  $('#rep-diff').innerHTML = '<b>Sizin formüllerinizle karşılaştırma</b><p>Su evresinin denklemleri (1–7) birebir aynı; 8. denkleme sıkıştırılmış havanın kütlesi eklendi. Simülasyon bunlara şunları ekler: ' + (S.cur.p.tube ? 'fırlatma tüpü evresi, ' : '') + 'su bittikten sonra kalan basınçlı havanın itkisi (hava evresi), 1 m\'lik rampa ve fırlatma açısı, üç boyutlu hareket ve rüzgâr, hücum açısıyla artan sürtünme, kanatçıkların ve varsa kanadın kuvvetleri, çalkalanan suyun kütlesi ve hareketi. Formül 10 rampadaki an için eksen boyunca yazıldı; uçuşta kuvvetler vektör olarak toplanır.</p>';
  // inputs
  const rows = [
    ['V_tank', 'Basınç haznesi (toplam) hacmi', nf(p.Vc, 2) + ' L'],
    ['V_su', 'Başlangıçtaki itki suyu', nf(p.Vc * p.wf, 2) + ' L'],
    ['V_hava0', 'Başlangıçtaki hava hacmi (V_tank − V_su)', nf(p.Vc * (1 - p.wf), 2) + ' L'],
    ['V_yük', 'Çalkalanma tankındaki su (tank %50 dolu)', nf(p.payload, 2) + ' L'],
    ['P_basılan', 'Pompayla basılan (gösterge) basınç', nf(p.P, 2) + ' bar'],
    ['P₀', 'Mutlak başlangıç basıncı (P_basılan + 1 atm)', nf(p.P + 1.01325, 3) + ' bar'],
    ['P_atm', 'Atmosfer basıncı', '1,01325 bar'],
    ['C_d', 'Nozul kayıp faktörü', nf(G.cd, 2)],
    ['D_nozul → A_n', 'Nozul çapı ve alanı, π · (D/2)²', `${nf(p.nozzle, 1)} mm → ${nf(G.An * 1e6, 0)} mm²`],
    ['D_şişe → A_b', 'Şişe çapı ve kesit alanı', `${nf(p.Dmm, 0)} mm → ${nf(G.Ab * 1e4, 1)} cm²`],
    ['C_D', 'Hava sürtünme katsayısı (gövde + kanatçık)', `${nf(p.CD, 2)} + ${nf(G.fin.cd, 2)} = ${nf(G.CD0, 2)}`],
    ['m_boş', 'Susuz kütle (yapı + kanatçık' + (G.wing ? ' + kanat ve kuyruk' : '') + ')', nf(G.mDry * 1000, 0) + ' g'],
    ['m_yük', 'Çalkalanan su kütlesi', nf(G.mPay * 1000, 0) + ' g'],
    ['m_su', 'İtki suyu kütlesi', nf(G.mProp * 1000, 0) + ' g'],
    ['m_toplam', 'Kalkış kütlesi', nf(G.TOW * 1000, 0) + ' g'],
    ['ρ_su · ρ_hava', 'Su ve hava yoğunluğu', '1000 · 1,2 kg/m³'],
    ['g · γ', 'Yerçekimi ivmesi · havanın özgül ısı oranı', '9,81 m/s² · 1,4'],
    ['θ', 'Fırlatma açısı (yataydan); rampa 1 m', nf(p.angle, 1) + '°'],
    ['Δt', 'Zaman adımı (itki · sonrası)', '0,5 ms · 10 ms'],
    ['b_levha · ζ', 'Bölme levhası (tank kesitine oranı) · sönüm oranı', `%${nf(p.baffle * 100, 1)} · ${nf(r.zeta, 2)}`],
    ['Kanatçık', 'Tip: açıklık × kök genişliği · kütle', `${G.fin.name}: ${nf(G.fin.span * 100, 1)} × ${nf(G.fin.chord * 100, 0)} cm · ${nf(G.fin.mass * 1000, 0)} g`],
    ['Tüp · paraşüt', 'Fırlatma tüpü (0,30 m) · tepeden 0,4 s sonra açılan paraşüt', `${p.tube ? 'var' : 'yok'} · ${p.chute ? 'var' : 'yok'}`],
    ['v_rüzgâr', 'Rüzgâr: ileri · yan (±%25 dalgalı)', `${nf(p.wind, 1)} · ${nf(p.windZ || 0, 1)} m/s`],
    ['L', 'Roket boyu', nf(G.L, 2) + ' m'],
  ];
  if (G.wing) rows.push(['b × c', 'Kanat açıklığı × genişliği', `${nf(p.wspan, 2)} m × ${nf(p.wchord * 100, 0)} cm`], ['S · AR', 'Kanat alanı · en-boy oranı', `${nf(G.Sw, 3)} m² · ${nf(G.AR, 1)}`], ['i', 'Kanat montaj açısı', nf(p.winc, 1) + '°']);
  $('#rep-in').innerHTML = rows.map(([a, b, c]) => `<tr><td class="sy">${a}</td><td>${b}</td><td class="num-c">${c}</td></tr>`).join('');
  // time table
  const step = $('#rep-step').value, times = [];
  if (step[0] === 'p') { let d = step === 'p2' ? 0.002 : 0.005; if (tB / d > 400) d = tB / 400; for (let i = 0; i * d < tB - 1e-6; i++) times.push(i * d); times.push(tB); }
  else { const d = step === 'f10' ? 0.1 : 0.25; for (let i = 0; i * d < r.time - 1e-6; i++) times.push(i * d); times.push(r.time); }
  const phName = (e, t) => t >= r.time - 1e-6 && r.left ? (r.truncated ? 'Kesildi' : 'İniş') : e.ph === 0 ? 'Tüp' : e.ph === 1 ? 'Su' : e.ph === 2 ? 'Hava' : e.chute ? 'Paraşüt' : (G.wing && e.wd > 0.5 && t > r.tApo) ? 'Süzülme' : t < r.tApo ? 'Tırmanış' : 'Düşüş';
  $('#rep-tt').innerHTML = times.map((t) => {
    const e = recAt(r, t), vsu = e.ph === 1 ? G.cd * Math.sqrt(2 * Math.max(0, e.P) * 1e5 / C.rhoW) : null;
    return `<tr><td>${nf(t, 3)}</td><td>${phName(e, t)}</td><td class="num-c">${nf(Math.max(0, e.y), 2)}</td><td class="num-c">${nf(Math.hypot(e.x, e.z), 2)}</td><td class="num-c">${nf(e.v, 2)}</td><td class="num-c">${nf(e.ax, 1)}</td><td class="num-c">${nf(e.T, 1)}</td><td class="num-c">${nf(e.P, 2)}</td><td class="num-c">${nf(e.Vw, 0)}</td><td class="num-c">${vsu == null ? '–' : nf(vsu, 1)}</td><td class="num-c">${nf(e.m, 3)}</td></tr>`;
  }).join('');
}

// ------------------------------------------------------------------ charts
const charts = {};
const cursorPlugin = {
  id: 'cursor',
  afterDatasetsDraw(chart) {
    const o = chart.options.plugins.cursor; if (!o || !o.on || !S.cur) return;
    const xv = S.t * (o.unit || 1);
    const xs = chart.scales.x; if (xv < xs.min || xv > xs.max) return;
    const x = xs.getPixelForValue(xv), a = chart.chartArea, c = chart.ctx;
    c.save(); c.strokeStyle = css('--ink'); c.globalAlpha = 0.55; c.lineWidth = 1.2; c.setLineDash([3, 3]);
    c.beginPath(); c.moveTo(x, a.top); c.lineTo(x, a.bottom); c.stroke(); c.restore();
  },
};
function chartBase(unit) {
  const ink = css('--ink-2'), grid = css('--line');
  return {
    responsive: true, maintainAspectRatio: false, animation: false, parsing: false, normalized: true,
    interaction: { mode: 'nearest', axis: 'x', intersect: false },
    plugins: {
      cursor: { on: true, unit: unit || 1 },
      legend: { labels: { color: ink, boxWidth: 14, boxHeight: 2, font: { family: css('--f-body'), size: 12 } } },
      tooltip: { backgroundColor: css('--panel'), titleColor: css('--ink'), bodyColor: css('--ink'), borderColor: grid, borderWidth: 1, titleFont: { family: css('--f-mono') }, bodyFont: { family: css('--f-mono') },
        callbacks: { title: (it) => it.length ? nf(it[0].parsed.x, unit === 1000 ? 0 : 2) + (unit === 1000 ? ' ms' : ' s') : '' } },
    },
    scales: {
      x: { type: 'linear', grid: { color: grid }, border: { color: grid }, ticks: { color: ink, font: { family: css('--f-mono'), size: 11 }, maxTicksLimit: 7 }, title: { display: true, color: ink, font: { family: css('--f-body'), size: 12 } } },
      y: { grid: { color: grid }, border: { color: grid }, ticks: { color: ink, font: { family: css('--f-mono'), size: 11 }, maxTicksLimit: 6 }, title: { display: true, color: ink, font: { family: css('--f-body'), size: 12 } } },
    },
  };
}
function mk(id, cfg) {
  if (!window.Chart) return null;
  if (charts[id]) charts[id].destroy();
  cfg.plugins = [cursorPlugin];
  charts[id] = new Chart(document.getElementById(id), cfg);
  return charts[id];
}
function noChart() {
  $$('.chart-box, .sweep-box').forEach((b) => { if (!b.querySelector('p')) b.innerHTML = '<p class="small">Grafik kütüphanesi yüklenemedi. İnternet bağlantısını kontrol edin.</p>'; });
}
function redrawCursor() { Object.values(charts).forEach((c) => c && c.draw()); }
function renderCharts() {
  if (!window.Chart) { noChart(); return; }
  Chart.defaults.locale = 'tr-TR';
  dirty.charts = false;
  const r = S.cur, ref = S.ref;
  const pEnd = (r.tAirEnd || 0) + 0.002;
  const tMax = Math.ceil(Math.max(r.tAirEnd || 0.1, ref.tAirEnd || 0) * 1.08 * 1000 / 25) * 25;
  const pw = (k, f) => r.rec.filter((e) => e.t <= pEnd && (k === null || e.ph === k || (k === 2 && e.ph === 3))).map((e) => ({ x: e.t * 1000, y: f(e) }));
  // thrust
  const o1 = chartBase(1000); o1.scales.x.title.text = 'zaman (ms)'; o1.scales.y.title.text = 'itki (N)'; o1.scales.x.min = 0; o1.scales.x.max = tMax; o1.scales.y.min = 0;
  const ds1 = [];
  const tube = pw(0, (e) => e.T); if (tube.length) ds1.push({ label: 'Tüp', data: tube, borderColor: css('--tube'), backgroundColor: css('--panel-2'), fill: 'origin', pointRadius: 0, borderWidth: 2 });
  ds1.push({ label: 'Su evresi', data: pw(1, (e) => e.T), borderColor: css('--water'), backgroundColor: css('--water-soft'), fill: 'origin', pointRadius: 0, borderWidth: 2 });
  const air = pw(2, (e) => e.T); if (air.length && r.tWaterEnd) air.unshift({ x: r.tWaterEnd * 1000, y: 0 }, { x: r.tWaterEnd * 1000, y: air[0].y });
  ds1.push({ label: 'Hava evresi', data: air, borderColor: css('--air'), backgroundColor: css('--warn-bg'), fill: 'origin', pointRadius: 0, borderWidth: 2 });
  mk('ch-thrust', { type: 'line', data: { datasets: ds1 }, options: o1 });
  // pressure
  const o2 = chartBase(1000); o2.scales.x.title.text = 'zaman (ms)'; o2.scales.y.title.text = 'basınç (bar, gösterge)'; o2.scales.x.min = 0; o2.scales.x.max = tMax; o2.scales.y.min = 0;
  const prs = (rr) => rr.rec.filter((e) => e.t <= (rr.tAirEnd || 0) + 0.002).map((e) => ({ x: e.t * 1000, y: e.P }));
  mk('ch-press', { type: 'line', data: { datasets: [
    { label: 'Bu tasarım', data: prs(r), borderColor: css('--accent'), pointRadius: 0, borderWidth: 2 },
    { label: 'Karşılaştırma', data: prs(ref), borderColor: css('--muted'), borderDash: [5, 5], pointRadius: 0, borderWidth: 2 },
  ] }, options: o2 });
  const tAll = Math.ceil(Math.max(r.time, ref.time) * 2) / 2;
  const series = (rr, f, every = 2) => rr.rec.filter((e, i) => i % every === 0).map((e) => ({ x: e.t, y: f(e) }));
  const two = (f, label) => [
    { label: 'Bu tasarım', data: series(r, f), borderColor: css('--accent'), pointRadius: 0, borderWidth: 2 },
    { label: 'Karşılaştırma', data: series(ref, f), borderColor: css('--muted'), borderDash: [5, 5], pointRadius: 0, borderWidth: 2 },
  ];
  const o3 = chartBase(1); o3.scales.x.title.text = 'zaman (s)'; o3.scales.y.title.text = 'irtifa (m)'; o3.scales.x.min = 0; o3.scales.x.max = tAll; o3.scales.y.min = 0;
  mk('ch-alt', { type: 'line', data: { datasets: two((e) => Math.max(0, e.y)) }, options: o3 });
  const o4 = chartBase(1); o4.scales.x.title.text = 'zaman (s)'; o4.scales.y.title.text = 'hız (m/s)'; o4.scales.x.min = 0; o4.scales.x.max = tAll; o4.scales.y.min = 0;
  mk('ch-vel', { type: 'line', data: { datasets: two((e) => e.v) }, options: o4 });
  const o5 = chartBase(1); o5.scales.x.title.text = 'zaman (s)'; o5.scales.y.title.text = 'yana kayma (%)'; o5.scales.x.min = 0; o5.scales.x.max = tAll; o5.scales.y.min = -100; o5.scales.y.max = 100; o5.scales.y.ticks.stepSize = 50;
  const sl = (rr) => rr.rec.filter((e) => !e.chute).map((e) => ({ x: e.t, y: e.s * 100 }));
  const slz = r.rec.filter((e) => !e.chute).map((e) => ({ x: e.t, y: e.sz * 100 }));
  const hasZ = slz.some((q) => Math.abs(q.y) > 1);
  mk('ch-slosh', { type: 'line', data: { datasets: [
    { label: 'Bu tasarım', data: sl(r), borderColor: css('--water'), pointRadius: 0, borderWidth: 1.6 },
    { label: 'Karşılaştırma', data: sl(ref), borderColor: css('--muted'), borderDash: [4, 4], pointRadius: 0, borderWidth: 1.6 },
  ].concat(hasZ ? [{ label: 'Yan düzlem', data: slz, borderColor: css('--accent-2'), pointRadius: 0, borderWidth: 1.4 }] : []) }, options: o5 });
  const o6 = chartBase(1); o6.scales.x.title.text = 'zaman (s)'; o6.scales.y.title.text = 'burun açısı, yataydan (°)'; o6.scales.x.min = 0; o6.scales.x.max = tAll; o6.scales.y.min = -90; o6.scales.y.max = 90; o6.scales.y.ticks.stepSize = 45;
  const pit = (rr) => rr.rec.filter((e, i) => !e.chute && i % 2 === 0).map((e) => ({ x: e.t, y: Math.asin(clamp(e.by, -1, 1)) * 180 / Math.PI }));
  mk('ch-pitch', { type: 'line', data: { datasets: [
    { label: 'Bu tasarım', data: pit(r), borderColor: css('--accent'), pointRadius: 0, borderWidth: 2, spanGaps: false },
    { label: 'Karşılaştırma', data: pit(ref), borderColor: css('--muted'), borderDash: [5, 5], pointRadius: 0, borderWidth: 2, spanGaps: false },
  ] }, options: o6 });
  const pwr = (rr, f) => rr.rec.filter((e) => e.t <= (rr.tAirEnd || 0) + 0.004).map((e) => ({ x: e.t * 1000, y: f(e, rr) }));
  const twoP = (f) => [
    { label: 'Bu tasarım', data: pwr(r, f), borderColor: css('--accent'), pointRadius: 0, borderWidth: 2 },
    { label: 'Karşılaştırma', data: pwr(ref, f), borderColor: css('--muted'), borderDash: [5, 5], pointRadius: 0, borderWidth: 2 },
  ];
  const o7 = chartBase(1000); o7.scales.x.title.text = 'zaman (ms)'; o7.scales.y.title.text = 'eksenel ivme (g)'; o7.scales.x.min = 0; o7.scales.x.max = tMax;
  mk('ch-acc', { type: 'line', data: { datasets: twoP((e) => e.ax) }, options: o7 });
  const o8 = chartBase(1000); o8.scales.x.title.text = 'zaman (ms)'; o8.scales.y.title.text = 'kalan itki suyu (ml)'; o8.scales.x.min = 0; o8.scales.x.max = tMax; o8.scales.y.min = 0;
  mk('ch-water', { type: 'line', data: { datasets: twoP((e) => e.Vw) }, options: o8 });
  const o9 = chartBase(1000); o9.scales.x.title.text = 'zaman (ms)'; o9.scales.y.title.text = 'v_su (m/s)'; o9.scales.x.min = 0; o9.scales.x.max = tMax; o9.scales.y.min = 0;
  mk('ch-ve', { type: 'line', data: { datasets: twoP((e, rr) => (e.ph === 1 ? rr.G.cd * Math.sqrt(2 * Math.max(0, e.P) * 1e5 / C.rhoW) : null)) }, options: o9 });
}

// ------------------------------------------------------------------ sweep
function range(a, b, s) { const o = []; for (let v = a; v <= b + 1e-9; v += s) o.push(Math.round(v * 1000) / 1000); return o; }
const SWEEPS = {
  wf: { label: 'İtki suyu (haznenin %\'si)', vals: () => range(5, 60, 2.5), set: (p, v) => (p.wf = v / 100), cur: (p) => p.wf * 100, fmt: (v, b) => '%' + nf(v, 1) + (b ? ' · ' + nf(v / 100 * b.Vc, 2) + ' L' : '') },
  P: { label: 'Başlangıç basıncı (bar)', vals: () => range(2, 12, 0.5), set: (p, v) => (p.P = v), cur: (p) => p.P, fmt: (v) => nf(v, 1) + ' bar' },
  dry: { label: 'Boş yapı kütlesi (g)', vals: () => range(200, 1500, 50), set: (p, v) => (p.dry = v / 1000), cur: (p) => p.dry * 1000, fmt: (v) => nf(v, 0) + ' g' },
  payload: { label: 'Çalkalanan su (L)', vals: () => range(0.25, 3, 0.125), set: (p, v) => (p.payload = v), cur: (p) => p.payload, fmt: (v) => nf(v, 2) + ' L' },
  Vc: { label: 'Basınç haznesi (L), itki suyu aynı', vals: () => range(1, 6, 0.25), set: (p, v) => { const w = p.Vc * p.wf; p.Vc = v; p.wf = Math.min(w, 0.95 * v) / v; }, cur: (p) => p.Vc, fmt: (v) => nf(v, 2) + ' L' },
  nozzle: { label: 'Nozul çapı (mm)', vals: () => range(6, 30, 1), set: (p, v) => (p.nozzle = v), cur: (p) => p.nozzle, fmt: (v) => nf(v, 0) + ' mm' },
  cd: { label: 'Nozul kayıp faktörü C_d', vals: () => range(0.8, 1, 0.01), set: (p, v) => (p.cd = v), cur: (p) => p.cd, fmt: (v) => nf(v, 2) },
  Dmm: { label: 'Şişe çapı (mm)', vals: () => range(70, 150, 5), set: (p, v) => (p.Dmm = v), cur: (p) => p.Dmm, fmt: (v) => nf(v, 0) + ' mm' },
  CD: { label: 'Hava sürtünme katsayısı C_D', vals: () => range(0.2, 0.9, 0.05), set: (p, v) => (p.CD = v), cur: (p) => p.CD, fmt: (v) => nf(v, 2) },
  wspan: { label: 'Kanat açıklığı (m)', vals: () => range(0.4, 2.6, 0.1), set: (p, v) => { p.wspan = v; if (p.wing !== 'sabit') p.wing = 'acilir'; }, cur: (p) => p.wspan, fmt: (v) => nf(v, 1) + ' m' },
  wchord: { label: 'Kanat genişliği, veter (cm)', vals: () => range(6, 35, 1), set: (p, v) => { p.wchord = v / 100; if (p.wing !== 'sabit') p.wing = 'acilir'; }, cur: (p) => p.wchord * 100, fmt: (v) => nf(v, 0) + ' cm' },
  winc: { label: 'Kanat montaj açısı (°)', vals: () => range(-2, 12, 0.5), set: (p, v) => { p.winc = v; if (p.wing !== 'sabit') p.wing = 'acilir'; }, cur: (p) => p.winc, fmt: (v) => nf(v, 1) + '°' },
  baffle: { label: 'Bölme levhası oranı (%)', vals: () => range(0, 80, 5), set: (p, v) => (p.baffle = v / 100), cur: (p) => p.baffle * 100, fmt: (v) => '%' + nf(v, 0) },
  wind: { label: 'Rüzgâr (m/s)', vals: () => range(-6, 6, 0.5), set: (p, v) => (p.wind = v), cur: (p) => p.wind, fmt: (v) => nf(v, 1) + ' m/s' },
  windZ: { label: 'Yan rüzgâr (m/s)', vals: () => range(-6, 6, 0.5), set: (p, v) => (p.windZ = v), cur: (p) => p.windZ, fmt: (v) => nf(v, 1) + ' m/s' },
  angle: { label: 'Fırlatma açısı (°)', vals: () => range(80, 90, 0.5), set: (p, v) => (p.angle = v), cur: (p) => p.angle, fmt: (v) => nf(v, 1) + '°' },
};
function renderSweep() {
  if (S.tab !== 'tarama') return;
  if (!window.Chart) { noChart(); return; }
  const key = $('#sweep-var').value, SW = SWEEPS[key], base = S.cur.p;
  const ok = [], bad = [], nofly = [];
  let best = null;
  const vals = SW.vals();
  $('#sweep-best').textContent = 'Hesaplanıyor…';
  sliced(vals, (v, i) => {
    const p = Object.assign({}, base); SW.set(p, v);
    const r = simulate(p, false);
    const pt = { x: v, y: r.score };
    (!r.left ? nofly : rulesOk(p, r) ? ok : bad).push(pt);
    if (allOk(p, r) && (!best || r.score > best.y)) best = pt;
    if (i % 4 === 0) $('#sweep-best').textContent = `Hesaplanıyor… ${i} / ${vals.length}`;
  }, () => drawSweep(SW, base, ok, bad, nofly, best));
}
function drawSweep(SW, base, ok, bad, nofly, best) {
  if (S.tab !== 'tarama') return;
  const o = chartBase(1);
  o.plugins.cursor = { on: false };
  o.scales.x.title.text = SW.label; o.scales.y.title.text = 'uçuş puanı'; o.scales.y.min = 0;
  o.interaction = { mode: 'nearest', intersect: true };
  o.plugins.legend.labels.usePointStyle = true; o.plugins.legend.labels.boxWidth = 10; o.plugins.legend.labels.boxHeight = 10;
  o.plugins.tooltip.callbacks = { title: () => '', label: (c) => `${SW.fmt(c.parsed.x, base)} → ${nf(c.parsed.y, 2)} puan` };
  o.onClick = (e, els, chart) => {
    if (!els.length) return;
    const pt = chart.data.datasets[els[0].datasetIndex].data[els[0].index];
    if (!pt) return;
    const q = Object.assign({}, base); SW.set(q, pt.x);
    writeControls(toCtl(q)); recompute(); showTab('tarama');
  };
  const curX = SW.cur(base);
  const ds = [
    { type: 'line', label: 'Kurallara uygun', data: ok, borderColor: css('--accent'), backgroundColor: css('--accent'), pointRadius: 3.5, pointHoverRadius: 6, borderWidth: 2, showLine: true },
    { type: 'scatter', label: 'Kural dışı', data: bad, borderColor: css('--bad'), backgroundColor: css('--bad'), pointStyle: 'crossRot', pointRadius: 6, pointHoverRadius: 8, borderWidth: 2 },
  ];
  if (nofly.length) ds.push({ type: 'scatter', label: 'Rampadan ayrılmıyor', data: nofly, borderColor: css('--warn'), backgroundColor: css('--warn'), pointStyle: 'triangle', pointRadius: 5, pointHoverRadius: 7, borderWidth: 1 });
  if (best) ds.push({ type: 'scatter', label: 'En iyi uygun nokta', data: [best], borderColor: css('--good'), backgroundColor: css('--panel'), pointRadius: 8, pointHoverRadius: 9, borderWidth: 3 });
  const curP = Object.assign({}, base); SW.set(curP, curX);
  const curY = (curP.wing === base.wing) ? S.cur.score : simulate(curP, false).score;
  ds.push({ type: 'scatter', label: curP.wing === base.wing ? 'Şu anki ayar' : 'Şu anki ayar (açılır kanatla)', data: [{ x: curX, y: curY }], borderColor: css('--ink'), backgroundColor: css('--ink'), pointStyle: 'rectRot', pointRadius: 6, pointHoverRadius: 8 });
  mk('ch-sweep', { type: 'line', data: { datasets: ds }, options: o });
  $('#sweep-best').textContent = best ? `En iyi uygun nokta: ${SW.fmt(best.x, base)} → ${nf(best.y, 2)} puan` : 'Bu aralıkta kurallara uygun ve uçan nokta yok.';
  dirty.sweep = false;
}

// ------------------------------------------------------------------ rocket tab
function placeLabels(items, obstacles) {
  const rows = [26, 74], placed = [];
  const width = (it) => Math.max(it.t1.length * 8.3, (it.t2 || '').length * 6.9);
  items.slice().sort((a, b) => a.x - b.x).forEach((it) => {
    const w = width(it);
    let done = null;
    for (const row of [0, 1]) {
      for (const anchor of ['start', 'end']) {
        const x0 = anchor === 'start' ? it.x + 12 : it.x - 12 - w, x1 = x0 + w;
        if (x0 < 4 || x1 > 996) continue;
        if (placed.some((q) => q.row === row && !(x1 < q.x0 - 12 || x0 > q.x1 + 12))) continue;
        if (row === 1 && placed.some((q) => q.row === 0 && q.lx > x0 - 6 && q.lx < x1 + 6)) continue;
        if (row === 1 && obstacles.some((o) => !(x1 < o[0] || x0 > o[1]) || (it.x > o[0] - 4 && it.x < o[1] + 4))) continue;
        if (row === 0 && placed.some((q) => q.row === 1 && it.x > q.x0 - 6 && it.x < q.x1 + 6)) continue;
        done = { row, anchor, x0, x1, lx: it.x, it }; break;
      }
      if (done) break;
    }
    if (!done) { const x0 = Math.min(it.x + 12, 996 - w); done = { row: 0, anchor: 'start', x0, x1: x0 + w, lx: it.x, it }; }
    placed.push(done);
  });
  return placed.map((q) => Object.assign({}, q, { y: rows[q.row] }));
}
function drawAnatomy(p, G) {
  const s = Math.min(586, 905 / Math.max(G.L, 1.5)), x0 = 60, cy = 168, R = G.R * s;
  const nose = x0 + G.L * s;
  const X = (xn) => nose - xn * s;
  const top = cy - R, bot = cy + R;
  const parts = [];
  const fs = G.fin.span * s, fc = G.fin.chord * s;
  const fx1 = X(G.L - 0.005), fx0 = fx1 + fc;
  parts.push(`<path class="svg-fin" d="M${fx0},${top} L${fx1},${top} L${fx1 - 4},${top - fs} L${fx1 + fc * 0.45},${top - fs} Z"/>`);
  parts.push(`<path class="svg-fin" d="M${fx0},${bot} L${fx1},${bot} L${fx1 - 4},${bot + fs} L${fx1 + fc * 0.45},${bot + fs} Z"/>`);
  const nx = X(G.L);
  parts.push(`<path class="svg-nozzle" d="M${nx + 14},${cy - 13} L${nx},${cy - 9} L${nx},${cy + 9} L${nx + 14},${cy + 13} Z"/>`);
  parts.push(`<rect class="svg-body" x="${X(G.L) + 10}" y="${top + 8}" width="${X(G.Ln + G.Lp + G.Lc) - X(G.L) - 10}" height="${2 * R - 16}" rx="4"/>`);
  const c0 = X(G.Ln + G.Lp + G.Lc), c1 = X(G.Ln + G.Lp);
  parts.push(`<rect class="svg-body" x="${c0}" y="${top}" width="${c1 - c0}" height="${2 * R}" rx="10"/>`);
  const wl = (G.Lc * p.wf) * s;
  parts.push(`<rect class="svg-water" x="${c0 + 1.5}" y="${top + 1.5}" width="${Math.max(0, wl - 1.5)}" height="${2 * R - 3}" rx="9"/>`);
  parts.push(`<line class="svg-water-line" x1="${c0 + wl}" y1="${top + 3}" x2="${c0 + wl}" y2="${bot - 3}"/>`);
  const p0 = X(G.Ln + G.Lp), p1 = X(G.Ln);
  parts.push(`<rect class="svg-body" x="${p0}" y="${top}" width="${p1 - p0}" height="${2 * R}" rx="10"/>`);
  let wave = `M${p0 + 2},${cy}`;
  for (let i = 0; i <= 24; i++) { const xx = p0 + 2 + (p1 - p0 - 4) * i / 24; wave += ` L${xx.toFixed(1)},${(cy + Math.sin(i / 24 * Math.PI * 3) * 3).toFixed(1)}`; }
  parts.push(`<path class="svg-water" d="${wave} L${p1 - 2},${bot - 2} L${p0 + 2},${bot - 2} Z"/>`);
  parts.push(`<path class="svg-water-line" d="${wave}"/>`);
  if (p.baffle > 0) {
    const stub = R * (1 - Math.sqrt(Math.max(0, 1 - p.baffle)));
    [1 / 3, 2 / 3].forEach((f) => {
      const bx = p0 + (p1 - p0) * f;
      parts.push(`<line class="svg-baffle" x1="${bx}" y1="${top + 2}" x2="${bx}" y2="${top + 2 + stub}"/>`);
      parts.push(`<line class="svg-baffle" x1="${bx}" y1="${bot - 2}" x2="${bx}" y2="${bot - 2 - stub}"/>`);
    });
  }
  const n0 = X(G.Ln);
  parts.push(`<path class="svg-body" d="M${n0},${top} C${n0 + 40},${top} ${nose - 12},${cy - 12} ${nose},${cy} C${nose - 12},${cy + 12} ${n0 + 40},${bot} ${n0},${bot} Z"/>`);
  if (G.wing) {
    // root section of a high wing, sitting on top of the body (see the top view below for the planform)
    const le = X(G.xW - 0.25 * G.chord), te = X(G.xW + 0.75 * G.chord), th = Math.max(5, 0.12 * G.chord * s), wy = top - 1;
    parts.push(`<path class="svg-wing" d="M${le},${wy} C${le - 6},${wy - th} ${te + 20},${wy - th * 0.6} ${te},${wy} C${te + 20},${wy + th * 0.15} ${le - 6},${wy + th * 0.3} ${le},${wy} Z"/>`);
  }
  const items = [
    { x: (n0 + nose) / 2, t1: 'Burun', t2: 'altimetre kapsülü' },
    { x: (p0 + p1) / 2, t1: 'Çalkalanma tankı', t2: `${nf(p.payload, 2)} L su · bölme %${nf(p.baffle * 100, 0)}` },
    { x: (c0 + c1) / 2, t1: 'Basınç haznesi', t2: `${nf(p.Vc, 1)} L · ${nf(p.Vc * p.wf, 2)} L itki suyu` },
    { x: fx1 + fc * 0.3, t1: 'Kanatçıklar', t2: `nozul ${nf(p.nozzle, p.nozzle % 1 ? 1 : 0)} mm · Rectus 26` },
  ];
  placeLabels(items, [[fx1 - 8, fx0 + 8]]).forEach((q) => {
    const tx = q.anchor === 'start' ? q.x0 : q.x1;
    const yTop = (q.it.t1 === 'Kanatçıklar') ? top - fs - 4 : top - 4;
    parts.push(`<path class="svg-lead" d="M${q.lx},${yTop} L${q.lx},${q.y + 8} L${q.anchor === 'start' ? q.x0 - 4 : q.x1 + 4},${q.y + 8}"/>`);
    parts.push(`<text class="svg-label" x="${tx}" y="${q.y}" text-anchor="${q.anchor}">${q.it.t1}</text>`);
    parts.push(`<text class="svg-sub" x="${tx}" y="${q.y + 24}" text-anchor="${q.anchor}">${q.it.t2}</text>`);
  });
  const cgx = X(G.xcg0), cg1x = X(G.xcg1), cpx = X(G.xNP0), my = bot + Math.max(30, fs + 14);
  [cgx, cpx].forEach((xx) => parts.push(`<line class="svg-dim" x1="${xx}" y1="${bot + 4}" x2="${xx}" y2="${my - 8}"/>`));
  parts.push(`<circle class="svg-cg" cx="${cgx}" cy="${my}" r="7"/>`);
  parts.push(`<circle cx="${cg1x}" cy="${my}" r="5" fill="none" stroke="var(--ink)" stroke-width="1.5" stroke-dasharray="2 2"/>`);
  parts.push(`<circle class="svg-cp" cx="${cpx}" cy="${my}" r="7"/>`);
  // marker labels under the markers, in up to two rows so they never overlap
  const labs = [{ x: cgx, t: 'AM dolu' }, { x: cg1x, t: 'AM boş' }, { x: cpx, t: G.wing === 'sabit' ? 'BM (kanatla)' : G.wing ? 'BM (kuyrukla)' : 'BM' }].sort((a, c) => a.x - c.x);
  const rowsEnd = [-1e9, -1e9]; let used = 1;
  labs.forEach((q) => {
    const w = q.t.length * 7.9, x0 = q.x - w / 2;
    const row = x0 > rowsEnd[0] + 10 ? 0 : 1; used = Math.max(used, row + 1);
    rowsEnd[row] = Math.max(rowsEnd[row], x0 + w);
    parts.push(`<text class="svg-dim-t" x="${q.x}" y="${my + 26 + row * 18}" text-anchor="middle">${q.t}</text>`);
  });
  const dy = Math.max(318, my + 30 + used * 18 + 12);
  $('#anatomy').setAttribute('viewBox', `0 0 1000 ${Math.round(dy + 26)}`);
  parts.push(`<path class="svg-dim" d="M${X(G.L)},${dy - 6} V${dy + 6} M${X(G.L)},${dy} H${nose} M${nose},${dy - 6} V${dy + 6}"/>`);
  parts.push(`<text class="svg-dim-t" x="${(X(G.L) + nose) / 2}" y="${dy - 8}" text-anchor="middle">boy ${nf(G.L, 2)} m</text>`);
  const lim = x0 + 1.5 * s;
  parts.push(`<line class="svg-limit" x1="${lim}" y1="${top - 10}" x2="${lim}" y2="${dy + 8}"/>`);
  parts.push(`<text class="svg-limit-t" x="${lim - 6}" y="${dy + 16}" text-anchor="end">sınır 1,5 m</text>`);
  $('#anatomy').innerHTML = parts.join('');
  drawPlanform(p, G);
}
// top view: nose up, span across; wing and tail to scale, with the 2,25 m span limit
function drawPlanform(p, G) {
  const s = 330, cx = 500, y0 = 46, Rb = Math.max(G.R * s, 6), Lpx = G.L * s;
  const Y = (xn) => y0 + xn * s;
  const yEnd = Y(G.L), H = Math.round(yEnd + 92);
  const svg = $('#planform'); if (!svg) return;
  svg.setAttribute('viewBox', `0 0 1000 ${H}`);
  const out = [];
  const lim = 1.125 * s;
  [-1, 1].forEach((sd) => out.push(`<line class="svg-limit" x1="${cx + sd * lim}" y1="20" x2="${cx + sd * lim}" y2="${yEnd + 34}"/>`));
  out.push(`<text class="svg-limit-t" x="${cx + lim - 8}" y="30" text-anchor="end">kanat açıklığı sınırı 2,25 m</text>`);
  // fins seen from above (the horizontal pair)
  const fs = G.fin.span * s, fc = G.fin.chord * s, ft = Y(G.L - 0.005);
  [-1, 1].forEach((sd) => out.push(`<path class="svg-fin" d="M${cx + sd * Rb},${ft - fc} L${cx + sd * Rb},${ft} L${cx + sd * (Rb + fs)},${ft + 4} L${cx + sd * (Rb + fs)},${ft - fc * 0.55} Z"/>`));
  if (G.wing) {
    // horizontal tail (aspect ratio 4) at the tail arm, and the main wing at the empty CG
    const tb = Math.sqrt(4 * G.St) * s, tcv = Math.sqrt(G.St / 4) * s, tLe = Y(G.xT) - 0.25 * tcv;
    [-1, 1].forEach((sd) => out.push(`<path class="svg-tail" d="M${cx + sd * Rb},${tLe} L${cx + sd * tb / 2},${tLe + tcv * 0.15} L${cx + sd * tb / 2},${tLe + tcv * 0.9} L${cx + sd * Rb},${tLe + tcv} Z"/>`));
    const c = G.chord * s, le = Y(G.xW) - 0.25 * c, hb = p.wspan / 2 * s;
    [-1, 1].forEach((sd) => out.push(`<path class="svg-wing${hb > lim ? ' over' : ''}" d="M${cx + sd * Rb},${le} L${cx + sd * hb},${le + 0.1 * c} L${cx + sd * hb},${le + 0.85 * c} L${cx + sd * Rb},${le + c} Z"/>`));
    if (G.wing === 'acilir') {
      const pack = Math.min(p.wspan / 2, 0.4 * G.L) * s, pw = Math.min(c, 2 * G.R * s) * 0.6;
      [-1, 1].forEach((sd) => out.push(`<rect class="svg-fold" x="${sd > 0 ? cx + Rb : cx - Rb - pw}" y="${le}" width="${pw}" height="${pack}" rx="3"/>`));
    }
    // span and chord dimensions
    const dy = yEnd + 54;
    out.push(`<path class="svg-dim" d="M${cx - hb},${dy - 6} V${dy + 6} M${cx - hb},${dy} H${cx + hb} M${cx + hb},${dy - 6} V${dy + 6}"/>`);
    out.push(`<text class="svg-dim-t" x="${cx}" y="${dy - 8}" text-anchor="middle">kanat açıklığı b = ${nf(p.wspan, 2)} m</text>`);
    const tx = cx + hb + 10 > 960 ? cx - hb - 10 : cx + hb + 10;
    out.push(`<text class="svg-dim-t" x="${tx}" y="${le + c / 2 + 4}" text-anchor="${tx > cx ? 'start' : 'end'}">c = ${nf(p.wchord * 100, 0)} cm</text>`);
  }
  // body and nose
  const n1 = Y(G.Ln);
  out.push(`<path class="svg-body" d="M${cx - Rb},${n1} C${cx - Rb},${y0 + 24} ${cx - 8},${y0 + 6} ${cx},${y0} C${cx + 8},${y0 + 6} ${cx + Rb},${y0 + 24} ${cx + Rb},${n1} Z"/>`);
  out.push(`<rect class="svg-body" x="${cx - Rb}" y="${n1}" width="${2 * Rb}" height="${yEnd - n1}" rx="3"/>`);
  out.push(`<rect class="svg-nozzle" x="${cx - Rb * 0.4}" y="${yEnd}" width="${Rb * 0.8}" height="7"/>`);
  out.push(`<circle class="svg-cg" cx="${cx}" cy="${Y(G.xcg1)}" r="5"/>`);
  if (G.wing) {
    const c = G.chord * s, le = Y(G.xW) - 0.25 * c;
    out.push(`<text class="svg-label" x="${cx + lim - 8}" y="${le - 30}" text-anchor="end">${G.wing === 'acilir' ? 'Açılır kanat' : 'Sabit kanat'}</text>`);
    out.push(`<text class="svg-sub" x="${cx + lim - 8}" y="${le - 12}" text-anchor="end">boş roketin ağırlık merkezi hizasında</text>`);
    if (G.wing === 'acilir') out.push(`<text class="svg-sub" x="${cx - Rb - Math.min(c, 2 * G.R * s) * 0.6 - 10}" y="${le + c + 22}" text-anchor="end">kesikli: katlı hâli (tırmanışta)</text>`);
    const tY = Y(G.xT);
    out.push(`<text class="svg-sub" x="${cx - Math.sqrt(4 * G.St) * s / 2 - 10}" y="${tY + 4}" text-anchor="end">yatay kuyruk ${nf(G.St * 1e4, 0)} cm²</text>`);
  } else {
    out.push(`<text class="svg-sub" x="${cx + 60}" y="${Y(G.L * 0.45)}">Kanat yok: yalnızca kanatçıklar.</text>`);
    out.push(`<text class="svg-sub" x="${cx + 60}" y="${Y(G.L * 0.45) + 20}">Simülasyon panelinde "Kanat tipi"ni Sabit veya Açılır yapın.</text>`);
  }
  svg.innerHTML = out.join('');
  const info = $('#plan-info');
  if (info) info.innerHTML = G.wing
    ? `<span>Kanat alanı <b>${nf(G.Sw, 3)} m²</b></span><span>En-boy oranı <b>${nf(G.AR, 1)}</b></span><span>Yatay kuyruk <b>${nf(G.St * 1e4, 0)} cm²</b></span><span>Kanat + kuyruk${G.wing === 'acilir' ? ' + mekanizma' : ''} <b>${nf(G.mWing * 1000, 0)} g</b></span><span>Montaj açısı <b>${nf(p.winc, 1)}°</b></span>${G.wing === 'acilir' ? '<span>Açılma <b>tepe + 0,15 s</b></span>' : ''}`
    : `<span>Kanat yok</span>`;
}
function renderMass() {
  const r = S.cur, G = r.G, p = r.p;
  const seg = [
    ['Yapı', p.dry, 'var(--c-struct)'],
    ['Kanatçıklar', G.fin.mass, 'var(--c-fin)'],
    ['Çalkalanan su', G.mPay, 'var(--c-pay)'],
    ['İtki suyu', G.mProp, 'var(--c-prop)'],
  ];
  if (G.wing) seg.splice(2, 0, [G.wing === 'acilir' ? 'Kanat, kuyruk, açma mek.' : 'Kanat ve kuyruk', G.mWing, 'var(--accent-2)']);
  const scale = Math.max(5.5, G.TOW * 1.05);
  $('#mass').innerHTML = `<div class="mass-bar">${seg.map(([k, v, c]) => `<i style="width:${(v / scale * 100).toFixed(2)}%;background:${c}" title="${k}: ${nf(v, 2)} kg"></i>`).join('')}</div>
    <div class="mass-limit" style="left:${(5 / scale * 100).toFixed(2)}%"><span>5 kg sınırı</span></div>`;
  $('#mass-legend').innerHTML = seg.map(([k, v, c]) => `<span><i style="background:${c}"></i>${k} <b>${nf(v, 2)} kg</b></span>`).join('') + `<span>Toplam <b>${nf(G.TOW, 2)} kg</b></span>`;
  const chipSM = (v) => v < 0 ? ['no', 'kararsız'] : v < 1 ? ['no', 'zayıf'] : v <= 2.5 ? ['ok', 'iyi'] : ['mid', 'çok kararlı'];
  const chipTW = (v) => v >= 5 ? ['ok', 'yeterli'] : v >= 3 ? ['mid', 'sınırda'] : ['no', 'düşük'];
  const sm0 = chipSM(G.sm0), sm1 = chipSM(G.sm1), tw = chipTW(G.tw);
  const fn = Math.sqrt((1.841 * C.g / G.R) * Math.tanh(1.841 * ((p.payload / 1000) / (Math.PI * G.R * G.R)) / G.R)) / (2 * Math.PI);
  const stats = [
    ['Kararlılık payı (kalkış)', nf(G.sm0, 1) + ' kalibre', sm0, 'AM ile BM arası, gövde çapı cinsinden. 1–2 kalibre yaygın hedef.' + (G.wing === 'sabit' ? ' BM hesabına açık kanat ve kuyruk dahil.' : G.wing ? ' BM hesabına yatay kuyruk dahil (kanat katlı).' : '') + (G.sm0 < 0 ? ' Negatif: AM, BM\'nin arkasında; roket kalkışta kararsız.' : '')],
    ['Kararlılık payı (su bitince)', nf(G.sm1, 1) + ' kalibre', sm1, 'İtki suyu bitince AM öne kayar, pay büyür.'],
    ['İtki / ağırlık (kalkış)', nf(G.tw, 1) + ' : 1', tw, 'Sunumdaki kural: ~5:1, en az 3:1.'],
    p.tube ? ['Kalkıştaki itki', nf(G.T0, 0) + ' N', null, `Fırlatma tüpü evresi: F = ΔP · A_tüp. Su evresi tüpten çıkınca ${nf(G.Tw0 * (r.Pstart - C.Patm) / (p.P * 1e5), 0)} N ile başlar.`]
      : ['Kalkıştaki itki', nf(G.T0, 0) + ' N', null, 'F = ṁ · v_su = 2 · C_d² · A_n · ΔP; nozul alanı ve basınçla doğru orantılı.'],
    ['Çalkalanma doğal frekansı (yerde)', nf(fn, 1) + ' Hz', null, 'İtkide ivmeyle artar, süzülmede düşer.'],
    ['Bölme levhası sönümü', 'ζ = ' + nf(r.zeta, 2), null, 'Levhalar suyun salınım enerjisini emer; hareket eden kütle ' + nf(r.mS * 1000, 0) + ' g.'],
  ];
  if (G.wing === 'acilir') stats.splice(2, 0, ['Kararlılık payı (kanat açık)', nf(G.smW, 1) + ' kalibre', chipSM(G.smW), 'Süzülmede: su bitmiş, kanat ve kuyruk açık. Kanat AM\'ye yakın olduğu için pay küçülür.']);
  $('#stats').innerHTML = stats.map(([k, v, ch, d]) => `<div class="stat"><span class="k">${k}</span><span class="v">${v} ${ch ? `<span class="chip ${ch[0]}">${ch[1]}</span>` : ''}</span><p>${d}</p></div>`).join('');
}

// ================================================================== actions
function applyParams(v, opts = {}) {
  writeControls(Object.assign({}, DEFAULTS, v));
  recompute({ t: opts.t });
}
function resetAll() {
  pause();
  applyParams({}, { t: 0 });
  S.ref = S.cur; S.refName = 'İKARUS-1 (örnek)';
  $('#sweep-var').value = 'wf'; $('#speed').value = '1'; $('#slowmo').checked = true;
  renderToolbar();
  dirty.charts = true;
  showTab('ucus');
  setT(0); setPlayUI();
}

function wire() {
  $('#nav-konu').addEventListener('click', () => showView('konu'));
  $('#nav-sim').addEventListener('click', () => showView('sim'));
  $$('.subtabs button').forEach((b) => b.addEventListener('click', () => showTab(b.dataset.tab)));
  $('#controls').addEventListener('input', (e) => {
    const t = e.target;
    if (t.type === 'range') { const k = t.id.slice(2); setVal(k, +t.value, 'range'); paintRange(t); recompute(); }
    else if (t.classList.contains('num')) { clearTimeout(t._tm); t._tm = setTimeout(() => commitNum(t, false), 600); }
  });
  $('#controls').addEventListener('change', (e) => {
    const t = e.target;
    if (t.type === 'radio' || t.type === 'checkbox') recompute();
    else if (t.classList.contains('num')) { clearTimeout(t._tm); commitNum(t, true); }
  });
  $('#controls').addEventListener('keydown', (e) => {
    const el = e.target; if (!el.classList || !el.classList.contains('num')) return;
    if (e.key === 'Enter') { e.preventDefault(); clearTimeout(el._tm); commitNum(el, true); return; }
    if (e.key === 'ArrowUp' || e.key === 'ArrowDown') {
      e.preventDefault();
      const k = el.id.slice(2), rg = $('#c-' + k), st = rg ? +rg.step || 0.1 : 0.1;
      setVal(k, VAL[k] + (e.key === 'ArrowUp' ? st : -st) * (e.shiftKey ? 10 : 1), 'num'); showNum(k, true); recompute();
    }
  });
  $('#controls').addEventListener('focusout', (e) => { if (e.target.classList && e.target.classList.contains('num')) showNum(e.target.id.slice(2)); });
  $$('[data-preset]').forEach((b) => b.addEventListener('click', () => {
    const live = S.tab === 'ucus' || S.tab === 'grafik';
    applyParams(PRESETS[b.dataset.preset], live ? { t: 0 } : {});
    if (live) play();
  }));
  $('#play').addEventListener('click', () => (S.playing ? pause() : play()));
  $('#rewind').addEventListener('click', () => { pause(); setT(0); setPlayUI(); });
  $('#reset').addEventListener('click', resetAll);
  $('#jump-ctl').addEventListener('click', () => $('#controls').scrollIntoView({ behavior: 'smooth', block: 'start' }));
  $('#pin').addEventListener('click', () => { S.ref = S.cur; S.refName = 'sabitlenen tasarım (' + nf(S.cur.score, 2) + ')'; renderToolbar(); dirty.charts = true; renderTab(true); });
  $('#scrub').addEventListener('input', (e) => { pause(); const v = +e.target.value; setT(v >= +e.target.max ? S.tEnd : v / 1000, true); setPlayUI(); });
  $('#events').addEventListener('click', (e) => { const b = e.target.closest('button[data-t]'); if (!b) return; pause(); setT(+b.dataset.t); setPlayUI(); });
  $('#sweep-var').addEventListener('change', renderSweep);
  $('#rep-step').addEventListener('change', renderReport);
  // guide links
  document.addEventListener('click', (e) => {
    const g = e.target.closest('[data-goto]');
    if (g) { e.preventDefault(); showTopic(g.dataset.goto, true); return; }
    const gs = e.target.closest('[data-gosim]');
    if (gs) { showTab(gs.dataset.gosim); window.scrollTo({ top: 0 }); return; }
    const tr = e.target.closest('.try');
    if (tr) {
      const tab = tr.dataset.sim || 'ucus';
      if (tr.dataset.try) {
        S.ref = simulate(fromCtl(DEFAULTS), true); S.refName = 'İKARUS-1 (örnek)';
        applyParams(JSON.parse(tr.dataset.try), { t: 0 });
      }
      showTab(tab); window.scrollTo({ top: 0 });
      if (tr.dataset.try && tab === 'ucus') play();
    }
  });
  let rt = null;
  window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(() => { if (S.view === 'sim' && S.tab === 'ucus') { if (S.dim === '3') size3(); else sizeCanvases(); drawFrame(); } }, 80); });
  $$('input[name=dim]').forEach((el) => el.addEventListener('change', () => setDim(el.value)));
  $$('[data-cam3]').forEach((b) => b.addEventListener('click', () => { V3.cam = b.dataset.cam3; V3.yaw = 0; V3.pitch = 0; V3.zoom = 1; store.set('ikarus-cam3', V3.cam); syncCam3(); drawFrame(); }));
  $('#pb-play').addEventListener('click', () => (S.playing ? pause() : play()));
  $$('[data-cam]').forEach((b) => b.addEventListener('click', () => { CAM.mode = b.dataset.cam; store.set('ikarus-cam', CAM.mode); syncCam(); drawFrame(); }));
  document.addEventListener('keydown', (e) => {
    if (S.view !== 'sim' || e.target.matches('input, select, textarea, button, a, summary')) return;
    if (e.code === 'Space') { e.preventDefault(); S.playing ? pause() : play(); }
  });
}

function commitNum(el, final) {
  const k = el.id.slice(2), c = NUMS[k];
  const txt = String(el.value).trim().replace(/\s/g, '').replace(',', '.');
  let v = /^[-+]?(\d+\.?\d*|\.\d+)$/.test(txt) ? parseFloat(txt) : NaN;
  if (!isFinite(v)) { if (final) showNum(k, true); return; }
  v /= c.scale;
  const cl = clampN(v, c.min, c.max);
  if (Math.abs(cl - VAL[k]) < 1e-12) { if (final && cl !== v) { const a = document.activeElement; el.blur(); showNum(k); if (a === el) el.focus(); } return; }
  setVal(k, cl, 'num');
  if (final && cl !== v) { el.value = fmtNum(k, cl); }
  recompute();
}
function syncCam() { $$('[data-cam]').forEach((b) => b.setAttribute('aria-pressed', b.dataset.cam === CAM.mode ? 'true' : 'false')); }
function init() {
  makeSprites();
  syncCam();
  syncCam3();
  buildToc();
  wire();
  writeControls(DEFAULTS);
  S.ref = simulate(fromCtl(DEFAULTS), true);
  recompute({ t: 0 });
  setPlayUI();
  setDim(store.get('ikarus-dim', '2'));
  route();
  window.addEventListener('hashchange', route);
}
function route() {
  const h = (location.hash || '').replace('#', '');
  if (h.startsWith('sim-') && ['ucus', 'puan', 'rapor', 'grafik', 'roket', 'tarama'].includes(h.slice(4))) { showTab(h.slice(4)); }
  else if (h === 'sim') showTab('ucus');
  else showTopic(topicIds.includes(h) ? h : S.topic || 'genel');
}
init();
})();
