#!/usr/bin/env wolframscript
(* cylinder_solver.wls — cylinders through 5 points via NSolve
   Following Lichtblau (ICMS 2006): parametrize axis as {y=ax+b, z=cx+d}.
   Gives 6 complex solutions generically; count real ones (Im≈0, rsqr>0).
   Key theorem: no real cylinder iff one point is inside convex hull of other 4.
*)

(* ── Perpendicular vector from axis to point ─────────────────────────── *)
(* vec = direction (1,a,c), offset = axis point (0,b,d) *)
perpVec[pt_, vec_, offset_] := pt - offset - Projection[pt - offset, vec, Dot]

(* ── Find all 6 cylinder solutions for 5 points ─────────────────────── *)
(* Returns NSolve solution list (may include complex solutions)          *)
solveCylinders[pts5_] := Module[{vec, offset, perps, exprs},
  vec    = {1, a, c};
  offset = {0, b, d};
  perps  = Table[perpVec[pts5[[k]], vec, offset], {k, 5}];
  (* Clear denominator |vec|² from each Projection — makes it polynomial *)
  exprs  = (Numerator @ Together[#.# - rsqr]) & /@ perps;
  Quiet @ NSolve[Thread[exprs == 0], {a, b, c, d, rsqr},
                 WorkingPrecision -> 30]
]

(* ── Filter: keep real solutions with positive radius squared ─────────── *)
realSols[sols_] := Select[sols,
  Max[Abs[Im[{a, b, c, d, rsqr} /. #]]] < 1.*^-6 &&
  Re[rsqr /. #] > 0 &]

(* ── Best-fit plane normal e3 (for comparison) ───────────────────────── *)
planeNormal[pts5_] := Module[{cen, ctrd, M, eig, ord, e1, e2},
  cen  = Mean[pts5];
  ctrd = # - cen & /@ pts5;
  M    = Transpose[ctrd].ctrd / 5.;
  eig  = Eigensystem[N[M]];
  ord  = Ordering[eig[[1]], All, Greater];
  e1   = eig[[2, ord[[1]]]];
  e2   = eig[[2, ord[[2]]]];
  Normalize[Cross[e1, e2]]
]

(* ── Planarity RMS ────────────────────────────────────────────────────── *)
planarityRms[pts5_] := Module[{cen, e3, e1, e2, ctrd, M, eig, ord},
  cen  = Mean[pts5];
  ctrd = # - cen & /@ pts5;
  M    = Transpose[ctrd].ctrd / 5.;
  eig  = Eigensystem[N[M]];
  ord  = Ordering[eig[[1]], All, Greater];
  e1   = eig[[2, ord[[1]]]];
  e2   = eig[[2, ord[[2]]]];
  Sqrt[Mean @ Table[Norm[(pts5[[i]]-cen) - ((pts5[[i]]-cen).e1)*e1 -
                          ((pts5[[i]]-cen).e2)*e2]^2, {i,5}]]
]

(* ── Knot error: how far do pts5 lie from the cylinder? ──────────────── *)
cylKnotErr[pts5_, sol_] := Module[{av, bv, cv, dv, rv, vec, off, perps},
  av = Re[a /. sol]; bv = Re[b /. sol];
  cv = Re[c /. sol]; dv = Re[d /. sol];
  rv = Sqrt[Max[0., Re[rsqr /. sol]]];
  vec = N[{1, av, cv}]; off = N[{0, bv, dv}];
  perps = Table[
    Norm[pts5[[k]] - off - (((pts5[[k]]-off).vec)/(vec.vec))*vec] - rv,
    {k, 5}];
  Sqrt[Mean[perps^2]]
]

(* ── Unroll: compute φ values for 5 pts on a cylinder solution ───────── *)
unrollPhi[pts5_, sol_] := Module[
  {av, bv, cv, dv, rv, uhat, vhat, what, tmp, cen, phi},
  av = Re[a/.sol]; bv = Re[b/.sol];
  cv = Re[c/.sol]; dv = Re[d/.sol];
  rv = Sqrt[Max[0., Re[rsqr/.sol]]];
  If[rv < 1.*^-10, Return[{0.,0.,0.,0.,0.}]];
  uhat = Normalize[N[{1, av, cv}]];
  tmp  = If[Abs[uhat[[1]]] < 0.9, {1.,0.,0.}, {0.,1.,0.}];
  vhat = Normalize[tmp - (tmp.uhat)*uhat];
  what = Cross[uhat, vhat];
  cen  = N[{0., bv, dv}];
  phi  = Table[
    Module[{p = pts5[[k]] - cen, pp},
      pp = p - (p.uhat)*uhat;
      ArcTan[pp.vhat, pp.what]], {k, 5}];
  (* Unwrap *)
  Do[Module[{d = phi[[k]] - phi[[k-1]]},
    While[d >  Pi, d -= 2*Pi; phi[[k]] -= 2*Pi];
    While[d < -Pi, d += 2*Pi; phi[[k]] += 2*Pi]],
    {k, 2, 5}];
  phi
]

(* ── φ-monotonicity and arc span for one solution ─────────────────────── *)
phiStats[pts5_, sol_] := Module[{phi, diffs, mono, span},
  phi   = unrollPhi[pts5, sol];
  diffs = Differences[phi];
  mono  = (Min[diffs] > 1.*^-6) || (Max[diffs] < -1.*^-6);
  span  = Abs[phi[[5]] - phi[[1]]];   (* Δφ in radians *)
  {mono, span}
]

(* ── Analyse one 5-point window ──────────────────────────────────────── *)
(* ── Fit planar conic through 5 unrolled (u,v) pts, return coeffs ──────── *)
(* Conic: A u²+B uv+C v²+D u+E v+F = 0, normalised so that |[A..F]|=1   *)
fitConic2D[uvPts_] := Module[{M, svd, coeffs},
  M = Table[With[{u=uvPts[[k,1]], v=uvPts[[k,2]]},
    {u^2, u*v, v^2, u, v, 1}], {k, 5}];
  svd = SingularValueDecomposition[N[M]];
  coeffs = svd[[3]][[All, -1]];  (* last right singular vector = null vector *)
  coeffs / Norm[coeffs]
]

(* ── Evaluate conic at given u, return v (solving Au²+Bu v+Cv²+Du+Ev+F=0) ─ *)
(* Returns the real root closest to vHint *)
conicV[coeffs_, u_, vHint_] := Module[{A,B,C,D,E,F,disc,roots},
  {A,B,C,D,E,F} = coeffs;
  (* Cv² + (Bu+E)v + (Au²+Du+F) = 0 *)
  If[Abs[C] < 1.*^-10,
    If[Abs[B*u+E] < 1.*^-10, Return[vHint]];
    Return[-(A*u^2+D*u+F)/(B*u+E)]];
  disc = (B*u+E)^2 - 4*C*(A*u^2+D*u+F);
  If[disc < 0, Return[vHint]];
  roots = {(-(B*u+E) + Sqrt[disc])/(2*C),
           (-(B*u+E) - Sqrt[disc])/(2*C)};
  (* Pick closest to hint *)
  roots[[First @ Ordering[Abs[roots - vHint]]]]
]

(* ── Full pipeline: correct comparison ────────────────────────────────── *)
(* For each intermediate t, project the TRUE point onto the cylinder,     *)
(* get its (φ_true, z_true), evaluate the conic there, re-roll, compare.  *)
(* This correctly accounts for the nonlinear t→φ mapping.                 *)
pipelineErr[pts5_, tvs5_, sol_, trueCurve_] := Module[
  {av, bv, cv, dv, rv, uhat, vhat, what, cen, phi5, z5, uvPts,
   coeffs, errs, tMid, pTrue, pp, phiTrue, zTrue, zConic, pCyl},
  av = Re[a/.sol]; bv = Re[b/.sol];
  cv = Re[c/.sol]; dv = Re[d/.sol];
  rv = Sqrt[Max[0., Re[rsqr/.sol]]];
  uhat = Normalize[N[{1, av, cv}]];
  vhat = Normalize[# - (#.uhat)*uhat & @ If[Abs[uhat[[1]]]<0.9,{1.,0.,0.},{0.,1.,0.}]];
  what = Cross[uhat, vhat];
  cen  = N[{0., bv, dv}];
  (* Unroll 5 knots *)
  phi5  = unrollPhi[pts5, sol];
  z5    = Table[(pts5[[k]]-cen).uhat, {k,5}];
  uvPts = Transpose[{rv*phi5, z5}];
  coeffs = fitConic2D[uvPts];
  (* For each of 4 midpoints + 8 sub-midpoints (12 total), project true
     curve point onto cylinder, then evaluate the fitted conic there *)
  errs = Flatten @ Table[
    Table[
      tMid   = tvs5[[k]] + s*(tvs5[[k+1]]-tvs5[[k]]);
      pTrue  = N[trueCurve[tMid]];
      (* Project pTrue perp to axis → get φ and z on cylinder surface *)
      pp     = pTrue - cen;
      zTrue  = pp.uhat;                      (* axial coord *)
      With[{pp2 = pp - zTrue*uhat},          (* radial part *)
        phiTrue = ArcTan[pp2.vhat, pp2.what]
      ];
      (* Unwrap phiTrue to be near the knot interval *)
      phiTrue += 2*Pi*Round[(phi5[[k]]-phiTrue)/(2*Pi)];
      zConic = conicV[coeffs, rv*phiTrue, (z5[[k]]+z5[[k+1]])/2.];
      (* Re-roll the conic point back to 3D *)
      pCyl   = cen + zConic*uhat + rv*Cos[phiTrue]*vhat + rv*Sin[phiTrue]*what;
      Norm[pCyl - pTrue],
      {s, {0.25, 0.5, 0.75}}],
    {k, 1, 4}];
  Max[errs]
]

(* ── Compare correct vs wrong cylinder ────────────────────────────────── *)
compareWindow[pts5_, tvs5_, label_String, trueCurve_] := Module[
  {planRms, e3, sols, real, nReal, stats, monoIdx, nMono, spans,
   sortedIdx, bestSol, worstSol, errBest, errWrong, fn},

  planRms = planarityRms[pts5];
  e3      = planeNormal[pts5];
  sols    = solveCylinders[pts5];
  real    = realSols[sols];
  nReal   = Length[real];

  fn[x_] := If[x == 0., "0", Module[{me = MantissaExponent[N[x]]},
    ToString[Round[me[[1]]*1000]/1000.] <> "e" <> ToString[me[[2]]-1]]];

  If[nReal < 2,
    Print[StringPadRight[label,8] <> "nReal=" <> ToString[nReal] <>
      "  [can't compare — need ≥2]"];
    Return[{}]];

  stats     = phiStats[pts5, #] & /@ real;
  spans     = stats[[All, 2]];
  sortedIdx = Ordering[spans];            (* ascending span *)
  bestSol   = real[[ sortedIdx[[1]] ]];  (* smallest span = "correct" *)
  worstSol  = real[[ sortedIdx[[-1]] ]]; (* largest span  = "wrong"   *)

  errBest  = pipelineErr[pts5, tvs5, bestSol,  trueCurve];
  errWrong = pipelineErr[pts5, tvs5, worstSol, trueCurve];

  Print[StringPadRight[label,8] <>
    "spans=" <> ToString[Round[N[spans[[sortedIdx]]]*180./Pi, 0.1]] <> "deg" <>
    "  err_correct=" <> fn[errBest] <>
    "  err_wrong=" <> fn[errWrong] <>
    "  ratio=" <> fn[errWrong/(errBest+1.*^-20)]];
]

analyseWindow[pts5_, label_String] := Module[
  {planRms, e3, sols, real, nReal, stats, monoIdx, nMono,
   spans, bestIdx, bestSol, knot, tilt},

  planRms = planarityRms[pts5];
  e3      = planeNormal[pts5];
  sols    = solveCylinders[pts5];
  real    = realSols[sols];
  nReal   = Length[real];

  If[nReal > 0,
    (* Compute {mono, span} for each real solution *)
    stats   = phiStats[pts5, #] & /@ real;
    (* Monotone solutions *)
    monoIdx = Select[Range[nReal], stats[[#,1]] &];
    nMono   = Length[monoIdx];
    (* Among monotone, pick smallest Δφ span — the "tightest local arc" *)
    bestIdx = If[nMono > 0,
      monoIdx[[First @ Ordering[stats[[monoIdx, 2]]]]],
      First @ Ordering[stats[[All, 2]]]   (* fallback: smallest span anyway *)
    ];
    bestSol = real[[bestIdx]];
    spans   = stats[[All, 2]];   (* Δφ for all real solutions *)
    knot    = cylKnotErr[pts5, bestSol];
    tilt    = ArcCos[Clip[Abs[
      Normalize[N[{1, Re[a/.bestSol], Re[c/.bestSol]}]] . e3
    ], {0.,1.}]] * 180./Pi,
    (* No real cylinder *)
    nMono = 0; spans = {}; knot = planRms; tilt = 0.
  ];

  fn[x_] := If[x == 0., "0", Module[{me = MantissaExponent[N[x]]},
    ToString[Round[me[[1]]*1000]/1000.] <> "e" <> ToString[me[[2]]-1]]];
  fnd[x_] := ToString[Round[N[x] * 180./Pi, 0.1]];   (* radians → degrees *)

  Print[StringPadRight[label, 8] <>
    "nReal=" <> ToString[nReal] <>
    "  nMono=" <> ToString[nMono] <>
    If[nReal > 0,
      "  spans(deg)=" <> ToString[Round[N[spans]*180./Pi, 0.1]] <>
      "  knot=" <> fn[knot] <>
      "  plan/knot=" <> fn[planRms/knot],
      "  [NO REAL CYL — use planar fallback]"]
  ];

  {nReal, nMono}
]

(* ═══════════════════════════════════════════════════════════════════════ *)

Print["=== Test A: Helix  r(t)=(cos 2t, sin 2t, 0.5t),  n=12 ==="];
Print["    exact axis=(0,0,1); expect nReal=6, knot→machine-eps"];
Module[{n=12, tv, hPts, pts5},
  tv   = Table[2.*Pi*i/(n-1), {i, 0, n-1}];
  hPts = {Cos[2.*#], Sin[2.*#], 0.5*#} & /@ tv;
  nrs = Table[
    pts5 = hPts[[i;;i+4]];
    analyseWindow[pts5, "win "<>ToString[i-1]],
    {i, 1, n-4}];
  Print["  {nReal, nMono} per window: " <> ToString[nrs]];
];

Print[""];
Print["=== Test B: Twisted cubic  r(t)=(t,t²,t³),  n=12 ==="];
Print["    Correct=smallest span, Wrong=largest span; err at 4 midpoints"];
Module[{n=12, tv, tcPts, pts5, tvs5, tc},
  tc    = Function[t, {t, t^2, t^3}];
  tv    = Table[N[(i-1)/(n-1)], {i, 1, n}];
  tcPts = tc /@ tv;
  Do[
    pts5 = tcPts[[i;;i+4]];
    tvs5 = tv[[i;;i+4]];
    compareWindow[pts5, tvs5, "win "<>ToString[i-1], tc],
    {i, 1, n-4}]
];
