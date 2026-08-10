(* Verify the stable analytic inversions used by the Uryu F-space solvers. *)
ClearAll[fe, fj, ee, delta, radialDenominator, angularFactor, lambda1,
  lambda2, fmax, ratio, q];

equatorialDelta = 2 fe ee/(1 + Sqrt[1 + 4 fe^2 ee]);
equatorialCheck = FullSimplify[
  equatorialDelta/(ee - equatorialDelta^2) - fe,
  Assumptions -> fe > 0 && ee > 0 &&
    Element[{fe, ee}, Reals]
  ];

localDelta = 2 fj radialDenominator/
  (angularFactor + Sqrt[angularFactor^2 +
    4 fj^2 angularFactor radialDenominator]);
localCheck = FullSimplify[
  localDelta angularFactor/
    (radialDenominator - localDelta^2 angularFactor) - fj,
  Assumptions -> fj > 0 && radialDenominator > 0 &&
    angularFactor > 0 &&
    Element[{fj, radialDenominator, angularFactor}, Reals]
  ];

uryuNumerator = lambda2 (ratio fmax)^q - lambda1 fmax^q;
boundaryCheck = FullSimplify[
  uryuNumerator /. ratio -> (lambda1/lambda2)^(1/q),
  Assumptions -> lambda1 > 0 && lambda2 > 0 && fmax > 0 &&
    Element[q, PositiveIntegers]
  ];
monotonicCheck = FullSimplify[
  D[uryuNumerator, ratio] > 0,
  Assumptions -> lambda1 > 0 && lambda2 > 0 && fmax > 0 &&
    ratio > 0 && Element[q, PositiveIntegers]
  ];

checks = {
  equatorialCheck === 0,
  localCheck === 0,
  boundaryCheck === 0,
  monotonicCheck === True
  };

If[And @@ checks,
  Print["URYU_MOMENTUM_INVERSION=PASS"],
  Print["URYU_MOMENTUM_INVERSION=FAIL"];
  Print[{equatorialCheck, localCheck, boundaryCheck, monotonicCheck}];
  Exit[1]
  ];
