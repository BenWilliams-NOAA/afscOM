/*
 * Selectivity functions
 */

// logistic selectivity using a50 and a95
template <class Type>
Type logis_a50_a95(int age, Type s50, Type a95) {
  Type slx = 1. / (1. + exp(-log(Type(19)) * ((age - a50) / (a95 - a50))));
  return slx;
}

// logistic selectivity function from ABL code
template <class Type>
Type logis(int a, Type a50, Type delta) {
  Type slx = 1. / (1. + exp(-log(19.0) * ((a+1) - a50) / delta));
  return slx;
}

