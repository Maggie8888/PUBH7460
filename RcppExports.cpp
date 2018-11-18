
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
// rcpp_distance
NumericMatrix rcpp_distance(NumericMatrix l, int m, int N);
RcppExport SEXP FastGP_rcpp_distance(SEXP lSEXP, SEXP mSEXP, SEXP NSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericMatrix >::type l(lSEXP);
  Rcpp::traits::input_parameter< int >::type m(mSEXP);
  Rcpp::traits::input_parameter< int >::type N(NSEXP);
  __result = Rcpp::wrap(rcpp_distance(l, m, N));
  return __result;
  END_RCPP
}
// durbin
NumericVector durbin(NumericVector r, int N);
RcppExport SEXP FastGP_durbin(SEXP rSEXP, SEXP NSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
  Rcpp::traits::input_parameter< int >::type N(NSEXP);
  __result = Rcpp::wrap(durbin(r, N));
  return __result;
  END_RCPP
}
// trench
NumericVector trench(NumericVector r, NumericVector y, int N);
RcppExport SEXP FastGP_trench(SEXP rSEXP, SEXP ySEXP, SEXP NSEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
  Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
  Rcpp::traits::input_parameter< int >::type N(NSEXP);
  __result = Rcpp::wrap(trench(r, y, N));
  return __result;
  END_RCPP
}
// rcppeigen_invert_matrix
Eigen::MatrixXd rcppeigen_invert_matrix(const Eigen::Map<Eigen::MatrixXd>& A);
RcppExport SEXP FastGP_rcppeigen_invert_matrix(SEXP ASEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type A(ASEXP);
  __result = Rcpp::wrap(rcppeigen_invert_matrix(A));
  return __result;
  END_RCPP
}
// rcppeigen_get_diag
Eigen::MatrixXd rcppeigen_get_diag(const Eigen::Map<Eigen::MatrixXd>& A);
RcppExport SEXP FastGP_rcppeigen_get_diag(SEXP ASEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type A(ASEXP);
  __result = Rcpp::wrap(rcppeigen_get_diag(A));
  return __result;
  END_RCPP
}
// rcppeigen_get_det
double rcppeigen_get_det(const Eigen::Map<Eigen::MatrixXd>& A);
RcppExport SEXP FastGP_rcppeigen_get_det(SEXP ASEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type A(ASEXP);
  __result = Rcpp::wrap(rcppeigen_get_det(A));
  return __result;
  END_RCPP
}
// rcppeigen_get_chol
Eigen::MatrixXd rcppeigen_get_chol(const Eigen::Map<Eigen::MatrixXd>& A);
RcppExport SEXP FastGP_rcppeigen_get_chol(SEXP ASEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type A(ASEXP);
  __result = Rcpp::wrap(rcppeigen_get_chol(A));
  return __result;
  END_RCPP
}
// rcppeigen_get_chol_stable
Eigen::MatrixXd rcppeigen_get_chol_stable(const Eigen::Map<Eigen::MatrixXd>& A);
RcppExport SEXP FastGP_rcppeigen_get_chol_stable(SEXP ASEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type A(ASEXP);
  __result = Rcpp::wrap(rcppeigen_get_chol_stable(A));
  return __result;
  END_RCPP
}
// rcppeigen_get_chol_diag
Eigen::VectorXd rcppeigen_get_chol_diag(const Eigen::Map<Eigen::MatrixXd>& A);
RcppExport SEXP FastGP_rcppeigen_get_chol_diag(SEXP ASEXP) {
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type A(ASEXP);
  __result = Rcpp::wrap(rcppeigen_get_chol_diag(A));
  return __result;
  END_RCPP
}
