#include "dochi2.h"

const double unit = 1e-39;
double chi2(const TMatrixD &cov, const TH1D &data, const TH1D &mc) {
  int dimension = cov.GetNrows();

  TMatrixD diff_matrix{};
  diff_matrix.ResizeTo(dimension, 1);
  for (int i = 0; i < dimension; i++) {
    diff_matrix(i, 0) =
        (data.GetBinContent(i + 1) - mc.GetBinContent(i + 1)) / unit;
  }
  auto diff_matrix_T = diff_matrix;
  diff_matrix_T.T();

  auto cov_inv = cov * (1 / (unit * unit));
  cov_inv.Invert();

  auto chi2 = diff_matrix_T * cov_inv * diff_matrix;
  return chi2(0, 0);
}