#include "dochi2.h"
#include <TH1.h>
#include <cassert>
#include <format>
#include <iostream>
#include <ranges>

double chi2(const TMatrixDSym &cov, const TH1D &data, const TH1D &mc,
            double unit) {
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

double chi2(const TMatrixT<double> &cov, const TH1D &data, const TH1D &mc,
            double unit) {
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

std::tuple<TH1D, double> to_shape(const TH1D &input) {
  auto ret = std::make_tuple(input, input.Integral("WIDTH"));
  auto &&[shape, norm] = ret;
  shape.SetName(
      std::format("{}_shape_only", input.GetName()).c_str()); // set name
  shape.Scale(1 / norm);
  return ret;
}

std::tuple<TMatrixDSym, TH1D> ExractShapeOnly(const TMatrixDSym &cov_in,
                                              const TH1D &data_in) {
  auto ret = std::make_tuple(cov_in, data_in);
  auto &&[new_cov, new_data] = ret;
  new_data.SetName(std::format("{}_shape_only", data_in.GetName()).c_str());
  int nbins = data_in.GetNbinsX();

  if (nbins != cov_in.GetNrows()) {
    std::cerr << "Error: nbins != cov_in.GetNrows()" << std::endl;
    throw std::runtime_error("Error: bins don't match");
  }

  auto data_norm = data_in.Integral("WIDTH");

  auto axis_iter = std::ranges::iota_view(0, nbins);

  double cov_sum{}, data_sum{};
  for (auto i : axis_iter) {
    data_sum += data_in.GetBinContent(i + 1);
    for (auto j : axis_iter) {
      cov_sum += cov_in[i][j];
    }
  }

  for (auto i : axis_iter) {
    auto data_bin_i = data_in.GetBinContent(i + 1);
    for (auto j : axis_iter) {
      auto data_bin_j = data_in.GetBinContent(j + 1);
      double norm_term =
          data_bin_i * data_bin_j / data_sum / data_sum * cov_sum;
      double mix_sum_1{}, mix_sum_2{};
      for (auto k : axis_iter) {
        mix_sum_1 = cov_in[k][j];
        mix_sum_2 = cov_in[i][k];
      }
      double mix_term1 =
          data_bin_i *
          (mix_sum_1 / data_sum - cov_sum * data_bin_j / data_sum / data_sum);
      double mix_term2 =
          data_bin_j *
          (mix_sum_2 / data_sum - cov_sum * data_bin_i / data_sum / data_sum);
      new_cov[i][j] = cov_in[i][j] - norm_term - mix_term1 - mix_term2;
    }
  }

  new_data.Scale(1 / data_norm);
  new_cov *= 1. / (data_norm * data_norm);

  for (auto i : axis_iter) {
    const auto new_bin_error = std::sqrt(new_cov[i][i]);
    new_data.SetBinError(i + 1, new_bin_error);
    // std::cout << std::format("HIST: {}, BIN: {}, Content: {}, Error: {}\n",
    //                          new_data.GetName(), i,
    //                          new_data.GetBinContent(i + 1), new_bin_error);
  }

  return ret;
}

std::tuple<double, TH1D /*shape only prediction*/, TH1D /*shape only data*/,
           double /*norm*/>
do_chi2_shape_only(const TMatrixDSym &cov, const TH1D &data,
                   const TH1D &prediction) {
  auto &&[shape_only_cov, shape_only_data] = ExractShapeOnly(cov, data);
  auto [shape_only_prediction, norm] = to_shape(prediction);
  auto chi2_val =
      chi2(shape_only_cov, shape_only_data, shape_only_prediction, 1e-2);
  return std::make_tuple(chi2_val, shape_only_prediction, shape_only_data,
                         norm);
}