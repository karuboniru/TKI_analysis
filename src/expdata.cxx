#include <TFile.h>
#include <TH2.h>
#include <expdata.h>

#include <TH1.h>
#include <TMatrixD.h>
#include <array>
#include <dochi2.h>
#include <iostream>

// pi0 channel: https://arxiv.org/src/2002.05812v4/anc/SupplementalMaterial2.txt
// 0pi channel: https://arxiv.org/src/1805.05486v4/anc/SupplementalMaterial2.txt

namespace {
const std::string pizero_path{PIZERO_DATA};
const std::string zeropi_path{ZEROPI_DATA};
} // namespace

namespace MINERvA_TKI::pi0 {
TFile pizero_file(pizero_path.c_str(), "READONLY");
}

namespace MINERvA_TKI::ZeroPi {
TFile zeropi_file(zeropi_path.c_str(), "READONLY");
}

namespace MINERvA_TKI::pi0::IApN {
std::array<double, dimension + 1> get_binning() {
  std::array<double, dimension + 1> bin_edges{
      0.0000,   55.0000,  110.0000, 165.0000, 220.0000, 275.0000, 330.0000,
      385.0000, 440.0000, 495.0000, 560.0000, 655.0000, 810.0000};
  for (auto &&i : bin_edges) {
    i /= 1000.;
  }
  return bin_edges;
}
TH1D get_hist_impl() {
  auto bin_edges = get_binning();
  // 10$^{-42}$ cm$^2$/MeV/c/nucleon -> 10$^{-39}$ cm$^2$/GeV/c/nucleon
  // to 10$^{-38}$ cm$^2$/GeV/c/nucleon
  // factor: 10^4 / 10^3 = 10
  // std::array<double, 12> exp_value{0.6, 1.3, 2.1, 2.8, 1.4, 1.2,
  //                                  1.1, 1.2, 1.3, 1.2, 0.8, 0.7};
  auto hist_from_root =
      (TH1D *)((pizero_file.Get<TList>("neutronmomentum"))->At(0));

  TH1D data("data_IApN_pi0", "MINERvA#kern[0.25]{#pi^{0}} Data", dimension,
            bin_edges.data());

  auto cov = get_cov();
  for (size_t i = 0; i < dimension; i++) {
    data.SetBinContent(i + 1,
                       hist_from_root->GetBinContent(i + 1) / (1e-38 / 1e3));
    // data.SetBinError(i + 1, std::sqrt(cov[i][i]));
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }

  return data;
}
const TH1D &get_hist() {
  static TH1D hist = get_hist_impl();
  return hist;
}
TMatrixDSym get_cov_impl() {
  auto cov_from_root =
      (TMatrixD *)((pizero_file.Get<TList>("neutronmomentum"))->At(2));
  TMatrixDSym cov_matrix(dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      cov_matrix[i][j] = (*cov_from_root)[i + 1][j + 1] / 1e-38 / 1e-38 * 1e6;
    }
  }
  return cov_matrix;
}
const TMatrixDSym &get_cov() {
  static TMatrixDSym cov = get_cov_impl();
  return cov;
}

double do_chi2(TH1 *hist) {
  auto &data = get_hist();
  return chi2(get_cov(), data, *((TH1D *)hist), 1e-6);
}

} // namespace MINERvA_TKI::pi0::IApN

namespace MINERvA_TKI::pi0::dalphat {
std::array<double, dimension + 1> get_binning() {
  std::array<double, dimension + 1> bin_edges{
      0.0000,   20.0000,  40.0000,  60.0000,  80.0000,
      100.0000, 120.0000, 140.0000, 160.0000, 180.0000};
  return bin_edges;
}

TH1D get_hist_impl() {
  auto bin_edges = get_binning();
  auto hist_from_root = (TH1D *)((pizero_file.Get<TList>("dalphat"))->At(0));
  TH1D data("data_dalphat_pi0", "MINERvA#kern[0.25]{#pi^{0}} Data", dimension,
            bin_edges.data());

  auto cov = get_cov();
  for (size_t i = 0; i < dimension; i++) {
    // data.SetBinContent(i + 1, exp_value[i] / 1e4);
    data.SetBinContent(i + 1, hist_from_root->GetBinContent(i + 1) / (1e-38));
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }

  return data;
}
const TH1D &get_hist() {
  static TH1D hist = get_hist_impl();
  return hist;
}
TMatrixDSym get_cov_impl() {
  auto cov_from_root = (TMatrixD *)((pizero_file.Get<TList>("dalphat"))->At(2));
  TMatrixDSym cov_matrix(dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      cov_matrix[i][j] = (*cov_from_root)[i + 1][j + 1] / 1e-38 / 1e-38;
    }
  }
  // cov_matrix *= 1e-18;
  return cov_matrix;
}
const TMatrixDSym &get_cov() {
  static TMatrixDSym cov = get_cov_impl();
  return cov;
}
double do_chi2(TH1 *hist) {
  auto &data = get_hist();
  return chi2(get_cov(), data, *((TH1D *)hist), 1e-12);
}
} // namespace MINERvA_TKI::pi0::dalphat

namespace MINERvA_TKI::ZeroPi::IApN {
std::array<double, dimension + 1> get_binning() {
  std::array<double, dimension + 1> bin_edges{
      0.0000, 0.0250, 0.0500, 0.0750, 0.1000, 0.1250, 0.1500, 0.1750, 0.2000,
      0.2250, 0.2500, 0.2750, 0.3000, 0.3500, 0.4000, 0.4500, 0.5000, 0.5500,
      0.6000, 0.6500, 0.7000, 0.8000, 1.0000, 1.2000, 2.0000};
  return bin_edges;
}
TH1D get_hist_impl() {
  auto bin_edges = get_binning();
  auto hist_from_root =
      (TH1D *)((zeropi_file.Get<TList>("neutronmomentum"))->At(0));

  TH1D data("data_IApN_0pi", "MINERvA 0#pi Data", dimension,
            bin_edges.data());

  auto cov = get_cov();
  for (size_t i = 0; i < dimension; i++) {
    // data.SetBinContent(i + 1, exp_value[i] / 1e3);
    data.SetBinContent(i + 1, hist_from_root->GetBinContent(i + 2) / (1e-38));
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }

  return data;
}
const TH1D &get_hist() {
  static TH1D hist = get_hist_impl();
  return hist;
}
TMatrixDSym get_cov_impl() {
  auto cov_from_root =
      (TMatrixD *)((zeropi_file.Get<TList>("neutronmomentum"))->At(2));
  TMatrixDSym cov_matrix(dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      cov_matrix[i][j] = (*cov_from_root)[i + 2][j + 2] / 1e-38 / 1e-38;
    }
  }
  return cov_matrix;
}
const TMatrixDSym &get_cov() {
  static TMatrixDSym cov = get_cov_impl();
  return cov;
}
double do_chi2(TH1 *hist) {
  auto &data = get_hist();

  return chi2(get_cov(), data, *((TH1D *)hist), 1e-6);
}
} // namespace MINERvA_TKI::ZeroPi::IApN

namespace MINERvA_TKI::ZeroPi::dalphat {
std::array<double, dimension + 1> get_binning() {
  std::array<double, dimension + 1> bin_edges{
      0.0000,   20.0000,  40.0000,  60.0000,  80.0000,  100.0000, 120.0000,
      130.0000, 140.0000, 150.0000, 160.0000, 170.0000, 180.0000};
  return bin_edges;
}

TH1D get_hist_impl() {
  auto bin_edges = get_binning();
  auto hist_from_root = (TH1D *)((zeropi_file.Get<TList>("dalphat"))->At(0));
  TH1D data("data_dalphat_0pi", "MINERvA 0#pi Data", dimension,
            bin_edges.data());

  auto cov = MINERvA_TKI::ZeroPi::dalphat::get_cov();
  for (size_t i = 0; i < dimension; i++) {
    // data.SetBinContent(i + 1, exp_value[i] / 1e4);
    data.SetBinContent(i + 1, hist_from_root->GetBinContent(i + 2) / (1e-38));
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }

  return data;
}
const TH1D &get_hist() {
  static TH1D hist = get_hist_impl();
  return hist;
}
TMatrixDSym get_cov_impl() {
  auto cov_from_root = (TMatrixD *)((zeropi_file.Get<TList>("dalphat"))->At(2));
  TMatrixDSym cov_matrix(dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      cov_matrix[i][j] = (*cov_from_root)[i + 2][j + 2] / 1e-38 / 1e-38;
    }
  }
  return cov_matrix;
}
const TMatrixDSym &get_cov() {
  static TMatrixDSym cov = get_cov_impl();
  return cov;
}
double do_chi2(TH1 *hist) {
  auto &data = get_hist();

  return chi2(get_cov(), data, *((TH1D *)hist), 1e-6);
}

} // namespace MINERvA_TKI::ZeroPi::dalphat

constexpr double rag_to_deg = 180. / M_PI;

namespace T2K_STK {
TFile T2K_STK_dat_file(T2K_STK_DATA, "READONLY");
std::array<double, dimension + 1> get_binning() {
  std::array<double, dimension + 1> bin_edges{0,    .47,  1.02, 1.54, 1.98,
                                              2.34, 2.64, 2.89, M_PI};
  for (auto &&i : bin_edges) {
    i *= rag_to_deg;
  }
  return bin_edges;
}

TMatrixDSym get_cov_impl() {
  auto cov_from_root = T2K_STK_dat_file.Get<TH2D>("Covariance_Matrix");
  TMatrixDSym cov_matrix(dimension);
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      cov_matrix[i][j] =
          cov_from_root->GetBinContent(i + 1, j + 1) / rag_to_deg / rag_to_deg;
    }
  }
  return cov_matrix;
}
const TMatrixDSym &get_cov() {
  static TMatrixDSym cov = get_cov_impl();
  return cov;
}

TH1D get_hist_impl() {
  auto data = T2K_STK_dat_file.Get<TH1D>("Result");
  auto bin_edges = get_binning();
  TH1D data_hist("data_IApN_T2K", "T2K 0#pi Data", dimension,
                 bin_edges.data());
  for (size_t i = 0; i < dimension; i++) {
    data_hist.SetBinContent(i + 1,
                            data->GetBinContent(i + 1) / 1e-38 / rag_to_deg);
    data_hist.SetBinError(i + 1, data->GetBinError(i + 1) / 1e-38 / rag_to_deg);
  }
  return data_hist;
}
const TH1D &get_hist() {
  static TH1D hist = get_hist_impl();
  return hist;
}

double do_chi2(TH1 *hist) {
  auto &data = get_hist();
  return chi2(get_cov(), data, *((TH1D *)hist), 1e-6);
}
} // namespace T2K_STK

namespace MicroBooNE {
namespace pi0_momentum {
std::array<double, 9> get_binning() {
  return {0, .1, .15, .2, .3, .4, .5, .6, .799};
}

TH1D get_hist_impl() {
  auto bin_edges = get_binning();
  TH1D data("data_pion_momentum", "MicroBooNE CC#pi^{0} Data", dimension,
            bin_edges.data());
  std::array<double, 8> exp_value{0.66, 1.92, 2.44, 1.15,
                                  0.73, 0.34, 0.18, 0.09};
  auto &cov = get_cov();
  for (size_t i = 0; i < dimension; i++) {
    data.SetBinContent(i + 1, exp_value[i] * 1e-1);
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }
  return data;
}
const TH1D &get_hist() {
  static TH1D hist = get_hist_impl();
  return hist;
}

TMatrixDSym get_cov_impl() {
  constexpr std::array<std::array<double, dimension>, dimension> matrix{
      {{0.032, 0.015, 0.005, 0.006, 0.008, 0.007, 0.006, 0.004},
       {0.015, 0.023, 0.021, 0.011, 0.008, 0.009, 0.007, 0.004},
       {0.005, 0.021, 0.031, 0.024, 0.013, 0.009, 0.008, 0.005},
       {0.006, 0.011, 0.024, 0.035, 0.023, 0.009, 0.006, 0.006},
       {0.008, 0.008, 0.013, 0.023, 0.027, 0.014, 0.003, -0.000},
       {0.007, 0.009, 0.009, 0.009, 0.014, 0.016, 0.007, -0.002},
       {0.006, 0.007, 0.008, 0.006, 0.003, 0.007, 0.012, 0.009},
       {0.004, 0.004, 0.005, 0.006, -0.000, -0.002, 0.009, 0.017}}};
  TMatrixDSym cov_matrix(dimension);
  auto &&pion_mom_bin_edges = get_binning();
  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      cov_matrix[i][j] = matrix[i][j] / 1e4 /
                         (pion_mom_bin_edges[j + 1] - pion_mom_bin_edges[j]) /
                         (pion_mom_bin_edges[i + 1] - pion_mom_bin_edges[i]);
    }
  }
  return cov_matrix;
}
const TMatrixDSym &get_cov() {
  static TMatrixDSym cov = get_cov_impl();
  return cov;
}

TMatrixDSym get_smear_impl() {
  constexpr std::array<std::array<double, dimension>, dimension> matrix{
      {{0.750, 0.320, 0.010, -0.027, -0.002, 0.014, 0.026, 0.024},
       {0.273, 0.393, 0.273, 0.025, -0.019, 0.012, 0.016, 0.012},
       {0.009, 0.282, 0.427, 0.229, 0.024, -0.014, 0.003, -0.004},
       {-0.061, 0.001, 0.266, 0.525, 0.284, -0.019, -0.084, -0.023},
       {-0.030, -0.040, 0.026, 0.259, 0.429, 0.239, -0.008, -0.069},
       {0.014, 0.011, -0.019, 0.006, 0.270, 0.477, 0.287, 0.016},
       {0.014, -0.001, -0.023, -0.070, 0.008, 0.273, 0.455, 0.313},
       {-0.004, -0.026, -0.057, -0.076, -0.132, -0.058, 0.270, 0.429}}};
  TMatrixDSym smear_matrix(dimension);
  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      smear_matrix[i][j] = matrix[i][j];
    }
  }
  return smear_matrix;
}
const TMatrixDSym &get_smear() {
  static TMatrixDSym smear = get_smear_impl();
  return smear;
}

double do_chi2(TH1 *hist_smeared) {
  auto &data = get_hist();
  return chi2(get_cov(), data, *((TH1D *)hist_smeared), 1e-6);
}

} // namespace pi0_momentum

namespace pi0_angular {
std::array<double, dimension + 1> get_binning() {
  return {-1.0, -0.5, -.25, 0, .25, .5, .75, 1.};
}

TH1D get_hist_impl() {
  auto bin_edges = get_binning();
  TH1D data("data_pion_angular", "MicroBooNE CC#pi^{0} Data", dimension,
            bin_edges.data());
  std::array<double, 7> exp_value{.1, .17, .20, .19, .23, .40, .49};
  auto &cov = get_cov();
  for (size_t i = 0; i < dimension; i++) {
    data.SetBinContent(i + 1, exp_value[i] * 1e-1);
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }
  return data;
}
const TH1D &get_hist() {
  static TH1D hist = get_hist_impl();
  return hist;
}

TMatrixDSym get_cov_impl() {
  constexpr std::array<std::array<double, dimension>, dimension> matrix{
      {{0.031, -0.004, 0.007, 0.013, 0.007, 0.019, 0.024},
       {-0.004, 0.026, -0.006, -0.002, 0.013, 0.004, 0.005},
       {0.007, -0.006, 0.028, 0.004, -0.010, 0.013, 0.013},
       {0.013, -0.002, 0.004, 0.023, 0.004, 0.006, 0.025},
       {0.007, 0.013, -0.010, 0.004, 0.035, 0.015, 0.003},
       {0.019, 0.004, 0.013, 0.006, 0.015, 0.050, 0.030},
       {0.024, 0.005, 0.013, 0.025, 0.003, 0.030, 0.085}}};
  TMatrixDSym cov_matrix(dimension);
  auto &&pion_mom_bin_edges = get_binning();
  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      cov_matrix[i][j] = matrix[i][j] / 1e4 /
                         (pion_mom_bin_edges[j + 1] - pion_mom_bin_edges[j]) /
                         (pion_mom_bin_edges[i + 1] - pion_mom_bin_edges[i]);
    }
  }
  return cov_matrix;
}
const TMatrixDSym &get_cov() {
  static TMatrixDSym cov = get_cov_impl();
  return cov;
}

TMatrixDSym get_smear_impl() {
  constexpr std::array<std::array<double, dimension>, dimension> matrix{
      {{0.872, 0.064, 0.000, 0.069, 0.036, 0.072, 0.061},
       {0.047, 0.792, 0.133, -0.063, 0.081, -0.037, -0.061},
       {-0.014, 0.118, 0.714, 0.168, -0.142, 0.010, -0.008},
       {0.047, -0.054, 0.185, 0.660, 0.166, -0.090, 0.075},
       {-0.007, 0.061, -0.132, 0.135, 0.640, 0.147, -0.125},
       {0.027, -0.020, 0.026, -0.096, 0.203, 0.721, 0.147},
       {-0.004, -0.053, -0.017, 0.044, -0.135, 0.084, 0.729}}};
  TMatrixDSym smear_matrix(dimension);
  for (int i = 0; i < dimension; ++i) {
    for (int j = 0; j < dimension; ++j) {
      smear_matrix[i][j] = matrix[i][j];
    }
  }
  return smear_matrix;
}
const TMatrixDSym &get_smear() {
  static TMatrixDSym smear = get_smear_impl();
  return smear;
}

double do_chi2(TH1 *hist_smeared) {
  auto &data = get_hist();
  return chi2(get_cov(), data, *((TH1D *)hist_smeared), 1e-6);
}

} // namespace pi0_angular
} // namespace MicroBooNE
