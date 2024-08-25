#include <TFile.h>
#include <expdata.h>

#include <TH1.h>
#include <TMatrixD.h>
#include <array>
#include <dochi2.h>
#include <iostream>

// pi0 channel: https://arxiv.org/src/2002.05812v4/anc/SupplementalMaterial2.txt
// 0pi channel: https://arxiv.org/src/1805.05486v4/anc/SupplementalMaterial2.txt

const std::string pizero_path{
    "/var/home/yan/neutrino/tkidata/PIZEROTKI_MINERvA.root"};
const std::string zeropi_path{
    "/var/home/yan/neutrino/tkidata/PhysRevD.101.092001.root"};
TFile pizero_file(pizero_path.c_str(), "READONLY");
TFile zeropi_file(zeropi_path.c_str(), "READONLY");

std::array<double, 13> get_binning_IApN_pi0() {
  std::array<double, 13> bin_edges{
      0.0000,   55.0000,  110.0000, 165.0000, 220.0000, 275.0000, 330.0000,
      385.0000, 440.0000, 495.0000, 560.0000, 655.0000, 810.0000};
  for (auto &&i : bin_edges) {
    i /= 1000.;
  }
  return bin_edges;
}

std::array<double, 25> get_binning_IApN_0pi() {
  std::array<double, 25> bin_edges{
      0.0000, 0.0250, 0.0500, 0.0750, 0.1000, 0.1250, 0.1500, 0.1750, 0.2000,
      0.2250, 0.2500, 0.2750, 0.3000, 0.3500, 0.4000, 0.4500, 0.5000, 0.5500,
      0.6000, 0.6500, 0.7000, 0.8000, 1.0000, 1.2000, 2.0000};
  return bin_edges;
}

std::array<double, 10> get_binning_dalphat_pi0() {
  std::array<double, 10> bin_edges{0.0000,   20.0000,  40.0000,  60.0000,
                                   80.0000,  100.0000, 120.0000, 140.0000,
                                   160.0000, 180.0000};
  return bin_edges;
}

std::array<double, 13> get_binning_dalphat_0pi() {
  std::array<double, 13> bin_edges{
      0.0000,   20.0000,  40.0000,  60.0000,  80.0000,  100.0000, 120.0000,
      130.0000, 140.0000, 150.0000, 160.0000, 170.0000, 180.0000};
  return bin_edges;
}

TH1D get_IApN_hist_pi0() {
  auto bin_edges = get_binning_IApN_pi0();
  // 10$^{-42}$ cm$^2$/MeV/c/nucleon -> 10$^{-39}$ cm$^2$/GeV/c/nucleon
  // to 10$^{-38}$ cm$^2$/GeV/c/nucleon
  // factor: 10^4 / 10^3 = 10
  // std::array<double, 12> exp_value{0.6, 1.3, 2.1, 2.8, 1.4, 1.2,
  //                                  1.1, 1.2, 1.3, 1.2, 0.8, 0.7};
  auto hist_from_root =
      (TH1D *)(((TList *)pizero_file.Get("neutronmomentum"))->At(0));

  TH1D data("data_IApN_pi0", "MINERvA #pi^{0} Data", bin_edges.size() - 1,
            bin_edges.data());

  auto cov = get_cov_IApN_pi0();
  for (size_t i = 0; i < bin_edges.size() - 1; i++) {
    data.SetBinContent(i + 1,
                       hist_from_root->GetBinContent(i + 1) / (1e-38 / 1e3));
    // data.SetBinError(i + 1, std::sqrt(cov[i][i]));
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }

  return data;
}

TH1D get_dalphat_hist_pi0() {
  auto bin_edges = get_binning_dalphat_pi0();
  // std::array<double, 9>
  // exp_value{4.5, 4.5, 4.3, 4.1, 4.3, 4.7, 5.9, 8.2, 10.2};
  auto hist_from_root = (TH1D *)(((TList *)pizero_file.Get("dalphat"))->At(0));
  TH1D data("data_dalphat_pi0", "MINERvA #pi^{0} Data", bin_edges.size() - 1,
            bin_edges.data());

  auto cov = get_cov_dalphat_pi0();
  for (size_t i = 0; i < bin_edges.size() - 1; i++) {
    // data.SetBinContent(i + 1, exp_value[i] / 1e4);
    data.SetBinContent(i + 1, hist_from_root->GetBinContent(i + 1) / (1e-38));
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }

  return data;
}

TH1D get_IApN_hist_0pi() {
  auto bin_edges = get_binning_IApN_0pi();
  // 10$^{-41}$ cm$^2$/GeV/c/nucleon
  // to 10$^{-38}$ cm$^2$/GeV/c/nucleon
  // factor: 10^3
  // std::array<double, 24> exp_value{35.0,  106.3, 185.2, 294.7, 371.7, 413.7,
  //                                  479.4, 500.9, 381.9, 302.2, 265.1, 258.5,
  //                                  267.9, 264.3, 195.2, 165.3, 158.5, 163.7,
  //                                  124.6, 114.4, 88.4,  46.0,  16.2,  1.6};

  auto hist_from_root =
      (TH1D *)(((TList *)zeropi_file.Get("neutronmomentum"))->At(0));

  TH1D data("data_IApN_0pi", "MINERvA 0#pi Data", bin_edges.size() - 1,
            bin_edges.data());

  auto cov = get_cov_IApN_0pi();
  for (size_t i = 0; i < bin_edges.size() - 1; i++) {
    // data.SetBinContent(i + 1, exp_value[i] / 1e3);
    data.SetBinContent(i + 1, hist_from_root->GetBinContent(i + 2) / (1e-38));
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }

  return data;
}

TH1D get_dalphat_hist_0pi() {
  auto bin_edges = get_binning_dalphat_0pi();
  // 10$^{-42}$ cm$^2$/deg
  // to 10$^{-38}$ cm$^2$/deg
  // factor: 10^3
  // std::array<double, 12> exp_value{7.8,  6.5,  7.3,  7.0,  8.0,  9.4,
  //                                  11.6, 14.1, 16.2, 19.4, 19.0, 16.9};
  auto hist_from_root = (TH1D *)(((TList *)zeropi_file.Get("dalphat"))->At(0));
  TH1D data("data_dalphat_0pi", "MINERvA 0#pi Data", bin_edges.size() - 1,
            bin_edges.data());

  auto cov = get_cov_dalphat_0pi();
  for (size_t i = 0; i < bin_edges.size() - 1; i++) {
    // data.SetBinContent(i + 1, exp_value[i] / 1e4);
    data.SetBinContent(i + 1, hist_from_root->GetBinContent(i + 2) / (1e-38));
    data.SetBinError(i + 1, std::sqrt(cov[i][i]));
  }

  return data;
}

double do_chi2_dalphat_pi0(TH1 *hist) {
  auto data = get_dalphat_hist_pi0();

  // for (size_t i{}; i < 12; i++) {
  //   std::cout << i << " th: " << cov_matrix[i][i] << " ";
  // }
  std::cout << std::endl;
  return chi2(get_cov_dalphat_pi0(), data, *((TH1D *)hist), 1e-12);
}

double do_chi2_IApN_pi0(TH1 *hist) {
  auto data = get_IApN_hist_pi0();

  // for (size_t i{}; i < 12; i++) {
  //   std::cout << i << " th: " << cov_matrix[i][i] << " ";
  // }
  std::cout << std::endl;
  return chi2(get_cov_IApN_pi0(), data, *((TH1D *)hist), 1e-6);
}

double do_chi2_IApN_0pi(TH1 *hist) {
  auto data = get_IApN_hist_0pi();

  std::cout << std::endl;
  return chi2(get_cov_IApN_0pi(), data, *((TH1D *)hist), 1e-6);
}

double do_chi2_dalphat_0pi(TH1 *hist) {
  auto data = get_dalphat_hist_0pi();

  std::cout << std::endl;
  return chi2(get_cov_dalphat_0pi(), data, *((TH1D *)hist), 1e-6);
}

TMatrixT<double> get_cov_dalphat_pi0() {
  // unit: 10$^{-47}$ cm$^2$/degree/nucleon
  // to: 10$^{-38}$ cm$^2$/degree/nucleon
  // scale: 10^9
  // const double matrix_element_decomp[9][9]{
  //     {101615, 78825, 47088, 23060, 11941, 944, -12437, -40048, -68951},
  //     {0, 38045, 54912, 34525, 30544, 12673, -945, -29507, -15789},
  //     {0, 0, 40833, 43038, 18788, 4871, 7983, 1325, -12945},
  //     {0, 0, 0, 50742, 55464, 31852, 13218, 19686, 11099},
  //     {0, 0, 0, 0, 47766, 59664, 50954, 8604, -45693},
  //     {0, 0, 0, 0, 0, 67023, 57869, 40189, 43947},
  //     {0, 0, 0, 0, 0, 0, 89938, 126346, 65696},
  //     {0, 0, 0, 0, 0, 0, 0, 108226, 151395},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 188252}};

  // TMatrixD cov_matrix_decomp(9, 9, &matrix_element_decomp[0][0]);

  // TMatrixD cov_matrix_transpose(cov_matrix_decomp);
  // cov_matrix_transpose.T();

  // auto cov_matrix = cov_matrix_transpose * cov_matrix_decomp;
  auto cov_from_root =
      (TMatrixD *)(((TList *)pizero_file.Get("dalphat"))->At(2));
  TMatrixD cov_matrix(9, 9);
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      cov_matrix[i][j] = (*cov_from_root)[i + 1][j + 1] / 1e-38 / 1e-38;
    }
  }
  // cov_matrix *= 1e-18;
  return cov_matrix;
}

TMatrixT<double> get_cov_dalphat_0pi() {
  // unit: 10$^{-46}$ cm$^2$/degree/nucleon
  // to: 10$^{-38}$ cm$^2$/degree/nucleon
  // scale: 10^8
  // const double matrix_element_decomp[12][12]{
  //     {15269, 9974, 7822, 4981, 5516, 4194, 5258, 4269, 3698, -2120, -4297,
  //      -6432},
  //     {0, 5998, 5956, 2869, 3617, 4716, 5519, 6098, 5834, 6523, 4562, 1998},
  //     {0, 0, 4646, 3945, 3042, 3442, 3941, 4653, 5042, 6819, 7571, 5498},
  //     {0, 0, 0, 5426, 4916, 4502, 5518, 7390, 9387, 12328, 12145, 13560},
  //     {0, 0, 0, 0, 5744, 4895, 3732, 3383, 3211, 4532, 3936, 3785},
  //     {0, 0, 0, 0, 0, 6887, 6703, 5879, 6785, 10476, 11622, 12575},
  //     {0, 0, 0, 0, 0, 0, 8719, 8842, 5503, 5707, 7951, 9020},
  //     {0, 0, 0, 0, 0, 0, 0, 10246, 10523, 8504, 11259, 10493},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 11181, 10052, 6936, 7588},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 14122, 10590, 7447},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15727, 11176},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16615}};

  // TMatrixD cov_matrix_decomp(12, 12, &matrix_element_decomp[0][0]);

  // TMatrixD cov_matrix_transpose(cov_matrix_decomp);
  // cov_matrix_transpose.T();

  // auto cov_matrix = cov_matrix_transpose * cov_matrix_decomp;
  // cov_matrix *= 1e-16;
  auto cov_from_root =
      (TMatrixD *)(((TList *)zeropi_file.Get("dalphat"))->At(2));
  TMatrixD cov_matrix(12, 12);
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      cov_matrix[i][j] = (*cov_from_root)[i + 2][j + 2] / 1e-38 / 1e-38;
    }
  }
  return cov_matrix;
}

TMatrixT<double> get_cov_IApN_0pi() {
  // unit: 10$^{-45}$ cm$^2$/GeV/c/nucleon
  // to: 10$^{-38}$ cm$^2$/GeV/c/nucleon
  // scale: 1e7
  // const double matrix_element_decomp[24][24]{
  //     {78769,  174422, 209919, 206290, 194890, 189406, 197725, 196553,
  //      139566, 109827, 116924, 102307, 85769,  66058,  41484,  47191,
  //      58203,  64875,  49373,  51054,  39645,  19856,  7572,   595},
  //     {0,      95257, 149500, 203630, 204702, 163306, 146870, 147512,
  //      101451, 87780, 86803,  88762,  77572,  47036,  35022,  33242,
  //      32945,  54116, 34377,  31908,  32224,  18982,  10262,  1330},
  //     {0,      0,      118376, 241382, 268375, 211168, 193736, 199308,
  //      143288, 118292, 90699,  94524,  101612, 69481,  55901,  57900,
  //      56178,  58814,  50271,  42849,  32875,  21295,  6855,   1051},
  //     {0,      0,     0,     217825, 328272, 258097, 180467, 174276,
  //      143946, 97893, 86490, 91738,  74965,  27092,  18254,  30748,
  //      36208,  50368, 40092, 36269,  35485,  21970,  9939,   1617},
  //     {0,      0,     0,     0,     178339, 296728, 309956, 291245,
  //      177210, 78767, 52640, 63145, 96850,  84434,  73314,  49559,
  //      48932,  46691, 26961, 20524, 12686,  2936,   -578,   -362},
  //     {0,      0,     0,     0,     0,     227369, 315481, 351720,
  //      224802, 80760, 32107, 43326, 77748, 106859, 81674,  53402,
  //      54078,  35547, 26637, 24185, 3461,  -5834,  -5987,  -1816},
  //     {0,      0,      0,     0,     0,     0,     220511, 359171,
  //      249826, 138874, 57582, 69681, 85290, 66140, 71528,  54537,
  //      61506,  52029,  37429, 36181, 24915, 13120, 2561,   220},
  //     {0,      0,      0,      0,      0,      0,     0,     276052,
  //      271536, 107865, 5516,   29187,  49169,  33609, 12840, 3330,
  //      -5640,  -8742,  -11069, -15027, -11514, -9147, -3826, -789},
  //     {0,      0,      0,      0,     0,     0,      0,     0,
  //      195784, 228959, 160652, 82617, 75452, 102961, 93896, 70822,
  //      63451,  48654,  40953,  40130, 20214, 5942,   -1205, -225},
  //     {0,     0,      0,      0,      0,      0,     0,     0,
  //      0,     176248, 210358, 164827, 125543, 98976, 93750, 61460,
  //      48282, 43683,  31388,  16967,  8548,   7822,  2697,  731},
  //     {0,     0,     0,      0,      0,     0,     0,     0,
  //      0,     0,     181146, 178319, 82289, 66379, 51427, 40503,
  //      21714, 18544, 14485,  14330,  5596,  -1935, -878,  157},
  //     {0,     0,     0,     0,      0,      0,     0,     0,
  //      0,     0,     0,     224805, 203754, 97867, 85404, 65175,
  //      40901, 29407, 26259, 19349,  8993,   177,   -1543, -89},
  //     {0,     0,     0,     0,     0,      0,      0,      0,
  //      0,     0,     0,     0,     213632, 191298, 127437, 86440,
  //      73810, 55123, 38454, 25320, 4526,   2718,   -3315,  83},
  //     {0,     0,     0,     0,     0,     0,      0,      0,
  //      0,     0,     0,     0,     0,     229394, 114497, 52171,
  //      28745, 31108, 25990, 19741, 14142, -165,   -1443,  -187},
  //     {0, 0, 0,      0,     0,     0,     0,     0,     0,    0,    0,    0,
  //      0, 0, 175591, 80134, 38065, 26100, 19301, 15450, 3305, 5811, -572,
  //      95},
  //     {0, 0, 0, 0,      0,     0,    0,     0,     0,     0,    0,    0,
  //      0, 0, 0, 128593, 62977, 8276, 16154, 16167, 14883, 7778, 1404, 405},
  //     {0, 0, 0, 0, 0,      0,     0,    0,     0,     0,    0,    0,
  //      0, 0, 0, 0, 126065, 78202, 9931, 10381, 11863, 8304, -153, -91},
  //     {0, 0, 0, 0, 0, 0,      0,     0,     0,     0,     0,    0,
  //      0, 0, 0, 0, 0, 128285, 65333, 14707, 17227, 16746, 8112, 1332},
  //     {0, 0, 0, 0, 0, 0, 0,      0,     0,    0,    0,    0,
  //      0, 0, 0, 0, 0, 0, 103993, 72726, 8846, 5920, 2227, 449},
  //     {0, 0, 0, 0, 0, 0, 0, 0,     0,     0,     0,    0,
  //      0, 0, 0, 0, 0, 0, 0, 97941, 61163, 12648, 7220, 870},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0,     0,     0,    0,
  //      0, 0, 0, 0, 0, 0, 0, 0, 84729, 26256, 6933, 1400},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0,     0,
  //      0, 0, 0, 0, 0, 0, 0, 0, 0, 61448, 11919, 1362},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     0,
  //      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 37477, 2167},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  //      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7435}};

  // TMatrixD cov_matrix_decomp(24, 24, &matrix_element_decomp[0][0]);

  // TMatrixD cov_matrix_transpose(cov_matrix_decomp);
  // cov_matrix_transpose.T();

  // auto cov_matrix = cov_matrix_transpose * cov_matrix_decomp;
  // cov_matrix *= 1e-14;
  auto cov_from_root =
      (TMatrixD *)(((TList *)zeropi_file.Get("neutronmomentum"))->At(2));
  TMatrixD cov_matrix(24, 24);
  for (int i = 0; i < 24; i++) {
    for (int j = 0; j < 24; j++) {
      cov_matrix[i][j] = (*cov_from_root)[i + 2][j + 2] / 1e-38 / 1e-38;
    }
  }
  return cov_matrix;
}

TMatrixT<double> get_cov_IApN_pi0() {
  // unit: 10$^{-47}$ cm$^2$/MeV/c/nucleon
  // to: 10$^{-38}$ cm$^2$/GeV/c/nucleon
  // scale: 10^9 / 10e3 = 10^6
  // const double matrix_element_decomp[12][12]{
  //     {12714, 16637, 7407, 2005, 1829, 1418, -822, -2963, -1775, -2892,
  //     -3822,
  //      -6706},
  //     {0, 15026, 18731, 13567, 1839, 781, -455, -1971, -366, -4464, -2025,
  //      -6245},
  //     {0, 0, 21091, 31461, 10158, 4100, 2114, 1556, 667, -1750, -2939, -635},
  //     {0, 0, 0, 25067, 18613, 10259, 5284, 949, 1992, 4051, 4716, -4554},
  //     {0, 0, 0, 0, 17620, 17591, 8716, 1303, 1944, -1240, -1948, -5097},
  //     {0, 0, 0, 0, 0, 16880, 21777, 17320, 4843, 184, 1668, 5120},
  //     {0, 0, 0, 0, 0, 0, 16465, 15579, 13441, 7898, 6227, -1467},
  //     {0, 0, 0, 0, 0, 0, 0, 18517, 19628, 13906, 4144, 5067},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 26267, 16989, 8748, -1683},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 25355, 21168, 14467},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18649, 7763},
  //     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25157}};

  // TMatrixD cov_matrix_decomp(12, 12, &matrix_element_decomp[0][0]);

  // TMatrixD cov_matrix_transpose(cov_matrix_decomp);
  // cov_matrix_transpose.T();

  // auto cov_matrix = cov_matrix_transpose * cov_matrix_decomp;
  // cov_matrix *= 1e-12;
  auto cov_from_root =
      (TMatrixD *)(((TList *)pizero_file.Get("neutronmomentum"))->At(2));
  TMatrixD cov_matrix(12, 12);
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      cov_matrix[i][j] = (*cov_from_root)[i + 1][j + 1] / 1e-38 / 1e-38 * 1e6;
    }
  }
  return cov_matrix;
}
