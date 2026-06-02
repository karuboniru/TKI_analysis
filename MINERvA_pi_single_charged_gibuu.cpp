#include "EvtTracker2event.h"
#include "dochi2.h"
#include "plottools.hxx"
#include "tkievent.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <TFile.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>
#include <boost/program_options.hpp>

#include <cmath>
#include <memory>
#include <print>
#include <sstream>
#include <string>
#include <vector>

ROOT::RDF::RNode vars_define(ROOT::RDF::RNode df) {
  return df
      .Define("Q2",
              [](const TLorentzVector &InitNeutrino,
                 const TLorentzVector &PrimaryLepton) {
                auto q0 = InitNeutrino - PrimaryLepton;
                return -q0.Mag2();
              },
              {"InitNeutrino", "PrimaryLepton"})
      .Define("W",
              [](const TLorentzVector &InitNeutrino,
                 const TLorentzVector &PrimaryLepton,
                 const TLorentzVector &InitNucleon) {
                auto had_system = InitNucleon + InitNeutrino - PrimaryLepton;
                return had_system.M();
              },
              {"InitNeutrino", "PrimaryLepton", "InitNucleon"});
}

struct TargetResult {
  TH1D h_Tk{};
  TH1D h_theta{};
  size_t n_runs{};
  size_t n_events{};
};

TargetResult process_target(const std::vector<std::string> &files,
                            const std::array<double, 8> &Tk_bin_edges,
                            const std::array<double, 15> &theta_bin_edges) {
  ROOT::RDataFrame input("out_tree", files);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(input);
#endif
  auto d = vars_define(TrackerPrepare(input));
  auto n_events = input.Count();

  auto d_signal =
      d.Filter([](double W) { return W < 1.4; }, {"W"}, "W < 1.4 GeV")
          .Filter(
              [](const NeutrinoEvent &event) {
                return event.count_post(211) + event.count_post(-211) == 1;
              },
              {"EventRecord"}, "single charged pion")
          .Define("pion_p4",
                  [](const NeutrinoEvent &event) {
                    for (auto &&[id, particle] : event.post_range(211))
                      return particle;
                    for (auto &&[id, particle] : event.post_range(-211))
                      return particle;
                    return TLorentzVector{};
                  },
                  {"EventRecord"})
          .Define("pion_Tk",
                  [](const TLorentzVector &pion) {
                    return (pion.E() - pion.M()) * 1000.;
                  },
                  {"pion_p4"})
          .Define("pion_theta",
                  [](const TLorentzVector &pion) {
                    return pion.Theta() * 180. / M_PI;
                  },
                  {"pion_p4"});

  auto h_Tk_tmp =
      d_signal.Histo1D({"h_Tk", "h_Tk", 7, Tk_bin_edges.data()}, "pion_Tk",
                       "weight");
  auto h_theta_tmp = d_signal.Histo1D(
      {"h_theta", "h_theta", 14, theta_bin_edges.data()}, "pion_theta",
      "weight");

  TargetResult result;
  result.h_Tk = *(TH1D *)h_Tk_tmp->Clone("h_Tk");
  result.h_theta = *(TH1D *)h_theta_tmp->Clone("h_theta");
  result.n_runs = files.size();
  result.n_events = *n_events;
  return result;
}

std::pair<TH1D, TMatrixDSym>
build_Tk_data(const std::array<double, 8> &Tk_bin_edges) {
  constexpr double Tk_bin_minerva[]{1.2, 1.3, 1.2, 0.9, 0.8, 0.7, 0.4};
  constexpr double Tk_bin_minerva_error[]{0.26, 0.21, 0.18, 0.17,
                                          0.18, 0.17, 0.22};

  TH1D h_data("h_Tk_data", "h_Tk_data", 7, Tk_bin_edges.data());
  for (int i = 0; i < 7; ++i) {
    h_data.SetBinContent(i + 1, Tk_bin_minerva[i] * 1e-3);
    h_data.SetBinError(i + 1, Tk_bin_minerva_error[i] * Tk_bin_minerva[i] *
                                  1e-3);
  }

  constexpr double cov_TK[7][7]{
      {1.0000, 0.7178, 0.6969, 0.6691, 0.6598, 0.5719, 0.5313},
      {0.7178, 1.0000, 0.8679, 0.8090, 0.7962, 0.6995, 0.6617},
      {0.6969, 0.8679, 1.0000, 0.8266, 0.8137, 0.7249, 0.6626},
      {0.6691, 0.8090, 0.8266, 1.0000, 0.8823, 0.8166, 0.7830},
      {0.6598, 0.7962, 0.8137, 0.8823, 1.0000, 0.8295, 0.8024},
      {0.5719, 0.6995, 0.7249, 0.8166, 0.8295, 1.0000, 0.8985},
      {0.5313, 0.6617, 0.6626, 0.7830, 0.8024, 0.8985, 1.0000}};

  TMatrixDSym cov(7);
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      cov[i][j] = cov_TK[i][j] * Tk_bin_minerva_error[i] *
                  Tk_bin_minerva_error[j] * Tk_bin_minerva[i] *
                  Tk_bin_minerva[j] * 1e-6;
  return {h_data, cov};
}

std::pair<TH1D, TMatrixDSym>
build_theta_data(const std::array<double, 15> &theta_bin_edges) {
  constexpr double theta_bin_minerva[]{1.2, 2.3, 2.8, 3.8, 3.6, 3.0, 2.3,
                                       2.2, 1.7, 1.0, 0.8, 0.6, 0.4, 0.2};
  constexpr double theta_bin_minerva_error[]{0.23, 0.21, 0.20, 0.20, 0.20,
                                             0.20, 0.21, 0.20, 0.19, 0.21,
                                             0.19, 0.19, 0.21, 0.26};

  TH1D h_data("h_theta_data", "h_theta_data", 14, theta_bin_edges.data());
  for (int i = 0; i < 14; ++i) {
    h_data.SetBinContent(i + 1, theta_bin_minerva[i] * 1e-3);
    h_data.SetBinError(i + 1,
                       theta_bin_minerva_error[i] * theta_bin_minerva[i] *
                           1e-3);
  }

  constexpr double cov_theta[14][14]{
      {1.0000, 0.8411, 0.8344, 0.8561, 0.8561, 0.8437, 0.8104, 0.7787, 0.8562,
       0.7050, 0.7002, 0.6942, 0.5834, 0.4134},
      {0.8411, 1.0000, 0.8468, 0.8617, 0.8606, 0.8438, 0.8254, 0.7975, 0.8330,
       0.6916, 0.6979, 0.6959, 0.6005, 0.4242},
      {0.8344, 0.8468, 1.0000, 0.8735, 0.8608, 0.8527, 0.8377, 0.8141, 0.8475,
       0.6781, 0.6873, 0.6856, 0.5911, 0.4354},
      {0.8561, 0.8617, 0.8735, 1.0000, 0.8897, 0.8810, 0.8574, 0.8376, 0.8573,
       0.7041, 0.7238, 0.7180, 0.6152, 0.4589},
      {0.8561, 0.8606, 0.8608, 0.8897, 1.0000, 0.8924, 0.8637, 0.8317, 0.8607,
       0.7207, 0.7273, 0.7261, 0.6309, 0.4594},
      {0.8437, 0.8438, 0.8527, 0.8810, 0.8924, 1.0000, 0.8564, 0.8264, 0.8589,
       0.7070, 0.7171, 0.7169, 0.6203, 0.4633},
      {0.8104, 0.8254, 0.8377, 0.8574, 0.8637, 0.8564, 1.0000, 0.8160, 0.8358,
       0.7020, 0.7147, 0.7176, 0.6325, 0.4621},
      {0.7787, 0.7975, 0.8141, 0.8376, 0.8317, 0.8264, 0.8160, 1.0000, 0.8095,
       0.7047, 0.7304, 0.7186, 0.6330, 0.5031},
      {0.8562, 0.8330, 0.8475, 0.8573, 0.8607, 0.8589, 0.8358, 0.8095, 1.0000,
       0.7619, 0.7509, 0.7447, 0.6452, 0.4745},
      {0.7050, 0.6916, 0.6781, 0.7041, 0.7207, 0.7070, 0.7020, 0.7047, 0.7619,
       1.0000, 0.7444, 0.7208, 0.6442, 0.4579},
      {0.7002, 0.6979, 0.6873, 0.7238, 0.7273, 0.7171, 0.7147, 0.7304, 0.7509,
       0.7444, 1.0000, 0.7534, 0.6670, 0.4832},
      {0.6942, 0.6959, 0.6856, 0.7180, 0.7261, 0.7169, 0.7176, 0.7186, 0.7447,
       0.7208, 0.7534, 1.0000, 0.6555, 0.4646},
      {0.5834, 0.6005, 0.5911, 0.6152, 0.6309, 0.6203, 0.6325, 0.6330, 0.6452,
       0.6442, 0.6670, 0.6555, 1.0000, 0.4535},
      {0.4134, 0.4242, 0.4354, 0.4589, 0.4594, 0.4633, 0.4621, 0.5031, 0.4745,
       0.4579, 0.4832, 0.4646, 0.4535, 1.0000}};

  TMatrixDSym cov(14);
  for (int i = 0; i < 14; ++i)
    for (int j = 0; j < 14; ++j)
      cov[i][j] = cov_theta[i][j] * theta_bin_minerva_error[i] *
                  theta_bin_minerva_error[j] * theta_bin_minerva[j] *
                  theta_bin_minerva[j] * 1e-6;
  return {h_data, cov};
}

int main(int argc, char *argv[]) {
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()("input-carbon", po::value<std::vector<std::string>>(),
                     "GiBUU ROOT files for carbon target")(
      "input-hydrogen", po::value<std::vector<std::string>>(),
      "GiBUU ROOT files for hydrogen target")(
      "output", po::value<std::string>()->default_value("output_pi_charged.root"),
      "Output ROOT file")("help", "produce help message");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help") ||
      (!vm.count("input-carbon") && !vm.count("input-hydrogen"))) {
    std::cout << desc << std::endl;
    return 1;
  }

  TH1::AddDirectory(false);
  ROOT::EnableImplicitMT();

  constexpr std::array<double, 8> Tk_bin_edges{35., 55., 75., 100.,  125.,
                                                150, 200, 350};
  constexpr std::array<double, 15> theta_bin_edges{
      0, 15, 22, 29, 36, 43, 50, 57, 72, 108, 130, 140, 150, 165, 180};
  constexpr int Tk_nbins = 7;
  constexpr int theta_nbins = 14;
  constexpr double C_nucleon = 12.;
  constexpr double H_nucleon = 1.;
  constexpr double CH_total_nucleon = C_nucleon + H_nucleon;

  auto [h_Tk_data, cov_Tk] = build_Tk_data(Tk_bin_edges);
  auto [h_theta_data, cov_theta] = build_theta_data(theta_bin_edges);

  std::vector<std::string> carbon_files, hydrogen_files;
  if (vm.count("input-carbon"))
    carbon_files = vm["input-carbon"].as<std::vector<std::string>>();
  if (vm.count("input-hydrogen"))
    hydrogen_files = vm["input-hydrogen"].as<std::vector<std::string>>();

  TargetResult result_C{}, result_H{};

  if (!carbon_files.empty()) {
    result_C = process_target(carbon_files, Tk_bin_edges, theta_bin_edges);
    std::println("Carbon: {} events (pre-cut) from {} files, n_runs={}",
                 result_C.n_events, carbon_files.size(), result_C.n_runs);
  }
  if (!hydrogen_files.empty()) {
    result_H = process_target(hydrogen_files, Tk_bin_edges, theta_bin_edges);
    std::println("Hydrogen: {} events (pre-cut) from {} files, n_runs={}",
                 result_H.n_events, hydrogen_files.size(), result_H.n_runs);
  }

  if (carbon_files.empty())
    std::println("WARNING: No carbon input");
  if (hydrogen_files.empty())
    std::println("WARNING: No hydrogen input");

  auto build_CH_hist = [&](int nbins, const double *bins, const char *name,
                           bool is_Tk) {
    TH1D h(name, name, nbins, bins);
    h.SetDirectory(nullptr);
    h.Sumw2(false);
    const int n_runs_C = carbon_files.size();
    const int n_runs_H = hydrogen_files.size();
    for (int i = 1; i <= nbins; ++i) {
      double sum = 0;
      if (n_runs_C > 0) {
        const auto &src = is_Tk ? result_C.h_Tk : result_C.h_theta;
        sum += C_nucleon * src.GetBinContent(i) / n_runs_C;
      }
      if (n_runs_H > 0) {
        const auto &src = is_Tk ? result_H.h_Tk : result_H.h_theta;
        sum += H_nucleon * src.GetBinContent(i) / n_runs_H;
      }
      h.SetBinContent(i, sum / CH_total_nucleon / 10.);
    }
    return h;
  };

  auto h_Tk_CH =
      build_CH_hist(Tk_nbins, Tk_bin_edges.data(), "h_Tk_CH", true);
  auto h_theta_CH =
      build_CH_hist(theta_nbins, theta_bin_edges.data(), "h_theta_CH", false);

  h_Tk_CH.Scale(1., "WIDTH");
  h_theta_CH.Scale(1., "WIDTH");

  auto chi2_Tk = chi2(cov_Tk, h_Tk_data, h_Tk_CH, 1e-2);
  auto chi2_theta = chi2(cov_theta, h_theta_data, h_theta_CH, 1e-2);

  std::println("chi2_Tk = {}", chi2_Tk);
  std::println("chi2_theta = {}", chi2_theta);

  auto output_file = std::make_unique<TFile>(
      vm["output"].as<std::string>().c_str(), "RECREATE");
  h_Tk_CH.Write("h_Tk_CH");
  h_theta_CH.Write("h_theta_CH");
  h_Tk_data.Write("h_Tk_data");
  h_theta_data.Write("h_theta_data");

  auto make_gi_label = [](double chi2, int nbins) {
    auto tex = std::make_unique<TLatex>();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    std::stringstream ss;
    ss << "GiBUU  #chi^{2}/NDF = " << std::round(chi2) << "/" << nbins;
    tex->SetTitle(ss.str().c_str());
    return tex;
  };

  auto gi_Tk = make_gi_label(chi2_Tk, Tk_nbins);
  auto gi_theta = make_gi_label(chi2_theta, theta_nbins);

  do_plot({&h_Tk_data, (TH1 *)&h_Tk_CH, gi_Tk.get()},
          "Tk_pi_charged",
          "d#sigma/dT_{#pi} (10^{-38} cm^{2}/MeV/nucleon)",
          "Pion kinetic energy (MeV)",
          {0.50, 0.65, 0.80, 0.85}, 350., "MINERvA #pi^{#pm} Data", "HIST C",
          0.050, {.top = 0.04, .bottom = 0.12});

  do_plot({&h_theta_data, (TH1 *)&h_theta_CH, gi_theta.get()},
          "theta_pi_charged",
          "d#sigma/d#theta_{#pi} (10^{-38} cm^{2}/degree/nucleon)",
          "Pion angle (degree)",
          {0.14, 0.65, 0.45, 0.85}, 0., "MINERvA #pi^{#pm} Data", "HIST C",
          0.040, {.top = 0.06, .bottom = 0.11});

  output_file->Close();

  return 0;
}
