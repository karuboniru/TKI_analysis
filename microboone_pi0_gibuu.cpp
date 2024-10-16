#include "EvtTracker2event.h"
#include "MINERvA_tki_cut.h"
#include "dochi2.h"
#include "expdata.h"
#include "plottools.hxx"
#include "tkievent.h"
#include "tkigeneral.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <TH1.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMatrixDfwd.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>
#include <algorithm>
#include <boost/program_options.hpp>

#include <cmath>
#include <print>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

plot_data get_info(const std::string &varname) {
  if (varname == "dalphat") {
    return plot_data{.bins = 30,
                     .xmin = 0.,
                     .xmax = 180.,
                     .name = "#delta#it{#alpha}_{T} (degree)",
                     .ytitle = "d#sigma/d#delta#it{#alpha}_{T} (10^{#minus 38} "
                               "cm^{2}/degree/nucleon)"};
  }
  if (varname == "IApN") {
    return plot_data{.bins = 30,
                     .xmin = 0.,
                     .xmax = 0.8,
                     .name = "#it{p}_{N} (GeV/#it{c})",
                     .ytitle = "d#sigma/d#it{p}_{N} (10^{#minus 38} "
                               "cm^{2}/GeV/#it{c}/nucleon)"};
  }
  if (varname == "dpL") {
    return plot_data{
        .bins = 30, .xmin = -1., .xmax = 1., .name = "#delta p_{L}"};
  }
  if (varname == "W") {
    return plot_data{.bins = 60,
                     .xmin = 0.7,
                     .xmax = 2.8,
                     .name = "#it{W} (GeV)",
                     .ytitle =
                         "d#sigma/d#it{W} (10^{#minus 38} cm^{2}/GeV/nucleon)",
                     .ymax_0pi = 0.15};
  }

  return {};
}

template <typename T, size_t N>
ROOT::RDF::RResultPtr<TH1> make_plots(T &&df_in, std::array<double, N> bins,
                                      std::string variable,
                                      const std::string &postfix = "") {
  auto plot_data = get_info(variable);
  auto plot_name = postfix.empty() ? variable : (variable + "_" + postfix);
  return df_in.template Histo1D<double>(
      {plot_name.c_str(), variable.c_str(), bins.size() - 1, bins.data()},
      variable, "weight");
}

template <typename T>
ROOT::RDF::RResultPtr<TH1> make_plots(T &&df_in, std::string variable,
                                      const std::string &postfix = "") {
  auto plot_data = get_info(variable);
  auto plot_name = postfix.empty() ? variable : (variable + "_" + postfix);
  return df_in.template Histo1D<double>({plot_name.c_str(), variable.c_str(),
                                         plot_data.bins, plot_data.xmin,
                                         plot_data.xmax},
                                        variable, "weight");
}

TH1D smear_prediction(TH1D hist, const TMatrixDSym &smear) {
  auto n_bins = hist.GetNbinsX();
  auto smear_hist = hist;
  smear_hist.Reset();
  smear_hist.SetName(std::format("{}_smeared", hist.GetName()).c_str());
  // Convention from MicroBooNE: Last bin contains overflow
  hist.SetBinContent(n_bins, hist.GetBinContent(n_bins) +
                                 hist.GetBinContent(n_bins + 1));
  for (int hist_bin_id = 1; hist_bin_id <= n_bins; ++hist_bin_id) {
    double sum = 0;
    for (int true_bin_id = 1; true_bin_id <= n_bins; ++true_bin_id) {
      sum += smear[true_bin_id - 1][hist_bin_id - 1] *
             hist.GetBinContent(true_bin_id);
    }
    smear_hist.SetBinContent(hist_bin_id, sum);
  }
  return smear_hist;
}

const std::vector<std::tuple<std::string, std::function<bool(int)>>>
    list_channel_channeldef{
        {"2#pi BG", [](int id) { return id == 37; }},
        {"1#pi BG", [](int id) { return id == 32 || id == 33; }},
        {"2p2h", [](int id) { return id == 35; }},
        {"QE", [](int id) { return id == 1; }},
        {"RES", [](int id) { return id >= 2 && id <= 31; }},
        {"DIS", [](int id) { return id == 34; }},
    };

template <typename T, size_t N = 0>
auto plot_channels(T df_in, std::string variable,
                   std::array<double, N> bins = std::array<double, 0>{}) {
  return list_channel_channeldef | std::views::enumerate |
         std::views::transform([&](auto &&name_id_count) {
           auto &[count, name_id] = name_id_count;
           auto &[name, id] = name_id;
           if constexpr (N != 0)
             return std::make_tuple(name,
                                    make_plots(df_in.Filter(id, {"channel"}),
                                               bins, variable, name + "_0pi"),
                                    count * 2 - (count > 0) - (count > 1));
           else
             return std::make_tuple(name,
                                    make_plots(df_in.Filter(id, {"channel"}),
                                               variable, name + "_0pi"),
                                    count * 2 - (count > 0) - (count > 1));
         }) |
         std::ranges::to<std::vector>();
}

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()("input-files", po::value<std::vector<std::string>>(),
                     "Input files")(
      "add-text", po::value<std::string>()->default_value(""),
      "Additional text")("help", "produce help message");
  po::positional_options_description p;
  p.add("input-files", -1);
  po::variables_map vm;
  po::store(
      po::command_line_parser(argc, argv).options(desc).positional(p).run(),
      vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-files")) {
    std::cout << desc << std::endl;
    return 1;
  }

  auto &&files = vm["input-files"].as<std::vector<std::string>>();
  auto n_runs = files.size();
  auto additional_text = vm["add-text"].as<std::string>();

  TH1::AddDirectory(false);
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame input("out_tree", files);

  // GiBUU data preparation
  auto d = TrackerPrepare(input);
  ROOT::RDF::Experimental::AddProgressBar(d);

  auto data_cut =
      d.Filter(
           [](const NeutrinoEvent &e) {
             // exactly one muon and one pion
             if (e.count_post(13) != 1 || e.count_post(111) != 1) {
               return false;
             }
             // and no other particles than muon, pion, proton, neutron
             return std::ranges::all_of(e.get_ids_post(), [](int pdg) {
               switch (pdg) {
               case 13:
               case 111:
               case 2212:
               case 2112:
                 return true;
               default:
                 return false;
               }
             });
           },
           {"EventRecord"}, "topocut")
          .Define("pion_mom",
                  [](const NeutrinoEvent &e) {
                    return e.post_range(111).begin()->second;
                  },
                  {"EventRecord"})
          // .Define(
          //     "pion_Tk",
          //     [](TLorentzVector &pion) { return (pion.E() - pion.M()) * 1e3;
          //     },
          //     {"pion_mom"})
          .Define("pion_p", [](TLorentzVector &pion) { return pion.P(); },
                  {"pion_mom"})
          // .Define("pion_angle",
          //         [](TLorentzVector &pion) {
          //           //  auto pion_dir = pion.Vect().Unit();
          //           auto pion_angle = pion.Vect().Angle(TVector3{0, 0, 1.});
          //           // convert to degrees
          //           pion_angle *= 180 / M_PI;
          //           return pion_angle;
          //         },
          //         {"pion_mom"})
          .Define("costh_pi",
                  [](TLorentzVector &pion) { return pion.Vect().CosTheta(); },
                  {"pion_mom"})
          .Define("W",
                  [](const TLorentzVector &InitNeutrino,
                     const TLorentzVector &PrimaryLepton,
                     const TLorentzVector &InitNucleon) {
                    auto had_system =
                        InitNucleon + InitNeutrino - PrimaryLepton;
                    return had_system.M();
                  },
                  {"InitNeutrino", "PrimaryLepton", "InitNucleon"})
      // .Define("proton",
      //         [](NeutrinoEvent &e) {
      //           ROOT::RVec<TLorentzVector> p;
      //           if (e.count_out(2212)) {
      //             p.push_back(e.get_leading(2212));
      //           }
      //           return p;
      //         },
      //         {"EventRecord"})
      // .Define("protonangle",
      //         [](ROOT::RVec<TLorentzVector> &p) {
      //           ROOT::RVec<double> var;
      //           for (auto &k : p) {
      //             var.push_back(k.Vect().Angle(TVector3{0, 0, 1.}) *
      //                           (180 / M_PI));
      //           }
      //           return var;
      //         },
      //         {"proton"})
      // .Define("protonmomentum",
      //         [](ROOT::RVec<TLorentzVector> &p) {
      //           ROOT::RVec<double> var;
      //           for (auto &k : p) {
      //             var.push_back(k.P());
      //           }
      //           return var;
      //         },
      //         {"proton"})
      ;

  auto data_cut_report = data_cut.Report();

  auto pi0_momentum_binning = MicroBooNE::pi0_momentum::get_binning();
  auto pi0_angular_binning = MicroBooNE::pi0_angular::get_binning();

  auto pi0_plot_list = plot_channels(data_cut, "pion_p", pi0_momentum_binning);
  auto pi0_angular_plot_list =
      plot_channels(data_cut, "costh_pi", pi0_angular_binning);

  auto pi0_momentum_hist = data_cut.Histo1D(
      {"pi0_momentum_hist", "", MicroBooNE::pi0_momentum::dimension,
       pi0_momentum_binning.data()},
      "pion_p", "weight");
  auto pi0_angular_hist = data_cut.Histo1D({"pi0_angular_hist", "",
                                            MicroBooNE::pi0_angular::dimension,
                                            pi0_angular_binning.data()},
                                           "costh_pi", "weight");

  auto vars = std::to_array({"W"});
  auto var_plots = vars | std::views::transform([&](std::string name) {
                     return plot_channels(data_cut, std::move(name));
                   }) |
                   std::ranges::to<std::vector>();

  ////////////////////////

  std::cout << "Finished data preparation\n";

  ////////////////////////
  auto pi0_momentum_hist_smeared = smear_prediction(
      pi0_momentum_hist.GetValue(), MicroBooNE::pi0_momentum::get_smear());
  auto pi0_angular_hist_smeared = smear_prediction(
      pi0_angular_hist.GetValue(), MicroBooNE::pi0_angular::get_smear());
  {
    auto scale_factor = 1. / n_runs / 10;
    pi0_momentum_hist_smeared.Scale(scale_factor, "WIDTH");
    pi0_angular_hist_smeared.Scale(scale_factor, "WIDTH");
  }

  auto xsec_inted = pi0_angular_hist_smeared.Integral("WIDTH");

  std::ranges::for_each(var_plots, [&](auto &&tup) {
    std::ranges::for_each(tup, [&](auto &&tup1) {
      // auto &&[_, hist, __] = tup1;
      auto &&hist = std::get<1>(tup1);
      hist->Scale(1. / n_runs / 10, "WIDTH");
    });
  });

  auto vars_stack = var_plots | std::views::transform([&](auto &&list) {
                      return build_stack_from_list(list, 0.05 * xsec_inted);
                    }) |
                    std::ranges::to<std::vector>();

  // auto scale_list = [&](auto &&list) {
  //   for (auto &&[name, hist, _] : list) {
  //     hist->Scale(1. / n_runs / 10, "WIDTH");
  //   }
  // };

  auto pi0_plot_list_smeared =
      pi0_plot_list | std::views::transform([&](auto &&tup) {
        auto &&[name, hist, count] = tup;
        auto hist_smeared = std::make_unique<TH1D>(smear_prediction(
            *((TH1D *)hist.GetPtr()), MicroBooNE::pi0_momentum::get_smear()));
        hist_smeared->Scale(1. / n_runs / 10, "WIDTH");
        return std::make_tuple(name, std::move(hist_smeared), count);
      }) |
      std::ranges::to<std::vector>();

  auto pi0_angular_plot_list_smeared =
      pi0_angular_plot_list | std::views::transform([&](auto &&tup) {
        auto &&[name, hist, count] = tup;
        auto hist_smeared = std::make_unique<TH1D>(smear_prediction(
            *((TH1D *)hist.GetPtr()), MicroBooNE::pi0_angular::get_smear()));
        hist_smeared->Scale(1. / n_runs / 10, "WIDTH");
        return std::make_tuple(name, std::move(hist_smeared), count);
      }) |
      std::ranges::to<std::vector>();

  // scale_list(pi0_plot_list);
  // scale_list(pi0_angular_plot_list);

  auto &&[stack_momentum, legend_momentum] =
      build_stack_from_list(pi0_plot_list_smeared, 0.01 * xsec_inted);
  auto &&[stack_angular, legend_angular] =
      build_stack_from_list(pi0_angular_plot_list_smeared, 0.01 * xsec_inted);

  auto pi0_momentum_chi2 =
      MicroBooNE::pi0_momentum::do_chi2(&pi0_momentum_hist_smeared);
  auto pi0_angular_chi2 =
      MicroBooNE::pi0_angular::do_chi2(&pi0_angular_hist_smeared);

  std::println("Chi2 Results: \n"
               " - Momentum: {}\n"
               " - Angular: {}\n",
               pi0_momentum_chi2, pi0_angular_chi2);
  data_cut_report->Print();
  auto hist_momentum = MicroBooNE::pi0_momentum::get_hist();
  auto hist_angular = MicroBooNE::pi0_angular::get_hist();
  auto form_legend = [&](TH1 *hist, double chi2) -> std::string {
    std::stringstream ss;
    ss << "#chi^{2}/NDF = " << round(chi2) << "/" << hist->GetNbinsX();
    return ss.str();
  };
  std::unique_ptr<TLatex> latex;
  if (!additional_text.empty()) {
    latex = std::make_unique<TLatex>(0.65, 0.5, additional_text.c_str());
    latex->SetNDC();
  }
  do_plot({&hist_momentum, &pi0_momentum_hist_smeared, &stack_momentum,
           &legend_angular, latex.get()},
          "mom",
          "d#sigma/d#it{p}_{#pi^{0}} (#times 10^{#minus 38} "
          "cm^{2}/GeV/#it{c}/nucleon)",
          "#it{p}_{#pi^{0}} (GeV/#it{c})", {.65, .45, .95, .9}, 0.,
          form_legend(&hist_momentum, pi0_momentum_chi2), "HIST", 0,
          {.top = 0.04, .bottom = 0.12});

  do_plot({&hist_angular, &pi0_angular_hist_smeared, &stack_angular,
           &legend_angular, latex.get()},
          "ang",
          "d#sigma/dcos #theta_{#pi^{0}} (#times 10^{#minus 38} "
          "cm^{2}/rad/nucleon)",
          "cos #theta_{#pi^{0}} (rad)", {.15, .55, .55, .85}, 0.,
          form_legend(&hist_angular, pi0_angular_chi2), "HIST", 0,
          {.top = 0.06, .bottom = 0.12});

  for (auto &&[name, tup] : std::views::zip(vars, vars_stack)) {
    auto &&[stack, legend] = tup;
    auto plotobj = get_info(name);
    do_plot({&stack, &legend, latex.get()}, name, plotobj.ytitle, plotobj.name,
            {.7, .55, .9, .95}, plotobj.xmax, "MicroBooNE CC #pi^{0}", "HIST",
            plotobj.ymax_0pi, {.top = 0.04, .bottom = 0.12});
  }

  return 0;
}