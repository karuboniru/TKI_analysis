#include "EvtTracker2event.h"
#include "MINERvA_tki_cut.h"
#include "dochi2.h"
#include "plottools.hxx"
#include "tkigeneral.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMatrixDfwd.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>
#include <boost/program_options.hpp>

#include <nlohmann/json.hpp>

#include <cmath>
#include <functional>
#include <memory>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "expdata.h"

// this is general channel definition
const std::vector<std::tuple<std::string, std::function<bool(int)>>>
    list_channel_channeldef{
        {"QE", [](int id) { return id == 1; }},
        {"2p2h", [](int id) { return id == 35; }},
        {"RES", [](int id) { return id >= 2 && id <= 31; }},
        {"DIS", [](int id) { return id == 34; }},
        {"1#pi BG", [](int id) { return id == 32 || id == 33; }},
        {"2#pi BG", [](int id) { return id == 37; }},
    };

// pion count definition for pi0 channel
const std::vector<std::tuple<std::string, std::function<bool(size_t)>>>
    list_pion_cut{{"1pi0", [](size_t i) { return i == 1; }},
                  {"Mpi0", [](size_t i) { return i > 1; }}};

std::string pion_pretty_name(std::string s) {
  if (s == "1pi0") {
    return "1#pi^{#lower[0.2]{0}}";
  }
  if (s == "Mpi0") {
    return "multi-#pi^{#lower[0.2]{0}}";
  }
  return s;
}

plot_data get_info(std::string varname) {
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
  if (varname == "p_cos_theta") {
    return plot_data{
        .bins = 30, .xmin = -1., .xmax = 1., .name = "cos#theta_{p}"};
  }
  if (varname == "p_theta") {
    return plot_data{
        .bins = 30, .xmin = 0., .xmax = 180., .name = "#theta_{p}"};
  }
  if (varname == "dtl") {
    return plot_data{
        .bins = 30, .xmin = -1., .xmax = 1., .name = "cos #delta#alpha_{L}"};
  }
  if (varname == "Q2") {
    return plot_data{.bins = 50,
                     .xmin = 0.,
                     .xmax = 5.,
                     .name = "#it{Q}^{2} (GeV^{#kern[0.3]{2}})",
                     .ytitle = "d#sigma/d#it{Q}^{2} (10^{#minus 38} "
                               "cm^{2}/GeV^{#kern[0.3]{2}}/nucleon)",
                     .ymax_0pi = 0.35,
                     .ymax_pi0 = 0.15};
  }
  if (varname == "xBj") {
    return plot_data{.bins = 50,
                     .xmin = 0.,
                     .xmax = 2.0,
                     .name = "#it{x}_{Bj}",
                     .ytitle =
                         "d#sigma/d#it{x}_{Bj} (10^{#minus 38} cm^{2}/nucleon)",
                     .ymax_0pi = .6,
                     .ymax_pi0 = .3};
  }
  if (varname == "W") {
    return plot_data{.bins = 60,
                     .xmin = 0.7,
                     .xmax = 4.0,
                     .name = "#it{W} (GeV)",
                     .ytitle =
                         "d#sigma/d#it{W} (10^{#minus 38} cm^{2}/GeV/nucleon)",
                     .ymax_0pi = 0.5,
                     .ymax_pi0 = .35,
                     .xmax_0pi = 2.,
                     .xmax_pi0 = 4.};
  }
  if (varname == "p_mu") {
    return plot_data{
        .bins = 30,
        .xmin = 0.,
        .xmax = 12.,
        .name = "#it{p}_{#mu} (GeV/#it{c})",
        .ytitle =
            "d#sigma/d#it{p}_{#mu} (10^{#minus 38} cm^{2}/GeV/#it{c}/nucleon)"};
  }
  if (varname == "theta_mu") {
    return plot_data{
        .bins = 30,
        .xmin = 0.,
        .xmax = 60.,
        .name = "#theta_{#mu} (degree)",
        .ytitle =
            "d#sigma/d#theta_{#mu} (10^{#minus 38} cm^{2}/degree/nucleon)"};
  }
  if (varname == "p_p") {
    return plot_data{
        .bins = 30,
        .xmin = 0.,
        .xmax = 20.,
        .name = "#it{p}_{p} (GeV/#it{c})",
        .ytitle =
            "d#sigma/d#it{p}_{p} (10^{#minus 38} cm^{2}/GeV/#it{c}/nucleon)"};
  }
  if (varname == "theta_p") {
    return plot_data{
        .bins = 30,
        .xmin = 0.,
        .xmax = 180.,
        .name = "#theta_{p} (degree)",
        .ytitle = "d#sigma/d#theta_{p} (10^{#minus 38} cm^{2}/degree/nucleon)"};
  }
  if (varname == "p_pion") {
    return plot_data{
        .bins = 30,
        .xmin = 0.,
        .xmax = 10.,
        .name = "#it{p}_{#pi} (GeV/#it{c})",
        .ytitle =
            "d#sigma/d#it{p}_{#pi} (10^{#minus 38} cm^{2}/GeV/#it{c}/nucleon)"};
  }
  if (varname == "theta_pion") {
    return plot_data{
        .bins = 30,
        .xmin = 0.,
        .xmax = 180.,
        .name = "#theta_{#pi} (degree)",
        .ytitle =
            "d#sigma/d#theta_{#pi} (10^{#minus 38} cm^{2}/degree/nucleon)"};
  }
  throw std::runtime_error("Unknown variable name");
}

template <typename T, size_t N>
ROOT::RDF::RResultPtr<TH1> make_plots(T &&df_in, std::array<double, N> bins,
                                      std::string variable,
                                      std::string postfix = "") {
  auto plot_data = get_info(variable);
  auto plot_name = postfix.empty() ? variable : (variable + "_" + postfix);
  return df_in.template Histo1D<double>(
      {plot_name.c_str(), variable.c_str(), bins.size() - 1, bins.data()},
      variable, "weight");
}

template <typename T>
ROOT::RDF::RResultPtr<TH1> make_plots(T &&df_in, std::string variable,
                                      std::string postfix = "") {
  auto plot_data = get_info(variable);
  auto plot_name = postfix.empty() ? variable : (variable + "_" + postfix);
  return df_in.template Histo1D<double>({plot_name.c_str(), variable.c_str(),
                                         plot_data.bins, plot_data.xmin,
                                         plot_data.xmax},
                                        variable, "weight");
}

template <typename T, size_t N = 0>
auto plot_channels_0pi(T &&df_in, std::string variable,
                       std::array<double, N> bins = std::array<double, 0>{}) {
  return list_channel_channeldef | std::views::enumerate |
         std::views::transform([&](auto &&name_id_count) {
           auto &[count, name_id] = name_id_count;
           auto &[name, id] = name_id;
           if constexpr (N != 0)
             return std::make_tuple(name,
                                    make_plots(df_in.Filter(id, {"channel"}),
                                               bins, variable, name + "_0pi"),
                                    count * 2);
           else
             return std::make_tuple(name,
                                    make_plots(df_in.Filter(id, {"channel"}),
                                               variable, name + "_0pi"),
                                    count * 2);
         }) |
         std::ranges::to<std::vector>();
}

template <typename T, size_t N = 0>
auto plot_channels_pi0(T &&df_in, std::string variable,
                       std::array<double, N> bins = std::array<double, 0>{},
                       const std::string &add_name = {}) {
  return std::views::cartesian_product(list_channel_channeldef |
                                           std::views::enumerate,
                                       list_pion_cut | std::views::enumerate) |
         std::views::transform([&](auto &&tup) {
           auto &&[count1_name_id, count2_pion_cut] = tup;

           auto &&[count1, name_id] = count1_name_id;
           auto &&[count2, pion_cut] = count2_pion_cut;

           auto &&[name, id] = name_id;
           auto &&[pion_name, cut] = pion_cut;
           if constexpr (N != 0)
             return std::make_tuple(
                 name + add_name + " " + pion_pretty_name(pion_name),
                 make_plots(df_in.Filter(id, {"channel"}).Filter(cut, {"npi0"}),
                            bins, variable, name + "_" + pion_name),
                 count1 * 2 + (count2));
           else
             return std::make_tuple(
                 name + add_name + " " + pion_pretty_name(pion_name),
                 make_plots(df_in.Filter(id, {"channel"}).Filter(cut, {"npi0"}),
                            variable, name + "_" + pion_name),
                 count1 * 2 + (count2));
         }) |
         std::ranges::to<std::vector>();
}

// constexpr std::array<int, 12> col{1014, 1003, 1002, 1009, 1010, 1007,
//                                   1012, 1005, 1011, 1008, 1014, 1017};
// auto col = GetColorArray(8);

int main(int argc, char *argv[]) {
  namespace po = boost::program_options;
  IniColorCBminerva2pibg();
  po::options_description desc("Options");
  desc.add_options()("input-files", po::value<std::vector<std::string>>(),
                     "Input files") //
      ("add-label-var", po::value<char>()->default_value(0),
       "tag to be added to vars") //
      ("add-label-tki", po::value<char>()->default_value(0),
       "tag to be added to TKI vars")
      //                (
      // "run-tag", po::value<std::string>()->default_value(""),
      // "Run tag")
      ("text-0pi", po::value<std::string>()->default_value(""),
       "Additional text") //
      ("text-pi0", po::value<std::string>()->default_value(""),
       "Additional text") //
      ("b", po::value<double>(),
       "binding energy_to_use (MeV)")("help", "produce help message");
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
  // auto additional_text = vm["add-text"].as<std::string>();
  auto text_0pi = vm["text-0pi"].as<std::string>();
  auto text_pi0 = vm["text-pi0"].as<std::string>();
  auto tag_tki_base = vm["add-label-tki"].as<char>();
  auto tag_var_base = vm["add-label-var"].as<char>();
  std::optional<double> binding_energy_to_use;
  if (vm.count("b")) {
    binding_energy_to_use = vm["b"].as<double>() / 1E3;
    std::cout << "Use binding energy " << binding_energy_to_use.value() * 1E3
              << " MeV" << std::endl;
  }
  std::cout << "Processing " << n_runs << " files" << std::endl;

  TH1::AddDirectory(false);
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame input("out_tree", files);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(input);
#endif
  auto d = TrackerPrepare(input);
  auto count = input.Count();
  auto channel_merged = d.Aggregate(
      [](std::set<int> &c_set, int channel) { c_set.insert(channel); },
      [](std::set<int> c_set1, const std::set<int> &c_set2) {
        c_set1.insert(c_set2.begin(), c_set2.end());
        return c_set1;
      },
      "channel", std::set<int>{});

  auto doTKI_pi0 = [=](ROOT::RDF::RNode d) {
    return CommonVariableDefinePI0(DoTKICut_MINERvA_pi0(std::move(d)), binding_energy_to_use)
        .Define("npi0",
                [](ROOT::RVec<TLorentzVector> &pions) { return pions.size(); },
                {"good_pion"})
        .Define("dalphat", [](TKIVars &vars) { return vars.dalphat; },
                {"TKIVars"})
        .Define("IApN", [](TKIVars &vars) { return vars.IApN; }, {"TKIVars"});
  };

  auto doTKI_0pi = [=](ROOT::RDF::RNode d) {
    return CommonVariableDefine0PI(DoTKICut_MINERvA0PI(std::move(d)), binding_energy_to_use)
        .Define("dalphat", [](TKIVars &vars) { return vars.dalphat; },
                {"TKIVars"})
        .Define("IApN", [](TKIVars &vars) { return vars.IApN; }, {"TKIVars"});
  };

  auto rdf_pi0_after_cut = doTKI_pi0(d);

  auto rdf_0pi_after_cut = doTKI_0pi(d);

  std::list<ROOT::RDF::RResultPtr<TH1>> plots{};

  // plots for full data
  auto pred_all_IApN_pi0 = plots.emplace_back(make_plots(
      rdf_pi0_after_cut, MINERvA_TKI::pi0::IApN::get_binning(), "IApN", "pi0"));
  auto pred_all_IApN_0pi = plots.emplace_back(
      make_plots(rdf_0pi_after_cut, MINERvA_TKI::ZeroPi::IApN::get_binning(),
                 "IApN", "0pi"));
  auto pred_all_dalphat_pi0 = plots.emplace_back(
      make_plots(rdf_pi0_after_cut, MINERvA_TKI::pi0::dalphat::get_binning(),
                 "dalphat", "pi0"));
  auto pred_all_dalphat_0pi = plots.emplace_back(
      make_plots(rdf_0pi_after_cut, MINERvA_TKI::ZeroPi::dalphat::get_binning(),
                 "dalphat", "0pi"));

  // per channel stacked plots
  // pn 0pi
  auto plots_0pi_IApN = plot_channels_0pi(
      rdf_0pi_after_cut, "IApN", MINERvA_TKI::ZeroPi::IApN::get_binning());

  // dat 0pi
  auto plots_0pi_dalphat =
      plot_channels_0pi(rdf_0pi_after_cut, "dalphat",
                        MINERvA_TKI::ZeroPi::dalphat::get_binning());

  // pn pi0
  auto plots_pi0_IApN = plot_channels_pi0(
      rdf_pi0_after_cut, "IApN", MINERvA_TKI::pi0::IApN::get_binning());

  // dat pi0
  auto plots_pi0_dalphat = plot_channels_pi0(
      rdf_pi0_after_cut, "dalphat", MINERvA_TKI::pi0::dalphat::get_binning());

  // auto plots_pi0_W = get_channel_list_pi0(rdf_pi0_after_cut, "W");
  std::list<std::string> vars_0pi{"W",        "Q2",  "xBj",    "p_mu",
                                  "theta_mu", "p_p", "theta_p"};
  std::list<std::string> vars_pi0{"W",       "Q2",       "xBj",
                                  "p_mu",    "theta_mu", "p_p",
                                  "theta_p", "p_pion",   "theta_pion"};
  auto var_plots_0pi =
      vars_0pi | std::views::transform([&](std::string name) {
        return plot_channels_0pi(rdf_0pi_after_cut, std::move(name));
      }) |
      std::ranges::to<std::vector>();
  auto &plot_list_W = var_plots_0pi[0];
  // re-sort the list following:
  // 2p2h, QE, RES, DIS, 1pi, 2pi
  std::ranges::sort(plot_list_W, [](auto &t1, auto &t2) -> bool {
    auto name1 = std::get<0>(t1);
    auto name2 = std::get<0>(t2);
    auto to_var = [](const std::string &x) -> short {
      if (x.contains("2p2h"))
        return 0;
      if (x.contains("QE"))
        return 1;
      if (x.contains("RES"))
        return 2;
      if (x.contains("DIS"))
        return 3;
      if (x.contains("1#pi"))
        return 4;
      if (x.contains("2#pi"))
        return 5;
      return -1;
    };
    return to_var(name1) < to_var(name2);
  });

  auto var_plots_pi0 =
      vars_pi0 | std::views::transform([&](std::string name) {
        return plot_channels_pi0(rdf_pi0_after_cut, std::move(name));
      }) |
      std::ranges::to<std::vector>();
  auto add_to_plots_vars = [&](auto &&list) {
    for (auto &&list_entry : list) {
      for (auto &&plot : list_entry | std::views::values) {
        plots.emplace_back(plot);
      }
    }
  };
  add_to_plots_vars(var_plots_0pi);
  add_to_plots_vars(var_plots_pi0);

  auto add_to_plots = [&](auto &&list) {
    for (auto &&plot : list | std::views::values) {
      plots.emplace_back(plot);
    }
  };
  add_to_plots(plots_0pi_IApN);
  add_to_plots(plots_0pi_dalphat);
  add_to_plots(plots_pi0_IApN);
  add_to_plots(plots_pi0_dalphat);

  auto df_pi0_report = rdf_pi0_after_cut.Report();
  auto df_0pi_report = rdf_0pi_after_cut.Report();

  auto hist2d =
      rdf_0pi_after_cut
          .Define(
              "InitNucleonMass",
              [](const TLorentzVector &InitNucleon) { return InitNucleon.M(); },
              {"InitNucleon"})
          .Histo2D({"W_vs_InitNucleonMass", "W_vs_InitNucleonMass;W;M", 60, 0.8,
                    1.2, 60, 0.8, 1.2},
                   "W", "InitNucleonMass", "weight");
  auto hist2d_another =
      rdf_0pi_after_cut
          .Define("deltam2",
                  [](const TLorentzVector &InitNeutrino,
                     const TLorentzVector &PrimaryLepton,
                     const TLorentzVector &InitNucleon) {
                    auto had_sys = InitNucleon + InitNeutrino - PrimaryLepton;
                    auto W = had_sys.M();
                    auto M = InitNucleon.M();
                    return W - M;
                  },
                  {"InitNeutrino", "PrimaryLepton", "InitNucleon"})
          .Define("pdotq",
                  [](const TLorentzVector &InitNeutrino,
                     const TLorentzVector &PrimaryLepton,
                     const TLorentzVector &InitNucleon) {
                    auto q = InitNeutrino - PrimaryLepton;
                    const auto &p = InitNucleon;
                    return p.Dot(q);
                  },
                  {"InitNeutrino", "PrimaryLepton", "InitNucleon"})
          .Histo2D({"deltam2_vs_pdotq", "deltam2_vs_pdotq;deltam2;pdotq", 60,
                    -0.5, 0.5, 60, 0, 5},
                   "deltam2", "pdotq", "weight");
  auto df_high_dat = rdf_pi0_after_cut.Filter(
      [](double dat) { return dat > 90; }, {"dalphat"});
  auto hist_q0_q2 = rdf_pi0_after_cut.Histo2D(
      {"Q2_vs_q0", "Q2_vs_q0;Q^{2};q_{0}", 50, 0, 5, 50, 0, 5}, "Q2", "q0",
      "weight");
  auto hist_q0_q2_high_dat = df_high_dat.Histo2D(
      {"Q2_vs_q0_high_dat", "Q2_vs_q0;Q^{2};q_{0}", 50, 0, 5, 50, 0, 5}, "Q2",
      "q0", "weight");
  auto hist_W = rdf_pi0_after_cut.Histo1D(
      {"W_all", "W;W (GeV);d#sigma/dW (10^{#minus 38} cm^{2}/GeV/nucleon)", 60,
       0.7, 4.0},
      "W", "weight");
  auto hist_W_high_dat = df_high_dat.Histo1D(
      {"W_high_dat",
       "W_high_dat;W (GeV);d#sigma/dW (10^{#minus 38} cm^{2}/GeV/nucleon)", 60,
       0.7, 4.0},
      "W", "weight");

  auto plot_W_pi0 =
      plot_channels_pi0(rdf_pi0_after_cut, "W", std::array<double, 0>{}, "all");
  auto plot_W_pi0_high_dat =
      plot_channels_pi0(df_high_dat, "W", std::array<double, 0>{}, "high_dat");
  // merge 2 lists
  plot_W_pi0.insert(plot_W_pi0.end(), plot_W_pi0_high_dat.begin(),
                    plot_W_pi0_high_dat.end());

  rdf_pi0_after_cut.Filter("channel>= 2 && channel <= 31")
      .Snapshot("check", "check.root",
                {"dalphat", "W", "Q2", "q0", "StdHepN", "StdHepPdg",
                 "StdHepStatus", "weight", "channel", "InitNucleon",
                 "InitNeutrino", "PrimaryLepton", "full_hadron"});

  auto file = std::make_unique<TFile>("dtl_all.root", "RECREATE");

  for (auto &&plot : plot_W_pi0 | std::views::elements<1>) {
    plot->Scale((12. / 13.) / n_runs / 10, "width");
    file->Add(plot.GetPtr());
  }

  hist2d->Scale((12. / 13.) / n_runs / 10, "width");
  file->Add(hist2d.GetPtr());
  hist2d_another->Scale((12. / 13.) / n_runs / 10, "width");
  file->Add(hist2d_another.GetPtr());
  hist_q0_q2->Scale((12. / 13.) / n_runs / 10, "width");
  file->Add(hist_q0_q2.GetPtr());

  hist_W->Scale((12. / 13.) / n_runs / 10, "width");
  file->Add(hist_W.GetPtr());

  hist_q0_q2_high_dat->Scale((12. / 13.) / n_runs / 10, "width");
  file->Add(hist_q0_q2_high_dat.GetPtr());

  hist_W_high_dat->Scale((12. / 13.) / n_runs / 10, "width");
  file->Add(hist_W_high_dat.GetPtr());

  for (auto &&plot : plots) {
    plot->Scale((12. / 13.) / n_runs / 10, "width");
    file->Add(plot.GetPtr());
  }

  std::ranges::for_each(var_plots_0pi[0] |
                            std::views::filter([](auto &tup) -> bool {
                              return std::get<0>(tup) == "QE";
                            }),
                        [](auto &hist) {
                          auto &[name, hist_ptr, count] = hist;
                          name = name + " (#times 1/3)";
                          hist_ptr->Scale(1 / 3.);
                        });

  constexpr double threshold_frac_0pi = .5e-2;
  constexpr double threshold_frac_pi0 = .5e-2;

  legend_conf conf_0pi{
      // .force_include = {"QE", "RES", "DIS"}
  };

  legend_conf conf_pi0{
      .force_include = {"QE", "RES", "DIS"},
      .force_exclude = {"2p2h"},
  };

  auto xsecint_0pi = pred_all_dalphat_0pi->Integral("WIDTH");
  auto xsecint_pi0 = pred_all_dalphat_pi0->Integral("WIDTH");
  auto [plots_0pi_IApN_stack, plots_0pi_IApN_leg] = build_stack_from_list(
      plots_0pi_IApN, pred_all_IApN_0pi->Integral("WIDTH") * threshold_frac_0pi,
      conf_0pi);
  auto [plots_0pi_dalphat_stack, plots_0pi_dalphat_leg] = build_stack_from_list(
      plots_0pi_dalphat,
      pred_all_dalphat_0pi->Integral("WIDTH") * threshold_frac_0pi, conf_0pi);

  auto [plots_pi0_IApN_stack, plots_pi0_IApN_leg] = build_stack_from_list(
      plots_pi0_IApN, pred_all_IApN_pi0->Integral("WIDTH") * threshold_frac_pi0,
      conf_pi0);
  auto [plots_pi0_dalphat_stack, plots_pi0_dalphat_leg] = build_stack_from_list(
      plots_pi0_dalphat,
      pred_all_dalphat_pi0->Integral("WIDTH") * threshold_frac_pi0, conf_pi0);

  auto stacked_vars_0pi =
      var_plots_0pi | std::views::transform([&](auto &&list) {
        return build_stack_from_list(list, xsecint_0pi * threshold_frac_0pi,
                                     conf_0pi);
      }) |
      std::ranges::to<std::vector>();
  auto stacked_vars_pi0 =
      var_plots_pi0 | std::views::transform([&](auto &&list) {
        return build_stack_from_list(list, xsecint_pi0 * threshold_frac_pi0,
                                     conf_pi0);
      }) |
      std::ranges::to<std::vector>();

  auto exp_data_hist_IApN_pi0 = MINERvA_TKI::pi0::IApN::get_hist();
  auto exp_data_hist_IApN_0pi = MINERvA_TKI::ZeroPi::IApN::get_hist();
  auto exp_data_hist_dalphat_pi0 = MINERvA_TKI::pi0::dalphat::get_hist();
  auto exp_data_hist_dalphat_0pi = MINERvA_TKI::ZeroPi::dalphat::get_hist();

  file->Add(&exp_data_hist_IApN_pi0);
  file->Add(&exp_data_hist_IApN_0pi);
  file->Add(&exp_data_hist_dalphat_pi0);
  file->Add(&exp_data_hist_dalphat_0pi);

  auto chi2_IApN_pi0 =
      MINERvA_TKI::pi0::IApN::do_chi2(pred_all_IApN_pi0.GetPtr());
  auto chi2_IApN_0pi =
      MINERvA_TKI::ZeroPi::IApN::do_chi2(pred_all_IApN_0pi.GetPtr());
  auto chi2_dalphat_pi0 =
      MINERvA_TKI::pi0::dalphat::do_chi2(pred_all_dalphat_pi0.GetPtr());
  auto chi2_dalphat_0pi =
      MINERvA_TKI::ZeroPi::dalphat::do_chi2(pred_all_dalphat_0pi.GetPtr());

  auto &&[shape_only_chi2_IApN_pi0, shape_only_prediction_IApN_pi0,
          shape_only_data_IApN_pi0, norm_IApN_pi0] =
      do_chi2_shape_only(MINERvA_TKI::pi0::IApN::get_cov(),
                         MINERvA_TKI::pi0::IApN::get_hist(),
                         *(dynamic_cast<TH1D *>(pred_all_IApN_pi0.GetPtr())));
  file->Add(&shape_only_data_IApN_pi0);
  file->Add(&shape_only_prediction_IApN_pi0);
  auto stack_shape_IApN_pi0 =
      scale_stack(plots_pi0_IApN_stack, 1 / norm_IApN_pi0);

  auto &&[shape_only_chi2_IApN_0pi, shape_only_prediction_IApN_0pi,
          shape_only_data_IApN_0pi, norm_IApN_0pi] =
      do_chi2_shape_only(MINERvA_TKI::ZeroPi::IApN::get_cov(),
                         MINERvA_TKI::ZeroPi::IApN::get_hist(),
                         *(dynamic_cast<TH1D *>(pred_all_IApN_0pi.GetPtr())));
  file->Add(&shape_only_data_IApN_0pi);
  file->Add(&shape_only_prediction_IApN_0pi);
  auto stack_shape_IApN_0pi =
      scale_stack(plots_0pi_IApN_stack, 1 / norm_IApN_0pi);

  auto &&[shape_only_chi2_dalphat_pi0, shape_only_prediction_dalphat_pi0,
          shape_only_data_dalphat_pi0, norm_dalphat_pi0] =
      do_chi2_shape_only(
          MINERvA_TKI::pi0::dalphat::get_cov(),
          MINERvA_TKI::pi0::dalphat::get_hist(),
          *(dynamic_cast<TH1D *>(pred_all_dalphat_pi0.GetPtr())));
  file->Add(&shape_only_data_dalphat_pi0);
  file->Add(&shape_only_prediction_dalphat_pi0);
  auto stack_shape_dalphat_pi0 =
      scale_stack(plots_pi0_dalphat_stack, 1 / norm_dalphat_pi0);

  auto &&[shape_only_chi2_dalphat_0pi, shape_only_prediction_dalphat_0pi,
          shape_only_data_dalphat_0pi, norm_dalphat_0pi] =
      do_chi2_shape_only(
          MINERvA_TKI::ZeroPi::dalphat::get_cov(),
          MINERvA_TKI::ZeroPi::dalphat::get_hist(),
          *(dynamic_cast<TH1D *>(pred_all_dalphat_0pi.GetPtr())));
  file->Add(&shape_only_data_dalphat_0pi);
  file->Add(&shape_only_prediction_dalphat_0pi);
  auto stack_shape_dalphat_0pi =
      scale_stack(plots_0pi_dalphat_stack, 1 / norm_dalphat_0pi);

  std::cout << "chi2_IApN_pi0 " << chi2_IApN_pi0 << std::endl;
  std::cout << "chi2_IApN_0pi " << chi2_IApN_0pi << std::endl;
  std::cout << "chi2_dalphat_pi0 " << chi2_dalphat_pi0 << std::endl;
  std::cout << "chi2_dalphat_0pi " << chi2_dalphat_0pi << std::endl;

  std::cout << "shape_only_chi2_IApN_pi0 " << shape_only_chi2_IApN_pi0
            << std::endl;
  std::cout << "shape_only_chi2_IApN_0pi " << shape_only_chi2_IApN_0pi
            << std::endl;
  std::cout << "shape_only_chi2_dalphat_pi0 " << shape_only_chi2_dalphat_pi0
            << std::endl;
  std::cout << "shape_only_chi2_dalphat_0pi " << shape_only_chi2_dalphat_0pi
            << std::endl;

  nlohmann::json json_out;
  json_out["chi2_IApN_pi0"] = chi2_IApN_pi0;
  json_out["chi2_IApN_0pi"] = chi2_IApN_0pi;
  json_out["chi2_dalphat_pi0"] = chi2_dalphat_pi0;
  json_out["chi2_dalphat_0pi"] = chi2_dalphat_0pi;
  json_out["shape_only_chi2_IApN_pi0"] = shape_only_chi2_IApN_pi0;
  json_out["shape_only_chi2_IApN_0pi"] = shape_only_chi2_IApN_0pi;
  json_out["shape_only_chi2_dalphat_pi0"] = shape_only_chi2_dalphat_pi0;
  json_out["shape_only_chi2_dalphat_0pi"] = shape_only_chi2_dalphat_0pi;
  json_out["xsecint_0pi"] = xsecint_0pi;
  json_out["xsecint_pi0"] = xsecint_pi0;
  json_out["xsecint_pn_0pi"] = pred_all_IApN_pi0->Integral("WIDTH");
  json_out["xsecint_dat_0pi"] = pred_all_dalphat_0pi->Integral("WIDTH");
  json_out["xsecint_pn_pi0"] = pred_all_IApN_pi0->Integral("WIDTH");
  json_out["xsecint_dat_pi0"] = pred_all_dalphat_pi0->Integral("WIDTH");
  std::fstream json_file("chi2.json", std::ios::out);
  json_file << json_out.dump(2) << std::endl;

  // if (!channel_merged->count(37)) {
  //   std::cout << "No 2pi BG channel found, marking legend as w/o 2#piBG"
  //             << std::endl;
  //   runtag += "w/o 2#piBG ";
  // }

  auto form_legend = [&](TH1 *hist, double chi2) -> std::string {
    // return std::string{"GiBUU"} + " #chi^{2}/NDF = " + std::to_string(chi2);
    std::stringstream ss;
    ss << "#chi^{2}/NDF = " << round(chi2) << "/" << hist->GetNbinsX();
    return ss.str();
  };

  auto build_add_text = [&](char tag, std::string additional_text) {
    std::unique_ptr<TLatex> latex;
    if (!additional_text.empty()) {
      if (tag)
        latex = std::make_unique<TLatex>(
            0.65, 0.5, std::format("({}) {}", tag, additional_text).c_str());
      else
        latex = std::make_unique<TLatex>(0.65, 0.5, additional_text.c_str());
      latex->SetNDC();
    }
    return latex;
  };

  // first entry always data, then anything else
  do_plot({&exp_data_hist_IApN_pi0, &plots_pi0_IApN_stack,
           &plots_pi0_dalphat_leg,
           build_add_text(tag_tki_base + 2, text_pi0).get()},
          "IApN_pi0", get_info("IApN").ytitle, get_info("IApN").name,
          {0.65, 0.45, 0.95, 0.93}, 0.8,
          form_legend(&exp_data_hist_IApN_pi0, chi2_IApN_pi0), "HISTC", 0,
          {.top = 0.02, .bottom = 0.125});

  do_plot({&shape_only_data_IApN_pi0, &stack_shape_IApN_pi0,
           &plots_pi0_IApN_leg,
           build_add_text(tag_tki_base + 2, text_pi0).get()},
          "IApN_pi0_shape", "shape", get_info("IApN").name,
          {0.65, 0.45, 0.95, 0.93}, 0.8,
          form_legend(&exp_data_hist_IApN_pi0, shape_only_chi2_IApN_pi0),
          "HISTC", 0, {.top = 0.02, .bottom = 0.125});

  do_plot({&exp_data_hist_IApN_0pi, &plots_0pi_IApN_stack, &plots_0pi_IApN_leg,
           build_add_text(tag_tki_base + 2, text_0pi).get()},
          "IApN_0pi", get_info("IApN").ytitle, get_info("IApN").name,
          {0.65, 0.55, 0.95, 0.93}, 0.8,
          form_legend(&exp_data_hist_IApN_0pi, chi2_IApN_0pi), "HISTC", 0,
          {.top = 0.02, .bottom = 0.125});

  do_plot({&shape_only_data_IApN_0pi, &stack_shape_IApN_0pi,
           &plots_0pi_IApN_leg,
           build_add_text(tag_tki_base + 2, text_pi0).get()},
          "IApN_0pi_shape", "shape", get_info("IApN").name,
          {0.65, 0.55, 0.95, 0.93}, 0.8,
          form_legend(&exp_data_hist_IApN_pi0, shape_only_chi2_IApN_pi0),
          "HISTC", 0, {.top = 0.02, .bottom = 0.125});

  do_plot({&exp_data_hist_dalphat_pi0, &plots_pi0_dalphat_stack,
           &plots_pi0_dalphat_leg,
           build_add_text(tag_tki_base, text_pi0).get()},
          "dalphat_pi0", get_info("dalphat").ytitle, get_info("dalphat").name,
          {0.15, 0.45, 0.4, 0.9}, 0.,
          form_legend(&exp_data_hist_dalphat_pi0, chi2_dalphat_pi0), "HISTC", 0,
          {.top = 0.06, .bottom = 0.11});

  do_plot({&shape_only_data_dalphat_pi0, &stack_shape_dalphat_pi0,
           &plots_pi0_dalphat_leg,
           build_add_text(tag_tki_base, text_pi0).get()},
          "dalphat_pi0_shape", "shape", get_info("dalphat").name,
          {0.2, 0.45, 0.5, 0.9}, 0.,
          form_legend(&exp_data_hist_dalphat_pi0, shape_only_chi2_dalphat_pi0),
          "HISTC", 0, {.top = 0.06, .bottom = 0.11});

  do_plot({&exp_data_hist_dalphat_0pi, &plots_0pi_dalphat_stack,
           &plots_0pi_dalphat_leg,
           build_add_text(tag_tki_base, text_0pi).get()},
          "dalphat_0pi", get_info("dalphat").ytitle, get_info("dalphat").name,
          {0.2, 0.6, 0.5, 0.9}, 0.,
          form_legend(&exp_data_hist_dalphat_0pi, chi2_dalphat_0pi), "HISTC", 0,
          {.top = 0.06, .bottom = 0.11});

  do_plot({&shape_only_data_dalphat_0pi, &stack_shape_dalphat_0pi,
           &plots_0pi_dalphat_leg,
           build_add_text(tag_tki_base, text_0pi).get()},
          "dalphat_0pi_shape", "shape", get_info("dalphat").name,
          {0.2, 0.5, 0.5, 0.9}, 0.,
          form_legend(&exp_data_hist_dalphat_0pi, shape_only_chi2_dalphat_0pi),
          "HISTC", 0, {.top = 0.06, .bottom = 0.11});

  for (auto &&[name, list] : std::views::zip(vars_0pi, stacked_vars_0pi)) {
    auto &&[stack, leg] = list;
    auto plot_ent = get_info(name);
    auto xmax = plot_ent.xmax_0pi == 0. ? plot_ent.xmax : plot_ent.xmax_0pi;
    do_plot({&stack, &leg, build_add_text(tag_var_base, text_0pi).get()},
            name + "_0pi", plot_ent.ytitle, plot_ent.name,
            {0.65, 0.45, 0.95, 0.93}, xmax, "MINERvA 0#pi", "HIST",
            plot_ent.ymax_0pi,
            {.top = plot_ent.ymax_0pi < 0.1 ? 0.06 : 0.02, .bottom = 0.11});
  }
  for (auto &&[name, list] : std::views::zip(vars_pi0, stacked_vars_pi0)) {
    auto &&[stack, leg] = list;
    auto plot_ent = get_info(name);
    auto xmax = plot_ent.xmax_pi0 == 0. ? plot_ent.xmax : plot_ent.xmax_pi0;
    do_plot({&stack, &leg, build_add_text(tag_var_base, text_pi0).get()},
            name + "_pi0", plot_ent.ytitle, plot_ent.name,
            {0.65, 0.45, 0.95, 0.93}, xmax, "MINERvA#kern[0.25]{#pi^{0}}",
            "HIST", plot_ent.ymax_pi0,
            {.top = plot_ent.ymax_pi0 < 0.1 ? 0.06 : 0.02, .bottom = 0.11});
  }
  /////////////////////
  auto plots_pi0_IApN_list_raw =
      plots_pi0_IApN | std::views::transform([](auto &&tup) -> plot_ptr_t {
        auto ptr = std::get<1>(tup);
        ptr->SetFillStyle(0);
        return ptr;
      }) |
      std::ranges::to<std::vector>();
  plots_pi0_IApN_list_raw.emplace_back(&plots_pi0_dalphat_leg);
  auto tex_pi0 = build_add_text(tag_tki_base, text_pi0);
  auto tex_0pi = build_add_text(tag_tki_base, text_0pi);
  plots_pi0_IApN_list_raw.emplace_back(tex_pi0.get());
  do_plot(plots_pi0_IApN_list_raw, "IApN_pi0_nostack", get_info("IApN").ytitle,
          get_info("IApN").name, {0.65, 0.45, 0.95, 0.93}, 0.8,
          form_legend(&exp_data_hist_IApN_pi0, chi2_IApN_pi0), "HISTC", 0.2,
          {.top = 0.02, .bottom = 0.125});

  auto plots_0pi_IApN_list_raw =
      plots_0pi_IApN | std::views::transform([](auto &&tup) -> plot_ptr_t {
        auto ptr = std::get<1>(tup);
        ptr->SetFillStyle(0);
        return ptr;
      }) |
      std::ranges::to<std::vector>();
  plots_0pi_IApN_list_raw.emplace_back(&plots_0pi_IApN_leg);
  plots_0pi_IApN_list_raw.emplace_back(tex_0pi.get());
  do_plot(plots_0pi_IApN_list_raw, "IApN_0pi_nostack", get_info("IApN").ytitle,
          get_info("IApN").name, {0.65, 0.55, 0.95, 0.93}, 0.8,
          form_legend(&exp_data_hist_IApN_0pi, chi2_IApN_0pi), "HISTC", 0.4,
          {.top = 0.02, .bottom = 0.125});
  /////////////////////

  file->Write();
  file->Close();

  for (auto channel : *channel_merged) {
    if (std::ranges::none_of(list_channel_channeldef | std::views::values,
                             [channel](auto &&id) { return id(channel); })) {
      std::cout << "WARNING!!! Unknown channel " << channel << std::endl;
    }
  }
  std::cout << std::endl;

  std::cout << "0pi Cut Summary:" << std::endl;
  df_0pi_report->Print();
  std::cout << "pi0 Cut Summary:" << std::endl;
  df_pi0_report->Print();
  return 0;
}
