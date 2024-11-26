#include "EvtTracker2event.h"
#include "MINERvA_tki_cut.h"
#include "T2K_tki_cut.h"
#include "dochi2.h"
#include "plottools.hxx"
#include "tkigeneral.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <THStack.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMatrixDfwd.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>
#include <boost/program_options.hpp>

#include <cmath>
#include <cstdlib>
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
        {"2#pi BG", [](int id) { return id == 37; }},
        {"2p2h", [](int id) { return id == 35; }},
        {"QE", [](int id) { return id == 1; }},
        {"RES", [](int id) { return id >= 2 && id <= 33; }},
        {"DIS", [](int id) { return id == 34; }},
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

plot_data get_info(const std::string &varname) {
  if (varname == "dalphat") {
    return plot_data{30, 0., 180., "#delta#it{#alpha}_{T} (degree)",
                     "d#sigma/d#delta#it{#alpha}_{T} (10^{#minus 38} "
                     "cm^{2}/degree/nucleon)"};
  }
  if (varname == "IApN") {
    return plot_data{30, 0., 0.8, "#it{p}_{N} (GeV/#it{c})",
                     "d#sigma/d#it{p}_{N} (10^{#minus 38} "
                     "cm^{2}/GeV/#it{c}/nucleon)"};
  }
  if (varname == "dpL") {
    return plot_data{30, -1., 1., "#delta p_{L}"};
  }
  if (varname == "p_cos_theta") {
    return plot_data{30, -1., 1., "cos#theta_{p}"};
  }
  if (varname == "p_theta") {
    return plot_data{30, 0., 180., "#theta_{p}"};
  }
  if (varname == "dtl") {
    return plot_data{30, -1., 1., "cos #delta#alpha_{L}"};
  }
  if (varname == "Q2") {
    return plot_data{50,
                     0.,
                     5.,
                     "#it{Q}^{2} (GeV^{#kern[0.3]{2}})",
                     "d#sigma/d#it{Q}^{2} (10^{#minus 38} "
                     "cm^{2}/GeV^{#kern[0.3]{2}}/nucleon)",
                     0.3,
                     0.15};
  }
  if (varname == "xBj") {
    return plot_data{50,
                     0.,
                     2.0,
                     "#it{x}_{Bj}",
                     "d#sigma/d#it{x}_{Bj} (10^{#minus 38} cm^{2}/nucleon)",
                     .85,
                     .25};
  }
  if (varname == "W") {
    return plot_data{60,
                     0.7,
                     1.8,
                     "#it{W} (GeV)",
                     "d#sigma/d#it{W} (10^{#minus 38} cm^{2}/GeV/nucleon)",
                     4.5,
                     .2};
  }
  if (varname == "p_mu") {
    return plot_data{
        30, 0., 12., "#it{p}_{#mu} (GeV/#it{c})",
        "d#sigma/d#it{p}_{#mu} (10^{#minus 38} cm^{2}/GeV/#it{c}/nucleon)"};
  }
  if (varname == "theta_mu") {
    return plot_data{
        30, 0., 60., "#theta_{#mu} (degree)",
        "d#sigma/d#theta_{#mu} (10^{#minus 38} cm^{2}/degree/nucleon)"};
  }
  if (varname == "p_p") {
    return plot_data{
        30, 0., 20., "#it{p}_{p} (GeV/#it{c})",
        "d#sigma/d#it{p}_{p} (10^{#minus 38} cm^{2}/GeV/#it{c}/nucleon)"};
  }
  if (varname == "theta_p") {
    return plot_data{
        30, 0., 180., "#theta_{p} (degree)",
        "d#sigma/d#theta_{p} (10^{#minus 38} cm^{2}/degree/nucleon)"};
  }
  throw std::runtime_error("Unknown variable name");
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
                       std::array<double, N> bins = std::array<double, 0>{}) {
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
                 name + " " + pion_pretty_name(pion_name),
                 make_plots(df_in.Filter(id, {"channel"}).Filter(cut, {"npi0"}),
                            bins, variable, name + "_" + pion_name),
                 count1 * 2 + count2);
           else
             return std::make_tuple(
                 name + " " + pion_pretty_name(pion_name),
                 make_plots(df_in.Filter(id, {"channel"}).Filter(cut, {"npi0"}),
                            variable, name + "_" + pion_name),
                 count1 * 2 + count2);
         }) |
         std::ranges::to<std::vector>();
}

int main(int argc, char *argv[]) {
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options() //
      ("input-files", po::value<std::vector<std::string>>(),
       "Input files") //
      // ("add-label", po::value<std::string>()->default_value(""),
      //  "Additional label")
      //                (
      // "run-tag", po::value<std::string>()->default_value(""),
      // "Run tag")
      ("add-text", po::value<std::string>()->default_value(""),
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

  auto doTKI_0pi = [](ROOT::RDF::RNode d) {
    return CommonVariableDefine0PI(DoTKICut_T2K_STK(std::move(d)))
        .Define("dalphat", [](TKIVars &vars) { return vars.dalphat; },
                {"TKIVars"})
        // .Define("dpL", [](TKIVars &vars) { return vars.dpL; }, {"TKIVars"})
        .Define("IApN", [](TKIVars &vars) { return vars.IApN; }, {"TKIVars"});
  };

  auto rdf_0pi_after_cut = doTKI_0pi(d);

  std::list<ROOT::RDF::RResultPtr<TH1>> plots{};

  auto pred_all_dalphat_0pi = plots.emplace_back(
      make_plots(rdf_0pi_after_cut, T2K_STK::get_binning(), "dalphat", "0pi"));

  // per channel stacked plots

  // dat 0pi
  auto plots_0pi_dalphat =
      plot_channels_0pi(rdf_0pi_after_cut, "dalphat", T2K_STK::get_binning());

  std::list<std::string> vars{"W",        "Q2",  "xBj",    "p_mu",
                              "theta_mu", "p_p", "theta_p"};
  auto var_plots_0pi =
      vars | std::views::transform([&](std::string name) {
        return plot_channels_0pi(rdf_0pi_after_cut, std::move(name));
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
  // add_to_plots_vars(var_plots_pi0);

  auto add_to_plots = [&](auto &&list) {
    for (auto &&plot : list | std::views::values) {
      plots.emplace_back(plot);
    }
  };
  add_to_plots(plots_0pi_dalphat);

  auto df_0pi_report = rdf_0pi_after_cut.Report();

  auto file = std::make_unique<TFile>("dtl_all.root", "RECREATE");
  for (auto &&plot : plots) {
    plot->Scale((12. / 13.) / n_runs / 10, "width");
    file->Add(plot.GetPtr());
  }

  constexpr double threshold_frac_0pi = 1.5e-2;

  auto xsecint_0pi = pred_all_dalphat_0pi->Integral("WIDTH");

  auto [plots_0pi_dalphat_stack, plots_0pi_dalphat_leg] = build_stack_from_list(
      plots_0pi_dalphat,
      pred_all_dalphat_0pi->Integral("WIDTH") * threshold_frac_0pi);

  auto stacked_vars_0pi =
      var_plots_0pi | std::views::transform([&](auto &&list) {
        return build_stack_from_list(list, xsecint_0pi * threshold_frac_0pi);
      }) |
      std::ranges::to<std::vector>();

  auto exp_data_hist_dalphat_0pi = T2K_STK::get_hist();

  file->Add(&exp_data_hist_dalphat_0pi);

  auto chi2_dalphat_0pi = T2K_STK::do_chi2(pred_all_dalphat_0pi.GetPtr());

  // std::cout << "chi2_IApN_0pi " << chi2_IApN_0pi << std::endl;
  std::cout << "chi2_dalphat_0pi " << chi2_dalphat_0pi << '\n';
  auto form_legend = [&](TH1 *hist, double chi2) -> std::string {
    // return std::string{"GiBUU"} + " #chi^{2}/NDF = " + std::to_string(chi2);
    std::stringstream ss;
    ss << "#chi^{2}/NDF = " << round(chi2) << "/" << hist->GetNbinsX();
    return ss.str();
  };

  // first entry always data, then anything else
  // do_plot({&exp_data_hist_IApN_pi0, pred_all_IApN_pi0, &plots_pi0_IApN_stack,
  //          &plots_pi0_IApN_leg},
  //         "IApN_pi0", get_info("IApN").ytitle, get_info("IApN").name,
  //         {0.5, 0.55, 0.8, 0.9}, 0.8,
  //         form_legend(&exp_data_hist_IApN_pi0, chi2_IApN_pi0));

  // do_plot({&exp_data_hist_IApN_0pi, pred_all_IApN_0pi, &plots_0pi_IApN_stack,
  //          &plots_0pi_IApN_leg},
  //         "IApN_0pi", get_info("IApN").ytitle, get_info("IApN").name,
  //         {0.5, 0.6, 0.8, 0.9}, 0.8,
  //         form_legend(&exp_data_hist_IApN_0pi, chi2_IApN_0pi));

  // do_plot({&exp_data_hist_dalphat_pi0, pred_all_dalphat_pi0,
  //          &plots_pi0_dalphat_stack, &plots_pi0_dalphat_leg},
  //         "dalphat_pi0", get_info("dalphat").ytitle,
  //         get_info("dalphat").name, {0.2, 0.5, 0.5, 0.9}, 0.,
  //         form_legend(&exp_data_hist_dalphat_pi0, chi2_dalphat_pi0));
  std::unique_ptr<TLatex> latex;
  if (!additional_text.empty()) {
    latex = std::make_unique<TLatex>(0.65, 0.5, additional_text.c_str());
    latex->SetNDC();
  }

  do_plot({&exp_data_hist_dalphat_0pi, pred_all_dalphat_0pi,
           &plots_0pi_dalphat_stack, &plots_0pi_dalphat_leg, latex.get()},
          "dalphat_0pi", get_info("dalphat").ytitle, get_info("dalphat").name,
          {0.15, 0.6, 0.5, 0.9}, 0.,
          form_legend(&exp_data_hist_dalphat_0pi, chi2_dalphat_0pi), "HISTC",
          2e-3, {.top = 0.06, .bottom = 0.12});

  for (auto &&[name, list] : std::views::zip(vars, stacked_vars_0pi)) {
    auto &&[stack, leg] = list;
    auto plot_ent = get_info(name);
    do_plot({&stack, &leg, latex.get()}, name + "_0pi", plot_ent.ytitle,
            plot_ent.name, {0.75, 0.55, 0.95, 0.9}, 0., "T2K", "HIST",
            plot_ent.ymax_0pi, {.top = 0.06, .bottom = 0.12});
  }
  // for (auto &&[name, list] : std::views::zip(vars, stacked_vars_pi0)) {
  //   auto &&[stack, leg] = list;
  //   auto plot_ent = get_info(name);
  //   do_plot({&stack, &leg}, name + "_pi0", plot_ent.ytitle, plot_ent.name,
  //           {0.55, 0.55, 0.85, 0.9}, 0.,
  //           "GiBUU 23p3 MINERvA#pi^{0} " + runtag, "HIST",
  //           plot_ent.ymax_pi0);
  // }

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
  // std::cout << "pi0 Cut Summary:" << std::endl;
  // df_pi0_report->Print();
  return 0;
}
