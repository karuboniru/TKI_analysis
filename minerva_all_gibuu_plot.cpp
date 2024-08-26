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
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMatrixDfwd.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>

#include <cmath>
#include <functional>
#include <memory>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "expdata.h"

using plot_ptr_t =
    std::variant<TH1 *, TLegend *, THStack *, ROOT::RDF::RResultPtr<TH1D>,
                 ROOT::RDF::RResultPtr<TH1>>;

// this is general channel definition
const std::vector<std::tuple<std::string, std::function<bool(int)>>>
    list_channel_channeldef{{"2#pi BG", [](int id) { return id == 37; }},
                            {"QE", [](int id) { return id == 1; }},
                            {"RES", [](int id) { return id >= 2 && id <= 33; }},
                            {"DIS", [](int id) { return id == 34; }},
                            {"2p2h", [](int id) { return id == 35; }}};

// pion count definition for pi0 channel
const std::vector<std::tuple<std::string, std::function<bool(size_t)>>>
    list_pion_cut{{"1pi0", [](size_t i) { return i == 1; }},
                  {"Mpi0", [](size_t i) { return i > 1; }}};

std::string pion_pretty_name(std::string s) {
  if (s == "1pi0") {
    return "1 #pi^{0}";
  }
  if (s == "Mpi0") {
    return "multi- #pi^{0}";
  }
  return s;
}

struct plot_data {
  int bins;
  double xmin, xmax;
  std::string name;
  std::string ytitle;
  double ymax_0pi{}, ymax_pi0{};
};

plot_data get_info(std::string varname) {
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
    return plot_data{30, -1., 1., "cos #theta_{p}"};
  }
  if (varname == "p_theta") {
    return plot_data{30, 0., 180., "#theta_{p}"};
  }
  if (varname == "dtl") {
    return plot_data{30, -1., 1., "cos #delta#alpha_{L}"};
  }
  if (varname == "Q2") {
    return plot_data{
        50,
        0.,
        5.,
        "#it{Q}^{2} (GeV^{ 2})",
        "d#sigma/d#it{Q}^{2} (10^{#minus 38} cm^{2}/GeV^{ 2}/nucleon)",
        0.3,
        0.15};
  }
  if (varname == "xBj") {
    return plot_data{50,
                     0.,
                     2.0,
                     "#it{x}_{Bj}",
                     "d#sigma/d#it{x}_{Bj} (10^{#minus 38} cm^{2}/nucleon)",
                     .5,
                     .2};
  }
  if (varname == "W") {
    return plot_data{60,
                     0.8,
                     4.0,
                     "#it{W} (GeV)",
                     "d#sigma/d#it{W} (10^{#minus 38} cm^{2}/GeV/nucleon)",
                     2,
                     .15};
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

// constexpr std::array<int, 12> col{1014, 1003, 1002, 1009, 1010, 1007,
//                                   1012, 1005, 1011, 1008, 1014, 1017};
// auto col = GetColorArray(8);

std::tuple<THStack, TLegend> build_stack_from_list(
    std::vector<std::tuple<std::string, ROOT::RDF::RResultPtr<TH1>, long>> list,
    double threshold) {
  IniColorCB();
  std::tuple<THStack, TLegend> tup;
  auto &&[stack, legend] = tup;
  size_t skipped = 0;
  for (auto &&[name, plot, color_diff] : list) {
    // auto &&[name, plot] = plot_entry;
    plot->SetLineColor(fgkColorBase + 1 + color_diff);
    plot->SetFillColor(fgkColorBase + 1 + color_diff);
    // plot->SetFillStyle(1001);

    stack.Add(plot.GetPtr());
    if (plot->Integral("WIDTH") > threshold)
      legend.AddEntry(plot.GetPtr(), name.c_str(), "f");
    else
      skipped++;
  }
  return tup;
}

void do_plot(std::vector<plot_ptr_t> plot_ptrs_list, std::string filename,
             std::string_view ytitle, std::string_view xtitle,
             std::array<double, 4> legend_pos, double xmax = 0.,
             std::string legend_head = "", std::string histopt = "HISTC",
             double ymax = 0.) {
  auto c = getCanvas();
  PadSetup(c.get());
  for (auto [id, obj] : plot_ptrs_list | std::views::enumerate) {
    // std::string hist_draw_opt = id == 0 ? "E0 X1" : "HISTC same";
    std::string hist_draw_opt = id == 0 ? histopt : (histopt + " same");
    std::visit(
        [&](auto &&arg) {
          if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, TH1 *> ||
                        std::is_same_v<std::decay_t<decltype(arg)>,
                                       ROOT::RDF::RResultPtr<TH1>>) {
            ResetStyle(arg);
            arg->SetStats(0);
            arg->GetXaxis()->SetTitle(xtitle.data());
            arg->GetYaxis()->SetTitle(ytitle.data());
            arg->SetMinimum(0);
            if (xmax > 0.) {
              arg->GetXaxis()->SetRangeUser(0., xmax);
            }
            if (ymax) {
              arg->SetMaximum(ymax);
            }
            if (id == 0) {
              arg->SetFillColor(kBlack);
              arg->SetLineColor(kBlack);
              arg->SetLineWidth(2);
              arg->SetMarkerStyle(20);
            }
            if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, TH1 *>)
              arg->Draw("E0 X1");
            else
              arg->Draw(hist_draw_opt.c_str());
          } else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>,
                                              TLegend *>) {
            ResetStyle(arg);
            arg->SetX1(legend_pos[0]);
            arg->SetY1(legend_pos[1]);
            arg->SetX2(legend_pos[2]);
            arg->SetY2(legend_pos[3]);
            if (std::holds_alternative<TH1 *>(plot_ptrs_list[0])) {
              auto datahist = std::get<TH1 *>(plot_ptrs_list[0]);
              arg->AddEntry(datahist, datahist->GetTitle(), "lep");
            }
            if (!legend_head.empty()) {
              arg->SetHeader(legend_head.data());
            }
            arg->Draw();
          } else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>,
                                              THStack *>) {
            if (id == 0) {
              auto last = dynamic_cast<TH1 *>(arg->GetStack()->Last());
              ResetStyle(last);
              last->Draw();
              last->GetXaxis()->SetTitle(xtitle.data());
              last->GetYaxis()->SetTitle(ytitle.data());
              if (ymax) {
                last->SetMaximum(ymax);
              }
              last->Draw(hist_draw_opt.c_str());
            }
            arg->Draw((hist_draw_opt + " SAME").c_str());
          }
        },
        obj);
  }
  std::visit(
      [&](auto &&obj) {
        if constexpr (std::is_same_v<std::decay_t<decltype(obj)>, TH1 *>)
          obj->Draw("E0 X1 SAME");
      },
      plot_ptrs_list[0]);
  c->SaveAs((filename + ".pdf").c_str());
  c->SaveAs((filename + ".eps").c_str());
  c->SaveAs((filename + ".svg").c_str());
}

int main(int argc, char *argv[]) {
  TH1::AddDirectory(false);
  ROOT::EnableImplicitMT();
  std::vector<std::string> files{};
  size_t n_runs{};
  for (int i = 1; i < argc; i++) {
    files.push_back(argv[i]);
    n_runs++;
  }
  ROOT::RDataFrame input("out_tree", files);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(input);
#endif
  auto d = TrackerPrepare(input);
  auto count = input.Count();
  auto channel_merged = d.Aggregate(
      [](std::set<int> &c_set, int channel) { c_set.insert(channel); },
      [](std::set<int> c_set1, const std::set<int> c_set2) {
        c_set1.insert(c_set2.begin(), c_set2.end());
        return c_set1;
      },
      "channel", std::set<int>{});

  auto doTKI_pi0 = [](ROOT::RDF::RNode d) {
    return CommonVariableDefinePI0(DoTKICut_MINERvA_pi0(d))
        .Define("npi0",
                [](ROOT::RVec<TLorentzVector> &pions) { return pions.size(); },
                {"good_pion"})
        // .Define("dtl", [](TKIVars &vars) { return vars.dpL / vars.IApN; },
        //         {"TKIVars"})
        // .Define("realpL",
        //         [](TLorentzVector InitNucleon) { return InitNucleon.Z(); },
        //         {"InitNucleon"})
        // .Define(
        //     "realpn",
        //     [](TLorentzVector InitNucleon) { return InitNucleon.CosTheta();
        //     },
        //     {"InitNucleon"})
        // .Define("p_cos_theta",
        //         [](const TLorentzVector &full_hadron) {
        //           return full_hadron.CosTheta();
        //         },
        //         {"leading_proton"})
        // .Define("p_theta",
        //         [](const TLorentzVector &full_hadron) {
        //           return full_hadron.Theta() / M_PI * 180.;
        //         },
        //         {"leading_proton"})
        .Define("dalphat", [](TKIVars &vars) { return vars.dalphat; },
                {"TKIVars"})
        // .Define("dpL", [](TKIVars &vars) { return vars.dpL; }, {"TKIVars"})
        .Define("IApN", [](TKIVars &vars) { return vars.IApN; }, {"TKIVars"});
  };

  auto doTKI_0pi = [](ROOT::RDF::RNode d) {
    return CommonVariableDefine0PI(DoTKICut_MINERvA0PI(d))
        // .Define("dtl", [](TKIVars &vars) { return vars.dpL / vars.IApN; },
        //         {"TKIVars"})
        // .Define("realpL",
        //         [](TLorentzVector InitNucleon) { return InitNucleon.Z(); },
        //         {"InitNucleon"})
        // .Define(
        //     "realpn",
        //     [](TLorentzVector InitNucleon) { return InitNucleon.CosTheta();
        //     },
        //     {"InitNucleon"})
        // .Define("p_cos_theta",
        //         [](const TLorentzVector &full_hadron) {
        //           return full_hadron.CosTheta();
        //         },
        //         {"full_hadron"})
        // .Define("p_theta",
        //         [](const TLorentzVector &full_hadron) {
        //           return full_hadron.Theta() / M_PI * 180.;
        //         },
        //         {"full_hadron"})
        .Define("dalphat", [](TKIVars &vars) { return vars.dalphat; },
                {"TKIVars"})
        // .Define("dpL", [](TKIVars &vars) { return vars.dpL; }, {"TKIVars"})
        .Define("IApN", [](TKIVars &vars) { return vars.IApN; }, {"TKIVars"});
  };

  auto rdf_pi0_after_cut = doTKI_pi0(d);

  auto rdf_0pi_after_cut = doTKI_0pi(d);

  std::list<ROOT::RDF::RResultPtr<TH1>> plots{};

  // plots for full data
  auto pred_all_IApN_pi0 = plots.emplace_back(
      make_plots(rdf_pi0_after_cut, get_binning_IApN_pi0(), "IApN", "pi0"));
  auto pred_all_IApN_0pi = plots.emplace_back(
      make_plots(rdf_0pi_after_cut, get_binning_IApN_0pi(), "IApN", "0pi"));
  auto pred_all_dalphat_pi0 = plots.emplace_back(make_plots(
      rdf_pi0_after_cut, get_binning_dalphat_pi0(), "dalphat", "pi0"));
  auto pred_all_dalphat_0pi = plots.emplace_back(make_plots(
      rdf_0pi_after_cut, get_binning_dalphat_0pi(), "dalphat", "0pi"));

  // per channel stacked plots
  // pn 0pi
  auto plots_0pi_IApN =
      plot_channels_0pi(rdf_0pi_after_cut, "IApN", get_binning_IApN_0pi());

  // dat 0pi
  auto plots_0pi_dalphat = plot_channels_0pi(rdf_0pi_after_cut, "dalphat",
                                             get_binning_dalphat_0pi());

  // pn pi0
  auto plots_pi0_IApN =
      plot_channels_pi0(rdf_pi0_after_cut, "IApN", get_binning_IApN_pi0());

  // dat pi0
  auto plots_pi0_dalphat = plot_channels_pi0(rdf_pi0_after_cut, "dalphat",
                                             get_binning_dalphat_pi0());

  // auto plots_pi0_W = get_channel_list_pi0(rdf_pi0_after_cut, "W");
  std::list<std::string> vars{"W", "Q2", "xBj"};
  auto var_plots_0pi = vars | std::views::transform([&](std::string name) {
                         return plot_channels_0pi(rdf_0pi_after_cut, name);
                       }) |
                       std::ranges::to<std::vector>();
  auto var_plots_pi0 = vars | std::views::transform([&](std::string name) {
                         return plot_channels_pi0(rdf_pi0_after_cut, name);
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
                    auto p = InitNucleon;
                    return p.Dot(q);
                  },
                  {"InitNeutrino", "PrimaryLepton", "InitNucleon"})
          .Histo2D({"deltam2_vs_pdotq", "deltam2_vs_pdotq;deltam2;pdotq", 60,
                    -0.5, 0.5, 60, 0, 5},
                   "deltam2", "pdotq", "weight");

  auto file = std::make_unique<TFile>("dtl_all.root", "RECREATE");
  file->Add(hist2d.GetPtr());
  file->Add(hist2d_another.GetPtr());
  for (auto &&plot : plots) {
    plot->Scale((12. / 13.) / n_runs / 10, "width");
    file->Add(plot.GetPtr());
  }

  constexpr double threshold_frac_0pi = 1.5e-2;
  constexpr double threshold_frac_pi0 = .5e-2;

  auto xsecint_0pi = pred_all_dalphat_0pi->Integral("WIDTH");
  auto xsecint_pi0 = pred_all_dalphat_pi0->Integral("WIDTH");
  auto [plots_0pi_IApN_stack, plots_0pi_IApN_leg] = build_stack_from_list(
      plots_0pi_IApN,
      pred_all_IApN_0pi->Integral("WIDTH") * threshold_frac_0pi);
  auto [plots_0pi_dalphat_stack, plots_0pi_dalphat_leg] = build_stack_from_list(
      plots_0pi_dalphat,
      pred_all_dalphat_0pi->Integral("WIDTH") * threshold_frac_0pi);

  auto [plots_pi0_IApN_stack, plots_pi0_IApN_leg] = build_stack_from_list(
      plots_pi0_IApN,
      pred_all_IApN_pi0->Integral("WIDTH") * threshold_frac_pi0);
  auto [plots_pi0_dalphat_stack, plots_pi0_dalphat_leg] = build_stack_from_list(
      plots_pi0_dalphat,
      pred_all_dalphat_pi0->Integral("WIDTH") * threshold_frac_pi0);

  auto stacked_vars_0pi =
      var_plots_0pi | std::views::transform([&](auto &&list) {
        return build_stack_from_list(list, xsecint_0pi * threshold_frac_0pi);
      }) |
      std::ranges::to<std::vector>();
  auto stacked_vars_pi0 =
      var_plots_pi0 | std::views::transform([&](auto &&list) {
        return build_stack_from_list(list, xsecint_pi0 * threshold_frac_pi0);
      }) |
      std::ranges::to<std::vector>();

  auto exp_data_hist_IApN_pi0 = get_IApN_hist_pi0();
  auto exp_data_hist_IApN_0pi = get_IApN_hist_0pi();
  auto exp_data_hist_dalphat_pi0 = get_dalphat_hist_pi0();
  auto exp_data_hist_dalphat_0pi = get_dalphat_hist_0pi();

  file->Add(&exp_data_hist_IApN_pi0);
  file->Add(&exp_data_hist_IApN_0pi);
  file->Add(&exp_data_hist_dalphat_pi0);
  file->Add(&exp_data_hist_dalphat_0pi);

  auto chi2_IApN_pi0 = do_chi2_IApN_pi0(pred_all_IApN_pi0.GetPtr());
  auto chi2_IApN_0pi = do_chi2_IApN_0pi(pred_all_IApN_0pi.GetPtr());
  auto chi2_dalphat_pi0 = do_chi2_dalphat_pi0(pred_all_dalphat_pi0.GetPtr());
  auto chi2_dalphat_0pi = do_chi2_dalphat_0pi(pred_all_dalphat_0pi.GetPtr());

  std::cout << "chi2_IApN_pi0 " << chi2_IApN_pi0 << std::endl;
  std::cout << "chi2_IApN_0pi " << chi2_IApN_0pi << std::endl;
  std::cout << "chi2_dalphat_pi0 " << chi2_dalphat_pi0 << std::endl;
  std::cout << "chi2_dalphat_0pi " << chi2_dalphat_0pi << std::endl;

  auto form_legend = [](TH1 *hist, double chi2) -> std::string {
    // return std::string{"GiBUU"} + " #chi^{2}/NDF = " + std::to_string(chi2);
    std::stringstream ss;
    ss << "GiBUU "
       << " #chi^{2}/NDF = " << round(chi2) << " / " << hist->GetNbinsX();
    return ss.str();
  };

  // first entry always data, then anything else
  do_plot({&exp_data_hist_IApN_pi0, pred_all_IApN_pi0, &plots_pi0_IApN_stack,
           &plots_pi0_IApN_leg},
          "IApN_pi0", get_info("IApN").ytitle, get_info("IApN").name,
          {0.6, 0.55, 0.9, 0.9}, 0.8,
          form_legend(&exp_data_hist_IApN_pi0, chi2_IApN_pi0));

  do_plot({&exp_data_hist_IApN_0pi, pred_all_IApN_0pi, &plots_0pi_IApN_stack,
           &plots_0pi_IApN_leg},
          "IApN_0pi", get_info("IApN").ytitle, get_info("IApN").name,
          {0.6, 0.6, 0.9, 0.9}, 0.8,
          form_legend(&exp_data_hist_IApN_0pi, chi2_IApN_0pi));

  do_plot({&exp_data_hist_dalphat_pi0, pred_all_dalphat_pi0,
           &plots_pi0_dalphat_stack, &plots_pi0_dalphat_leg},
          "dalphat_pi0", get_info("dalphat").ytitle, get_info("dalphat").name,
          {0.2, 0.5, 0.5, 0.9}, 0.,
          form_legend(&exp_data_hist_dalphat_pi0, chi2_dalphat_pi0));

  do_plot({&exp_data_hist_dalphat_0pi, pred_all_dalphat_0pi,
           &plots_0pi_dalphat_stack, &plots_0pi_dalphat_leg},
          "dalphat_0pi", get_info("dalphat").ytitle, get_info("dalphat").name,
          {0.2, 0.6, 0.5, 0.9}, 0.,
          form_legend(&exp_data_hist_dalphat_0pi, chi2_dalphat_0pi));

  for (auto &&[name, list] : std::views::zip(vars, stacked_vars_0pi)) {
    auto &&[stack, leg] = list;
    auto plot_ent = get_info(name);
    do_plot({&stack, &leg}, name + "_0pi", plot_ent.ytitle, plot_ent.name,
            {0.6, 0.55, 0.9, 0.9}, 0., "GiBUU MINERvA 0#pi", "HIST",
            plot_ent.ymax_0pi);
  }
  for (auto &&[name, list] : std::views::zip(vars, stacked_vars_pi0)) {
    auto &&[stack, leg] = list;
    auto plot_ent = get_info(name);
    do_plot({&stack, &leg}, name + "_pi0", plot_ent.ytitle, plot_ent.name,
            {0.6, 0.55, 0.9, 0.9}, 0., "GiBUU MINERvA #pi^{0}", "HIST",
            plot_ent.ymax_pi0);
  }

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
