#include "MINERvA_tki_cut.h"
#include "NuWro2event.h"
#include "plottools.hxx"
#include "tkigeneral.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>

#include <filesystem>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

template <typename T>
std::vector<ROOT::RDF::RResultPtr<TH1>>
make_plots(T &&df_in, std::string variable = "dtl", std::string prefix = "") {
  auto muon_cut = df_in.Filter(
      [](const TLorentzVector &muon) {
        return muon.P() > 1.5 && muon.P() < 10.0 &&
               muon.Theta() < 20 * M_PI / 180;
      },
      {"muon_p"});
  auto proton_momentum_min_cut = muon_cut.Filter(
      [](const TLorentzVector &pproton) { return pproton.P() > 0.45; },
      {"full_hadron"});
  auto proton_momentum_min_max_cut = proton_momentum_min_cut.Filter(
      [](const TLorentzVector &pproton) { return pproton.P() < 1.2; },
      {"full_hadron"});
  auto proton_momentum_min_cut_theta_cut = proton_momentum_min_max_cut.Filter(
      [](const TLorentzVector &pproton) {
        return pproton.Theta() < 70 * M_PI / 180;
      },
      {"full_hadron"});
  auto make_DTL_plot = [&](auto &&df_in, std::string plot_name) {
    auto h = df_in.template Histo1D<double>(
        {plot_name.c_str(), plot_name.c_str(), 200, -1., 1.}, "dtl", "weight");
    return h;
  };
  std::vector<ROOT::RDF::RResultPtr<TH1>> plots{};
  plots.emplace_back(make_DTL_plot(df_in, prefix + variable + "_nocut"));
  plots.emplace_back(make_DTL_plot(muon_cut, prefix + variable + "_muon_cut"));
  plots.emplace_back(make_DTL_plot(
      proton_momentum_min_cut, prefix + variable + "_proton_momentum_min_cut"));
  plots.emplace_back(
      make_DTL_plot(proton_momentum_min_max_cut,
                    prefix + variable + "_proton_momentum_min_max_cut"));
  plots.emplace_back(
      make_DTL_plot(proton_momentum_min_cut_theta_cut,
                    prefix + variable + "_proton_momentum_min_cut_theta_cut"));
  return plots;
}

struct plot_data {
  int bins;
  double xmin, xmax;
  std::string name;
};

plot_data get_info(std::string varname) {
  if (varname == "dalphat") {
    return plot_data{80, 0., 180., "#delta#alpha_{T}"};
  }
  if (varname == "IApN") {
    return plot_data{80, 0., 0.8, "p_{N}"};
  }
  if (varname == "dpL") {
    return plot_data{80, -1., 1., "#delta p_{L}"};
  }
  if (varname == "p_cos_theta") {
    return plot_data{80, -1., 1., "cos #theta_{p}"};
  }
  if (varname == "p_theta") {
    return plot_data{80, 0., 180., "#theta_{p}"};
  }
  if (varname == "dtl") {
    return plot_data{80, -1., 1., "cos #delta#alpha_{L}"};
  }
  throw std::runtime_error("Unknown variable name");
}

std::vector<ROOT::RDF::RResultPtr<TH1>> make_2D_plots(ROOT::RDF::RNode df_in,
                                                      std::string prefix = "") {
  const std::array<std::string, 5> vars{"p_theta", "dalphat", "IApN", "dpL",
                                        "dtl"};
  std::vector<ROOT::RDF::RResultPtr<TH1>> plots{};
  for (size_t i{}; i < vars.size(); i++) {
    for (size_t j{i + 1}; j < vars.size(); j++) {
      auto var1 = vars[i];
      auto var2 = vars[j];
      auto info1 = get_info(var1);
      auto info2 = get_info(var2);
      auto h = df_in.Histo2D<double, double>(
          {(prefix + var1 + "_" + var2).c_str(),
           (prefix + ";" + info1.name + ";" + info2.name).c_str(), info1.bins,
           info1.xmin, info1.xmax, info2.bins, info2.xmin, info2.xmax},
          var1.c_str(), var2.c_str(), "weight");
      plots.emplace_back(h);
    }
  }
  return plots;
}

void plot_and_save_2d(TH2 *plot, std::string path_prefix = "2dplots/") {
  thread_local style::global_style style{};

  auto path_dir = std::filesystem::path(path_prefix + "test");
  auto dirname = path_dir.parent_path();
  if (!std::filesystem::exists(dirname)) {
    std::filesystem::create_directories(dirname);
  }
  auto canvas = getCanvas();
  ResetStyle(plot);
  // plot->SetMarkerSize(2.5);
  plot->Draw("colz");
  canvas->Print((path_prefix + plot->GetName() + ".pdf").c_str());
  canvas->Print((path_prefix + plot->GetName() + ".svg").c_str());
  canvas->Print((path_prefix + plot->GetName() + ".eps").c_str());
}

int main(int argc, char *argv[]) {
  gSystem->ResetSignal(kSigBus);
  gSystem->ResetSignal(kSigSegmentationViolation);
  gSystem->ResetSignal(kSigIllegalInstruction);
  TH1::AddDirectory(false);
  ROOT::EnableImplicitMT();

  std::vector<std::string> files{};
  for (int i = 1; i < argc; i++) {
    files.push_back(argv[i]);
  }
  auto d = NuWroPrepare(ROOT::RDataFrame("treeout", files));

  // auto weight = d.Mean("weight");
  auto count = d.Count();

  auto d_TKIResult =
      CommonVariableDefine0PI(DoTKICut_MINERvA0PI(d))
          .Define("dtl", [](TKIVars &vars) { return vars.dpL / vars.IApN; },
                  {"TKIVars"})
          .Define("realpL",
                  [](event &e) {
                    // return e.in[1].p4();
                    auto p4 = e.in[1].p4();
                    return p4.z;
                  },
                  {"e"})
          .Define("realpn",
                  [](event &e) {
                    // return e.in[1].p4();
                    auto p4 = e.in[1].p4();
                    return TVector3{p4.x, p4.y, p4.z}.CosTheta();
                  },
                  {"e"})
          .Define("outsize",
                  [](event &e) {
                    // return e.in[1].p4();
                    // auto p4 = [1].p4();
                    return e.out.size();
                  },
                  {"e"})
          .Define("p_cos_theta",
                  [](const TLorentzVector &full_hadron) {
                    return full_hadron.CosTheta();
                  },
                  {"full_hadron"})
          .Define("p_theta",
                  [](const TLorentzVector &full_hadron) {
                    return full_hadron.Theta() / M_PI * 180.;
                  },
                  {"full_hadron"})
          .Define("dalphat", [](TKIVars &vars) { return vars.dalphat; },
                  {"TKIVars"})
          .Define("dpL", [](TKIVars &vars) { return vars.dpL; }, {"TKIVars"})
          .Define("IApN", [](TKIVars &vars) { return vars.IApN; }, {"TKIVars"});

  auto count_after_cut = d_TKIResult.Count();

  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "RECREATE";
  opts.fLazy = true;

  auto action =
      d_TKIResult.Snapshot("tree", "output0pi.root",
                           {"TKIVars", "weight", "targetA", "targetZ",
                            "neutrino_p", "muon_p", "full_hadron", "dpl_alt",
                            "realpn", "factor", "realpL", "flag", "outsize"},
                           opts);

  // const double MINERvA_pi0_IApN_bin_edges[] = {
  //     0,     0.055, 0.11,  0.165, 0.22,  0.275, 0.33,
  //     0.385, 0.44,  0.495, 0.56,  0.655, 0.81};
  // auto h =
  //     d_TKIResult
  //         .Define("IApN", [](TKIVars &vars) { return vars.IApN; },
  //         {"TKIVars"}) .Histo1D<double>({"h", "h", 12,
  //         MINERvA_pi0_IApN_bin_edges}, "IApN",
  //                          "weight");
  // auto report = d_TKIResult.Report();
  // report->Print();
  // h->Scale((12. / 13.) / count.GetValue(), "width");
  // h->SaveAs("IApN.root");
  auto plots = make_plots(d_TKIResult);

  {
    auto plots_dat = make_plots(d_TKIResult, "dalphat");
    std::move(plots_dat.begin(), plots_dat.end(), std::back_inserter(plots));
  }
  {
    auto plots_pn = make_plots(d_TKIResult, "IApN");
    std::move(plots_pn.begin(), plots_pn.end(), std::back_inserter(plots));
  }
  {
    auto plots_QE = make_plots(
        d_TKIResult.Filter(
            [](TKIEvent &event) {
              return event.count_out(2212) + event.count_out(2112) == 1;
            },
            {"TKIEvent"}),
        "1p1h");
    auto plots_2p2h = make_plots(
        d_TKIResult.Filter(
            [](TKIEvent &event) {
              return event.count_out(2212) + event.count_out(2112) > 1;
            },
            {"TKIEvent"}),
        "2p2h");

    std::move(plots_QE.begin(), plots_QE.end(), std::back_inserter(plots));
    std::move(plots_2p2h.begin(), plots_2p2h.end(), std::back_inserter(plots));
  }

  {
    auto plots_QE = make_plots(
        d_TKIResult.Filter([](event &e) { return e.flag.qel; }, {"e"}), "dtl",
        "QE");
    // auto plots_QE_SRC = make_plots(
    //     d_TKIResult.Filter([](event &e) { return e.flag.qel; }, {"e"}),
    //     "QE");
    auto plots_2p2h = make_plots(
        d_TKIResult.Filter([](event &e) { return e.flag.mec; }, {"e"}), "dtl",
        "MEC");

    std::move(plots_QE.begin(), plots_QE.end(), std::back_inserter(plots));
    std::move(plots_2p2h.begin(), plots_2p2h.end(), std::back_inserter(plots));
  }
  {
    auto plot_2d = make_2D_plots(d_TKIResult, "2D");
    std::move(plot_2d.begin(), plot_2d.end(), std::back_inserter(plots));
  }

  {
    auto plot_2d = make_2D_plots(
        d_TKIResult.Filter([](event &e) { return e.flag.qel; }, {"e"}),
        "2DQEL");
    std::move(plot_2d.begin(), plot_2d.end(), std::back_inserter(plots));
  }

  auto file = std::make_unique<TFile>("dtl_0pi.root", "RECREATE");
  for (auto &&plot : plots) {
    plot->Scale((12. / 13.) / count.GetValue(), "width");
    file->Add(plot.GetPtr());
    auto ptr = dynamic_cast<TH2D *>(plot.GetPtr());
    if (ptr) {
      auto normalize_x = normalize_slice(ptr, true);
      auto normalize_y = normalize_slice(ptr, false);
      normalize_x->SetMinimum(0.);
      normalize_y->SetMinimum(0.);
      file->Add(normalize_x);
      file->Add(normalize_y);
      plot_and_save_2d(ptr);
      plot_and_save_2d(normalize_x);
      plot_and_save_2d(normalize_y);
    }
  }

  file->Write();
  file->Close();

  return 0;
}