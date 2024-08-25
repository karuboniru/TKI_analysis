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
#include <memory>
#include <string>
#include <vector>

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

template <typename T>
std::vector<ROOT::RDF::RResultPtr<TH1>>
make_plots_pi0(T &&df_in, std::string variable = "dtl",
               std::string prefix = "") {
  auto plot_data = get_info(variable);
  auto muon_cut = df_in.Filter(
      [](const TLorentzVector &muon) {
        return muon.P() > 1.5 && muon.P() < 20.0 &&
               muon.Theta() < 25 * M_PI / 180;
      },
      {"muon_p"});
  auto proton_momentum_min_cut = muon_cut.Filter(
      [](const TLorentzVector &pproton) { return pproton.P() > 0.45; },
      {"leading_proton"});
  auto make_DTL_plot = [&](auto &&df_in, std::string plot_name) {
    auto h =
        df_in.template Histo1D<double>({plot_name.c_str(), plot_name.c_str(),
                                        200, plot_data.xmin, plot_data.xmax},
                                       variable, "weight");
    return h;
  };
  std::vector<ROOT::RDF::RResultPtr<TH1>> plots{};
  plots.emplace_back(make_DTL_plot(df_in, prefix + variable + "_nocut"));
  plots.emplace_back(make_DTL_plot(muon_cut, prefix + variable + "_muon_cut"));
  plots.emplace_back(make_DTL_plot(
      proton_momentum_min_cut, prefix + variable + "_proton_momentum_min_cut"));
  return plots;
}

template <typename T>
std::vector<ROOT::RDF::RResultPtr<TH1>>
make_plots_0pi(T &&df_in, std::string variable = "dtl",
               std::string prefix = "") {
  auto plot_data = get_info(variable);
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
    auto h =
        df_in.template Histo1D<double>({plot_name.c_str(), plot_name.c_str(),
                                        200, plot_data.xmin, plot_data.xmax},
                                       variable, "weight");
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

void plot_and_save_2d(TH2 *plot, std::string path_prefix = "2dplots_all/") {
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
  ROOT::RDataFrame input("treeout", files);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(input);
#endif
  auto d = NuWroPrepare(input, true);
  auto d_nofsi = NuWroPrepare(input, false);

  auto count = input.Count();

  auto doTKI_pi0 = [](ROOT::RDF::RNode d) {
    return CommonVariableDefinePI0(DoTKICut_MINERvA_pi0(d))
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
                {"leading_proton"})
        .Define("p_theta",
                [](const TLorentzVector &full_hadron) {
                  return full_hadron.Theta() / M_PI * 180.;
                },
                {"leading_proton"})
        .Define("dalphat", [](TKIVars &vars) { return vars.dalphat; },
                {"TKIVars"})
        .Define("dpL", [](TKIVars &vars) { return vars.dpL; }, {"TKIVars"})
        .Define("IApN", [](TKIVars &vars) { return vars.IApN; }, {"TKIVars"});
  };

  auto doTKI_0pi = [](ROOT::RDF::RNode d) {
    return CommonVariableDefine0PI(DoTKICut_MINERvA0PI(d))
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
  };

  std::vector<ROOT::RDF::RResultPtr<TH1>> plots{};

  auto plot_pi0 = [&count, &plots](ROOT::RDF::RNode df_in,
                                   std::string plottag) {
    {
      auto plots_new = make_plots_pi0(df_in, "dtl", plottag);
      std::move(plots_new.begin(), plots_new.end(), std::back_inserter(plots));
    }
    {
      auto plots_new = make_plots_pi0(df_in, "dpL", plottag);
      std::move(plots_new.begin(), plots_new.end(), std::back_inserter(plots));
    }
    {
      auto plots_dat = make_plots_pi0(df_in, "dalphat", plottag);
      std::move(plots_dat.begin(), plots_dat.end(), std::back_inserter(plots));
    }
    {
      auto plots_pn = make_plots_pi0(df_in, "IApN", plottag);
      std::move(plots_pn.begin(), plots_pn.end(), std::back_inserter(plots));
    }
    {
      auto plot_2d = make_2D_plots(df_in, "2D" + plottag);
      std::move(plot_2d.begin(), plot_2d.end(), std::back_inserter(plots));
    }
    {
      auto plot_2d = make_2D_plots(
          df_in.Filter([](event &e) { return e.flag.mec; }, {"e"}),
          "2DMEC" + plottag);
      std::move(plot_2d.begin(), plot_2d.end(), std::back_inserter(plots));
    }
  };

  auto plot_0pi = [&count, &plots](ROOT::RDF::RNode df_in,
                                   std::string plottag) {
    {
      auto plots_new = make_plots_0pi(df_in, "dtl", plottag);
      std::move(plots_new.begin(), plots_new.end(), std::back_inserter(plots));
    }
    {
      auto plots_new = make_plots_0pi(df_in, "dpL", plottag);
      std::move(plots_new.begin(), plots_new.end(), std::back_inserter(plots));
    }
    {
      auto plots_dat = make_plots_0pi(df_in, "dalphat", plottag);
      std::move(plots_dat.begin(), plots_dat.end(), std::back_inserter(plots));
    }
    {
      auto plots_pn = make_plots_0pi(df_in, "IApN", plottag);
      std::move(plots_pn.begin(), plots_pn.end(), std::back_inserter(plots));
    }
    {
      auto plot_2d = make_2D_plots(df_in, "2D" + plottag);
      std::move(plot_2d.begin(), plot_2d.end(), std::back_inserter(plots));
    }
    {
      auto plot_2d = make_2D_plots(
          df_in.Filter([](event &e) { return e.flag.mec; }, {"e"}),
          "2DMEC" + plottag);
      std::move(plot_2d.begin(), plot_2d.end(), std::back_inserter(plots));
    }
  };
  plot_pi0(doTKI_pi0(d), "pi0_fsi");
  plot_pi0(doTKI_pi0(d_nofsi), "pi0_nofsi");
  plot_0pi(doTKI_0pi(d), "0pi_fsi");
  plot_0pi(doTKI_0pi(d_nofsi), "0pi_nofsi");

  auto file = std::make_unique<TFile>("dtl_all.root", "RECREATE");
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