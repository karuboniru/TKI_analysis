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
#include <TLorentzVector.h>
#include <TMatrixDfwd.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

auto get_IApN_hist_pi0() {
  std::array<double, 13> bin_edges{
      0.0000,   55.0000,  110.0000, 165.0000, 220.0000, 275.0000, 330.0000,
      385.0000, 440.0000, 495.0000, 560.0000, 655.0000, 810.0000};
  for (auto &&i : bin_edges) {
    i /= 1000.;
  }
  // 10$^{-42}$ cm$^2$/MeV/c/nucleon
  // to 10$^{-38}$ cm$^2$/GeV/c/nucleon
  // factor: 10^4 / 10^3 = 10
  std::array<double, 12> exp_value{0.6, 1.3, 2.1, 2.8, 1.4, 1.2,
                                   1.1, 1.2, 1.3, 1.2, 0.8, 0.7};

  TH1D data("data", "", 12, bin_edges.data());

  for (size_t i = 0; i < 12; i++) {
    data.SetBinContent(i + 1, exp_value[i] / 10.);
  }

  return data;
}

double do_chi2_IApN_pi0(TH1 *hist) {
  auto data = get_IApN_hist_pi0();
  // unit: 10$^{-47}$ cm$^2$/MeV/c/nucleon
  // to: 10$^{-38}$ cm$^2$/GeV/c/nucleon
  // scale: 10^9 / 10e3 = 10^6
  const double matrix_element_decomp[12][12]{
      {12714, 16637, 7407, 2005, 1829, 1418, -822, -2963, -1775, -2892, -3822,
       -6706},
      {0, 15026, 18731, 13567, 1839, 781, -455, -1971, -366, -4464, -2025,
       -6245},
      {0, 0, 21091, 31461, 10158, 4100, 2114, 1556, 667, -1750, -2939, -635},
      {0, 0, 0, 25067, 18613, 10259, 5284, 949, 1992, 4051, 4716, -4554},
      {0, 0, 0, 0, 17620, 17591, 8716, 1303, 1944, -1240, -1948, -5097},
      {0, 0, 0, 0, 0, 16880, 21777, 17320, 4843, 184, 1668, 5120},
      {0, 0, 0, 0, 0, 0, 16465, 15579, 13441, 7898, 6227, -1467},
      {0, 0, 0, 0, 0, 0, 0, 18517, 19628, 13906, 4144, 5067},
      {0, 0, 0, 0, 0, 0, 0, 0, 26267, 16989, 8748, -1683},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 25355, 21168, 14467},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18649, 7763},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25157}};

  TMatrixD cov_matrix_decomp(12, 12, &matrix_element_decomp[0][0]);

  TMatrixD cov_matrix_transpose(cov_matrix_decomp);
  cov_matrix_transpose.T();

  auto cov_matrix = cov_matrix_transpose * cov_matrix_decomp;
  cov_matrix *= 1e-12;
  for (size_t i{}; i < 12; i++) {
    std::cout << i << " th: " << cov_matrix[i][i] << " ";
  }
  std::cout << std::endl;
  return chi2(cov_matrix, data, *((TH1D *)hist), 1e-6);
}

struct plot_data {
  int bins;
  double xmin, xmax;
  std::string name;
};

plot_data get_info(std::string varname) {
  if (varname == "dalphat") {
    return plot_data{30, 0., 180., "#delta#alpha_{T}"};
  }
  if (varname == "IApN") {
    return plot_data{30, 0., 0.8, "p_{N}"};
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
  auto make_plot = [&](auto &&df_in, std::string plot_name) {
    auto h =
        df_in.template Histo1D<double>({plot_name.c_str(), plot_name.c_str(),
                                        200, plot_data.xmin, plot_data.xmax},
                                       variable, "weight");
    return h;
  };
  std::vector<ROOT::RDF::RResultPtr<TH1>> plots{};
  plots.emplace_back(make_plot(df_in, prefix + variable + "_nocut"));
  plots.emplace_back(make_plot(muon_cut, prefix + variable + "_muon_cut"));
  plots.emplace_back(make_plot(proton_momentum_min_cut,
                               prefix + variable + "_proton_momentum_min_cut"));
  if (variable == "IApN") {
    auto exphist = get_IApN_hist_pi0();
    exphist.SetName("expbin_IApN_pi0");
    plots.emplace_back(proton_momentum_min_cut.template Histo1D<double>(
        exphist, "IApN", "weight"));
  }
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
  size_t nruns{};
  for (int i = 1; i < argc; i++) {
    files.push_back(argv[i]);
    nruns++;
  }
  ROOT::RDataFrame input("out_tree", files);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(input);
#endif
  auto d = TrackerPrepare(input);
  // auto d_nofsi = NuWroPrepare(input, false);
  d.Snapshot("test", "test.root", {"InitNucleon", "InitNeutrino", "weight"});
  auto count = input.Count();

  auto doTKI_pi0 = [](ROOT::RDF::RNode d) {
    return CommonVariableDefinePI0(DoTKICut_MINERvA_pi0(d))
        .Define("dtl", [](TKIVars &vars) { return vars.dpL / vars.IApN; },
                {"TKIVars"})
        .Define("realpL",
                [](TLorentzVector InitNucleon) {
                  // return e.in[1].p4();
                  // auto p4 = e.in[1].p4();
                  // return p4.z;
                  return InitNucleon.Z();
                },
                {"InitNucleon"})
        .Define("realpn",
                [](TLorentzVector InitNucleon) {
                  // return e.in[1].p4();
                  // auto p4 = e.in[1].p4();
                  // return TVector3{p4.x, p4.y, p4.z}.CosTheta();
                  return InitNucleon.CosTheta();
                },
                {"InitNucleon"})
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
                [](TLorentzVector InitNucleon) {
                  // return e.in[1].p4();
                  // auto p4 = e.in[1].p4();
                  // return p4.z;
                  return InitNucleon.Z();
                },
                {"InitNucleon"})
        .Define("realpn",
                [](TLorentzVector InitNucleon) {
                  // return e.in[1].p4();
                  // auto p4 = e.in[1].p4();
                  // return TVector3{p4.x, p4.y, p4.z}.CosTheta();
                  return InitNucleon.CosTheta();
                },
                {"InitNucleon"})
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
          df_in.Filter("channel == 35 || channel == 36"), "2DMEC" + plottag);
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
          df_in.Filter("channel == 35 || channel == 36"), "2DMEC" + plottag);
      std::move(plot_2d.begin(), plot_2d.end(), std::back_inserter(plots));
    }
  };
  auto df_pi0 = doTKI_pi0(d);
  auto df_0pi = doTKI_0pi(d);
  auto df_pi0_report = df_pi0.Report();
  auto df_0pi_report = df_0pi.Report();

  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "RECREATE";
  opts.fLazy = true;

  plot_pi0(df_pi0, "pi0_fsi");
  plot_0pi(df_0pi, "0pi_fsi");

  auto file = std::make_unique<TFile>("dtl_all.root", "RECREATE");
  for (auto &&plot : plots) {
    plot->Scale((12. / 13.) / nruns / 10, "width");
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
    if (std::string_view{plot.GetPtr()->GetName()} == "expbin_IApN_pi0") {
      std::cout << "Chi2: " << do_chi2_IApN_pi0(plot.GetPtr()) << std::endl;
    }
  }
  file->Write();
  file->Close();

  std::cout << "0pi" << std::endl;
  df_0pi_report->Print();

  std::cout << "pi0" << std::endl;
  df_pi0_report->Print();
  return 0;
}
