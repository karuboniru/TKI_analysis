#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <TLorentzVector.h>

#include <TMath.h>
#include <TSystem.h>
#include <TVector3.h>
#include <event1.h>

#include <functional>
#include <string>
#include <vector>

#include "plottools.hxx"

ROOT::RDF::RResultPtr<TH1>
make_dal_plot(ROOT::RDF::RNode df_in, std::string plot_name,
              std::function<bool(const event &)> cut) {
  return df_in.Filter(cut, {"e"})
      .Histo1D<double>({plot_name.c_str(), plot_name.c_str(), 200, -1., 1.},
                       "realDAL", "weight");
}

ROOT::RDF::RResultPtr<TH1>
make_daT_plot(ROOT::RDF::RNode df_in, std::string plot_name,
              std::function<bool(const event &)> cut) {
  return df_in.Filter(cut, {"e"})
      .Histo1D<double>({plot_name.c_str(), plot_name.c_str(), 200, 0., 180.},
                       "realDAT", "weight");
}

ROOT::RDF::RResultPtr<TH1>
make_dalv_dat_plot(ROOT::RDF::RNode df_in, std::string plot_name,
                   std::function<bool(const event &)> cut) {
  return df_in.Filter(cut, {"e"})
      .Histo2D<double, double>(
          {plot_name.c_str(), plot_name.c_str(), 200, -1., 1., 200, 0., 180.},
          "realDAL", "realDAT", "weight");
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
  auto d =
      ROOT::RDataFrame("treeout", files)
          .Define("realDAL",
                  [](event &e) {
                    // return e.in[1].p4();
                    auto p4 = e.in[1].p4();
                    return TVector3{p4.x, p4.y, p4.z}.CosTheta();
                  },
                  {"e"})
          .Define("realDAT",
                  [](event &e) {
                    auto lepton_final_state = e.out[0];
                    TVector2 plT_r{-lepton_final_state.p4().x,
                                   -lepton_final_state.p4().y};
                    auto p4N = e.in[1].p4();
                    TVector2 dpt{p4N.x, p4N.y};
                    return TMath::ACos(plT_r * dpt / plT_r.Mod() / dpt.Mod()) /
                           TMath::Pi() * 180.;
                  },
                  {"e"});
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(d);
#endif

  std::vector<ROOT::RDF::RResultPtr<TH1>> plots{
      make_dal_plot(d, "dal_nocut", [](const event &e) { return true; }),
      make_dal_plot(d, "dal_qel", [](const event &e) { return e.flag.qel; }),
      make_dal_plot(d, "dal_res", [](const event &e) { return e.flag.res; }),
      make_dal_plot(d, "dal_mec", [](const event &e) { return e.flag.mec; }),
      make_dal_plot(d, "dal_dis", [](const event &e) { return e.flag.dis; }),
      make_daT_plot(d, "daT_nocut", [](const event &e) { return true; }),
      make_daT_plot(d, "daT_qel", [](const event &e) { return e.flag.qel; }),
      make_daT_plot(d, "daT_res", [](const event &e) { return e.flag.res; }),
      make_daT_plot(d, "daT_mec", [](const event &e) { return e.flag.mec; }),
      make_daT_plot(d, "daT_dis", [](const event &e) { return e.flag.dis; }),
      make_dalv_dat_plot(d, "dalvdat_nocut",
                         [](const event &e) { return true; }),
      make_dalv_dat_plot(d, "dalvdat_qel",
                         [](const event &e) { return e.flag.qel; }),
      make_dalv_dat_plot(d, "dalvdat_res",
                         [](const event &e) { return e.flag.res; }),
      make_dalv_dat_plot(d, "dalvdat_mec",
                         [](const event &e) { return e.flag.mec; }),
      make_dalv_dat_plot(d, "dalvdat_dis",
                         [](const event &e) { return e.flag.dis; }),
  };

  auto file = std::make_unique<TFile>("dal_true.root", "RECREATE");

  auto count = d.Count();

  for (auto &&plot : plots) {
    plot->Scale((12. / 13.) / count.GetValue(), "width");
    plot->SetMinimum(0.);
    file->Add(plot.GetPtr());
    auto ptr = dynamic_cast<TH2D *>(plot.GetPtr());
    if (ptr) {
      file->Add(normalize_slice(ptr, true));
      file->Add(normalize_slice(ptr, false));
    }
  }

  file->Write();
  file->Close();

  return 0;
}