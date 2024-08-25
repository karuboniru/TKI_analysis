#include "NuWro2event.h"
#include "T2K_tki_cut.h"
#include "dochi2.h"
#include "tkigeneral.h"

#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>

#include <memory>
#include <string>
#include <vector>

TH1D GetHist_PN() {
  const double pn_edges[]{0, .120, .240, .600, 1.500};
  TH1D hist("xsec_pn", "xsec_pn", 4, pn_edges);
  hist.SetBinContent(1, 3.57422e-43 * 1e3);
  hist.SetBinContent(2, 3.29035e-43 * 1e3);
  hist.SetBinContent(3, 9.70798e-44 * 1e3);
  hist.SetBinContent(4, 6.08949e-45 * 1e3);
  hist.SetBinError(1, sqrt(6.15822e-87 * 1e6));
  hist.SetBinError(2, sqrt(1.26014e-86 * 1e6));
  hist.SetBinError(3, sqrt(9.91557e-88 * 1e6));
  hist.SetBinError(4, sqrt(1.60775e-89 * 1e6));
  return hist;
}

TMatrixD GetCov_PN() {
  // TMatrixD ret(4,4);
  double matrix[4][4]{{6.15822e-87, -1.36217e-87, 4.07904e-88, 3.4919e-89},
                      {-1.36217e-87, 1.26014e-86, 3.11855e-88, 5.01091e-89},
                      {4.07904e-88, 3.11855e-88, 9.91557e-88, 2.42863e-89},
                      {3.4919e-89, 5.01091e-89, 2.42863e-89, 1.60775e-89}};
  TMatrixD ret(4, 4, &matrix[0][0]);
  ret = ret * 1e6;
  return ret;
}

TH1D GetHist_dpTT() {
  double bin_edges[]{-.700, -.300, -.100, .100, .300, .700};
  TH1D hist("xsec_dpTT", "xsec_dpTT", 5, bin_edges);
  hist.SetBinContent(1, 3.23935e-44 * 1e3);
  hist.SetBinContent(2, 7.75956e-44 * 1e3);
  hist.SetBinContent(3, 4.1688e-43 * 1e3);
  hist.SetBinContent(4, 6.88565e-44 * 1e3);
  hist.SetBinContent(5, 2.99973e-44 * 1e3);

  hist.SetBinError(1, sqrt(2.05641e-88 * 1e6));
  hist.SetBinError(2, sqrt(1.41872e-87 * 1e6));
  hist.SetBinError(3, sqrt(6.30918e-87 * 1e6));
  hist.SetBinError(4, sqrt(1.34319e-87 * 1e6));
  hist.SetBinError(5, sqrt(2.19258e-88 * 1e6));
  return hist;
}

TMatrixD GetCov_dpTT() {
  double data[5][5]{
      {2.05641e-88, -2.89315e-90, 2.99082e-88, 7.62505e-89, 4.53107e-89},
      {-2.89315e-90, 1.41872e-87, 3.06741e-88, 2.32161e-88, 5.20707e-89},
      {2.99082e-88, 3.06741e-88, 6.30918e-87, 9.27811e-89, 2.98112e-88},
      {7.62505e-89, 2.32161e-88, 9.27811e-89, 1.34319e-87, 1.62104e-89},
      {4.53107e-89, 5.20707e-89, 2.98112e-88, 1.62104e-89, 2.19258e-88}};
  TMatrixD ret(5, 5, &data[0][0]);
  ret = ret * 1e6;
  return ret;
}

TH1D GetHist_daT() {
  double bin_edges[]{0, 60, 120, 180};
  TH1D hist("xsec_daT", "xsec_daT", 3, bin_edges);
  hist.SetBinContent(1, 7.77984e-43);
  hist.SetBinContent(2, 7.34878e-43);
  hist.SetBinContent(3, 8.18796e-43);
  hist.SetBinError(1, sqrt(2.51148e-86));
  hist.SetBinError(2, sqrt(2.50636e-86));
  hist.SetBinError(3, sqrt(3.61188e-86));
  return hist;
}

TMatrixD GetCov_daT() {
  double data[3][3]{{2.51148e-86, 1.45016e-86, 1.26919e-86},
                    {1.45016e-86, 2.50636e-86, 1.81899e-86},
                    {1.26919e-86, 1.81899e-86, 3.61188e-86}};
  TMatrixD ret(3, 3, &data[0][0]);
  return ret;
}

int main(int argc, char *argv[]) {
  ROOT::EnableImplicitMT();

  std::vector<std::string> files{};
  for (int i = 1; i < argc; i++) {
    files.push_back(argv[i]);
  }
  auto d = NuWroPrepare(ROOT::RDataFrame("treeout", files));

  // auto weight = d.Mean("weight");
  auto count = d.Count();

  auto d_TKICut = CommonVariableDefinePI0(DoTKICut_T2K(d));
  auto count_after_cut = d_TKICut.Count();
  auto count_H =
      d_TKICut.Filter([](int Z) { return Z == 1; }, {"targetZ"}).Count();
  auto count_C =
      d_TKICut.Filter([](int Z) { return Z == 6; }, {"targetZ"}).Count();

  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "RECREATE";
  opts.fLazy = true;

  auto action =
      d_TKICut.Snapshot("tree", "output.root",
                        {"TKIVars", "weight", "targetA", "targetZ",
                         "neutrino_p", "muon_p", "full_hadron"},
                        opts);

  const double T2K_pi0_IApN_bin_edges[] = {0, 0.12, 0.24, 0.6, 1.5};
  const double T2K_pi0_dpTT_bin_edges[] = {-0.7, -0.3, -0.1, 0.1, 0.3, 0.7};
  const double T2K_pi0_dalphaT_bin_edges[] = {0, 60, 120, 180.};

  TH1::AddDirectory(kFALSE);
  auto h_IApN =
      d_TKICut
          .Define("IApN", [](TKIVars &vars) { return vars.IApN; }, {"TKIVars"})
          .Histo1D<double>({"IApN", "IApN", 4, T2K_pi0_IApN_bin_edges}, "IApN",
                           "weight");
  auto h_dpTT =
      d_TKICut
          .Define("dpTT", [](TKIVars &vars) { return vars.dpTT; }, {"TKIVars"})
          .Histo1D<double>({"dpTT", "dpTT", 5, T2K_pi0_dpTT_bin_edges}, "dpTT",
                           "weight");
  auto h_dalphaT =
      d_TKICut
          .Define("dalphaT", [](TKIVars &vars) { return vars.dalphat; },
                  {"TKIVars"})
          .Histo1D<double>({"dalphaT", "dalphaT", 3, T2K_pi0_dalphaT_bin_edges},
                           "dalphaT", "weight");

  auto report = d_TKICut.Report();
  report->Print();
  h_IApN->Scale(1. / count.GetValue(), "width");
  h_dpTT->Scale(1. / count.GetValue(), "width");
  h_dalphaT->Scale(1. / count.GetValue(), "width");

  auto chi2_IApN = chi2(GetCov_PN(), GetHist_PN(), *(h_IApN.GetPtr()));
  auto chi2_dpTT = chi2(GetCov_dpTT(), GetHist_dpTT(), *(h_dpTT.GetPtr()));
  auto chi2_dalphaT = chi2(GetCov_daT(), GetHist_daT(), *(h_dalphaT.GetPtr()));

  std::cout << "chi2_IApN: " << chi2_IApN << std::endl;
  std::cout << "chi2_dpTT: " << chi2_dpTT << std::endl;
  std::cout << "chi2_dalphaT: " << chi2_dalphaT << std::endl;

  auto file = std::make_unique<TFile>("T2K_pi0.root", "RECREATE");
  file->Add(h_IApN.GetPtr());
  file->Add(h_dpTT.GetPtr());
  file->Add(h_dalphaT.GetPtr());
  file->Write();
  file->Close();
  return 0;
}