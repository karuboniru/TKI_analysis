#include "MINERvA_tki_cut.h"
#include "NuWro2event.h"
#include "tkigeneral.h"

#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>

#include <string>
#include <vector>

int main(int argc, char *argv[]) {
  ROOT::EnableImplicitMT();

  std::vector<std::string> files{};
  for (int i = 1; i < argc; i++) {
    files.push_back(argv[i]);
  }
  auto d = NuWroPrepare(ROOT::RDataFrame("treeout", files));

  // auto weight = d.Mean("weight");
  auto count = d.Count();

  auto d_TKIResult = CommonVariableDefine(DoTKICut_MINERvA(d));

  auto count_after_cut = d_TKIResult.Count();

  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "RECREATE";
  opts.fLazy = true;

  auto action = d_TKIResult.Snapshot("tree", "output.root",
                                     {"TKIVars", "weight", "targetA", "targetZ",
                                      "neutrino_p", "muon_p", "full_hadron"},
                                     opts);

  const double MINERvA_pi0_IApN_bin_edges[] = {
      0,     0.055, 0.11,  0.165, 0.22,  0.275, 0.33,
      0.385, 0.44,  0.495, 0.56,  0.655, 0.81};
  auto h =
      d_TKIResult
          .Define("IApN", [](TKIVars &vars) { return vars.IApN; }, {"TKIVars"})
          .Histo1D<double>({"h", "h", 12, MINERvA_pi0_IApN_bin_edges}, "IApN",
                           "weight");
  auto report = d_TKIResult.Report();
  report->Print();
  h->Scale((12. / 13.) / count.GetValue(), "width");
  h->SaveAs("IApN.root");
  return 0;
}