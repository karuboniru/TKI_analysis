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
  auto d = NuWroPrepare(input, true)
               .Define("Q2",
                       [](const TLorentzVector &InitNeutrino,
                          const TLorentzVector &PrimaryLepton) {
                         auto q0 = InitNeutrino - PrimaryLepton;
                         return -q0.Mag2();
                       },
                       {"InitNeutrino", "PrimaryLepton"});
  auto count = d.Count();
  std::vector<ROOT::RDF::RResultPtr<TH1>> plots{};

  auto file = std::make_unique<TFile>("q2raw.root", "RECREATE");
  auto plot = d.Histo1D({"Q2", "Q^{2}", 100, 0, 10.}, "Q2", "weight");
  plot->Scale(1. / count.GetValue(), "width");
  file->Add(plot.GetPtr());
  file->Write();
  file->Close();

  return 0;
}