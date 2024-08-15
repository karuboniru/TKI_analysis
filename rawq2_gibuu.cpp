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
int main(int argc, char *argv[]) {
  gSystem->ResetSignal(kSigBus);
  gSystem->ResetSignal(kSigSegmentationViolation);
  gSystem->ResetSignal(kSigIllegalInstruction);
  TH1::AddDirectory(false);
  ROOT::EnableImplicitMT();
  std::vector<std::string> files{};
  size_t nruns{};
  for (int i = 1; i < argc; i++) {
    nruns++;
    files.push_back(argv[i]);
  }
  ROOT::RDataFrame input("out_tree", files);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(input);
#endif
  auto d =
      TrackerPrepare(input).Define("Q2",
                                   [](const TLorentzVector &InitNeutrino,
                                      const TLorentzVector &PrimaryLepton) {
                                     auto q0 = InitNeutrino - PrimaryLepton;
                                     return -q0.Mag2();
                                   },
                                   {"InitNeutrino", "PrimaryLepton"});

  std::vector<ROOT::RDF::RResultPtr<TH1>> plots{};

  auto file = std::make_unique<TFile>("q2raw.root", "RECREATE");
  auto plot = d.Histo1D({"Q2", "Q^{2}", 100, 0, 10.}, "Q2", "weight");
  plot->Scale(1. / nruns / 10., "width");
  file->Add(plot.GetPtr());
  file->Write();
  file->Close();

  return 0;
}