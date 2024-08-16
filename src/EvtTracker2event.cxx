#include "EvtTracker2event.h"
#include "tkievent.h"
#include <TLorentzVector.h>

ROOT::RDF::RNode TrackerPrepare(ROOT::RDF::RNode df) {
  return df
      .Define("TKIEvent",
              [](int StdHepN, const ROOT::RVec<int> &StdHepPdg,
                 const ROOT::RVec<int> &StdHepStatus,
                 const ROOT::RVec<double> &StdHepP4_) {
                TKIEvent e{};
                for (int i = 0; i < StdHepN; i++) {
                  TLorentzVector p4{StdHepP4_[i * 4], StdHepP4_[i * 4 + 1],
                                    StdHepP4_[i * 4 + 2], StdHepP4_[i * 4 + 3]};
                  switch (StdHepStatus[i]) {
                  case 0:
                    e.add_in(StdHepPdg[i], p4);
                    break;
                  case 1:
                    e.add_post(StdHepPdg[i], p4);
                    e.add_out(StdHepPdg[i], p4);
                    break;
                  default:
                    throw;
                  }
                }
                return e;
              },
              {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4"})
      .Define("targetA", []() { return 12; }, {})
      .Define("targetZ", []() { return 6; }, {})
      .Define("InitNucleon",
              [](const ROOT::RVec<double> &StdHepP4_) {
                return TLorentzVector{StdHepP4_[4], StdHepP4_[5], StdHepP4_[6],
                                      StdHepP4_[7]};
              },
              {"StdHepP4"})
      .Define("InitNeutrino",
              [](const ROOT::RVec<double> &StdHepP4_) {
                return TLorentzVector{StdHepP4_[0], StdHepP4_[1], StdHepP4_[2],
                                      StdHepP4_[3]};
              },
              {"StdHepP4"})
      .Define("PrimaryLepton",
              [](const ROOT::RVec<double> &StdHepP4_) {
                return TLorentzVector{StdHepP4_[8], StdHepP4_[9], StdHepP4_[10],
                                      StdHepP4_[11]};
              },
              {"StdHepP4"})
      .Filter([](double weight) { return abs(weight) <= 1; }, {"weight"},
              "sane weight");
}
