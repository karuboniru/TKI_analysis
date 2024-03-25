#include "MINERvA_tki_cut.h"
#include "T2K_tki_cut.h"
#include "tkievent.h"
#include <TDatabasePDG.h>
#include <TLorentzVector.h>

bool acceptEventMINERvAPi0(const TKIEvent &e) {
  for (auto &&id : e.get_ids_post()) {
    auto absid = std::abs(id);
    if ((absid != 13) && (absid != 2212) && (absid != 2112) && (absid != 111) &&
        ((absid > 99 && absid < 1000) || (absid == 22 || absid == 11) ||
         (absid > 3000 && absid < 5000) || (absid == 2103) ||
         (absid == 2203))) {
      return false;
    }
  }

  return true;
}

TLorentzVector get_full_hadron_MINERvAPi0(const TKIEvent &e) {
  return e.get_leading(2212) + e.get_leading(111);
}

ROOT::RDF::RNode DoTKICut_MINERvA(ROOT::RDF::RNode df) {
  return df
      .Define("good_muon",
              [](TKIEvent &event) {
                ROOT::RVec<TLorentzVector> muons{};
                for (auto &&[id, particle] : event.post_range(13)) {
                  auto p = particle.P();
                  auto theta = particle.Theta();
                  if (p > 1.5 && p < 20.0 && theta < 25 * M_PI / 180) {
                    muons.push_back(particle);
                  }
                }
                return muons;
              },
              {"TKIEvent"})
      .Filter(
          [](ROOT::RVec<TLorentzVector> &muons) { return muons.size() == 1; },
          {"good_muon"}, "1 mu-")
      .Define("good_proton",
              [](TKIEvent &event) {
                ROOT::RVec<TLorentzVector> protons{};
                for (auto &&[id, particle] : event.post_range(2212)) {
                  auto p = particle.P();
                  auto theta = particle.Theta();
                  if (p > 0.450) {
                    protons.push_back(particle);
                  }
                }
                return protons;
              },
              {"TKIEvent"})
      .Filter(
          [](ROOT::RVec<TLorentzVector> &protons) {
            return protons.size() >= 1;
          },
          {"good_proton"}, "at least one proton")
      .Define("good_pion",
              [](TKIEvent &event) {
                ROOT::RVec<TLorentzVector> pions{};
                for (auto &&[id, particle] : event.post_range(111)) {
                  pions.push_back(particle);
                }
                return pions;
              },
              {"TKIEvent"}) // all pi0(s) are good
      .Filter(
          [](ROOT::RVec<TLorentzVector> &pions) { return pions.size() >= 1; },
          {"good_pion"}, "1+ pi0")
      .Filter(acceptEventMINERvAPi0, {"TKIEvent"}, "no other meson");
}
