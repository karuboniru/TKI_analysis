#include "T2K_tki_cut.h"
#include "tkievent.h"
#include <ROOT/RVec.hxx>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>

// CC1π+Xp
bool acceptEventT2K(const NeutrinoEvent &e) {
  // no Meson
  for (auto &&id : e.get_ids_post()) {
    auto absid = std::abs(id);
    if ((absid != 211) &&
        ((absid > 99 && absid < 1000) || (absid == 22 || absid == 11) ||
         (absid > 3000 && absid < 5000) || (absid == 2103) ||
         (absid == 2203))) {
      return false;
    }
    if (id == -211) {
      return false;
    }
    // if (id != 211 && id != 2212 && id != 13 && id != 2112) {
    //   return false;
    // }
  }
  // if (e.count_post(211)!=1) return false;
  return true;
}

ROOT::RDF::RNode DoTKICut_T2K(ROOT::RDF::RNode df) {
  return df
      .Define("good_muon",
              [](NeutrinoEvent &event) {
                ROOT::RVec<TLorentzVector> muons{};
                for (auto &&[id, particle] : event.post_range(13)) {
                  auto p = particle.P();
                  auto theta = particle.Theta();
                  if (p > 0.250 && p < 7.000 && theta < 70 * M_PI / 180) {
                    muons.push_back(particle);
                  }
                }
                return muons;
              },
              {"EventRecord"})
      .Filter(
          [](ROOT::RVec<TLorentzVector> &muons) { return muons.size() == 1; },
          {"good_muon"}, "1 mu-")
      .Define("good_proton",
              [](NeutrinoEvent &event) {
                ROOT::RVec<TLorentzVector> protons{};
                for (auto &&[id, particle] : event.post_range(2212)) {
                  auto p = particle.P();
                  auto theta = particle.Theta();
                  if (p > 0.450 && p < 1.200 && theta < 70 * M_PI / 180) {
                    protons.push_back(particle);
                  }
                }
                return protons;
              },
              {"EventRecord"})
      .Filter(
          [](ROOT::RVec<TLorentzVector> &protons) {
            return protons.size() >= 1;
          },
          {"good_proton"}, "at least one proton")
      .Define("good_pion",
              [](NeutrinoEvent &event) {
                ROOT::RVec<TLorentzVector> pions{};
                for (auto &&[id, particle] : event.post_range(211)) {
                  auto p = particle.P();
                  auto theta = particle.Theta();
                  if (p > 0.150 && p < 1.200 && theta < 70 * M_PI / 180) {
                    pions.push_back(particle);
                  }
                }
                return pions;
              },
              {"EventRecord"})
      .Filter(
          [](ROOT::RVec<TLorentzVector> &pions) { return pions.size() == 1; },
          {"good_pion"}, "1 pi+")
      .Filter(acceptEventT2K, {"EventRecord"}, "no other meson");
}

bool acceptEventT2K_STK(const NeutrinoEvent &e) {
  // no Meson
  for (auto &&id : e.get_ids_post()) {
    auto absid = std::abs(id);
    if (((absid > 99 && absid < 1000) || (absid == 22 || absid == 11) ||
         (absid > 3000 && absid < 5000) || (absid == 2103) ||
         (absid == 2203))) {
      return false;
    }
  }
  // if (e.count_post(211)!=1) return false;
  return true;
}

ROOT::RDF::RNode DoTKICut_T2K_STK(ROOT::RDF::RNode df) {
  return df
      .Define("good_muon",
              [](NeutrinoEvent &event) {
                ROOT::RVec<TLorentzVector> muons{};
                for (auto &&[id, particle] : event.post_range(13)) {
                  auto p = particle.P();
                  auto ctheta = particle.CosTheta();
                  if (p > 0.250 && ctheta > -0.6) {
                    muons.push_back(particle);
                  }
                }
                return muons;
              },
              {"EventRecord"})
      .Filter(
          [](ROOT::RVec<TLorentzVector> &muons) { return muons.size() == 1; },
          {"good_muon"}, "1 mu-")
      .Filter([](NeutrinoEvent &event) { return event.count_post(2212) != 0; },
              {"EventRecord"})
      .Define("leading_proton_before_cut",
              [](NeutrinoEvent &event) { return event.get_leading(2212); },
              {"EventRecord"})
      .Filter(
          [](const TLorentzVector &proton) {
            auto p = proton.P();
            auto ctheta = proton.CosTheta();
            return p > 0.450 && p < 1.0 && ctheta > 0.4;
          },
          {"leading_proton_before_cut"})

      .Define("good_proton",
              [](const TLorentzVector &proton) {
                return ROOT::RVec<TLorentzVector>{proton};
              },
              {"leading_proton_before_cut"})
      .Filter(acceptEventT2K_STK, {"EventRecord"}, "no other meson");
}