#include "MINERvA_tki_cut.h"
#include "T2K_tki_cut.h"
#include "tkievent.h"
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TRandom.h>

bool acceptEventMINERvAPi0(const NeutrinoEvent &e) {
  for (auto &&id : e.get_ids_post()) {
    auto absid = std::abs(id);
    if ((id != 13) && (id != 2212) && (id != 2112) && (absid != 111) &&
        ((absid > 99 && absid < 1000) || (absid == 22 || absid == 11) ||
         (absid > 3000 && absid < 5000) || (absid == 2103) ||
         (absid == 2203))) {
      return false;
    }
  }

  return true;
}

bool acceptEventMINERvA0PI(const NeutrinoEvent &e) {
  for (auto &&id : e.get_ids_post()) {
    auto absid = std::abs(id);
    if ((id != 13) && (id != 2212) && (id != 2112) &&
        ((absid > 99 && absid < 1000) || (absid == 22 || absid == 11) ||
         (absid > 3000 && absid < 5000) || (absid == 2103) ||
         (absid == 2203))) {
      return false;
    }
  }

  return true;
}

ROOT::RDF::RNode DoTKICut_MINERvA_pi0(ROOT::RDF::RNode df) {
  return df
      .Define("good_muon",
              [](NeutrinoEvent &event) {
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
                  if (p > 0.450) {
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
                for (auto &&[id, particle] : event.post_range(111)) {
                  pions.push_back(particle);
                }
                return pions;
              },
              {"EventRecord"}) // all pi0(s) are good
      .Filter(
          [](ROOT::RVec<TLorentzVector> &pions) { return pions.size() >= 1; },
          {"good_pion"}, "1+ pi0")
      .Filter(acceptEventMINERvAPi0, {"EventRecord"}, "no other meson");
}

ROOT::RDF::RNode DoTKICut_MINERvA0PI(ROOT::RDF::RNode df) {
  return df
      .Define("good_muon",
              [](NeutrinoEvent &event) {
                ROOT::RVec<TLorentzVector> muons{};
                for (auto &&[id, particle] : event.post_range(13)) {
                  auto p = particle.P();
                  auto theta = particle.Theta();
                  if (p > 1.5 && p < 10.0 && theta < 20 * M_PI / 180) {
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
                  if (p > 0.450 && p < 1.2 && theta < 70 * M_PI / 180) {
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
      .Filter(acceptEventMINERvA0PI, {"EventRecord"}, "no other meson");
}

ROOT::RDF::RNode vars_define(ROOT::RDF::RNode df) {
  return df
      .Define("Q2",
              [](const TLorentzVector &InitNeutrino,
                 const TLorentzVector &PrimaryLepton) {
                auto q0 = InitNeutrino - PrimaryLepton;
                return -q0.Mag2();
              },
              {"InitNeutrino", "PrimaryLepton"})
      .Define("q0",
              [](const TLorentzVector &InitNeutrino,
                 const TLorentzVector &PrimaryLepton) {
                auto q0 = InitNeutrino.E() - PrimaryLepton.E();
                return q0;
              },
              {"InitNeutrino", "PrimaryLepton"})
      .Define("W",
              [](const TLorentzVector &InitNeutrino,
                 const TLorentzVector &PrimaryLepton,
                 const TLorentzVector &InitNucleon) {
                auto had_system = InitNucleon + InitNeutrino - PrimaryLepton;
                return had_system.M();
              },
              {"InitNeutrino", "PrimaryLepton", "InitNucleon"})
      .Define("xBj",
              [](const TLorentzVector &InitNeutrino,
                 const TLorentzVector &PrimaryLepton,
                 const TLorentzVector &InitNucleon) {
                auto q = InitNeutrino - PrimaryLepton;
                return -q.Mag2() / (2 * InitNucleon.Dot(q));
              },
              {"InitNeutrino", "PrimaryLepton", "InitNucleon"})
      .Define(
          "p_mu",
          [](const TLorentzVector &PrimaryLepton) { return PrimaryLepton.P(); },
          {"PrimaryLepton"})
      .Define("theta_mu",
              [](const TLorentzVector &PrimaryLepton) {
                return PrimaryLepton.Theta() / M_PI * 180.;
              },
              {"PrimaryLepton"})
      .Define("p_p",
              [](const TLorentzVector &leading_proton) {
                return leading_proton.P();
              },
              {"leading_proton"})
      .Define("theta_p",
              [](const TLorentzVector &leading_proton) {
                return leading_proton.Theta() / M_PI * 180.;
              },
              {"leading_proton"});
}

ROOT::RDF::RNode CommonVariableDefinePI0(ROOT::RDF::RNode df,
                                         std::optional<double> b) {
  return vars_define(
      df.Define("leading_proton",
                [](ROOT::RVec<TLorentzVector> &protons) {
                  TLorentzVector leading_proton;
                  for (auto &&proton : protons) {
                    if (proton.P() > leading_proton.P()) {
                      leading_proton = proton;
                    }
                  }
                  return leading_proton;
                },
                {"good_proton"})
          .Define("leading_pion",
                  [](ROOT::RVec<TLorentzVector> &pions) {
                    TLorentzVector leading_pion;
                    for (auto &&pion : pions) {
                      if (pion.P() > leading_pion.P()) {
                        leading_pion = pion;
                      }
                    }
                    return leading_pion;
                  },
                  {"good_pion"})
          .Define("full_hadron",
                  [](const TLorentzVector &leading_proton,
                     const TLorentzVector &leading_pion) {
                    return leading_proton + leading_pion;
                  },
                  {"leading_proton", "leading_pion"})
          .Alias("neutrino_p", "InitNeutrino")
          .Alias("muon_p", "PrimaryLepton")
          .Define("p_pion",
                  [](const TLorentzVector &leading_pion) {
                    return leading_pion.P();
                  },
                  {"leading_pion"})
          .Define("theta_pion",
                  [](const TLorentzVector &leading_pion) {
                    return leading_pion.Theta() / M_PI * 180.;
                  },
                  {"leading_pion"})
          .Define("TKIVars",
                  [=](NeutrinoEvent &e, const TLorentzVector &neutrino_p,
                      const TLorentzVector &muon_p,
                      const TLorentzVector &full_hadron) {
                    auto ret = getCommonTKI(12, 6, &neutrino_p, &muon_p,
                                            &full_hadron, b);
                    if (ret.dalphat == -999) {
                      ret.dalphat = gRandom->Uniform(0, 1) * 180;
                    }
                    return ret;
                  },
                  {"EventRecord", "neutrino_p", "muon_p", "full_hadron"}));
}

ROOT::RDF::RNode CommonVariableDefine0PI(ROOT::RDF::RNode df,
                                         std::optional<double> b) {
  return vars_define(
      df.Define("leading_proton",
                [](ROOT::RVec<TLorentzVector> &protons) {
                  TLorentzVector leading_proton;
                  for (auto &&proton : protons) {
                    if (proton.P() > leading_proton.P()) {
                      leading_proton = proton;
                    }
                  }
                  return leading_proton;
                },
                {"good_proton"})
          .Alias("full_hadron", "leading_proton")
          .Alias("neutrino_p", "InitNeutrino")
          .Alias("muon_p", "PrimaryLepton")
          .Define("TKIVars",
                  [=](NeutrinoEvent &e, TLorentzVector &neutrino_p,
                     TLorentzVector &muon_p, TLorentzVector &full_hadron) {
                    auto ret =
                        getCommonTKI(12, 6, &neutrino_p, &muon_p, &full_hadron, b);
                    if (ret.dalphat == -999) {
                      ret.dalphat = gRandom->Uniform(0, 1) * 180;
                    }
                    return ret;
                  },
                  {"EventRecord", "neutrino_p", "muon_p", "leading_proton"})
          .Define("dpl_alt",
                  [](const TLorentzVector &muon, const TLorentzVector &hadron) {
                    return getdpLMassless(muon, hadron);
                  },
                  {"muon_p", "leading_proton"})
          .Define("factor",
                  [](const TLorentzVector &muon, const TLorentzVector &hadron) {
                    return get_factor_pdv(muon, hadron);
                  },
                  {"muon_p", "leading_proton"}));
}
