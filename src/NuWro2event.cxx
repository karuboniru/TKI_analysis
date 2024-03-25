#include "NuWro2event.h"
#include "tkigeneral.h"

#include <TRandom.h>
#include <TRandom3.h>

TKIEvent NuWro2event(event &e1) {
  TKIEvent e2;
  for (const auto &p : e1.in) {
    e2.add_in(p.pdg, TLorentzVector{p.x / 1000., p.y / 1000., p.z / 1000.,
                                    p.t / 1000.});
  }
  for (const auto &p : e1.out) {
    e2.add_out(p.pdg, TLorentzVector{p.x / 1000., p.y / 1000., p.z / 1000.,
                                     p.t / 1000.});
  }
  for (const auto &p : e1.post) {
    e2.add_post(p.pdg, TLorentzVector{p.x / 1000., p.y / 1000., p.z / 1000.,
                                      p.t / 1000.});
  }
  return e2;
}

ROOT::RDF::RNode NuWroPrepare(ROOT::RDF::RNode df) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(df);
#endif
  return df.Define("TKIEvent", NuWro2event, {"e"})
      .Define("targetA",
              [](event &e) { return e.par.nucleus_n + e.par.nucleus_p; }, {"e"})
      .Define("targetZ", [](event &e) { return e.par.nucleus_p; }, {"e"})
      // .Define("weight", [](event &e) { return e.weight; }, {"e"})
      ;
}

ROOT::RDF::RNode CommonVariableDefine(ROOT::RDF::RNode df) {
  return df
      .Define("leading_proton",
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
              [](TLorentzVector &leading_proton, TLorentzVector &leading_pion) {
                return leading_proton + leading_pion;
              },
              {"leading_proton", "leading_pion"})
      .Define("neutrino_p",
              [](event &e) {
                return TLorentzVector{e.in[0].x / 1000., e.in[0].y / 1000.,
                                      e.in[0].z / 1000., e.in[0].t / 1000.};
              },
              {"e"})
      .Define("muon_p",
              [](event &e) {
                return TLorentzVector{e.post[0].x / 1000., e.post[0].y / 1000.,
                                      e.post[0].z / 1000., e.post[0].t / 1000.};
              },
              {"e"})
      .Define("TKIVars",
              [](TKIEvent &e, TLorentzVector &neutrino_p,
                 TLorentzVector &muon_p, TLorentzVector &full_hadron) {
                auto ret = getCommonTKI(12, 6, &neutrino_p, &muon_p, &full_hadron);
                if (ret.dalphat == -999){
                  ret.dalphat = gRandom->Uniform(0, 1)* 180;
                }
                return ret;
              },
              {"TKIEvent", "neutrino_p", "muon_p", "full_hadron"});
}