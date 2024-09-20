#include "NuWro2event.h"
#include "tkigeneral.h"

#include <TLorentzVector.h>
#include <TRandom.h>
#include <TRandom3.h>

NeutrinoEvent NuWro2event_nofsi(event &e1) {
  NeutrinoEvent e2;
  for (const auto &p : e1.in) {
    e2.add_in(p.pdg, TLorentzVector{p.x / 1000., p.y / 1000., p.z / 1000.,
                                    p.t / 1000.});
  }
  for (const auto &p : e1.out) {
    e2.add_out(p.pdg, TLorentzVector{p.x / 1000., p.y / 1000., p.z / 1000.,
                                     p.t / 1000.});
  }
  for (const auto &p : e1.out) {
    e2.add_post(p.pdg, TLorentzVector{p.x / 1000., p.y / 1000., p.z / 1000.,
                                      p.t / 1000.});
  }
  return e2;
}

NeutrinoEvent NuWro2event(event &e1) {
  NeutrinoEvent e2;
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

ROOT::RDF::RNode NuWroPrepare(ROOT::RDF::RNode df, bool fsi) {
  // #if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  //   ROOT::RDF::Experimental::AddProgressBar(df);
  // #endif
  return df.Define("EventRecord", fsi ? NuWro2event : NuWro2event_nofsi, {"e"})
      .Define("targetA",
              [](event &e) { return e.par.nucleus_n + e.par.nucleus_p; }, {"e"})
      .Define("targetZ", [](event &e) { return e.par.nucleus_p; }, {"e"})
      .Define("InitNucleon",
              [](event &e) {
                const auto &p = e.in[1].p4();
                return TLorentzVector{p.x / 1000., p.y / 1000., p.z / 1000.,
                                      p.t / 1000.};
              },
              {"e"})
      .Define("InitNeutrino",
              [](event &e) {
                const auto &p = e.in[0].p4();
                return TLorentzVector{p.x / 1000., p.y / 1000., p.z / 1000.,
                                      p.t / 1000.};
              },
              {"e"})
      .Define("PrimaryLepton",
              [](event &e) {
                const auto &p = e.out[0].p4();
                return TLorentzVector{p.x / 1000., p.y / 1000., p.z / 1000.,
                                      p.t / 1000.};
              },
              {"e"});
}
