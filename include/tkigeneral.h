#pragma once

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <optional>

struct TKIVars {
  double dalphat, dphit, dpt, dpTT, beamCalcP, IApN, recoilM, recoilP, dpL,
      deriv;
};

class TLorentzVector;

TKIVars getCommonTKI(const int targetA, const int targetZ,
                     const TLorentzVector *tmp4pBeam,
                     const TLorentzVector *tmp4pScatter,
                     const TLorentzVector *tmp4pRecoil, std::optional<double> b = std::nullopt);

double getdpLMassless(TLorentzVector pmu, TLorentzVector p_hadron);
double get_factor_pdv(TLorentzVector pmu, TLorentzVector p_hadron);

ROOT::RDF::RNode CommonVariableDefinePI0(ROOT::RDF::RNode df);
ROOT::RDF::RNode CommonVariableDefine0PI(ROOT::RDF::RNode df);
