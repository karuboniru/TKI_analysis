#pragma once

struct TKIVars {
  double dalphat, dphit, dpt, dpTT, beamCalcP, IApN, recoilM, recoilP;
};

class TLorentzVector;

TKIVars getCommonTKI(const int targetA, const int targetZ,
                     const TLorentzVector *tmp4pBeam,
                     const TLorentzVector *tmp4pScatter,
                     const TLorentzVector *tmp4pRecoil);

