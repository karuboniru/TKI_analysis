#pragma once

#include <TMatrixD.h>
#include <TH1D.h>

double chi2(const TMatrixD &cov, const TH1D &data, const TH1D &mc, double unit = 1e-39);