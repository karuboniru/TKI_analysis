#pragma once

#include <TH1D.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>

double chi2(const TMatrixDSym &cov, const TH1D &data, const TH1D &mc,
            double unit = 1e-39);
double chi2(const TMatrixT<double> &cov, const TH1D &data, const TH1D &mc,
            double unit = 1e-39);

// TH1D to_shape(const TH1D &input);

std::tuple<TMatrixDSym, TH1D> ExractShapeOnly(const TMatrixDSym &cov_in,
                                              const TH1D &data_in);

std::tuple<double, TH1D /*shape only prediction*/, TH1D /*shape only data*/,
           double /*norm*/>
do_chi2_shape_only(const TMatrixDSym &, const TH1D &, const TH1D &);