#pragma once
#include <TH1.h>
#include <TMatrixD.h>
#include <array>

std::array<double, 13> get_binning_IApN_pi0();
std::array<double, 25> get_binning_IApN_0pi() ;
std::array<double, 10> get_binning_dalphat_pi0();
std::array<double, 13> get_binning_dalphat_0pi();

TH1D get_IApN_hist_pi0();
TH1D get_dalphat_hist_pi0();
TH1D get_IApN_hist_0pi();
TH1D get_dalphat_hist_0pi();
double do_chi2_dalphat_pi0(TH1 *hist);
double do_chi2_IApN_pi0(TH1 *hist);
double do_chi2_IApN_0pi(TH1 *hist);
double do_chi2_dalphat_0pi(TH1 *hist) ;

TMatrixT<double> get_cov_dalphat_pi0();
TMatrixT<double> get_cov_IApN_pi0();
TMatrixT<double> get_cov_dalphat_0pi();
TMatrixT<double> get_cov_IApN_0pi();
