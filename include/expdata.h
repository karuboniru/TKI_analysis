#pragma once
#include <TH1.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <array>

namespace MINERvA_TKI {
namespace pi0 {
namespace IApN {
constexpr size_t dimension = 12;
std::array<double, dimension + 1> get_binning();
const TH1D &get_hist();
double do_chi2(TH1 *hist);
const TMatrixDSym &get_cov();
} // namespace IApN

namespace dalphat {
constexpr size_t dimension = 9;
std::array<double, dimension + 1> get_binning();
const TH1D &get_hist();
double do_chi2(TH1 *hist);
const TMatrixDSym &get_cov();
} // namespace dalphat
} // namespace pi0

namespace ZeroPi {
namespace IApN {
constexpr size_t dimension = 24;
std::array<double, dimension + 1> get_binning();
const TH1D &get_hist();
double do_chi2(TH1 *hist);
const TMatrixDSym &get_cov();
} // namespace IApN

namespace dalphat {
constexpr size_t dimension = 12;
std::array<double, dimension + 1> get_binning();
const TH1D &get_hist();
double do_chi2(TH1 *hist);
const TMatrixDSym &get_cov();
} // namespace dalphat
} // namespace ZeroPi
} // namespace MINERvA_TKI

namespace T2K_STK {
constexpr size_t dimension = 8;
const TH1D &get_hist();
const TMatrixDSym &get_cov();
std::array<double, dimension + 1> get_binning();
double do_chi2(TH1 *hist);
} // namespace T2K_STK

namespace MicroBooNE {
namespace pi0_momentum {
constexpr size_t dimension = 8;
std::array<double, dimension + 1> get_binning();
const TH1D &get_hist();
double do_chi2(TH1 *hist_smeared);
const TMatrixDSym &get_cov();
const TMatrixDSym &get_smear();
} // namespace pi0_momentum

namespace pi0_angular {
constexpr size_t dimension = 7;
std::array<double, dimension + 1> get_binning();
const TH1D &get_hist();
double do_chi2(TH1 *hist_smeared);
const TMatrixDSym &get_cov();
const TMatrixDSym &get_smear();
} // namespace pi0_angular
} // namespace MicroBooNE