#include "ROOT/RVec.hxx"
#include <ROOT/RDataFrame.hxx>
#include <TF1.h>
#include <TFile.h>
#include <TSpline.h>

#include <array>
#include <fstream>
#include <generator>
#include <vector>

constexpr auto e_bins = std::to_array(
    {0.2, 0.54, 0.705, 0.805, 0.92, 1.05, 1.2, 1.375, 1.57, 2.05, 4.0});
constexpr auto E_center = std::to_array(
    {0.3818, 0.622, 0.7546, 0.8615, 0.9833, 1.122, 1.282, 1.463, 1.735, 2.619});

std::generator<std::pair<double, double>> iter_flux_file(std::string file) {
  std::ifstream infile(file);
  if (!infile.is_open()) {
    throw std::runtime_error("Could not open flux file: " + file);
  }
  for (std::string line; std::getline(infile, line);) {
    if (line.empty() || line[0] == '#')
      continue; // Skip empty lines and comments
    std::istringstream iss(line);
    double energy{}, flux{};
    if (!(iss >> energy >> flux)) {
      throw std::runtime_error("Error reading line: " + line);
    }
    co_yield std::make_pair(energy, flux);
  }
}

auto get_flux_spline() {
  const std::string map_file = BASEPATH "/data/microboone.flux";
  std::vector<double> energies{0}, fluxes{0};
  try {
    for (const auto &[energy, flux] : iter_flux_file(map_file)) {
      energies.push_back(energy);
      fluxes.push_back(flux);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error("Error reading flux file: " +
                             std::string(e.what()));
  }
  if (energies.empty() || fluxes.empty()) {
    throw std::runtime_error("No data found in flux file: " + map_file);
  }
  return TSpline5("flux_spline", energies.data(), fluxes.data(),
                  energies.size());
}

int main(int argc, char **argv) {
  if (argc == 1) {
    std::cerr << "Usage: " << argv[0] << " <input_files...>\n";
    return 1;
  }
  std::vector<std::string> input_files{argv + 1, argv + argc};
  auto n_files = input_files.size();
  auto df = ROOT::RDataFrame{"out_tree", input_files}.Define(
      "enu", [](const ROOT::RVecD &v) { return v[3]; }, {"StdHepP4"});
  auto hist = df.Histo1D({"enu_hist", "Neutrino Energy Distribution",
                          e_bins.size() - 1, e_bins.data()},
                         "enu", "weight");
  // auto total_weight = df.Sum("weight");
  hist->Scale(1. / n_files / 10.);

  auto flux_spline = get_flux_spline();
  TF1 func{"flux_func",
           [&flux_spline](const double *x, const double *p) {
             return flux_spline.Eval(x[0]);
           },
           0, 6, 0};
  auto flux_int = func.Integral(0, 6);
  TH1D flux_hist_binned{"flux_hist", "Flux Distribution", e_bins.size() - 1,
                        e_bins.data()};
  for (size_t i = 0; i < e_bins.size() - 1; ++i) {
    int bin_id = i + 1; // ROOT bins are 1-indexed
    double bin_low = flux_hist_binned.GetBinLowEdge(bin_id);
    double bin_high = flux_hist_binned.GetBinLowEdge(bin_id + 1);
    auto flux_int_bin = func.Integral(bin_low, bin_high, 1e-8);
    flux_hist_binned.SetBinContent(bin_id, flux_int_bin / flux_int);
    flux_hist_binned.SetBinError(bin_id, 0);
  }

  hist->Divide(&flux_hist_binned);
  for (int i = 0; i < hist->GetNbinsX(); ++i) {
    double bin_content = hist->GetBinContent(i + 1);
    hist->SetBinContent(i + 1, bin_content / E_center[i]);
  }
  hist->SaveAs("enu_flux_hist.root");
}
