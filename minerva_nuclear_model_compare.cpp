#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>

#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct observable {
  std::string histogram;
  std::string data_histogram;
  std::string label;
  std::string x_title;
  std::string y_title;
};

const std::vector<observable> observables{
    {.histogram = "IApN_pi0",
     .data_histogram = "data_IApN_pi0",
     .label = "MINERvA #pi^{0}",
     .x_title = "p_{N} (GeV/c)",
     .y_title = "d#sigma/dp_{N}"},
    {.histogram = "IApN_0pi",
     .data_histogram = "data_IApN_0pi",
     .label = "MINERvA 0#pi",
     .x_title = "p_{N} (GeV/c)",
     .y_title = "d#sigma/dp_{N}"},
    {.histogram = "dalphat_pi0",
     .data_histogram = "data_dalphat_pi0",
     .label = "MINERvA #pi^{0}",
     .x_title = "#delta#alpha_{T} (degree)",
     .y_title = "d#sigma/d#delta#alpha_{T}"},
    {.histogram = "dalphat_0pi",
     .data_histogram = "data_dalphat_0pi",
     .label = "MINERvA 0#pi",
     .x_title = "#delta#alpha_{T} (degree)",
     .y_title = "d#sigma/d#delta#alpha_{T}"},
};

nlohmann::json read_json(const std::filesystem::path &path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Cannot open JSON file: " + path.string());
  }
  nlohmann::json value;
  input >> value;
  return value;
}

std::unique_ptr<TH1> clone_histogram(TFile &file, const std::string &name,
                                     const std::string &clone_name,
                                     bool required = true) {
  auto *source = dynamic_cast<TH1 *>(file.Get(name.c_str()));
  if (!source) {
    if (required) {
      throw std::runtime_error("Missing histogram '" + name + "' in " +
                               file.GetName());
    }
    return nullptr;
  }
  auto result = std::unique_ptr<TH1>(
      dynamic_cast<TH1 *>(source->Clone(clone_name.c_str())));
  if (!result) {
    throw std::runtime_error("Failed to clone histogram '" + name + "'");
  }
  result->SetDirectory(nullptr);
  return result;
}

void check_compatible(const TH1 &default_hist, const TH1 &cda_hist,
                      const std::string &name) {
  if (default_hist.GetNbinsX() != cda_hist.GetNbinsX()) {
    throw std::runtime_error("Incompatible bin counts for " + name);
  }
  for (int bin = 1; bin <= default_hist.GetNbinsX() + 1; ++bin) {
    if (std::abs(default_hist.GetXaxis()->GetBinLowEdge(bin) -
                 cda_hist.GetXaxis()->GetBinLowEdge(bin)) > 1e-12) {
      throw std::runtime_error("Incompatible bin edges for " + name);
    }
  }
}

std::unique_ptr<TH1> make_ratio(const TH1 &numerator, const TH1 &denominator,
                                const std::string &name) {
  auto ratio =
      std::unique_ptr<TH1>(dynamic_cast<TH1 *>(numerator.Clone(name.c_str())));
  ratio->SetDirectory(nullptr);
  ratio->Reset("ICES");
  if (ratio->GetSumw2N() == 0) {
    ratio->Sumw2();
  }
  for (int bin = 1; bin <= numerator.GetNbinsX(); ++bin) {
    const double num = numerator.GetBinContent(bin);
    const double den = denominator.GetBinContent(bin);
    if (den == 0.) {
      ratio->SetBinContent(bin, 0.);
      ratio->SetBinError(bin, 0.);
      continue;
    }
    const double value = num / den;
    const double num_error = numerator.GetBinError(bin);
    const double den_error = denominator.GetBinError(bin);
    const double variance = std::pow(num_error / den, 2) +
                            std::pow(num * den_error / std::pow(den, 2), 2);
    ratio->SetBinContent(bin, value);
    ratio->SetBinError(bin, std::sqrt(variance));
  }
  return ratio;
}

nlohmann::json histogram_summary(const TH1 &default_hist, const TH1 &cda_hist,
                                 const TH1 &ratio_hist) {
  nlohmann::json result;
  double default_error{};
  double cda_error{};
  const double default_integral = default_hist.IntegralAndError(
      1, default_hist.GetNbinsX(), default_error, "width");
  const double cda_integral =
      cda_hist.IntegralAndError(1, cda_hist.GetNbinsX(), cda_error, "width");
  result["integral"] = {{"default", default_integral},
                        {"default_error", default_error},
                        {"use_cda", cda_integral},
                        {"use_cda_error", cda_error}};
  if (default_integral != 0.) {
    result["integral"]["ratio"] = cda_integral / default_integral;
    result["integral"]["percent_change"] =
        100. * (cda_integral / default_integral - 1.);
  } else {
    result["integral"]["ratio"] = nullptr;
    result["integral"]["percent_change"] = nullptr;
  }

  result["bins"] = nlohmann::json::array();
  for (int bin = 1; bin <= default_hist.GetNbinsX(); ++bin) {
    const double default_value = default_hist.GetBinContent(bin);
    nlohmann::json bin_value{
        {"low", default_hist.GetXaxis()->GetBinLowEdge(bin)},
        {"high", default_hist.GetXaxis()->GetBinUpEdge(bin)},
        {"default", default_value},
        {"default_error", default_hist.GetBinError(bin)},
        {"use_cda", cda_hist.GetBinContent(bin)},
        {"use_cda_error", cda_hist.GetBinError(bin)},
    };
    if (default_value != 0.) {
      bin_value["ratio"] = ratio_hist.GetBinContent(bin);
      bin_value["ratio_error"] = ratio_hist.GetBinError(bin);
      bin_value["percent_change"] = 100. * (ratio_hist.GetBinContent(bin) - 1.);
    } else {
      bin_value["ratio"] = nullptr;
      bin_value["ratio_error"] = nullptr;
      bin_value["percent_change"] = nullptr;
    }
    result["bins"].push_back(std::move(bin_value));
  }
  return result;
}

void draw_comparison(const observable &info, TH1 &default_hist, TH1 &cda_hist,
                     TH1 *data_hist, TH1 &ratio_hist,
                     const std::filesystem::path &output_dir) {
  TCanvas canvas((info.histogram + "_canvas").c_str(), "", 800, 800);
  TPad upper("upper", "", 0., 0.30, 1., 1.);
  TPad lower("lower", "", 0., 0., 1., 0.30);
  upper.SetBottomMargin(0.02);
  upper.SetLeftMargin(0.13);
  upper.SetRightMargin(0.04);
  lower.SetTopMargin(0.02);
  lower.SetBottomMargin(0.34);
  lower.SetLeftMargin(0.13);
  lower.SetRightMargin(0.04);
  upper.Draw();
  lower.Draw();

  upper.cd();
  default_hist.SetTitle("");
  default_hist.SetLineColor(kBlue + 1);
  default_hist.SetLineWidth(3);
  default_hist.SetFillStyle(0);
  cda_hist.SetLineColor(kRed + 1);
  cda_hist.SetLineWidth(3);
  cda_hist.SetLineStyle(2);
  cda_hist.SetFillStyle(0);
  default_hist.GetYaxis()->SetTitle(info.y_title.c_str());
  default_hist.GetYaxis()->SetTitleSize(0.055);
  default_hist.GetYaxis()->SetTitleOffset(1.10);
  default_hist.GetYaxis()->SetLabelSize(0.045);
  default_hist.GetXaxis()->SetLabelSize(0.);
  double maximum = std::max(default_hist.GetMaximum(), cda_hist.GetMaximum());
  if (data_hist) {
    maximum = std::max(maximum, data_hist->GetMaximum());
  }
  default_hist.SetMaximum(maximum > 0. ? maximum * 1.35 : 1.);
  default_hist.SetMinimum(0.);
  default_hist.Draw("HIST E");
  cda_hist.Draw("HIST E SAME");

  if (data_hist) {
    data_hist->SetMarkerStyle(20);
    data_hist->SetMarkerSize(0.9);
    data_hist->SetMarkerColor(kBlack);
    data_hist->SetLineColor(kBlack);
    data_hist->Draw("E1 SAME");
  }

  TLegend legend(0.57, 0.68, 0.94, 0.91);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  if (data_hist) {
    legend.AddEntry(data_hist, "MINERvA data", "lep");
  }
  legend.AddEntry(&default_hist, "GiBUU useCdA=F", "le");
  legend.AddEntry(&cda_hist, "GiBUU useCdA=T", "le");
  legend.Draw();

  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.05);
  label.DrawLatex(0.16, 0.91, info.label.c_str());

  lower.cd();
  ratio_hist.SetTitle("");
  ratio_hist.SetLineColor(kRed + 1);
  ratio_hist.SetMarkerColor(kRed + 1);
  ratio_hist.SetMarkerStyle(20);
  ratio_hist.SetMarkerSize(0.7);
  ratio_hist.GetYaxis()->SetTitle("T / F");
  ratio_hist.GetXaxis()->SetTitle(info.x_title.c_str());
  ratio_hist.GetYaxis()->SetNdivisions(505);
  ratio_hist.GetYaxis()->SetTitleSize(0.12);
  ratio_hist.GetYaxis()->SetTitleOffset(0.48);
  ratio_hist.GetYaxis()->SetLabelSize(0.10);
  ratio_hist.GetXaxis()->SetTitleSize(0.13);
  ratio_hist.GetXaxis()->SetTitleOffset(1.05);
  ratio_hist.GetXaxis()->SetLabelSize(0.10);

  double ratio_min = 1.;
  double ratio_max = 1.;
  for (int bin = 1; bin <= ratio_hist.GetNbinsX(); ++bin) {
    if (default_hist.GetBinContent(bin) == 0.) {
      continue;
    }
    ratio_min = std::min(ratio_min, ratio_hist.GetBinContent(bin));
    ratio_max = std::max(ratio_max, ratio_hist.GetBinContent(bin));
  }
  const double margin = std::max(0.08, 0.15 * (ratio_max - ratio_min));
  ratio_hist.SetMinimum(std::max(0., ratio_min - margin));
  ratio_hist.SetMaximum(ratio_max + margin);
  ratio_hist.Draw("E1");
  TLine unity(ratio_hist.GetXaxis()->GetXmin(), 1.,
              ratio_hist.GetXaxis()->GetXmax(), 1.);
  unity.SetLineStyle(2);
  unity.SetLineColor(kGray + 2);
  unity.Draw();

  canvas.cd();
  for (const std::string extension : {"pdf", "svg", "eps"}) {
    canvas.SaveAs((output_dir / (info.histogram + "." + extension)).c_str());
  }
}

} // namespace

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description options("Options");
  options.add_options()("help", "Show this help message")(
      "default-root", po::value<std::string>()->required(),
      "dtl_all.root produced with useCdA=F")(
      "cda-root", po::value<std::string>()->required(),
      "dtl_all.root produced with useCdA=T")(
      "default-json", po::value<std::string>()->required(),
      "chi2.json produced with useCdA=F")("cda-json",
                                          po::value<std::string>()->required(),
                                          "chi2.json produced with useCdA=T")(
      "output-dir", po::value<std::string>()->required(),
      "Output directory for plots, ROOT, and JSON summaries");

  po::variables_map values;
  try {
    po::store(po::parse_command_line(argc, argv, options), values);
    if (values.count("help")) {
      std::cout << options << '\n';
      return 0;
    }
    po::notify(values);
  } catch (const po::error &error) {
    std::cerr << error.what() << "\n\n" << options << '\n';
    return 1;
  }

  try {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    const auto default_root_path =
        std::filesystem::absolute(values["default-root"].as<std::string>());
    const auto cda_root_path =
        std::filesystem::absolute(values["cda-root"].as<std::string>());
    const auto default_json_path =
        std::filesystem::absolute(values["default-json"].as<std::string>());
    const auto cda_json_path =
        std::filesystem::absolute(values["cda-json"].as<std::string>());
    const auto output_dir =
        std::filesystem::absolute(values["output-dir"].as<std::string>());
    std::filesystem::create_directories(output_dir);

    TFile default_file(default_root_path.c_str(), "READ");
    TFile cda_file(cda_root_path.c_str(), "READ");
    if (default_file.IsZombie() || cda_file.IsZombie()) {
      throw std::runtime_error("Cannot open one or both ROOT inputs");
    }

    TFile comparison_file((output_dir / "comparison.root").c_str(), "RECREATE");
    if (comparison_file.IsZombie()) {
      throw std::runtime_error("Cannot create comparison.root");
    }

    nlohmann::json summary{
        {"inputs",
         {{"default_root", default_root_path.string()},
          {"use_cda_root", cda_root_path.string()},
          {"default_json", default_json_path.string()},
          {"use_cda_json", cda_json_path.string()}}},
        {"observables", nlohmann::json::object()},
    };

    for (const auto &info : observables) {
      auto default_hist = clone_histogram(default_file, info.histogram,
                                          info.histogram + "_default");
      auto cda_hist = clone_histogram(cda_file, info.histogram,
                                      info.histogram + "_use_cda");
      auto data_hist =
          clone_histogram(default_file, info.data_histogram,
                          info.data_histogram + "_comparison", false);
      check_compatible(*default_hist, *cda_hist, info.histogram);
      auto ratio_hist = make_ratio(*cda_hist, *default_hist,
                                   info.histogram + "_ratio_T_over_F");

      summary["observables"][info.histogram] =
          histogram_summary(*default_hist, *cda_hist, *ratio_hist);
      draw_comparison(info, *default_hist, *cda_hist, data_hist.get(),
                      *ratio_hist, output_dir);

      comparison_file.cd();
      default_hist->Write();
      cda_hist->Write();
      ratio_hist->Write();
      if (data_hist) {
        data_hist->Write();
      }
    }

    const auto default_chi2 = read_json(default_json_path);
    const auto cda_chi2 = read_json(cda_json_path);
    summary["chi2"] = {{"default", default_chi2}, {"use_cda", cda_chi2}};
    summary["chi2"]["delta"] = nlohmann::json::object();
    for (auto iterator = default_chi2.begin(); iterator != default_chi2.end();
         ++iterator) {
      if (iterator.value().is_number() && cda_chi2.contains(iterator.key()) &&
          cda_chi2.at(iterator.key()).is_number()) {
        summary["chi2"]["delta"][iterator.key()] =
            cda_chi2.at(iterator.key()).get<double>() -
            iterator.value().get<double>();
      }
    }

    comparison_file.Write();
    comparison_file.Close();
    std::ofstream json_output(output_dir / "comparison.json");
    json_output << summary.dump(2) << '\n';
    if (!json_output) {
      throw std::runtime_error("Failed to write comparison.json");
    }
  } catch (const std::exception &error) {
    std::cerr << "Comparison failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
