#include "dochi2.h"
#include "expdata.h"

#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>

#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

struct observable {
  std::string key;
  std::string channel;
  std::string variable;
  std::string histogram;
  const TH1D &(*data)();
  const TMatrixDSym &(*covariance)();
  double (*absolute_chi2)(TH1 *);
};

const std::vector<observable> observables{
    {.key = "IApN_pi0",
     .channel = "pi0",
     .variable = "IApN",
     .histogram = "IApN_pi0",
     .data = &MINERvA_TKI::pi0::IApN::get_hist,
     .covariance = &MINERvA_TKI::pi0::IApN::get_cov,
     .absolute_chi2 = &MINERvA_TKI::pi0::IApN::do_chi2},
    {.key = "IApN_0pi",
     .channel = "0pi",
     .variable = "IApN",
     .histogram = "IApN_0pi",
     .data = &MINERvA_TKI::ZeroPi::IApN::get_hist,
     .covariance = &MINERvA_TKI::ZeroPi::IApN::get_cov,
     .absolute_chi2 = &MINERvA_TKI::ZeroPi::IApN::do_chi2},
    {.key = "dalphat_pi0",
     .channel = "pi0",
     .variable = "dalphat",
     .histogram = "dalphat_pi0",
     .data = &MINERvA_TKI::pi0::dalphat::get_hist,
     .covariance = &MINERvA_TKI::pi0::dalphat::get_cov,
     .absolute_chi2 = &MINERvA_TKI::pi0::dalphat::do_chi2},
    {.key = "dalphat_0pi",
     .channel = "0pi",
     .variable = "dalphat",
     .histogram = "dalphat_0pi",
     .data = &MINERvA_TKI::ZeroPi::dalphat::get_hist,
     .covariance = &MINERvA_TKI::ZeroPi::dalphat::get_cov,
     .absolute_chi2 = &MINERvA_TKI::ZeroPi::dalphat::do_chi2},
};

struct model_input {
  std::string label;
  std::filesystem::path root_file;
};

model_input parse_model(const std::string &value) {
  const auto separator = value.find('=');
  if (separator == std::string::npos || separator == 0 ||
      separator + 1 == value.size()) {
    throw std::runtime_error("Expected --model LABEL=ROOT_FILE, got: " + value);
  }
  return {.label = value.substr(0, separator),
          .root_file = std::filesystem::absolute(value.substr(separator + 1))};
}

double checked_value(double value, const std::string &description) {
  if (!std::isfinite(value)) {
    throw std::runtime_error("Non-finite value for " + description);
  }
  return value;
}

nlohmann::json calculate_model(const model_input &model) {
  TFile input(model.root_file.c_str(), "READ");
  if (input.IsZombie()) {
    throw std::runtime_error("Cannot open ROOT input: " +
                             model.root_file.string());
  }

  nlohmann::json result{
      {"input_root", model.root_file.string()},
      {"observables", nlohmann::json::object()},
  };
  for (const auto &info : observables) {
    TH1D *prediction{};
    input.GetObject(info.histogram.c_str(), prediction);
    if (!prediction) {
      throw std::runtime_error("Missing TH1D '" + info.histogram + "' in " +
                               model.root_file.string());
    }

    const auto &data = info.data();
    if (prediction->GetNbinsX() != data.GetNbinsX()) {
      throw std::runtime_error("Prediction/data bin mismatch for " + info.key);
    }
    const int absolute_ndf = prediction->GetNbinsX();
    const int shape_only_ndf = absolute_ndf - 1;
    const double absolute =
        checked_value(info.absolute_chi2(prediction),
                      model.label + "/" + info.key + "/absolute chi2");
    const auto shape_result =
        do_chi2_shape_only(info.covariance(), data, *prediction);
    const double shape_only =
        checked_value(std::get<0>(shape_result),
                      model.label + "/" + info.key + "/shape-only chi2");
    const double prediction_norm = std::get<3>(shape_result);

    result["observables"][info.key] = {
        {"channel", info.channel},
        {"variable", info.variable},
        {"histogram", info.histogram},
        {"absolute",
         {{"chi2", absolute},
          {"ndf", absolute_ndf},
          {"chi2_per_ndf", absolute / absolute_ndf}}},
        {"shape_only",
         {{"chi2", shape_only},
          {"ndf", shape_only_ndf},
          {"chi2_per_ndf", shape_only / shape_only_ndf}}},
        {"integral",
         {{"prediction",
           checked_value(prediction_norm,
                         model.label + "/" + info.key + "/integral")},
          {"data", data.Integral("WIDTH")}}},
    };
  }
  return result;
}

nlohmann::json compare_models(const std::string &reference_label,
                              const nlohmann::json &reference,
                              const std::string &candidate_label,
                              const nlohmann::json &candidate) {
  nlohmann::json result{
      {"reference", reference_label},
      {"candidate", candidate_label},
      {"delta_convention", "candidate_minus_reference"},
      {"observables", nlohmann::json::object()},
  };
  for (const auto &info : observables) {
    const auto &reference_observable = reference.at("observables").at(info.key);
    const auto &candidate_observable = candidate.at("observables").at(info.key);
    auto &comparison = result["observables"][info.key];
    for (const std::string fit : {"absolute", "shape_only"}) {
      const double reference_chi2 =
          reference_observable.at(fit).at("chi2").get<double>();
      const double candidate_chi2 =
          candidate_observable.at(fit).at("chi2").get<double>();
      const double reference_per_ndf =
          reference_observable.at(fit).at("chi2_per_ndf").get<double>();
      const double candidate_per_ndf =
          candidate_observable.at(fit).at("chi2_per_ndf").get<double>();
      comparison[fit] = {
          {"reference_chi2", reference_chi2},
          {"candidate_chi2", candidate_chi2},
          {"delta_chi2", candidate_chi2 - reference_chi2},
          {"reference_chi2_per_ndf", reference_per_ndf},
          {"candidate_chi2_per_ndf", candidate_per_ndf},
          {"delta_chi2_per_ndf", candidate_per_ndf - reference_per_ndf},
          {"lower_chi2_model",
           candidate_chi2 < reference_chi2 ? candidate_label : reference_label},
      };
    }
  }
  return result;
}

void write_json(const std::filesystem::path &output_path,
                const nlohmann::json &value) {
  if (output_path.has_parent_path()) {
    std::filesystem::create_directories(output_path.parent_path());
  }
  auto temporary_path = output_path;
  temporary_path += ".tmp";
  {
    std::ofstream output(temporary_path);
    output << value.dump(2) << '\n';
    if (!output) {
      throw std::runtime_error("Failed to write JSON output: " +
                               temporary_path.string());
    }
  }
  std::filesystem::rename(temporary_path, output_path);
}

} // namespace

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description options("Options");
  options.add_options()("help", "Show this help message")(
      "model", po::value<std::vector<std::string>>()->composing()->required(),
      "Model input as LABEL=ROOT_FILE; repeat for each model")(
      "reference", po::value<std::string>(),
      "Reference model label (defaults to the first --model)")(
      "output", po::value<std::string>()->required(), "Destination JSON file");

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
    std::vector<model_input> models;
    std::set<std::string> labels;
    for (const auto &value : values["model"].as<std::vector<std::string>>()) {
      auto model = parse_model(value);
      if (!labels.insert(model.label).second) {
        throw std::runtime_error("Duplicate model label: " + model.label);
      }
      models.push_back(std::move(model));
    }
    if (models.empty()) {
      throw std::runtime_error("At least one --model is required");
    }

    const std::string reference_label =
        values.count("reference") ? values["reference"].as<std::string>()
                                  : models.front().label;
    if (!labels.contains(reference_label)) {
      throw std::runtime_error("Unknown reference model: " + reference_label);
    }

    nlohmann::json output{
        {"schema_version", 1},
        {"program", "minerva_tki_chi2"},
        {"reference", reference_label},
        {"models", nlohmann::json::object()},
        {"comparisons", nlohmann::json::object()},
        {"notes",
         {"Each observable uses its published MINERvA covariance matrix. "
          "No combined chi2 is reported because cross-observable covariance "
          "matrices are unavailable."}},
    };
    for (const auto &model : models) {
      output["models"][model.label] = calculate_model(model);
    }
    const auto &reference = output["models"].at(reference_label);
    for (const auto &model : models) {
      if (model.label == reference_label) {
        continue;
      }
      output["comparisons"][model.label] =
          compare_models(reference_label, reference, model.label,
                         output["models"].at(model.label));
    }

    const auto output_path =
        std::filesystem::absolute(values["output"].as<std::string>());
    write_json(output_path, output);
    std::cout << "Wrote MINERvA TKI chi2 comparison to " << output_path << '\n';
  } catch (const std::exception &error) {
    std::cerr << "MINERvA TKI chi2 calculation failed: " << error.what()
              << '\n';
    return 1;
  }
  return 0;
}
