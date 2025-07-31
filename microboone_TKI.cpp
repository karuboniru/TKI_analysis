#include "EvtTracker2event.h"
#include "MINERvA_tki_cut.h"
#include "dochi2.h"
#include "expdata.h"
#include "plottools.hxx"
#include "tkievent.h"
#include "tkigeneral.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <TH1.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMatrixDfwd.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>
#include <algorithm>
#include <boost/program_options.hpp>

#include <nlohmann/json.hpp>

#include <cmath>
#include <print>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

const auto bin_edges =
    std::to_array({0.0, 22.0, 44.0, 66.0, 88.0, 110.0, 145.0, 180.0});

// Cross sections and uncertainties
const auto cross_sections =
    std::to_array({0.047295775, 0.044962975, 0.044549171, 0.050441444,
                   0.066811717, 0.078273462, 0.090827721});

const auto uncertainties =
    std::to_array({0.0069867236, 0.0058600051, 0.0056673462, 0.0070400151,
                   0.008404137, 0.0087227652, 0.01011356});

const std::vector<std::tuple<std::string, std::function<bool(int)>>>
    list_channel_channeldef{
        {"QE", [](int id) { return id == 1; }},
        {"2p2h", [](int id) { return id == 35; }},
        {"RES", [](int id) { return id >= 2 && id <= 31; }},
        {"DIS", [](int id) { return id == 34; }},
        {"1#pi BG", [](int id) { return id == 32 || id == 33; }},
        {"2#pi BG", [](int id) { return id == 37; }},
    };

TH1D make_exp_hist() {
  TH1D exp_hist("expbin", "MicroBooNE CC0#pi;dalpha_{t} [deg];Events",
                bin_edges.size() - 1, bin_edges.data());
  for (size_t i = 0; i < cross_sections.size(); ++i) {
    exp_hist.SetBinContent(i + 1, cross_sections[i] / 40.);
    exp_hist.SetBinError(i + 1, uncertainties[i] / 40.);
  }
  return exp_hist;
}

int main(int argc, char **argv) {
  IniColorCBminerva2pibg();

  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()("input-files", po::value<std::vector<std::string>>(),
                     "Input files")(
      "add-text", po::value<std::string>()->default_value(""),
      "Additional text")("help", "produce help message");
  po::positional_options_description p;
  p.add("input-files", -1);
  po::variables_map vm;
  po::store(
      po::command_line_parser(argc, argv).options(desc).positional(p).run(),
      vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-files")) {
    std::cout << desc << std::endl;
    return 1;
  }

  auto &&files = vm["input-files"].as<std::vector<std::string>>();
  auto n_runs = files.size();
  auto additional_text = vm["add-text"].as<std::string>();

  TH1::AddDirectory(false);
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame input("out_tree", files);

  // GiBUU data preparation
  auto d = TrackerPrepare(input);
  ROOT::RDF::Experimental::AddProgressBar(d);

  auto signal =
      d.Filter(
           [](const NeutrinoEvent &event) {
             // all final state particles are muon, charged pion, proton,
             // neutron existance of any other particle leads to rejection
             if (std::ranges::any_of(event.get_ids_post(), [](int pdg) {
                   switch (pdg) {
                   case 13:   // muon
                   case 2212: // proton
                   case 2112: // neutron
                   case 211:  // pion+
                   case -211: // pion-
                     return false;
                   default:
                     return true; // other particles are not allowed
                   }
                 })) {
               return false; // reject events with other particles
             }

             // a final state muon at 0.1 - 1.2 GeV/c
             if (event.count_post(13) != 1) {
               return false; // reject events with no muon or more than one muon
             }

             auto muon = event.post_range(13).begin()->second;
             if (muon.P() < 0.1 || muon.P() > 1.2) {
               return false; // reject events with muon outside of the range
             }

             auto n_proton_in_range = std::ranges::distance(
                 event.post_range(2212) | std::views::values |
                 std::views::filter(
                     [](auto &&vec) { return vec.P() > 0.3 && vec.P() < 1.; }));
             if (n_proton_in_range != 1)
               return false; // require exactly one proton in the range

             auto selection_above_70_MeV =
                 std::views::values |
                 std::views::filter([](auto &&vec) { return vec.P() > 0.07; });
             auto n_charged_pion_above_70 =
                 // this is when I want the concat view...
                 std::ranges::distance(event.post_range(211) |
                                       selection_above_70_MeV) +
                 std::ranges::distance(event.post_range(-211) |
                                       selection_above_70_MeV);

             return n_charged_pion_above_70 ==
                    0; // require no charged pions above 70 MeV/c, final cut.
           },
           {"EventRecord"}, "Particle Type Cut")
          .Define("proton",
                  [](const NeutrinoEvent &event) {
                    return *(event.post_range(2212) | std::views::values |
                             std::views::filter([](auto &&vec) {
                               return vec.P() > 0.3 && vec.P() < 1.;
                             })).begin();
                  },
                  {"EventRecord"})
          .Define(
              "muon",
              [](const NeutrinoEvent &event) { return event.get_leading(13); },
              {"EventRecord"})
          .Define(
              "TKI",
              [](const TLorentzVector &neutrino_p, const TLorentzVector &muon_p,
                 const TLorentzVector &full_hadron) {
                auto ret =
                    getCommonTKI(40, 18, &neutrino_p, &muon_p, &full_hadron);
                if (ret.dalphat == -999) {
                  exit(1);
                }
                return ret;
              },
              {"InitNeutrino", "muon", "proton"})
          .Define("dalphat", [](const TKIVars &tki) { return tki.dalphat; },
                  {"TKI"});
  auto cutreport = signal.Report();

  auto dalphat =
      signal.Histo1D({"dalphat", "dalpha_{t};dalpha_{t} [deg];Events",
                      bin_edges.size() - 1, bin_edges.data()},
                     "dalphat", "weight");
  auto per_channel_dalphat =
      list_channel_channeldef | std::views::enumerate |
      std::views::transform([&](const auto &id_channel) {
        auto &&[id, channel] = id_channel;
        auto &[name, channel_def] = channel;
        return std::make_tuple(
            name,
            ROOT::RDF::RResultPtr<TH1>(
                signal.Filter(channel_def, {"channel"})
                    .Histo1D({("dalphat_" + name).c_str(),
                              "dalpha_{t};dalpha_{t} [deg];Events",
                              bin_edges.size() - 1, bin_edges.data()},
                             "dalphat", "weight")),
            id * 2);
      }) |
      std::ranges::to<std::vector>();

  for (auto &plot : per_channel_dalphat | std::views::elements<1>) {
    plot->Scale(1. / n_runs / 10, "width");
  }
  dalphat->Scale(1. / n_runs / 10, "width");

  auto &&[stack, legend] = build_stack_from_list(
      per_channel_dalphat, dalphat->Integral("width") * 0.015, {});
  auto exp_hist = make_exp_hist();

  double chi2{};
  for (size_t i = 0; i < bin_edges.size() - 1; ++i) {
    auto pred = dalphat->GetBinContent(i + 1);
    auto exp = cross_sections[i]/40.;
    auto unc = uncertainties[i]/40.;
    chi2 += std::pow((pred - exp) / unc, 2);
  }

  std::unique_ptr<TLatex> latex;
  if (!additional_text.empty()) {
    latex = std::make_unique<TLatex>(0.65, 0.5, additional_text.c_str());
    latex->SetNDC();
  }

  do_plot(
      {&exp_hist, dalphat, &stack, &legend, latex.get()}, "dalphat",
      "d#sigma/d#delta#it{#alpha}_{T} (10^{#minus 38} cm^{2}/degree/nucleon)",
      "#delta#it{#alpha}_{T} (degree)", {0.15, 0.6, 0.5, 0.9}, 0.,
      std::format("#chi^{{2}}/NDF = {:.0f}/7", chi2), "HISTC", 3e-3,
      {.top = 0.06, .bottom = 0.12});

  std::println("Chi2: {:.3f}", chi2);
  auto file = std::make_unique<TFile>("tki.root", "RECREATE");
  file->Add(dalphat.GetPtr());
  exp_hist.SetName("exphist");
  file->Add(&exp_hist);
  file->Write();
  file->Close();

  nlohmann::json j;
  j["chi2"] = chi2;

  std::ofstream out("tki.json");
  if (out.is_open()) {
    out << j.dump(4);
    out.close();
  } else {
    std::println("Could not open tki.json for writing");
  }
  std::println("Finished processing {} files", n_runs);

  cutreport->Print();
  std::cout << "Done." << std::endl;

  return 0;
}
