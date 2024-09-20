#include "EvtTracker2event.h"
#include "MINERvA_tki_cut.h"
#include "plottools.hxx"
#include "tkievent.h"
#include "tkigeneral.h"

#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <THStack.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMatrixDfwd.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>
#include <boost/program_options.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

const std::array<
    std::tuple<std::string, std::function<bool(const NeutrinoEvent &)>>, 6>
    pion_channels{
        std::tuple<std::string, std::function<bool(const NeutrinoEvent &)>>{
            "2#kern[0.2]{#pi}^{+}",
            [](const NeutrinoEvent &event) {
              return event.count_post(211) == 2 &&
                     event.count_post(-211) == 0 && event.count_post(111) == 0;
            }},
        {"1#kern[0.2]{#pi}^{+}1#kern[0.2]{#pi}^{0}",
         [](const NeutrinoEvent &event) {
           return event.count_post(211) == 1 && event.count_post(-211) == 0 &&
                  event.count_post(111) == 1;
         }},
        {"1#kern[0.2]{#pi}^{+}1#kern[0.2]{#pi}^{-}",
         [](const NeutrinoEvent &event) {
           return event.count_post(211) == 1 && event.count_post(-211) == 1 &&
                  event.count_post(111) == 0;
         }},
        {"2#kern[0.2]{#pi}^{0}",
         [](const NeutrinoEvent &event) {
           return event.count_post(211) == 0 && event.count_post(-211) == 0 &&
                  event.count_post(111) == 2;
         }},
        {"1#kern[0.2]{#pi}^{0}1#kern[0.2]{#pi}^{-}",
         [](const NeutrinoEvent &event) {
           return event.count_post(211) == 0 && event.count_post(-211) == 1 &&
                  event.count_post(111) == 1;
         }},
        {"2#kern[0.2]{#pi}^{-}", [](const NeutrinoEvent &event) {
           return event.count_post(211) == 0 && event.count_post(-211) == 2 &&
                  event.count_post(111) == 0;
         }}};

const std::array<std::tuple<std::string, std::function<bool(int)>>, 2>
    interaction_cut{
        std::tuple<std::string, std::function<bool(int)>>{
            "2#kern[0.2]{#pi}BG", [](int c) { return c == 37; }},
        std::tuple<std::string, std::function<bool(int)>>{
            "non-2#kern[0.2]{#pi}BG", [](int c) { return c != 37; }}};

ROOT::RDF::RNode vars_define(ROOT::RDF::RNode df) {
  return df
      .Define("Q2",
              [](const TLorentzVector &InitNeutrino,
                 const TLorentzVector &PrimaryLepton) {
                auto q0 = InitNeutrino - PrimaryLepton;
                return -q0.Mag2();
              },
              {"InitNeutrino", "PrimaryLepton"})
      .Define("W",
              [](const TLorentzVector &InitNeutrino,
                 const TLorentzVector &PrimaryLepton,
                 const TLorentzVector &InitNucleon) {
                auto had_system = InitNucleon + InitNeutrino - PrimaryLepton;
                return had_system.M();
              },
              {"InitNeutrino", "PrimaryLepton", "InitNucleon"})
      .Define("xBj",
              [](const TLorentzVector &InitNeutrino,
                 const TLorentzVector &PrimaryLepton,
                 const TLorentzVector &InitNucleon) {
                auto q = InitNeutrino - PrimaryLepton;
                return -q.Mag2() / (2 * InitNucleon.Dot(q));
              },
              {"InitNeutrino", "PrimaryLepton", "InitNucleon"});
}

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()("input-files", po::value<std::vector<std::string>>(),
                     "Input files")(
      "run-tag", po::value<std::string>()->default_value(""),
      "Additional text to legend")("help", "produce help message");
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
  auto nruns = files.size();
  auto &&runtag = vm["run-tag"].as<std::string>();

  TH1::AddDirectory(false);
  IniColorCB2pibg();
  ROOT::EnableImplicitMT();
  // std::vector<std::string> files{};
  // std::string runtag = argv[1];
  // // std::string leg_head =
  // for (int i = 2; i < argc; i++) {
  //   files.emplace_back(argv[i]);
  // }
  // auto nruns = files.size();
  ROOT::RDataFrame input("out_tree", files);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 30, 0)
  ROOT::RDF::Experimental::AddProgressBar(input);
#endif
  auto d = vars_define(TrackerPrepare(input));
  auto count = input.Count();
  auto d_pions = d.Filter(
      [](const NeutrinoEvent &event) {
        return std::ranges::all_of(event.get_ids_post(), [](int id) {
          return id == 211 || id == -211 || id == 111 || id == 2112 ||
                 id == 2212 || abs(id) == 13;
        });
      },
      {"EventRecord"}, "No other particles");
  auto channel_merged_2pi = d_pions.Aggregate(
      [](std::set<int> &c_set, int channel) { c_set.insert(channel); },
      [](std::set<int> c_set1, const std::set<int> &c_set2) {
        c_set1.insert(c_set2.begin(), c_set2.end());
        return c_set1;
      },
      "channel", std::set<int>{});
  auto d_2pibg =
      d.Filter([](int channel) { return channel == 37; }, {"channel"});
  constexpr size_t nbins = 100;
  auto df_single_pion = d_pions.Filter(
      [](const NeutrinoEvent &event) {
        return event.count_post(211) + event.count_post(-211) +
                   event.count_post(111) ==
               1;
      },
      {"EventRecord"}, "1#kern[0.2]{#pi}^{+} only");
  auto plot_single_pion = df_single_pion.Histo1D(
      {"1#kern[0.2]{#pi}", "1#kern[0.2]{#pi}", nbins, 0.8, 4.}, "W", "weight");
  auto plot_single_pion_bg =
      df_single_pion
          .Filter([](int channel) { return channel == 32 || channel == 33; },
                  {"channel"})
          .Histo1D(
              {"1#kern[0.2]{#pi} BG", "1#kern[0.2]{#pi} BG", nbins, 0.8, 4.},
              "W", "weight");
  auto plot_list =
      std::views::cartesian_product(interaction_cut, pion_channels) |
      std::views::transform([&d_pions, nbins](auto &&tup) {
        auto &&[interaction_cut, channel_defintion] = tup;
        auto &&[channel_name, channel_def] = channel_defintion;
        auto &&[interaction_name, interaction_def] = interaction_cut;
        auto name = std::format("{:<14} {}", interaction_name, channel_name);
        return std::make_tuple(
            name, ROOT::RDF::RResultPtr<TH1>(
                      d_pions.Filter(channel_def, {"EventRecord"})
                          .Filter(interaction_def, {"channel"})
                          .Histo1D({name.c_str(), name.c_str(), nbins, 0.8, 4.},
                                   "W", "weight")));
      }) |
      std::ranges::to<std::vector>() | std::views::filter([](auto &&tup) {
        auto &&[name, hist] = tup;
        return hist->Integral("WIDTH") > 0;
      }) |
      std::ranges::to<std::vector>();
  auto plot_list_2pibg = plot_list | std::views::filter([](auto &&plot) {
                           auto &&[name, _] = plot;
                           return !name.starts_with("non");
                         }) |
                         std::ranges::to<std::vector>();

  for (auto &&hist : plot_list | std::views::values) {
    hist->Scale(1. / nruns / 50, "WIDTH");
    hist->SetMinimum(0);
    std::cout << std::format("For Channel {:15}, the xsec is: {:3f} e-38 cm^2",
                             hist->GetName(), hist->Integral("WIDTH"))
              << std::endl;
  }
  plot_single_pion->Scale(1. / nruns / 50, "WIDTH");
  plot_single_pion_bg->Scale(1. / nruns / 50, "WIDTH");

  std::cout << std::format("For Channel {:15}, the xsec is: {:3f} e-38 cm^2",
                           plot_single_pion->GetName(),
                           plot_single_pion->Integral("WIDTH"))
            << std::endl;
  plot_single_pion->SetLineColor(kBlack);
  plot_single_pion->SetLineWidth(3);
  plot_single_pion_bg->SetLineColor(kRed);
  plot_single_pion_bg->SetLineWidth(2);
  std::string x{"#it{W} (GeV)"},
      y{"d#sigma/d#it{W} (10^{#minus 38} cm^{2}/GeV/nucleon)"};
  auto [stack, legend] = build_stack_from_list(plot_list, -1);

  legend.AddEntry(plot_single_pion.GetPtr());
  legend.AddEntry(plot_single_pion_bg.GetPtr());
  auto ymax = plot_single_pion->GetMaximum() * 1.1;
  do_plot({&stack, plot_single_pion, plot_single_pion_bg, &legend},
          "output_2pi", y, x, {0.7, 0.6, 0.9, 0.9}, 4.0, runtag, "HIST", ymax);
  auto [stack_2pibg, legend_2pibg] = build_stack_from_list(plot_list_2pibg, -1);
  do_plot({&stack_2pibg, &legend_2pibg}, "output_2pibg", y, x,
          {0.7, 0.6, 0.9, 0.9}, 4.0, runtag, "HIST", 0);

  std::cout << "channel that I see for 2pi: ";
  for (auto &&channel : *channel_merged_2pi) {
    std::cout << channel << " ";
  }
  std::cout << std::endl;
  return 0;
}