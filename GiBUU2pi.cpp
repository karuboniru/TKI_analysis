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
#include <fstream>
#include <functional>
#include <print>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

// const std::array<
//     std::tuple<std::string, std::function<bool(const NeutrinoEvent &)>>, 6>
//     pion_channels{
//         std::tuple<std::string, std::function<bool(const NeutrinoEvent &)>>{
//             "2#pi^{+}", // +2
//             [](const NeutrinoEvent &event) {
//               return event.count_post(211) == 2 &&
//                      event.count_post(-211) == 0 && event.count_post(111) ==
//                      0;
//             }},
//         {"1#pi^{0}", // +1
//          [](const NeutrinoEvent &event) {
//            return event.count_post(211) == 1 && event.count_post(-211) == 0
//            &&
//                   event.count_post(111) == 1;
//          }},
//         {"1#pi^{-}", // 0
//          [](const NeutrinoEvent &event) {
//            return event.count_post(211) == 1 && event.count_post(-211) == 1
//            &&
//                   event.count_post(111) == 0;
//          }},
//         {"2#pi^{0}", // 0
//          [](const NeutrinoEvent &event) {
//            return event.count_post(211) == 0 && event.count_post(-211) == 0
//            &&
//                   event.count_post(111) == 2;
//          }},
//         {"1#pi^{-}", // -1
//          [](const NeutrinoEvent &event) {
//            return event.count_post(211) == 0 && event.count_post(-211) == 1
//            &&
//                   event.count_post(111) == 1;
//          }},
//         {"2#pi^{-}", // -2
//          [](const NeutrinoEvent &event) {
//            return event.count_post(211) == 0 && event.count_post(-211) == 2
//            &&
//                   event.count_post(111) == 0;
//          }}};

const std::array<
    std::tuple<std::string, std::function<bool(const NeutrinoEvent &)>>, 1>
    pion_channels{
        std::tuple<std::string, std::function<bool(const NeutrinoEvent &)>>{
            " ", // +2
            [](const NeutrinoEvent &event) {
              return event.count_post(211) + event.count_post(-211) +
                         event.count_post(111) ==
                     2;
            }}};

const std::array<std::tuple<std::string, std::function<bool(int, double)>>, 3>
    interaction_cut{
        // std::tuple<std::string, std::function<bool(int)>>{
        //     "2#pi non-BG", [](int c) { return c != 37; }},
        std::tuple<std::string, std::function<bool(int, double)>>{
            "2#pi RES",
            [](int c, double) { return c >= 2 && c <= 31; }},
        // std::tuple<std::string, std::function<bool(int, double)>>{
        //     "2#pi SIS",
        //     [](int c, double W) { return c == 34 && W <= 3; }},
        std::tuple<std::string, std::function<bool(int, double)>>{
            "2#pi DIS", [](int c, double W) { return c == 34; }},
        std::tuple<std::string, std::function<bool(int, double)>>{
            "2#pi BG", [](int c, double) { return c == 37; }},
    };

const auto idlist =
    std::to_array<int>({2, 3, 4, 7, 10, 16, 31, 32, 33, 34, 37});

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
                     "Input files")("run-tag",
                                    po::value<std::string>()->default_value(""),
                                    "Additional text to legend")(
      "ymax", po::value<double>()->default_value(0), "Maximum y value")(
      "cross-section", po::value<std::string>()->default_value("xsec.txt"),
      "Cross section file")("help", "produce help message");
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
  auto ymax = vm["ymax"].as<double>();
  auto xsecfile =
      std::ofstream(vm["cross-section"].as<std::string>(), std::ios::trunc);

  TH1::AddDirectory(false);
  // IniColorCB2pibg();
  ROOT::EnableImplicitMT();
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
  auto weight_sum_channel = d.Define("tup",
                                     [](int channel, double weight) {
                                       return std::make_tuple(channel, weight);
                                     },
                                     {"channel", "weight"})
                                .Aggregate(
                                    [](std::map<int, double> &map,
                                       const std::tuple<int, double> &var) {
                                      auto &&[channel, weight] = var;
                                      map[channel] += weight;
                                    },
                                    [](std::map<int, double> map1,
                                       const std::map<int, double> &map2) {
                                      for (auto &&[channel, weight] : map2) {
                                        map1[channel] += weight;
                                      }
                                      return map1;
                                    },
                                    "tup", std::map<int, double>{});
  // auto weight_sum_by_channel =
  //     idlist | std::views::transform([&](int x) {
  //       return d.Filter([x](int channel) { return channel == x; },
  //       {"channel"})
  //           .Sum<double>("weight");
  //     }) |
  //     std::ranges::to<std::vector>();
  auto channel_merged = d.Aggregate(
      [](std::set<int> &c_set, int channel) { c_set.insert(channel); },
      [](std::set<int> c_set1, const std::set<int> &c_set2) {
        c_set1.insert(c_set2.begin(), c_set2.end());
        return c_set1;
      },
      "channel", std::set<int>{});
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
      {"EventRecord"}, "1#pi^{+} only");
  auto plot_single_pion = df_single_pion.Histo1D(
      {"1#pi", nbins, 0.8, 4.}, "W", "weight");
  auto plot_single_pion_res =
      df_single_pion
          .Filter([](int channel) { return channel >= 2 && channel <= 31; },
                  {"channel"})
          .Histo1D(
              {"1#pi RES", nbins, 0.8, 4.},
              "W", "weight");
  auto plot_single_pion_non_bg =
      df_single_pion
          .Filter(
              [](int channel) {
                return channel == 34 || (channel >= 2 && channel <= 31);
              },
              {"channel"})
          .Histo1D(
              {"1#kern[0.2]{#pi} DIS", "1#kern[0.2]{#pi} DIS", nbins, 0.8, 4.},
              "W", "weight");
  auto plot_list =
      std::views::cartesian_product(interaction_cut, pion_channels) |
      std::views::enumerate |
      std::views::transform([&d_pions, nbins](auto &&tup_id) {
        auto &&[idx, tup] = tup_id;
        auto &&[interaction_cut, channel_defintion] = tup;
        auto &&[channel_name, channel_def] = channel_defintion;
        auto &&[interaction_name, interaction_def] = interaction_cut;
        auto name = std::format("{:<14} {}", interaction_name, channel_name);
        return std::make_tuple(
            name,
            ROOT::RDF::RResultPtr<TH1>(
                d_pions.Filter(channel_def, {"EventRecord"})
                    .Filter(interaction_def, {"channel", "W"})
                    .Histo1D({name.c_str(), name.c_str(), nbins, 0.8, 4.}, "W",
                             "weight")),
            idx);
      }) |
      std::ranges::to<std::vector>() | std::views::filter([](auto &&tup) {
        // auto &&[name, hist, col] = tup;
        auto &hist = std::get<1>(tup);
        return hist->Integral("WIDTH") > 0;
      }) |
      std::ranges::to<std::vector>();
  auto plot_list_2pibg = plot_list | std::views::filter([](auto &&plot) {
                           //  auto &&[name, _unused1, _unused2] = plot;
                           auto &name = std::get<0>(plot);
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
  // plot_single_pion->Scale(1. / nruns / 50, "WIDTH");
  // plot_single_pion_bg->Scale(1. / nruns / 50, "WIDTH");
  std::ranges::for_each(
      std::to_array(
          {plot_single_pion, plot_single_pion_res, plot_single_pion_non_bg}),
      [&](auto &&hist) { hist->Scale(1. / nruns / 50, "WIDTH"); });

  std::cout << std::format("For Channel {:15}, the xsec is: {:3f} e-38 cm^2",
                           plot_single_pion->GetName(),
                           plot_single_pion->Integral("WIDTH"))
            << std::endl;

  plot_single_pion->SetLineColor(kBlack);
  plot_single_pion->SetLineWidth(4);
  plot_single_pion_non_bg->SetLineColor(kRed);
  plot_single_pion_non_bg->SetLineWidth(2);
  plot_single_pion_res->SetLineColor(kViolet);
  plot_single_pion_res->SetLineWidth(2);
  std::string x{"#it{W} (GeV)"},
      y{"d#sigma/d#it{W} (10^{#minus 38} cm^{2}/GeV/nucleon)"};
  auto [stack, stacklegend] = build_stack_from_list(plot_list, -1);
  TLegend legend(0.7, 0.6, 0.9, 0.9);

  legend.AddEntry(plot_single_pion.GetPtr());
  legend.AddEntry(plot_single_pion_non_bg.GetPtr());
  legend.AddEntry(plot_single_pion_res.GetPtr());
  for (auto &&[name, hist, _] : plot_list) {
    legend.AddEntry(hist.GetPtr());
  }
  // auto ymax = plot_single_pion->GetMaximum() * 1.1;
  do_plot({&stack, plot_single_pion, plot_single_pion_non_bg,
           plot_single_pion_res, &legend},
          "output_2pi", y, x, {0.65, 0.5, 0.9, 0.95}, 4.0, runtag, "HIST", ymax,
          {.top = 0.015, .bottom = 0.11});

  // auto [stack_2pibg, legend_2pibg] = build_stack_from_list(plot_list_2pibg,
  // -1); do_plot({&stack_2pibg, &legend_2pibg}, "output_2pibg", y, x,
  //         {0.7, 0.6, 0.9, 0.9}, 4.0, runtag, "HIST", 0);

  std::cout << "channel that I see for 2pi: ";
  for (auto &&channel : *channel_merged_2pi) {
    std::cout << channel << " ";
  }
  std::cout << std::endl;

  std::cout << "channel that I see for all: ";
  for (auto &&channel : *channel_merged) {
    std::cout << channel << " ";
  }
  std::cout << std::endl;

  // for (auto &&[channel, weight] :
  //      std::views::zip(idlist, weight_sum_by_channel)) {
  //   if (weight.GetValue() == 0) {
  //     continue;
  //   }
  //   std::println("Channel {:^3} has weight sum: {:>8}", channel,
  //                std::round(weight.GetValue()));
  //   std::println(xsecfile, "{} {}", channel, weight.GetValue());
  // }

  for (auto &&[channel, weight] : weight_sum_channel) {
    if (weight == 0) {
      continue;
    }
    std::println("Channel {:^3} has weight sum: {:>8.4f}", channel,
                 weight / nruns);
    std::println(xsecfile, "{} {}", channel, weight / nruns);
  }

  return 0;
}