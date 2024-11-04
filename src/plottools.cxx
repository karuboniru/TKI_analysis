#include "plottools.hxx"
#include <TNamed.h>
#include <algorithm>
#include <ranges>

THStack scale_stack(THStack &stack, double scale) {
  THStack stack_new;
  for (auto &&hist : *stack.GetHists()) {
    auto hist_new = (TH1 *)hist->Clone();
    hist_new->Scale(scale);
    stack_new.Add(hist_new);
  }
  return stack_new;
}
namespace {
static bool kset = false;
}

void IniColorCB() {
  if (kset) {
    printf("style::IniColorCB arleady set\n");
    return;
  } else {
    printf("style::IniColorCB creating new color\n");
  }

  // http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
  Int_t id = fgkColorBase + 1;
  // new TColor(id++, 0. / 255., 0. / 255., 0. / 255., "CB1_Black", 1.0);
  // 2pibg
  // new TColor(id++, 0. / 255., 73. / 255., 73. / 255., "CB2_Forest", 1.0);
  // new TColor(id++, 0. / 255., 146. / 255., 146. / 255., "CB3_Teal", 1.0);
  // new TColor(id++, 0. / 255., 219. / 255., 219. / 255., "CB3_Teal", 1.0);

  // new TColor(id++, 0. / 255., 146. / 255., 146. / 255., "CB3_Teal", 1.0);
  // QE
  new TColor(id++, 36. / 255., 255. / 255., 36. / 255., "CB14_DayGleen", 1.0);
  new TColor(id++, 0. / 255., 146. / 255., 146. / 255., "CB3_Teal", 1.0);
  // 2p2h
  new TColor(id++, 182. / 255., 109. / 255., 255. / 255., "CB8_Lilac", 1.0);
  // 2p2h-multipi: unused!
  new TColor(id++, 255. / 255., 109. / 255., 182. / 255., "dummy", 1.0);
  // RES
  // new TColor(id++, 146. / 255., 0. / 255., 0. / 255., "CB11_Maroon", 1.0);
  // new TColor(id++, 146. / 255., 73. / 255., 0. / 255., "CB12_Tan", 1.0);
  new TColor(id++, 182. / 255., 219. / 255., 255. / 255., "CB10_SpaceWolves",
             1.0);
  new TColor(id++, 0. / 255., 109. / 255., 219. / 255., "CB7_RoyalBlue", 1.0);
  // DIS
  new TColor(id++, 255. / 255., 182. / 255., 119. / 255., "CB5_BabyPink", 1.0);
  new TColor(id++, 146. / 255., 0. / 255., 0. / 255., "CB11_Maroon", 1.0);
  // 2pibg
  new TColor(id++, 0. / 255., 73. / 255., 73. / 255., "CB2_Forest", 1.0);
  new TColor(id++, 35. / 255., 141. / 255., 88. / 255., "CB9_BlueGrey", 1.0);
  // new TColor(id++, 182. / 255., 219. / 255., 255. / 255., "CB10_SpaceWolves",
  //  1.0);

  new TColor(id++, 219. / 255., 209. / 255., 0. / 255., "CB13_Orange", 1.0);
  // new TColor(id++, 36. / 255., 255. / 255., 36. / 255.,
  // "CB14_DayGleen", 1.0);
  new TColor(id++, 255. / 255., 255. / 255., 109. / 255., "CB15_SunFlower",
             1.0);

  kset = true;
}
// void IniColorCB2pibg() {
//   if (kset) {
//     printf("style::IniColorCB arleady set\n");
//     return;
//   } else {
//     printf("style::IniColorCB creating new color\n");
//   }

//   //
//   http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
//   Int_t id = fgkColorBase + 1;
//   // new TColor(id++, 0. / 255., 0. / 255., 0. / 255., "CB1_Black", 1.0);
//   // new TColor(id++, 0. / 255., 73. / 255., 73. / 255., "CB2_Forest", 1.0);
//   // new TColor(id++, 0. / 255., 146. / 255., 146. / 255., "CB3_Teal", 1.0);
//   // new TColor(id++, 0. / 255., 219. / 255., 219. / 255., "CB3_Teal", 1.0);
//   // 2pibg
//   new TColor(id++, 35. / 255., 141. / 255., 88. / 255., "CB9_BlueGrey", 1.0);
//   new TColor(id++, 0. / 255., 73. / 255., 73. / 255., "CB2_Forest", 1.0);
//   // 1pibg
//   new TColor(id++, 182. / 255., 219. / 255., 255. / 255., "CB10_SpaceWolves",
//              1.0);
//   new TColor(id++, 0. / 255., 109. / 255., 219. / 255.,
//   "CB7_RoyalBlue", 1.0);
//   // new TColor(id++, 109. / 255., 182. / 255., 255. / 255.,
//   // "CB9_BlueGrey", 1.0);
//   // new TColor(id++, 0. / 255., 146. / 255., 146. / 255., "CB3_Teal", 1.0);
//   // 2p2h
//   new TColor(id++, 255. / 255., 182. / 255., 119. / 255.,
//   "CB5_BabyPink", 1.0); new TColor(id++, 146. / 255., 0. / 255., 0. / 255.,
//   "CB11_Maroon", 1.0);
//   // QE
//   new TColor(id++, 36. / 255., 255. / 255., 36. / 255.,
//   "CB14_DayGleen", 1.0); new TColor(id++, 0. / 255., 146. / 255., 146. /
//   255., "CB3_Teal", 1.0);
//   // RES
//   // new TColor(id++, 146. / 255., 0. / 255., 0. / 255., "CB11_Maroon", 1.0);
//   // new TColor(id++, 146. / 255., 73. / 255., 0. / 255., "CB12_Tan", 1.0);

//   // DIS
//   new TColor(id++, 255. / 255., 182. / 255., 119. / 255.,
//   "CB5_BabyPink", 1.0); new TColor(id++, 146. / 255., 0. / 255., 0. / 255.,
//   "CB11_Maroon", 1.0);
//   // 2p2h
//   new TColor(id++, 182. / 255., 109. / 255., 255. / 255., "CB8_Lilac", 1.0);
//   // 2p2h-multipi: unused!
//   new TColor(id++, 255. / 255., 109. / 255., 182. / 255., "dummy", 1.0);

//   // new TColor(id++, 182. / 255., 219. / 255., 255. / 255.,
//   "CB10_SpaceWolves",
//   //  1.0);

//   // new TColor(id++, 219. / 255., 209. / 255., 0. / 255.,
//   "CB13_Orange", 1.0);
//   // new TColor(id++, 36. / 255., 255. / 255., 36. / 255.,
//   // "CB14_DayGleen", 1.0); new TColor(id++, 255. / 255., 255. / 255., 109. /
//   // 255., "CB15_SunFlower",
//   //            1.0);

//   kset = true;
// }

std::tuple<THStack, TLegend> build_stack_from_list(
    std::vector<std::tuple<std::string, ROOT::RDF::RResultPtr<TH1>, long>> list,
    double threshold, const legend_conf &conf) {
  IniColorCB();
  std::tuple<THStack, TLegend> tup;
  auto &&[stack, legend] = tup;
  for (auto &&[name, plot, color_diff] : list) {
    // auto &&[name, plot] = plot_entry;
    plot->SetLineColor(fgkColorBase + 1 + color_diff);
    plot->SetFillColor(fgkColorBase + 1 + color_diff);

    stack.Add(plot.GetPtr());
  }
  for (auto &&[name, plot, _] : list) {
    if ((std::ranges::any_of(conf.force_include,
                             [&](auto &&str) { return name.contains(str); }) ||
         plot->Integral("WIDTH") > threshold) &&
        (std::ranges::none_of(conf.force_exclude, [&](auto &&str) {
          return name.contains(str);
        }))) {
      legend.AddEntry(plot.GetPtr(), name.c_str(), "f");
    }
  }
  return tup;
}

// std::tuple<THStack, TLegend> build_stack_from_list(
//     std::vector<
//         std::tuple<std::string, ROOT::RDF::RResultPtr<TH1>, long, short>>
//         list,
//     double threshold, const legend_conf &conf) {
//   IniColorCB();
//   std::tuple<THStack, TLegend> tup;
//   auto &&[stack, legend] = tup;
//   for (auto &&[name, plot, color_diff, fill_attr] : list) {
//     // auto &&[name, plot] = plot_entry;
//     plot->SetLineColor(fgkColorBase + 1 + color_diff);
//     plot->SetFillColor(fgkColorBase + 1 + color_diff);
//     plot->SetFillStyle(fill_attr);
//     // plot->SetFillStyle(1001);

//     stack.Add(plot.GetPtr());
//     // if (plot->Integral("WIDTH") > threshold)
//     //   legend.AddEntry(plot.GetPtr(), name.c_str(), "f");
//     // else
//     //   skipped++;
//   }
//   for (auto &&[name, plot, _, _u] : list) {
//     if ((std::ranges::any_of(conf.force_include,
//                              [&](auto &&str) { return name.contains(str); })
//                              ||
//          plot->Integral("WIDTH") > threshold) &&
//         (std::ranges::none_of(conf.force_exclude, [&](auto &&str) {
//           return name.contains(str);
//         }))) {
//       legend.AddEntry(plot.GetPtr(), name.c_str(), "f");
//     }
//   }
//   return tup;
// }

std::tuple<THStack, TLegend> build_stack_from_list(
    std::vector<std::tuple<std::string, std::unique_ptr<TH1D>, long>> &list,
    double threshold, const legend_conf &conf) {
  IniColorCB();
  std::tuple<THStack, TLegend> tup;
  auto &&[stack, legend] = tup;
  for (auto &&[name, plot, color_diff] : list) {
    // auto &&[name, plot] = plot_entry;
    plot->SetLineColor(fgkColorBase + 1 + color_diff);
    plot->SetFillColor(fgkColorBase + 1 + color_diff);
    // plot->SetFillStyle(1001);

    stack.Add(plot.get());
    // if (plot->Integral("WIDTH") > threshold)
    //   legend.AddEntry(plot.GetPtr(), name.c_str(), "f");
    // else
    //   skipped++;
  }
  for (auto &&[name, plot, _] : list) {
    if ((std::ranges::any_of(conf.force_include,
                             [&](auto &&str) { return name.contains(str); }) ||
         plot->Integral("WIDTH") > threshold) &&
        (std::ranges::none_of(conf.force_exclude, [&](auto &&str) {
          return name.contains(str);
        }))) {
      legend.AddEntry(plot.get(), name.c_str(), "f");
    }
  }
  return tup;
}

std::tuple<THStack, TLegend> build_stack_from_list(
    std::vector<std::tuple<std::string, ROOT::RDF::RResultPtr<TH1>>> list,
    double threshold, const legend_conf &conf) {
  IniColorCB();
  std::tuple<THStack, TLegend> tup;
  auto &&[stack, legend] = tup;
  // size_t skipped = 0;
  for (auto &&[id, obj] : list | std::ranges::views::enumerate) {
    auto &&[name, plot] = obj;
    const auto color_diff = id + 1;
    // auto &&[name, plot] = plot_entry;
    plot->SetLineColor(fgkColorBase + 1 + color_diff);
    plot->SetFillColor(fgkColorBase + 1 + color_diff);
    // plot->SetFillStyle(1001);

    stack.Add(plot.GetPtr());
    // if (plot->Integral("WIDTH") > threshold)
    //   legend.AddEntry(plot.GetPtr(), name.c_str(), "f");
    // else
    //   skipped++;
  }
  for (auto &&[name, plot] : list) {
    if ((std::ranges::any_of(conf.force_include,
                             [&](auto &&str) { return name.contains(str); }) ||
         plot->Integral("WIDTH") > threshold) &&
        (std::ranges::none_of(conf.force_exclude, [&](auto &&str) {
          return name.contains(str);
        }))) {
      legend.AddEntry(plot.GetPtr(), name.c_str(), "f");
    }
  }
  return tup;
}

void do_plot(std::vector<plot_ptr_t> plot_ptrs_list,
             const std::string &filename, std::string_view ytitle,
             std::string_view xtitle, std::array<double, 4> legend_pos,
             double xmax, std::string legend_head, std::string histopt,
             double ymax, override_margin m) {
  auto plot = [](TNamed &obj, std::string_view option) {
    std::cout << std::format("Drawing obj {:5} with option {:5}", obj.GetName(),
                             option)
              << std::endl;
    obj.Draw(option.data());
  };

  auto c = getCanvas();
  PadSetup(c.get());
  if (m.left)
    c->SetLeftMargin(m.left);
  if (m.right)
    c->SetRightMargin(m.right);
  if (m.top)
    c->SetTopMargin(m.top);
  if (m.bottom)
    c->SetBottomMargin(m.bottom);
  std::unique_ptr<TLegend> leg_copy;
  for (auto [id, obj] : plot_ptrs_list | std::views::enumerate) {
    // std::string hist_draw_opt = id == 0 ? "E0 X1" : "HISTC same";
    std::string hist_draw_opt = id == 0 ? histopt : (histopt + " same");
    std::visit(
        [&](auto &&arg) {
          if (!arg) // ignore nullptrs
            return;
          if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, TH1 *> ||
                        std::is_same_v<std::decay_t<decltype(arg)>,
                                       ROOT::RDF::RResultPtr<TH1>> ||
                        std::is_same_v<std::decay_t<decltype(arg)>,
                                       ROOT::RDF::RResultPtr<TH1D>>) {
            if (id == 0)
              ResetStyle(arg);
            arg->SetStats(0);
            arg->GetXaxis()->SetTitle(xtitle.data());
            arg->GetYaxis()->SetTitle(ytitle.data());
            arg->SetMinimum(0);
            if (xmax > 0.) {
              arg->GetXaxis()->SetRangeUser(arg->GetXaxis()->GetXmin(), xmax);
            }
            if (ymax) {
              arg->SetMaximum(ymax);
            }
            if (id == 0) {
              arg->SetFillColor(kBlack);
              arg->SetLineColor(kBlack);
              arg->SetLineWidth(2);
              arg->SetMarkerStyle(20);
              if (!ymax)
                arg->SetMaximum(arg->GetMaximum() * 1.5);
            }
            if (std::is_same_v<std::decay_t<decltype(arg)>, TH1 *> && id == 0) {
              if (id == 0)
                // arg->Draw("E0 X1");
                plot(*arg, "E0 X1");
              else
                // arg->Draw("E0 X1 SAME");
                plot(*arg, "E0 X1 SAME");
            } else
              // arg->Draw(hist_draw_opt.c_str());
              plot(*arg, hist_draw_opt);
          } else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>,
                                              TLegend *>) {
            leg_copy = std::make_unique<TLegend>(*arg);
            ResetStyle(leg_copy);
            leg_copy->SetX1(legend_pos[0]);
            leg_copy->SetY1(legend_pos[1]);
            leg_copy->SetX2(legend_pos[2]);
            leg_copy->SetY2(legend_pos[3]);
            if (std::holds_alternative<TH1 *>(plot_ptrs_list[0])) {
              auto datahist = std::get<TH1 *>(plot_ptrs_list[0]);
              leg_copy->AddEntry(datahist, datahist->GetTitle(), "lep");
            }
            if (!legend_head.empty()) {
              leg_copy->SetHeader(legend_head.data());
            }
            leg_copy->Draw();
            // plot(*leg_copy, "");
          } else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>,
                                              THStack *>) {
            TObjArray *stack = arg->GetStack();
            auto size = stack->GetSize();
            auto last = dynamic_cast<TH1 *>(stack->Last());
            if (id == 0) {
              last->Draw();
              if (xmax > 0.) {
                last->GetXaxis()->SetRangeUser(last->GetXaxis()->GetXmin(),
                                               xmax);
              }
              last->GetXaxis()->SetTitle(xtitle.data());
              last->GetYaxis()->SetTitle(ytitle.data());
              if (ymax != 0.) {
                last->SetMaximum(ymax);
              }
              // last->Draw(histopt.c_str());
              plot(*last, histopt);
            }
            for (int i = size - 1; i >= 0; i--) {
              auto stack_obj = dynamic_cast<TH1 *>(stack->At(i));
              stack_obj->GetXaxis()->SetTitle(xtitle.data());
              stack_obj->GetYaxis()->SetTitle(ytitle.data());
              if (ymax != 0.) {
                stack_obj->SetMaximum(ymax);
              } else {
                stack_obj->SetMaximum(stack_obj->GetMaximum() * 1.3);
              }
              ResetStyle(stack_obj);
              stack_obj->SetFillColor(stack_obj->GetLineColor());
              stack_obj->SetFillStyle(1001);
              // stack_obj->Draw((histopt + " SAME").c_str());
              plot(*stack_obj, (histopt + " SAME"));
            }
            last->Draw("AXIS SAME");
          } else if constexpr (std::is_same_v<std::decay_t<decltype(arg)>,
                                              TLatex *>) {
            if (legend_pos[0] < 0.3) {
              arg->SetX(legend_pos[2] - 0.05);
              arg->SetY(legend_pos[3] - 0.05);
            } else {
              arg->SetX(legend_pos[0] - 0.45);
              arg->SetY(legend_pos[3] - 0.05);
            }
            // arg->Draw();
            plot(*arg, "");
          } else {
            throw std::runtime_error("Unknown type in plot_ptrs_list");
          }
        },
        obj);
  }
  std::visit(
      [&](auto &&obj) {
        if constexpr (std::is_same_v<std::decay_t<decltype(obj)>, TH1 *>)
          obj->Draw("E0 X1 SAME");
      },
      plot_ptrs_list[0]);
  c->SaveAs((filename + ".pdf").c_str());
  c->SaveAs((filename + ".eps").c_str());
  c->SaveAs((filename + ".svg").c_str());
}

void IniColorCBminerva2pibg() {
  if (kset) {
    printf("style::IniColorCB arleady set\n");
    return;
  } else {
    printf("style::IniColorCB creating new color\n");
  }

  // http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
  Int_t id = fgkColorBase + 1;
  // new TColor(id++, 0. / 255., 0. / 255., 0. / 255., "CB1_Black", 1.0);
  // 2pibg
  // new TColor(id++, 0. / 255., 73. / 255., 73. / 255., "CB2_Forest", 1.0);
  // new TColor(id++, 0. / 255., 146. / 255., 146. / 255., "CB3_Teal", 1.0);
  // new TColor(id++, 0. / 255., 219. / 255., 219. / 255., "CB3_Teal", 1.0);

  // 1

  // new TColor(id++, 146. / 255., 0. / 255., 0. / 255., "CB11_Maroon", 1.0);
  // QE
  new TColor(id++, 36. / 255., 255. / 255., 36. / 255., "CB14_DayGleen", 1.0);
  new TColor(id++, 0. / 255., 146. / 255., 146. / 255., "CB3_Teal", 1.0);
  // 2p2h
  new TColor(id++, 182. / 255., 109. / 255., 255. / 255., "CB8_Lilac", 1.0);
  // 2p2h-multipi: unused!
  new TColor(id++, 255. / 255., 109. / 255., 182. / 255., "dummy", 1.0);

  // RES
  new TColor(id++, 182. / 255., 219. / 255., 255. / 255., "CB10_SpaceWolves",
             1.0);
  new TColor(id++, 0. / 255., 109. / 255., 219. / 255., "CB7_RoyalBlue", 1.0);

  // DIS
  new TColor(id++, 255. / 255., 182. / 255., 119. / 255., "CB5_BabyPink", 1.0);
  new TColor(id++, 146. / 255., 0. / 255., 0. / 255., "CB11_Maroon", 1.0);
  // 1pibg
  new TColor(id++, 219. / 255., 209. / 255., 0. / 255., "CB8_Lilac_dark", 1.0);
  new TColor(id++, 182. / 255., 109. / 255., 255. / 255., "CB8_Lilac", 1.0);
  // 2pibg
  new TColor(id++, 0. / 255., 73. / 255., 73. / 255., "CB2_Forest", 1.0);
  new TColor(id++, 35. / 255., 141. / 255., 88. / 255., "CB9_BlueGrey", 1.0);
  kset = true;
}
