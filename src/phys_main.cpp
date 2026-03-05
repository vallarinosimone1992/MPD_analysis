#include "mpd_common.h"

#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLine.h>
#include <TParameter.h>
#include <TStyle.h>

#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Options {
  std::string input;
  std::string output;
  int avgEvents = 0;
  double pitch = 0.0;
  double stripCenter = -999.0;
};

void Usage(const char *prog)
{
  std::cout
      << "Usage: " << prog << " -i reco.root [options]\n\n"
      << "Options:\n"
      << "  -i  Input reco ROOT file (uClu, uEvt)\n"
      << "  -o  Output PDF file (default: output/<input>_phys.pdf)\n"
      << "  -A  Number of avg event profiles to plot from uEvt (0 disables)\n"
      << "  -P  Strip pitch in mm (override config/analysis.yaml)\n"
      << "  -c  Strip center (override config/analysis.yaml; negative = auto)\n"
      << "  -h  Show this help\n";
}

int RunPhys(const Options &opt)
{
  const TString inFile = opt.input.c_str();
  const TString outPdf = opt.output.c_str();

  LoadLegacyDictionaries();

  double pitch = opt.pitch;
  double stripCenter = opt.stripCenter;
  LoadConfigParams(pitch, stripCenter);
  if (pitch <= 0.0) {
    pitch = 1.0;
  }

  if (gSystem->AccessPathName(inFile)) {
    std::cerr << "[ERROR] Reco file not found: " << inFile << "\n";
    return 1;
  }

  TString outDir = gSystem->DirName(outPdf);
  if (outDir.Length() == 0) {
    outDir = ".";
  }
  gSystem->mkdir(outDir.Data(), true);

  TFile fIn(inFile, "READ");
  if (fIn.IsZombie()) {
    std::cerr << "[ERROR] Cannot open file: " << inFile << "\n";
    return 1;
  }

  double nsigma = 3.0;
  if (auto *p = dynamic_cast<TParameter<double>*>(fIn.Get("nsigma"))) {
    nsigma = p->GetVal();
  }

  TTreeReader tr("uClu", &fIn);
  if (!tr.GetTree()) {
    std::cerr << "[ERROR] Tree uClu not found in: " << inFile << "\n";
    return 1;
  }

  struct Range {
    bool set = false;
    double min = 0.0;
    double max = 0.0;
    void Update(double v) {
      if (!set) {
        min = max = v;
        set = true;
        return;
      }
      if (v < min) {
        min = v;
      }
      if (v > max) {
        max = v;
      }
    }
  };

  auto padRange = [](const Range &r, double defMin, double defMax) {
    if (!r.set) {
      return std::pair<double, double>(defMin, defMax);
    }
    double lo = r.min;
    double hi = r.max;
    if (lo == hi) {
      lo -= 1.0;
      hi += 1.0;
    }
    const double pad = 0.05 * (hi - lo);
    return std::pair<double, double>(lo - pad, hi + pad);
  };

  Range chargeAllR, charge0R, charge1R;
  Range posAllR, pos0R, pos1R;
  Range posStripR;

  TTreeReaderValue<Int_t> nclu(tr, "nclu");
  TTreeReaderValue<std::vector<short>> plane(tr, "plane");
  TTreeReaderValue<std::vector<short>> size(tr, "size");
  TTreeReaderValue<std::vector<double>> charge(tr, "charge");
  TTreeReaderValue<std::vector<double>> pos_mm(tr, "pos_mm");
  TTreeReaderValue<std::vector<double>> pos_strip(tr, "pos_strip");

  while (tr.Next()) {
    const size_t n = std::min({plane->size(), size->size(), charge->size(),
                               pos_mm->size(), pos_strip->size()});
    for (size_t i = 0; i < n; ++i) {
      const short pl = (*plane)[i];
      const double q = (*charge)[i];
      const double pos = (*pos_mm)[i];
      const double posS = (*pos_strip)[i];
      chargeAllR.Update(q);
      posAllR.Update(pos);
      posStripR.Update(posS);
      if (pl == 0) {
        charge0R.Update(q);
        pos0R.Update(pos);
      } else if (pl == 1) {
        charge1R.Update(q);
        pos1R.Update(pos);
      }
    }
  }

  const auto charge0Range = padRange(charge0R, 0.0, 2000.0);
  const auto charge1Range = padRange(charge1R, 0.0, 2000.0);
  const auto chargeAllRange = padRange(chargeAllR, 0.0, 2000.0);
  const auto pos0Range = padRange(pos0R, -50.0, 50.0);
  const auto pos1Range = padRange(pos1R, -50.0, 50.0);
  const auto posAllRange = padRange(posAllR, -50.0, 50.0);
  const auto posStripRange = padRange(posStripR, 0.0, 255.0);

  TTreeReader tr2("uClu", &fIn);
  TTreeReaderValue<Int_t> nclu2(tr2, "nclu");
  TTreeReaderValue<std::vector<short>> plane2(tr2, "plane");
  TTreeReaderValue<std::vector<short>> size2(tr2, "size");
  TTreeReaderValue<std::vector<double>> charge2(tr2, "charge");
  TTreeReaderValue<std::vector<double>> pos_mm2(tr2, "pos_mm");
  TTreeReaderValue<std::vector<double>> pos_strip2(tr2, "pos_strip");

  TH1D hNclu("hNclu", "Clusters per event;N;counts", 21, -0.5, 20.5);

  TH1D hSize0("hSize0", "Cluster size X (plane 0);size (strips);counts", 31, -0.5, 30.5);
  TH1D hSize1("hSize1", "Cluster size Y (plane 1);size (strips);counts", 31, -0.5, 30.5);
  TH1D hSizeAll("hSizeAll", "Cluster size all planes;size (strips);counts", 31, -0.5, 30.5);

  TH1D hCharge0("hCharge0", "Cluster charge X (plane 0);Q (ADC);counts", 200,
                charge0Range.first, charge0Range.second);
  TH1D hCharge1("hCharge1", "Cluster charge Y (plane 1);Q (ADC);counts", 200,
                charge1Range.first, charge1Range.second);
  TH1D hChargeAll("hChargeAll", "Cluster charge all planes;Q (ADC);counts", 200,
                  chargeAllRange.first, chargeAllRange.second);

  TH1D hPos0("hPos0", "Cluster position X (plane 0);pos (mm);counts", 200,
             pos0Range.first, pos0Range.second);
  TH1D hPos1("hPos1", "Cluster position Y (plane 1);pos (mm);counts", 200,
             pos1Range.first, pos1Range.second);
  TH1D hPosAll("hPosAll", "Cluster position all planes;pos (mm);counts", 200,
               posAllRange.first, posAllRange.second);

  TH1D h1PosStrip("h1PosStrip",
                  "Most significant cluster per event: position (strip);pos_strip;counts",
                  200, posStripRange.first, posStripRange.second);
  TH1D h1PosMm("h1PosMm",
               "Most significant cluster per event: position (mm);pos (mm);counts",
               200, posAllRange.first, posAllRange.second);
  TH1D h1Charge("h1Charge",
                "Most significant cluster per event: cluster charge (significance);Q (ADC);counts",
                200, chargeAllRange.first, chargeAllRange.second);

  TH2D hNcluVsSize("hNcluVsSize",
                   "Clusters per event vs cluster size;N clusters;cluster size (strips)",
                   21, -0.5, 20.5, 31, -0.5, 30.5);
  hNcluVsSize.SetCanExtend(TH1::kAllAxes);

  while (tr2.Next()) {
    hNclu.Fill(*nclu2);
    const size_t n = std::min({plane2->size(), size2->size(), charge2->size(),
                               pos_mm2->size(), pos_strip2->size()});
    for (size_t i = 0; i < n; ++i) {
      const short pl = (*plane2)[i];
      const int sz = (*size2)[i];
      const double q = (*charge2)[i];
      const double pos = (*pos_mm2)[i];

      hSizeAll.Fill(sz);
      hChargeAll.Fill(q);
      hPosAll.Fill(pos);
      hNcluVsSize.Fill(*nclu2, sz);

      if (pl == 0) {
        hSize0.Fill(sz);
        hCharge0.Fill(q);
        hPos0.Fill(pos);
      } else if (pl == 1) {
        hSize1.Fill(sz);
        hCharge1.Fill(q);
        hPos1.Fill(pos);
      }
    }

    if (n > 0) {
      size_t imax = 0;
      double qmax = (*charge2)[0];
      for (size_t i = 1; i < n; ++i) {
        if ((*charge2)[i] > qmax) {
          qmax = (*charge2)[i];
          imax = i;
        }
      }
      h1PosStrip.Fill((*pos_strip2)[imax]);
      h1PosMm.Fill((*pos_mm2)[imax]);
      h1Charge.Fill(qmax);
    }
  }

  gStyle->SetOptStat(0);
  TCanvas c("c", "c", 1600, 900);
  c.Print(Form("%s[", outPdf.Data()));

  c.Clear();
  c.Divide(2, 3);
  c.cd(1);
  hSize0.Draw("hist");
  c.cd(2);
  hSize1.Draw("hist");
  c.cd(3);
  gPad->SetLogy(true);
  hCharge0.Draw("hist");
  c.cd(4);
  gPad->SetLogy(true);
  hCharge1.Draw("hist");
  c.cd(5);
  gPad->SetLogy(true);
  hPos0.Draw("hist");
  c.cd(6);
  gPad->SetLogy(true);
  hPos1.Draw("hist");
  c.Print(outPdf.Data());

  c.Clear();
  c.Divide(2, 2);
  c.cd(1);
  gPad->SetLogy(false);
  hNclu.Draw("hist");
  c.cd(2);
  gPad->SetLogy(false);
  hSizeAll.Draw("hist");
  c.cd(3);
  gPad->SetLogy(true);
  hChargeAll.Draw("hist");
  c.cd(4);
  gPad->SetLogy(true);
  hPosAll.Draw("hist");
  c.Print(outPdf.Data());

  c.Clear();
  c.cd(1);
  gPad->SetLogy(false);
  hNcluVsSize.Draw("colz");
  c.Print(outPdf.Data());

  c.Clear();
  c.Divide(3, 1);
  c.cd(1);
  h1PosStrip.Draw("hist");
  c.cd(2);
  h1PosMm.Draw("hist");
  c.cd(3);
  h1Charge.Draw("hist");
  c.Print(outPdf.Data());

  std::vector<TH1D*> hPed(6, nullptr);
  for (int s = 0; s < 6; ++s) {
    hPed[s] = (TH1D *)fIn.Get(Form("hPed_%d", s));
  }

  if (opt.avgEvents > 0) {
    bool havePedContent = false;
    c.Clear();
    c.Divide(2, 3);
    for (int s = 0; s < 6; ++s) {
      c.cd(s + 1);
      if (hPed[s] != nullptr) {
        havePedContent = true;
        hPed[s]->SetLineColor(kBlue + 1);
        hPed[s]->SetLineWidth(1);
        hPed[s]->Draw("hist");
      }
    }
    if (havePedContent) {
      c.Print(outPdf.Data());
    }

    bool havePed = false;
    std::vector<std::unique_ptr<TH1D>> thrHists;
    thrHists.reserve(6);
    c.Clear();
    c.Divide(2, 3);
    for (int s = 0; s < 6; ++s) {
      c.cd(s + 1);
      if (hPed[s] != nullptr) {
        havePed = true;
        const int nBins = hPed[s]->GetNbinsX();
        const double xMin = hPed[s]->GetXaxis()->GetXmin();
        const double xMax = hPed[s]->GetXaxis()->GetXmax();
        auto hThr = std::make_unique<TH1D>(
            Form("hThrSample_%d", s),
            Form("Good strip threshold sample %d;Strip;%.1f#sigma (ADC)", s, nsigma),
            nBins, xMin, xMax);
        double maxVal = 0.0;
        for (int i = 1; i <= nBins; ++i) {
          const double thr = nsigma * hPed[s]->GetBinError(i);
          hThr->SetBinContent(i, thr);
          if (thr > maxVal) {
            maxVal = thr;
          }
        }
        if (maxVal <= 0.0) {
          maxVal = 1.0;
        }
        hThr->SetLineColor(kGreen + 2);
        hThr->SetLineWidth(1);
        hThr->SetMaximum(maxVal * 1.2);
        hThr->SetMinimum(0.0);
        hThr->Draw("hist");
        thrHists.push_back(std::move(hThr));
      }
    }
    if (havePed) {
      c.Print(outPdf.Data());
    }
  }

  TTree *tRwell = (TTree *)fIn.Get("uRwell");
  Int_t rEvt = 0;
  Int_t rNhit = 0;
  std::vector<short> *rStrip = nullptr;
  std::vector<short> *rSample = nullptr;
  std::vector<short> *rPlane = nullptr;
  std::vector<short> *rModule = nullptr;
  std::vector<double> *rAdc = nullptr;
  bool hasRwell = (tRwell != nullptr);
  if (hasRwell) {
    tRwell->SetBranchAddress("evt", &rEvt);
    tRwell->SetBranchAddress("nhit", &rNhit);
    tRwell->SetBranchAddress("strip", &rStrip);
    tRwell->SetBranchAddress("sample", &rSample);
    tRwell->SetBranchAddress("plane", &rPlane);
    tRwell->SetBranchAddress("module", &rModule);
    tRwell->SetBranchAddress("adc", &rAdc);
    tRwell->BuildIndex("evt");
  }

  if (opt.avgEvents > 0) {
    TTreeReader trEvt("uEvt", &fIn);
    if (trEvt.GetTree()) {
      TTreeReaderValue<Int_t> evt(trEvt, "evt");
      TTreeReaderValue<Short_t> planeE(trEvt, "plane");
      TTreeReaderValue<Short_t> moduleE(trEvt, "module");
      TTreeReaderValue<Int_t> sfirst(trEvt, "sfirst");
      TTreeReaderValue<Int_t> sdelta(trEvt, "sdelta");
      TTreeReaderValue<std::vector<double>> avgStrip(trEvt, "avg_strip");

      int shown = 0;
      while (trEvt.Next() && shown < opt.avgEvents) {
        const int sf = *sfirst;
        const int sd = *sdelta;
        const double usedCenter = (stripCenter < 0.0) ? (sf + 0.5 * (sd - 1)) : stripCenter;
        const double xMin = (sf - usedCenter) * pitch;
        const double xMax = (sf + sd - usedCenter) * pitch;

        const char *planeTag = (*planeE == 0) ? "X" : ((*planeE == 1) ? "Y" : "plane");
        if (hasRwell) {
          Long64_t entry = tRwell->GetEntryWithIndex(*evt);
          if (entry >= 0) {
            std::vector<std::vector<double>> ev(6, std::vector<double>(sd, 0.0));
            const size_t n = std::min({rStrip->size(), rSample->size(), rPlane->size(),
                                       rModule->size(), rAdc->size(),
                                       static_cast<size_t>(rNhit)});
            for (size_t i = 0; i < n; ++i) {
              if ((*rPlane)[i] != *planeE || (*rModule)[i] != *moduleE) {
                continue;
              }
              const int s = (*rSample)[i];
              if (s < 0 || s >= 6) {
                continue;
              }
              const int idx = (*rStrip)[i] - sf;
              if (idx < 0 || idx >= sd) {
                continue;
              }
              ev[s][idx] += (*rAdc)[i];
            }

            c.Clear();
            c.Divide(2, 3);
            std::vector<std::unique_ptr<TH1D>> hRaw, hSub, hThr;
            std::vector<std::unique_ptr<THStack>> stacks;
            std::unique_ptr<TLegend> leg;
            hRaw.reserve(6);
            hSub.reserve(6);
            hThr.reserve(6);
            stacks.reserve(6);
            for (int s = 0; s < 6; ++s) {
              const TString tName = Form("Evt %d %s plane %d mod %d sample %d;Strip;ADC",
                                          *evt, planeTag, *planeE, *moduleE, s);
              auto hS = std::make_unique<TH1D>(Form("hSub_%d_s%d", *evt, s),
                                               tName, sd, sf, sf + sd);
              auto hR = std::make_unique<TH1D>(Form("hRaw_%d_s%d", *evt, s),
                                               tName, sd, sf, sf + sd);
              auto hT = std::make_unique<TH1D>(Form("hThr_%d_s%d", *evt, s),
                                               "threshold", sd, sf, sf + sd);
              double maxVal = 0.0;
              for (int idx = 0; idx < sd; ++idx) {
                const double sub = ev[s][idx];
                double ped = 0.0;
                double thr = 0.0;
                if (hPed[s] != nullptr) {
                  const int bin = hPed[s]->FindBin(sf + idx);
                  ped = hPed[s]->GetBinContent(bin);
                  thr = nsigma * hPed[s]->GetBinError(bin);
                }
                const double raw = sub + ped;
                hS->SetBinContent(idx + 1, sub);
                hR->SetBinContent(idx + 1, raw);
                hT->SetBinContent(idx + 1, thr);
                maxVal = std::max(maxVal, std::max(raw, std::max(sub, thr)));
              }
              if (maxVal <= 0.0) {
                maxVal = 1.0;
              }
              hS->SetLineColor(kBlue + 1);
              hR->SetLineColor(kRed + 1);
              hT->SetLineColor(kGreen + 2);
              hT->SetLineWidth(1);
              hS->SetMaximum(maxVal * 1.2);
              hS->SetMinimum(std::min(0.0, -0.05 * maxVal));

              auto st = std::make_unique<THStack>(Form("st_%d_s%d", *evt, s), tName);
              st->Add(hR.get());
              st->Add(hS.get());

              c.cd(s + 1);
              st->Draw("nostack hist");
              hT->Draw("hist same");
              if (s == 0) {
                leg = std::make_unique<TLegend>(0.65, 0.7, 0.88, 0.88);
                leg->AddEntry(hR.get(), "Raw", "l");
                leg->AddEntry(hS.get(), "Ped-sub", "l");
                leg->AddEntry(hT.get(), Form("%.1f#sigma", nsigma), "l");
                leg->Draw();
              }

              hSub.push_back(std::move(hS));
              hRaw.push_back(std::move(hR));
              hThr.push_back(std::move(hT));
              stacks.push_back(std::move(st));
            }
            c.Print(outPdf.Data());
          }
        }

        c.Clear();
        TH1D hStrip("hAvgStrip",
                    Form("Evt %d %s plane %d mod %d;Strip;ADC", *evt, planeTag, *planeE, *moduleE),
                    sd, sf, sf + sd);
        TH1D hMm("hAvgMm",
                 Form("Evt %d %s plane %d mod %d;pos (mm);ADC", *evt, planeTag, *planeE, *moduleE),
                 sd, xMin, xMax);

        const size_t n = std::min(avgStrip->size(), static_cast<size_t>(sd));
        for (size_t i = 0; i < n; ++i) {
          hStrip.SetBinContent(static_cast<int>(i) + 1, (*avgStrip)[i]);
          hMm.SetBinContent(static_cast<int>(i) + 1, (*avgStrip)[i]);
        }

        c.Divide(2, 1);
        c.cd(1);
        hStrip.Draw("hist");
        c.cd(2);
        hMm.Draw("hist");
        c.Print(outPdf.Data());

        shown++;
      }
    }
  }

  c.Print(Form("%s]", outPdf.Data()));
  std::cout << "Output PDF: " << outPdf << "\n";
  return 0;
}

} // namespace

int main(int argc, char **argv)
{
  Options opt;

  const struct option longopts[] = {
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, "hi:o:A:P:c:", longopts, nullptr)) != -1) {
    switch (c) {
      case 'i':
        opt.input = optarg;
        break;
      case 'o':
        opt.output = optarg;
        break;
      case 'A':
        opt.avgEvents = std::stoi(optarg);
        break;
      case 'P':
        opt.pitch = std::stod(optarg);
        break;
      case 'c':
        opt.stripCenter = std::stod(optarg);
        break;
      case 'h':
        Usage(argv[0]);
        return 0;
      default:
        Usage(argv[0]);
        return 1;
    }
  }

  const TString suite = ResolveSuite();
  if (suite.Length() == 0) {
    std::cerr << "MPD_SUITE is not set. From MPD_dev root: "
                 "export MPD_SUITE=\"$(cd devel/MPD_analysis && pwd)\"\n";
    return 1;
  }

  const std::string recoDir = std::string(suite.Data()) + "/../RECO_DATA";
  const std::string outDir = std::string(suite.Data()) + "/../output";

  if (opt.input.empty()) {
    opt.input = recoDir + "/run_0254.dat_apv.root";
  }

  if (opt.output.empty()) {
    const std::string base = BaseName(opt.input);
    std::string pdfName;
    if (base.size() >= 5 && base.substr(base.size() - 5) == ".root") {
      pdfName = ReplaceExtension(base, "_phys.pdf");
    } else {
      pdfName = base + "_phys.pdf";
    }
    opt.output = outDir + "/" + pdfName;
  } else if (opt.output.find('/') == std::string::npos) {
    opt.output = outDir + "/" + opt.output;
  }

  gSystem->mkdir(outDir.c_str(), true);

  return RunPhys(opt);
}
