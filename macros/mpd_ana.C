// Monolithic MPD analysis macro (RDataFrame-based)
// Usage:
//   root -l -b -q '$MPD_SUITE/macros/mpd_ana.C(1234, 1200)'
//   root -l -b -q '$MPD_SUITE/macros/mpd_ana.C(1234, 1200, "$MPD_SUITE/../DATA/out")'
//   root -l -b -q '$MPD_SUITE/macros/mpd_ana.C("$MPD_SUITE/../DATA/out/run_0297.dat_apv.root", "$MPD_SUITE/../DATA/out/run_0295.dat_apv.root", "$MPD_SUITE/../output/run_0297_mpdAna.root")'
//
// Notes:
// - This is a first-approximation, inspired by geAna/siAna.
// - Event selection counts hits above nsigma*pedError (approximation of geAna goodEvent).
// - Output TTree "uRwell" contains ped-subtracted ADCs.

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using ROOT::VecOps::RVec;

namespace {

TString ResolveSuite()
{
  const char *env = gSystem->Getenv("MPD_SUITE");
  if (env != nullptr && TString(env).Length() > 0) {
    TString s(env);
    if (s.EndsWith("/")) {
      s.Chop();
    }
    return s;
  }
  return "";
}

void LoadConfigParams(double &pitch, double &stripCenter)
{
  const TString suite = ResolveSuite();
  if (suite.Length() == 0) {
    return;
  }
  const TString cfgPath = suite + "/config/analysis.yaml";
  std::ifstream in(cfgPath.Data());
  if (!in) {
    return;
  }
  std::string line;
  double cfgPitch = 0.0;
  double cfgCenter = std::numeric_limits<double>::quiet_NaN();
  while (std::getline(in, line)) {
    auto hash = line.find('#');
    if (hash != std::string::npos) {
      line.erase(hash);
    }
    auto trim = [](std::string &s) {
      const char *ws = " \t\r\n";
      s.erase(0, s.find_first_not_of(ws));
      s.erase(s.find_last_not_of(ws) + 1);
    };
    trim(line);
    if (line.empty()) {
      continue;
    }
    auto parseValue = [&](const std::string &key, double &out) {
      auto pos = line.find(key);
      if (pos == std::string::npos) {
        return;
      }
      auto colon = line.find(':', pos + key.size());
      if (colon == std::string::npos) {
        return;
      }
      std::string val = line.substr(colon + 1);
      trim(val);
      if (!val.empty()) {
        out = std::stod(val);
      }
    };
    parseValue("pitch_mm", cfgPitch);
    parseValue("strip_center", cfgCenter);
  }
  if (pitch <= 0.0 && cfgPitch > 0.0) {
    pitch = cfgPitch;
  }
  if (stripCenter < -900.0 && !std::isnan(cfgCenter)) {
    stripCenter = cfgCenter;
  }
}

std::vector<ROOT::RDF::RResultPtr<TProfile>> BuildPedestals(ROOT::RDF::RNode df,
                                                           int nsamples,
                                                           int sfirst,
                                                           int sdelta)
{
  std::vector<ROOT::RDF::RResultPtr<TProfile>> profs;
  profs.reserve(nsamples);

  const int slast = sfirst + sdelta;

  for (int i = 0; i < nsamples; ++i) {
    const TString maskName = Form("mask_s%d", i);
    const TString stripName = Form("strip_s%d", i);
    const TString adcName = Form("adc_s%d", i);

    auto df_i = df.Define(maskName.Data(),
                          Form("(ft.sample==%d) && (ft.strip>=%d) && (ft.strip<%d)",
                               i, sfirst, slast))
                    .Define(stripName.Data(), Form("ft.strip[%s]", maskName.Data()))
                    .Define(adcName.Data(), Form("ft.adc[%s]", maskName.Data()));

    profs.push_back(df_i.Profile1D(
        {Form("hPed_%d", i), Form("Pedestal sample %d;Strip;ADC", i),
         sdelta, static_cast<double>(sfirst), static_cast<double>(slast)},
        stripName.Data(),
        adcName.Data()));
  }

  return profs;
}

void LoadLegacyDictionaries()
{
  // Optional: load legacy dictionaries to silence MPDdb/APVdb/DAQpar errors.
  const TString suite = ResolveSuite();
  if (suite.Length() == 0) {
    return;
  }
  const TString legacyDir = suite + "/src/legacy";

  const TString dqdb = legacyDir + "/DQdb.h";

  if (!gSystem->AccessPathName(dqdb)) {
    gROOT->ProcessLine(Form(".L %s+", dqdb.Data()));
  }
}

int RunAnalysis(const TString& sigFile,
                const TString& pedFile,
                const TString& outFile,
                double nsigma,
                int sfirst,
                int sdelta,
                int minGood,
                double pitch,
                double stripCenter,
                int plotEvents,
                bool makePlots,
                const TString& logFile)
{
  LoadLegacyDictionaries();
  ROOT::EnableImplicitMT(0);

  LoadConfigParams(pitch, stripCenter);

  if (gSystem->AccessPathName(sigFile)) {
    std::cerr << "[ERROR] Signal file not found: " << sigFile << "\n";
    return 1;
  }
  if (gSystem->AccessPathName(pedFile)) {
    std::cerr << "[ERROR] Pedestal file not found: " << pedFile << "\n";
    return 1;
  }

  const int nsamples = 6;

  const double stripCenterInput = stripCenter;
  if (stripCenter < 0) {
    ROOT::RDataFrame dfTmp("fT", sigFile.Data());
    auto minStrip = *dfTmp.Define("min_strip", "ROOT::VecOps::Min(ft.strip)").Min("min_strip");
    auto maxStrip = *dfTmp.Define("max_strip", "ROOT::VecOps::Max(ft.strip)").Max("max_strip");
    stripCenter = 0.5 * (static_cast<double>(minStrip) + static_cast<double>(maxStrip));
  }

  ROOT::RDataFrame dfPed("fT", pedFile.Data());
  auto pedProfiles = BuildPedestals(dfPed, nsamples, sfirst, sdelta);
  for (auto& p : pedProfiles) {
    p.GetValue();
  }

  std::vector<TProfile*> pedPtr(nsamples, nullptr);
  for (int i = 0; i < nsamples; ++i) {
    pedPtr[i] = pedProfiles[i].GetPtr();
  }

  std::vector<std::vector<double>> pedMean(nsamples, std::vector<double>(sdelta, 0.0));
  std::vector<std::vector<double>> pedErr(nsamples, std::vector<double>(sdelta, 0.0));
  for (int is = 0; is < nsamples; ++is) {
    if (pedPtr[is] == nullptr) {
      continue;
    }
    for (int idx = 0; idx < sdelta; ++idx) {
      const int strip = sfirst + idx;
      const int bin = pedPtr[is]->FindBin(strip);
      pedMean[is][idx] = pedPtr[is]->GetBinContent(bin);
      pedErr[is][idx] = pedPtr[is]->GetBinError(bin);
    }
  }

  auto take_n_short = [](const RVec<short>& v, int nhit) {
    const size_t n = std::min<size_t>(v.size(), static_cast<size_t>(nhit));
    return v[ROOT::VecOps::Range(0, n)];
  };

  auto compute_pos = [pitch, stripCenter](const RVec<short>& strip) {
    RVec<double> out(strip.size());
    for (size_t i = 0; i < strip.size(); ++i) {
      out[i] = (static_cast<double>(strip[i]) - stripCenter) * pitch;
    }
    return out;
  };

  auto compute_x = [pitch, stripCenter](const RVec<short>& strip,
                                        const RVec<short>& plane) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t n = std::min(strip.size(), plane.size());
    RVec<double> x(n, nan);
    for (size_t i = 0; i < n; ++i) {
      if (plane[i] == 0) {
        x[i] = (static_cast<double>(strip[i]) - stripCenter) * pitch;
      }
    }
    return x;
  };

  auto compute_y = [pitch, stripCenter](const RVec<short>& strip,
                                        const RVec<short>& plane) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const size_t n = std::min(strip.size(), plane.size());
    RVec<double> y(n, nan);
    for (size_t i = 0; i < n; ++i) {
      if (plane[i] == 1) {
        y[i] = (static_cast<double>(strip[i]) - stripCenter) * pitch;
      }
    }
    return y;
  };

  auto compute_adc_sub = [pedMean, sfirst, sdelta](const RVec<short>& strip,
                                                   const RVec<short>& sample,
                                                   const RVec<double>& adc,
                                                   int nhit) {
    const size_t n = std::min({adc.size(), strip.size(), sample.size(),
                               static_cast<size_t>(nhit)});
    RVec<double> out(n);
    for (size_t i = 0; i < n; ++i) {
      const int s = static_cast<int>(sample[i]);
      const int idx = static_cast<int>(strip[i]) - sfirst;
      if (s < 0 || s >= static_cast<int>(pedMean.size()) ||
          idx < 0 || idx >= sdelta) {
        out[i] = adc[i];
      } else {
        out[i] = adc[i] - pedMean[s][idx];
      }
    }
    return out;
  };

  auto compute_q_strip = [sfirst, sdelta](const RVec<short>& strip,
                                          const RVec<short>& plane,
                                          const RVec<short>& module,
                                          const RVec<double>& adc,
                                          int nhit) {
    const size_t n = std::min({adc.size(), strip.size(), plane.size(), module.size(),
                               static_cast<size_t>(nhit)});
    std::map<std::pair<short, short>, std::vector<double>> qsum;
    for (size_t i = 0; i < n; ++i) {
      const int idx = static_cast<int>(strip[i]) - sfirst;
      if (idx < 0 || idx >= sdelta) {
        continue;
      }
      const auto key = std::make_pair(plane[i], module[i]);
      auto &vec = qsum[key];
      if (vec.empty()) {
        vec.assign(sdelta, 0.0);
      }
      vec[idx] += adc[i];
    }
    RVec<double> out(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      const int idx = static_cast<int>(strip[i]) - sfirst;
      if (idx < 0 || idx >= sdelta) {
        out[i] = 0.0;
      } else {
        const auto key = std::make_pair(plane[i], module[i]);
        auto it = qsum.find(key);
        if (it != qsum.end()) {
          out[i] = it->second[idx];
        } else {
          out[i] = 0.0;
        }
      }
    }
    return out;
  };

  auto is_good = [pedMean, pedErr, nsigma, minGood, sfirst, sdelta](const RVec<short>& strip,
                                                                   const RVec<short>& sample,
                                                                   const RVec<double>& adc,
                                                                   int nhit) {
    const int nsamples = static_cast<int>(pedMean.size());
    std::vector<std::vector<double>> ev(nsamples, std::vector<double>(sdelta, 0.0));

    const size_t n = std::min<size_t>(adc.size(), static_cast<size_t>(nhit));
    for (size_t i = 0; i < n; ++i) {
      const int s = static_cast<int>(sample[i]);
      const int idx = static_cast<int>(strip[i]) - sfirst;
      if (s < 0 || s >= nsamples) {
        continue;
      }
      if (idx < 0 || idx >= sdelta) {
        continue;
      }
      ev[s][idx] += adc[i];
    }

    int goodBins = 0;
    for (int s = 0; s < nsamples; ++s) {
      for (int idx = 0; idx < sdelta; ++idx) {
        const double perr = pedErr[s][idx];
        if (perr <= 0) {
          continue;
        }
        if ((ev[s][idx] - pedMean[s][idx]) > nsigma * perr) {
          ++goodBins;
        }
      }
    }
    return goodBins >= minGood;
  };

  ROOT::RDataFrame dfSig("fT", sigFile.Data());
  auto nTotal = dfSig.Count();
  auto dfBase = dfSig.Define("nhit", "ft.nhit")
                     .Define("strip", take_n_short, {"ft.strip", "nhit"})
                     .Define("sample", take_n_short, {"ft.sample", "nhit"})
                     .Define("plane", take_n_short, {"ft.plane", "nhit"})
                     .Define("module", take_n_short, {"ft.module", "nhit"})
                     .Define("pos_mm", compute_pos, {"strip"})
                     .Define("x_mm", compute_x, {"strip", "plane"})
                     .Define("y_mm", compute_y, {"strip", "plane"})
                     .Define("adc_sub", compute_adc_sub,
                             {"ft.strip", "ft.sample", "ft.adc", "nhit"})
                     .Define("q_strip", compute_q_strip, {"strip", "plane", "module", "adc_sub", "nhit"});

  ROOT::RDF::RNode df = dfBase;
  if (minGood > 0) {
    df = dfBase.Define("is_good", is_good,
                       {"ft.strip", "ft.sample", "ft.adc", "ft.nhit"})
               .Filter("is_good");
  }
  auto nSelected = df.Count();

  auto dfOut = df.Define("evt", "ft.evt")
                   .Define("adc", "adc_sub");

  TString outDir = gSystem->DirName(outFile);
  if (outDir.Length() == 0) {
    outDir = ".";
  }
  gSystem->mkdir(outDir.Data(), true);

  dfOut.Snapshot("uRwell", outFile.Data(),
                 {"evt", "nhit", "module", "plane", "strip", "sample", "adc", "q_strip", "pos_mm", "x_mm", "y_mm"});

  long recoEvents = 0;
  {
    TFile fOut(outFile, "UPDATE");
    TTree *tIn = (TTree*)fOut.Get("uRwell");
    if (tIn != nullptr) {
      TTreeReader tr(tIn);
      TTreeReaderValue<Int_t> evt(tr, "evt");
      TTreeReaderValue<Int_t> nhit(tr, "nhit");
      TTreeReaderValue<RVec<short>> strip(tr, "strip");
      TTreeReaderValue<RVec<short>> plane(tr, "plane");
      TTreeReaderValue<RVec<short>> module(tr, "module");
      TTreeReaderValue<RVec<double>> adc(tr, "adc");

      int oEvt = 0;
      int oNclu = 0;
      std::vector<short> oPlane, oModule, oStripMin, oStripMax, oSize;
      std::vector<double> oCharge, oPosStrip, oPosMm;

      TTree tClu("uClu", "cluster data");
      tClu.Branch("evt", &oEvt, "evt/I");
      tClu.Branch("nclu", &oNclu, "nclu/I");
      tClu.Branch("plane", &oPlane);
      tClu.Branch("module", &oModule);
      tClu.Branch("strip_min", &oStripMin);
      tClu.Branch("strip_max", &oStripMax);
      tClu.Branch("size", &oSize);
      tClu.Branch("charge", &oCharge);
      tClu.Branch("pos_strip", &oPosStrip);
      tClu.Branch("pos_mm", &oPosMm);

      int eEvt = 0;
      short ePlane = 0;
      short eModule = 0;
      int eSfirst = sfirst;
      int eSdelta = sdelta;
      std::vector<double> eAvgStrip;

      TTree tEvt("uEvt", "event reconstructed distribution");
      tEvt.Branch("evt", &eEvt, "evt/I");
      tEvt.Branch("plane", &ePlane, "plane/S");
      tEvt.Branch("module", &eModule, "module/S");
      tEvt.Branch("sfirst", &eSfirst, "sfirst/I");
      tEvt.Branch("sdelta", &eSdelta, "sdelta/I");
      tEvt.Branch("avg_strip", &eAvgStrip);

      while (tr.Next()) {
        oPlane.clear();
        oModule.clear();
        oStripMin.clear();
        oStripMax.clear();
        oSize.clear();
        oCharge.clear();
        oPosStrip.clear();
        oPosMm.clear();

        std::map<std::pair<short, short>, std::vector<double>> qsum;

        const size_t n = std::min<size_t>(adc->size(), static_cast<size_t>(*nhit));
        for (size_t i = 0; i < n; ++i) {
          const int idx = static_cast<int>((*strip)[i]) - sfirst;
          if (idx < 0 || idx >= sdelta) {
            continue;
          }
          const short pl = (i < plane->size()) ? (*plane)[i] : -1;
          const short mod = (i < module->size()) ? (*module)[i] : -1;
          const auto key = std::make_pair(pl, mod);
          auto &vec = qsum[key];
          if (vec.empty()) {
            vec.assign(sdelta, 0.0);
          }
          vec[idx] += (*adc)[i];
        }

        eEvt = *evt;
        for (const auto &kv : qsum) {
          ePlane = kv.first.first;
          eModule = kv.first.second;
          eAvgStrip = kv.second;
          for (double &v : eAvgStrip) {
            v /= static_cast<double>(nsamples);
          }
          tEvt.Fill();
        }

        oEvt = *evt;
        oNclu = 0;

        for (const auto &kv : qsum) {
          const short pl = kv.first.first;
          const short mod = kv.first.second;
          const auto &qvec = kv.second;
          bool in = false;
          int start = -1;
          int last = -1;
          int size = 0;
          double sumQ = 0.0;
          double sumPos = 0.0;
          for (int idx = 0; idx < sdelta; ++idx) {
            const double qv = qvec[idx];
            if (qv > 0.0) {
              if (!in) {
                in = true;
                start = idx;
                size = 0;
                sumQ = 0.0;
                sumPos = 0.0;
              }
              sumQ += qv;
              sumPos += qv * (sfirst + idx);
              last = idx;
              size++;
            } else if (in) {
              const double posStrip = (sumQ > 0.0) ? (sumPos / sumQ) : (sfirst + start);
              const double posMm = (posStrip - stripCenter) * pitch;
              oPlane.push_back(pl);
              oModule.push_back(mod);
              oStripMin.push_back(sfirst + start);
              oStripMax.push_back(sfirst + last);
              oSize.push_back(size);
              oCharge.push_back(sumQ);
              oPosStrip.push_back(posStrip);
              oPosMm.push_back(posMm);
              oNclu++;
              in = false;
            }
          }
          if (in) {
            const double posStrip = (sumQ > 0.0) ? (sumPos / sumQ) : (sfirst + start);
            const double posMm = (posStrip - stripCenter) * pitch;
            oPlane.push_back(pl);
            oModule.push_back(mod);
            oStripMin.push_back(sfirst + start);
            oStripMax.push_back(sfirst + last);
            oSize.push_back(size);
            oCharge.push_back(sumQ);
            oPosStrip.push_back(posStrip);
            oPosMm.push_back(posMm);
            oNclu++;
          }
        }

        tClu.Fill();
        recoEvents++;
      }

      tClu.Write("", TObject::kOverwrite);
      tEvt.Write("", TObject::kOverwrite);
    }
    fOut.Close();
  }

  const auto totalEvents = nTotal.GetValue();
  const auto selectedEvents = nSelected.GetValue();

  TString logPath = logFile;
  if (logPath.Length() == 0) {
    const TString suite = ResolveSuite();
    if (suite.Length() > 0) {
      TString base = gSystem->BaseName(sigFile);
      if (base.EndsWith(".root")) {
        base.ReplaceAll(".root", ".log");
      } else {
        base += ".log";
      }
      logPath = suite + "/logs/" + base;
    }
  }
  if (logPath.Length() > 0) {
    TString logDir = gSystem->DirName(logPath);
    if (logDir.Length() > 0) {
      gSystem->mkdir(logDir.Data(), true);
    }
    std::ofstream logf(logPath.Data());
    if (logf) {
      logf << "MPD_analysis log\n";
      logf << "input_signal=" << sigFile << "\n";
      logf << "input_pedestal=" << pedFile << "\n";
      logf << "output_root=" << outFile << "\n";
      logf << "nsigma=" << nsigma << "\n";
      logf << "sfirst=" << sfirst << "\n";
      logf << "sdelta=" << sdelta << "\n";
      logf << "minGood=" << minGood << "\n";
      logf << "pitch_mm=" << pitch << "\n";
      logf << "strip_center_input=" << stripCenterInput << "\n";
      logf << "strip_center_used=" << stripCenter << "\n";
      logf << "plot_events=" << plotEvents << "\n";
      logf << "make_plots=" << (makePlots ? "true" : "false") << "\n";
      logf << "events_total=" << totalEvents << "\n";
      logf << "events_selected=" << selectedEvents << "\n";
      logf << "events_reco=" << recoEvents << "\n";
    }
  }

  if (makePlots) {
    gStyle->SetOptStat(0);

    TString pdfOut = outFile;
    if (pdfOut.EndsWith(".root")) {
      pdfOut.ReplaceAll(".root", ".pdf");
    } else {
      pdfOut += ".pdf";
    }

    std::vector<TH1D*> hPed(nsamples, nullptr);
    for (int s = 0; s < nsamples; ++s) {
      hPed[s] = new TH1D(Form("hPedPlot_%d", s),
                         Form("Pedestal sample %d;Strip;ADC", s),
                         sdelta, sfirst, sfirst + sdelta);
      for (int idx = 0; idx < sdelta; ++idx) {
        hPed[s]->SetBinContent(idx + 1, pedMean[s][idx]);
      }
    }

    TCanvas c("c", "c", 1600, 900);
    c.Print(Form("%s[", pdfOut.Data()));
    c.Divide(2, 3);
    for (int s = 0; s < nsamples; ++s) {
      c.cd(s + 1);
      hPed[s]->Draw("hist");
    }
    c.Print(pdfOut.Data());

    if (plotEvents > 0) {
      std::vector<TH1D*> hGoodCount(nsamples, nullptr);
      std::vector<TH1D*> hGoodDist(nsamples, nullptr);
      for (int s = 0; s < nsamples; ++s) {
        hGoodCount[s] = new TH1D(Form("hGoodCount_%d", s),
                                 Form("N strips > %.0f #sigma (sample %d);N;counts", nsigma, s),
                                 21, -0.5, 20.5);
        hGoodDist[s] = new TH1D(Form("hGoodStrip_%d", s),
                                Form("Good strip distribution (sample %d);Strip;counts", s),
                                sdelta, sfirst, sfirst + sdelta);
      }

      TFile fsig(sigFile, "read");
      TTreeReader reader("fT", &fsig);
      TTreeReaderValue<Int_t> evt(reader, "ft.evt");
      TTreeReaderValue<Int_t> nhit(reader, "ft.nhit");
      TTreeReaderArray<Short_t> strip(reader, "ft.strip");
      TTreeReaderArray<Short_t> sample(reader, "ft.sample");
      TTreeReaderArray<Double_t> adc(reader, "ft.adc");

      TLatex label;
      label.SetNDC(true);

      int evCount = 0;
      while (reader.Next() && evCount < plotEvents) {
        std::vector<std::vector<double>> ev(nsamples, std::vector<double>(sdelta, 0.0));
        for (int i = 0; i < *nhit; ++i) {
          const int s = sample[i];
          const int idx = strip[i] - sfirst;
          if (s < 0 || s >= nsamples) {
            continue;
          }
          if (idx < 0 || idx >= sdelta) {
            continue;
          }
          ev[s][idx] += adc[i];
        }

        std::vector<TH1D*> hRaw(nsamples, nullptr);
        std::vector<TH1D*> hSub(nsamples, nullptr);

        TLegend *leg = nullptr;
        for (int s = 0; s < nsamples; ++s) {
          hRaw[s] = new TH1D(Form("hRaw_%d_%d", *evt, s),
                             Form("Evt %d sample %d;Strip;ADC", *evt, s),
                             sdelta, sfirst, sfirst + sdelta);
          for (int idx = 0; idx < sdelta; ++idx) {
            hRaw[s]->SetBinContent(idx + 1, ev[s][idx]);
          }
          hSub[s] = (TH1D*)hRaw[s]->Clone(Form("hSub_%d_%d", *evt, s));
          for (int idx = 0; idx < sdelta; ++idx) {
            hSub[s]->SetBinContent(idx + 1, ev[s][idx] - pedMean[s][idx]);
          }

          int goodCount = 0;
          for (int idx = 0; idx < sdelta; ++idx) {
            if (pedErr[s][idx] <= 0) {
              continue;
            }
            if ((ev[s][idx] - pedMean[s][idx]) > nsigma * pedErr[s][idx]) {
              ++goodCount;
              hGoodDist[s]->Fill(sfirst + idx);
            }
          }
          hGoodCount[s]->Fill(goodCount);

          c.cd(s + 1);
          hRaw[s]->SetLineColor(kRed + 1);
          hSub[s]->SetLineColor(kBlue + 1);
          hRaw[s]->Draw("hist");
          hSub[s]->Draw("hist same");
          if (s == 0) {
            leg = new TLegend(0.7, 0.7, 0.9, 0.9);
            leg->AddEntry(hSub[s], "Ped sub", "l");
            leg->AddEntry(hRaw[s], "Raw", "l");
            leg->Draw();
          }
          label.DrawLatex(0.2, 0.75, Form("Over %.0f #sigma = %d", nsigma, goodCount));
        }

        c.Print(pdfOut.Data());

        if (leg != nullptr) {
          delete leg;
        }

        for (int s = 0; s < nsamples; ++s) {
          delete hRaw[s];
          delete hSub[s];
        }

        ++evCount;
      }

      c.Clear();
      c.Divide(2, 3);
      for (int s = 0; s < nsamples; ++s) {
        c.cd(s + 1);
        hGoodCount[s]->Draw("hist");
      }
      c.Print(pdfOut.Data());

      c.Clear();
      c.Divide(2, 3);
      for (int s = 0; s < nsamples; ++s) {
        c.cd(s + 1);
        hGoodDist[s]->Draw("hist");
      }
      c.Print(pdfOut.Data());

      for (int s = 0; s < nsamples; ++s) {
        delete hGoodCount[s];
        delete hGoodDist[s];
      }
    }

    c.Print(Form("%s]", pdfOut.Data()));
    for (int s = 0; s < nsamples; ++s) {
      delete hPed[s];
    }
  }

  std::cout << "Output ROOT: " << outFile << "\n";
  if (logPath.Length() > 0) {
    std::cout << "Log: " << logPath << "\n";
  }
  return 0;
}

} // namespace

int mpd_ana(const TString& sigFile,
            const TString& pedFile,
            const TString& outFile,
            double nsigma = 3.0,
            int sfirst = 0,
            int sdelta = 256,
            int minGood = 6,
            double pitch = 0.0,
            double stripCenter = -999.0,
            int plotEvents = 0,
            bool makePlots = true,
            const TString& logFile = "")
{
  return RunAnalysis(sigFile, pedFile, outFile, nsigma, sfirst, sdelta, minGood, pitch, stripCenter, plotEvents, makePlots, logFile);
}

int mpd_ana(int srun,
            int prun,
            TString spath = "",
            double nsigma = 3.0,
            int sfirst = 0,
            int sdelta = 256,
            int minGood = 6,
            double pitch = 0.0,
            double stripCenter = -999.0,
            int plotEvents = 0,
            bool makePlots = true,
            const TString& logFile = "")
{
  if (spath.Length() == 0) {
    const TString suite = ResolveSuite();
    if (suite.Length() > 0) {
      spath = suite + "/../DATA/out";
    }
  }
  const TString sigFile = Form("%s/run_%04d.dat_apv.root", spath.Data(), srun);
  const TString pedFile = Form("%s/run_%04d.dat_apv.root", spath.Data(), prun);
  TString outFile;
  {
    const TString suite = ResolveSuite();
    if (suite.Length() > 0) {
      outFile = Form("%s/output/run_%04d.dat_apv.root", suite.Data(), srun);
    } else {
      outFile = Form("run_%04d.dat_apv.root", srun);
    }
  }
  return RunAnalysis(sigFile, pedFile, outFile, nsigma, sfirst, sdelta, minGood, pitch, stripCenter, plotEvents, makePlots, logFile);
}
