#include "mpd_common.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Rtypes.h>
#include <TFile.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>

#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using ROOT::VecOps::RVec;

namespace {

struct Options {
  std::string input;
  std::string pedestal;
  std::string output;
  std::string log;
  double nsigma = 3.0;
  int sfirst = 0;
  int sdelta = 256;
  int minGood = 6;
  double pitch = 0.0;
  double stripCenter = -999.0;
};

void Usage(const char *prog)
{
  std::cout
      << "Usage: " << prog << " -i signal.root -p pedestal.root [options]\n\n"
      << "Options:\n"
      << "  -i  Signal input ROOT file (ft.* branches)\n"
      << "  -p  Pedestal ROOT file (ft.* branches)\n"
      << "  -o  Output reco ROOT file (default: RECO_DATA/<input basename>)\n"
      << "  -l  Log file (default: MPD_SUITE/logs/<input basename>.log)\n"
      << "  -n  Threshold in units of pedestal RMS (nsigma)\n"
      << "  -s  First strip index used for pedestal binning\n"
      << "  -d  Number of strips used for pedestal binning\n"
      << "  -m  Minimum number of bins above threshold (0 disables filtering)\n"
      << "  -P  Strip pitch in mm (override config/analysis.yaml)\n"
      << "  -c  Strip center (override config/analysis.yaml; negative = auto)\n"
      << "  -h  Show this help\n";
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

TH1D *ProfileToTH1D(const TProfile *prof, const TString &name)
{
  if (!prof) {
    return nullptr;
  }
  const int nBins = prof->GetNbinsX();
  const double xMin = prof->GetXaxis()->GetXmin();
  const double xMax = prof->GetXaxis()->GetXmax();
  TH1D *hist = new TH1D(name, prof->GetTitle(), nBins, xMin, xMax);
  for (int i = 1; i <= nBins; ++i) {
    hist->SetBinContent(i, prof->GetBinContent(i));
    hist->SetBinError(i, prof->GetBinError(i));
  }
  return hist;
}

int RunReco(const Options &opt)
{
  const TString sigFile = opt.input.c_str();
  const TString pedFile = opt.pedestal.c_str();
  const TString outFile = opt.output.c_str();
  const TString logFile = opt.log.c_str();

  LoadLegacyDictionaries();
  ROOT::EnableImplicitMT(0);

  double pitch = opt.pitch;
  double stripCenter = opt.stripCenter;
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
  auto pedProfiles = BuildPedestals(dfPed, nsamples, opt.sfirst, opt.sdelta);
  for (auto &p : pedProfiles) {
    p.GetValue();
  }

  std::vector<TProfile *> pedPtr(nsamples, nullptr);
  for (int i = 0; i < nsamples; ++i) {
    pedPtr[i] = pedProfiles[i].GetPtr();
    if (pedPtr[i] != nullptr) {
      pedPtr[i]->SetErrorOption("s"); // Use RMS (not error on the mean)
    }
  }

  std::vector<std::vector<double>> pedMean(nsamples, std::vector<double>(opt.sdelta, 0.0));
  std::vector<std::vector<double>> pedErr(nsamples, std::vector<double>(opt.sdelta, 0.0));
  for (int is = 0; is < nsamples; ++is) {
    if (pedPtr[is] == nullptr) {
      continue;
    }
    for (int idx = 0; idx < opt.sdelta; ++idx) {
      const int strip = opt.sfirst + idx;
      const int bin = pedPtr[is]->FindBin(strip);
      pedMean[is][idx] = pedPtr[is]->GetBinContent(bin);
      pedErr[is][idx] = pedPtr[is]->GetBinError(bin);
    }
  }

  auto take_n_short = [](const RVec<short> &v, int nhit) {
    const size_t n = std::min<size_t>(v.size(), static_cast<size_t>(nhit));
    return v[ROOT::VecOps::Range(0, n)];
  };

  auto compute_pos = [pitch, stripCenter](const RVec<short> &strip) {
    RVec<double> out(strip.size());
    for (size_t i = 0; i < strip.size(); ++i) {
      out[i] = (static_cast<double>(strip[i]) - stripCenter) * pitch;
    }
    return out;
  };

  auto compute_x = [pitch, stripCenter](const RVec<short> &strip,
                                        const RVec<short> &plane) {
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

  auto compute_y = [pitch, stripCenter](const RVec<short> &strip,
                                        const RVec<short> &plane) {
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

  auto compute_adc_sub = [pedMean, sfirst = opt.sfirst, sdelta = opt.sdelta](
                              const RVec<short> &strip,
                              const RVec<short> &sample,
                              const RVec<double> &adc,
                              int nhit) {
    const size_t n = std::min({adc.size(), strip.size(), sample.size(),
                               static_cast<size_t>(nhit)});
    RVec<double> out(n);
    for (size_t i = 0; i < n; ++i) {
      const int s = static_cast<int>(sample[i]);
      const int idx = static_cast<int>(strip[i]) - sfirst;
      if (s < 0 || s >= static_cast<int>(pedMean.size()) || idx < 0 || idx >= sdelta) {
        out[i] = adc[i];
      } else {
        out[i] = adc[i] - pedMean[s][idx];
      }
    }
    return out;
  };

  auto compute_q_strip = [sfirst = opt.sfirst, sdelta = opt.sdelta](
                            const RVec<short> &strip,
                            const RVec<short> &plane,
                            const RVec<short> &module,
                            const RVec<double> &adc,
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
        }
      }
    }
    return out;
  };

  auto is_good = [pedMean, pedErr, opt](const RVec<short> &strip,
                                        const RVec<short> &sample,
                                        const RVec<double> &adc,
                                        int nhit) {
    const int nsamples = static_cast<int>(pedMean.size());
    std::vector<std::vector<double>> ev(nsamples, std::vector<double>(opt.sdelta, 0.0));

    const size_t n = std::min<size_t>(adc.size(), static_cast<size_t>(nhit));
    for (size_t i = 0; i < n; ++i) {
      const int s = static_cast<int>(sample[i]);
      const int idx = static_cast<int>(strip[i]) - opt.sfirst;
      if (s < 0 || s >= nsamples) {
        continue;
      }
      if (idx < 0 || idx >= opt.sdelta) {
        continue;
      }
      ev[s][idx] += adc[i];
    }

    int goodBins = 0;
    for (int idx = 0; idx < opt.sdelta; ++idx) {
      bool allGood = true;
      for (int s = 0; s < nsamples; ++s) {
        const double perr = pedErr[s][idx];
        if (perr <= 0) {
          allGood = false;
          break;
        }
        if ((ev[s][idx] - pedMean[s][idx]) <= opt.nsigma * perr) {
          allGood = false;
          break;
        }
      }
      if (allGood) {
        ++goodBins;
      }
    }
    return goodBins >= opt.minGood;
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
                     .Define("q_strip", compute_q_strip,
                             {"strip", "plane", "module", "adc_sub", "nhit"});

  ROOT::RDF::RNode df = dfBase;
  if (opt.minGood > 0) {
    df = dfBase.Define("is_good", is_good, {"ft.strip", "ft.sample", "ft.adc", "ft.nhit"})
             .Filter("is_good");
  }
  auto nSelected = df.Count();

  auto dfOut = df.Define("evt", "ft.evt").Define("adc", "adc_sub");

  TString outDir = gSystem->DirName(outFile);
  if (outDir.Length() == 0) {
    outDir = ".";
  }
  gSystem->mkdir(outDir.Data(), true);

  dfOut.Snapshot("uRwell", outFile.Data(),
                 {"evt", "nhit", "module", "plane", "strip", "sample", "adc", "q_strip",
                  "pos_mm", "x_mm", "y_mm"});

  long recoEvents = 0;
  {
    TFile fOut(outFile, "UPDATE");

    for (int s = 0; s < nsamples; ++s) {
      TH1D *hPed = ProfileToTH1D(pedPtr[s], Form("hPed_%d", s));
      if (hPed != nullptr) {
        hPed->Write("", TObject::kOverwrite);
        delete hPed;
      }
    }
    TParameter<double> pN("nsigma", opt.nsigma);
    pN.Write("nsigma", TObject::kOverwrite);

    TTree *tIn = (TTree *)fOut.Get("uRwell");
    if (tIn != nullptr) {
      TTreeReader tr(tIn);
      TTreeReaderValue<Int_t> evt(tr, "evt");
      TTreeReaderValue<Int_t> nhit(tr, "nhit");
      TTreeReaderValue<RVec<short>> strip(tr, "strip");
      TTreeReaderValue<RVec<short>> sample(tr, "sample");
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
      int eSfirst = opt.sfirst;
      int eSdelta = opt.sdelta;
      std::vector<double> eAvgStrip;
      std::vector<double> eGoodStrip;

      TTree tEvt("uEvt", "event reconstructed distribution");
      tEvt.Branch("evt", &eEvt, "evt/I");
      tEvt.Branch("plane", &ePlane, "plane/S");
      tEvt.Branch("module", &eModule, "module/S");
      tEvt.Branch("sfirst", &eSfirst, "sfirst/I");
      tEvt.Branch("sdelta", &eSdelta, "sdelta/I");
      tEvt.Branch("avg_strip", &eAvgStrip);
      tEvt.Branch("good_strip", &eGoodStrip);

      struct PlaneData {
        std::vector<double> qsum;
        std::vector<std::vector<double>> ev;
      };

      while (tr.Next()) {
        oPlane.clear();
        oModule.clear();
        oStripMin.clear();
        oStripMax.clear();
        oSize.clear();
        oCharge.clear();
        oPosStrip.clear();
        oPosMm.clear();

        std::map<std::pair<short, short>, PlaneData> pdata;

        const size_t n = std::min({adc->size(), strip->size(), sample->size(), plane->size(), module->size(),
                                   static_cast<size_t>(*nhit)});
        for (size_t i = 0; i < n; ++i) {
          const int idx = static_cast<int>((*strip)[i]) - opt.sfirst;
          if (idx < 0 || idx >= opt.sdelta) {
            continue;
          }
          const int s = (*sample)[i];
          const short pl = (*plane)[i];
          const short mod = (*module)[i];
          const auto key = std::make_pair(pl, mod);
          auto &pd = pdata[key];
          if (pd.qsum.empty()) {
            pd.qsum.assign(opt.sdelta, 0.0);
            pd.ev.assign(nsamples, std::vector<double>(opt.sdelta, 0.0));
          }
          pd.qsum[idx] += (*adc)[i];
          if (s >= 0 && s < nsamples) {
            pd.ev[s][idx] += (*adc)[i];
          }
        }

        eEvt = *evt;
        for (const auto &kv : pdata) {
          ePlane = kv.first.first;
          eModule = kv.first.second;
          const auto &pd = kv.second;
          eAvgStrip = pd.qsum;
          for (double &v : eAvgStrip) {
            v /= static_cast<double>(nsamples);
          }
          eGoodStrip.assign(opt.sdelta, 0.0);
          for (int idx = 0; idx < opt.sdelta; ++idx) {
            bool good = true;
            for (int s = 0; s < nsamples; ++s) {
              if (pedErr[s][idx] <= 0.0) {
                good = false;
                break;
              }
              if (pd.ev[s][idx] <= opt.nsigma * pedErr[s][idx]) {
                good = false;
                break;
              }
            }
            if (good) {
              eGoodStrip[idx] = eAvgStrip[idx];
            }
          }
          tEvt.Fill();
        }

        oEvt = *evt;
        oNclu = 0;

        for (const auto &kv : pdata) {
          const short pl = kv.first.first;
          const short mod = kv.first.second;
          const auto &pd = kv.second;
          const auto &qvec = pd.qsum;
          std::vector<char> goodMask(opt.sdelta, 0);
          for (int idx = 0; idx < opt.sdelta; ++idx) {
            bool good = true;
            for (int s = 0; s < nsamples; ++s) {
              if (pedErr[s][idx] <= 0.0) {
                good = false;
                break;
              }
              if (pd.ev[s][idx] <= opt.nsigma * pedErr[s][idx]) {
                good = false;
                break;
              }
            }
            if (good) {
              goodMask[idx] = 1;
            }
          }
          bool in = false;
          int start = -1;
          int last = -1;
          int size = 0;
          double sumQ = 0.0;
          double sumPos = 0.0;
          for (int idx = 0; idx < opt.sdelta; ++idx) {
            if (goodMask[idx]) {
              const double qv = qvec[idx];
              if (!in) {
                in = true;
                start = idx;
                size = 0;
                sumQ = 0.0;
                sumPos = 0.0;
              }
              sumQ += qv;
              sumPos += qv * (opt.sfirst + idx);
              last = idx;
              size++;
            } else if (in) {
              const double posStrip = (sumQ > 0.0) ? (sumPos / sumQ) : (opt.sfirst + start);
              const double posMm = (posStrip - stripCenter) * pitch;
              oPlane.push_back(pl);
              oModule.push_back(mod);
              oStripMin.push_back(opt.sfirst + start);
              oStripMax.push_back(opt.sfirst + last);
              oSize.push_back(size);
              oCharge.push_back(sumQ);
              oPosStrip.push_back(posStrip);
              oPosMm.push_back(posMm);
              oNclu++;
              in = false;
            }
          }
          if (in) {
            const double posStrip = (sumQ > 0.0) ? (sumPos / sumQ) : (opt.sfirst + start);
            const double posMm = (posStrip - stripCenter) * pitch;
            oPlane.push_back(pl);
            oModule.push_back(mod);
            oStripMin.push_back(opt.sfirst + start);
            oStripMax.push_back(opt.sfirst + last);
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

  if (!logFile.IsNull()) {
    TString logDir = gSystem->DirName(logFile);
    if (logDir.Length() > 0) {
      gSystem->mkdir(logDir.Data(), true);
    }
    std::ofstream logf(logFile.Data());
    if (logf) {
      logf << "MPD reco log\n";
      logf << "input_signal=" << sigFile << "\n";
      logf << "input_pedestal=" << pedFile << "\n";
      logf << "output_root=" << outFile << "\n";
      logf << "nsigma=" << opt.nsigma << "\n";
      logf << "sfirst=" << opt.sfirst << "\n";
      logf << "sdelta=" << opt.sdelta << "\n";
      logf << "minGood=" << opt.minGood << "\n";
      logf << "pitch_mm=" << pitch << "\n";
      logf << "strip_center_input=" << stripCenterInput << "\n";
      logf << "strip_center_used=" << stripCenter << "\n";
      logf << "events_total=" << totalEvents << "\n";
      logf << "events_selected=" << selectedEvents << "\n";
      logf << "events_reco=" << recoEvents << "\n";
    }
  }

  std::cout << "Output RECO: " << outFile << "\n";
  if (!logFile.IsNull()) {
    std::cout << "Log: " << logFile << "\n";
  }
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
  while ((c = getopt_long(argc, argv, "hi:p:o:l:n:s:d:m:P:c:", longopts, nullptr)) != -1) {
    switch (c) {
      case 'i':
        opt.input = optarg;
        break;
      case 'p':
        opt.pedestal = optarg;
        break;
      case 'o':
        opt.output = optarg;
        break;
      case 'l':
        opt.log = optarg;
        break;
      case 'n':
        opt.nsigma = std::stod(optarg);
        break;
      case 's':
        opt.sfirst = std::stoi(optarg);
        break;
      case 'd':
        opt.sdelta = std::stoi(optarg);
        break;
      case 'm':
        opt.minGood = std::stoi(optarg);
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
    std::cerr << "MPD_SUITE is not set. Example: export MPD_SUITE=/path/to/MPD_dev/devel/MPD_analysis\n";
    return 1;
  }

  const std::string localDataDir = std::string(suite.Data()) + "/data";
  const std::string dataDir = std::string(suite.Data()) + "/../DATA/out";
  const std::string recoDir = std::string(suite.Data()) + "/../RECO_DATA";
  const std::string logDir = std::string(suite.Data()) + "/logs";

  auto exists = [](const std::string &path) {
    return gSystem->AccessPathName(path.c_str()) == false;
  };

  const std::string localSig = localDataDir + "/run_0254.dat_apv.root";
  const std::string localPed = localDataDir + "/run_0253.dat_apv.root";
  const bool haveLocal = exists(localSig) && exists(localPed);

  if (opt.input.empty()) {
    opt.input = haveLocal ? localSig : (dataDir + "/run_0254.dat_apv.root");
  }
  if (opt.pedestal.empty()) {
    opt.pedestal = haveLocal ? localPed : (dataDir + "/run_0253.dat_apv.root");
  }

  if (opt.output.empty()) {
    opt.output = recoDir + "/" + BaseName(opt.input);
  } else if (opt.output.find('/') == std::string::npos) {
    opt.output = recoDir + "/" + opt.output;
  }

  if (opt.log.empty()) {
    const std::string base = BaseName(opt.input);
    opt.log = logDir + "/" + ReplaceExtension(base, ".log");
  } else if (opt.log.find('/') == std::string::npos) {
    opt.log = logDir + "/" + opt.log;
  }

  gSystem->mkdir(recoDir.c_str(), true);
  gSystem->mkdir(logDir.c_str(), true);

  return RunReco(opt);
}
