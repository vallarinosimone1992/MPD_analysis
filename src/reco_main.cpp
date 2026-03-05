#include "mpd_common.h"
#include "mpd_rwell_reco.h"

#include <getopt.h>
#include <iostream>
#include <string>

namespace {

void Usage(const char *prog)
{
  std::cout
      << "Usage: " << prog << " -i signal.root -p pedestal.root [options]\n\n"
      << "Options:\n"
      << "  -i  Signal input ROOT file (ft.* branches)\n"
      << "  -p  Pedestal ROOT file (ft.* branches)\n"
      << "  -o  Output reco ROOT file (default: $MPD_SUITE/../RECO_DATA/<input basename>, derived from input)\n"
      << "  -l  Log file (default: $MPD_SUITE/logs/<input basename>.log)\n"
      << "  -n  Threshold in units of pedestal RMS (nsigma)\n"
      << "  -s  First strip index used for pedestal binning\n"
      << "  -d  Number of strips used for pedestal binning\n"
      << "  -m  Minimum number of bins above threshold (0 disables filtering)\n"
      << "  -P  Strip pitch in mm (override config/analysis.yaml)\n"
      << "  -c  Strip center (override config/analysis.yaml; negative = auto)\n"
      << "  -h  Show this help\n";
}

} // namespace

int main(int argc, char **argv)
{
  MpdRwellRecoOptions opt;

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
    std::cerr << "MPD_SUITE is not set. From MPD_analysis directory: "
                 "export MPD_SUITE=\"$PWD\"\n";
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

  MpdRwellReco reco(std::move(opt));
  return reco.Run();
}
