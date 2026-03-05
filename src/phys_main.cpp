#include "mpd_common.h"
#include "mpd_rwell_ana.h"

#include <getopt.h>
#include <iostream>
#include <string>

namespace {

void Usage(const char *prog)
{
  std::cout
      << "Usage: " << prog << " -i reco.root [options]\n\n"
      << "Options:\n"
      << "  -i  Input reco ROOT file (uClu, uEvt)\n"
      << "  -o  Output PDF file (default: $MPD_SUITE/../output/<input>_phys.pdf, derived from input)\n"
      << "  -A  Number of avg event profiles to plot from uEvt (0 disables)\n"
      << "  -P  Strip pitch in mm (override config/analysis.yaml)\n"
      << "  -c  Strip center (override config/analysis.yaml; negative = auto)\n"
      << "  -h  Show this help\n";
}

} // namespace

int main(int argc, char **argv)
{
  MpdRwellAnaOptions opt;

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
    std::cerr << "MPD_SUITE is not set. From MPD_analysis directory: "
                 "export MPD_SUITE=\"$PWD\"\n";
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

  MpdRwellAna ana(std::move(opt));
  return ana.Run();
}
