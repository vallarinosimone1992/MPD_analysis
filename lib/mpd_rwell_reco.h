#pragma once

#include <string>

struct MpdRwellRecoOptions {
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

class MpdRwellReco {
public:
  explicit MpdRwellReco(MpdRwellRecoOptions opt);
  int Run();

private:
  MpdRwellRecoOptions opt_;
};
