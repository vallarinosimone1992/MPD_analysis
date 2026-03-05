#pragma once

#include <string>

struct MpdRwellAnaOptions {
  std::string input;
  std::string output;
  int avgEvents = 0;
  double pitch = 0.0;
  double stripCenter = -999.0;
};

class MpdRwellAna {
public:
  explicit MpdRwellAna(MpdRwellAnaOptions opt);
  int Run();

private:
  MpdRwellAnaOptions opt_;
};
