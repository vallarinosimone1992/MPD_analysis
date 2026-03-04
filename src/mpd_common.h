#pragma once

#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <string>

inline TString ResolveSuite()
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

inline void LoadLegacyDictionaries()
{
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

inline void LoadConfigParams(double &pitch, double &stripCenter)
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

inline std::string BaseName(const std::string &path)
{
  const auto pos = path.find_last_of("/");
  if (pos == std::string::npos) {
    return path;
  }
  return path.substr(pos + 1);
}

inline std::string ReplaceExtension(const std::string &name,
                                    const std::string &ext)
{
  const auto pos = name.rfind('.');
  if (pos == std::string::npos) {
    return name + ext;
  }
  return name.substr(0, pos) + ext;
}
