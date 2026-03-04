# MPD_analysis

MPD analysis suite based on ROOT/RDataFrame. Two separate steps:
- **Reco**: pedestal subtraction and event/cluster building.
- **Ana**: physics plots from the reconstructed TTree(s).

## Signal vs pedestal
- Signal: decoded ROOT files with signal events `run_XXXX.dat_apv.root`.
- Pedestal: decoded ROOT files with pedestal `run_YYYY.dat_apv.root`.
- Reco produces a ROOT file with `TTree` `uRwell`, `uClu`, `uEvt`.
- Ana produces a PDF with plots.

## Structure
- `bin/` compiled executables (`reco`, `ana`) + environment setup script
- `src/` C++ sources
- `config/` configuration (pitch/strip_center)
- `macros/` legacy macros (not used)
- `../RECO_DATA/` reconstructed ROOT files (sibling of MPD_analysis)
- `../output/` PDFs from analysis (sibling of MPD_analysis)
- `src/legacy/` minimal headers for ROOT dictionaries (MPDdb/APVdb/DAQpar)
- `logs/` reconstruction logs

## Setup
Set the environment variable (must point to `devel/MPD_analysis`):
```bash
export MPD_SUITE=/path/to/MPD_dev/devel/MPD_analysis
```
Or use the setup script (recommended):
```bash
source devel/MPD_analysis/bin/setup_mpd_env.sh
```

Decoded data are expected in:
- `$MPD_SUITE/../DATA/out/`

### Included test dataset
For quick tests the repo includes a small dataset in:
- `$MPD_SUITE/data/` (`run_0254` signal, `run_0253` pedestal)

If present, these files are used as defaults by `reco`.

### Data symlink
The path `$MPD_SUITE/../DATA` must point to the MPD data directory. In this repo
there is already a symlink `devel/DATA` to `old_software/MPD_Eth/gdaq_mpd_eth/data`.
To recreate it:
```bash
ln -s /path/to/MPD_dev/old_software/MPD_Eth/gdaq_mpd_eth/data \
  /path/to/MPD_dev/devel/DATA
```

## Build (CMake)
Build the executables:
```bash
cd devel/MPD_analysis
cmake -S . -B build
cmake --build build -j
```

## Quick run
Reco with defaults (uses bundled test data if present):
```bash
./devel/MPD_analysis/bin/reco
```

Analysis (plots from `uClu`):
```bash
./devel/MPD_analysis/bin/ana -i $MPD_SUITE/../RECO_DATA/run_0254.dat_apv.root
```

Default outputs:
- RECO ROOT: `$MPD_SUITE/../RECO_DATA/<input basename>`
- PDF: `$MPD_SUITE/../output/<input basename>_phys.pdf`
- LOG: `$MPD_SUITE/logs/<input basename>.log`

## Recommended usage (signal + pedestal)
1. Pick a signal run and a pedestal run (typically the previous run).
2. Run reconstruction:
```bash
./devel/MPD_analysis/bin/reco \
  -i $MPD_SUITE/../DATA/out/run_0254.dat_apv.root \
  -p $MPD_SUITE/../DATA/out/run_0253.dat_apv.root \
  -o run_0254.dat_apv.root
```
3. Run analysis:
```bash
./devel/MPD_analysis/bin/ana \
  -i $MPD_SUITE/../RECO_DATA/run_0254.dat_apv.root \
  -A 5
```
4. Open the PDF in `devel/output/`.

## Executables (bin/)
- `reco` - reconstruction (pedestal subtraction + reconstructed TTree)
- `ana` - physics analysis and plots
- `setup_mpd_env.sh` - MPD_SUITE setup + standard directories

Main options (reco):
- `-i` signal input ROOT
- `-p` pedestal ROOT
- `-o` output ROOT (if relative, it is written to `../RECO_DATA/`)
- `-l` log file (if relative, it is written to `logs/`)
- `-n` threshold in RMS units (nsigma)
- `-s` first strip
- `-d` number of strips
- `-m` minimum number of bins above threshold (0 disables filtering)
- `-P` strip pitch in mm (override config)
- `-c` strip center (override config; negative = auto)

Main options (ana):
- `-i` input reco ROOT (uClu)
- `-o` output PDF (if relative, it is written to `../output/`)
- `-A` number of average event profiles (uEvt) appended at the end of the PDF
- `-P` strip pitch in mm (override config)
- `-c` strip center (override config; negative = auto)

## Sources (src/)
- `reco_main.cpp` - reconstruction.
- `phys_main.cpp` - analysis from `uClu` (+ optional mean profiles from `uEvt`).

## Operational notes
- Default strip range: `0..255` (`-s 0 -d 256`).
- Pitch/strip center are read from `config/analysis.yaml` (override with `-P`/`-c`).
- If `strip_center < 0` the strip center is estimated from min/max strip in the run.
- Plane convention: `plane=0 -> X`, `plane=1 -> Y` (from `ft.plane` in the decoded data).

## Definitions (good strip / good bins)
For each event and each sample `s`, build an array of strips `idx` (bins):
- `ev[s][idx]` = sum of ADC for sample `s` and strip `idx` within `sfirst..sfirst+sdelta-1`.
- `pedMean[s][idx]` = pedestal mean (TProfile on pedestal run).
- `pedErr[s][idx]` = pedestal RMS (TProfile error option `"s"`).

**Good strip** (all samples): strip `idx` is good if it is above threshold in **all 6 samples**:
```
for all s in [0..5]:
  pedErr[s][idx] > 0  AND  (ev[s][idx] - pedMean[s][idx]) > nsigma * pedErr[s][idx]
```

**Number of good strips**: count of strips that satisfy the condition above (all samples).

**Event selection (good bins)**: count good strips (global `goodBins`).
Event is accepted if `goodBins >= minGood` (default `minGood=6`).

## Definitions (Q_strip / clustering)
For each event and each plane/module:
- `Q_strip` = simple sum of ped-sub ADC over the 6 samples for the same strip.
- **Good strip**: strip is good only if it is above threshold in all 6 samples.
- **Cluster**: contiguous good strips (1 or more).

For each cluster:
- `charge` = sum of `Q_strip` in the cluster.
- `significance` = cluster charge (same as `charge`).
- `pos_strip` = weighted mean `sum(Q_strip * strip) / sum(Q_strip)`.
- `pos_mm` = `(pos_strip - strip_center) * pitch_mm`.

Tree `uClu` contains: `evt, nclu, plane, module, strip_min, strip_max, size, charge, pos_strip, pos_mm`.

## Definitions (event reconstruction)
For each event and each plane/module:
- `avg_strip[idx]` = mean of the 6 samples after pedestal subtraction.
- `good_strip[idx]` = `avg_strip[idx]` only for good strips, zero otherwise.

Tree `uEvt` contains:
- `evt, plane, module, sfirst, sdelta, avg_strip, good_strip`.

## Legacy classes
To avoid errors when ROOT reads UserInfo (MPDdb/APVdb/DAQpar), we keep minimal
headers in `src/legacy/` and load the dictionary at runtime. The analysis remains
independent from `old_software/`.
