# Studies of Jet Quenching using Boosted W and Z Bosons Decaying Hadronically in Heavy-Ion Collisions

Master's Thesis — Francesca Zavatti  
Sapienza Università di Roma, 2024/2025

---

## What this repository contains

This repository contains the simulation and analysis code used for my Master's thesis on time-dependent jet quenching in heavy-ion collisions using boosted W and Z bosons.

The code is organised as follows:

- `WZgen.cc`, main simulation: generates pp events with Pythia8, embeds thermal background (OO or PbPb), applies a phenomenological quenching model, and computes jet substructure observables. Output: ROOT files with one TTree per run.
- `creating.C` — reads the ROOT files produced by WZgen.cc and fills 2D histograms of SoftDrop groomed mass vs N2, split by pT bin and event type (W, Z, QCD). Output: `WZ_histograms4bins.root`.
- `fit1.C`, Phase 1 fit: fits the QCD background with a 4th-degree polynomial (signal fixed to zero).
- `fit2.C`, Phase 2 fit: fits the W mass peak with a Crystal Ball function (background fixed from Phase 1).
- `fit3.C`, Phase 3 fit: fits the Z mass peak with a Crystal Ball function (background fixed from Phase 1).
- `fit4.C`, Phase 4 fit: full W+Z+background fit with all shape parameters fixed from Phases 1-3. Extracts the reconstructed mass peaks as a function of pT bin.
- `fitsPDF.C`, Toy MC study: generates pseudo-experiments from the fitted function and refits them to estimate the statistical uncertainty on the mass peaks.

---

## Dependencies

- Pythia8
- FastJet + contrib (N-subjettiness, SoftDrop, ConstituentSubtractor, EnergyCorrelator)
- ROOT 6
- `libMyLib.so` (Shared library compiled from `myLib.cc` for plotting utilities)

---

## How to run

### WZgen.cc

Compile with the local script (environment-specific — adapt include and library paths if needed):
```bash
./compile_pythia.sh WZgen
```

Run:
```bash
./WZgen <tune> <seed> <nEvts> <jetR> <massflag> <ns>
```

The main flag is `massflag`: use `6` for the signal process (HardQCD + W/Z+parton). The other values are for QCD background studies only.

Before compiling, set the medium parameters at the top of the file.
For OO: L_medium = TMath::Power(16, 1./3.), quenching_normal = 0.03, qhat = 1.0 GeV²/fm.
For PbPb: replace 16 with 208, quenching_normal = 0.10, qhat = 4.0 GeV²/fm.
Three quenching scenarios are available inside the particle loop — uncomment the desired one (search for `SCENARIO` in the file).


### Plotting Library Compilation
Before running the analysis and plotting macros, compile the custom shared utility library:
```bash
g++ -shared -fPIC $(root-config --cflags --libs) myLib.cc -o libMyLib.so
```s

### creating.C

Edit the default arguments at the top of the file (input file list and output path), then run:
```bash
root -l creating.C
```

### fit1.C — fit4.C

Set `EOS_BASE` at the top of each file to your output directory, then run in order for each pT bin (`pt150_350`, `pt350_500`, `pt500_700`, `pt700_1000`):

```bash
EOS=/your/path
# Phase 1
root -l -q "fit1.C(\"$EOS/WZ_histograms4bins.root\",\"WZfit_phase1\",\"pt150_350\")"
# Phase 2 (needs Phase 1 output)
root -l -q "fit2.C(\"$EOS/WZ_histograms4bins.root\",\"$EOS/figuresB/WZfit_phase1_pt150_350/phase1_results.root\",\"WZfit_phase2\",\"pt150_350\")"
# Phase 3 (needs Phase 1 output)
root -l -q "fit3.C(\"$EOS/WZ_histograms4bins.root\",\"$EOS/figuresB/WZfit_phase1_pt150_350/phase1_results.root\",\"WZfit_phase3\",\"pt150_350\")"
# Phase 4 (needs Phases 1-3 output, "auto" reconstructs paths from EOS_BASE)
root -l -q "fit4.C(\"$EOS/WZ_histograms4bins.root\",\"auto\",\"auto\",\"auto\",\"WZfit_phase4\",\"pt150_350\")"
```
Repeat for each pT bin.

### fitsPDF.C

Replace the parameter arrays at the top of the file with the results from fit1.C–fit4.C for your sample. Set `nEventsList` to the number of events in each pT bin after the N2 < 0.18 cut. Then run:
```bash
root -l fitsPDF.C
```

---

## Reference

L. Apolinário, J. G. Milhano, G. P. Salam, C. A. Salgado, *Probing the time structure of the quark-gluon plasma with top quarks*, Phys. Rev. Lett. 120 (2018) 232301. [arXiv:1711.01161](https://arxiv.org/abs/1711.01161)
