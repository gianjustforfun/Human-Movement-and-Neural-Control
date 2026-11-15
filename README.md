# Experiment E — Personalized gait-phase calibration

This project uses **one `main.m`** and six ordinary MATLAB functions. It reads only the original Camargo (2021) and Lhoste (2026) datasets. All outputs are written inside `Experiment_E/results/`. The datasets are never modified.

## Files

```text
Project/
├── DATASET - CAMARGO (2021)/
├── DATASET - LHOSTE (2026)/
└── Experiment_E/
    ├── main.m                  Run this file
    ├── readCamargo.m           Camargo markers, gait cycles, 204 causal features
    ├── readLhoste.m            Lhoste controller cycles and continuous histories
    ├── trainPhaseModel.m       Camargo-only circular phase Ridge
    ├── evaluatePhase.m         Original vs temporally rescaled phase
    ├── predictTransitions.m    E08/E11/E12 transition models and ablation
    ├── trainLhosteModel.m      Optional held-out Lhoste-only phase comparison
    ├── README.md
    └── results/                Created automatically
        └── cache/              Subject data and phase predictions
```

## Running

Open `main.m` and click **Run**. The dataset folders must be beside `Experiment_E`, as shown above. The project needs the original MATLAB files and dataset CSV/MAT files; no outputs from other experiments are used.

Flags near the top of `main.m`:

- `redo = false`: rebuild the subject caches from the original datasets when set to `true`.
- `redoModel = false`: retrain the Camargo phase model when set to `true`.
- `redoPhase = false`: regenerate Lhoste phase predictions when set to `true`.
- `redoTransitions = false`: recompute the transition models and event scores when set to `true`.
- `runLhosteComparison = true`: run the optional Lhoste-only phase comparison.
- `redoLhoste = false`: retrain/re-evaluate this additional comparison when set to `true`.

With all flags `false`, existing verified caches and results are reused. After updating `predictTransitions.m`, use `redoTransitions = true` once if transition-result CSVs already exist. Rebuilding an earlier stage also triggers all dependent downstream analyses.

## Scientific methods

**Camargo preprocessing (E01, E02, E04).** Marker resets (`diff < -50`) define heel strikes and toe-offs. Valid hip/knee gait cycles are assembled into continuous within-trial histories. A 1,000-ms history at 50 Hz gives 51 points per channel and 204 causal features: hip angle, knee angle, hip velocity, and knee velocity. Phase derived from completed cycles is used only as a training/validation label.

**Camargo model (E05).** Deterministic, subject-disjoint 17/5 split (`rng(574)`). Ridge models predict cosine and sine of circular phase. The original selection compares 500- and 1,000-ms histories and four regularization strengths. Only sampled Camargo training windows are fitted; all held-out windows are evaluated with equal subject weight. The final model is refitted on all 22 Camargo subjects without using Lhoste.

**Lhoste preprocessing and calibration (E03, E06).** Controller `3→0` identifies cycle boundaries; `1→2` is the reference stance-to-swing controller event, **not** measured anatomical toe-off. Kinematics are resampled to 200 Hz. The first **K=5 or K=10 complete cycles for each leg** determine personalized cycle duration, phase threshold and temporal scale. Later cycles are the evaluation set. Segment indices preserve chronology across recording resets.

**Phase transfer (E07).** Applies the fixed Camargo model with original or temporally rescaled past/current kinematic histories. Angular velocities receive the matching time-scale transformation. K5 and K10 are assessed separately, and also on identical `post10_common` anchors. Group results average subject-level metrics, not pooled windows.

**Controller transition prediction (E08, E11, E12).** Five methods are compared on identical eligible anchors: `Original`, `TimeRescaled`, `E11_A_PhaseRate`, `E12_PhaseOnlyRidge` (3 phase/timing features), and `E11_B_KinematicRidge` (the same 3 features plus 4 hip/knee kinematics). Both Ridge models use precisely the same K-cycle calibration examples, residual time-to-event target, lambda=2, standardization procedure and correction clamp (±0.35 s). Method A uses a causal 160-ms phase-rate lookback. Corrections apply only when Method A predicts at most 1 s remaining.

The alert horizons are 100, 200 and 300 ms. Alerts require two consecutive positive 20-ms anchors and a refractory interval of `max(0.3 s, 0.5 × calibrated cycle duration)`. Events must occur 20 ms or more after the alarm, within the horizon, and can match only one alarm. Device exposure counts overlapping left/right time bins once. Controller events after calibration are scoring references **only**, never inference inputs. The results remain exploratory and are not a validated closed-loop controller.

**Optional Lhoste-only comparison.** `trainLhosteModel.m` fits the same circular phase Ridge form used in E05 using only Lhoste: one mixed-condition model from 14 training subjects, tested on five entirely unseen subjects. The 14/5 split is reproducible (`rng(574)`). The history is fixed at 1,000 ms (204 causal features), and lambda is the value selected on Camargo, so no parameters are tuned using the five Lhoste test subjects. Training samples at most 180 equally spaced windows per continuous Lhoste run. The controller-defined full-cycle phase is an **offline supervised label** during training, not a predictor input. This phase reference is not anatomically measured gait phase.

The three methods are compared on exactly the **same post-K5 E07 anchors** of the five test subjects: `Camargo Original`, `Camargo + K5`, and `Lhoste-only`. Each subject contributes equally to the reported group MAE. Lhoste-only is trained on many cycles from **other** exoskeleton users, whereas Camargo + K5 uses five calibration cycles per leg from the current test user. The data requirements are therefore different. Do not directly compare the new five-subject averages against the existing E07 19-subject averages.

The offline E09/E10 oracles are intentionally excluded from the final executable.

## Outputs

`results/` includes:

- `phase_model.mat`, `phase_validation_models.csv`, `phase_validation_subjects.csv`
- `calibration_summary.csv`, `phase_subject_metrics.csv`, `phase_group_summary.csv`
- `transition_subject_metrics.csv`, `transition_group_summary.csv`
- `transition_training_QC.csv`, `transition_ablation_H200.csv`
- Three original presentation figures: `phase_calibration.png`, `transition_recall.png`, `transition_ablation.png`
- Optional comparison: `lhoste_only_model.mat`, `lhoste_only_subject_metrics.csv`, `lhoste_only_group_summary.csv`, `lhoste_only_comparison.png`
- `cache/`: per-subject source-derived arrays, plus 38 phase prediction MAT files

Main phase-transfer and transition group metrics give each of the 19 Lhoste subjects equal weight. The optional Lhoste-only comparison averages the five held-out subject scores with equal weight. H200/K5 transition comparisons use common post-10 evaluation anchors. The full transition CSV also includes K10 and horizons 100/300 ms.

## Verification status

**Verified locally from the MATLAB logs provided on October 25, 2026:**

| Check | Refactored result |
|---|---:|
| Camargo subjects/trials/cycles | 22 / 680 / 14,942 |
| Camargo causal windows | 835,934 |
| Lhoste subjects, SM/Transparent cycles | 19 / 1,467 / 1,497 |
| K5/K10 subject-condition calibrations | 76 |
| Camargo phase validation MAE | 2.38 phase pp |
| Camargo estimates within ±10 pp | 98.4% |
| Selected history, lambda | 1,000 ms / 0.0001 |
| Final Camargo training samples | 122,400 |
| SM K5 phase MAE (original → rescaled) | 12.39 → 7.36 pp |
| Transparent K5 phase MAE (original → rescaled) | 14.60 → 8.37 pp |
| SM post10_common K5/K10 rescaled phase MAE | 7.27 / 7.18 pp |
| Transparent post10_common K5/K10 rescaled phase MAE | 8.25 / 8.29 pp |

**Transition prediction verified locally, October 25, 2026.** The final MATLAB execution completed without reported errors. The uploaded `transition_subject_metrics.csv`, `transition_group_summary.csv`, `transition_training_QC.csv`, and `transition_ablation_H200.csv` reproduce the main E08/E11/E12 reference metrics to the reported precision. The group-level phase and transition scores were independently recomputed from subject rows. There are 76 personalized calibration fits (19 subjects × 2 conditions × 2 K values), 2,280 transition subject-score rows and 120 group-score rows. K5/K10 common-anchor comparisons and the five methods' H200 matched event counts are internally consistent.

### Final transition performance (K5, H200, post10_common)

| Method | SM recall | Transparent recall |
| --- | ---: | ---: |
| Original | 9.53% | 5.49% |
| TimeRescaled | 19.07% | 11.00% |
| E11-A causal phase rate | 20.91% | 10.80% |
| 3-feature personalized Ridge | 26.64% | 13.83% |
| 7-feature personalized Ridge | **31.37%** | **19.03%** |

The final seven-feature Ridge achieves **32.61% precision and 25.93 false alarms/device-minute for SM** and **21.42% precision and 24.61 false alarms/device-minute for Transparent**, on the same H200/K5 post-10 anchors. These absolute results do not support reliable autonomous exoskeleton control.

### Controlled kinematic feature ablation (H200, matched post10_common)

| Condition | K cycles per leg | Recall gain (7 minus 3 features) | Extra false alarms/device-minute |
| --- | ---: | ---: | ---: |
| SM | 5 | +4.73 pp | +3.09 |
| SM | 10 | +5.69 pp | +2.71 |
| Transparent | 5 | +5.19 pp | +2.23 |
| Transparent | 10 | +1.70 pp | +1.14 |

The improvement in H200 mean recall is accompanied by increased false alarms in all four comparisons. It is not uniform across subjects: with SM/K5, the full Ridge has better recall in only 8 of 19 subjects (5 worse, 6 tied). At H300, the added kinematics **reduce** Transparent mean recall by 1.36 pp (K5) and 0.69 pp (K10).

**Scope of verification:** Matching counts, selected settings, reported metrics, and internal CSV consistency are verified. This is not a byte-for-byte comparison of all feature matrices, coefficients or time-series predictions against the archived E00–E12 MAT files. The Lhoste 1→2 reference is controller-defined, not independently measured anatomical toe-off. The same Lhoste cohort was explored during model development, so this is retrospective, exploratory evidence, not an independent test of deployment performance.

## Lhoste-only comparison: verification and results (October 25, 2026)

The optional `trainLhosteModel.m` experiment completed locally without reported errors. The uploaded `lhoste_only_model.mat`, `lhoste_only_subject_metrics.csv`, and `lhoste_only_group_summary.csv` were checked for model settings, held-out subjects, sample counts, and aggregation consistency.

- **Training subjects (14):** AB68, AB69, AB71, AB72, AB75, AB76, AB78, AB79, AB80, AB81, HL02, HL03, HL04, HL05.
- **Held-out subjects (5):** AB70, AB73, AB77, AB82, HL01.
- **Training:** 49,680 Lhoste windows, 204 kinematic features, 1,000-ms causal history and lambda = 0.0001. The saved coefficient matrix is 204 × 2 with finite values.
- **Paired evaluation:** post-K5 E07 anchors: 33,182 SM anchors and 41,926 Transparent anchors, identical across the three methods. These anchor totals are **not** used to weight the group average: each of the five subjects contributes equally.

| Method | SM circular phase MAE | Transparent circular phase MAE | SM within ±10 pp | Transparent within ±10 pp |
| --- | ---: | ---: | ---: | ---: |
| Camargo Original | 10.52 pp | 12.70 pp | 52.96% | 47.58% |
| Camargo + K5 | **5.99 pp** | 8.85 pp | **81.24%** | 67.73% |
| Lhoste-only | 6.49 pp | **7.98 pp** | 78.95% | **73.35%** |

For SM, Camargo + K5 had a 0.50-pp lower mean MAE than Lhoste-only and performed better in three of five test subjects. For Transparent, Lhoste-only had a 0.87-pp lower mean MAE and performed better in four of five test subjects. Relative to Camargo Original, Camargo + K5 reduced mean phase MAE by 43.1% in SM and 30.3% in Transparent on these five subjects.

**Interpretation:** A brief, individual temporal calibration makes a natural-walking-trained model competitive with direct target-domain training in this exploratory comparison. This does **not** establish statistical equivalence or superiority with five held-out subjects. The strategies also use different training information: Lhoste-only is trained on other exoskeleton users (49,680 windows), while Camargo + K5 uses a model fitted to 122,400 sampled Camargo windows and adapts with five cycles per leg from each test user. The Lhoste-only training selection is performed per continuous run; the Camargo training selection is performed per trial. No p-values or out-of-cohort clinical generalization claims are made. The controller-defined full-cycle phase is an offline label, not an independently measured anatomical phase.

All original E01–E12 results remain unchanged. The comparison contributes a secondary research result suitable for the Discussion, rather than requiring additional transition-prediction analyses.
