# Gait Phase Estimation

Personalized gait-phase estimation and advance stance-to-swing transition prediction using natural walking and lower-limb exoskeleton kinematics.

This MATLAB project investigates whether a model trained on **Camargo (2021)** natural walking data can be transferred to the **Lhoste (2026)** exoskeleton dataset through short subject-specific temporal calibration.

The analysis also compares direct training on exoskeleton data and evaluates the contribution of additional kinematic features to controller-transition prediction.

## Research Questions

- Can a gait-phase model trained on natural walking generalize to exoskeleton walking?
- Can a short individual calibration improve phase estimation without retraining the original model?
- Is personalized temporal calibration competitive with training directly on exoskeleton data?
- Can causal kinematic information improve advance stance-to-swing transition prediction?

The study considers two exoskeleton control conditions: **State-Machine (SM)** and **Transparent**.

## MATLAB Files

| File | Purpose |
| --- | --- |
| `main.m` | Runs the complete analysis |
| `readCamargo.m` | Reads Camargo gait cycles and constructs causal kinematic features |
| `readLhoste.m` | Reads Lhoste controller states, cycles, and continuous kinematic histories |
| `trainPhaseModel.m` | Trains and validates the Camargo circular phase model |
| `evaluatePhase.m` | Evaluates original and temporally calibrated phase estimation |
| `predictTransitions.m` | Predicts controller transitions and compares personalized models |
| `trainLhosteModel.m` | Trains a Lhoste-only model and compares it on unseen subjects |

Only **`main.m`** needs to be executed. All other `.m` files are MATLAB functions called by the main script.

## Setup and Execution

Download the original Camargo and Lhoste datasets separately.

The dataset folders must be located one level above `Gait_Phase_Estimation`, as shown below.

```text
Human-Movement-and-Neural-Control/
├── Gait_Phase_Estimation/
│   ├── main.m
│   ├── readCamargo.m
│   ├── readLhoste.m
│   ├── trainPhaseModel.m
│   ├── evaluatePhase.m
│   ├── predictTransitions.m
│   ├── trainLhosteModel.m
│   ├── README.md
│   └── results/
│       └── cache/
├── DATASET - CAMARGO (2021)/
└── DATASET - LHOSTE (2026)/
```

Open `Gait_Phase_Estimation/main.m` in MATLAB and click **Run**.

The first execution prepares the datasets, trains the models, evaluates phase estimation, and generates the transition-prediction results.

Intermediate data are saved in `results/cache/` to avoid repeating computationally expensive operations.

The following settings control which analyses are recomputed:

| Setting | Purpose |
| --- | --- |
| `redo` | Rebuild the dataset caches |
| `redoModel` | Retrain the Camargo phase model |
| `redoPhase` | Recalculate phase predictions |
| `redoTransitions` | Recalculate transition-prediction results |
| `runLhosteComparison` | Enable the additional Lhoste-only comparison |
| `redoLhoste` | Retrain and reevaluate the Lhoste-only model |

All settings can normally remain `false`, except `runLhosteComparison = true` when the additional comparison is desired.

Missing results are generated automatically. Existing results can be reused on subsequent executions.

The original datasets are **not included** in this repository and are never modified.

## Methods

### 1. Data Preparation

**Camargo:** Heel strikes and toe-offs are identified from gait-marker signals. Valid cycles are assembled into continuous kinematic histories containing hip and knee joint angles and angular velocities.

**Lhoste:** Gait cycles are defined using controller transitions `3→0`, while the controller transition `1→2` is used as the stance-to-swing reference.

Lhoste joint kinematics are resampled at 200 Hz, and recording chronology is preserved across separate acquisition segments.

The `1→2` event represents a **controller-defined transition**, not an independently measured anatomical toe-off.

### 2. Camargo Gait-Phase Model

A Circular Phase Ridge Regression model estimates gait phase from four kinematic signals:

- Hip joint angle
- Knee joint angle
- Hip angular velocity
- Knee angular velocity

Features are constructed causally using past and current measurements.

| Parameter | Value |
| --- | --- |
| Feature sampling frequency | 50 Hz |
| Selected history window | 1,000 ms |
| Features | 204 |
| Regression | Circular Phase Ridge |
| Selected regularization λ | 0.0001 |
| Training/validation split | 17/5 subjects |

The model predicts the sine and cosine of gait phase, allowing phase to be reconstructed while accounting for its circular nature.

Training and validation subjects are strictly separated. Hyperparameters are selected using Camargo data only.

The final model is refitted using all 22 Camargo subjects, without using Lhoste data.

### 3. Personalized Temporal Calibration

The trained Camargo model is applied to Lhoste kinematics under two approaches:

**Original:** The model is transferred directly without temporal adjustment.

**TimeRescaled:** Kinematic history windows are individually rescaled according to the participant's gait-cycle duration.

The first **5 or 10 complete cycles per leg** are used for calibration.

The personalized temporal scale is computed as:

\[
s = \frac{\operatorname{median}(T_{\mathrm{Lhoste,cal}})}
{\operatorname{median}(T_{\mathrm{Camargo}})}
\]

Angular velocities are transformed consistently with the temporal rescaling.

The Camargo model coefficients remain fixed, and performance is evaluated exclusively on subsequent cycles.

Both calibration lengths are also compared on common post-calibration samples to ensure a paired evaluation.

### 4. Advance Controller-Transition Prediction

The objective is to predict the upcoming controller transition `1→2` using only past and current kinematic measurements.

Five approaches are compared:

| Method | Description |
| --- | --- |
| Original | Prediction using the uncalibrated Camargo phase |
| TimeRescaled | Prediction using temporally calibrated phase |
| Phase Rate | Estimated phase and causal phase-advancement rate |
| 3-feature Ridge | Personalized time-to-transition correction using phase and timing features |
| 7-feature Ridge | Same correction with four additional hip/knee kinematic features |

The personalized regressions use only the initial 5 or 10 calibration cycles per leg.

Both Ridge models use identical calibration samples, regression targets, and evaluation anchors.

Prediction horizons of **100, 200, and 300 ms** are investigated.

Performance is quantified using event recall, precision, and false alarms per device-minute.

Events occurring after calibration are used exclusively as offline references for scoring. They are never provided as predictive inputs.

### 5. Lhoste-Only Training Comparison

An additional Circular Phase Ridge model is trained directly on Lhoste kinematics.

| Parameter | Value |
| --- | --- |
| Training subjects | 14 |
| Unseen test subjects | 5 |
| Training windows | 49,680 |
| Features | 204 |
| History | 1,000 ms |
| Ridge λ | 0.0001 |
| Training conditions | SM and Transparent |

The model is trained on both exoskeleton conditions and evaluated separately for each condition.

The five test participants are excluded entirely from model training.

Three approaches are compared on identical post-K5 evaluation samples:

1. Camargo Original
2. Camargo + K5 calibration
3. Lhoste-only

The comparison evaluates the benefit of personalized transfer learning relative to training directly on exoskeleton data.

## Results

### Dataset Preparation

| Dataset | Subjects | Valid cycles |
| --- | ---: | ---: |
| Camargo | 22 | 14,942 |
| Lhoste SM | 19 | 1,467 |
| Lhoste Transparent | 19 | 1,497 |

Camargo processing produced **835,934 causal kinematic windows** across 680 trials.

### Camargo Model Validation

| Metric | Result |
| --- | ---: |
| Circular Phase MAE | **2.38 percentage points** |
| Predictions within ±10 pp | **98.4%** |
| Selected history | 1,000 ms |
| Selected λ | 0.0001 |

These values were obtained using subject-disjoint validation.

### Temporal Calibration

Circular phase MAE after five calibration cycles per leg:

| Condition | Original | TimeRescaled | Relative reduction |
| --- | ---: | ---: | ---: |
| State-Machine | 12.39 pp | **7.36 pp** | 40.6% |
| Transparent | 14.60 pp | **8.37 pp** | 42.7% |

The results demonstrate that short temporal calibration substantially improves phase estimation without retraining the Camargo model.

Increasing calibration from 5 to 10 cycles per leg produced only small additional changes, with no consistent improvement across both conditions.

### Lhoste-Only Comparison

Circular phase MAE on five unseen Lhoste subjects:

| Model | SM | Transparent |
| --- | ---: | ---: |
| Camargo Original | 10.52 pp | 12.70 pp |
| Camargo + K5 | **5.99 pp** | 8.85 pp |
| Lhoste-only | 6.49 pp | **7.98 pp** |

Camargo + K5 performed slightly better in SM, while Lhoste-only performed better in Transparent.

The results suggest that a short subject-specific calibration can make a model trained on natural walking competitive with direct exoskeleton-domain training.

The comparison is exploratory and does not establish statistical equivalence between the approaches.

### Advance Transition Prediction

Event recall with K=5 and a 200-ms prediction horizon:

| Method | SM recall | Transparent recall |
| --- | ---: | ---: |
| Original | 9.53% | 5.49% |
| TimeRescaled | 19.07% | 11.00% |
| Phase Rate | 20.91% | 10.80% |
| 3-feature Ridge | 26.64% | 13.83% |
| 7-feature Ridge | **31.37%** | **19.03%** |

Personalized Ridge correction improves transition prediction, particularly when additional joint kinematics are included.

However, the best model still produces:

| Metric | SM | Transparent |
| --- | ---: | ---: |
| Precision | 32.61% | 21.42% |
| False alarms/device-minute | 25.93 | 24.61 |

These performances are insufficient for reliable autonomous exoskeleton control.

### Kinematic Feature Ablation

The contribution of four additional kinematic features was evaluated by comparing the personalized 3-feature and 7-feature Ridge models.

At a 200-ms prediction horizon:

| Condition | Calibration | Recall improvement | Additional false alarms/min |
| --- | --- | ---: | ---: |
| SM | K5 | +4.73 pp | +3.09 |
| SM | K10 | +5.69 pp | +2.71 |
| Transparent | K5 | +5.19 pp | +2.23 |
| Transparent | K10 | +1.70 pp | +1.14 |

The additional kinematic features improve mean recall at 200 ms, but also increase false alarms.

The benefit is not consistent across all participants or prediction horizons.

## Generated Outputs

All results are saved automatically inside `results/`.

| Output | Contents |
| --- | --- |
| `phase_model.mat` | Trained Camargo circular phase model |
| `phase_validation_models.csv` | Model-selection and validation results |
| `phase_validation_subjects.csv` | Subject-level validation metrics |
| `calibration_summary.csv` | Individual K5/K10 calibration parameters |
| `phase_subject_metrics.csv` | Subject-level phase-estimation results |
| `phase_group_summary.csv` | Group phase-estimation metrics |
| `phase_calibration.png` | Original vs calibrated phase estimation |
| `transition_subject_metrics.csv` | Subject-level transition-prediction metrics |
| `transition_group_summary.csv` | Group transition-prediction results |
| `transition_training_QC.csv` | Personalized training summaries |
| `transition_ablation_H200.csv` | Kinematic feature-ablation results |
| `transition_recall.png` | Recall comparison across prediction methods |
| `transition_ablation.png` | Recall improvement and false-alarm trade-off |
| `lhoste_only_model.mat` | Trained Lhoste-only phase model |
| `lhoste_only_subject_metrics.csv` | Held-out subject comparison |
| `lhoste_only_group_summary.csv` | Group Lhoste-only comparison |
| `lhoste_only_comparison.png` | Comparison of the three training strategies |

The `results/cache/` directory contains intermediate data and individual phase predictions used to accelerate subsequent MATLAB executions.

Cached files and large intermediate outputs do not need to be uploaded to GitHub.

## Interpretation and Limitations

- Natural walking and exoskeleton walking belong to different acquisition domains, with different kinematics and experimental conditions.
- The Lhoste `1→2` transition is a controller-defined event and must not be interpreted as independently measured anatomical toe-off.
- Temporal calibration improves phase estimation, but does not guarantee accurate prediction of future transitions.
- The Lhoste-only comparison uses five unseen test subjects and does not establish statistical equivalence or superiority.
- The Camargo and Lhoste-only approaches have different training-data requirements.
- Additional kinematic features improve mean recall at 200 ms but also generate more false alarms.
- Transition prediction is evaluated retrospectively and has not been validated in a real-time closed-loop exoskeleton controller.
- The reported group metrics give equal weight to individual subjects rather than pooling all observations.
- No claims of statistical significance or clinical effectiveness are made.

## Main Findings

**Personalized temporal calibration substantially improves the transfer of gait-phase estimation from natural walking to exoskeleton walking.**

Using only five calibration cycles per leg, the Camargo-trained model achieves phase-estimation performance competitive with a model trained directly on exoskeleton data.

Additional personalized kinematic information improves advance controller-transition recall, but high false-alarm rates remain a major limitation.

These findings support the potential of short individual calibration for gait-phase estimation while highlighting the challenges of reliable prospective transition prediction.
