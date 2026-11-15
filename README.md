# Human Movement and Neural Control

**Biomechanical Analysis and Personalized Gait-Phase Estimation in Natural and Exoskeleton Walking**

MATLAB project developed for **BME 574 – Human Movement and Neural Control** at the University of Illinois Chicago (UIC).

This project investigates how lower-limb exoskeleton control modes influence human gait biomechanics and whether models trained on natural walking can be adapted to estimate gait phase and anticipate controller transitions during exoskeleton-assisted walking.

The analysis combines two complementary components:

1. **Gait Biomechanics:** Comparison of natural walking, State-Machine (SM), and Transparent exoskeleton control.
2. **Gait Phase Estimation:** Transfer of a natural-walking-trained model to exoskeleton walking, personalized temporal calibration, and advance transition prediction.

## Research Objectives

The project addresses three main questions:

- How do gait timing and joint kinematics differ between natural walking and different exoskeleton control modes?
- Can a gait-phase model trained on natural walking be transferred to exoskeleton walking using short subject-specific calibration?
- Can personalized phase and kinematic information predict upcoming stance-to-swing controller transitions?

The broader objective is to investigate how biomechanical analysis and personalized modeling can contribute to adaptive exoskeleton control strategies.

## Datasets

Two publicly described research datasets are used:

| Dataset | Walking condition | Subjects | Valid gait cycles |
| --- | --- | ---: | ---: |
| Camargo (2021) | Natural walking | 22 | 14,942 |
| Lhoste (2026) | State-Machine | 19 | 1,467 |
| Lhoste (2026) | Transparent | 19 | 1,497 |

The datasets contain lower-limb kinematic information, including hip and knee joint angles and angular velocities.

The Lhoste dataset also provides exoskeleton controller states used to identify gait-cycle boundaries and controller-defined stance-to-swing transitions.

**The original datasets are not included in this repository.** They must be obtained separately from their respective sources.

## Repository Structure

```text
Human-Movement-and-Neural-Control/
│
├── Gait_Biomechanics/
│   ├── main.m
│   ├── readCamargo.m
│   ├── readLhoste.m
│   ├── eventAlignedAnalysis.m
│   └── README.md
│
├── Gait_Phase_Estimation/
│   ├── main.m
│   ├── readCamargo.m
│   ├── readLhoste.m
│   ├── trainPhaseModel.m
│   ├── evaluatePhase.m
│   ├── predictTransitions.m
│   ├── trainLhosteModel.m
│   └── README.md
│
└── README.md
```

Each analysis is **self-contained** and has its own `main.m` script.

Both components read the original datasets independently and do not require intermediate outputs from one another.

For implementation details, see:

- [Gait Biomechanics](Gait_Biomechanics/README.md)
- [Gait Phase Estimation](Gait_Phase_Estimation/README.md)

## 1. Gait Biomechanics

This component investigates differences in gait patterns between natural walking and two exoskeleton control modes.

The analysis includes:

- Gait-cycle duration and stance/swing timing.
- Hip and knee joint angles.
- Joint angular velocities and range of motion (ROM).
- Differences between State-Machine and Transparent control.
- Event-aligned kinematics during the 300 ms preceding the controller-defined stance-to-swing transition.

### Main Findings

Compared with State-Machine walking, Transparent control showed:

| Biomechanical measure | Transparent − SM |
| --- | ---: |
| Gait-cycle duration | +0.59 s |
| Hip ROM | +13.82° |
| Knee ROM | −12.06° |
| Stance duration | −2.09 percentage points |

These findings highlight distinct biomechanical patterns associated with the two exoskeleton control modes.

The comparison between SM and Transparent is paired within the same 19 participants. Comparisons with natural walking are descriptive because the datasets were collected using different experimental protocols.

## 2. Gait Phase Estimation

This component investigates the transfer of gait-phase models from natural walking to exoskeleton-assisted walking.

A **Circular Phase Ridge Regression** model is trained on Camargo kinematics using causal hip and knee angle and angular-velocity features.

The analysis includes:

- Subject-disjoint training and validation.
- Personalized temporal calibration using 5 or 10 gait cycles per leg.
- Comparison between direct model transfer and temporally calibrated estimation.
- Lhoste-only model training on separate participants.
- Advance controller-transition prediction at 100, 200, and 300 ms.
- Comparison between phase-based and additional kinematic prediction features.

### Main Findings

The Camargo phase model achieved a validation mean absolute error (MAE) of **2.38 phase percentage points**.

After five calibration cycles per leg, temporal calibration substantially improved phase estimation:

| Condition | Original MAE | Calibrated MAE |
| --- | ---: | ---: |
| State-Machine | 12.39 pp | **7.36 pp** |
| Transparent | 14.60 pp | **8.37 pp** |

A separate comparison on five unseen Lhoste subjects found that the calibrated Camargo model performed competitively with a model trained directly on exoskeleton data.

For advance controller-transition prediction at a 200-ms horizon, the personalized model with additional kinematic features achieved:

| Metric | State-Machine | Transparent |
| --- | ---: | ---: |
| Event recall | 31.37% | 19.03% |
| False alarms/device-minute | 25.93 | 24.61 |

Although personalized modeling improved event recall, the high false-alarm rates remain a major limitation for practical closed-loop control.

## Running the Project

The analyses are implemented in MATLAB.

1. Clone or download this repository.
2. Obtain the original Camargo (2021) and Lhoste (2026) datasets.
3. Place the dataset folders alongside the two analysis folders, using the names expected by the MATLAB scripts:

```text
Human-Movement-and-Neural-Control/
├── Gait_Biomechanics/
├── Gait_Phase_Estimation/
├── DATASET - CAMARGO (2021)/
└── DATASET - LHOSTE (2026)/
```

4. Open the desired `main.m` in MATLAB.
5. Run the script.

Each main script calls the required functions and saves generated outputs in its own `results/` directory.

Intermediate results may be cached to avoid repeating expensive computations. The original datasets remain unchanged.

The individual README files provide additional information about execution settings, data processing, and output files.

## Scientific Interpretation

The two components address different aspects of the same problem.

Biomechanical analysis identifies differences in gait timing and joint motion associated with exoskeleton control strategies.

Gait-phase estimation investigates whether a model trained on natural walking can be transferred to these different movement conditions through individual calibration.

Together, the findings suggest that **short subject-specific temporal calibration can substantially improve phase estimation across walking domains**, while accurate advance prediction of controller transitions remains challenging.

The study is exploratory and subject to several limitations:

- Natural walking and exoskeleton walking originate from different datasets and experimental protocols.
- Controller-defined stance-to-swing transitions are not independently measured anatomical toe-off events.
- Transition prediction was evaluated offline, not in a real-time closed-loop controller.
- High false-alarm rates limit the practical reliability of advance transition prediction.
- No claims of clinical effectiveness or statistical significance are made.

## References

The analysis is based on the following research datasets:

- **Camargo et al. (2021):** Lower-limb biomechanics during natural walking and other locomotor activities.
- **Lhoste et al. (2026):** Lower-limb robotic exoskeleton walking under different control modes.

Detailed methodology and dataset descriptions are provided in the respective analysis folders.

## Project Information

**Course:** BME 574 – Human Movement and Neural Control

**Institution:** University of Illinois Chicago

**Language:** MATLAB

**Project type:** Biomechanical analysis and personalized gait modeling

**Author:** Gianluigi D'Antonio

---

*This project was developed for academic research and educational purposes. The results are exploratory and are not intended for direct clinical or autonomous exoskeleton control applications.*
