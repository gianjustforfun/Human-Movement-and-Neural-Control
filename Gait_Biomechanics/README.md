# Gait Biomechanics

Biomechanical comparison of **natural walking** and two lower-limb exoskeleton control modes: **State-Machine (SM)** and **Transparent**.

This self-contained MATLAB analysis uses the original **Camargo (2021)** and **Lhoste (2026)** datasets to examine gait-cycle timing, hip and knee kinematics, and joint motion before a controller-defined stance-to-swing transition.

## Research Questions

- How do gait timing and hip/knee kinematics differ across natural walking, SM, and Transparent walking?
- How do SM and Transparent kinematics differ during the 300 ms preceding the controller's stance-to-swing transition?

Comparisons between Camargo and Lhoste are **descriptive** because the datasets differ in participants and experimental protocols. The SM–Transparent comparison is paired within the same 19 Lhoste participants.

## MATLAB Files

| File | Purpose |
| --- | --- |
| `main.m` | Runs data preparation, biomechanical summaries, and plots |
| `readCamargo.m` | Reads Camargo gait markers and inverse-kinematics data |
| `readLhoste.m` | Reads Lhoste joint kinematics and controller states |
| `eventAlignedAnalysis.m` | Analyzes the 300 ms before controller transition `1→2` |

Only **`main.m`** needs to be run. The other `.m` files are functions called by it.

## Setup and Execution

Download the original datasets separately. Keep the folder names shown below, because `main.m` looks for them **one level above** `Gait_Biomechanics`.

```text
Human-Movement-and-Neural-Control/
├── Gait_Biomechanics/
│   ├── main.m
│   ├── readCamargo.m
│   ├── readLhoste.m
│   ├── eventAlignedAnalysis.m
│   ├── README.md
│   └── results/
├── DATASET - CAMARGO (2021)/
└── DATASET - LHOSTE (2026)/
```

Open `Gait_Biomechanics/main.m` in MATLAB and click **Run**. The first run reads the source files and saves subject-level summaries in `results/cache/`. Later runs reuse these cache files. To rebuild the analysis from the original datasets, set `redo = true` near the beginning of `main.m`.

The analysis uses standard MATLAB functions, including `boxchart` (MATLAB R2020a or later). No machine-learning toolbox is required.

The raw datasets are **not included** in this repository and are never modified by the scripts.

## Methods

### Gait-Cycle Analysis

**Camargo:** Heel strikes and toe-offs are extracted from `gcLeft` and `gcRight` marker signals. Valid cycles run from one heel strike to the next and include exactly one toe-off. Hip and knee angles come from inverse-kinematics data.

**Lhoste:** A gait cycle is defined by successive controller transitions `3→0`, with the expected intermediate state sequence `0→1→2→3`. Transition `1→2` separates controller-defined stance and swing. Joint angles are resampled at 200 Hz, preserving the recorded sign conventions.

For both datasets, joint angular velocities are computed with a causal five-sample backward difference (25 ms). Valid cycles are interpolated to 101 points over 0–100% of the gait cycle.

Scalar measures are summarized within each subject before group means are calculated, so every subject receives equal weight.

### Pre-Transition Kinematics

`eventAlignedAnalysis.m` examines hip and knee angles and angular velocities from **−300 to −5 ms** before each valid Lhoste controller `1→2` transition, at 5-ms intervals. The event at 0 ms is only an alignment reference, not an input sample.

Events are averaged within each leg, the two legs contribute equally within each subject, and the group mean gives equal weight to the 19 participants.

The shaded regions in the event-aligned plots show ±1 between-subject standard deviation. Differences are reported as **Transparent minus SM** for paired participants.

## Results

### Gait Timing and Range of Motion

| Condition | Subjects | Valid cycles | Cycle duration (s) | Stance (%) | Hip ROM (°) | Knee ROM (°) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Natural walking (Camargo) | 22 | 14,942 | 1.16 | 57.63 | 45.38 | 59.81 |
| State-Machine | 19 | 1,467 | 2.92 | 62.21 | 46.71 | 68.81 |
| Transparent | 19 | 1,497 | 3.51 | 60.12 | 60.53 | 56.74 |

Values are group means of **subject-level medians**.

Within the Lhoste cohort, Transparent relative to SM showed:

- **+0.59 s** cycle duration (longer in 18/19 participants).
- **+13.82°** hip ROM (greater in 18/19 participants).
- **−12.06°** knee ROM (lower in 16/19 participants).
- **−2.09 percentage points** stance duration.

### Kinematics Preceding the Controller Transition

Across the paired Lhoste participants:

- Hip angle at −300 ms was **6.17° higher** in Transparent than SM (19/19 participants).
- Knee angle change from −300 to −5 ms was **0.93° higher** in Transparent (17/19 participants).

These are descriptive differences in controller-aligned kinematics, not evidence of accurate prospective event prediction.

## Generated Outputs

The scripts create `results/` automatically.

| Output | Contents |
| --- | --- |
| `subject_summary.csv` | Timing and ROM by subject and condition |
| `group_summary.csv` | Equally weighted group summaries |
| `01_cycle_timing.png` | Subject-level timing boxplots |
| `02_joint_angles.png` | Hip and knee angle profiles over the gait cycle |
| `03_angular_velocities.png` | Hip and knee angular-velocity profiles |
| `04_transparent_minus_SM.png` | Differences between group-average kinematic profiles |
| `event_aligned/01_event_aligned_kinematics.png` | Pre-transition mean kinematics ±1 SD |
| `event_aligned/02_paired_differences.png` | Paired participant-level differences |
| `event_aligned/subject_metrics.csv` | Pre-transition measurements by participant and condition |
| `event_aligned/paired_differences.csv` | Paired SM–Transparent comparisons |
| `event_aligned/group_metrics.csv` | Group-level pre-transition summaries |

## Interpretation and Limitations

- The Camargo and Lhoste cohorts, acquisition protocols, and event definitions are different. Cross-dataset comparisons should not be interpreted as causal effects of exoskeleton use.
- The Lhoste `1→2` transition is a **controller-state event**, not independently measured anatomical toe-off.
- Signed angular velocity should not be confused with absolute movement speed.
- Retrospective event alignment describes motion before known controller transitions; it does not demonstrate real-time prediction.
- The analysis is descriptive. No hypothesis tests or claims of statistical significance are made.

**Main finding:** SM and Transparent walking show distinct temporal and joint-kinematic patterns, including measurable differences immediately before controller-defined stance-to-swing transitions.
