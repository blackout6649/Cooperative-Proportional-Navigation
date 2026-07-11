# Cooperative Proportional Navigation (CPN) — 3D Monte Carlo Simulation

Final project for the Principles of Guidance & Navigation course.

## Members
Ilan Elzam
Erik Gurfinkel
Tamar Itzhaki
Margarita Birchanski

## Overview

This project simulates a **cooperative proportional navigation (CPN)** guidance law for multiple missiles engaging a stationary target in 3D space. A Monte Carlo (MC) framework evaluates robustness under randomized initial conditions, wind gusts, measurement noise, and system lags.


## File Structure

| File | Description |
|------|-------------|
| `CPN_3D.m` | Core 3D CPN guidance simulation function |
| `CPN_2D.m` | 2D variant (reference/legacy) |
| `run_mc_3d.m` | Monte Carlo driver script — runs N scenarios and saves results |
| `display_mc_3d_stats.m` | Loads results and generates all statistical plots |
| `plot_trajectory_3d.m` | 3D trajectory visualization with time-color gradient |
| `save_current_results_v73.m` | Utility to re-save workspace results in MAT v7.3 |
| `MC parameters.md` | Documentation of MC input ranges and output fields |
| `project_main_code.m` | Top-level entry point / scratch script |

## How to Run

### 1. Run Monte Carlo

```matlab
run_mc_3d
```

This executes `N_mc` scenarios (configurable at top of file), saves `mc_results_3d.mat`.

### 2. Display Results

```matlab
display_mc_3d_stats
```

Generates all statistical plots. At the end, a pop-up dialog offers to save all figures (`.fig` + `.png`) to a user-selected directory.

### 3. Single-Run Trajectory

```matlab
% Set N_mc = 1 in run_mc_3d.m, then run it
% Or call CPN_3D directly and pass output to:
plot_trajectory_3d(data)
```

## Simulation Parameters

- **Missiles**: 3 cooperative missiles
- **Navigation constant**: N = 3
- **Cooperative gain**: K = 40
- **Acceleration limit**: 35 g
- **Hit radius**: 1 m
- **Time step**: 1 ms
- **Max simulation time**: 30 s

## Monte Carlo Randomization

Each scenario randomly samples:
- Mechanical, filter, and seeker lag time constants
- Missile initial positions (annulus around target)
- Missile speeds and heading angles
- Wind gust parameters (Gauss–Markov model)
- Position and velocity measurement noise

See `MC parameters.md` for full ranges.

## Key Outputs

- Hit probability (all missiles must hit)
- Miss distance distribution
- Impact time and first-to-last hit spread
- Time-varying N_CPN, epsilon, acceleration commands
- Hit angle distributions (polar histograms)
- Wind gust profiles

## Dependencies

- MATLAB (R2020b or later recommended for `exportgraphics`)
- No external toolboxes required
