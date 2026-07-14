# Real-Time 3D Trajectory Visualization Guide

## Overview
A new real-time 3D trajectory plotting function has been created to visualize missile trajectories as an animation, showing the missiles moving through space like a video.

## Files Created/Modified

### New File
- **`plot_trajectory_3d_realtime.m`** - Real-time animation function for 3D trajectories

### Modified File
- **`run_mc_3d.m`** - Added visualization options and conditional plotting logic

## How to Use

### Quick Start
To enable real-time visualization in `run_mc_3d.m`:

1. **Set `N_mc = 1`** - The real-time animation only works with a single run
2. **Set `plot_realtime = true`** - Enable the real-time plotting flag
3. **Run the script** - You'll see the animated trajectory

### Example Configuration

In `run_mc_3d.m`, near the top of the file, you'll find these parameters:

```matlab
%% Fixed scenario parameters
N_mc            = 1;        % <-- Set to 1 for real-time visualization

%% VISUALIZATION OPTIONS
plot_realtime   = true;             % Enable real-time animation
rt_playback_speed = 1;              % Playback speed (1 = real-time, 0.5 = half, 2 = double)
rt_save_video   = false;            % Set to true to save animation as AVI file
```

### Playback Speed Options

| Value | Speed |
|-------|-------|
| 0.5 | Half-speed (slower, better for detail observation) |
| 1.0 | Real-time (default, simulates actual elapsed time) |
| 2.0 | Double-speed (faster viewing) |

### Saving the Animation

To save the visualization as a video file:

```matlab
rt_save_video = true;
```

The video will be saved as `CPN_3D_RT_animation.avi` in your working directory.

## Features

### Animation Display Shows
- **Missile Trails** - Colored lines showing the path traveled
- **Current Positions** - Colored dots showing current missile positions
- **Target** - Large red star at the origin
- **Real-time Elapsed Time** - Updated in the plot title
- **Missile Labels** - M1, M2, M3, etc. for identification

### Interactivity
- **Rotate View** - Click and drag the mouse to rotate the 3D view
- **Zoom** - Use mouse scroll wheel to zoom in/out
- **Pan** - Use middle mouse button to pan

## Function Signature

```matlab
plot_trajectory_3d_realtime(data, playback_speed, save_video)

Inputs:
  data            - Simulation data struct from CPN_3D (required)
  playback_speed  - Playback speed factor (optional, default: 1)
  save_video      - Boolean to save video (optional, default: false)
```

## Example Workflows

### Scenario 1: Quick Visual Check
```matlab
N_mc = 1;
plot_realtime = true;
rt_playback_speed = 1;
rt_save_video = false;
% Run run_mc_3d.m
```

### Scenario 2: Detailed Inspection (Slow Motion)
```matlab
N_mc = 1;
plot_realtime = true;
rt_playback_speed = 0.25;      % Quarter-speed for detailed observation
rt_save_video = false;
% Run run_mc_3d.m
```

### Scenario 3: Save Presentation Video
```matlab
N_mc = 1;
plot_realtime = true;
rt_playback_speed = 0.5;       % Half-speed for smooth video
rt_save_video = true;          % Save as AVI file
% Run run_mc_3d.m
% Video saved as: CPN_3D_RT_animation.avi
```

## Notes

- Real-time animation only activates when `N_mc = 1` (single run mode)
- If `plot_realtime = false` and `N_mc = 1`, the original static 3D plot displays instead
- The animation respects the frame rate display (~30 fps) regardless of simulation time step
- You can stop the animation anytime by pressing Ctrl+C in MATLAB

## Troubleshooting

### Animation plays too fast/slow
Adjust `rt_playback_speed`:
- Too fast? Decrease value (e.g., 0.5)
- Too slow? Increase value (e.g., 2.0)

### Video file not created
- Ensure `rt_save_video = true`
- Check that the working directory is writable
- Video file should appear as `CPN_3D_RT_animation.avi`

### Can't rotate the view
- Make sure you click on the 3D plot area first
- `rotate3d` is automatically enabled - click and drag with mouse
