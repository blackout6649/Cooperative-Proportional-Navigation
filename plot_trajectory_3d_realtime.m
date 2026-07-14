%% GUIDANCE AND NAVIGATION - REAL-TIME 3D TRAJECTORY VISUALIZATION
% This function plots missile trajectories in real-time as an animation,
% showing missiles moving in 3D space like a video.
% 
% Usage:
%   plot_trajectory_3d_realtime(data, playback_speed, save_video, video_fps, video_quality)
%
% Inputs:
%   data            - struct containing simulation data (from CPN_3D)
%                     Required fields: RM, RT, t_vec, m
%   playback_speed  - (optional) playback speed factor. 1 = real-time, 0.5 = half-speed, 2 = double-speed
%                     Default: 1 (real-time relative to simulation)
%   save_video      - (optional) if true, saves animation as video file. Default: false
%   video_fps       - (optional) video output frame rate in fps. Default: 15 (good compromise)
%   video_quality   - (optional) 'low' (720p), 'medium' (1080p), 'high' (1440p). Default: 'medium'

function plot_trajectory_3d_realtime(data, playback_speed, save_video, video_fps, video_quality)

if nargin < 2
    playback_speed = 1;  % real-time by default
end
if nargin < 3
    save_video = false;
end
if nargin < 4
    video_fps = 15;  % 15 fps is good compromise for file size
end
if nargin < 5
    video_quality = 'medium';  % 'low' (720p), 'medium' (1080p), 'high' (1440p)
end

% Video resolution settings
video_resolutions = struct(...
    'low', [1280, 720], ...
    'medium', [1920, 1080], ...
    'high', [2560, 1440]);

if isfield(video_resolutions, video_quality)
    video_resolution = video_resolutions.(video_quality);
else
    video_resolution = video_resolutions.medium;
    warning('Unknown quality "%s". Using medium.', video_quality);
end

% Extract data
RM    = data.RM;
RT    = data.RT;
t_vec = data.t_vec;
m     = data.m;
dt    = t_vec(2) - t_vec(1);  % time step

% Target position (stationary)
xT = RT(1, 1);
yT = RT(2, 1);
zT = RT(3, 1);

% Create figure
fig = figure('Name', 'Real-Time 3D Trajectory Visualization', 'Color', 'w', 'NumberTitle', 'off');

% If saving video, set figure size to match video resolution to preserve aspect ratio
if save_video
    set(fig, 'Units', 'pixels', 'Position', [100, 100, video_resolution(1), video_resolution(2)]);
end

ax  = axes(fig);
hold(ax, 'on')
axis(ax, 'equal')
grid(ax, 'on')
ax.GridAlpha      = 0.25;
ax.MinorGridAlpha = 0.10;
ax.Box            = 'on';
ax.LineWidth      = 1.0;
ax.FontSize       = 11;

% Determine axis limits based on trajectory extent
x_all = [];
y_all = [];
z_all = [];
for i = 1:m
    x_all = [x_all, RM{i}(1, :)];
    y_all = [y_all, RM{i}(2, :)];
    z_all = [z_all, RM{i}(3, :)];
end

x_lim = [min([x_all, xT]) - 500, max([x_all, xT]) + 500];
y_lim = [min([y_all, yT]) - 500, max([y_all, yT]) + 500];
z_lim = [min([z_all, zT]) - 500, max([z_all, zT]) + 500];

xlim(ax, x_lim)
ylim(ax, y_lim)
zlim(ax, z_lim)

% Labels and title
xlabel(ax, 'x [m]', 'FontSize', 12)
ylabel(ax, 'y [m]', 'FontSize', 12)
zlabel(ax, 'z [m] (vertical)', 'FontSize', 12)
title(ax, 'Real-Time 3D CPN Trajectory Animation', sprintf('%d cooperative missiles', m), ...
    'FontSize', 14, 'FontWeight', 'bold')

% View angle and rotate3d
view(ax, 45, 25)
camproj(ax, 'perspective')
rotate3d(ax, 'on')

% Color scheme for missiles
colors = [
    0.2 0.4 0.8;    % blue
    0.8 0.2 0.4;    % red
    0.2 0.8 0.4;    % green
];
if m > 3
    % Generate additional colors if more than 3 missiles
    for i = 4:m
        colors = [colors; hsv2rgb([mod(i/m, 1), 0.7, 0.9])];
    end
end

% Initialize plot objects for each missile
missile_lines = gobjects(m, 1);
missile_points = gobjects(m, 1);
missile_trails = gobjects(m, 1);

for i = 1:m
    % Line for trail (will be updated each frame)
    missile_trails(i) = plot3(ax, nan, nan, nan, ...
        'Color', colors(i, :), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    
    % Scatter point for current missile position
    missile_points(i) = scatter3(ax, nan, nan, nan, 100, colors(i, :), ...
        'filled', 'o', 'HandleVisibility', 'off', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Label for missile
    text(ax, nan, nan, nan, sprintf('M%d', i), 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
end

% Target marker
h_tgt = plot3(ax, xT, yT, zT, 'p', 'MarkerSize', 20, ...
    'MarkerFaceColor', [0.85 0.10 0.10], 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2, 'DisplayName', 'Target');

% Legend
h_path_proxy = plot3(ax, nan, nan, nan, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 2.5, ...
    'DisplayName', 'Missile trail');
h_missile_proxy = scatter3(ax, nan, nan, nan, 100, [0.5 0.5 0.5], 'filled', 'o', ...
    'DisplayName', 'Missile position');
legend(ax, [h_path_proxy, h_missile_proxy, h_tgt], ...
    {'Missile trail', 'Missile position', 'Target'}, ...
    'Location', 'northeast', 'FontSize', 10, 'Box', 'on')

n_sim = length(t_vec);

% Setup video writer if requested
video_writer = [];
if save_video
    % Use MP4 (MPEG-4) with H.264 compression - much smaller file size
    video_filename = 'CPN_3D_RT_animation.mp4';
    video_writer = VideoWriter(video_filename, 'MPEG-4');
    video_writer.FrameRate = video_fps;
    video_writer.Quality = 95;  % Quality 0-100, higher = better but larger file
    open(video_writer);
    fprintf('Recording animation to: %s\n', video_filename);
    fprintf('  Video FPS: %d (for %.1f sec video, ~%d frames will be captured)\n', ...
        video_fps, t_vec(end)/playback_speed, round(t_vec(end)/playback_speed * video_fps));
    fprintf('  Resolution: %dx%d (%s)\n', video_resolution(1), video_resolution(2), video_quality);
end

% Timing for smooth playback
tic;
frame_count = 0;
% Only sample frames at the rate we'll actually record them
frame_skip = max(1, round(n_sim / (t_vec(end)/playback_speed * video_fps)));

% Display information
fprintf('Starting real-time animation...\n');
fprintf('Playback speed: %.2fx\n', playback_speed);
fprintf('Simulation duration: %.2f seconds\n', t_vec(end));
fprintf('Use mouse to rotate view. Press Ctrl+C to stop.\n\n');

%% Animation loop
for t = 1:frame_skip:n_sim
    % Update missile positions and trails
    for i = 1:m
        % Current position
        x_t = RM{i}(1, t);
        y_t = RM{i}(2, t);
        z_t = RM{i}(3, t);
        
        % Update scatter point (current position)
        set(missile_points(i), 'XData', x_t, 'YData', y_t, 'ZData', z_t);
        
        % Trail from start to current time
        trail_x = RM{i}(1, 1:t);
        trail_y = RM{i}(2, 1:t);
        trail_z = RM{i}(3, 1:t);
        
        set(missile_trails(i), 'XData', trail_x, 'YData', trail_y, 'ZData', trail_z);
    end
    
    % Update title with current time
    current_time = t_vec(t);
    title(ax, sprintf('Real-Time 3D CPN Trajectory Animation (t = %.2f s)', current_time), ...
        sprintf('%d cooperative missiles', m), 'FontSize', 14, 'FontWeight', 'bold')
    
    % Redraw
    drawnow;
    
    % Frame timing control
    frame_count = frame_count + 1;
    elapsed_time = toc;
    sim_time = current_time / playback_speed;
    
    % Control playback speed
    if elapsed_time < sim_time
        pause(0.001);  % Small pause to avoid busy-waiting
    end
    
    % Save frame to video if requested
    if ~isempty(video_writer)
        frame = getframe(fig);
        writeVideo(video_writer, frame);
    end
end

if ~isempty(video_writer)
    close(video_writer);
    fprintf('\nAnimation saved to: %s\n', video_filename);
end

fprintf('Animation complete!\n');
hold(ax, 'off')

end
