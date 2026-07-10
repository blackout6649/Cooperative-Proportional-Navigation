%% GUIDANCE AND NAVIGATION - FINAL PROJECT NUMERICAL SIMULATION

%% 3-D trajectory plot with time color gradient
function plot_trajectory_3d(data, save_fig)
if nargin < 2
    save_fig = 0;
end

RM    = data.RM;
RT    = data.RT;
t_vec = data.t_vec;
m     = data.m;

% target is stationary -> single point (x-y horizontal plane, z vertical)
xT = RT(1, 1);
yT = RT(2, 1);
zT = RT(3, 1);

fig = figure('Color', 'w');
ax  = axes(fig);
hold(ax, 'on')

% ---- missile trajectories as time-colored gradient lines ----
h_init = gobjects(1, 1);
for i = 1:m
    x = RM{i}(1, :);
    y = RM{i}(2, :);
    z = RM{i}(3, :);
    c = t_vec;

    % colored line: color encodes simulation time along the path
    surface(ax, [x; x], [y; y], [z; z], [c; c], ...
        'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.5, ...
        'HandleVisibility', 'off');

    % initial position marker
    hi = plot3(ax, x(1), y(1), z(1), 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.0, 'HandleVisibility', 'off');
    if i == 1
        h_init = hi;
    end
end

% ---- stationary target: single marker ----
h_tgt = plot3(ax, xT, yT, zT, 'p', 'MarkerSize', 20, ...
    'MarkerFaceColor', [0.85 0.10 0.10], 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2, 'HandleVisibility', 'off');

% ---- colorbar (time legend) ----
% single-hue gradient (light -> dark blue) encoding time
n_shades = 256;
c_lo     = [0.78 0.86 0.98];   % light blue (t start)
c_hi     = [0.03 0.12 0.45];   % dark blue  (t end)
time_cmap = [linspace(c_lo(1), c_hi(1), n_shades)', ...
             linspace(c_lo(2), c_hi(2), n_shades)', ...
             linspace(c_lo(3), c_hi(3), n_shades)'];
colormap(ax, time_cmap)
cb = colorbar(ax);
cb.Label.String     = 'time [s]';
cb.Label.FontSize   = 11;
cb.Label.FontWeight = 'bold';
if t_vec(end) > t_vec(1)
    clim(ax, [t_vec(1), t_vec(end)])
end

% ---- legend via proxy handles ----
h_path_proxy = plot3(ax, nan, nan, nan, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 2.5);
legend(ax, [h_path_proxy, h_init, h_tgt], ...
    {'Missile path (color = time)', 'Missile initial position', 'Target'}, ...
    'Location', 'northeast', 'FontSize', 10, 'Box', 'on')

% ---- formatting ----
xlabel(ax, 'x [m]', 'FontSize', 12)
ylabel(ax, 'y [m]', 'FontSize', 12)
zlabel(ax, 'z [m]  (vertical)', 'FontSize', 12)
title(ax, '3D CPN Trajectories', sprintf('%d cooperative missiles', m), ...
    'FontSize', 14, 'FontWeight', 'bold')

grid(ax, 'on')
ax.GridAlpha      = 0.25;
ax.MinorGridAlpha = 0.10;
ax.Box            = 'on';
ax.LineWidth      = 1.0;
ax.FontSize       = 11;
axis(ax, 'equal')
axis(ax, 'tight')
view(ax, 45, 25)
camproj(ax, 'perspective')
rotate3d(ax, 'on')
hold(ax, 'off')

if save_fig
    set(fig, 'units', 'pix', 'pos', [0, 0, 1920, 1080])
    saveas(fig, "3D_GNC_sim_" + m + "_missiles" + ".png")
end

end

% inputs:
    % data      - output struct from CPN_3D (uses data.RM, data.RT,
    %             data.t_vec and data.m)
    % save_fig  - save figure as png, 1 = save, 0 = no save (default 0)

% notes:
    % x-y is the horizontal plane, z is the vertical axis (gravity direction).
    % The missile path color encodes simulation time (see colorbar).
