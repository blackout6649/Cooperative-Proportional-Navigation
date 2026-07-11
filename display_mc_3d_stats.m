%% DISPLAY MONTE CARLO STATISTICS (3D CPN)
% Plots mean, std x1, std x3, and min/max vs time across all MC runs.

clear; clc; close all;

if ~isfile('mc_results_3d.mat')
    error('mc_results_3d.mat not found. Run run_mc_3d.m first.');
end

S = load('mc_results_3d.mat', 'results');
results = S.results;

fprintf('MC runs: %d\n', results.N_mc);
fprintf('All-missiles hit probability: %.3f\n', results.hit_probability);

m = size(results.N_CPN_ts, 1);

% Common time grid (all runs share dt and start at t = 0).
dt = infer_dt(results.t_run_vec);
n_max = max_length(results.t_run_vec);
t_grid = (0:n_max-1) * dt;

%% N_CPN vs time (per missile - subplots)
figure('Name', 'N_{CPN} - Per Missile');
for i = 1:m
    subplot(m, 1, i);
    X = stack_cell_row(results.N_CPN_ts, i, n_max);
    plot_time_stats_ax(gca, t_grid, X, sprintf('N_{CPN} - Missile %d', i), 'N_{CPN}');
end

%% epsilon vs time (per missile - subplots)
figure('Name', '\epsilon - Per Missile');
for i = 1:m
    subplot(m, 1, i);
    X = stack_cell_row(results.epsilon_ts, i, n_max);
    plot_time_stats_ax(gca, t_grid, X, sprintf('\\epsilon - Missile %d', i), '\epsilon [sec]');
end

%% N_CPN vs time (all missiles combined)
X = stack_cell_all(results.N_CPN_ts, n_max);
plot_time_stats(t_grid, X, 'N_{CPN} - All Missiles Combined', 'N_{CPN}');

%% epsilon vs time (all missiles combined)
X = stack_cell_all(results.epsilon_ts, n_max);
plot_time_stats(t_grid, X, '\epsilon - All Missiles Combined', '\epsilon [sec]');

%% v_gust vs time (per axis)
axis_names = {'x', 'y', 'z'};
for ax = 1:3
    X = stack_gust_axis(results.v_gust_ts, ax, n_max);
    plot_time_stats(t_grid, X, sprintf('v_{gust} - %s axis', axis_names{ax}), ...
        'v_{gust} [m/s]');
end

%% a_com vs time (single figure: all missiles, all axes, in g)
X = stack_cell_all_axes(results.a_com_ts, n_max) ./ 9.8;
plot_time_stats(t_grid, X, 'a_{com} - All Missiles and Axes Combined', 'a_{com} [g]', false, results.A_max);

%% a_uncapped vs time (single figure: all missiles, all axes, in g)
X = stack_cell_all_axes(results.a_uncapped_ts, n_max) ./ 9.8;
plot_time_stats(t_grid, X, 'a_{uncapped} - All Missiles and Axes Combined', 'a_{uncapped} [g]', false, results.A_max);

%% Miss-distance histogram
plot_miss_distance_hist(results.miss_distance_min);

%% Time of impact histogram
plot_impact_time_hist(results.impact_time, results.N_mc);

%% First-vs-last missile hit-time difference histogram
plot_hit_spread_hist(results, m);

%% Hit angle histograms
plot_hit_angle_histograms(results.hit_angle, results.N_mc);

%% Save all figures (optional, via pop-up)
save_all_figures_prompt();

%% Local helper functions
function dt = infer_dt(t_run_vec)
    dt = 1;
    for k = 1:numel(t_run_vec)
        t = t_run_vec{k};
        if numel(t) >= 2
            dt = t(2) - t(1);
            return;
        end
    end
end

function n = max_length(t_run_vec)
    n = 0;
    for k = 1:numel(t_run_vec)
        n = max(n, numel(t_run_vec{k}));
    end
    if n == 0; n = 1; end
end

function X = stack_cell_row(C, row_idx, n_max)
% Build [n_runs x n_max] matrix; shorter runs padded with NaN.

    n_runs = size(C, 2);
    X = nan(n_runs, n_max);
    for k = 1:n_runs
        v = C{row_idx, k};
        if ~isempty(v)
            v = v(:).';
            X(k, 1:numel(v)) = v;
        end
    end
end

function X = stack_gust_axis(v_gust_ts, ax, n_max)
% Build [n_runs x n_max] matrix for one gust axis; padded with NaN.

    n_runs = numel(v_gust_ts);
    X = nan(n_runs, n_max);
    for k = 1:n_runs
        v = v_gust_ts{k};
        if ~isempty(v)
            row = v(ax, :);
            X(k, 1:numel(row)) = row;
        end
    end
end

function X = stack_cell_row_axis(C, row_idx, ax, n_max)
% Build [n_runs x n_max] matrix from one missile row and one vector axis.

    n_runs = size(C, 2);
    X = nan(n_runs, n_max);
    for k = 1:n_runs
        v = C{row_idx, k};
        if ~isempty(v)
            row = v(ax, :);
            X(k, 1:numel(row)) = row;
        end
    end
end

function plot_time_stats(t, X, ttl, ylab, show_three_sigma, accel_limit)
% Plot mean, mean +/- 1 std, optional mean +/- 3 std, and min/max vs time.

    if nargin < 5
        show_three_sigma = true;
    end
    if nargin < 6
        accel_limit = [];
    end

    mu    = mean(X, 1, 'omitnan');
    s     = std(X, 0, 1, 'omitnan');
    x_min = min(X, [], 1, 'omitnan');
    x_max = max(X, [], 1, 'omitnan');

    valid = any(~isnan(X), 1);
    t = t(valid);
    mu = mu(valid); s = s(valid);
    x_min = x_min(valid); x_max = x_max(valid);

    figure('Name', ttl);
    hold on; grid on;

    % min/max envelope
    fill_between(t, x_min, x_max, [0.90 0.90 0.90], 'min/max');

    if show_three_sigma
        % +/- 3 sigma band
        fill_between(t, mu - 3*s, mu + 3*s, [0.75 0.83 0.95], '\pm 3\sigma');
    end

    % +/- 1 sigma band
    fill_between(t, mu - s, mu + s, [0.55 0.70 0.90], '\pm 1\sigma');

    if ~isempty(accel_limit) && isfinite(accel_limit)
        yline(accel_limit, 'r--', 'LineWidth', 1.4, 'DisplayName', '+a_{max}');
        yline(-accel_limit, 'r--', 'LineWidth', 1.4, 'DisplayName', '-a_{max}');
    end

    plot(t, mu, 'k-', 'LineWidth', 1.8, 'DisplayName', 'mean');

    title(ttl);
    xlabel('t [s]');
    ylabel(ylab);
    legend('Location', 'best');
    hold off;
end

function plot_time_stats_ax(ax, t, X, ttl, ylab, show_three_sigma, accel_limit)
% Same as plot_time_stats but draws into a given axes handle (for subplots).

    if nargin < 6
        show_three_sigma = true;
    end
    if nargin < 7
        accel_limit = [];
    end

    mu    = mean(X, 1, 'omitnan');
    s     = std(X, 0, 1, 'omitnan');
    x_min = min(X, [], 1, 'omitnan');
    x_max = max(X, [], 1, 'omitnan');

    valid = any(~isnan(X), 1);
    t = t(valid);
    mu = mu(valid); s = s(valid);
    x_min = x_min(valid); x_max = x_max(valid);

    hold(ax, 'on'); grid(ax, 'on');

    % min/max envelope
    fill_between_ax(ax, t, x_min, x_max, [0.90 0.90 0.90], 'min/max');

    if show_three_sigma
        fill_between_ax(ax, t, mu - 3*s, mu + 3*s, [0.75 0.83 0.95], '\pm 3\sigma');
    end

    fill_between_ax(ax, t, mu - s, mu + s, [0.55 0.70 0.90], '\pm 1\sigma');

    if ~isempty(accel_limit) && isfinite(accel_limit)
        yline(ax, accel_limit, 'r--', 'LineWidth', 1.4, 'DisplayName', '+a_{max}');
        yline(ax, -accel_limit, 'r--', 'LineWidth', 1.4, 'DisplayName', '-a_{max}');
    end

    plot(ax, t, mu, 'k-', 'LineWidth', 1.8, 'DisplayName', 'mean');

    title(ax, ttl);
    xlabel(ax, 't [s]');
    ylabel(ax, ylab);
    legend(ax, 'Location', 'best');
    hold(ax, 'off');
end

function fill_between_ax(ax, t, y_lo, y_hi, rgb, name)
% Shaded band into a specific axes handle.

    t = t(:).';
    y_lo = y_lo(:).';
    y_hi = y_hi(:).';

    xx = [t, fliplr(t)];
    yy = [y_lo, fliplr(y_hi)];

    fill(ax, xx, yy, rgb, 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
        'DisplayName', name);
end

function fill_between(t, y_lo, y_hi, rgb, name)
% Shaded band between y_lo and y_hi.

    t = t(:).';
    y_lo = y_lo(:).';
    y_hi = y_hi(:).';

    xx = [t, fliplr(t)];
    yy = [y_lo, fliplr(y_hi)];

    fill(xx, yy, rgb, 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
        'DisplayName', name);
end

function X = stack_cell_all(C, n_max)
% Build [(n_runs*m) x n_max] matrix pooling all missiles; padded with NaN.

    m = size(C, 1);
    blocks = cell(m, 1);
    for i = 1:m
        blocks{i} = stack_cell_row(C, i, n_max);
    end
    X = vertcat(blocks{:});
end

function X = stack_cell_all_axis(C, ax, n_max)
% Build [(n_runs*m) x n_max] matrix pooling all missiles for one vector axis.

    m = size(C, 1);
    blocks = cell(m, 1);
    for i = 1:m
        blocks{i} = stack_cell_row_axis(C, i, ax, n_max);
    end
    X = vertcat(blocks{:});
end

function X = stack_cell_all_axes(C, n_max)
% Build [(n_runs*m*3) x n_max] matrix pooling all missiles and all axes.

    m = size(C, 1);
    blocks = cell(m, 1);
    for i = 1:m
        row_x = stack_cell_row_axis(C, i, 1, n_max);
        row_y = stack_cell_row_axis(C, i, 2, n_max);
        row_z = stack_cell_row_axis(C, i, 3, n_max);
        blocks{i} = [row_x; row_y; row_z];
    end
    X = vertcat(blocks{:});
end

function X = stack_cell_all_magnitudes(C, n_max)
% Build [(n_runs*m) x n_max] matrix of vector magnitudes pooled over all missiles.

    m = size(C, 1);
    blocks = cell(m, 1);
    for i = 1:m
        n_runs = size(C, 2);
        X_i = nan(n_runs, n_max);
        for k = 1:n_runs
            v = C{i, k};
            if ~isempty(v)
                mag = sqrt(sum(v.^2, 1));
                X_i(k, 1:numel(mag)) = mag;
            end
        end
        blocks{i} = X_i;
    end
    X = vertcat(blocks{:});
end

function plot_miss_distance_hist(miss_distance_min)
% Histogram of closest-approach miss distance across all MC runs.

    miss_distance_min = miss_distance_min(:);

    figure('Name', 'Miss Distance Histogram');
    histogram(miss_distance_min);
    grid on;

    mu = mean(miss_distance_min);
    s  = std(miss_distance_min);

    title(sprintf('Miss Distance Distribution (mean = %.4g, std = %.4g)', mu, s));
    xlabel('miss\_distance\_min [m]');
    ylabel('count');
end

function plot_impact_time_hist(impact_time, n_mc)
% Histogram of impact times across hit runs.

    impact_time = impact_time(~isnan(impact_time));
    n_bins = max(5, round(n_mc / 10));

    figure('Name', 'Impact Time Histogram');

    if isempty(impact_time)
        histogram(nan);
        title('Impact Time Distribution (no hits)');
    else
        histogram(impact_time, n_bins);
        mu = mean(impact_time);
        s  = std(impact_time);
        title(sprintf('Impact Time Distribution (mean = %.4g, std = %.4g)', mu, s));
    end

    grid on;
    xlabel('time of impact [s]');
    ylabel('count');
end

function plot_hit_spread_hist(results, m)
% Histogram of per-scenario hit-time spread: last missile hit minus first missile hit.

    if ~isfield(results, 't_hit')
        warning(['Cannot plot first-last hit-time histogram: results.t_hit is missing. ', ...
                 'Need per-missile hit times as an m x N_mc matrix in results.t_hit.']);
        return;
    end

    t_hit = results.t_hit;
    if ~ismatrix(t_hit) || size(t_hit, 1) ~= m
        warning('Cannot plot first-last hit-time histogram: results.t_hit must be m x N_mc.');
        return;
    end

    all_hit_runs = all(~isnan(t_hit), 1);
    spread = max(t_hit(:, all_hit_runs), [], 1) - min(t_hit(:, all_hit_runs), [], 1);
    n_bins = max(5, round(results.N_mc / 10));

    figure('Name', 'First-to-Last Hit Time Difference Histogram');

    if isempty(spread)
        histogram(nan);
        title('First-to-Last Hit Time Difference (no all-hit runs)');
    else
        histogram(spread, n_bins);
        mu = mean(spread);
        s  = std(spread);
        title(sprintf('First-to-Last Hit Time Difference (mean = %.4g, std = %.4g)', mu, s));
    end

    grid on;
    xlabel('\Delta t_{hit} = t_{last} - t_{first} [s]');
    ylabel('count of scenarios');
end

function plot_hit_angle_histograms(hit_angle, n_mc)
% Polar histograms of hit angles [psi, theta, phi] for each missile and combined.

    angle_names = {'\psi', '\theta', '\phi'};
    edges_deg = 0:30:360;
    m = size(hit_angle, 1);

    for i = 1:m
        figure('Name', sprintf('Hit Angle Polar Histogram - Missile %d', i));

        for j = 1:3
            subplot(3, 1, j);
            x = squeeze(hit_angle(i, :, j));
            x = rad2deg(x(~isnan(x)));
            x = wrap_angle_deg(x);

            if isempty(x)
                polarhistogram('BinEdges', deg2rad(edges_deg), 'BinCounts', zeros(1, numel(edges_deg) - 1));
                title(sprintf('%s Missile %d (no hits)', angle_names{j}, i));
            else
                polarhistogram(deg2rad(x), 'BinEdges', deg2rad(edges_deg));
                title(sprintf('%s Missile %d (mean = %.4g deg, std = %.4g deg)', ...
                    angle_names{j}, i, mean(x), std(x)));
            end
        end
    end

    figure('Name', 'Hit Angle Polar Histogram - All Missiles Combined');
    for j = 1:3
        subplot(3, 1, j);
        x = reshape(hit_angle(:, :, j), [], 1);
        x = rad2deg(x(~isnan(x)));
        x = wrap_angle_deg(x);

        if isempty(x)
            polarhistogram('BinEdges', deg2rad(edges_deg), 'BinCounts', zeros(1, numel(edges_deg) - 1));
            title(sprintf('%s All Missiles (no hits)', angle_names{j}));
        else
            polarhistogram(deg2rad(x), 'BinEdges', deg2rad(edges_deg));
            title(sprintf('%s All Missiles (mean = %.4g deg, std = %.4g deg)', ...
                angle_names{j}, mean(x), std(x)));
        end
    end
end

function x = wrap_angle_deg(x)
% Wrap angles in degrees to [0, 360).

    x = mod(x, 360);
end

function save_all_figures_prompt()
% Pop-up dialog asking whether to save all open figures, and if yes, where.

    answer = questdlg('Save all figures to disk?', 'Save Figures', 'Yes', 'No', 'No');
    if ~strcmp(answer, 'Yes')
        return;
    end

    save_dir = uigetdir(pwd, 'Select folder to save figures');
    if isequal(save_dir, 0)
        return;  % user cancelled
    end

    figs = findall(0, 'Type', 'figure');
    if isempty(figs)
        msgbox('No figures to save.', 'Info');
        return;
    end

    for i = 1:numel(figs)
        fig = figs(i);
        fig_name = get(fig, 'Name');
        if isempty(fig_name)
            fig_name = sprintf('figure_%d', fig.Number);
        end
        % sanitize filename
        fig_name = regexprep(fig_name, '[\\/:*?"<>|]', '_');

        % save as .fig and .png
        savefig(fig, fullfile(save_dir, [fig_name '.fig']));
        exportgraphics(fig, fullfile(save_dir, [fig_name '.png']), 'Resolution', 200);
    end

    msgbox(sprintf('Saved %d figures to:\n%s', numel(figs), save_dir), 'Done');
end
