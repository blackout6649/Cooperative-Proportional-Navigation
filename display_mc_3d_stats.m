%% DISPLAY MONTE CARLO STATISTICS (3D CPN)
% Plots mean, std x1, std x3, and min/max vs time across all MC runs.

clear; clc; close all;

if ~isfile('mc_results_3d.mat')
    error('mc_results_3d.mat not found. Run run_mc_3d.m first.');
end

S = load('mc_results_3d.mat', 'results');
results = S.results;

fprintf('MC runs: %d\n', results.N_mc);
fprintf('Hit probability: %.6f\n', results.hit_probability);

m = size(results.N_CPN_ts, 1);

% Common time grid (all runs share dt and start at t = 0).
dt = infer_dt(results.t_run_vec);
n_max = max_length(results.t_run_vec);
t_grid = (0:n_max-1) * dt;

%% N_CPN vs time (per missile)
for i = 1:m
    X = stack_cell_row(results.N_CPN_ts, i, n_max);
    plot_time_stats(t_grid, X, sprintf('N_{CPN} - Missile %d', i), 'N_{CPN}');
end

%% epsilon vs time (per missile)
for i = 1:m
    X = stack_cell_row(results.epsilon_ts, i, n_max);
    plot_time_stats(t_grid, X, sprintf('\\epsilon - Missile %d', i), '\epsilon [sec]');
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

%% R_hit histogram
plot_R_hit_hist(results.R_hit);

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

function plot_time_stats(t, X, ttl, ylab)
% Plot mean, mean +/- 1 std, mean +/- 3 std, and min/max vs time.

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

    % +/- 3 sigma band
    fill_between(t, mu - 3*s, mu + 3*s, [0.75 0.83 0.95], '\pm 3\sigma');

    % +/- 1 sigma band
    fill_between(t, mu - s, mu + s, [0.55 0.70 0.90], '\pm 1\sigma');

    plot(t, mu, 'k-', 'LineWidth', 1.8, 'DisplayName', 'mean');

    title(ttl);
    xlabel('t [s]');
    ylabel(ylab);
    legend('Location', 'best');
    hold off;
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

function plot_R_hit_hist(R_hit)
% Histogram of R_hit across all MC runs.

    R_hit = R_hit(:);

    figure('Name', 'R_{hit} Histogram');
    histogram(R_hit);
    grid on;

    mu = mean(R_hit);
    s  = std(R_hit);

    title(sprintf('R_{hit} Distribution (mean = %.4g, std = %.4g)', mu, s));
    xlabel('R_{hit} [m]');
    ylabel('count');
end
