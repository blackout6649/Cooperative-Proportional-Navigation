%% GUIDANCE AND NAVIGATION - MONTE CARLO (3D CPN)
% This script runs MC trials for CPN_3D with method (A) wind modeling.
% Wind is modeled as zero-mean Gaussian velocity samples with sampled std.

clear; clc;

%% Fixed scenario parameters
rng(1); % reproducible MC sequence

N_mc            = 10;

m_3D            = 3;
N_3D            = 3;
K_gain_3D       = 40;
A_max_3D        = 100;
R_hit_3D        = 1;

tend_3D         = 30;
dt_3D           = 0.001;
t_vec           = 0:dt_3D:tend_3D;
n_sim           = length(t_vec);

T_x0_3D         = [0, 0, 2500];

% Gaussian gust std ranges for MC sampling [m/s] in x,y,z.
gust_sigma_min  = [2; 2; 1];
gust_sigma_max  = [12; 12; 5];

%% Outputs tracked per run
hit_vec         = false(1, N_mc);
R_hit_vec       = zeros(1, N_mc);
N_CPN_end       = zeros(m_3D, N_mc);
r_min_end       = zeros(m_3D, N_mc);
t_end_vec       = zeros(1, N_mc);

N_CPN_ts        = cell(m_3D, N_mc);
epsilon_ts      = cell(m_3D, N_mc);
v_gust_ts       = cell(1, N_mc);
t_run_vec       = cell(1, N_mc);

tau_m_vec       = zeros(1, N_mc);
tau_f_vec       = zeros(1, N_mc);
tau_sk_vec      = zeros(1, N_mc);

wb = waitbar(0, 'Running Monte Carlo simulation...');

%% Monte Carlo loop
for k = 1:N_mc

    % Random inputs from MC parameters.md
    tau_m_3D = sample_uniform(0.005, 0.05);
    tau_f_3D = sample_uniform(0.01, 0.1);

    tau_sk_3D = sample_uniform(0.01, 0.2);

    M_x0_vec_3D     = cell([m_3D 1]);
    M_V_vec_3D      = cell([m_3D 1]);
    M_gamma0_vec_3D = cell([m_3D 1]);

    for i = 1:m_3D
        M_x0_vec_3D{i}     = sample_uniform(-500 .* ones(1, 3), 500 .* ones(1, 3));
        M_V_vec_3D{i}      = sample_uniform(100, 150);
        M_gamma0_vec_3D{i} = sample_uniform(60 .* ones(1, 3), 120 .* ones(1, 3));
    end

    % Additional random input: gust realization for this run.
    sigma_gust = sample_uniform(gust_sigma_min, gust_sigma_max);
    v_gust = generate_normal_gust(sigma_gust, n_sim);

    data = CPN_3D(m_3D, N_3D, K_gain_3D, A_max_3D, R_hit_3D, tau_m_3D, ...
        tau_f_3D, tau_sk_3D, tend_3D, dt_3D, T_x0_3D, M_x0_vec_3D, ...
        M_V_vec_3D, M_gamma0_vec_3D, v_gust, 0, 0);

    hit_vec(k)   = logical(data.hit);
    R_hit_vec(k) = data.R_hit;
    t_end_vec(k) = data.t_vec(end);
    t_run_vec{k} = data.t_vec;

    n_k = length(data.t_vec);
    v_gust_ts{k} = data.v_gust(:, 1:n_k);

    for i = 1:m_3D
        N_CPN_end(i, k) = data.N_CPN{i}(end);
        r_min_end(i, k) = min(data.r_go{i});
        N_CPN_ts{i, k} = data.N_CPN{i};
        epsilon_ts{i, k} = data.epsilon{i};
    end

    tau_m_vec(k)  = tau_m_3D;
    tau_f_vec(k)  = tau_f_3D;
    tau_sk_vec(k) = tau_sk_3D;

    waitbar(k / N_mc, wb, sprintf('Running Monte Carlo simulation... %d/%d', k, N_mc));
end

if isgraphics(wb)
    close(wb);
end

%% Aggregate MC results
results = struct;
results.N_mc         = N_mc;
results.hit          = hit_vec;
results.R_hit        = R_hit_vec;
results.N_CPN_end    = N_CPN_end;
results.N_CPN_ts     = N_CPN_ts;
results.epsilon_ts   = epsilon_ts;
results.v_gust_ts    = v_gust_ts;
results.t_run_vec    = t_run_vec;
results.r_min_end    = r_min_end;
results.t_end        = t_end_vec;
results.tau_m        = tau_m_vec;
results.tau_f        = tau_f_vec;
results.tau_sk       = tau_sk_vec;

results.P_hit        = mean(hit_vec);
results.hit_probability = results.P_hit;
results.mean_r_min   = mean(r_min_end, 2);
results.mean_N_CPN   = mean(N_CPN_end, 2);

save('mc_results_3d.mat', 'results');

fprintf('MC runs: %d\n', N_mc);
fprintf('Hit probability: %.4f\n', results.P_hit);
fprintf('Saved results to mc_results_3d.mat\n');
for i = 1:m_3D
    fprintf('Missile %d: mean min range = %.3f m, mean final N_CPN = %.3f\n', ...
        i, results.mean_r_min(i), results.mean_N_CPN(i));
end

%% Local helper functions
function x = sample_uniform(a, b)
    x = a + (b - a) .* rand(size(a));
end

function w = generate_normal_gust(sigma, n_sim)
% Generate 3D gust profile with zero-mean Gaussian samples.
% sigma: 3x1 standard deviation [m/s] for x,y,z.

    sigma = sigma(:);
    w = sigma .* randn(3, n_sim);
end
