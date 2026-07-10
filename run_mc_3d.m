%% GUIDANCE AND NAVIGATION - MONTE CARLO (3D CPN)
% This script runs MC trials for CPN_3D with method (A) wind modeling.
% Wind is modeled as a Gauss-Markov process with sampled std and length scale.

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

% Gauss-Markov gust parameter ranges for MC sampling.
% sigma_gust is the steady-state gust std [m/s] in x,y,z.
% L_gust is the turbulence length scale [m] in x,y,z.
gust_sigma_min  = [2; 2; 1];
gust_sigma_max  = [12; 12; 5];
gust_L_min      = [100; 100; 50];
gust_L_max      = [400; 400; 200];

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
    L_gust = sample_uniform(gust_L_min, gust_L_max);

    V_ref = mean(cell2mat(M_V_vec_3D));
    v_gust = generate_gauss_markov_gust(sigma_gust, L_gust, V_ref, dt_3D, n_sim);

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

%% Local helper functions
function x = sample_uniform(a, b)
    x = a + (b - a) .* rand(size(a));
end

function w = generate_gauss_markov_gust(sigma, Lc, V_ref, dt, n_sim)
% Generate 3D gust profile using an axis-wise first-order Gauss-Markov process.
% sigma: 3x1 steady-state std [m/s], Lc: 3x1 length scale [m].

    sigma = sigma(:);
    Lc = Lc(:);

    V_ref = max(V_ref, 1e-6);
    tau = max(Lc ./ V_ref, 1e-6);

    beta = exp(-dt ./ tau);
    q = sigma .* sqrt(1 - beta.^2);

    w = zeros(3, n_sim);
    w(:, 1) = sigma .* randn(3, 1);

    for t = 1:(n_sim - 1)
        w(:, t+1) = beta .* w(:, t) + q .* randn(3, 1);
    end
end
