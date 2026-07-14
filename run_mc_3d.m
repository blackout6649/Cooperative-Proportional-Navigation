%% GUIDANCE AND NAVIGATION - MONTE CARLO (3D CPN)
% This script runs MC trials for CPN_3D 
% Wind is modeled as a Gauss-Markov process with sampled std and length scale.

clear; clc;

%% Fixed scenario parameters
%rng(1); % reproducible MC sequence

N_mc            = 1000;

m_3D            = 3;
N_3D            = 3;
K_gain_3D       = 40;
A_max_3D        = 35;
R_hit_3D        = 1;

tend_3D         = 30;
dt_3D           = 0.001;
t_vec           = 0:dt_3D:tend_3D;
n_sim           = length(t_vec);

T_x0_3D         = [0, 0, 0];

%% VISUALIZATION OPTIONS
% Set to true to enable real-time 3D trajectory animation (only works with N_mc = 1)
plot_realtime   = false;
% Playback speed: 1 = real-time, 0.5 = half-speed, 2 = double-speed
rt_playback_speed = 1;
% Save animation to video file (only if plot_realtime = true)
rt_save_video   = false;
% Video output frame rate in fps (lower = smaller file). Default: 15 fps
% For 30 sec video: 15 fps = 450 frames, 10 fps = 300 frames, 5 fps = 150 frames
rt_video_fps    = 15;
% Video quality/resolution: 'low' (720p), 'medium' (1080p), 'high' (1440p)
rt_video_quality = 'medium';

% Gauss-Markov gust parameter ranges for MC sampling.
% sigma_gust is the steady-state gust std [m/s] in x,y,z.
% L_gust is the turbulence length scale [m] in x,y,z.
gust_sigma_min  = [2; 2; 1];
gust_sigma_max  = [12; 12; 5];
gust_L_min      = [100; 100; 50];
gust_L_max      = [400; 400; 200];

% Measurement noise std ranges for RM [m] and VM [m/s].
sigma_RM_min    = [0; 0; 0];
sigma_RM_max    = [0.3; 0.3; 0.3];
sigma_VM_min    = [0; 0; 0];
sigma_VM_max    = [0.05; 0.05; 0.05];

%% Outputs tracked per run
hit_vec         = false(1, N_mc);
N_CPN_end       = zeros(m_3D, N_mc);
r_min_end       = zeros(m_3D, N_mc);
miss_distance_min_vec = nan(1, N_mc);
t_end_vec       = zeros(1, N_mc);

N_CPN_ts        = cell(m_3D, N_mc);
epsilon_ts      = cell(m_3D, N_mc);
a_com_ts        = cell(m_3D, N_mc);
a_uncapped_ts   = cell(m_3D, N_mc);
v_gust_ts       = cell(1, N_mc);
t_run_vec       = cell(1, N_mc);

tau_m_vec       = zeros(1, N_mc);
tau_f_vec       = zeros(1, N_mc);
tau_sk_vec      = zeros(1, N_mc);
impact_time_vec = nan(1, N_mc);
t_hit_mat       = nan(m_3D, N_mc);
hit_flags_mat   = false(m_3D, N_mc);
t_hit_first_vec = nan(1, N_mc);
t_hit_last_vec  = nan(1, N_mc);
t_hit_spread_vec = nan(1, N_mc);
hit_angle_vec   = nan(m_3D, N_mc, 3);
sigma_RM_vec    = zeros(3, N_mc);
sigma_VM_vec    = zeros(3, N_mc);

wb = waitbar(0, 'Running Monte Carlo simulation...');

%% Monte Carlo loop
for k = 1:N_mc

    % Random inputs from MC parameters.md
    tau_m_3D = sample_uniform(0.3, 1);
    tau_f_3D = sample_uniform(0.01, 0.1);

    tau_sk_3D = sample_uniform(0.01, 0.2);

    M_x0_vec_3D     = cell([m_3D 1]);
    M_V_vec_3D      = cell([m_3D 1]);
    M_gamma0_vec_3D = cell([m_3D 1]);

    % initial positions: sample each missile an annulus (R_min <= r <= R_max) on the x-y plane (z = 0) so they start at least R_min from the target.
    R_max_3D = 3000;
    R_min_3D  = 2000;

    for i = 1:m_3D
        ang0              = 2 * pi * (i - 1 + rand);
        r0                = sqrt(R_min_3D^2 + rand * (R_max_3D^2 - R_min_3D^2));
        M_x0_vec_3D{i}     = [r0 * cos(ang0), r0 * sin(ang0), 0];
        M_V_vec_3D{i}      = sample_uniform(100, 150);
        % heading [phi, theta, psi]: elevation theta ~ 90 deg (straight up) +-10 deg
        M_gamma0_vec_3D{i} = [sample_uniform(60, 120), ...
                              sample_uniform(80, 100), ...
                              sample_uniform(60, 120)];
    end

    % Additional random input: gust realization for this run.
    sigma_gust = sample_uniform(gust_sigma_min, gust_sigma_max);
    L_gust = sample_uniform(gust_L_min, gust_L_max);
    sigma_RM = sample_uniform(sigma_RM_min, sigma_RM_max);
    sigma_VM = sample_uniform(sigma_VM_min, sigma_VM_max);

    V_ref = mean(cell2mat(M_V_vec_3D));
    v_gust = generate_gauss_markov_gust(sigma_gust, L_gust, V_ref, dt_3D, n_sim);


    data = CPN_3D(m_3D, N_3D, K_gain_3D, A_max_3D, R_hit_3D, tau_m_3D, ...
        tau_f_3D, tau_sk_3D, tend_3D, dt_3D, T_x0_3D, M_x0_vec_3D, ...
        M_V_vec_3D, M_gamma0_vec_3D, v_gust, sigma_RM, sigma_VM);

    if isfield(data, 'hit_flags')
        run_success = all(logical(data.hit_flags));
    else
        run_success = logical(data.hit);
    end

    hit_vec(k)   = run_success;
    t_end_vec(k) = data.t_vec(end);
    t_run_vec{k} = data.t_vec;
    sigma_RM_vec(:, k) = sigma_RM(:);
    sigma_VM_vec(:, k) = sigma_VM(:);

    if run_success
        impact_time_vec(k) = data.t_vec(end);
    end

    if isfield(data, 'hit_flags')
        hit_flags_mat(:, k) = logical(data.hit_flags(:));
    end

    if isfield(data, 't_hit')
        t_hit_k = data.t_hit(:);
        t_hit_mat(:, k) = t_hit_k;

        valid_hit_k = ~isnan(t_hit_k);
        if any(valid_hit_k)
            t_hit_first_vec(k) = min(t_hit_k(valid_hit_k));
            t_hit_last_vec(k) = max(t_hit_k(valid_hit_k));
            t_hit_spread_vec(k) = t_hit_last_vec(k) - t_hit_first_vec(k);
        end
    end

    n_k = length(data.t_vec);
    v_gust_ts{k} = data.v_gust(:, 1:n_k);

    miss_distance_min_run = inf;

    for i = 1:m_3D
        N_CPN_end(i, k) = data.N_CPN{i}(end);
        r_min_end(i, k) = min(data.r_go{i});
        N_CPN_ts{i, k} = data.N_CPN{i};
        epsilon_ts{i, k} = data.epsilon{i};
        a_com_ts{i, k} = data.a_com{i};
        a_uncapped_ts{i, k} = data.a_uncapped{i};

        r_go_true_i = sqrt(sum((data.RT - data.RM{i}) .^ 2, 1));
        miss_distance_min_run = min(miss_distance_min_run, min(r_go_true_i));

        if run_success
            gamma_hit = data.gammaM{i}(:, end);
            hit_angle_vec(i, k, :) = [gamma_hit(3), gamma_hit(2), gamma_hit(1)];
        end
    end

    miss_distance_min_vec(k) = miss_distance_min_run;

    tau_m_vec(k)  = tau_m_3D;
    tau_f_vec(k)  = tau_f_3D;
    tau_sk_vec(k) = tau_sk_3D;

    waitbar(k / N_mc, wb, sprintf('Running Monte Carlo simulation... %d/%d', k, N_mc));
end

if isgraphics(wb)
    close(wb);
end

% for a single run, show the 3D trajectory of that run
if N_mc == 1
    if plot_realtime
        % Real-time animation
        fprintf('\n========================================\n');
        fprintf('Plotting REAL-TIME 3D Trajectory Animation\n');
        fprintf('========================================\n');
        plot_trajectory_3d_realtime(data, rt_playback_speed, rt_save_video, rt_video_fps, rt_video_quality);
    else
        % Static 3D plot
        plot_trajectory_3d(data);
    end
end

%% Aggregate MC results
results = struct;
results.N_mc         = N_mc;
results.hit          = hit_vec;
results.miss_distance_min = miss_distance_min_vec;
results.N_CPN_end    = N_CPN_end;
results.N_CPN_ts     = N_CPN_ts;
results.epsilon_ts   = epsilon_ts;
results.a_com_ts     = a_com_ts;
results.a_uncapped_ts = a_uncapped_ts;
results.v_gust_ts    = v_gust_ts;
results.t_run_vec    = t_run_vec;
results.r_min_end    = r_min_end;
results.t_end        = t_end_vec;
results.impact_time  = impact_time_vec;
results.t_hit        = t_hit_mat;
results.hit_flags_missile = hit_flags_mat;
results.t_hit_first  = t_hit_first_vec;
results.t_hit_last   = t_hit_last_vec;
results.t_hit_spread = t_hit_spread_vec;
results.hit_angle    = hit_angle_vec;
results.tau_m        = tau_m_vec;
results.tau_f        = tau_f_vec;
results.tau_sk       = tau_sk_vec;
results.sigma_RM     = sigma_RM_vec;
results.sigma_VM     = sigma_VM_vec;
results.A_max        = A_max_3D;

results.P_hit        = mean(hit_vec);
results.hit_probability = results.P_hit;
results.mean_r_min   = mean(r_min_end, 2);
results.mean_N_CPN   = mean(N_CPN_end, 2);

save('mc_results_3d.mat', 'results', '-v7.3');

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
