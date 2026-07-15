## MC SETUP

N_mc = 1000
m_3D = 3
N_3D = 3
K_gain_3D = 40
A_max_3D = 35 [g]
R_hit_3D = 1 [m]
tend_3D = 30 [sec]
dt_3D = 0.001 [sec]

## RANDOM INPUTS PER SCENARIO

tau_m = U(0.3, 1.0)            # mechanical lag time constant [sec]
tau_f = U(0.01, 0.1)              # filter lag time constant [sec]
tau_sk = U(0.01, 0.2)             # seeker lag time constant [sec]

M_x0_vec_3D{i} = [r0*cos(ang0), r0*sin(ang0), 0]    # missile initial position [m]
ang0 = 2*pi*(i - 1 + U(0,1))                         # angular sample (per missile)
r0 = sqrt(R_min_3D^2 + U(0,1)*(R_init_3D^2 - R_min_3D^2))
R_min_3D = 2000
R_init_3D = 3000

M_V_vec_3D{i} = U(100, 150)       # missile speed [m/sec]

M_gamma0_vec_3D{i} = [theta0, psi0]                  # initial heading [deg]
theta0 = U(80, 100)               # elevation around 90 deg (+/-10 deg)
psi0 = U(60, 120)

sigma_gust = U^3([2, 2, 1], [12, 12, 5])            # gust std per axis [m/sec]
L_gust = U^3([100, 100, 50], [400, 400, 200])       # Gauss-Markov length scale per axis [m]
v_gust(t) = first-order Gauss-Markov process         # 3 x n_sim gust profile [m/sec]

sigma_RM = U^3([0, 0, 0], [0.3, 0.3, 0.3])          # RM measurement-noise std [m]
sigma_VM = U^3([0, 0, 0], [0.05, 0.05, 0.05])       # VM measurement-noise std [m/sec]

## SUCCESS CRITERION

hit = 1 only if all missiles hit in the scenario
hit = 0 otherwise

## RESULTS FIELDS (mc_results_3d.mat -> results)

N_mc                      # number of MC scenarios
hit                       # 1 x N_mc, all-missiles-hit flag per scenario
hit_probability, P_hit    # mean(hit)

miss_distance_min         # 1 x N_mc, closest-approach distance per scenario [m]
r_min_end                 # m x N_mc, per-missile minimum r_go [m]

N_CPN_end                 # m x N_mc, terminal N_CPN per missile
mean_N_CPN                # m x 1, mean terminal N_CPN over MC
mean_r_min                # m x 1, mean of r_min_end over MC

N_CPN_ts                  # m x N_mc cell, N_CPN time-series
epsilon_ts                # m x N_mc cell, epsilon time-series [sec]
a_com_ts                  # m x N_mc cell, commanded acceleration [m/sec^2]
a_uncapped_ts             # m x N_mc cell, uncapped acceleration [m/sec^2]
v_gust_ts                 # 1 x N_mc cell, gust profile per scenario [m/sec]
t_run_vec                 # 1 x N_mc cell, actual time vector per scenario [sec]

t_end                     # 1 x N_mc, simulation end time per scenario [sec]
impact_time               # 1 x N_mc, impact time for all-hit scenarios, NaN otherwise [sec]

t_hit                     # m x N_mc, per-missile hit time per scenario, NaN if missile did not hit [sec]
hit_flags_missile         # m x N_mc logical, per-missile hit status per scenario
t_hit_first               # 1 x N_mc, first missile hit time per scenario [sec]
t_hit_last                # 1 x N_mc, last missile hit time per scenario [sec]
t_hit_spread              # 1 x N_mc, t_hit_last - t_hit_first per scenario [sec]

hit_angle                 # m x N_mc x 2, final [psi, theta] on all-hit scenarios [rad]

tau_m, tau_f, tau_sk      # 1 x N_mc, sampled lag constants
sigma_RM, sigma_VM        # 3 x N_mc, sampled measurement-noise std per axis

A_max                     # scalar, acceleration cap [g]
