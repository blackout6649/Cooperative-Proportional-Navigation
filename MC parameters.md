## INPUT

tau_m = U(0.005, 0.05)            # mechanical lag time constant [sec]
tau_f = U(0.01, 0.1)              # filter lag time constant [sec]
tau_sk = U(0.01, 0.2)             # seeker lag time constant [sec]

M_x0_vec_3D{i} = [r0*cos(ang0), r0*sin(ang0), 0]    # missile initial position [m]
ang0 = 2*pi*(i - 1 + U(0,1))                         # angular sample (per missile)
r0 = sqrt(R_min_3D^2 + U(0,1)*(R_init_3D^2 - R_min_3D^2))
R_min_3D = 9000
R_init_3D = 15000

M_V_vec_3D{i} = U(200, 300)       # missile speed [m/sec]

M_gamma0_vec_3D{i} = [phi0, theta0, psi0]            # initial heading [deg]
phi0 = U(60, 120)
theta0 = U(80, 100)               # elevation around 90 deg (+/-10 deg)
psi0 = U(60, 120)

## ADDITIONAL

sigma_gust = U^3([2, 2, 1], [12, 12, 5])            # gust std per axis [m/sec]
L_gust = U^3([100, 100, 50], [400, 400, 200])       # Gauss-Markov length scale per axis [m]

v_gust(t) = first-order Gauss-Markov process         # 3 x n_sim gust profile [m/sec]

sigma_RM = U^3([0, 0, 0], [1, 1, 1])                # RM measurement-noise std [m]
sigma_VM = U^3([0, 0, 0], [0.01, 0.01, 0.01])       # VM measurement-noise std [m/sec]

## SUCCESS CRITERION

hit = 1 only if all missiles hit in the run
hit = 0 otherwise

## OUTPUT

N_CPN(t)                 # time-varying navigation gain; mx1 cell of time-series
epsilon(t)               # cooperative time-to-go error; mx1 cell of time-series
a_com(t)                 # commanded acceleration; mx1 cell of [3 x n_sim] [m/sec^2]
a_uncapped(t)            # uncapped guidance acceleration; mx1 cell of [3 x n_sim] [m/sec^2]
v_gust(t)                # gust time-series used in the run; 3 x n_sim

miss_distance_min        # closest-approach distance per run [m]
hit                      # run success flag (all missiles must hit)
P_hit                    # probability of successful runs over MC

t_impact                 # impact time for successful runs [sec]
hit_angle                # final hit attitude [psi, theta, phi] from gammaM(:, end) [rad]
