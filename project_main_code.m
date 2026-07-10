%% GUIDANCE AND NAVIGATION - FINAL PROJECT NUMERICAL SIMULATION

%% 2-D TEST SIMULATION - random initial missile conditions
gvis_2D         = 1;
save_fig_2D     = 0;

m_2D            = 3;
N_2D            = 5;
K_gain_2D       = 20;
A_max_2D        = 30;
tau_2D          = 0;
R_hit_2D        = 0.05;

tend_2D         = 30;
dt_2D           = 0.001;

T_x0_2D         = [0, 2500]; % stationary target

M_x0_vec_2D     = zeros(m_2D, 2);
M_V_vec_2D      = zeros(m_2D, 1);
M_gamma0_vec_2D = zeros(m_2D, 1);

% random starting conditions
    for i = 1:m_2D
        M_x0_vec_2D(i, :)   = -500 + 1000 .* [rand, rand];
        M_V_vec_2D(i)       = 100 + 50 * rand;
        M_gamma0_vec_2D(i)  = 60 + 60 * rand;
    end

data_2D_test = CPN_2D(m_2D, N_2D, K_gain_2D, A_max_2D, R_hit_2D, tau_2D,...
    tend_2D, dt_2D, T_x0_2D, M_x0_vec_2D, M_V_vec_2D, M_gamma0_vec_2D, ...
    gvis_2D, save_fig_2D);

%% 3-D TEST SIMULATION - random initial missile conditions
gvis_3D         = 1;
save_fig_3D     = 0;

m_3D            = 3;
N_3D            = 3;
K_gain_3D       = 40;
A_max_3D        = 100;
R_hit_3D        = 1;
tau_m_3D        = 0;
tau_f_3D        = 0;
tau_sk_3D       = 0;
v_gust          = [100; 50; 20];

tend_3D         = 30;
dt_3D           = 0.001;

T_x0_3D         = [0, 0, 2500];

M_x0_vec_3D     = cell([m_3D 1]);
M_V_vec_3D      = cell([m_3D 1]);
M_gamma0_vec_3D = cell([m_3D 1]);

% random starting conditions

    for i = 1:m_3D
        M_x0_vec_3D{i}     = -500 + 1000 .* [rand, rand, rand];
        M_V_vec_3D{i}      = 100 + 50 * rand;
        M_gamma0_vec_3D{i} = 60 + 60 * [rand, rand, rand];
    end

data_3D_test = CPN_3D(m_3D, N_3D, K_gain_3D, A_max_3D, R_hit_3D, tau_m_3D, ...
    tau_f_3D, tau_sk_3D, tend_3D, dt_3D, T_x0_3D, M_x0_vec_3D, M_V_vec_3D, ...
    M_gamma0_vec_3D, v_gust, gvis_3D, save_fig_3D);