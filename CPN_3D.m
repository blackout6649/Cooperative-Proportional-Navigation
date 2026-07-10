%% GUIDANCE AND NAVIGATION - FINAL PROJECT NUMERICAL SIMULATION

%% 3-D GNC simulation function
function data = CPN_3D(m, N, K_gain, A_max, R_hit, tau_m, tau_f, tau_sk, tend, dt, T_x0, M_x0_vec, M_V_vec, M_gamma0_vec, v_gust, gvis, save_fig, sigma_RM, sigma_VM)
if nargin < 19
    sigma_VM = 0;
end
if nargin < 18
    sigma_RM = 0;
end
if nargin == 16
    save_fig = 0;
elseif nargin == 15
    save_fig = 0;
    gvis     = 0;
end

%% preallcation
t_vec = 0:dt:tend;
n_sim = length(t_vec);
a_max = 9.8 * A_max;
K = K_gain;
hit = 0;

if isvector(v_gust)
    v_gust = repmat(v_gust(:), 1, n_sim);
elseif size(v_gust, 1) ~= 3 || size(v_gust, 2) ~= n_sim
    error('v_gust must be 3x1 or 3xn_sim.');
end

sigma_RM = expand_noise_std(sigma_RM);
sigma_VM = expand_noise_std(sigma_VM);

VMtot       = M_V_vec;
RT          = zeros(3, n_sim); RT(:, 1) = T_x0;
RM          = cell([m 1]);
gammaM      = cell([m 1]);

r_rel       = cell([m 1]);
VM          = cell([m 1]);
V_rel       = cell([m 1]);
r_go        = cell([m 1]);
r_go_dot    = cell([m 1]);
t_go        = cell([m 1]);
sigma       = cell([m 1]);
epsilon     = cell([m 1]);
omega       = cell([m 1]);
omega_m     = cell([m 1]);
omega_f     = cell([m 1]);
N_CPN       = cell([m 1]);
a_com       = cell([m 1]);
a_ideal     = cell([m 1]);
w_com       = cell([m 1]);
mean_t_go   = zeros(1, n_sim);
VT          = zeros(3, n_sim);


for i = 1:m
    gammaM{i}       = zeros(3, n_sim);
    gammaM{i}(:, 1) = deg2rad(M_gamma0_vec{i});
    
    RM{i}           = zeros(3, n_sim);
    RM{i}(:, 1)     = M_x0_vec{i};

    r_rel       {i} = zeros(3, n_sim);
    VM          {i} = zeros(3, n_sim);
    V_rel       {i} = zeros(3, n_sim);
    omega       {i} = zeros(3, n_sim);
    omega_m     {i} = zeros(3, n_sim);
    omega_f     {i} = zeros(3, n_sim);
    a_com       {i} = zeros(3, n_sim);
    a_ideal     {i} = zeros(3, n_sim);
    w_com       {i} = zeros(3, n_sim);
    r_go        {i} = zeros(1, n_sim);
    r_go_dot    {i} = zeros(1, n_sim);
    t_go        {i} = zeros(1, n_sim);
    N_CPN       {i} = zeros(1, n_sim);
    sigma       {i} = zeros(1, n_sim);
    epsilon     {i} = zeros(1, n_sim);
end

% w_max = a_max .* VMtot;
% VT = T_V0;
% aT = T_a;
% gammaT = zeros(1, n_sim); gammaT(1)    = deg2rad(T_gamma0);

%%  main simulation loop
for t = 1:n_sim

for i = 1:m
    phi     = gammaM{i}(1, t); % roll
    theta   = gammaM{i}(2, t); % pitch
    psi     = gammaM{i}(3, t); % yaw

    r_rel_true      = RT(:, t) - RM{i}(:, t);
    r_go_true       = sqrt(sum(r_rel_true .^ 2));

    RM_meas         = RM{i}(:, t) + sigma_RM .* randn(3, 1);
    r_rel{i}(:, t)  = RT(:, t) - RM_meas;
    r_go{i}(t)      = max(sqrt(sum(r_rel{i}(:, t) .^ 2)), eps);
if t == 1
    VM{i}(:, t)     = VMtot{i} .* [cos(theta)*cos(psi),...
                                   cos(theta)*sin(psi),...
                                   sin(theta)];
end

    VM_meas         = VM{i}(:, t) + sigma_VM .* randn(3, 1);

    V_rel{i}(:, t)  = VT(:, t) - (VM_meas + v_gust(:, t));
    r_go_dot{i}(t)  = dot(r_rel{i}(:, t), V_rel{i}(:, t)) / r_go{i}(t);
    omega{i}(:, t)  = cross(r_rel{i}(:, t), V_rel{i}(:, t)) / r_go{i}(t);

    if tau_sk > 0 && t < n_sim
        omega_m{i}(:, t+1)  = omega_m{i}(:, t) + (omega{i}(:, t) - omega_m{i}(:, t));
    elseif tau_sk == 0
        omega_m{i}(:, t)    = omega{i}(:, t);
    end

    if tau_f > 0 && t < n_sim
        omega_f{i}(:, t+1)  = omega_f{i}(:, t) + (omega_m{i}(:, t) - omega_f{i}(:, t));
    elseif tau_f == 0
        omega_f{i}(:, t)    = omega_m{i}(:, t);
    end

    vm_meas_norm = max(norm(VM_meas), eps);
    cos_sigma = dot(r_rel{i}(:, t), VM_meas) / (r_go{i}(t) * vm_meas_norm);
    sigma{i}(t) = acos(max(-1, min(1, cos_sigma)));

    t_go{i}(t)  = (1 + ( sigma{i}(t)^2 / (2*(2*N-1)) )) * (r_go{i}(:, t) / VMtot{i});
    
    mean_t_go(t) = mean_t_go(t) + t_go{i}(t) / (m-1);
end

if t == 1
    mean_t0 = 0;
    mean_r0 = 0;

    for i = 1:m
        mean_t0 = mean_t0 + t_go{i}(t) / m;
        mean_r0 = mean_r0 + r_go{i}(t) / m;
    end

    K = K_gain / (mean_r0 * mean_t0);
end

for i = 1:m
    epsilon{i}(t) = mean_t_go(t) - (m ./ (m-1)) .* t_go{i}(t);

    N_CPN{i}(t) = N * (1 - K * r_go{i}(t) * epsilon{i}(t));

    a = N_CPN{i}(t) .* cross(omega_f{i}(:, t), VM{i}(:, t));
    a_ideal{i}(:, t) = sign(a) .* min(abs(a), a_max);
    if tau_m > 0 && t < n_sim
        a_com{i}(:, t+1)  = a_com{i}(:, t) + (a_ideal{i}(:, t) - a_com{i}(:, t));
    elseif tau_m == 0
        a_com{i}(:, t)    = a_ideal{i}(:, t);
    end
    
    w_com{i}(:, t) = -cross(VM{i}(:, t), a_com{i}(:, t)) / (VMtot{i}^2);

    if t < n_sim
    RM{i}(:, t+1)   = RM{i}(:, t) + dt .* (VM{i}(:, t) + v_gust(:, t));
    RT(:, t+1)      = RT(:, t);
    
    gammaM{i}(:, t+1)  = gammaM{i}(:, t) + dt .* w_com{i}(:, t);
    V = (VM{i}(:, t) + a_com{i}(:, t) .* dt);
    VM{i}(:, t+1) = V .* VMtot{i}/sqrt(sum(V.^2));
    end

    if r_go_true <= R_hit
        hit = 1;
        break
    end
end
    if hit == 1; break; end
end

if hit == 1
    t_vec           = t_vec     (:, 1:t);
    RT              = RT        (:, 1:t);
    VT              = VT        (:, 1:t);
    mean_t_go       = mean_t_go (:, 1:t);

    for i = 1:m
    r_rel       {i} = r_rel     {i}(:, 1:t);
    RM          {i} = RM        {i}(:, 1:t);
    gammaM      {i} = gammaM    {i}(:, 1:t);
    VM          {i} = VM        {i}(:, 1:t);
    V_rel       {i} = V_rel     {i}(:, 1:t);
    omega       {i} = omega     {i}(:, 1:t);
    omega_m     {i} = omega_m   {i}(:, 1:t);
    omega_f     {i} = omega_f   {i}(:, 1:t);
    a_com       {i} = a_com     {i}(:, 1:t);
    a_ideal     {i} = a_ideal   {i}(:, 1:t);
    w_com       {i} = w_com     {i}(:, 1:t);
    r_go        {i} = r_go      {i}(:, 1:t);
    r_go_dot    {i} = r_go_dot  {i}(:, 1:t);
    t_go        {i} = t_go      {i}(:, 1:t);
    N_CPN       {i} = N_CPN     {i}(:, 1:t);
    sigma       {i} = sigma     {i}(:, 1:t);
    epsilon     {i} = epsilon   {i}(:, 1:t);
    end
end


%% graphing

xT = RT(1, :);
yT = RT(2, :);
zT = RT(3, :);

xM = zeros(m, length(t_vec));
yM = zeros(m, length(t_vec));
zM = zeros(m, length(t_vec));

for i = 1:m
    xM(i, :) = RM{i}(1, :);
    yM(i, :) = RM{i}(2, :);
    zM(i, :) = RM{i}(3, :);
end
marker_step = 5000;

xM_marker = [xM(:, 1:marker_step:end), xM(:, end)];
yM_marker = [yM(:, 1:marker_step:end), yM(:, end)];
zM_marker = [zM(:, 1:marker_step:end), zM(:, end)];
xT_marker = [xT(:, 1:marker_step:end), xT(:, end)];
yT_marker = [yT(:, 1:marker_step:end), yT(:, end)];
zT_marker = [zT(:, 1:marker_step:end), zT(:, end)];

min_x = min([min(xM), xT]);
min_y = min([min(yM), yT]);
min_z = min([min(zM), zT]);
max_x = max([max(xM), xT]);
max_y = max([max(yM), yT]);
max_z = max([max(zM), zT]);

% figure(Visible=gvis)
if gvis || save_fig
figure(Visible=gvis)

    plot3(xM(1, :), yM(1, :), zM(1, :), 'r-', LineWidth=1.5, DisplayName="Missile paths")
hold on
    
        plot3([xM_marker(1, :); xT_marker], [yM_marker(1, :); yT_marker], [zM_marker(1, :); zT_marker], ...
             'k--', HandleVisibility='off')
        
        plot3(xM_marker(1, :), yM_marker(1, :), zM_marker(1, :), 'or', LineWidth=2, ...
             MarkerSize=7, MarkerFaceColor='auto', HandleVisibility='off')
    
    if m > 1; for i = 2:m
    
        plot3(xM(i, :), yM(i, :), zM(i, :), 'r-', LineWidth=1.5, HandleVisibility='off')
        
            plot3([xM_marker(i, :); xT_marker], [yM_marker(i, :); yT_marker], [zM_marker(i, :); zT_marker], ...
                 'k--', HandleVisibility='off')
            
            plot3(xM_marker(i, :), yM_marker(i, :), zM_marker(i, :), 'or', LineWidth=2, ...
                 MarkerSize=7, MarkerFaceColor='auto', HandleVisibility='off')
    
    end; end
    
    plot3(xT, yT, zT, 'ob-.', MarkerIndices=[1, length(xT)], MarkerSize=7, ...
        MarkerFaceColor='auto', LineWidth=2.5, HandleVisibility='off')
    
        plot3(xT_marker, yT_marker, zT_marker, 'ob', LineWidth=2, ...
            MarkerSize=7, MarkerFaceColor='auto', DisplayName="Target")
    
    title("3D CPN Simulation", m + " Cooprative missiles")
    
    legend(autoupdate="on", Location="southeast")
    xlabel("x [m]")
    ylabel("z [m]")
    xlim([min_x - 100, max_x + 100])
    ylim([min_y - 100, max_y + 100])
    zlim([min_z - 100, max_z + 100])
    grid on
    hold off

if save_fig
    set(gcf,'units','pix','pos',[0,0,1920,1080])
    saveas(gcf, "3D_GNC_sim_" + m + "_missiles" + ".png")
end
end

%% data logging
data = struct;

data.t_vec       = t_vec     ;
data.RT          = RT        ;
% data.VT          = VT        ;
% data.mean_t_go   = mean_t_go ;
data.r_rel       = r_rel     ;
data.RM          = RM        ;
data.gammaM      = gammaM    ;
data.VM          = VM        ;
data.V_rel       = V_rel     ;
data.omega       = omega     ;
data.omega_m     = omega_m   ;
data.omega_f     = omega_f   ;
data.a_com       = a_com     ;
data.a_ideal     = a_ideal   ;
% data.w_com       = w_com     ;
data.r_go        = r_go      ;
data.r_go_dot    = r_go_dot  ;
data.t_go        = t_go      ;
data.N_CPN       = N_CPN     ;
data.sigma       = sigma     ;
data.epsilon     = epsilon   ;
data.hit         = hit       ;
data.R_hit       = R_hit     ;
data.K           = K         ;
data.m           = m         ;
data.N           = N         ;
data.A_max       = A_max     ;
data.VMtot       = VMtot     ;
data.v_gust      = v_gust    ;


end

% inputs:
    % m             - number of cooprative missiles (m > 1)
    % N             - parallel navigation constant
    % K_gain        - gain for cooprative navigation
    % A_max         - max acceleration allowed for missiles, in [G]
    % tau_m         - mechanical lag time constant, in [sec]
    % tau_f         - filter lag time constant, in [sec]
    % tau_sk        - seeker lag time constant, in [sec]
    % R_hit         - hit detection radius, in [m]
    % tend          - end time if no hit detected, in [sec]
    % dt            - simulation time step, in [sec] 
    % T_x0          - target initial position [x0, y0, z0]; in [m]

    % M_x0_vec      - missile initial positions; mx1 cell; in [m]
                        % M_x0_vec{i} = [x0, y0, z0] for missile i

    % M_V_vec       - missile velocities; mx1 cell; in [m/sec]
                        % M_V_vec{i} = V_tot for missile i
                        
    % M_gamma0_vec  - missile initial headings; mx1 cell; in [deg]
                        % M_gamma0_vec{i} = [phi0, theta0, psi0] for missile i

    % v_gust        - wind gust velocity profile [m/sec]
                        % either 3x1 (constant gust) or 3xn_sim (time-varying)

    % gvis          - show graph of missile paths, 1 = show, 0 = no show
    % save_fig      - save graph of missile paths, 1 = save, 0 = no save
    % sigma_RM      - RM measurement noise std, scalar or 3x1, in [m]
    % sigma_VM      - VM measurement noise std, scalar or 3x1, in [m/sec]


% output:
    % data - struct containing the output values

    % data.t_vec       - simulation time vector, in [sec]
    % data.RT          - target location over time, in [m]
    % data.r_rel       - target location relative to missile over time, in [m]
    % data.RM          - missiles' location over time, in [m]
    % data.gammaM      - missiles' heading angles over time [phi; theta; psi], in [rad]
    % data.VM          - missiles' velocity over time, in [m/sec]
    % data.V_rel       - target velocity relative to missile over time, in [m]
    % data.omega       - ideal LOS rate of change over time, in [rad/sec]
    % data.omega_m     - measured LOS rate of change over time, in [rad/sec]
    % data.omega_f     - measured & filtered LOS rate of change over time, in [rad/sec]
    % data.a_ideal     - ideal acceleration command for each missile, in [m/sec^2]
    % data.a_com       - true acceleration command for each missile, in [m/sec^2]
    % data.r_go        - missiles' closing distance over time, in [m]
    % data.r_go_dot    - missiles' closing velocity over time, in [m/sec]
    % data.t_go        - estimated time to go for each missile, in [sec]
    % data.N_CPN       - time-varying navigation gain
    % data.sigma       - spatial heading error, in [rad]
    % data.epsilon     - relative ttg error for the cooprative control
    % data.hit         - simulation hit flag, 1 = hit, 0 = no hit
    % data.R_hit       - hit detection radius, in [m]
    % data.K           - cooprative navigation constant
    % data.m           - number of missiles
    % data.N           - parallel navigation constant
    % data.A_max       - max acceleration allowed for missiles, in [G]
    % data.VMtot       - missile total velocities, in [m/sec]
    % data.v_gust      - gust profile used in the simulation, in [m/sec]

function sigma = expand_noise_std(sigma)
    if isscalar(sigma)
        sigma = sigma .* ones(3, 1);
    else
        sigma = sigma(:);
    end

    if numel(sigma) ~= 3
        error('Measurement noise std must be scalar or 3x1.');
    end
end