%% GUIDANCE AND NAVIGATION - FINAL PROJECT NUMERICAL SIMULATION

%% 2-D GNC simulation function
function data = CPN_2D(m, N, K_gain, A_max, R_hit, tau, tend, dt, T_x0, M_x0_vec, M_V_vec, M_gamma0_vec, gvis, save_fig)
if nargin == 13
    save_fig = 0;
elseif nargin == 12
    save_fig = 0;
    gvis     = 0;
end
%% preallcation
t_vec = 0:dt:tend;
n_sim = length(t_vec);
a_max = 9.8 * A_max;
K = K_gain;
hit = 0;

xM = zeros(m, n_sim); xM(:, 1) = M_x0_vec(:, 1);
zM = zeros(m, n_sim); zM(:, 1) = M_x0_vec(:, 2);

VM = M_V_vec;

w_max = a_max .* VM;

xT = zeros(1, n_sim); xT(1) = T_x0(1);
zT = zeros(1, n_sim); zT(1) = T_x0(2);

% VT = T_V0;
% aT = T_a;

% gammaT = zeros(1, n_sim); gammaT(1)    = deg2rad(T_gamma0);
gammaM = zeros(m, n_sim); gammaM(:, 1) = deg2rad(M_gamma0_vec);

lambda      = zeros(m, n_sim);
sigma       = zeros(m, n_sim);
VMx         = zeros(m, n_sim);
VMz         = zeros(m, n_sim);
% VTx         = zeros(1, n_sim);
% VTz         = zeros(1, n_sim);
r_go        = zeros(m, n_sim);
r_go_dot    = zeros(m, n_sim);
lambda_dot  = zeros(m, n_sim);
t_go        = zeros(m, n_sim);
epsilon     = zeros(m, n_sim);
omega       = zeros(m, n_sim);
N_CPN       = zeros(m, n_sim);
a_com       = zeros(m, n_sim);
w_com       = zeros(m, n_sim);
w_ideal     = zeros(m, n_sim);
w           = zeros(m, n_sim);
min_r       = zeros(1, n_sim);



%%  main simulation loop
for i = 1:n_sim
t = t_vec(i);

VMx(:, i)   = VM .* cos(gammaM(:, i));
VMz(:, i)   = VM .* sin(gammaM(:, i));
% VTx(i)      = VT .* cos(gammaT(i));
% VTz(i)      = VT .* sin(gammaT(i));

lambda(:, i)        = atan2(zT(i)-zM(:, i), xT(i)-xM(:, i));
sigma(:, i)         = gammaM(:, i) - lambda(:, i);

r_go(:, i)          = sqrt( (zT(i)-zM(:, i)).^2 + (xT(i)-xM(:, i)).^2 );
r_go_dot(:, i)      = -1 .* VM .* cos(sigma(:, i));

lambda_dot(:, i)    = -1 .* VM .* sin(sigma(:, i)) ./ r_go(:, i);
t_go(:, i)          = (r_go(:, i) ./ VM) .* (1 + (sigma(:, i)).^2 ./ (2.*(2.*N-1)));

if i == 1
    K = K_gain / (mean(r_go(:, i)) * mean(t_go(:, i)));
end

epsilon(:, i)       = (1 ./ (m-1)) .* sum(t_go(:, i)) - (m ./ (m-1)) .* t_go(:, i);

omega(:, i)         = K .* r_go(:, i) .* epsilon(:, i);

N_CPN(:, i)         = N .* (1 - omega(:, i));

w(:, i)             = N_CPN(:, i) .* lambda_dot(:, i);


w_ideal(:, i) = sign(w(:, i)) .* min(abs(w(:, i)), w_max);

if tau > 0 && i > 1
    w_com(:, i) = w_com(:, i-1) + dt/tau * (w_ideal(:, i) - w_com(:, i-1));
elseif i > 1
    w_com(:, i) = w_ideal(:, i);
end

a_com(:, i) = VM .* w_com(:, i);

if i < n_sim

xM(:, i+1)      = xM(:, i) + dt .* VMx(:, i);
zM(:, i+1)      = zM(:, i) + dt .* VMz(:, i);
xT(i+1)         = xT(i); %    + dt .* VTx(i)   ;
zT(i+1)         = zT(i); %    + dt .* VTz(i)   ;

gammaM(:, i+1)  = gammaM(:, i) + dt .* w_com(:, i);

end

min_r(i) = min(r_go(:, i));

if min_r(i) <= R_hit
    hit = 1;

    t_vec       = t_vec     (:, 1:i);
    xM          = xM        (:, 1:i);
    zM          = zM        (:, 1:i);
    xT          = xT        (:, 1:i);
    zT          = zT        (:, 1:i);
    gammaM      = gammaM    (:, 1:i);
    lambda      = lambda    (:, 1:i);
    sigma       = sigma     (:, 1:i);
    VMx         = VMx       (:, 1:i);
    VMz         = VMz       (:, 1:i);
    r_go        = r_go      (:, 1:i);
    r_go_dot    = r_go_dot  (:, 1:i);
    lambda_dot  = lambda_dot(:, 1:i);
    t_go        = t_go      (:, 1:i);
    epsilon     = epsilon   (:, 1:i);
    omega       = omega     (:, 1:i);
    N_CPN       = N_CPN     (:, 1:i);
    a_com       = a_com     (:, 1:i);
    w_com       = w_com     (:, 1:i);
    w_ideal     = w_ideal   (:, 1:i);
    w           = w         (:, 1:i);
    min_r       = min_r     (:, 1:i);
    
    break
end

end

%% graphing
marker_step = 5000;

xM_marker = [xM(:, 1:marker_step:end), xM(:, end)];
zM_marker = [zM(:, 1:marker_step:end), zM(:, end)];
xT_marker = [xT(:, 1:marker_step:end), xT(:, end)];
zT_marker = [zT(:, 1:marker_step:end), zT(:, end)];

min_x = min([min(xM), xT]);
min_z = min([min(zM), zT]);
max_x = max([max(xM), xT]);
max_z = max([max(zM), zT]);


figure(Visible=gvis)
hold on

    plot(xM(1, :), zM(1, :), 'r-', LineWidth=1.5, DisplayName="Missile paths")
    
    plot([xM_marker(1, :); xT_marker], [zM_marker(1, :); zT_marker], ...
         'k--', HandleVisibility='off')
    
    plot(xM_marker(1, :), zM_marker(1, :), 'or', LineWidth=2, ...
         MarkerSize=7, MarkerFaceColor='auto', HandleVisibility='off')
    
    if m > 1; for i = 2:m
    
        plot(xM(i, :), zM(i, :), 'r-', LineWidth=1.5, HandleVisibility='off')
        
        plot([xM_marker(i, :); xT_marker], [zM_marker(i, :); zT_marker], ...
             'k--', HandleVisibility='off')
        
        plot(xM_marker(i, :), zM_marker(i, :), 'or', LineWidth=2, ...
             MarkerSize=7, MarkerFaceColor='auto', HandleVisibility='off')
    
    end; end
    
    plot(xT, zT, 'b-.', LineWidth=2.5, HandleVisibility='off')
    
    plot(xT_marker, zT_marker, 'ob', LineWidth=2, ...
         MarkerSize=7, MarkerFaceColor='auto', DisplayName="Target")
    
    title("2D CPN Simulation", m + " Cooprative missiles")
    
    legend(autoupdate="on", Location="southeast")
    xlabel("x [m]")
    ylabel("z [m]")
    xlim([min_x - 100, max_x + 100])
    ylim([min_z - 100, max_z + 100])
    grid on
    hold off

if save_fig
    set(gcf,'units','pix','pos',[0,0,1920,1080])
    saveas(gcf, "2D_GNC_sim_" + m + "_missiles" + ".png")
end

%% data logging
data = struct;

data.t_vec       = t_vec     ;
data.xM          = xM        ;
data.zM          = zM        ;
data.xT          = xT        ;
data.zT          = zT        ;
data.gammaM      = gammaM    ;
data.lambda      = lambda    ;
data.sigma       = sigma     ;
data.VMx         = VMx       ;
data.VMz         = VMz       ;
data.r_go        = r_go      ;
data.r_go_dot    = r_go_dot  ;
data.lambda_dot  = lambda_dot;
data.t_go        = t_go      ;
data.epsilon     = epsilon   ;
data.omega       = omega     ;
data.N_CPN       = N_CPN     ;
data.a_com       = a_com     ;
data.w_com       = w_com     ;
% data.w_ideal     = w_ideal   ;
% data.w           = w         ;
% data.min_r       = min_r     ;
data.hit         = hit       ;
data.R_hit       = R_hit     ;
data.K           = K         ;
data.min_miss    = min(min_r);
data.m           = m         ;
data.N           = N         ;
data.A_max       = A_max     ;
data.VM          = VM        ;


end

% inputs:
    % m             - number of cooprative missiles (m > 1)
    % N             - parallel navigation constant
    % K_gain        - gain for cooprative navigation
    % A_max         - max acceleration allowed for missiles, in [G]
    % tau           - time constant for first order dynamics, in [sec]
    % R_hit         - hit detection radius, in [m]
    % tend          - end time if no hit detected, in [sec]
    % dt            - simulation time step, in [sec] 
    % T_x0          - target initial position [x0, z0]; in [m]

    % M_x0_vec      - missile initial positions; mx2 array; in [m]
                        % M_x0_vec(i, :) = [x0, z0] for missile i

    % M_V_vec       - missile velocities; mx1 array; in [m/sec]
                        % M_V_vec(i, :) = V_tot for missile i
                        
    % M_gamma0_vec  - missile initial headings; mx1 array; in [deg]
                        % M_gamma0_vec(i, :) = heading angle for missile i

    % gvis          - show graph of missile paths, 1 = show, 0 = no show
    % save_fig      - save graph of missile paths, 1 = save, 0 = no save


% output:
    % data - struct containing the output values

    % data.t_vec       - simulation time vector, in [sec]
    % data.xM          - missiles' x values over time, in [m]
    % data.zM          - missiles' z values over time, in [m]
    % data.xT          - target x values over time, in [m]
    % data.zT          - target z values over time, in [m]
    % data.gammaM      - missiles' heading angles over time, in [rad]
    % data.lambda      - LOS angles over time, in [rad]
    % data.sigma       - missiles' heading errors over time, in [rad]
    % data.VMx         - missiles' velocity in x over time, in [m/sec]
    % data.VMz         - missiles' velocity in z over time, in [m/sec]
    % data.r_go        - missiles' closing distance over time, in [m]
    % data.r_go_dot    - missiles' closing velocity over time, in [m/sec]
    % data.lambda_dot  - LOS rate of change over time, in [rad/sec]
    % data.t_go        - estimated time to go for each missile, in [sec]
    % data.epsilon     - relative ttg error for the cooprative control
    % data.omega       - cooprative feedback term
    % data.N_CPN       - time-varying navigation gain
    % data.a_com       - acceleration command for each missile, in [m/sec^2]
    % data.w_com       - turn rate command for each missile, in [rad/sec]
    % data.hit         - simulation hit flag, 1 = hit, 0 = no hit
    % data.R_hit       - hit detection radius, in [m]
    % data.K           - cooprative navigation constant
    % data.min_miss    - closest any missile got to the target, in [m]
    % data.m           - number of missiles
    % data.N           - parallel navigation constant
    % data.A_max       - max acceleration allowed for missiles, in [G]
    % data.VM          - missile total velocities, in [m/sec]


