% =====================================================================
% 4-VSCMG (Variable Speed Control Moment Gyroscope) – Main Simulation
% =====================================================================
addpath("Functions")
addpath("Graphics")
clear; clc;
warning('off')

%% ---------------------------------------------------------
% Spacecraft & VSCMG inertias
%% ---------------------------------------------------------
I_spacecraft   = diag([86, 85, 113]);      % Spacecraft inertia [kg·m^2]
J_vscmg        = diag([0.13, 0.04, 0.03]); % VSCMG frame inertia about (s,t,g) axes
J_vscmg_table  = repmat([0.13, 0.04, 0.03], 4, 1); % (unused table form; kept for completeness)
I_wheel        = 0.1;                      % Wheel (rotor) inertia about s-axis [kg·m^2]
J_cmg_effect   = J_vscmg(1,1) - I_wheel;   % Effective CMG inertia about spin axis

% Wheel spin speeds (one per VSCMG), row vector for convenience
wheelOmega = [14.4 14.4 14.4 14.4];        % [rad/s]

% Initial gimbal angles and rates
gimbalAngle0     = [0 0 90 -90] * pi/180;  % [rad]
gimbalRate0      = [0 0 0 0];              % [rad/s]

% Pyramid half-angle of the VSCMG layout (typical for 4-CMG pyramids)
pyramidAngle = 54.75 * pi/180;             % [rad]

%% ---------------------------------------------------------
% Build initial gimbal frames (spin s, transverse t, gimbal g) in body {B}
% For each VSCMG k, we define unit vectors g_s(k), g_t(k), g_g(k) in {B}.
%% ---------------------------------------------------------
% VSCMG 1
g_g1 = [ cos(pyramidAngle); 0;  sin(pyramidAngle)];
g_s1 = [ 0; 1; 0];
g_t1 = cross(g_g1, g_s1);

% VSCMG 2
g_g2 = [-cos(pyramidAngle); 0;  sin(pyramidAngle)];
g_s2 = [ 0; -1; 0];
g_t2 = cross(g_g2, g_s2);

% VSCMG 3
g_g3 = [ 0;  cos(pyramidAngle); sin(pyramidAngle)];
g_s3 = [ 1; 0; 0];
g_t3 = cross(g_g3, g_s3);

% VSCMG 4
g_g4 = [ 0; -cos(pyramidAngle); sin(pyramidAngle)];
g_s4 = [-1; 0; 0];
g_t4 = cross(g_g4, g_s4);

% Stack the 4 triads (each column is one VSCMG)
Gs0 = [g_s1, g_s2, g_s3, g_s4];            % spin axes in {B}
Gt0 = [g_t1, g_t2, g_t3, g_t4];            % transverse axes in {B}
Gg0 = [g_g1, g_g2, g_g3, g_g4];            % gimbal axes in {B}

% Build the 12x3 "stacked DCM" from gimbal frames to body:
%   For k=1..4, rows (3k-2:3k, :) = [g_s(k), g_t(k), g_g(k)]
GS = [Gs0(:,1); Gs0(:,2); Gs0(:,3); Gs0(:,4)];
GT = [Gt0(:,1); Gt0(:,2); Gt0(:,3); Gt0(:,4)];
GG = [Gg0(:,1); Gg0(:,2); Gg0(:,3); Gg0(:,4)];
DCM_bg0 = [GS, GT, GG];                    % 12x3 stacked orientation for all VSCMGs

%% ---------------------------------------------------------
% Initial spacecraft attitude & rates
%% ---------------------------------------------------------
sigmaBN   = [0.414; 0.3; 0.2];             % Attitude (MRPs) of body {B} w.r.t inertial {N}
omegaBN_B = [0.01; 0.05; -0.01];           % Body angular velocity {B} [rad/s]

L_ext = 0; u_s = 0; u_g = 0; u_body = 0;   
omegaBN_B_dot = [0; 0; 0];                 % Angular acceleration init

%% ---------------------------------------------------------
% Initialize  inertia
%% ---------------------------------------------------------
DCM_bg  = DCM_bg0;                          % Current stacked gimbal->body DCMs
DCM_bg1 = DCM_bg0( 1: 3, 1:3);
DCM_bg2 = DCM_bg0( 4: 6, 1:3);
DCM_bg3 = DCM_bg0( 7: 9, 1:3);
DCM_bg4 = DCM_bg0(10:12, 1:3);

gimbalAngle     = gimbalAngle0;             % Current gimbal angles
gimbalRate      = gimbalRate0;              % Current gimbal rates

% (Optional) attitude DCM from MRPs (not strictly needed below)
DCM_bn = ConvertAttitude(sigmaBN, 'MRP', 'DCM');

% Effective inertia including VSCMG frames (about body)
I_body = I_spacecraft ...
    + DCM_bg1*J_vscmg*DCM_bg1' ...
    + DCM_bg2*J_vscmg*DCM_bg2' ...
    + DCM_bg3*J_vscmg*DCM_bg3' ...
    + DCM_bg4*J_vscmg*DCM_bg4';

%% ---------------------------------------------------------
% Simulation settings and loggers
%% ---------------------------------------------------------
% Steering command memory (for numerical gimbal acceleration)
gimbalRateCmd      = [0 0 0 0];             % γ̇_cmd(k) from previous step
gimbalAccelCmd     = [0 0 0 0];             % γ̈_cmd(k) (approx via finite diff)

% Logs

%% ---------------------------------------------------------
% Preallocate logs
%% ---------------------------------------------------------
t_final = 200;                  % Simulation duration [s]
dt      = 0.005;                  % Time step [s] (define properly here)
time    = 0:dt:t_final;
Nsteps  = length(time);

% Sizes (adjust depending on your states)
log_sigmaBN = zeros(3, Nsteps);        % Attitude (MRPs)
log_wheel   = zeros(4, Nsteps);        % Wheel speeds [rad/s]
log_gimbal  = zeros(4, Nsteps);        % Gimbal angles [rad]
log_frames  = cell(1, Nsteps);         % DCM per step

wheelangles = zeros(Nsteps,4);         % Integrated wheel angles
wheelangle  = [0 0 0 0];
gimbalTmp   = [0 0 0 0];
gimbalAngles= zeros(4,Nsteps);

w_s = zeros(1,4);
w_t = zeros(1,4);
w_g = zeros(1,4);
wheelOmegaDot= zeros(1,4);

%% ---------------------------------------------------------
% Main integration loop
%% ---------------------------------------------------------
for k = 1:Nsteps
    t = time(k);

    % Keep MRPs in principal region
    sigmaBN = MRPswitch(sigmaBN);

    % Current attitude DCM
    DCM_bn = ConvertAttitude(sigmaBN,'MRP','DCM'); 

    % --------- Logging ---------
    log_sigmaBN(:,k) = sigmaBN;
    log_wheel(:,k)   = wheelOmega;
    log_gimbal(:,k)  = gimbalAngle;
    wheelangles(k,:) = wheelangle;
    log_frames{k}    = DCM_bg;

    % Express body rate in each gimbal frame
    for j = 1:4
        [w_s(j), w_t(j), w_g(j)] = omega_body2gimbal_frame(j, DCM_bg, omegaBN_B); 
    end

    % Body torque command from attitude controller (PD)
    u_body = control(I_body, t, dt, omegaBN_B, sigmaBN, I_wheel, DCM_bg, ...
                     w_s, w_t, w_g, wheelOmega);

    % Steering call for gimbal accel prediction
    [~, gimbalRateCmdPred] = Stering(DCM_bg, J_vscmg, wheelOmega, w_s, u_body);
    gimbalAccelCmd = (gimbalRateCmdPred - gimbalRateCmd) / dt;

    % Final steering call for wheel & gimbal commands
    [wheelOmegaDot_cmd, gimbalRateCmd] = Stering(DCM_bg, J_vscmg, wheelOmega, w_s, u_body);

    % Actuator models
    for j = 1:4
        wheelOmegaDot(j) = Wheel_Motor_Eqn(j, wheelOmegaDot_cmd);
        gimbalRate(j)    = Gimbal_Motor_Eqn(j, gimbalRate, gimbalRateCmd, gimbalAccelCmd);
    end

    % Rigid body dynamics with VSCMGs
    omegaBN_B_dot = omega_dot_eqn( ...
        I_body, I_wheel, J_vscmg, DCM_bg, ...
        w_s, w_t, w_g, ...
        wheelOmega, wheelOmegaDot, ...
        gimbalRate, gimbalAccelCmd, ...
        omegaBN_B);

    % Attitude kinematics
    sigmaBN_dot = dMRP(sigmaBN, omegaBN_B);

    % ---------------- Time update (Euler forward) ----------------
    gimbalAngle = gimbalAngle + gimbalRate * dt;
    wheelOmega  = wheelOmega  + wheelOmegaDot * dt;
    wheelangle  = wheelangle  + wheelOmega * dt;
    sigmaBN     = sigmaBN     + sigmaBN_dot * dt;
    omegaBN_B   = omegaBN_B   + omegaBN_B_dot * dt;

    % Update gimbal frame DCMs
    DCM_bg1 = gimbal_frame_vectors(1, gimbalAngle0, DCM_bg0, gimbalAngle);
    DCM_bg2 = gimbal_frame_vectors(2, gimbalAngle0, DCM_bg0, gimbalAngle);
    DCM_bg3 = gimbal_frame_vectors(3, gimbalAngle0, DCM_bg0, gimbalAngle);
    DCM_bg4 = gimbal_frame_vectors(4, gimbalAngle0, DCM_bg0, gimbalAngle);
    DCM_bg  = [DCM_bg1; DCM_bg2; DCM_bg3; DCM_bg4];

    % Refresh inertia with new gimbal orientations
    I_body = I_spacecraft ...
        + DCM_bg1*J_vscmg*DCM_bg1' ...
        + DCM_bg2*J_vscmg*DCM_bg2' ...
        + DCM_bg3*J_vscmg*DCM_bg3' ...
        + DCM_bg4*J_vscmg*DCM_bg4';
end


%% ---------------------------------------------------------
% Plots
%% ---------------------------------------------------------
subplot(1,3,1)
plot(time, log_sigmaBN)
title('Attitude (MRP)')
xlabel('Time [s]'); ylabel('\sigma')
grid on

subplot(1,3,2)
plot(time, log_gimbal * 180/pi)
title('Gimbal Angles')
xlabel('Time [s]'); ylabel('Angle [deg]')
grid on

subplot(1,3,3)
plot(time, log_wheel)
title('Wheel Speeds')
xlabel('Time [s]'); ylabel('\Omega [rad/s]')
grid on

function BG = gimbal_frame_vectors(i, gamma_0, BG00, gamma)
% gimbal_frame_vectors
% ---------------------
% Computes the rotated gimbal frame vectors for the i-th gimbal.
%
% Inputs:
%   i       - index of the gimbal (1, 2, ...)
%   gamma_0 - vector of initial gimbal angles [rad]
%   BG00    - stacked initial gimbal frames (each is 3x3, stacked vertically)
%   gamma   - vector of current gimbal angles [rad]
%
% Output:
%   BG      - [3x3] matrix of gimbal frame vectors after rotation
%             (columns are g_s, g_t, g_g)

    % Extract initial frame for i-th gimbal
    BG0 = BG00((i-1)*3+1:(i-1)*3+3, 1:3);

    % Initial basis vectors
    gs0 = BG0(:,1);
    gt0 = BG0(:,2);
    gg0 = BG0(:,3);

    % Rotate about gimbal axis (gg0)
    dgamma = gamma(i) - gamma_0(i);
    g_s = cos(dgamma) * gs0 + sin(dgamma) * gt0;
    g_t = -sin(dgamma) * gs0 + cos(dgamma) * gt0;
    g_g = gg0;

    % Updated gimbal frame
    BG = [g_s, g_t, g_g];
end


function [ws, wt, wg] = omega_body2gimbal_frame(i, DCM_bg_all, omega_BN_B)
% omega_body2gimbal_frame
% -----------------------------------------
% Map the body angular velocity vector into the i-th VSCMG gimbal frame.
%
% Inputs
%   i            : VSCMG index (1..4)
%   DCM_bg_all   : 12x3 block-stacked DCMs from each gimbal frame to body
%                  (each 3x3 block is [g_s, g_t, g_g] in body axes)
%   omega_BN_B   : 3x1 body angular velocity {B} [rad/s]
%
% Outputs
%   ws, wt, wg   : components of body rate expressed in the gimbal frame:
%                  spin-axis (s), transverse-axis (t), and gimbal-axis (g)

% Extract the i-th 3x3 DCM (gimbal->body)
DCM_bg_i = DCM_bg_all((i-1)*3+1:(i-1)*3+3, 1:3);

% To express a body vector in the gimbal frame, use body->gimbal = DCM_bg_i'
omega_BN_G = DCM_bg_i' * omega_BN_B;

ws = omega_BN_G(1);
wt = omega_BN_G(2);
wg = omega_BN_G(3);
end
function omegaDot = Wheel_Motor_Eqn(i, OMEGA_dot_cmd)
% Wheel_Motor_Eqn
% -----------------------------------------
% (Simple wheel actuator model / placeholder)
% Maps the commanded wheel angular acceleration to the actual one.
%
% Inputs
%   i             : wheel index (1..4)  [not used here, kept for interface]
%   OMEGA_dot_cmd : 1x4 commanded wheel angular accelerations [rad/s^2]
%
% Output
%   omegaDot      : scalar angular acceleration for wheel i [rad/s^2]

omegaDot = OMEGA_dot_cmd(i);   % pass-through (no motor dynamics modeled)
end
function Gammadoubledot = Gimbal_Motor_Eqn(i, ~, gamma_dot_cmd, ~)
% Gimbal_Motor_Eqn
% -----------------------------------------
% (Simple gimbal actuator model / placeholder)
% Returns the commanded gimbal angular acceleration directly.
%
% Inputs
%   i              : gimbal index (1..4)  [not used here, kept for interface]
%   gamma_dot      : 1x4 current gimbal rates [rad/s] (unused in this simple model)
%   gamma_dot_cmd  : 1x4 commanded gimbal rates [rad/s]
%   gamma_ddot_cmd : 1x4 commanded gimbal accelerations [rad/s^2] (unused)
%
% Output
%   Gammadoubledot : scalar gimbal angular acceleration for gimbal i [rad/s^2]

Gammadoubledot = gamma_dot_cmd(i);   % pass-through (no motor dynamics modeled)
end

function omega_dot = omega_dot_eqn(I_B, I_wheel, J_vscmg, DCM_bg, w_s, w_t, w_g, ...
                                   OMEGA, OMEGA_DOT, gamma_dot, gamma_ddot, omega_BN_B)
% omega_dot_eqn
% -----------------------------------------
% Rigid-body rotational dynamics with 4 VSCMGs.
% Computes body angular acceleration given current state and actuator rates.
%
% Inputs
%   I_B         : 3x3 total spacecraft inertia in body frame (incl. VSCMG frames)
%   I_wheel     : scalar wheel inertia about spin axis
%   J_vscmg     : 3x3 inertia of the VSCMG gimbal frame
%   DCM_bg      : 12x3 stacked DCMs [g_s, g_t, g_g] for 4 VSCMGs (gimbal->body)
%   w_s, w_t, w_g : 1x4 components of body rate expressed in each gimbal frame
%   OMEGA       : 1x4 wheel spin rates [rad/s]
%   OMEGA_DOT   : 1x4 wheel spin accelerations [rad/s^2]
%   gamma_dot   : 1x4 gimbal rates [rad/s]
%   gamma_ddot  : 1x4 gimbal accelerations [rad/s^2]
%   omega_BN_B  : 3x1 body angular velocity [rad/s]
%
% Output
%   omega_dot   : 3x1 body angular acceleration [rad/s^2]

% Build direction matrices for all gimbals:
% Each column of GS, Gt, Gg is a 3x1 axis (in {B}) for one VSCMG.
GS = [DCM_bg(1:3,1),   DCM_bg(4:6,1),   DCM_bg(7:9,1),   DCM_bg(10:12,1)];
Gt = [DCM_bg(1:3,2),   DCM_bg(4:6,2),   DCM_bg(7:9,2),   DCM_bg(10:12,2)];
Gg = [DCM_bg(1:3,3),   DCM_bg(4:6,3),   DCM_bg(7:9,3),   DCM_bg(10:12,3)];

% Internal momentum-rate contributions in each gimbal frame:
%   \dot{H}_s, \dot{H}_t, \dot{H}_g for each VSCMG (stacked as 4x1 vectors).
% Derivation follows standard VSCMG momentum relations.
g_s_hat_part = [ ...
    J_vscmg(1,1)*gamma_dot(1)*w_t(1) + I_wheel*OMEGA_DOT(1) - (J_vscmg(2,2)-J_vscmg(3,3))*w_t(1)*gamma_dot(1); ...
    J_vscmg(1,1)*gamma_dot(2)*w_t(2) + I_wheel*OMEGA_DOT(2) - (J_vscmg(2,2)-J_vscmg(3,3))*w_t(2)*gamma_dot(2); ...
    J_vscmg(1,1)*gamma_dot(3)*w_t(3) + I_wheel*OMEGA_DOT(3) - (J_vscmg(2,2)-J_vscmg(3,3))*w_t(3)*gamma_dot(3); ...
    J_vscmg(1,1)*gamma_dot(4)*w_t(4) + I_wheel*OMEGA_DOT(4) - (J_vscmg(2,2)-J_vscmg(3,3))*w_t(4)*gamma_dot(4) ...
];

g_t_hat_part = [ ...
    (J_vscmg(1,1)*w_s(1) + I_wheel*OMEGA(1))*gamma_dot(1) - (J_vscmg(2,2)+J_vscmg(3,3))*w_s(1)*gamma_dot(1) + I_wheel*OMEGA(1)*w_g(1); ...
    (J_vscmg(1,1)*w_s(2) + I_wheel*OMEGA(2))*gamma_dot(2) - (J_vscmg(2,2)+J_vscmg(3,3))*w_s(2)*gamma_dot(2) + I_wheel*OMEGA(2)*w_g(2); ...
    (J_vscmg(1,1)*w_s(3) + I_wheel*OMEGA(3))*gamma_dot(3) - (J_vscmg(2,2)+J_vscmg(3,3))*w_s(3)*gamma_dot(3) + I_wheel*OMEGA(3)*w_g(3); ...
    (J_vscmg(1,1)*w_s(4) + I_wheel*OMEGA(4))*gamma_dot(4) - (J_vscmg(2,2)+J_vscmg(3,3))*w_s(4)*gamma_dot(4) + I_wheel*OMEGA(4)*w_g(4) ...
];

g_g_hat_part = [ ...
    J_vscmg(3,3)*gamma_ddot(1) - I_wheel*OMEGA(1)*w_t(1); ...
    J_vscmg(3,3)*gamma_ddot(2) - I_wheel*OMEGA(2)*w_t(2); ...
    J_vscmg(3,3)*gamma_ddot(3) - I_wheel*OMEGA(3)*w_t(3); ...
    J_vscmg(3,3)*gamma_ddot(4) - I_wheel*OMEGA(4)*w_t(4) ...
];

% Net internal torque in {B} from all VSCMGs:
Hdot_internal_B = - GS * g_s_hat_part - Gt * g_t_hat_part - Gg * g_g_hat_part;

% Euler’s rotational equation: I*omega_dot = -omega×(I*omega) + Hdot_internal
H_dot_total_B = -cross(omega_BN_B, I_B*omega_BN_B) + Hdot_internal_B;

omega_dot = I_B \ H_dot_total_B;   % numerically better than inv(I_B)*H_dot_total_B
end
% function y = w_sig_dot(x, u, L)
% % w_sig_dot
% % -----------------------------------------
% % Combined attitude kinematics (MRP) and rigid-body dynamics.
% % State x = [sigma(1:3); omega(1:3)]
% %
% % Inputs
% %   x : 6x1 state [sigma; omega]
% %   u : 3x1 control torque in body frame
% %   L : 3x1 external/disturbance torque in body frame
% %
% % Output
% %   y : 6x1 time derivative [sigma_dot; omega_dot]
% 
% sigma = x(1:3);
% omega = x(4:6);
% 
% % MRP kinematics: sigma_dot = 0.25 * B(sigma) * omega
% B_sigma = (1 - norm(sigma)^2) * eye(3) + 2*skew(sigma) + 2*(sigma*sigma');
% sigma_dot = 0.25 * B_sigma * omega;
% 
% % Simple rigid-body dynamics with diagonal inertia (example values)
% I_body = diag([100; 75; 80]);
% omega_dot = I_body \ (-skew(omega)*I_body*omega + u + L);
% 
% y = [sigma_dot; omega_dot];
% end
function u = control(~, ~, ~, omega_BN_B, sigma_BN, I_wheel, DCM_bg, ~, w_t, w_g, OMEGA)
% control
% -----------------------------------------
% Simple PD attitude controller in body frame.
%
% Inputs
%   I_body      : 3x3 spacecraft inertia
%   t, dt       : time and step (unused here; kept for extensibility)
%   omega_BN_B  : 3x1 body angular velocity
%   sigma_BN    : 3x1 MRP attitude error/state
%   I_wheel     : wheel inertia (unused here; reserved)
%   DCM_bg      : 12x3 gimbal->body DCM stack (used to form Gt, Gg)
%   w_s,w_t,w_g : 1x4 body-rate components in each gimbal frame
%   OMEGA       : 1x4 wheel speeds
%
% Output
%   u           : 3x1 control torque in body frame

% Axes collections
Gt = [DCM_bg(1:3,2),   DCM_bg(4:6,2),   DCM_bg(7:9,2),   DCM_bg(10:12,2)];
Gg = [DCM_bg(1:3,3),   DCM_bg(4:6,3),   DCM_bg(7:9,3),   DCM_bg(10:12,3)];

% (Optional) internal momentum coupling term (not used in torque here)
% sumL = Σ I_wheel*OMEGA_i*(w_g_i * Gt_i - w_t_i * Gg_i)
% Kept to show how wheel speed couples into apparent torques.
sumL = I_wheel * ( ...
      OMEGA(1) * ( w_g(1) * Gt(:,1) - w_t(1) * Gg(:,1) ) ...
    + OMEGA(2) * ( w_g(2) * Gt(:,2) - w_t(2) * Gg(:,2) ) ...
    + OMEGA(3) * ( w_g(3) * Gt(:,3) - w_t(3) * Gg(:,3) ) ...
    + OMEGA(4) * ( w_g(4) * Gt(:,4) - w_t(4) * Gg(:,4) ) );

% PD gains (tune as needed)
Kp = 1.7;                                        % attitude (MRP) gain
Kd = diag([13.13, 13.04, 15.08]);              % rate feedback gains

% Control torque (body frame)
u = -Kp * sigma_BN - Kd * omega_BN_B+sumL;

% If you want to include momentum coupling explicitly:
% u = -Kp*sigma_BN - Kd*omega_BN_B + sumL;
end
function [dOmegad, dgamad] = Stering(DCM_bg, J_vscmg, OMEGA, w_s, u)
% Stering (intentional spelling preserved for compatibility)
% -----------------------------------------
% Steering law that maps desired body torque u to actuator accelerations:
%   Q * [dOMEGA; dGAMMA] = -u   (weighted least-squares)
%
% Inputs
%   DCM_bg : 12x3 gimbal->body DCM stack for the 4 VSCMGs
%   J_vscmg: 3x3 VSCMG frame inertia
%   OMEGA  : 1x4 wheel spin rates [rad/s]
%   w_s    : 1x4 body-rate s-axis components in each gimbal frame
%   u      : 3x1 desired body torque
%
% Outputs
%   dOmegad : 1x4 wheel spin accelerations [rad/s^2]
%   dgamad  : 1x4 gimbal angular accelerations [rad/s^2]

% Axis collections per gimbal
GS = [DCM_bg(1:3,1),   DCM_bg(4:6,1),   DCM_bg(7:9,1),   DCM_bg(10:12,1)];
Gt = [DCM_bg(1:3,2),   DCM_bg(4:6,2),   DCM_bg(7:9,2),   DCM_bg(10:12,2)];

% Build the allocation matrix Q = [D0  D1] such that:
%   body torque ≈ -Q * [dOMEGA; dGAMMA]
% where
%   D0 columns (per wheel) ~ effect of wheel accel about s-axis
%   D1 columns (per gimbal) ~ effect of gimbal accel via transverse axis
D0 = [ GS(:,1)*J_vscmg(1,1), GS(:,2)*J_vscmg(1,1), GS(:,3)*J_vscmg(1,1), GS(:,4)*J_vscmg(1,1) ];
D1 = [ Gt(:,1)*J_vscmg(1,1)*(OMEGA(1)+w_s(1)), ...
       Gt(:,2)*J_vscmg(1,1)*(OMEGA(2)+w_s(2)), ...
       Gt(:,3)*J_vscmg(1,1)*(OMEGA(3)+w_s(3)), ...
       Gt(:,4)*J_vscmg(1,1)*(OMEGA(4)+w_s(4)) ];

Q = [D0, D1];   % 3x8

% Simple scalar weight that downweights near-singular configurations
de   = det(D1*D1');           % proxy for gimbal geometry richness
Ws   = 2*exp(-1e-9 * de);     % small when det is large (OK), near 2 when small
W    = diag([Ws, Ws, Ws, Ws, 1, 1, 1, 1]);   % 8x8

% Weighted least-squares solution: minimize || W^(1/2) * (Q*x + u) ||
% x = [dOMEGA; dGAMMA]
% NOTE: pinv(Q*W*Q') is numerically safer than inv(...)
x = - W * Q' * ((Q*W*Q') \ u);

dOmegad = x(1:4)';     % 1x4
dgamad  = x(5:8)';     % 1x4
end

