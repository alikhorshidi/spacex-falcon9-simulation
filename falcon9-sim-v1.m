% =========================================================================
% Project: SpaceX Falcon 9 Ascent Simulation (3D)
% Course: Fundamentals of Programming - Aerospace Engineering
% Instructor/Author: Dr. Ali Khorshidi Benam
%
% VIDEO TUTORIAL (Step-by-Step): https://www.youtube.com/watch?v=0PvsDQLXEQI
%
% Description: Numerical integration of rocket motion using Euler's Method.
% 
% DISCLAIMER: 
% This code is for EDUCATIONAL PURPOSES ONLY. The parameters, physics constants, 
% and simulation results are simplified approximations and may differ from 
% actual real-world telemetry data of the Falcon 9 vehicle.
%
% COMPATIBILITY:
% - Developed and Tested on: MATLAB 2025b
% - Minimum Requirement: MATLAB R2016b or later
% =========================================================================

clc; clear; close all;

%% 1. Simulation Setup
dt = 0.1;                   % Time Step (0.1 seconds per frame)
T_total = 300;              % Total Simulation Time (300 seconds)
time = 0:dt:T_total;        % Time Vector
N = length(time);           % Total number of calculation steps

% Pre-allocation (Prevents the code from running slowly)
h = zeros(1, N);            % Altitude (meters)
x = zeros(1, N);            % Downrange Distance (meters)
vx = zeros(1, N);           % Horizontal Velocity (m/s)
vz = zeros(1, N);           % Vertical Velocity (m/s)
m = zeros(1, N);            % Instantaneous Mass (kg)
pitch = zeros(1, N);        % Pitch Angle (degrees relative to horizon)

%% 2. Initial Conditions (at t=0)
h(1) = 0;                   % Starting on the ground
x(1) = 0;                   
vx(1) = 0; vz(1) = 0;       % Initial velocity is zero
m(1) = 549054;              % Falcon 9 Lift-off Mass (approx. kg)
pitch(1) = 90;              % Rocket starts vertically (90 degrees)

g0 = 9.81;                  % Standard Gravity (m/s^2)

%% 3. Main Physics Loop
disp('Simulation Started...');

for i = 1:N-1
    t_current = time(i);
    
    % --- A) Environmental Sensor (Atmosphere) ---
    % Function defined at the bottom of the script
    [rho, ~] = get_atmosphere(h(i));
    
    % --- B) Propulsion Status (Engine Controller) ---
    % Is the engine on? What is the thrust? Is staging occurring?
    [Thrust, m_dot, stage_mode] = get_engine_status(t_current, h(i));
    
    % --- C) Staging Event ---
    % If the computer commands staging (Mode 2), drop the dry mass of Stage 1
    if stage_mode == 2 && time(i-1) < 162
         m(i) = m(i) - 25600; % Subtract Stage 1 dry mass (approx. 25.6 tons)
    end
    
    % Update Mass: Current Mass - (Burn Rate * Time Step)
    m(i+1) = m(i) - (m_dot * dt);
    
    % --- D) Guidance & Autopilot (Gravity Turn) ---
    % Strategy: Go straight up initially, then tilt gradually after 1km altitude
    if h(i) > 1000 && pitch(i) > 20
        pitch(i+1) = pitch(i) - (0.04 * dt); % Slow pitch-over rate
    else
        pitch(i+1) = pitch(i);
    end
    
    % Convert angle to radians (since sin/cos require radians)
    theta_rad = deg2rad(pitch(i));
    
    % --- E) Aerodynamic Forces (Physics Engine) ---
    V_total = sqrt(vx(i)^2 + vz(i)^2); % Total Velocity Magnitude
    
    % Drag Calculation: D = 0.5 * rho * V^2 * Cd * A
    Area = 10.75; % Falcon 9 Cross-sectional Area (m^2)
    Cd = 0.3;     % Drag Coefficient
    Drag = 0.5 * rho * V_total^2 * Cd * Area;
    
    % Force Decomposition (Drag always opposes motion)
    if V_total > 0
        Dx = Drag * (vx(i) / V_total);
        Dz = Drag * (vz(i) / V_total);
    else
        Dx=0; Dz=0;
    end
    
    % Thrust Decomposition (Thrust is aligned with pitch angle)
    Tx = Thrust * cos(theta_rad);
    Tz = Thrust * sin(theta_rad);
    
    % Newton's Second Law: a = F_net / m
    ax = (Tx - Dx) / m(i);              % Horizontal Acceleration
    Weight = m(i) * g0;
    az = (Tz - Weight - Dz) / m(i);     % Vertical Acceleration
    
    % --- F) Numerical Integration (Euler's Method) ---
    % Next Velocity = Current Velocity + Acceleration * Time Step
    vx(i+1) = vx(i) + ax * dt;
    vz(i+1) = vz(i) + az * dt;
    
    % Next Position = Current Position + Velocity * Time Step
    x(i+1) = x(i) + vx(i+1) * dt;
    h(i+1) = h(i) + vz(i+1) * dt;
    
    % Ground Collision Check
    if h(i+1) < 0
        h(i+1) = 0; vz(i+1)=0; vx(i+1)=0;
    end
end

%% 4. Visualization
figure('Name', 'SpaceX Falcon 9 Telemetry', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);

% --- Plot 1: 3D Flight Profile ---
subplot(2, 2, [1 3]);
surface([x/1000; x/1000], [time; time], [h/1000; h/1000], ...
    [sqrt(vx.^2+vz.^2); sqrt(vx.^2+vz.^2)], ...
    'FaceColor', 'no', 'EdgeColor', 'interp', 'LineWidth', 2.5);
hold on;
% Draw Earth Surface
fill3([-10 100 100 -10], [0 0 300 300], [0 0 0 0], [0.8 0.9 0.8], 'EdgeAlpha', 0); 
grid on; view(50, 20); colorbar;
c = colorbar; c.Label.String = 'Velocity (m/s)';
xlabel('Downrange (km)'); ylabel('Time (s)'); zlabel('Altitude (km)');
title('Falcon 9 Ascent Profile (3D View)');

% --- Plot 2: Altitude & Velocity ---
subplot(2, 2, 2);
yyaxis left; plot(time, h/1000, 'b-', 'LineWidth', 2); ylabel('Altitude (km)');
yyaxis right; plot(time, sqrt(vx.^2 + vz.^2)/343, 'r-', 'LineWidth', 2); ylabel('Mach Number');
title('Flight Performance'); grid on; xline(162, '--k', 'MECO');

% --- Plot 3: Dynamic Pressure (Max-Q) ---
subplot(2, 2, 4);
rho_vec = 1.225 * exp(-h/8500);
Q = 0.5 .* rho_vec .* (vx.^2 + vz.^2);
plot(time, Q/1000, 'k-', 'LineWidth', 1.5);
title('Dynamic Pressure (Q)'); ylabel('Q (kPa)'); xlabel('Time (s)'); grid on;
[maxQ, idx] = max(Q); text(time(idx), maxQ/1000, ' \leftarrow Max-Q', 'Color', 'r', 'FontWeight', 'bold');

% =========================================================================
% Helper Functions
% =========================================================================

function [rho, P] = get_atmosphere(altitude)
    % Calculate Air Density based on Altitude (Isothermal Exponential Model)
    rho0 = 1.225;       % Sea Level Density (kg/m^3)
    H_scale = 8500;     % Scale Height (m)
    
    rho = rho0 * exp(-altitude / H_scale);
    P = 101325 * exp(-altitude / H_scale); % Pressure (Pa)
end

function [F, m_dot, mode] = get_engine_status(t, h)
    g0 = 9.81;
    MECO_TIME = 162;    % Main Engine Cut-Off Time (seconds)
    SEP_DELAY = 4;      % Separation Delay (seconds)
    
    if t <= MECO_TIME
        % --- Phase 1: Main Boost (Stage 1) ---
        mode = 1;
        F_sl = 7607000;   % Thrust at Sea Level (N)
        F_vac = 8227000;  % Thrust in Vacuum (N)
        
        % Thrust increases with altitude (up to 50km approximation)
        factor = min(h / 50000, 1); 
        F = F_sl + (F_vac - F_sl) * factor;
        
        % Isp (Specific Impulse) also increases
        Isp = 282 + (311 - 282) * factor; 
        
    elseif t > MECO_TIME && t <= MECO_TIME + SEP_DELAY
        % --- Phase 2: Coasting & Separation ---
        mode = 2; F = 0; Isp = 1; % Isp=1 avoids division by zero
        
    else
        % --- Phase 3: Stage 2 Ignition (Vacuum Engine) ---
        mode = 3;
        F = 934000; % Lower thrust but highly efficient
        Isp = 348;  
    end
    
    % Calculate Mass Flow Rate (Fuel Consumption)
    if F == 0
        m_dot = 0; 
    else
        m_dot = F / (Isp * g0); 
    end
end
