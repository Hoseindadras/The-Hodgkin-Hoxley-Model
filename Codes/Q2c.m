% Parameters
mem_capacitance = 1.0; % membrane capacitance, in uF/cm^2
max_gK = 36.0; % maximum conductance for potassium, in mS/cm^2
max_gNa = 120.0; % maximum conductance for sodium, in mS/cm^2
leakage_gL = 0.3; % leakage conductance, in mS/cm^2
reversal_pot_K = -77.0; % reversal potential for potassium, in mV
reversal_pot_Na = 50.0; % reversal potential for sodium, in mV
reversal_pot_L = -54.4; % leakage reversal potential, in mV

% Time parameters
time_range = [0 100]; % time range, in ms
max_step_size = 0.01; % desired maximum time step, in ms

% Alpha and Beta functions
alpha_m_func = @(V) 0.1 * (V + 40) ./ (1 - exp(-(V + 40) / 10));
beta_m_func = @(V) 4 * exp(-(V + 65) / 18);
alpha_h_func = @(V) 0.07 * exp(-(V + 65) / 20);
beta_h_func = @(V) 1 ./ (1 + exp(-(V + 35) / 10));
alpha_n_func = @(V) 0.01 * (V + 55) ./ (1 - exp(-(V + 55) / 10));
beta_n_func = @(V) 0.125 * exp(-(V + 65) / 80);

% Initial conditions
initial_V = -65.0;
initial_m = alpha_m_func(initial_V) / (alpha_m_func(initial_V) + beta_m_func(initial_V));
initial_h = alpha_h_func(initial_V) / (alpha_h_func(initial_V) + beta_h_func(initial_V));
initial_n = alpha_n_func(initial_V) / (alpha_n_func(initial_V) + beta_n_func(initial_V));
initial_conditions = [initial_V, initial_n, initial_m, initial_h];

% Applied current function for the given amplitude and duration
current_amplitude = 34; % in uA/cm^2
current_duration = 0.2; % in ms
applied_current_func = @(t) (t >= 10 & t <= (10 + current_duration)) * current_amplitude;

% Define solver options to limit the timestep
solver_options = odeset('MaxStep', max_step_size);

% ODE solver with applied current
[time, state_vars] = ode45(@(t, y) HodgkinHuxleyModel(t, y, mem_capacitance, max_gK, max_gNa, leakage_gL, reversal_pot_K, reversal_pot_Na, reversal_pot_L, applied_current_func(t), alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func), time_range, initial_conditions, solver_options);

% Extract variables
membrane_potential = state_vars(:, 1);
gating_var_n = state_vars(:, 2);
gating_var_m = state_vars(:, 3);
gating_var_h = state_vars(:, 4);

% Calculate ionic currents
potassium_current = max_gK * gating_var_n.^4 .* (membrane_potential - reversal_pot_K);
sodium_current = max_gNa * gating_var_m.^3 .* gating_var_h .* (membrane_potential - reversal_pot_Na);
leak_current = leakage_gL .* ( membrane_potential - reversal_pot_L);

% Plot ionic currents IK and INa
figure;
plot(time, sodium_current, 'r', time, potassium_current, 'b', time, leak_current, 'g');
legend('I_{Na}', 'I_{K}', 'I_{L}');
title('Sodium and Potassium Currents Over Time');
xlabel('Time (ms)');
ylabel('Current (\muA/cm^2)');
grid on;

% Hodgkin-Huxley model
function dy = HodgkinHuxleyModel(~, y, mem_capacitance, max_gK, max_gNa, leakage_gL, reversal_pot_K, reversal_pot_Na, reversal_pot_L, applied_current, alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func)
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);

    % Ionic currents
    IK = max_gK * n^4 * (V - reversal_pot_K);
    INa = max_gNa * m^3 * h * (V - reversal_pot_Na);
    IL = leakage_gL * (V - reversal_pot_L);

    % Differential equations
    dVdt = (applied_current - IK - INa - IL) / mem_capacitance;
    dndt = alpha_n_func(V) * (1 - n) - beta_n_func(V) * n;
    dmdt = alpha_m_func(V) * (1 - m) - beta_m_func(V) * m;
    dhdt = alpha_h_func(V) * (1 - h) - beta_h_func(V) * h;

    dy = [dVdt; dndt; dmdt; dhdt];
end
