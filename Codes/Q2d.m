% Parameters for membrane capacitance
mem_cap_std = 1.0; % standard membrane capacitance, in uF/cm^2
mem_cap_incr = 2.0; % increased membrane capacitance, in uF/cm^2

max_gK = 36.0; % maximum conductance for potassium, in mS/cm^2
max_gNa = 120.0; % maximum conductance for sodium, in mS/cm^2
leak_gL = 0.3; % leakage conductance, in mS/cm^2
rev_pot_K = -77.0; % reversal potential for potassium, in mV
rev_pot_Na = 50.0; % reversal potential for sodium, in mV
rev_pot_L = -54.4; % leakage reversal potential, in mV

% Time parameters
time_span = [0 100]; % time range, in ms
max_time_step = 0.01; % desired maximum time step, in ms

% Alpha and Beta functions
alpha_m_func = @(V) 0.1 * (V + 40) ./ (1 - exp(-(V + 40) / 10));
beta_m_func = @(V) 4 * exp(-(V + 65) / 18);
alpha_h_func = @(V) 0.07 * exp(-(V + 65) / 20);
beta_h_func = @(V) 1 ./ (1 + exp(-(V + 35) / 10));
alpha_n_func = @(V) 0.01 * (V + 55) ./ (1 - exp(-(V + 55) / 10));
beta_n_func = @(V) 0.125 * exp(-(V + 65) / 80);

% Initial conditions
init_V = -65.0;
init_m = alpha_m_func(init_V) / (alpha_m_func(init_V) + beta_m_func(init_V));
init_h = alpha_h_func(init_V) / (alpha_h_func(init_V) + beta_h_func(init_V));
init_n = alpha_n_func(init_V) / (alpha_n_func(init_V) + beta_n_func(init_V));
init_conditions = [init_V, init_n, init_m, init_h];

% Applied current function for the given amplitude and duration
current_amp1 = 32.70; % in uA/cm^2
current_duration = 0.2; % in ms

applied_current_func = @(t) ((t >= 10 & t <= (10 + current_duration)) * current_amp1);

% Define solver options to limit the timestep
solver_options = odeset('MaxStep', max_time_step);

% ODE solver with standard membrane capacitance
[time_std, state_std] = ode45(@(t, y) HodgkinHuxleyModel(t, y, mem_cap_std, max_gK, max_gNa, leak_gL, rev_pot_K, rev_pot_Na, rev_pot_L, applied_current_func(t), alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func), time_span, init_conditions, solver_options);

% ODE solver with increased membrane capacitance
[time_incr, state_incr] = ode45(@(t, y) HodgkinHuxleyModel(t, y, mem_cap_incr, max_gK, max_gNa, leak_gL, rev_pot_K, rev_pot_Na, rev_pot_L, applied_current_func(t), alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func), time_span, init_conditions, solver_options);

% Extract membrane potentials
mem_pot_std = state_std(:, 1);
mem_pot_incr = state_incr(:, 1);

% Plot membrane potentials for standard and increased capacitance
figure;
plot(time_std, mem_pot_std, 'b', time_incr, mem_pot_incr, 'r');
yyaxis right
plot(time_std, arrayfun(applied_current_func, time_std));
ylabel('Applied Current (\muA/cm^2)');
title(sprintf('Action Potential with Standard (%.2f \muF/cm^2) and Increased (%.2f \muF/cm^2) Membrane Capacitance', mem_cap_std, mem_cap_incr));
legend(['C_m = ', num2str(mem_cap_std), ' \muF/cm^2'], ['C_m = ', num2str(mem_cap_incr), ' \muF/cm^2'],'Applied Current');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
grid on;

% Hodgkin-Huxley model function
function dy = HodgkinHuxleyModel(~, y, mem_cap, max_gK, max_gNa, leak_gL, rev_pot_K, rev_pot_Na, rev_pot_L, applied_current, alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func)
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);

    % Ionic currents
    I_K = max_gK * n^4 * (V - rev_pot_K);
    I_Na = max_gNa * m^3 * h * (V - rev_pot_Na);
    I_L = leak_gL * (V - rev_pot_L);

    % Differential equations
    dVdt = (applied_current - I_K - I_Na - I_L) / mem_cap;
    dndt = alpha_n_func(V) * (1 - n) - beta_n_func(V) * n;
    dmdt = alpha_m_func(V) * (1 - m) - beta_m_func(V) * m;
    dhdt = alpha_h_func(V) * (1 - h) - beta_h_func(V) * h;

    dy = [dVdt; dndt; dmdt; dhdt];
end
