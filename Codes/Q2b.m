% Parameters
mem_cap = 1.0; % membrane capacitance, in uF/cm^2
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
initial_conditions = [init_V, init_n, init_m, init_h];

% Applied current function for the given amplitude and duration
current_amp = 32.70; % in uA/cm^2
current_dur = 0.2; % in ms
applied_current_func = @(t) (t >= 10 & t <= (10 + current_dur)) * current_amp;

% Define solver options to limit the timestep
solver_options = odeset('MaxStep', max_time_step);

% ODE solver with applied current
[time, states] = ode45(@(t, y) HodgkinHuxleyModel(t, y, mem_cap, max_gK, max_gNa, leak_gL, rev_pot_K, rev_pot_Na, rev_pot_L, applied_current_func(t), alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func), time_span, initial_conditions, solver_options);

% Extract gating variables and calculate conductances
V = states(:, 1);
n = states(:, 2);
m = states(:, 3);
h = states(:, 4);
conductance_K = max_gK * n.^4;
conductance_Na = max_gNa * m.^3 .* h;

% Calculate derivatives of gating variables
dn_dt = alpha_n_func(V) .* (1 - n) - beta_n_func(V) .* n;
dm_dt = alpha_m_func(V) .* (1 - m) - beta_m_func(V) .* m;
dh_dt = alpha_h_func(V) .* (1 - h) - beta_h_func(V) .* h;

% Plot membrane potential and input current
figure;
yyaxis left
plot(time, arrayfun(applied_current_func, time));
ylabel('Applied Current (\muA/cm^2)');
yyaxis right
plot(time, V);
title('Membrane Potential with 32.70 \muA/cm^2 Current for 0.2 ms');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
grid on;

% Plot conductances gNa and gK
figure;
yyaxis right
plot(time, arrayfun(applied_current_func, time));
ylabel('Applied Current (\muA/cm^2)');
yyaxis left
plot(time, conductance_Na, 'r', time, conductance_K, 'b');
legend('g_{Na}', 'g_{K}');
title('Sodium and Potassium Conductances Over Time');
xlabel('Time (ms)');
ylabel('Conductance (mS/cm^2)');
grid on;

% Plot gating variables m, n, and h
figure;
yyaxis right
plot(time, arrayfun(applied_current_func, time));
ylabel('Applied Current (\muA/cm^2)');
yyaxis left
plot(time, m, 'g', time, n, 'b', time, h, 'r');
legend('m', 'n', 'h');
title('Gating Variables Over Time');
xlabel('Time (ms)');
ylabel('Gating Variables');
grid on;

% Plot derivatives of gating variables and input current
figure;
yyaxis left
plot(time, dn_dt, 'b', time, dm_dt, 'g', time, dh_dt, 'r');
legend('dn/dt', 'dm/dt', 'dh/dt');
ylabel('Gating Variable Derivatives');
yyaxis right
plot(time, arrayfun(applied_current_func, time), 'k');
ylabel('Applied Current (\muA/cm^2)');
title('Gating Variable Derivatives and Applied Current Over Time');
xlabel('Time (ms)');
grid on;

% Hodgkin-Huxley model
function dy = HodgkinHuxleyModel(~, y, mem_cap, cond_K, cond_Na, leak_gL, rev_pot_K, rev_pot_Na, rev_pot_L, I_ext, alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func)
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);

    % Ionic currents
    I_K = cond_K * n^4 * (V - rev_pot_K);
    I_Na = cond_Na * m^3 * h * (V - rev_pot_Na);
    I_L = leak_gL * (V - rev_pot_L);

    % Differential equations
    dVdt = (I_ext - I_K - I_Na - I_L) / mem_cap;
    dndt = alpha_n_func(V) * (1 - n) - beta_n_func(V) * n;
    dmdt = alpha_m_func(V) * (1 - m) - beta_m_func(V) * m;
    dhdt = alpha_h_func(V) * (1 - h) - beta_h_func(V) * h;

    dy = [dVdt; dndt; dmdt; dhdt];
end
