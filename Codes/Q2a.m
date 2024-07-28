% Parameters
mem_cap = 1.0; % membrane capacitance, in uF/cm^2
cond_K = 36.0; % maximum conductance for potassium, in mS/cm^2
cond_Na = 120.0; % maximum conductance for sodium, in mS/cm^2
cond_L = 0.3; % leakage conductance, in mS/cm^2
rev_pot_K = -77.0; % reversal potential for potassium, in mV
rev_pot_Na = 50.0; % reversal potential for sodium, in mV
rev_pot_L = -54.4; % leakage reversal potential, in mV

% Time parameters
time_range = [0 100]; % time range, in ms
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
current_amp = 20; % in uA/cm^2
current_dur = 0.2; % in ms
applied_current_func = @(t) (t >= 10 & t <= (10 + current_dur)) * current_amp;

% Define solver options to limit the timestep
solver_options = odeset('MaxStep', max_time_step);

% ODE solver with applied current
[time, states] = ode45(@(t, y) HHModel(t, y, mem_cap, cond_K, cond_Na, cond_L, rev_pot_K, rev_pot_Na, rev_pot_L, applied_current_func(t), alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func), time_range, initial_conditions, solver_options);

% Plotting the results for the given amplitude
figure;
yyaxis left
plot(time, arrayfun(applied_current_func, time));
ylabel('Input Current (\muA/cm^2)');
yyaxis right
plot(time, states(:, 1));
title('Action Potential with 20 \muA/cm^2 Current for 0.2 ms');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
grid on;

% Finding the minimum amplitude of excitation
min_amp = 0;
spike_detected = false;
current_dur = 0.2; % in ms

while ~spike_detected
    min_amp = min_amp + 0.1; % Increment the amplitude
    applied_current_func = @(t) (t >= 10 & t <= (10 + current_dur)) * min_amp;

    % ODE solver with varying current
    [time, states] = ode45(@(t, y) HHModel(t, y, mem_cap, cond_K, cond_Na, cond_L, rev_pot_K, rev_pot_Na, rev_pot_L, applied_current_func(t), alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func), time_range, initial_conditions, solver_options);
    
    % Check for spike (threshold crossing)
    if any(states(:, 1) > 0) % If the membrane potential crosses 0 mV
        spike_detected = true;
    end
end

fprintf('The minimum amplitude of excitation for a spike is %.2f µA/cm^2\n', min_amp);

% Plotting the first spike for the minimum amplitude
figure;
yyaxis left
plot(time, arrayfun(applied_current_func, time));
ylabel('Input Current (\muA/cm^2)');
yyaxis right
plot(time, states(:, 1));
title(sprintf('First Spike with Minimum Amplitude %.2f µA/cm^2 for 0.2 ms', min_amp));
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
grid on;

% Hodgkin-Huxley model
function dy = HHModel(~, y, mem_cap, cond_K, cond_Na, cond_L, rev_pot_K, rev_pot_Na, rev_pot_L, I_ext, alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func)
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);

    % Ionic currents
    I_K = cond_K * n^4 * (V - rev_pot_K);
    I_Na = cond_Na * m^3 * h * (V - rev_pot_Na);
    I_L = cond_L * (V - rev_pot_L);

    % Differential equations
    dVdt = (I_ext - I_K - I_Na - I_L) / mem_cap;
    dndt = alpha_n_func(V) * (1 - n) - beta_n_func(V) * n;
    dmdt = alpha_m_func(V) * (1 - m) - beta_m_func(V) * m;
    dhdt = alpha_h_func(V) * (1 - h) - beta_h_func(V) * h;

    dy = [dVdt; dndt; dmdt; dhdt];
end
