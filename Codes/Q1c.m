% Parameters
mem_cap = 1.0; % membrane capacitance, in uF/cm^2
max_gK = 36.0; % maximum conductances, in mS/cm^2
max_gNa = 120.0; % maximum conductances, in mS/cm^2
leak_gL = 0.3; % leakage conductance, in mS/cm^2
rev_EK = -77.0; % reversal potentials, in mV
rev_ENa = 50.0; % reversal potentials, in mV
rev_EL = -54.4; % reversal potentials, in mV

% Time parameters
time_span = [0 100]; % time range, in ms
time_step = 0.001; % desired maximum time step, in ms

% Alpha and Beta functions
alpha_m_func = @(V) 0.1 * (V + 40) ./ (1 - exp(-(V + 40) / 10));
beta_m_func = @(V) 4 * exp(-(V + 65) / 18);
alpha_h_func = @(V) 0.07 * exp(-(V + 65) / 20);
beta_h_func = @(V) 1 ./ (1 + exp(-(V + 35) / 10));
alpha_n_func = @(V) 0.01 * (V + 55) ./ (1 - exp(-(V + 55) / 10));
beta_n_func = @(V) 0.125 * exp(-(V + 65) / 80);

% Define applied current function
applied_current_func = @(t, I_amp, I_dur) (t <= I_dur) * I_amp;

% Define solver options to limit the timestep
solver_options = odeset('MaxStep', time_step);

% Different initial conditions
initial_conditions_list = {
    [-65.0, 0.3177, 0.0529, 0.5961],
    [-70.0, 0.32, 0.15, 0.58],
    [-60.0, 0.65, 0.05, 0.6],
    [-60.0, 0.33, 0.06, 0.22], 
    [-80.0, 0.33, 0.06, 0.64]
};

current_amplitude = 10; % in uA/cm^2
current_duration = 10; % in ms

simulation_results = struct();

for ic_idx = 1:length(initial_conditions_list)
    initial_conditions = initial_conditions_list{ic_idx};
    % Display initial conditions
    fprintf('Simulating with initial conditions: V0 = %.1f, n0 = %.4f, m0 = %.4f, h0 = %.4f\n', ...
        initial_conditions(1), initial_conditions(2), initial_conditions(3), initial_conditions(4));
    
    % ODE solver with varying current
    [time, states] = ode45(@(t, y) HodgkinHuxleyModel(t, y, mem_cap, max_gK, max_gNa, leak_gL, rev_EK, rev_ENa, rev_EL, applied_current_func(t, current_amplitude, current_duration), alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func), time_span, initial_conditions, solver_options);

    % Store results
    simulation_results(ic_idx).initial_conditions = initial_conditions;
    simulation_results(ic_idx).time = time;
    simulation_results(ic_idx).voltage = states(:, 1);
end

% Plotting the results
figure;
for ic_idx = 1:length(initial_conditions_list)
    subplot(length(initial_conditions_list), 1, ic_idx);
    plot(simulation_results(ic_idx).time, simulation_results(ic_idx).voltage);
    title(sprintf('Initial Conditions: V0 = %.1f, n0 = %.4f, m0 = %.4f, h0 = %.4f', ...
        simulation_results(ic_idx).initial_conditions(1), simulation_results(ic_idx).initial_conditions(2), ...
        simulation_results(ic_idx).initial_conditions(3), simulation_results(ic_idx).initial_conditions(4)));
    xlabel('Time (ms)');
    ylabel('Membrane Potential (mV)');
    grid on;
end

% Hodgkin-Huxley model
function dy = HodgkinHuxleyModel(~, y, mem_cap, max_gK, max_gNa, leak_gL, rev_EK, rev_ENa, rev_EL, I_ext, alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func)
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);

    % Ionic currents
    I_K = max_gK * n^4 * (V - rev_EK);
    I_Na = max_gNa * m^3 * h * (V - rev_ENa);
    I_L = leak_gL * (V - rev_EL);

    % Differential equations
    dVdt = (I_ext - I_K - I_Na - I_L) / mem_cap;
    dndt = alpha_n_func(V) * (1 - n) - beta_n_func(V) * n;
    dmdt = alpha_m_func(V) * (1 - m) - beta_m_func(V) * m;
    dhdt = alpha_h_func(V) * (1 - h) - beta_h_func(V) * h;

    dy = [dVdt; dndt; dmdt; dhdt];
end
