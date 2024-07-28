% Parameters
membrane_capacitance = 1.0; % membrane capacitance, in uF/cm^2
max_gK = 36.0; % maximum conductances, in mS/cm^2
max_gNa = 120.0; % maximum conductances, in mS/cm^2
leakage_gL = 0.3; % leakage conductance, in mS/cm^2
reversal_EK = -77.0; % reversal potentials, in mV
reversal_ENa = 50.0; % reversal potentials, in mV
reversal_EL = -54.4; % reversal potentials, in mV

% Time parameters
time_range = [0 15]; % time range, in ms
time_step = 0.000001; % time step, in ms

% Alpha and Beta functions
alpha_m_func = @(voltage) 0.1 * (voltage + 40) ./ (1 - exp(-(voltage + 40) / 10));
beta_m_func = @(voltage) 4 * exp(-(voltage + 65) / 18);
alpha_h_func = @(voltage) 0.07 * exp(-(voltage + 65) / 20);
beta_h_func = @(voltage) 1 ./ (1 + exp(-(voltage + 35) / 10));
alpha_n_func = @(voltage) 0.01 * (voltage + 55) ./ (1 - exp(-(voltage + 55) / 10));
beta_n_func = @(voltage) 0.125 * exp(-(voltage + 65) / 80);

% Initial conditions
initial_voltage = -65.0;
initial_m = alpha_m_func(initial_voltage) / (alpha_m_func(initial_voltage) + beta_m_func(initial_voltage));
initial_h = alpha_h_func(initial_voltage) / (alpha_h_func(initial_voltage) + beta_h_func(initial_voltage));
initial_n = alpha_n_func(initial_voltage) / (alpha_n_func(initial_voltage) + beta_n_func(initial_voltage));
initial_conditions = [initial_voltage, initial_n, initial_m, initial_h];

% Define threshold detection function
detect_action_potential = @(voltage) any(voltage > 0); % Define threshold as 0 mV

% Find threshold current and store voltage responses
current_start = 0; % starting current, in uA/cm^2
current_end = 5; % ending current, in uA/cm^2
current_step = 0.1; % current step, in uA/cm^2

threshold_current = current_start;
found_action_potential = false;
voltage_responses = {};

for external_current = current_start:current_step:current_end
    % ODE solver
    [time, state_variables] = ode45(@(t, y) HodgkinHuxleyModel(t, y, membrane_capacitance, max_gK, max_gNa, leakage_gL, reversal_EK, reversal_ENa, reversal_EL, external_current, alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func), time_range, initial_conditions);
    
    % Store voltage response
    voltage_responses{end+1} = struct('external_current', external_current, 'time', time, 'voltage', state_variables(:, 1));
    
    % Check if action potential occurs
    if detect_action_potential(state_variables(:, 1)) && ~found_action_potential
        threshold_current = external_current;
        found_action_potential = true;
    end
end

if found_action_potential
    fprintf('Threshold current is: %.2f uA/cm^2\n', threshold_current);
else
    fprintf('No action potential elicited within the given range of currents.\n');
end

% Plotting the voltage responses for some excitations
figure;
hold on;
plot_colors = lines(5);
for i = 1:10:length(voltage_responses)
    plot_data = voltage_responses{i};
    plot(plot_data.time, plot_data.voltage, 'DisplayName', ['I_{ext} = ', num2str(plot_data.external_current), ' uA/cm^2']);
end
title('Voltage Response Over Time for Various External Currents');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
legend('show');
grid on;
hold off;

% Hodgkin-Huxley model
function dy = HodgkinHuxleyModel(~, y, membrane_capacitance, max_gK, max_gNa, leakage_gL, reversal_EK, reversal_ENa, reversal_EL, external_current, alpha_m_func, beta_m_func, alpha_h_func, beta_h_func, alpha_n_func, beta_n_func)
    voltage = y(1);
    n = y(2);
    m = y(3);
    h = y(4);

    % Ionic currents
    potassium_current = max_gK * n^4 * (voltage - reversal_EK);
    sodium_current = max_gNa * m^3 * h * (voltage - reversal_ENa);
    leakage_current = leakage_gL * (voltage - reversal_EL);

    % Differential equations
    dVdt = (external_current - potassium_current - sodium_current - leakage_current) / membrane_capacitance;
    dndt = alpha_n_func(voltage) * (1 - n) - beta_n_func(voltage) * n;
    dmdt = alpha_m_func(voltage) * (1 - m) - beta_m_func(voltage) * m;
    dhdt = alpha_h_func(voltage) * (1 - h) - beta_h_func(voltage) * h;

    dy = [dVdt; dndt; dmdt; dhdt];
end
