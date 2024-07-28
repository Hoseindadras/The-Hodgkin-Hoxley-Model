% Parameters
Cm = 1.0; % membrane capacitance, in uF/cm^2
gK = 36.0; % maximum conductances, in mS/cm^2
gNa = 120.0; % maximum conductances, in mS/cm^2
gL = 0.3; % leakage conductance, in mS/cm^2
EK = -77.0; % reversal potentials, in mV
ENa = 50.0; % reversal potentials, in mV
EL = -54.4; % reversal potentials, in mV

% Time parameters
tspan = [0 100]; % time range, in ms
dt = 0.01; % desired maximum time step, in ms

% Alpha and Beta functions
alpha_m = @(V) 0.1 * (V + 40) ./ (1 - exp(-(V + 40) / 10));
beta_m = @(V) 4 * exp(-(V + 65) / 18);
alpha_h = @(V) 0.07 * exp(-(V + 65) / 20);
beta_h = @(V) 1 ./ (1 + exp(-(V + 35) / 10));
alpha_n = @(V) 0.01 * (V + 55) ./ (1 - exp(-(V + 55) / 10));
beta_n = @(V) 0.125 * exp(-(V + 65) / 80);

% Initial conditions
V0 = -65.0;
m0 = alpha_m(V0) / (alpha_m(V0) + beta_m(V0));
h0 = alpha_h(V0) / (alpha_h(V0) + beta_h(V0));
n0 = alpha_n(V0) / (alpha_n(V0) + beta_n(V0));
y0 = [V0, n0, m0, h0];

% Define applied current function
I_ext_function = @(t, I_amp, I_dur) (t >= 10 & t <= (10 + I_dur)) * I_amp;

% Define solver options to limit the timestep
options = odeset('MaxStep', dt);

% Investigate effects of varying amplitude and duration of applied current
I_amplitudes = [5, 10, 15, 20]; % in uA/cm^2
I_durations = [5, 10, 20, 30, 40, 50]; % in ms

results = struct();

for amp_idx = 1:length(I_amplitudes)
    for dur_idx = 1:length(I_durations)
        I_amp = I_amplitudes(amp_idx);
        I_dur = I_durations(dur_idx);

        % Display input parameters
        fprintf('Simulating with I_amp = %d uA/cm^2 and I_dur = %d ms\n', I_amp, I_dur);

        % ODE solver with varying current
        [t, y] = ode45(@(t, y) HodgkinHuxleyModel(t, y, Cm, gK, gNa, gL, EK, ENa, EL, I_ext_function(t, I_amp, I_dur), alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n), tspan, y0, options);

        % Store results
        results(amp_idx, dur_idx).I_amp = I_amp;
        results(amp_idx, dur_idx).I_dur = I_dur;
        results(amp_idx, dur_idx).t = t;
        results(amp_idx, dur_idx).V = y(:, 1);
        results(amp_idx, dur_idx).firing_rate = compute_firing_rate(t, y(:, 1));
    end
end

% Plotting the results
figure;
for amp_idx = 1:length(I_amplitudes)
    for dur_idx = 1:length(I_durations)
        subplot(length(I_amplitudes), length(I_durations), (amp_idx-1)*length(I_durations) + dur_idx);
        plot(results(amp_idx, dur_idx).t, results(amp_idx, dur_idx).V);
        title(sprintf('A: %d uA/cm^2, D: %d ms', results(amp_idx, dur_idx).I_amp, results(amp_idx, dur_idx).I_dur));
        xlabel('Time (ms)');
        ylabel('Membrane Potential (mV)');
        grid on;
    end
end
% Hodgkin-Huxley model
function dy = HodgkinHuxleyModel(~, y, Cm, gK, gNa, gL, EK, ENa, EL, I_ext, alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n)
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);

    % Ionic currents
    IK = gK * n^4 * (V - EK);
    INa = gNa * m^3 * h * (V - ENa);
    IL = gL * (V - EL);

    % Differential equations
    dVdt = (I_ext - IK - INa - IL) / Cm;
    dndt = alpha_n(V) * (1 - n) - beta_n(V) * n;
    dmdt = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dhdt = alpha_h(V) * (1 - h) - beta_h(V) * h;

    dy = [dVdt; dndt; dmdt; dhdt];
end

% Compute firing rate
function rate = compute_firing_rate(t, V)
    % Detect spikes
    spikes = find(V(1:end-1) < 0 & V(2:end) >= 0);
    num_spikes = length(spikes);
    
    % Compute firing rate
    rate = num_spikes / (t(end) - t(1)) * 1000; % in Hz
end
