clear;
close all;
clc;

N = 100;
Rb = 100;
Tb = 1/Rb;
fs = 10000;
fc = 1000;
samples_per_bit = fs * Tb;

%% Generate Random Binary Data
data = randi([0 1], 1, N);
data_upsampled = repelem(data, samples_per_bit);

t = (0:length(data_upsampled)-1) / fs;

%% ASK Modulation
carrier = cos(2*pi*fc*t);
ask_signal = data_upsampled .* carrier;

%% Plot Transmitted Signal (Zoomed Time + Frequency)
num_bits_to_plot = 5;
samples_to_plot = samples_per_bit * num_bits_to_plot;

figure;
subplot(2,1,1);
plot(t(1:samples_to_plot), ask_signal(1:samples_to_plot));
title(['ASK Transmitted Signal (First ', num2str(num_bits_to_plot), ' Bits)']);
xlabel('Time [s]'); ylabel('Amplitude'); grid on;

% Frequency Domain
Y = fftshift(fft(ask_signal, 2^nextpow2(length(ask_signal))));
f = linspace(-fs/2, fs/2, length(Y));
subplot(2,1,2);
plot(f, abs(Y));
title('ASK Transmitted Signal (Frequency Domain)');
xlabel('Frequency [Hz]'); ylabel('|Y(f)|'); grid on;

%% Coherent Demodulation with Phase Offsets
phases_deg = [30, 60, 90];

for i = 1:length(phases_deg)
    phase_rad = deg2rad(phases_deg(i));
    osc = cos(2*pi*fc*t + phase_rad);

    received = ask_signal .* osc;

    % Low-pass filter using moving average
    lpf = ones(1, round(samples_per_bit)) / samples_per_bit;
    filtered = conv(received, lpf, 'same');

    % Downsample and make decisions
    sampled = filtered(round(samples_per_bit/2):samples_per_bit:end);
    recovered_bits = sampled > 0.25;

    % Bit error analysis
    num_errors = sum(recovered_bits ~= data);
    disp(['Phase = ', num2str(phases_deg(i)), '°: Bit Errors = ', num2str(num_errors), ' / ', num2str(N)]);

    % Plot Received and Filtered Signal (Zoomed)
    figure;
    subplot(2,1,1);
    plot(t(1:samples_to_plot), received(1:samples_to_plot));
    title(['Received Signal with Phase = ', num2str(phases_deg(i)), '° (Time Domain)']);
    xlabel('Time [s]'); ylabel('Amplitude'); grid on;

    subplot(2,1,2);
    plot(t(1:samples_to_plot), filtered(1:samples_to_plot));
    title(['Filtered Signal with Phase = ', num2str(phases_deg(i)), '°']);
    xlabel('Time [s]'); ylabel('Amplitude'); grid on;
end
