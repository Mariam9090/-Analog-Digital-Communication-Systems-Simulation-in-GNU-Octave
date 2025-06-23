clc; 
clear; 
close all;

n = 64;                         
fs = 1000;  
ts = 1/fs;
bitrate = 1;                   
Tb = 1/bitrate; 
N = ceil((Tb*n)/ts);
samples_per_bit = fs / bitrate;
t = 0 : ts : (N-1)*ts;   

bits = randi([0 1], 1, n);
disp('Random bitstream:');
disp(bits);

manchester = zeros(size(t));
unipolar = zeros(size(t));

for i = 1:n
    idx = (i-1)*samples_per_bit + 1 : i*samples_per_bit;
    half = floor(samples_per_bit / 2);

    % manchester
    if bits(i) == 1
        manchester(idx(1:half)) = 1;
        manchester(idx(half+1:end)) = -1;
    else
        manchester(idx(1:half)) = -1;
        manchester(idx(half+1:end)) = 1;
    end

    % unipolar NRZ
    if bits(i) == 1
        unipolar(idx) = 1;
    else
        unipolar(idx) = 0;
    end
end

%% plot manchester
figure;
subplot(2,1,1);
plot(t, manchester, 'b', 'LineWidth', 2);
title('Manchester Line Coding');
xlabel('Time (s)'); ylabel('Amplitude');
xlim([0 64]);
ylim([-1.5 1.5]);
grid on;

%% plot unipolar NRZ
subplot(2,1,2);
plot(t, unipolar, 'r', 'LineWidth', 2);
title('Unipolar NRZ Line Coding');
xlabel('Time (s)'); ylabel('Amplitude');
xlim([0 64]);
ylim([-0.5 1.5]);
grid on;

%% frequency domain
df = bitrate/n;
if(rem(N,2)==0) %% Even
  f = - (0.5*fs) : df : (0.5*fs-df) ; 
else %% Odd
  f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df) ;
end

M = abs(fftshift(fft(manchester))* ts);
U = abs(fftshift(fft(unipolar))* ts);

figure;
subplot(2,1,1);
plot(f, M, 'b');
xlim([-10 10]);
title('Spectrum of Manchester Code');
xlabel('Frequency (Hz)'); ylabel('|X(f)|'); grid on;

subplot(2,1,2);
plot(f, U, 'r');
xlim([-10 10]);
title('Spectrum of Unipolar NRZ Code');
xlabel('Frequency (Hz)'); ylabel('|X(f)|'); grid on;