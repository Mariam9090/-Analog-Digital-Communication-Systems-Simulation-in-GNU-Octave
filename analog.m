clc;
clear
close all
%% 1. plot the signal x(t)
fs = 100;
ts = 1/fs;
df = 0.01;
T_m = 1/df;
N = ceil(T_m/ts);
t = -((N-1)*ts)/2:ts:((N-1)*ts)/2;
x = zeros(size(t));

leftside = (t > -4) & (t < 0);
x(leftside) = t(leftside) + 5;
rightside = (t >= 0) & (t < 4);
x(rightside) = 5 - t(rightside);

figure;
plot(t, x, 'LineWidth',1.5);
xlabel('t (s)');
xlim([-5 5]);
ylabel('x(t)');
title('Time Domain Signal x(t)');
grid on;
%% 2. analytical expression fourier transform X(f) to plot later
if(rem(N,2)==0)
  f = - (0.5*fs) : df : (0.5*fs-df) ;
else %% Odd
  f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df) ;
end
X_analytic = 16*((sinc(4*f)).^2) + 8*sinc(8*f);
X_analytic(f == 0) = 24;

%% 3. calculate fourier transform
X = fftshift(fft(x))*ts;
figure;
plot(f, abs(X_analytic), 'LineWidth',1.5);
hold on;
plot(f, abs(X), '--','LineWidth',1);
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
legend('Analytical','Numerical');
title('Analytical vs. Numerical Spectra');
xlim([-2, 2]);
grid on;

%% 4. estimate the BW
Power = abs(X_analytic).^2;
Pmax = max(Power);
bandi = find(Power >= 0.05 * Pmax);
index = find(Power == Pmax);
BW_est = max(abs(f(bandi)));
disp(['Estimated BW of x(t): ', num2str(BW_est), ' Hz']);

%% 5. ideal LPF with BW = 1 Hz
H1 = abs(f) <= 1;
Y1 = X .* H1;
y1 = real(ifft(ifftshift(Y1) /ts));

figure;
plot(t, x, 'b','LineWidth',1.5);
hold on;
plot(t, y1, 'r--','LineWidth',1.5);
xlabel('Time t (s)');
xlim([-7 7]);
ylabel('Amplitude');
legend('Original','Filtered BW=1Hz');
title('LPF Output (BW=1 Hz)');
grid on;

%% 6. ideal LPF with BW = 0.3 Hz
H2 = abs(f) <= 0.3;
Y2 = X .* H2;
y2 = real(ifft(ifftshift(Y2) /ts));
figure;
plot(t, x, 'b','LineWidth',1.5);
hold on;
plot(t, real(y2(1:length(t))), 'r--','LineWidth',1.5);
xlabel('Time t (s)');
xlim([-10 10]);
ylabel('Amplitude');
legend('Original','Filtered BW=0.3Hz');
title('LPF Output (BW=0.3 Hz)');
grid on;
%% 7.1) define the signal m(t)
m = zeros(size(t));
fm = 0.5;
m(rightside) = cos(2*pi*fm*t(rightside));
figure;
plot(t, m, 'LineWidth',1.5);
xlabel('t (s)');
xlim([-1 5]);
ylabel('m(t)');
ylim([-1.5 1.5]);
title('Time Domain Signal m(t)'); grid on;
%% 7.2) analytical expression fourier transform M(f) to plot later
T_m = 4;
M_analytic = 2 * ( sinc((f - fm) * T_m) .* exp(-1j * pi * (f - fm) * T_m) + ...
    sinc((f + fm) * T_m) .* exp(-1j * pi * (f + fm) * T_m) );

%% 7.3) calculate fourier transform of m(t)
M = fftshift(fft(m))*ts;
figure;
plot(f, abs(M_analytic), 'LineWidth',1.5);
hold on;
plot(f, abs(M), '--','LineWidth',1);
xlabel('Frequency (Hz)');
ylabel('|m(f)|');
legend('Analytical','Numerical');
title('Analytical vs. Numerical Spectra');
xlim([-5,5]);
grid on;

%% 7.4) estimate the BW of m(t)
Power_m = abs(M_analytic).^2;
Pm = max(Power_m);
bandm = abs(f(Power_m>=0.05*Pm));
BW_m = max(bandm);
fprintf('Estimated BW of m(t): %.2f Hz\n',BW_m);

%% 8. FDM modulation (x(t) in DSB-SC, m(t) in SSB)
fc1 = 20;
c1 = cos(2*pi*fc1*t);
s1 = y1 .* c1;
fc2 = 23; %% fc1(20) + BW(1) + guard band(2)
c2 = 2 * cos(2*pi*fc2*t);
H = (f>fc2 & f<(fc2+2)) | (f<-fc2 & f>-(fc2+2)); %% USB
s2 = m .* c2;
S2 = fftshift(fft (s2)) *ts;
S2 = S2.*H;
s2 = real(ifft(ifftshift(S2) / ts));
s = s1 + s2;

figure;
plot(t, s, 'LineWidth',1.2);
xlabel('t (s)');
ylabel('s(t)');
title('Composite FDM Signal');
grid on;

S = fftshift(fft (s)) *ts;
figure;
plot(f, abs(S), 'LineWidth',1.2);
xlabel('Frequency (Hz)');
ylabel('|S(f)|');
title('Spectrum of FDM Signal');
xlim([-40,40]);
grid on;

%% 12. coherent demodulation
xr_unfiltered = 2 * s .* cos(2*pi*fc1*t);
Xr_unfiltered = fftshift(fft (xr_unfiltered)) *ts;
H_x = double(abs(f) <= (2));
Xr = H_x .* Xr_unfiltered;
x_rec = real(ifft(ifftshift(Xr) /ts));

mr_unfiltered = 2 * s .* cos(2*pi*fc2*t);
Mr_unfiltered = fftshift(fft (mr_unfiltered)) *ts;
H_m = double(abs(f) <= (2));
Mr = H_m .* Mr_unfiltered;
m_rec_full = real(ifft(ifftshift(Mr) /ts));
m_rec = m_rec_full(1:length(t));

figure;
subplot(2,1,1);
plot(t, x, 'b', t, x_rec, 'r--');
legend('x(t)','Recovered x');
xlabel('t (s)');
xlim([-6 6]);
ylim([-1 5.5]);
title('x(t) Recovery');
grid on;
subplot(2,1,2);
plot(t, m, 'b', t, m_rec, 'r--');
legend('m(t)','Recovered m');
xlabel('t (s)');
xlim([-1 5]);
ylim([-1.5 1.5]);
title('m(t) Recovery'); grid on;

