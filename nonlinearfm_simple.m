function [x_nl,ft] = nonlinearfm(x,fs,B)
% x : LFM waveform, fs : Sample Frequency (Hz), B : Bandwidth (Hz)
N = numel(x); T=N/fs; df = fs/N; dt = 1/fs;
f = [-fs/2:df:fs/2-df];
t = -T/2:dt:(T/2-dt);
% Compute Spectrum
xfft = fftshift(fft(x));
xfft = xfft(abs(f)<=(B/2));
% Compute Group Time Delay function
Tg = cumsum(abs(xfft).^2);
% Solve Boundary Conditions
c1 = T/(Tg(end)-Tg(1));
c2 = -T/2-c1*Tg(1);
Tg = c1*Tg+c2;
% Invert Tg using interpolation
ft = interp1(Tg, f(abs(f)<=(B/2)), t);
% Integrate frequency function to obtain phase
phi = 2*pi*cumsum(ft)/fs;
phi = phi+pi/2-phi(1);
% Compute nonlinear FM waveform
x_nl = exp(1i*phi);

