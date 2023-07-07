%% Parameters consideration
% Sample rate: 1000Hz
% Noise voltage: max 10uV (pp), respectively 3uV (RMS) with input short circuit
% Input voltage: +- 16mV
% Resolution: 16-bit with 0.5uV/LSB (2000 A/D units per mV)
% ADC range: +- 16.384mV
%
% Signal to noise ratio: 20log(0.1uV/0.1mV) = -60dB
% Cutoff frequency: 0.5Hz
% Stop frequency: 0.1Hz
%%
name = 's_0010_i.csv';
ECG = readtable(name);

signal_length = 38401;
Fs = 1000;
T = signal_length/Fs;

time = table2array(ECG(2:signal_length,2));
value = table2array(ECG(2:signal_length,3));
%% Spectrum
value_spectrum = fft(value);
F_domain = (0:signal_length-2)/T;
%% Plot original
figure(1), clf
% subplot(2,1,1)
% plot(time,value)
% set(gca,'xlim',[0 38])
% xlabel('Time (s)'), ylabel('Value (mV)')

hold on
subplot(2,1,1)
plot(time,value)
% set(gca,'xlim',[0 2])
xlabel('Time (s)'), ylabel('Value (mV)')
%% Plot spectrum
% figure(2)
hold on
subplot(2,1,2)
plot(F_domain,2*abs(value_spectrum)/signal_length)
set(gca,'xlim',[0 180],'ylim',[-0.01 Inf])
xlabel('Frequency (Hz)'), ylabel('Amplitude (mVpp)')
%% Median filter
hann_length = 3000;
hann_window = hann(hann_length)/sum(hann(hann_length));
value_base = conv(value, hann_window, 'same');
p = polyfit(time,value,20);
value_base_polyfit = polyval(p,time);

figure(3), clf
subplot(2,1,1)
plot(time,value)
hold on
plot(time,value_base, 'LineWidth', 4)
title('Estimate baseline noise', 'FontSize', 24)
% hold on
% plot(time,value_base_polyfit, 'LineWidth', 4)
legend('signal', 'window')
% set(gca, 'xlim',[0 5])
%% Baseline spectrum
Tb = 38.4;
N_freqb = (0:1600)/16;
F_domainb = N_freqb/Tb;

value_base_fft = fft(value_base);
% value_base_fft_b = fourier.spectrum(time,value_base,N_freqb,Tb/2);

% figure(4), clf
figure(3), hold on, subplot(2,1,2)
plot(F_domain, abs(value_base_fft)/signal_length)
% hold on
% plot(F_domainb, value_base_fft_b)
set(gca, 'xlim',[0 2])
xlabel('Frequency (Hz)'), ylabel('Amplitude (mVpp)')
title('Baseline spectrum (Window estimate)', 'FontSize', 24)
%% Butterworth filter (high pass)
cutoff_freq = 0.1;            % Cutoff frequency (Hz)
% passband_ripple = 0.1;      % Maximum passband ripple in dB

% Design the Butterworth high-pass filter
order = 3;  % Filter order
[b_base, a_base] = butter(order, cutoff_freq/(Fs/2), 'high');
% [z,p,k] = butter(order,cutoff_freq/(Fs/2),'high');

figure(2), clf
freqz(b_base,a_base,20*Fs,Fs)
% sos = zp2sos(z,p,k);
% fvtool(sos,'Analysis','freq')
set(gca,'xlim',[0 2])
title('Butterworth frequency response','FontSize',15)
%% Low pass with FIR
fcutoff = 150;
transw  = 0.2;
order   = round( 24*Fs/fcutoff );
npnts   = length(time);
hz = linspace(0,Fs/2,floor(npnts/2)+1);

shape   = [ 1 1 0 0 ];
frex    = [ 0 fcutoff fcutoff+fcutoff*transw Fs/2 ] / (Fs/2);

% filter kernel
filtkern = firls(order,frex,shape);

% its power spectrum
filtkernX = abs(fft(filtkern,npnts)).^2;

figure(5), clf
subplot(211), hold on
plot(frex*Fs/2,shape,'r','linew',1)
plot(hz,filtkernX(1:length(hz)),'k','linew',2)
set(gca,'xlim',[0 300])
xlabel('Frequency (Hz)'), ylabel('Gain')
title('Filter kernel spectrum')

subplot(212), hold on
plot(hz,10*log(filtkernX(1:length(hz))),'k','linew',2)
set(gca,'xlim',[0 300])
xlabel('Frequency (Hz)'), ylabel('Gain (dBm)')
title('Filter kernel spectrum')
%% Notch filter butterworth
[b_notch,a_notch] = butter(3,[49.5 50.5]/(Fs/2),'stop');
figure(7), clf
freqz(b_notch,a_notch,Fs*10,Fs)
set(gca,'xlim',[40 60])
title('50Hz notch filter','FontSize',16)
%% Low pass (100Hz with Window)
[b_low,a_low] = butter(3,100/(Fs/2),'low');
figure(8), clf
freqz(b_notch,a_notch,Fs*10,Fs)
% set(gca,'xlim',[40 60])
%% Let's filter
%% Section 1: Remove noise
filter_value = value;
filter_value = filter(b_base,a_base,filter_value);
filter_value = filter(b_notch,a_notch,filter_value);
filter_value = filtfilt(filtkern,1,filter_value);
% filter_value = conv(value, hann_window_100Hz, 'same');

filter_value_base = conv(filter_value,hann_window,'same');

figure(10), clf
subplot(211)
plot(time,filter_value)
% hold on
% plot(time,filter_value_base,'LineWidth',4)
set(gca,'xlim',[0 2.5])
title('Filtered')
% legend('signal','baseline')
xlabel('Time (s)'), ylabel('Value (mV)')

hold on, subplot(212)
plot(time,value)
% hold on
% plot(time,value_base,'LineWidth',4)
set(gca,'xlim',[0 2.5])
title('Original')
% legend('signal','baseline')
xlabel('Time (s)'), ylabel('Value (mV)')

%% Section 2: Spectrum
filter_value_spectrum = fft(filter_value);

figure(11), clf
subplot(212)
plot(F_domain,2*abs(value_spectrum)/signal_length)
set(gca,'xlim',[-1.5 100],'ylim',[-0.003 0.07])
xlabel('Frequency (Hz)'), ylabel('Amplitude (mVpp)')
title('Original','FontSize',15)

hold on
subplot(211)
plot(F_domain,2*abs(filter_value_spectrum)/signal_length)
set(gca,'xlim',[-1.5 100],'ylim',[-0.003 0.07])
xlabel('Frequency (Hz)'), ylabel('Amplitude (mVpp)')
title('Filtered','FontSize',15)