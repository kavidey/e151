 % genBodePlot.m
 % David E Olumese (dolumese@g.hmc.edu) | 28th Jan 2019
 % E151/E153 Course Development

%% variables
 % data file parameters
finname = './data/scope_13.csv';
foutname = './data/scope_12.csv';
fStartRow = 3; 
fStartCol = 0;

 % estimated number of system poles & zeros
np = 1;
nz = 0;

% maximum frequency in sweep
Fmin_Hz  = 10^5;     % [Hz]
Fmax_Hz  = 1.5*10^6; % [Hz]

 % bode plot figure bounds
mindB = -20;
maxdB = 10;
minDeg = -135;
maxDeg = 0;

%% Pull data from the file
M = csvread(finname, fStartRow, fStartCol);
tIn    = M(:,1);
inV  = M(:,2);
inV = inV(tIn > 0 & tIn < 10^(-3));
tIn = tIn(tIn > 0 & tIn < 10^(-3));
% triggerInV = M(:,3);

M = csvread(foutname, fStartRow, fStartCol);
tOut    = M(:,1);
outV  = M(:,2);
outV = outV(tOut > 0 & tOut < 10^(-3));
tOut = tOut(tOut > 0 & tOut < 10^(-3));
% triggerOutV = M(:,3);

Ts = tIn(2);  % [s] Sampling period
Fs = 1/Ts;  % [Hz] Sampling frequency
windw = []; % windowing function
%% Calculate Frequencies
FRange_Hz = logspace(log10(Fmin_Hz), log10(Fmax_Hz), size(tIn, 1));
Fmax      = Fmax_Hz*(2*pi); % [rad/s]
%% Plot data
figure(1)
plot(tIn, inV, tIn, outV);
% plot(tIn, triggerInV, tOut, triggerOutV);
legend('Input signal', 'Output signal')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Input and Output Signals from Frequency Sweep')

%% Make bode plot with estimated poles & zeros
 % generate an estimated system
data = iddata(outV, inV, Ts);
sys = tfest(data, np, nz);

 % plot the bode plot
figure(2)
H = bodeplot(sys);
setoptions(H, 'FreqUnits', 'Hz');
title('Estimated System Bode Plot (Based on Number of Poles and Zeros)')
grid on

%% Make bode plot using transfer function estimation
[txy, ft] = tfestimate(inV, outV, windw, [], FRange_Hz, Fs); % generate tf estimate

 % determine system parameters
A  = abs(txy);           A_dB = mag2db(A);
Ph = unwrap(angle(txy)); Ph_deg = 180/pi*Phz;
w  = 2*pi*ft;

 % plot bode (include only useful information [F0, Fmax])
figure(3);
subplot(2, 1, 1)
semilogx(ft, A_dB);
axis([0 Fmax_Hz mindB maxdB])
title('System Bode Plot from Estimated Transfer Function')
ylabel('Magnitude (dB)')
grid on

subplot(2, 1, 2)
semilogx(ft, Ph_deg);
set(gca, 'YTick', [-180 -90 -45 0 45 90 180])
axis([0 Fmax_Hz minDeg maxDeg])
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
grid on