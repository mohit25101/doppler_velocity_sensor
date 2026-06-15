clc;
clear;
close all;

%% Radar Parameters

fc = 10e9;          % 10 GHz carrier
A = 1;

%% Simulation Parameters

fs = 1e6;           % Sampling frequency
T = 0.01;           % Observation time

t = 0:1/fs:T;

%% Transmitted Signal

tx = A*cos(2*pi*fc*t);

%% Plot

figure;
plot(t(1:1000),tx(1:1000));
xlabel('Time (s)');
ylabel('Amplitude');
title('Transmitted Signal');
grid on;

%% Aircraft Parameters

V = 100;                % m/s

theta = 60*pi/180;      % beam angle

c = 3e8;                % speed of light

lambda = c/fc;

fd = (2*V*cos(theta))/lambda;

disp(['Doppler Frequency = ',num2str(fd),' Hz'])

%% Received Signal

rx = A*cos(2*pi*(fc+fd)*t);

figure;
plot(t(1:1000),rx(1:1000));

xlabel('Time');
ylabel('Amplitude');

title('Received Signal');
grid on;

%% Mixing

beat = tx .* rx;

figure;
plot(t(1:5000),beat(1:5000));

title('Mixer Output');
xlabel('Time');
ylabel('Amplitude');
grid on;

%% FFT

N = length(beat);

FFT_BEAT = fft(beat);

freq = (0:N-1)*(fs/N);

figure;

plot(freq,abs(FFT_BEAT));

xlim([0 10000]);

xlabel('Frequency (Hz)');
ylabel('Magnitude');

title('FFT Spectrum');
grid on;

%% Peak Detection

[~,index] = max(abs(FFT_BEAT));

fd_est = freq(index);

disp(['Estimated Doppler Frequency = ', ...
    num2str(fd_est),' Hz']);

%% Velocity Estimation

V_est = (fd_est*lambda)/(2*cos(theta));

disp(['Estimated Velocity = ', ...
    num2str(V_est),' m/s']);

%% AWGN

SNR = 10;

rx_noisy = awgn(rx,SNR,'measured');

%realistic noisy signal
figure;
plot(t(1:1000), rx_noisy(1:1000));
title('Realistic (Noisy) Received Signal');

beat_noisy = tx .* rx_noisy;
FFT_NOISY = fft(beat_noisy);

figure;

plot(freq,abs(FFT_NOISY));

xlim([0 10000]);

title('FFT with Noise');

xlabel('Frequency');

ylabel('Magnitude');

grid on;

%% Aircraft True Motion

Vx_true = 100;      % m/s
Vy_true = 20;       % m/s

%% Radar

fc = 10e9;
c = 3e8;
lambda = c/fc;

gamma = 60*pi/180;
beta = 45*pi/180;


T = 20;

dt = 0.1;

t = 0:dt:T;

Vx_true = 100 + 20*sin(0.2*t);

Vy_true = 15*cos(0.1*t);

Vz_true = 25*cos(0.05*t);
%% Altitude Profile

Altitude = zeros(size(t));

Altitude(1) = 1000;

for k = 2:length(t)

    Altitude(k) = Altitude(k-1) + ...
        Vz_true(k-1)*dt;

end

for k=1:length(t)

    fd1(k) = (2/lambda)*( ...
          Vx_true(k)*sin(gamma)*cos(beta)+ ...
          Vy_true(k)*sin(gamma)*sin(beta)+ ...
          Vz_true(k)*cos(gamma));

fd2(k) = (2/lambda)*( ...
          Vx_true(k)*sin(gamma)*cos(beta)- ...
          Vy_true(k)*sin(gamma)*sin(beta)+ ...
          Vz_true(k)*cos(gamma));

fd3(k) = (2/lambda)*( ...
         -Vx_true(k)*sin(gamma)*cos(beta)+ ...
          Vy_true(k)*sin(gamma)*sin(beta)+ ...
          Vz_true(k)*cos(gamma));

fd4(k) = (2/lambda)*( ...
         -Vx_true(k)*sin(gamma)*cos(beta)- ...
          Vy_true(k)*sin(gamma)*sin(beta)+ ...
          Vz_true(k)*cos(gamma));

end

for k=1:length(t)

    Vx_est(k)=A*(fd1(k)+fd2(k));

    Vy_est(k)=A*(fd1(k)-fd2(k));

end

figure;

plot(t,Vx_true,'LineWidth',2);
hold on;

plot(t,Vx_est,'--','LineWidth',2);

legend('True','Estimated');

xlabel('Time');
ylabel('Velocity');

title('Forward Velocity');
grid on;


%% Altitude Plot

figure;

plot(t,Altitude,'LineWidth',2);

grid on;

xlabel('Time (s)');
ylabel('Altitude (m)');

title('Aircraft Altitude Profile');

SNR_dB = 10;

%% Altitude Dependent SNR Parameters

SNR_ref = 20;

Altitude_ref = 1000;

for k = 1:length(t)

    signal_power1 = fd1(k)^2;
    noise_power1 = signal_power1/(10^(SNR_dB/10));
    sigma1 = sqrt(noise_power1);

    signal_power2 = fd2(k)^2;
    noise_power2 = signal_power2/(10^(SNR_dB/10));
    sigma2 = sqrt(noise_power2);

    % fd1_m(k) = fd1(k) + sigma1*randn;
    % fd2_m(k) = fd2(k) + sigma2*randn;
        fd1_true = fd1(k);
        fd2_true = fd2(k);
        fd3_true = fd3(k);
        fd4_true = fd4(k);
        SNR_dB = SNR_ref ...
            - 40*log10(Altitude(k)/Altitude_ref);

        SNR_hist(k) = SNR_dB;
        t_local = 0:1/fs:0.01;
        rx1 = exp(1j*2*pi*fd1_true*t_local);

rx2 = exp(1j*2*pi*fd2_true*t_local);

rx3 = exp(1j*2*pi*fd3_true*t_local);

rx4 = exp(1j*2*pi*fd4_true*t_local);

        rx1_noisy = awgn(rx1,SNR_dB,'measured');
        rx2_noisy = awgn(rx2,SNR_dB,'measured');
        rx3_noisy = awgn(rx3,SNR_dB,'measured');
        rx4_noisy = awgn(rx4,SNR_dB,'measured');

        FFT1 = fftshift(abs(fft(rx1_noisy)));

FFT2 = fftshift(abs(fft(rx2_noisy)));

FFT3 = fftshift(abs(fft(rx3_noisy)));

FFT4 = fftshift(abs(fft(rx4_noisy)));
        N = length(FFT1);

        f_axis = (-N/2:N/2-1)*(fs/N);
        [~,idx1] = max(FFT1);
        fd1_est = f_axis(idx1);

        [~,idx2] = max(FFT2);
        fd2_est = f_axis(idx2);

        [~,idx3] = max(FFT3);
        fd3_est = f_axis(idx3);

        [~,idx4] = max(FFT4);
        fd4_est = f_axis(idx4);

        fd1_fft(k) = fd1_est;
        fd2_fft(k) = fd2_est;
        fd3_fft(k)=fd3_est;
        fd4_fft(k)=fd4_est;
        if k==1

            fprintf('\n');

            fprintf('fd1_est = %.2f Hz\n',fd1_est);
            fprintf('fd2_est = %.2f Hz\n',fd2_est);
            fprintf('fd3_est = %.2f Hz\n',fd3_est);
            fprintf('fd4_est = %.2f Hz\n',fd4_est);

        end

        %% 3D Velocity Solver

       B = lambda/(8*sin(gamma)*cos(beta));

C = lambda/(8*sin(gamma)*sin(beta));

D = lambda/(8*cos(gamma));

Vx_meas(k)=B*( ...
              fd1_est + fd2_est ...
             -fd3_est - fd4_est );

Vy_meas(k)=C*( ...
              fd1_est - fd2_est ...
             +fd3_est - fd4_est );

Vz_meas(k)=D*( ...
              fd1_est + fd2_est ...
             +fd3_est + fd4_est );

end
fprintf('\nMeasured Vx (first sample) = %.2f\n',Vx_meas(1));
fprintf('Measured Vy (first sample) = %.2f\n',Vy_meas(1));
%% Dynamic SNR Plot

figure;

plot(t,SNR_hist,'LineWidth',2);

grid on;

xlabel('Time (s)');
ylabel('SNR (dB)');

title('Altitude Dependent SNR');


%EKF

F = eye(3);

H = eye(3);

xhat = [Vx_meas(1);
        Vy_meas(1);
        Vz_meas(1)];

P = eye(3);

Q = 0.5*eye(3);

R = 2*eye(3);


for k=1:length(t)

    z = [Vx_meas(k);
        Vy_meas(k);
        Vz_meas(k)];

    %% Prediction

    xpred = F*xhat;

    Ppred = F*P*F' + Q;

    %% Gain

    K = Ppred*H'/(H*Ppred*H'+R);

    %% Update

    xhat = xpred + K*(z-H*xpred);

    P = (eye(3)-K*H)*Ppred;

    Vx_EKF(k)=xhat(1);
    Vy_EKF(k)=xhat(2);
    Vz_EKF(k)=xhat(3);

end

%% Altitude Estimation

Altitude_EKF = zeros(size(t));

Altitude_EKF(1)=Altitude(1);

for k=2:length(t)

    Altitude_EKF(k)= ...
        Altitude_EKF(k-1)+ ...
        Vz_EKF(k-1)*dt;

end


%% Doppler Frequency Error

Error_fd1 = fd1 - fd1_fft;
Error_fd2 = fd2 - fd2_fft;

RMSE_fd1 = sqrt(mean(Error_fd1.^2));
RMSE_fd2 = sqrt(mean(Error_fd2.^2));

fprintf('\n');
fprintf('fd1 RMSE = %.3f Hz\n',RMSE_fd1);
fprintf('fd2 RMSE = %.3f Hz\n',RMSE_fd2);
figure;

plot(t,Error_fd1,'LineWidth',2);
grid on;

xlabel('Time (s)');
ylabel('Frequency Error (Hz)');
title('Beam 1 Doppler Error');

figure;

plot(t,Error_fd2,'LineWidth',2);
grid on;

xlabel('Time (s)');
ylabel('Frequency Error (Hz)');
title('Beam 2 Doppler Error');

%%%%%%%

figure;

plot(t,Vx_true,'k','LineWidth',2);
hold on;

plot(t,Vx_meas,'r.');

plot(t,Vx_EKF,'b','LineWidth',2);

legend('True','Measured','EKF');

grid on;

figure;

plot(t,Vy_true,'k','LineWidth',2);
hold on;

plot(t,Vy_meas,'r.');

plot(t,Vy_EKF,'b','LineWidth',2);

legend('True','Measured','EKF');

xlabel('Time (s)');
ylabel('Vy (m/s)');
title('Lateral Velocity Estimation');

grid on;

figure;

plot(t,Vz_true,'k','LineWidth',2);
hold on;

plot(t,Vz_meas,'r.');

plot(t,Vz_EKF,'b','LineWidth',2);

legend('True','Measured','EKF');

xlabel('Time (s)');
ylabel('Vz (m/s)');

title('Vertical Velocity Estimation');

grid on;

%% RMSE

RMSE_Vx_Meas = sqrt(mean((Vx_true - Vx_meas).^2));
RMSE_Vx_EKF  = sqrt(mean((Vx_true - Vx_EKF).^2));

RMSE_Vy_Meas = sqrt(mean((Vy_true - Vy_meas).^2));
RMSE_Vy_EKF  = sqrt(mean((Vy_true - Vy_EKF).^2));

fprintf('\n');
fprintf('Vx Measurement RMSE = %.3f m/s\n',RMSE_Vx_Meas);
fprintf('Vx EKF RMSE         = %.3f m/s\n',RMSE_Vx_EKF);

fprintf('\n');

fprintf('Vy Measurement RMSE = %.3f m/s\n',RMSE_Vy_Meas);
fprintf('Vy EKF RMSE         = %.3f m/s\n',RMSE_Vy_EKF);

RMSE_Vz_Meas = ...
    sqrt(mean((Vz_true-Vz_meas).^2));

RMSE_Vz_EKF = ...
    sqrt(mean((Vz_true-Vz_EKF).^2));

fprintf('\n');

fprintf('Vz Measurement RMSE = %.3f m/s\n', ...
    RMSE_Vz_Meas);

fprintf('Vz EKF RMSE = %.3f m/s\n', ...
    RMSE_Vz_EKF);

% Ground speed tracking
GS_true = sqrt(Vx_true.^2 + Vy_true.^2);

GS_EKF = sqrt(Vx_EKF.^2 + Vy_EKF.^2);

figure;

plot(t,GS_true,'k','LineWidth',2);
hold on;

plot(t,GS_EKF,'b','LineWidth',2);

legend('True','EKF');

xlabel('Time (s)');
ylabel('Ground Speed (m/s)');

title('Ground Speed Estimation');

grid on;

%%Altitude plot

figure;

plot(t,Altitude,'k','LineWidth',2);
hold on;

plot(t,Altitude_EKF,'b','LineWidth',2);

legend('True','Estimated');

xlabel('Time (s)');
ylabel('Altitude (m)');

title('Altitude Estimation');

grid on;

%Drift angle tracking
Drift_true = atan2d(Vy_true,Vx_true);

Drift_EKF = atan2d(Vy_EKF,Vx_EKF);

figure;

plot(t,Drift_true,'k','LineWidth',2);
hold on;

plot(t,Drift_EKF,'b','LineWidth',2);

legend('True','EKF');

xlabel('Time (s)');
ylabel('Drift Angle (deg)');

title('Drift Angle Estimation');

grid on;

%quantitative accuracy measurement
%% RMSE Analysis

RMSE_Vx_Meas = sqrt(mean((Vx_true - Vx_meas).^2));
RMSE_Vx_EKF  = sqrt(mean((Vx_true - Vx_EKF).^2));

RMSE_Vy_Meas = sqrt(mean((Vy_true - Vy_meas).^2));
RMSE_Vy_EKF  = sqrt(mean((Vy_true - Vy_EKF).^2));

GS_true = sqrt(Vx_true.^2 + Vy_true.^2);
GS_meas = sqrt(Vx_meas.^2 + Vy_meas.^2);
GS_EKF  = sqrt(Vx_EKF.^2 + Vy_EKF.^2);

RMSE_GS_Meas = sqrt(mean((GS_true - GS_meas).^2));
RMSE_GS_EKF  = sqrt(mean((GS_true - GS_EKF).^2));

fprintf('\n========== RMSE RESULTS ==========\n');

fprintf('Vx Measurement RMSE = %.3f m/s\n',RMSE_Vx_Meas);
fprintf('Vx EKF RMSE         = %.3f m/s\n\n',RMSE_Vx_EKF);

fprintf('Vy Measurement RMSE = %.3f m/s\n',RMSE_Vy_Meas);
fprintf('Vy EKF RMSE         = %.3f m/s\n\n',RMSE_Vy_EKF);

fprintf('Ground Speed Measurement RMSE = %.3f m/s\n',RMSE_GS_Meas);
fprintf('Ground Speed EKF RMSE         = %.3f m/s\n',RMSE_GS_EKF);


%Error plots
Error_Vx = Vx_true - Vx_EKF;
Error_Vy = Vy_true - Vy_EKF;

figure;
plot(t,Error_Vx,'LineWidth',2);
grid on;
xlabel('Time (s)');
ylabel('Error (m/s)');
title('Vx Estimation Error');

figure;
plot(t,Error_Vy,'LineWidth',2);
grid on;
xlabel('Time (s)');
ylabel('Error (m/s)');
title('Vy Estimation Error');

%mean error analysis
Mean_Vx_Error = mean(Error_Vx);
Mean_Vy_Error = mean(Error_Vy);

fprintf('\n');
fprintf('Mean Vx Error = %.4f m/s\n',Mean_Vx_Error);
fprintf('Mean Vy Error = %.4f m/s\n',Mean_Vy_Error);

%Maximum error analysis
Max_Vx_Error = max(abs(Error_Vx));
Max_Vy_Error = max(abs(Error_Vy));

fprintf('Max Vx Error = %.4f m/s\n',Max_Vx_Error);
fprintf('Max Vy Error = %.4f m/s\n',Max_Vy_Error);


