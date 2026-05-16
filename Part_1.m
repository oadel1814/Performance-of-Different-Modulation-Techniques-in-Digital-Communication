%% =========================================================
% Comparison of OOK, PRK, and BFSK Modulation Schemes
% Manual Simulation + MATLAB Built-in Functions
%% =========================================================
clear; clc; close all;

%% Parameters
N = 1e5;
SNR_dB = 0:3:30;
SNR_lin = 10.^(SNR_dB/10);

BER_OOK = zeros(size(SNR_dB));
BER_PRK = zeros(size(SNR_dB));
BER_BFSK = zeros(size(SNR_dB));

BER_OOK_builtin = zeros(size(SNR_dB));
BER_PRK_builtin = zeros(size(SNR_dB));
BER_BFSK_builtin = zeros(size(SNR_dB));

bits = randi([0 1],1,N);

%% Theoretical BER
BER_OOK_th = 0.5 * erfc(sqrt(SNR_lin/4));
BER_PRK_th = 0.5 * erfc(sqrt(SNR_lin));
BER_BFSK_th = 0.5 * erfc(sqrt(SNR_lin/2));

%% FSK Parameters
freqsep = 1;
nsamp = 8;

%% =========================================================
% SNR Loop
%% =========================================================
for k = 1:length(SNR_dB)
    SNR = SNR_lin(k);
    
    %% Manual OOK
    s_OOK = double(bits);
    Ps_OOK = mean(abs(s_OOK).^2);
    noise = sqrt(Ps_OOK/SNR) * randn(1,N);
    r_OOK = s_OOK + noise;
    det_OOK = r_OOK > 0.5;
    BER_OOK(k) = sum(det_OOK ~= bits)/N;
    
    %% Manual PRK
    s_PRK = 2*bits - 1;
    Ps_PRK = mean(abs(s_PRK).^2);
    noise = sqrt(Ps_PRK/SNR) * randn(1,N);
    r_PRK = s_PRK + noise;
    det_PRK = r_PRK > 0;
    BER_PRK(k) = sum(det_PRK ~= bits)/N;
    
    %% Manual BFSK
    s_BFSK = complex(double(bits==0), double(bits==1));
    Ps_BFSK = mean(abs(s_BFSK).^2);
    sigma = sqrt(Ps_BFSK/(2*SNR));
    noise = sigma * (randn(1,N) + 1i*randn(1,N));
    r_BFSK = s_BFSK + noise;
    det_BFSK = imag(r_BFSK) > real(r_BFSK);
    BER_BFSK(k) = sum(det_BFSK ~= bits)/N;
    
    %% Built-in OOK
    tx_OOK = pammod(bits,2);
    tx_OOK = (tx_OOK + 1)/2;
    Ps = mean(abs(tx_OOK).^2);
    noise = sqrt(Ps/SNR) * randn(1,N);
    rx_OOK = tx_OOK + noise;
    det_OOK_builtin = rx_OOK > 0.5;
    BER_OOK_builtin(k) = sum(det_OOK_builtin ~= bits)/N;
    
    %% Built-in PRK
    tx_PRK = pskmod(bits,2);
    Ps = mean(abs(tx_PRK).^2);
    noise = sqrt(Ps/SNR) * randn(1,N);
    rx_PRK = tx_PRK + noise;
    det_PRK_builtin = pskdemod(rx_PRK,2);
    BER_PRK_builtin(k) = sum(det_PRK_builtin ~= bits)/N;
    
    %% Built-in BFSK (Fixed)
    M = 2;
    tx_BFSK = fskmod(bits,M,freqsep,nsamp,nsamp);
    Ps = mean(abs(tx_BFSK).^2);
    SNR_sample = SNR/nsamp;
    sigma = sqrt(Ps / (2 * SNR_sample));
    noise = sigma * (randn(size(tx_BFSK)) + 1i*randn(size(tx_BFSK)));
    rx_BFSK = tx_BFSK + noise;
    det_BFSK_builtin = fskdemod(rx_BFSK,M,freqsep,nsamp,nsamp);
    BER_BFSK_builtin(k) = sum(det_BFSK_builtin ~= bits)/N;
end

%% =========================================================
% Figure 1 — Manual Simulation
%% =========================================================
figure('Name','Manual Simulation','NumberTitle','off');
semilogy(SNR_dB,BER_OOK,'ro-','LineWidth',1.8,'MarkerSize',6);
hold on;
semilogy(SNR_dB,BER_PRK,'bs-','LineWidth',1.8,'MarkerSize',6);
semilogy(SNR_dB,BER_BFSK,'g^-','LineWidth',1.8,'MarkerSize',6);
semilogy(SNR_dB,BER_OOK_th,'r--','LineWidth',1.2);
semilogy(SNR_dB,BER_PRK_th,'b--','LineWidth',1.2);
semilogy(SNR_dB,BER_BFSK_th,'g--','LineWidth',1.2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR — Manual Simulation');
legend('OOK sim','PRK sim','BFSK sim','OOK theory','PRK theory','BFSK theory','Location','southwest');
ylim([1e-6 1]);
xlim([0 30]);

%% =========================================================
% Figure 2 — MATLAB Built-in Functions
%% =========================================================
figure('Name','Built-in Functions','NumberTitle','off');
semilogy(SNR_dB,BER_OOK_builtin,'ro-','LineWidth',1.8,'MarkerSize',6);
hold on;
semilogy(SNR_dB,BER_PRK_builtin,'bs-','LineWidth',1.8,'MarkerSize',6);
semilogy(SNR_dB,BER_BFSK_builtin,'g^-','LineWidth',1.8,'MarkerSize',6);
semilogy(SNR_dB,BER_OOK_th,'r--','LineWidth',1.2);
semilogy(SNR_dB,BER_PRK_th,'b--','LineWidth',1.2);
semilogy(SNR_dB,BER_BFSK_th,'g--','LineWidth',1.2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR — MATLAB Built-in Functions vs Theory');
legend('OOK built-in','PRK built-in','BFSK built-in','OOK theory','PRK theory','BFSK theory','Location','southwest');
ylim([1e-6 1]);
xlim([0 30]);