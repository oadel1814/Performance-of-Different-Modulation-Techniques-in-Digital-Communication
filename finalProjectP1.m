%% Comparison of OOK, PRK, and BFSK Modulation Schemes
%  BER vs SNR performance (simulation only)
clear; clc; close all;

%% ---------- Simulation Parameters ----------
N         = 1e6;             % bits per SNR value
SNR_dB    = 0:3:60;          % SNR range in dB
SNR_lin   = 10.^(SNR_dB/10); % SNR (linear)

BER_OOK  = [];
BER_PRK  = [];
BER_BFSK = [];

%% ---------- Random Binary Data ----------
bits = randi([0 1], 1, N);   % same bits used for all schemes (fair comparison)

%% ---------- Loop over SNR ----------
for k = 1:length(SNR_dB)
    SNR = SNR_lin(k);

    % ===== OOK : 0 -> 0, 1 -> 1 =====
    s_OOK   = bits;
    Ps_OOK  = mean(abs(s_OOK).^2);     % signal power = 0.5
    Pn_OOK  = Ps_OOK / SNR;            % noise power from SNR definition
    sigma   = sqrt(Pn_OOK);            % real AWGN std
    noise   = sigma * randn(1, N);
    r_OOK   = s_OOK + noise;
    det_OOK = r_OOK > 0.5;             % optimum threshold = 0.5
    BER_OOK = [BER_OOK sum(det_OOK ~= bits)/N];

    % ===== PRK : 0 -> -1, 1 -> +1 =====
    s_PRK   = 2*bits - 1;
    Ps_PRK  = mean(abs(s_PRK).^2);     % signal power = 1
    Pn_PRK  = Ps_PRK / SNR;
    sigma   = sqrt(Pn_PRK);
    noise   = sigma * randn(1, N);
    r_PRK   = s_PRK + noise;
    det_PRK = r_PRK > 0;               % threshold = 0
    BER_PRK = [BER_PRK sum(det_PRK ~= bits)/N];

    % ===== BFSK : 0 -> 1 (in-phase), 1 -> j (quadrature) =====
    s_BFSK = zeros(1, N);
    s_BFSK(bits == 0) = 1;
    s_BFSK(bits == 1) = 1i;
    Ps_BFSK = mean(abs(s_BFSK).^2);    % signal power = 1
    Pn_BFSK = Ps_BFSK / SNR;
    sigma   = sqrt(Pn_BFSK/2);         % per-dimension std (complex noise)
    noise   = sigma * (randn(1, N) + 1i*randn(1, N));
    r_BFSK  = s_BFSK + noise;
    det_BFSK = imag(r_BFSK) > real(r_BFSK);   % decide closer orthogonal axis
    BER_BFSK = [BER_BFSK sum(det_BFSK ~= bits)/N];
end

%% ---------- Plot ----------
figure('Color','w');
semilogy(SNR_dB, BER_OOK,  'ro-', 'LineWidth', 1.6, 'MarkerSize', 7); hold on;
semilogy(SNR_dB, BER_PRK,  'bs-', 'LineWidth', 1.6, 'MarkerSize', 7);
semilogy(SNR_dB, BER_BFSK, 'g^-', 'LineWidth', 1.6, 'MarkerSize', 7);

grid on; box on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for OOK, PRK, and BFSK');
legend({'OOK','PRK','BFSK'}, 'Location','southwest');
ylim([1e-7 1]);
xlim([SNR_dB(1) SNR_dB(end)]);