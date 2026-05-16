% ================ PART III =================
%% Compare the Performance of BPSK, QPSK, 4-QAM, and 4-ASK
%  Noise added in signal-space domain using proper Eb/N0 normalization

clear; clc; close all;

%% ---------- Simulation Parameters ----------
N   = 1e6;          % number of bits
snr = 0:3:60;       % Eb/N0 range in dB

Rb  = 1000;         % bit rate
Tb  = 1/Rb;         % bit period

%% ---------- Common Bit Energy ----------
Eb = Tb/2;

%% ---------- Common Amplitude ----------
% Choose amplitude so BPSK symbol energy = Eb
A = sqrt(Eb);

%% ==========================================================
%% BPSK
%% Constellation: {-A, +A}
%% ==========================================================
fprintf('Simulating BPSK...\n');

bits_BPSK = randi([0 1], 1, N);

% Mapping
tx_BPSK = (2*bits_BPSK - 1) * A;

BER_BPSK = zeros(size(snr));

for i = 1:length(snr)

    SNR_lin = 10^(snr(i)/10);

    % AWGN std
    sigma = sqrt(Eb/(2*SNR_lin));

    % Channel
    rx = tx_BPSK + sigma*randn(1,N);

    % Detection
    rx_bits = rx >= 0;

    % BER
    BER_BPSK(i) = sum(bits_BPSK ~= rx_bits)/N;
end

%% ==========================================================
%% QPSK (Gray coded)
%% Constellation: (±A, ±A)
%% ==========================================================
fprintf('Simulating QPSK...\n');

bits_QPSK = randi([0 1], 1, N);

num_sym_QPSK = N/2;

bits_pairs = reshape(bits_QPSK,2,num_sym_QPSK);

% Gray-coded mapping
I_QPSK =  (2*bits_pairs(1,:) - 1) * A;
Q_QPSK = -(2*bits_pairs(2,:) - 1) * A;

BER_QPSK = zeros(size(snr));

for i = 1:length(snr)

    SNR_lin = 10^(snr(i)/10);

    sigma = sqrt(Eb/(2*SNR_lin));

    % AWGN channel
    I_rx = I_QPSK + sigma*randn(1,num_sym_QPSK);
    Q_rx = Q_QPSK + sigma*randn(1,num_sym_QPSK);

    % Detection
    rx_b1 = I_rx >= 0;
    rx_b2 = Q_rx <= 0;

    rx_bits = reshape([rx_b1; rx_b2],1,[]);

    % BER
    BER_QPSK(i) = sum(bits_QPSK ~= rx_bits)/N;
end

%% ==========================================================
%% 4-QAM
%% Standard rectangular 4-QAM
%% ==========================================================
fprintf('Simulating 4-QAM...\n');

bits_4QAM = randi([0 1],1,N);

num_sym_4QAM = N/2;

bits_pairs4 = reshape(bits_4QAM,2,num_sym_4QAM);

I_4QAM = (2*bits_pairs4(1,:) - 1) * A;
Q_4QAM = (2*bits_pairs4(2,:) - 1) * A;

BER_4QAM = zeros(size(snr));

for i = 1:length(snr)

    SNR_lin = 10^(snr(i)/10);

    sigma = sqrt(Eb/(2*SNR_lin));

    % AWGN channel
    I_rx4 = I_4QAM + sigma*randn(1,num_sym_4QAM);
    Q_rx4 = Q_4QAM + sigma*randn(1,num_sym_4QAM);

    % Detection
    rx_b1 = I_rx4 >= 0;
    rx_b2 = Q_rx4 >= 0;

    rx_bits4 = reshape([rx_b1; rx_b2],1,[]);

    % BER
    BER_4QAM(i) = sum(bits_4QAM ~= rx_bits4)/N;
end

%% ==========================================================
%% 4-ASK (Gray coded)
%% Levels: {-3A, -A, +A, +3A}
%%
%% Mapping:
%% 00 -> -3A
%% 01 -> -A
%% 11 -> +A
%% 10 -> +3A
%% ==========================================================
fprintf('Simulating 4-ASK...\n');

bits_4ASK = randi([0 1],1,N);

num_sym_4ASK = N/2;

bits_pairs_4ASK = reshape(bits_4ASK,2,num_sym_4ASK);

sym_idx = bits_pairs_4ASK(1,:)*2 + bits_pairs_4ASK(2,:);

tx_4ASK = zeros(1,num_sym_4ASK);

tx_4ASK(sym_idx == 0) = -3*A;   % 00
tx_4ASK(sym_idx == 1) = -1*A;   % 01
tx_4ASK(sym_idx == 3) = +1*A;   % 11
tx_4ASK(sym_idx == 2) = +3*A;   % 10

% Symbol energy
Es_4ASK = mean(tx_4ASK.^2);

% Bit energy
Eb_4ASK = Es_4ASK / 2;

BER_4ASK = zeros(size(snr));

for i = 1:length(snr)

    SNR_lin = 10^(snr(i)/10);

    sigma = sqrt(Eb_4ASK/(2*SNR_lin));

    % AWGN channel
    rx_4ASK = tx_4ASK + sigma*randn(1,num_sym_4ASK);

    rx_bits_4ASK = zeros(2,num_sym_4ASK);

    % Decision regions
    idx = rx_4ASK > 2*A;
    rx_bits_4ASK(:,idx) = repmat([1;0],1,sum(idx));

    idx = (rx_4ASK > 0) & (rx_4ASK <= 2*A);
    rx_bits_4ASK(:,idx) = repmat([1;1],1,sum(idx));

    idx = (rx_4ASK <= 0) & (rx_4ASK > -2*A);
    rx_bits_4ASK(:,idx) = repmat([0;1],1,sum(idx));

    idx = rx_4ASK <= -2*A;
    rx_bits_4ASK(:,idx) = repmat([0;0],1,sum(idx));

    rx_bits4ASK = reshape(rx_bits_4ASK,1,[]);

    % BER
    BER_4ASK(i) = sum(bits_4ASK ~= rx_bits4ASK)/N;
end

%% ==========================================================
%% Plot Results
%% ==========================================================
figure('Color','w');

semilogy(snr, BER_BPSK, 'b-o', ...
    'LineWidth',2,'MarkerSize',6);
hold on;

semilogy(snr, BER_QPSK, 'r-s', ...
    'LineWidth',2,'MarkerSize',6);

semilogy(snr, BER_4QAM, 'm-^', ...
    'LineWidth',2,'MarkerSize',6);

semilogy(snr, BER_4ASK, 'k--d', ...
    'LineWidth',2,'MarkerSize',6);

grid on;

xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');

title('BER vs E_b/N_0 for BPSK, QPSK, 4-QAM, and 4-ASK');

legend('BPSK', ...
       'QPSK', ...
       '4-QAM', ...
       '4-ASK', ...
       'Location','southwest');

xlim([snr(1) snr(end)]);
ylim([1e-7 1]);

fprintf('Done!\n');