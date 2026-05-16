% ================PART I=================
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

% ================PART II=================
% Constants
fs = 20000;
fc = 1000;
Rb = 1000;
Tb = 1/Rb;
% Input random bits --> 10^6 bits
bits = 1e6;
% Vector of random bits
vector = randi([0,1],1,bits);
% SNR in dB, incrementing by 3 dB from 0 to 60 dB
snr = 0:3:60;
% BASK
% Initialize BER vector to store results for each SNR
BER_BASK = zeros(1,length(snr));
% Samples per bit to represent each symbol
samples_per_bit = fs/Rb;
% Define orthonormal basis for ASK
t_bit = 0:1/fs:Tb-1/fs;
% Carrier Signal Generation
E_bit = Tb / 2;
A = sqrt(E_bit)/2;
% Symbol mask 0 --> -A and 1 --> A
symbol_mask = (2*vector - 1) * A;
phi_I_BASK = sqrt(2/Tb) * cos(2*pi*fc*t_bit);
% Transmission:
% Multiply the basis by the transmission vector
% symbol_mask is 1*1000000 and phi_I_BASK is 1*20
y_matrix_form = phi_I_BASK' .* symbol_mask;
% Reshape the transmission into a continuous frame
y = reshape(y_matrix_form, 1, []);

decision_boundary = 0;

for i = 1:length(snr)
    % Add noise to the input vector
    Rx_sequence = awgn(y, snr(i), 'measured');
    % Demodulate the transmitted waveform
    Rx_matrix = reshape(Rx_sequence, samples_per_bit, bits);
    % Multiply by the basis function and integrate
    I_received = sum(Rx_matrix .* phi_I_BASK') * (1/fs);
    noisy_vector = I_received >= decision_boundary;
    % Calculate the Bit Error Rate (BER) using the XOR operator
    errors = sum(xor(vector, noisy_vector));
    BER_BASK(i) = errors / bits;
end

% 4-ASK
% Initialize BER vector to store results for each SNR
BER_4ASK = zeros(1,length(snr));
% Signal transmitted --> 2 bits per symbol
Ts = 2 * Tb;
samples_per_symbol = Ts * fs;
% Group bits into pairs (as 4-ASK groups every two symbols and transmits
% them)
bits_grouped = reshape(vector, 2, bits/2);
% Carrier Signal Generation
E_bit = Tb / 2;
A = sqrt(E_bit)/2;
% Gray coding --> 0 then 1 then 3 then 2
symbol_index = bits_grouped(1,:)*2 + bits_grouped(2,:);
symbol_mask = zeros(1,bits/2);

% 00
symbol_mask(symbol_index == 0) = -3 * A;
% 01
symbol_mask(symbol_index == 1) = -1 * A;
% 11
symbol_mask(symbol_index == 3) = 1 * A;
% 10
symbol_mask(symbol_index == 2) = 3 * A;

% Define orthonormal basis for 4ASK
t_sym = 0:1/fs:Ts-1/fs;
phi_I_4ASK = sqrt(2/Ts) * cos(2*pi*fc*t_sym);

% Modulate the signal:
y_matrix_form = phi_I_4ASK' .* symbol_mask;
% Reshape the transmission into a continuous frame
y = reshape(y_matrix_form, 1, []);

for i = 1:length(snr)
    % Add noise to the input vector
    Rx_sequence = awgn(y, snr(i), 'measured');

    % Demodulate the transmitted signal:
    Rx_matrix = reshape(Rx_sequence, samples_per_symbol, bits/2);
    I_received = sum(Rx_matrix .* phi_I_4ASK') * (1/fs);

    % Remapping back to bits:
    rx_bits = zeros(2,bits/2);
    % Dividing the space into regions comprised of different decision
    % boundaries

    % Region 1 --> > 2A ==> Symbol would be 3A and the Bits are 10
    idx = I_received > 2*A;
    rx_bits(1,idx) = 1; rx_bits(2,idx) = 0;

    % Region 2 --> 0 to 2A ==> Symbol would be 1A and the Bits are 11
    idx = (I_received > 0) & (I_received <= 2*A);
    rx_bits(1,idx) = 1; rx_bits(2,idx) = 1;

    % Region 3 --> 0 down to -2A ==> Symbol would be -1A and the Bits are
    % 01
    idx = (I_received <= 0) & (I_received > -2*A);
    rx_bits(1,idx) = 0; rx_bits(2,idx) = 1;

    % Region 4 --> <= -2A ==> Symbol would be -3A and the Bits are 00
    idx = I_received <= -2*A;
    rx_bits(1,idx) = 0; rx_bits(2,idx) = 0;

    % Reshape received bits to 1D vector before calculating the BER
    rx_vector = reshape(rx_bits, 1, []);
    % Calculate the Bit Error Rate (BER) using the XOR operator
    errors = sum(xor(vector, rx_vector));
    BER_4ASK(i) = errors / bits;
end


% 8-ASK
% Initialize BER vector to store results for each SNR
BER_8ASK = zeros(1,length(snr));
% Number of bits will be changed in order to be able to divide by 3 without
% throwing any errors:
bits = 999996;
% Update vector of random bits
vector = randi([0,1],1,bits);
% Signal transmitted --> 2 bits per symbol
Ts = 3 * Tb;
samples_per_symbol = Ts * fs;
% Group bits into threes (as 8-ASK groups every three symbols and transmits
% them)
bits_grouped = reshape(vector, 3, bits/3);
% Carrier Signal Generation
E_bit = Tb / 2;
A = sqrt(E_bit)/2;
% Gray coding --> 0 --> 1 --> 3 --> 2 --> 6 --> 7 --> 5 --> 4
symbol_index = bits_grouped(1,:)*4 + bits_grouped(2,:)*2 + bits_grouped(3,:);
symbol_mask = zeros(1,bits/3);

% 000
symbol_mask(symbol_index == 0) = -7 * A;
% 001
symbol_mask(symbol_index == 1) = -5 * A;
% 011
symbol_mask(symbol_index == 3) = -3 * A;
% 010
symbol_mask(symbol_index == 2) = -1 * A;
% 110
symbol_mask(symbol_index == 6) = 1 * A;
% 111
symbol_mask(symbol_index == 7) = 3 * A;
% 101
symbol_mask(symbol_index == 5) = 5 * A;
% 100
symbol_mask(symbol_index == 4) = 7 * A;

% Define orthonormal basis for 8ASK
t_sym = 0:1/fs:Ts-1/fs;
phi_I_8ASK = sqrt(2/Ts) * cos(2*pi*fc*t_sym);

% Modulate the signal:
y_matrix_form = phi_I_8ASK' .* symbol_mask;
% Reshape the transmission into a continuous frame
y = reshape(y_matrix_form, 1, []);

for i = 1:length(snr)
    % Add noise to the input vector
    Rx_sequence = awgn(y, snr(i), 'measured');

    % Demodulate the transmitted signal:
    Rx_matrix = reshape(Rx_sequence, samples_per_symbol, bits/3);
    I_received = sum(Rx_matrix .* phi_I_8ASK') * (1/fs);

    % Remapping back to bits:
    rx_bits = zeros(3,bits/3);
    % Dividing the space into regions comprised of different decision
    % boundaries

    % Region 1 --> > 6A ==> Symbol would be 7A and the Bits are 100
    idx = I_received > 6*A;
    rx_bits(1,idx) = 1; rx_bits(2,idx) = 0; rx_bits(3,idx) = 0;

    % Region 2 --> 4A to 6A ==> Symbol would be 5A and the Bits are 101
    idx = (I_received > 4*A) & (I_received <= 6*A);
    rx_bits(1,idx) = 1; rx_bits(2,idx) = 0; rx_bits(3,idx) = 1;

    % Region 3 --> 2A to 4A ==> Symbol would be 3A and the Bits are 111
    idx = (I_received > 2*A) & (I_received <= 4*A);
    rx_bits(1,idx) = 1; rx_bits(2,idx) = 1; rx_bits(3,idx) = 1;

    % Region 4 --> 0 to 2A ==> Symbol would be 1A and the Bits are 110
    idx = (I_received > 0) & (I_received <= 2*A);
    rx_bits(1,idx) = 1; rx_bits(2,idx) = 1; rx_bits(3,idx) = 0;

    % Region 5 --> -2A to 0 ==> Symbol would be -1A and the Bits are 010
    idx = (I_received > -2*A) & (I_received <= 0);
    rx_bits(1,idx) = 0; rx_bits(2,idx) = 1; rx_bits(3,idx) = 0;

    % Region 6 --> -4A to -2A ==> Symbol would be -3A and the Bits are 011
    idx = (I_received > -4*A) & (I_received <= -2*A);
    rx_bits(1,idx) = 0; rx_bits(2,idx) = 1; rx_bits(3,idx) = 1;

    % Region 7 --> -6A to -4A ==> Symbol would be -5A and the Bits are 001
    idx = (I_received > -6*A) & (I_received <= -4*A);
    rx_bits(1,idx) = 0; rx_bits(2,idx) = 0; rx_bits(3,idx) = 1;

    % Region 7 --> -6A <= ==> Symbol would be -7A and the Bits are 000
    idx = (I_received <= -6*A);
    rx_bits(1,idx) = 0; rx_bits(2,idx) = 0; rx_bits(3,idx) = 0;

    % Reshape received bits to 1D vector before calculating the BER
    rx_vector = reshape(rx_bits, 1, []);
    % Calculate the Bit Error Rate (BER) using the XOR operator
    errors = sum(xor(vector, rx_vector));
    BER_8ASK(i) = errors / bits;
end

% Plot BER vs SNR
figure;
semilogy(snr, BER_8ASK, 'g-^', 'LineWidth', 2);
hold on;
semilogy(snr, BER_4ASK, 'b-o', 'LineWidth', 2);
hold on;
semilogy(snr, BER_BASK, 'r-s', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR -- 2, 4, 8 ASK Modulation');
legend('8ASK', '4ASK', 'BASK (2ASK)');

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