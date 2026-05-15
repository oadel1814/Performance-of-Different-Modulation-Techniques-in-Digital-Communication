% ================PART I=================
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