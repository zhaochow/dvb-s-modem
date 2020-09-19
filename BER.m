close all; clear; clc;

% Parameters
nb_symb = 5000;
ratio_EbN0_db = -6:0.1:25;
ber = zeros(length(ratio_EbN0_db),4);
for j = 1:4
    if j <= 2
        constel_name = 'pam';
        constel_nb = j;
    else
        constel_name = 'qam';
        constel_nb = 2*(j-1);
    end
nbps = constel_nb; % Number of bits per symbol
nb_bits = nb_symb*nbps;
bits_tx = randi(2,nb_bits,1)-1;
M = 4;
f_symb = 1e6; % symbols freq
fs = M*f_symb; % sampling freq
ts = (0:(nb_bits-1))/fs;
beta = 0.3;
T_symb = 1/f_symb;
nb_taps = 33;

%% TX
% Generate signal in baseband
symb_tx = mapping(bits_tx,nbps,constel_name);

% Create the halfroot Nyquist filter
[h_RRC,t] = half_Nyquist_filt(beta,T_symb,fs,nb_taps); % impulse response

% Convolution of the signal with the impulse response of the filter
signal_upSamp = upsample(symb_tx,M); % Up sampling by M
signal_tx = conv(signal_upSamp,h_RRC);

res = zeros(length(ratio_EbN0_db),1);
for i = 1:length(ratio_EbN0_db)
    [noise,N0] = Gaussian_noise(signal_tx,fs,ratio_EbN0_db(i),nb_bits);

    if j <= 2
        signal_rx = signal_tx + real(noise);
    else
        signal_rx = signal_tx + noise;
    end

%% RX
    % Convolution of the signal with the impulse response of the filter
    symb_rx = conv(signal_rx,h_RRC);
    symb_rx = symb_rx(nb_taps:end-(nb_taps-1));
    symb_rx = downsample(symb_rx,M); % Down sample by M
    
    bits_rx = demapping(symb_rx,nbps,constel_name);
    % Check result
    err = bits_tx - bits_rx;
    res(i) = sum(abs(err))/nb_bits;
end
ber(:,j) = res;
end

%% PLOT
figure
for i = 1:4
semilogy(ratio_EbN0_db,ber(:,i))
hold on
end
legend('2-PAM','4-PAM','16-QAM','64-QAM')
grid on

load('ber_pam2.mat')
semilogy(bers0_pam2.data{1,1},bers0_pam2.data{1,2},'--')
load('ber_pam4.mat')
semilogy(bers0_pam4.data{1,1},bers0_pam4.data{1,2},'--')
load('ber_qam16.mat')
semilogy(bers0_qam16.data{1,1},bers0_qam16.data{1,2},'--')
load('ber_qam64.mat')
semilogy(bers0_qam64.data{1,1},bers0_qam64.data{1,2},'--')
% legend('2-PAM (ideal)','4-PAM (ideal)','16-QAM (ideal)','64-QAM (ideal)')
ylim([1e-5 1])