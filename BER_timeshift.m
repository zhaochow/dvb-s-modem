close all; clear; clc;

% Parameters
nb_symb = 5e5;
nbps = 2; % Number of bits per symbol
constel = 'qam';
nb_bits = nb_symb*nbps;
bits_tx = randi(2,nb_bits,1)-1;
M = 50;
f_symb = 1e6; % symbols freq
fs = M*f_symb; % sampling freq
beta = 0.3;
T_symb = 1/f_symb;
nb_taps = 20*M+1;

t_shift = 1:5;
ratio_EbN0_db = -6:1:15;
ber = zeros(length(ratio_EbN0_db),length(t_shift));

%% TX
% Generate signal in baseband
symb_tx = mapping(bits_tx,nbps,constel);

% Create the halfroot Nyquist filter
[h_RRC,t] = half_Nyquist_filt(beta,T_symb,fs,nb_taps); % impulse response

% Convolution of the signal with the impulse response of the filter
signal_upSamp = upsample(symb_tx,M); % Up sampling by M
signal_tx = conv(signal_upSamp,h_RRC);

parfor j = 1:length(t_shift)
    res = zeros(length(ratio_EbN0_db),1);
    for i = 1:length(ratio_EbN0_db)
        [noise,N0] = Gaussian_noise(signal_tx,fs,ratio_EbN0_db(i),nb_bits);

%         signal_rx = signal_tx + real(noise); % For PAM
        signal_rx = signal_tx + noise;

        %% RX
        % Convolution of the signal with the impulse response of the filter
        signal_rx = conv(signal_rx,h_RRC);
        signal_rx_shift = signal_rx(nb_taps+t_shift(j):end-(nb_taps-1)+t_shift(j)); % Shifted samples
        signal_rx_shift = downsample(signal_rx_shift,M); % Down sample by M
        
        bits_rx_shift = demapping(signal_rx_shift,nbps,constel);
        err_shift = sum(abs(bits_tx - bits_rx_shift));
        res(i) = sum(abs(err_shift))/nb_bits;
    end
    ber(:,j) = res;
end

%%
figure
for j = 1:length(t_shift)
    semilogy(ratio_EbN0_db,ber(:,j))
    hold on
end
% legend('PAM 1 bit','PAM 2 bits','QAM 4 bits','QAM 6 bits')
grid on

% semilogy(bers0_qam4.data{1,1},bers0_qam4.data{1,2},'--')

% load('ber_pam2.mat')
% semilogy(bers0_pam2.data{1,1},bers0_pam2.data{1,2},'--')
% load('ber_pam4.mat')
% semilogy(bers0_pam4.data{1,1},bers0_pam4.data{1,2},'--')
% load('ber_qam16.mat')
% semilogy(bers0_qam16.data{1,1},bers0_qam16.data{1,2},'--')
% load('ber_qam64.mat')
% semilogy(bers0_qam64.data{1,1},bers0_qam64.data{1,2},'--')
ylim([1e-7 1])
% legend('0.02T','0.04T','0.06T','0.08T','0.1T','No time offset')