close all; clear; clc;

% Parameters
nb_symb = 1e5;
nbps = 4; % Number of bits per symbol
constel = 'qam';
nb_bits = nb_symb*nbps;
bits_tx = randi(2,nb_bits,1)-1;
M = 4;
f_symb = 1e6; % symbols freq
fs = M*f_symb; % sampling freq
beta = 0.3;
T_symb = 1/f_symb;
nb_taps = 20*M+1;

fc = 2e9;
deltaF = [2 10 30 40 50 80 100]*fc*1e-6;
phi = 0;
ratio_EbN0_db = -6:1:25;
ber = zeros(length(ratio_EbN0_db),length(deltaF));

%% TX
% Generate signal in baseband
symb_tx = mapping(bits_tx,nbps,constel);

% Create the halfroot Nyquist filter
[h_RRC,t_filt] = half_Nyquist_filt(beta,T_symb,fs,nb_taps); % impulse response

% Convolution of the signal with the impulse response of the filter
signal_upSamp = upsample(symb_tx,M); % Up sampling by M
signal_tx = conv(signal_upSamp,h_RRC);

for j = 1:length(deltaF)
    res = zeros(length(ratio_EbN0_db),1);
    for i = 1:length(ratio_EbN0_db)
        [noise,N0] = Gaussian_noise(signal_tx,fs,ratio_EbN0_db(i),nb_bits);

%         signal_rx = signal_tx + abs(noise);
        signal_rx = signal_tx + noise;

        %% RX
        % Convolution of the signal with the impulse response of the filter
        ts = ((0:length(signal_rx)-1)/fs)';
        signal_rx_CFO = signal_rx.*exp(1j*(2*pi*deltaF(j)*ts+phi)); % Add CFO
%         signal_rx_CFO = signal_rx;
        signal_rx_CFO = conv(signal_rx_CFO,h_RRC);
        ts = ((-(nb_taps-1)/2:length(signal_rx_CFO)-(nb_taps-1)/2-1)/fs)';
        signal_rx_CFO = signal_rx_CFO.*exp(-1j*(2*pi*deltaF(j)*ts+phi)); % Remove phase drift
        signal_rx_CFO = signal_rx_CFO(nb_taps:end-(nb_taps-1));
        signal_rx_CFO = downsample(signal_rx_CFO,M); % Down sample by M
        
%         signal_rx_CFO = real(signal_rx_CFO); % For PAM
        bits_rx_CFO = demapping(signal_rx_CFO,nbps,constel);
        err_CFO = sum(abs(bits_tx - bits_rx_CFO));
        res(i) = err_CFO/nb_bits;
    end
    ber(:,j) = res;
end

%%
figure
for j = 1:length(deltaF)
    semilogy(ratio_EbN0_db,ber(:,j))
    hold on
end
% legend('PAM 1 bit','PAM 2 bits','QAM 4 bits','QAM 6 bits')
grid on

% load('ber_pam2.mat')
% semilogy(bers0_pam2.data{1,1},bers0_pam2.data{1,2},'--')
% load('ber_pam4.mat')
% semilogy(bers0_pam4.data{1,1},bers0_pam4.data{1,2},'--')
load('ber_qam16.mat')
semilogy(bers0_qam16.data{1,1},bers0_qam16.data{1,2},'--')
% load('ber_qam64.mat')
% semilogy(bers0_qam64.data{1,1},bers0_qam64.data{1,2},'--')
ylim([1e-5 1])