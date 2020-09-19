close all; clear; clc;

% Parameters
nb_symb = 512;
nbps = 2; % Number of bits per symbol
constel = 'qam';
nb_data = nb_symb*nbps;
M = 50;
f_symb = 1e6; % symbols freq
fs = M*f_symb; % sampling freq
beta = 0.3;
T_symb = 1/f_symb;
nb_taps = 20*M+1;

% ---For CFO---------
fc = 2e9;
deltaF = fc*5e-5;
% deltaF = 0;
phi = 0;
% phi = 0.02:0.02:0.2;
% --------------------
t_shift = round(0.02*M);
% t_shift = 0;
nb_pilots = [40];
K = [8];
nb_it = 200;
ratio_EbN0_db = 0:2:4;
acqu_err = zeros(length(ratio_EbN0_db),length(nb_pilots));

for m = 1:length(K)
for j = 1:length(nb_pilots)
    nb_bits = nb_data + nb_pilots(j)*nbps;
    for i = 1:length(ratio_EbN0_db)
        tmp = zeros(1,nb_it);
        for k = 1:nb_it                % For computing the standard deviation
            pilot_bits = randi(2,nb_pilots(j)*nbps,1)-1;
            bits_tx = [randi(2,nb_data/2,1)-1; pilot_bits; randi(2,nb_data/2,1)-1];

            %% TX
            % Generate signal in baseband
            symb_tx = mapping(bits_tx,nbps,constel);
            pilot = symb_tx(nb_symb/2+1:nb_symb/2+nb_pilots(j));

            % Create the halfroot Nyquist filter
            [h_RRC,t] = half_Nyquist_filt(beta,T_symb,fs,nb_taps); % impulse response

            % Convolution of the signal with the impulse response of the filter
            signal_upSamp = upsample(symb_tx,M); % Up sampling by M
            signal_tx = conv(signal_upSamp,h_RRC);

            [noise,N0] = Gaussian_noise(signal_tx,fs,ratio_EbN0_db(i),nb_bits);

            if(strcmp(constel,'pam'))
                signal_rx = signal_tx + real(noise);
            else
                signal_rx = signal_tx + noise;
            end

            %% RX
            % Convolution of the signal with the impulse response of the filter
            
            ts = ((0:length(signal_rx)-1)/fs)';
            signal_rx_CFO = signal_rx.*exp(1j*(2*pi*deltaF*ts+phi));
%             signal_rx_CFO = signal_rx;
            
            signal_rx_CFO = conv(signal_rx_CFO,h_RRC);
            signal_rx_shift = signal_rx_CFO(nb_taps+t_shift:end-(nb_taps-1)+t_shift);
            signal_rx_shift = downsample(signal_rx_shift,M); % Down sample by M
            
            
            Dk = diffcorr(signal_rx_shift, pilot, K(m));
            n_est = get_pilot_pos(Dk);
            tmp(k) = n_est-(nb_symb/2);
%             delta_est = get_phase_diff(Dk, n_est,T_symb);
%             tmp(k) = delta_est - deltaF;
            
            if(strcmp(constel,'pam'))
                signal_rx_CFO = real(signal_rx_CFO); % For PAM
            end
            bits_rx_CFO = demapping(signal_rx_CFO,nbps,constel);
%             err_CFO = sum(abs(bits_tx - bits_rx_CFO));
%             res = sum(abs(err_CFO))/nb_bits;
        end
        acqu_err(i,j,m) = sqrt(var(tmp));
    end
end
end

% figure
% for j = 1:length(nb_pilots)
%     plot(ratio_EbN0_db,(acqu_err(:,j,1)))
%     ylim([0, 200]);
%     hold on
% end
% grid on

% figure
% for j = 1:length(K)
%     plot(ratio_EbN0_db,acqu_err(:,1,j))
%     ylim([0, 1]);
%     hold on
% end
% grid on

figure
for j = 1:length(K)
    plot(ratio_EbN0_db,10e5*(acqu_err(:,1,j)/fc))
    ylim([0, 1]);
    hold on
end
grid on

% legend('PAM 1 bit','PAM 2 bits','QAM 4 bits','QAM 6 bits')


% load('ber_pam2.mat')
% semilogy(bers0_pam2.data{1,1},bers0_pam2.data{1,2},'--')
% load('ber_pam4.mat')
% semilogy(bers0_pam4.data{1,1},bers0_pam4.data{1,2},'--')
% load('ber_qam16.mat')
% semilogy(bers0_qam16.data{1,1},bers0_qam16.data{1,2},'--')
% load('ber_qam64.mat')
% semilogy(bers0_qam64.data{1,1},bers0_qam64.data{1,2},'--')
% ylim([1e-5 1])