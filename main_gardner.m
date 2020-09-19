close all; clear; clc;

% Parameters
nb_symb = 1000;
nbps = 1; % Number of bits per symbol
constel = 'pam';
nb_bits = nb_symb*nbps;

M = 100;
f_symb = 1e6; % symbols freq
fs = M*f_symb; % sampling freq
beta = 0.3;
T_symb = 1/f_symb;
nb_taps = 20*M+1;

deltaF = [0, 2e9*1e-5, 2e9*5e-5];
t_shift = 30;
ratio_EbN0_db = 10;
err_weight = [0.01, 0.01, 0.01];
time_err = zeros(nb_symb,length(err_weight));
err_var = zeros(nb_symb,length(err_weight));
iter_max = 100;
%% TX
% Generate signal in baseband

for i = 1:length(err_weight)
    iter_err = zeros(nb_symb,iter_max);
    for j = 1:iter_max
    
    bits_tx = randi(2,nb_bits,1)-1;
    symb_tx = mapping(bits_tx,nbps,constel);

    % Create the halfroot Nyquist filter
    [h_RRC,t] = half_Nyquist_filt(beta,T_symb,fs,nb_taps); % impulse response
    % h_RRC = rcosdesign(beta,8,4,'sqrt');

    % Convolution of the signal with the impulse response of the filter
    % signal_upSamp = repelem(symb_tx,M); % Up sampling by M
    signal_upSamp = upsample(symb_tx,M); % Up sampling by M
    signal_tx = conv(signal_upSamp,h_RRC);

%     [noise,N0] = Gaussian_noise(signal_tx,fs,ratio_EbN0_db,nb_bits);

%     signal_rx = signal_tx + real(noise); % For PAM
%     signal_rx = signal_tx + noise;
    
    %% RX
    
    signal_rx = signal_tx;
    ts = ((0:length(signal_rx)-1)/fs)';
    signal_rx_CFO = signal_rx.*exp(1j*(2*pi*deltaF(i)*ts));   
    
    % Convolution of the signal with the impulse response of the filter
    signal_rx_CFO = conv(signal_rx_CFO,h_RRC);
    signal_rx_real = signal_rx_CFO(nb_taps:end-(nb_taps-1));
    idx_real = (1:M:length(signal_rx_real))';
    signal_rx_shift = signal_rx_CFO(nb_taps+t_shift:end-(nb_taps-1)+t_shift);
%     shift_err = shift_err + K * real(signal_rx(n-0.5)*(conj(signal_rx(n)) - ...
%         conj(signal_rx(n-1))));
%     a = downsample(signal_rx_real,M); % Down sample by M
    [signal_rx_corr,idx_corr,errcorr] = gardnercorr(signal_rx_shift,fs,f_symb,nb_symb,err_weight(i));
    idx_corr = idx_corr + t_shift;
    
    idx_err = idx_corr - idx_real;
    iter_err(:,j) = idx_err*1/fs;
    
%     bits_rx = demapping(a,nbps,constel);
%     err = sum(abs(bits_tx - bits_rx));
%     res(i) = sum(abs(err))/nb_bits;
    end
    err_var(:,i) = var(iter_err,0, 2);
    time_err(:,i) = mean(iter_err,2);
end

%%
colors = ["b","r","k","m","c","y","w"];
plots = [];
figure
for i = 1:length(err_weight)
    subset = 1:20:nb_symb;
    plot(subset,time_err(subset,i) + sqrt(err_var(subset,i)), colors(i)+"--");
    hold on
    h = plot(subset,time_err(subset,i), colors(i));
    plots(i) = h;
    hold on
    plot(subset,time_err(subset,i) - sqrt(err_var(subset,i)), colors(i)+"--");
end
legend(plots,'No CFO', '10 ppm', '50 ppm')
grid on

% windowSize = 10; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% y = filter(b,a,time_err);
% 
% figure
% for i = 1:length(err_weight)
%     plot([time_err(1:windowSize-1,i) ; y(windowSize:end,i)])
%     hold on
% end
% grid on
