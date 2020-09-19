function [bits_rx, symb_rx, err, N0] = channel_txrx(packet,ratio_EbN0_db)
    % Parameters
    nbps = 1; % Number of bits per symbol
    nb_bits = length(packet);
    nb_symb = nb_bits/nbps;
    bits_tx = packet;
    M = 4;
    f_symb = 1e6; % symbols freq
    fs = M*f_symb; % sampling freq
    beta = 0.3;
    T_symb = 1/f_symb;
    nb_taps = 33;

    %% TX
    % Generate signal in baseband
    symb_tx = mapping(bits_tx,nbps,'pam');

    % Create the halfroot Nyquist filter
    [h_RRC,t] = half_Nyquist_filt(beta,T_symb,fs,nb_taps); % impulse response

    % Convolution of the signal with the impulse response of the filter
    signal_upSamp = upsample(symb_tx,M); % Up sampling by M
    signal_tx = conv(signal_upSamp,h_RRC);

    [noise,N0] = Gaussian_noise(signal_tx,fs,ratio_EbN0_db,nb_bits);

    signal_rx = signal_tx + real(noise);
%     signal_rx = signal_tx + noise;

    %% RX
    % Convolution of the signal with the impulse response of the filter
    symb_rx = conv(signal_rx,h_RRC);
    symb_rx = symb_rx(nb_taps:end-(nb_taps-1));
    symb_rx = downsample(symb_rx,M); % Down sample by M

    bits_rx = demapping(symb_rx,nbps,'pam');
    % Check result
    err = sum(abs(bits_tx - bits_rx));
end

