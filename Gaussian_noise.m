function [noise,N0] = Gaussian_noise(s,fs,ratio_db,nb_bits)
% INPUTS: 
% - s: signal
% - fs: sampling frequency
% - ratio_db: ratio Eb/N0 in dB
% - nb_bits: number of bits
% 
% OUTPUTS: 
% - n: noise
% - N0: noise spectral density

% Creates n signal corresponding to the Gaussian noise

ratio = 10^(ratio_db/10);
e_b = 1/(length(s)*fs)*trapz(abs(s.^2)); % Energy per bit (upsampled)
e_b = e_b * length(s)/nb_bits; % Scale energy per bit (because of downsampling)
e_b = e_b/2; % power bandpass = 1/2 power complex envelope
N0 = e_b/ratio;
P_n = 2*fs*N0;
noise = sqrt(P_n/2)*(randn(length(s),1)+1i*randn(length(s),1));

end