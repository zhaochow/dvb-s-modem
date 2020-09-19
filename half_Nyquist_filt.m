function [h_RRC,t] = half_Nyquist_filt(beta,T,fs,nb_taps)

% INPUTS:
% - symb_tx: vector of input symbols
% - beta: roll-off factor
% - T: symbol period
% - fs: sampling frequency
% - nb_taps: number of taps of the filter (odd nb of taps only)
%
% OUTPUTS:
% - h_RRC: halfroot Nyquist filter
% - t: time vector

stepOffset = 1/nb_taps*fs;
highestFreq = stepOffset*(nb_taps-1)/2;
freqGrid = linspace(-highestFreq,highestFreq,nb_taps);
H_RC = zeros(nb_taps,1);

for i = 1:nb_taps
    if abs(freqGrid(i)) <= (1-beta)/(2*T)
        H_RC(i) = T;
    elseif abs(freqGrid(i)) <= (1+beta)/(2*T)
        H_RC(i) = T/2*(1+cos(pi*T/beta* ( abs(freqGrid(i))-((1-beta)/(2*T)) )));
    end
end

h_RC = ifft(ifftshift(H_RC),'symmetric');
H_RC_norm = H_RC/h_RC(1);

% plot(abs(H_RC_norm2))
H_RRC = sqrt(H_RC_norm);

delta_t = 1/fs;
t = (-(nb_taps-1)/2:(nb_taps-1)/2)*delta_t;

% figure
% stem(t,fftshift(h_RC/h_RC(1)))
% grid on;
h_RRC = ifft(ifftshift(H_RRC),'symmetric');
h_RRC = fftshift(h_RRC);
% figure
% plot(freqGrid,abs(H_RRC))
% figure
% stem(t,h_RRC)
