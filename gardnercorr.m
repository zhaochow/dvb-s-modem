function [s_corr,idx,err] = gardnercorr(s,fs,f_symb,nb_symb,K)

n = length(s);
M = fs/f_symb; % Upsampling factor (need to be even)


err = zeros(floor(n/M),1);
idx = zeros(floor(n/M),1);
s_corr = zeros(floor(n/M),1);
s_corr(1) = s(1);
idx(1) = 1;

i = 1;
index = M+1;
while index < n && i < nb_symb
    idx(i+1) = index;
    s_corr(i+1) = s(index);
    err(i+1) = err(i) + K*real(s(index - M/2)*(conj(s(index)) - conj(s(index-M))));
    index = (i+1)*M + 1 - round(err(i+1)*M);
    i = i + 1;
end

s_corr = s_corr(1:i);
idx = idx(1:i);
