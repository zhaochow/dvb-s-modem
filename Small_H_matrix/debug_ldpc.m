close all; clear; clc;

bits_tx = randi(2,100,1)-1;
[encoded_bits,H] = ldpc_encoder(bits_tx);
% encoded_bits'

a = randi(5,1,1);
bits_corr = encoded_bits;
for i = 1:length(a)
    bits_corr(a(i)) = not(bits_corr(a(i)));
end
% bits_corr'
err = encoded_bits - bits_corr;
e1 = sum(abs(err))
mod = 2*bits_corr - 1;

sigma = 0.7;
max_it = 3;

decoded_bits = ldpc_soft_decoder(mod,sigma,H,max_it);
% decoded_bits'
% demod = (decoded_bits + 1)/2;
err2 = encoded_bits - decoded_bits;
e2 = sum(abs(err2))