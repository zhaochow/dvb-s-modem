close all; clear; clc;

H = [1 1 0 1 1 0 0 1 0 0;
     0 1 1 0 1 1 1 0 0 0;
     0 0 0 1 0 0 0 1 1 1;
     1 1 0 0 0 1 1 0 1 0;
     0 0 1 0 0 1 0 1 0 1];

Ht = H';
Ht = Ht(1:5,1:5);

combili = abs(Ht\eye(5));

H_new = (H'*combili)';
H_new = mod(H_new,2);

I = H_new(1:5,1:5);
P = H_new(1:5,6:end);

G = [P' I];

%% Test
bits_tx = randi(2,100,1)-1;
[encoded_bits,H] = ldpc_encoder(bits_tx);
