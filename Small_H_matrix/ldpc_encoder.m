function [encoded_bits,H] = ldpc_encoder(bit_stream)

% INPUTS:
% - bit_stream: input bit stream
%
% OUTPUTS:
% - encoded_bits: encoded bit stream

H = [1     0     0     0     0     0     1     1     1     0;
     0     1     0     0     0     1     0     1     0     0;
     0     0     1     0     0     1     0     1     0     1;
     0     0     0     1     0     0     0     1     1     1;
     0     0     0     0     1     1     1     0     0     1];

G = [0 1 1 0 1 1 0 0 0 0;
     1 0 0 0 1 0 1 0 0 0;
     1 1 1 1 0 0 0 1 0 0;
     1 0 0 1 0 0 0 0 1 0;
     0 0 1 1 1 0 0 0 0 1;];

bits = bit_stream;

i = mod(length(bits),5);
if i ~= 0
    bits = [bits; zeros(5-i,1)];
end

n = length(bits);
encoded_bits = zeros(2*n,1);
for i = 1:n/5
    block = bits(5*i-4:5*i);
    encoded_bits(10*i-9:10*i) = block' * G;
end

encoded_bits = mod(encoded_bits,2);
