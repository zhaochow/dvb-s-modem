function decoded_bits = ldpc_soft_decoder(bit_stream,sigma,H,max_it)

% INPUTS:
% - bit_stream: input bit stream
% - sigma: noise power / 2
% - H: parity check matrix
% - max_it: maximum number of iterations
%
% OUTPUTS:
% - decoded_bits: decoded bit stream

phi = @(x) -log10(tanh(0.5*x));
% phi = @(x) log10((exp(x)+1)./(exp(x)-1));

bits = bit_stream;
var1 = var([bits(bits<0)+1; bits(bits>0)-1]);

[n,m] = size(bits);
if m > 1
    bits = reshape(bits,[m,1]);
end
decoded_bits = zeros(n,1);

% i for check nodes, j for variable nodes
[I,J] = size(H); % I check nodes, J variables nodes
L_Q = zeros(J,1); % Soft decision matrix
L_q = zeros(I,J); % q_ij = message from var j to check i
L_r = zeros(I,J); % r_ij = message from check i to var j
for k = 1:n/J
    c = bits(J*(k-1)+1:J*k);
    L_c = -2*c/(var1);
    % Initialization
    for i = 1:I % for the number of rows in H
        cols_1 = find(H(i,:)==1); % get column indexes with 1
        L_q(i,cols_1) = L_c(cols_1)';
    end
    it = 0;
    syndrome = sum(mod(c'*H',2));
    while syndrome ~= 0 && it < max_it
        % Computation of check nodes
        for i = 1:I
            cols_1 = find(H(i,:)==1);
            alphas = sign(L_q(i,cols_1));
            betas = abs(L_q(i,cols_1));
            for l = 1:length(cols_1) % for each variable node used in check node
                alphas_cpy = alphas;
                betas_cpy = betas;
                alphas_cpy(l) = [];
                betas_cpy(l) = [];
                prod_alpha = prod(alphas_cpy);
%                 L_r(i,cols_1(l)) = prod_alpha * phi(sum(phi(betas_cpy)));
                L_r(i,cols_1(l)) = prod_alpha * min(betas_cpy);
            end
        end
        % Computation of variable nodes
        for j = 1:J
            rows_1 = find(H(:,j)==1); % get row indexes with 1
            message_from_check_nodes = L_r(rows_1,j);     
            for l = 1:length(rows_1) % for each check node linked to variable node
                message_from_check_cpy = message_from_check_nodes;
                message_from_check_cpy(l) = [];
                L_q(rows_1(l),j) = L_c(j) + sum(message_from_check_cpy);
            end
            % Soft decision
            L_Q(j) = L_c(j) + sum(message_from_check_nodes);
            % Hard decision
            if L_Q(j) < 0
                c(j) = 1;
            else
                c(j) = 0;
            end
        end
        syndrome = sum(mod(c'*H',2));
        it = it + 1;
    end
    decoded_bits(J*(k-1)+1:J*k) = c;
end

if m > 1
    decoded_bits = reshape(decoded_bits,[1,m]);
end