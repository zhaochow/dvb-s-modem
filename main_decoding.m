clc; clear; close all

nb_bits = 128;
cr = 1/2;
NB_PACKETS = 1000;
SNR = -10:0.2:10;
K = [1 2 3 4 5];
H0 = makeLdpc(nb_bits, nb_bits/cr,0,1,3); % Create initial parity check matrix of size 128 x 256

infobits = (randi(2,nb_bits,NB_PACKETS)-1);

[checkbits, H] = makeParityChk(infobits, H0, 0); % Encode information bits (128 x number of packets) and generate final parity check matrix
encodedbits = [checkbits' infobits'];
BER = zeros(1+2*length(K),length(SNR));

for i = 1:length(SNR)
    for p = 1:NB_PACKETS
        % CAREFUL: USE BPSK IN CHANNEL_TXRX
        [encodedbits_err, symb_rx, errors, P_n] = channel_txrx(encodedbits(p,:)', SNR(i));
        encodedbits_err =  encodedbits_err';

        infobits_err = encodedbits_err(size(checkbits,1)+1:end);
        errors_info = sum(abs(infobits(:,p)' - infobits_err));
        error_rate = (errors_info/nb_bits);
        BER(1,i) = BER(1,i) + error_rate;
                
        %======HARD DECODING========
        for k = 1:length(K)
            decoded_bits = hard_decoder(encodedbits_err, H, K(k));
        
            infobits_decoded = decoded_bits(size(checkbits,1)+1:end);
            errors_decoded = sum(abs(infobits(:,p)' - infobits_decoded));

            error_rate_ldpc = (errors_decoded/nb_bits);
            BER(k+1,i) = BER(k+1,i) + error_rate_ldpc;
        end
        %======SOFT DECODING========
        for k = 1:length(K)
            decoded_bits = ldpc_soft_decoder(symb_rx,P_n/2,H,K(k));
            decoded_bits = decoded_bits';
        
            infobits_decoded = decoded_bits(size(checkbits,1)+1:end);
            errors_decoded = sum(abs(infobits(:,p)' - infobits_decoded));

            error_rate_ldpc = (errors_decoded/nb_bits);
            BER(k+1+length(K),i) = BER(k+1+length(K),i) + error_rate_ldpc;
        end
    end
    BER(:,i) = BER(:,i)/NB_PACKETS;
end

%%
figure
semilogy(SNR,BER(:,:))
hold on

load('ber_pam2.mat')
semilogy(bers0_pam2.data{1,1},bers0_pam2.data{1,2},'--')
% load('ber_pam4.mat')
% semilogy(bers0_pam4.data{1,1},bers0_pam4.data{1,2},'--')
% load('ber_qam16.mat')
% semilogy(bers0_qam16.data{1,1},bers0_qam16.data{1,2},'--')
% load('ber_qam64.mat')
% semilogy(bers0_qam64.data{1,1},bers0_qam64.data{1,2},'--')
% legend('2-PAM (ideal)','4-PAM (ideal)','16-QAM (ideal)','64-QAM (ideal)')
ylim([1e-5 1])
