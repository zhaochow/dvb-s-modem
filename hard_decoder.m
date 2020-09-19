function decoded_bits = hard_decoder(y, H, K)
%HARD_DECODER Summary of this function goes here
%   Detailed explanation goes here

f = zeros(1,size(H,1));
c = zeros(length(f),length(y));
c(1,:) = y;
c_n = ones(1,length(y));

for k = 1:K
%     for i = 1:length(f)
%         f(i) = mod(c(1,:)*H(i,:)',2);
%     end
    f = mod(c(1,:)*H',2);   % Faster
    for i = 1:size(c,2)
        for j = 1:length(f)
            if(H(j,i) == 1)
                c_n(i) = c_n(i) + 1;
                if(f(j) == 1)
                    c(c_n(i),i) = not(c(1,i));
                else
                    c(c_n(i),i) = c(1,i);
                end
            end
        end
        chck = sum(c(:,i))/c_n(i);   % Faster to calculate once
        if(chck > 0.5)
            c(1,i) = 1;
        elseif(chck == 0.5)
            c(1,i) = c(1,i);
        else
            c(1,i) = 0;
        end
        c_n(i) = 1;
    end
%     if(c(1,:)*H' == 0)
%         break;
%     end
    if(c(1,:) == y)
        break;          %Faster condition
    end
    c(2:end,:) = 0;
    y = c(1,:);
end

decoded_bits = y;

end

