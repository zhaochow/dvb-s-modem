function Dk = diffcorr(y, a, K)
    M = length(y); % Length frame
    N = length(a); % Length pilot
    y = [y;zeros(N,1)]; % Add padding
    Dk = zeros(M, K);
    for k = 1:K
        for n = 1:M
            Dk(n,k) = sum(conj(y(n+k+1:n+N)).*a(k+1:N).*y(n+1:n+N-k).*conj(a(1:N-k)));
        end
        Dk(:,k) = Dk(:,k)/(N-k);
    end
end