function delta_est = get_phase_diff(Dk, n_est,Tsymb)
    K = size(Dk,2);
    k_vect = 1:K;
    Dk_ang = angle(Dk(n_est,:))./(2*pi*k_vect*Tsymb);
    delta_est = (-1/K)*sum(Dk_ang);
    
end