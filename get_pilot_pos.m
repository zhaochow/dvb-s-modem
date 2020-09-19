function n_est = get_pilot_pos(Dk)
%   Return the estimated pilot position __in terms of MATLAB index__
    
    summed_Dk = sum(abs(Dk),2);
    [max_val, max_idx] = max(summed_Dk);
    n_est = max_idx;
end