function out = polar_code_rate_recovery(in,N,P)
    
    out = zeros(1,N);
    pos = setdiff(0:N-1,P,'stable');
    out(pos+1) = in;

end