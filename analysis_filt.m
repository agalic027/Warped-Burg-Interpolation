function [res, states] = analysis_filt(a, x, lambda)
% Analysis filter used to remember current states of known part of the
% signal
%
% States are later used as internal excitation of synthesis filter, in
% order to compute one-sided extrapolation

    res = 0;
    states = zeros(numel(a), 1);
    wrpd = x;
    for i = 1 : length(a)
        res = res + a(i) * wrpd;
        states(i) = wrpd(end);
        wrpd = ap_delay(wrpd, lambda);  
    end
end

