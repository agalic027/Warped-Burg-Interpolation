function extrapolated = warped_burg_extrapolate(data, G, p, lambda)
% Perform one-sided extrapolation using warped Burg's method
%
% extrapolated = warped_burg_extrapolate(data, G, p, lambda)
%   Computes frequency-warped extrapolation on G zero samples following
%   data using frequency-warped linear prediction

    data = data(:);
    N = numel(data);
    extrap = cat(1, data, zeros(G, 1));
    extrap(1 : N) = data;
    wb = warburg(data, p, lambda);
    [~, zf] = analysis_filt(wb, data, lambda);   
    extrap(N + 1 : end) = synthesis_filt(wb, extrap(N + 1 : end),...
        lambda, zf);
    extrapolated = extrap(N + 1 : end);
end

