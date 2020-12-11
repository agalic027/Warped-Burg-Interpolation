function interpolated =  warped_burg_method(y, fs, unknown, p, lambda, alpha)
% Provide full interpolation based on warped Burg's method
%   interpolated  =   WARPED_BURG_METHOD(y, fs, unknown, p, lambda, alpha, debug)
%    
% Input Arguments:
%   y - signal to be interpolated, including the unknown samples
%   fs - signal sample frequency. If unimportant, set to 1
%   unknown - two element vector containing sample number of last known
%   sample in first known part and last unkown sample in the gap
%   p - model order for linear predictive filter
%   lambda - warping factor in range from -1 to 1 (Recommended value is 0.723).
%   See WARBURG for more details. 
%   alpha - steepness of crossfade window (Recommended value 3)
%
%   A. Galic 2020
%
%   Original paper:
%   Esquef, Paulo Antonio & Välimäki, Vesa & Roth, Kari & Kauppinen, Ismo. (2003).
%   Interpolation Of Long Gaps In Audio Signals Using The Warped Burg's Method. 
%
%   Functions ap_delay and genwfiltercoef are part of third party matlab
%   toolbox WarpTB provided here:
%   http://www.acoustics.hut.fi/software/warp/ 
%   and are not modified in any way
%

    % argument check
    if (abs(lambda) > 1)
        error('Lambda outside of [-1,1] range')
    end

    % assert input data as a column vector
    y = y(:);
    
    % compute number of unknown samples
    t = (0 : length(y) - 1) / fs;
    G = sum(t > unknown(1) & t<= unknown(2));
    
    % find left-side extrapolation
    y_left = y(t <= unknown(1));
    left_extrapolation = warped_burg_extrapolate ...
        (y_left, G, p, lambda);

    % find right-side extrapolation accounting for necessary time reversal
    y_right = y(t > unknown(2));
    right_extrapolation =  warped_burg_extrapolate ...
        (flipud(y_right), G, p, lambda);
    right_extrapolation =  flipud(right_extrapolation);

    % compute interpolation by crossfading both sided extrapolations
    extrap = crossfade_window(left_extrapolation, ...
        right_extrapolation, alpha);
    interpolated = cat(1, y_left, extrap, y_right);

end

