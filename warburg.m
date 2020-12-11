function [a, e] = warburg(x, order, lambda)
% Frequency-warped Burg method for computing AR coefficients
% [a, e] = WARBURG(x, order, lambda)
%
% Lambda needs to be between -1 and 1. Positive values provide better
% resolution in lower frequency range, negative values provide better
% resolution in higher frequency range. Running this function with lambda
% set to 0 is the same as running built-in arburg function.

    % argument check
    if (abs(lambda) > 1)
        error('Lambda outside of [-1,1] range')
    end
    
    % assert x as a column vector
    x = x(:);
    N = length(x);
    
    % initiate error and first AR coefficient
    e = (1 / N) * (x' * x);
    a = 1;
    
    % modified levinson durbin recursion
    f_old = x;
    b_old = x;
    for i = 1 : order
        
        %Truncate error vectors,warp b instead of delay
        f = f_old(2 : end);
        b_old = ap_delay(b_old, lambda);
        b = b_old(2 : end);
        
        %Compute reflection coefficients
        k = (-2 * b' * f) / (f' * f + b' * b);
        
        % Compute LPC coeffs and prediction error via levinson
        a = [a; 0] + k * flipud(conj([a; 0]));
        e = (1 - k' * k) * e;
        
        %Update prediction errors
        f_old = f + k * b;
        b_old = b + conj(k) * f;
        
    end
    a = a.';

end

