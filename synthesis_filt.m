function res = synthesis_filt(a, x, lambda, init_state)
% Synthesize samples on unknown part of the signal, using internal
% excitation provided by analysis filter.

    [sigma, gain] = genwfiltercoef(a, lambda);
    res = zeros(size(x));
    aux_states = zeros(numel(a), 1);
    init_state = cat(1, init_state, zeros(numel(a) - numel(init_state), 1));
    
    
    for i = 1 : length(x)
        aux_states(1) = gain * (x(i) - sigma' * init_state);
        for j = 2 : numel(a)
            aux_states(j) = init_state(j - 1) + ...
                lambda * (init_state(j) - aux_states(j - 1));
        end
        res(i) = aux_states(1);
        init_state = aux_states;
    end

end

function [aw,g]=genwfiltercoef(a,lambda)

% INPUTS
% a: coefficients of the denominator (original ARMA)
% lambda: coefficient of the all pass filter (frequency warping factor)
% lambda must be bettwen -1 and 1 

% OUTPUTS
% aw: modified coefficients of the feedback loop (denominator) 
% g: direct gain given by 1-lambda*aw(1)= 1/aw(0) = 1/sigma(0)

a(1)=[]; % removes the first coefficient of the monic polynomial
N=length(a); % number of feedback coefficients  
%  Note that if there are N alpha (a) feedback coefficients, there will be N+1 
% sigma (aw) feedback coefficients. 
if N>0 % if the structure has feedback (or is IIR)
    aw=zeros(N+1,1); % defines an empty vector for the modified feedback parameters 
    S=zeros(N,1); % auxiliary vector
    aw(N+1)=lambda*a(N); S(N)=a(N); % computes the N+1 sigma coeff
    for k=N:-1:2
        S(k-1)=a(k-1)-lambda*S(k);
        aw(k)=lambda*S(k-1)+S(k);  % computes the kth sigma coeff from N to 2
    end
    aw(1)=S(1); % computes the first sigma coeff.
    g=1/(1-lambda*S(1)); % computes the direct gain: g=1/sigma(0)     
elseif N==0 % if the structure is FIR
    aw=[];
    g=1;
end
end