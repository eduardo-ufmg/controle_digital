function [x, pH, xc, pHc] = simrk_pH(x0, u1, u2, h, Ts, params, Kas)
% simrk_pH - Simulates the pH dynamics using Runge-Kutta 4th order method
%
% Syntax: [x, pH, xc, pHc] = simrk_pH(x0, u1, u2, h, Ts, params, Kas)
%
% Inputs:
%   x0     - Initial state structure
%   u1     - Input flow rate 1
%   u2     - Input flow rate 2
%   h      - Step size
%   Ts     - Total simulation time
%   params - Structure containing parameters
%   Kas    - Structure containing acid dissociation constants
%
% Outputs:
%   x      - State structure at the final time step
%   pH     - pH value at the final time step
%   xc     - State structures at each time step
%   pHc    - pH values at each time step

    xc = repmat(x0, 1, Ts/h);
    pHc = zeros(1, Ts/h);
    xc(1) = x0;
    pHc(1) = computepH(x0, Kas);

    for i = 2:Ts/h
        xc(i) = rkpH(xc(i-1), u1, u2, h, params);
        pHc(i) = computepH(xc(i), Kas);
    end
    
    x = xc(end);
    pH = pHc(end);

end


function ph = computepH(xc, Kas)
    r = roots([1 (Kas.ka1 - xc.wa) (Kas.ka1*Kas.ka2-Kas.ka1 * xc.wa-Kas.kw-Kas.ka1 * xc.wb) -(Kas.ka1*Kas.kw+Kas.ka1*Kas.ka2 * xc.wa+2*Kas.ka1*Kas.ka2 * xc.wb) -Kas.ka1*Kas.ka2*Kas.kw]);
    ph = -log10(r(r > 0 & r < 1));
end
