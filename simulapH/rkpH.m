function x = rkpH(x0, ux, uy, h, params)
% rkpH - Runge-Kutta 4th order method for solving ODEs
%
% Syntax: x = rkpH(x0, ux, uy, h, params)
%
% Inputs:
%   x0     - Initial state structure
%   ux     - Input flow rate 1
%   uy     - Input flow rate 2
%   h      - Step size
%   params - Structure containing parameters
%
% Outputs:
%   x      - State structure at the next time step

  k1 = dvpH(x0, ux, uy, params);

  x_temp = structfun(@(field, k) field + 0.5 * h * k, x0, k1, 'UniformOutput', false);
  k2 = dvpH(x_temp, ux, uy, params);

  x_temp = structfun(@(field, k) field + 0.5 * h * k, x0, k2, 'UniformOutput', false);
  k3 = dvpH(x_temp, ux, uy, params);
  
  x_temp = structfun(@(field, k) field + h * k, x0, k3, 'UniformOutput', false);
  k4 = dvpH(x_temp, ux, uy, params);

  x = structfun(@(field, k1, k2, k3, k4) field + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4), x0, k1, k2, k3, k4, 'UniformOutput', false);
end