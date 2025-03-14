function xdot = dvpH(x, ux, uy, params)
% dvpH - Computes the time derivative of the state vector x
%
% Syntax: xdot = dvpH(x, ux, uy, params)
%
% Inputs:
%   x      - State vector
%   ux     - Input flow rate 1
%   uy     - Input flow rate 2
%   params - Structure containing parameters
%
% Outputs:
%   xdot   - Time derivative of the state vector x

  dx.rh = (ux + params.flow.Q2 + uy - (params.constants.c * (sqrt(x.rh - params.geometry.h0)))) / params.geometry.Ar;

  dx.wa = ((ux * (params.concentration.wa1 - x.wa)) + (params.flow.Q2 * (params.concentration.wa2 - x.wa)) + (uy * (params.concentration.wa3 - x.wa))) / params.geometry.Vr;

  dx.wb = ((ux * (params.concentration.wb1 - x.wb)) + (params.flow.Q2 * (params.concentration.wb2 - x.wb)) + (uy * (params.concentration.wb3 - x.wb))) / params.geometry.Vr;

  dx.hta = -ux / params.geometry.Ata;

  dx.htt = -params.flow.Q2 / params.geometry.Att;

  dx.htb = -uy / params.geometry.Atb;

  dx.htc = params.constants.c / params.geometry.Atc;

  xdot = dx';
  
end
