function [input] = WaxmanSmits1968(input,k,kk)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Waxman and Smits equation for oil-bearing shaly sands, See also Butler
% and Knight (1998), 
input.elec_cond_WS(k,kk)=((input.Sw(k)^input.sat_exp)./input.Ff(k,kk))*(input.cond_wat);
input.surf_cond_WS(k,kk)=((input.Sw(k)^(input.sat_exp-1))./input.Ff(k,kk))*(((input.Beta_ws)*input.Qv(k,kk)));

input.cond_ws(k,kk)=input.elec_cond_WS(k,kk)+input.surf_cond_WS(k,kk); % S/m

end

