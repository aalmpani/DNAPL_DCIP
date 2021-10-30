function [input] = Revil2013(input,k,kk)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% input.elec_cond_Revil(k,kk)=(((input.Sw(k)^input.sat_exp)*(input.por_total^input.cem_exp_mix))*input.cond_wat);
input.elec_cond_Revil(k,kk)=(1/input.Ff(k,kk))*(input.Sw(k)^input.sat_exp)*(input.cond_wat);
input.surf_cond_Revil(k,kk)=(1/(input.Ff(k,kk))*(input.Sw(k)^(input.sat_exp-1))*(input.Beta_revil_low)*input.Qv(k,kk));
% Low-frequency conductivity
input.cond_low_revil(k,kk)= input.elec_cond_Revil(k,kk) + input.surf_cond_Revil(k,kk) ;
% High-frequency conductivity
% input.cond_high_revil(k,kk)= input.elec_cond_Revil(k,kk) + input.surf_cond_high_Revil(k,kk) ;
% Chargeability by Qi et al., (2018) Eq., A-5
input.M(k,kk)=((input.dens_gr*input.lamda*input.CECmix(k,kk))/(input.cond_wat*input.Sw(k)+input.dens_gr*input.Beta_revil*input.CECmix(k,kk)));
% Normalized chargeability
input.Mn(k,kk)=(input.Sw(k)^(input.sat_exp-1))*(input.por_total^(input.cem_exp_mix-1))*input.dens_gr*input.lamda*input.CECmix(k,kk);
% input.Mn(k,kk)=input.M(k,kk)/(1./input.cond_low_revil(k,kk));


end

