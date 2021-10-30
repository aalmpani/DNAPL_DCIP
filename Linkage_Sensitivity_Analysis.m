%***************************************************************************************************
% This code compares the different published models.
% Waxman and Smits (1968)
% Berg (2007)
% Revil (2013a)
%***************************************************************************************************
clc
clear all
close all

% If Comp=1 then comparison between Revil and WS surface conductivites. Otherwise, if Comp=0 comparison of
% conductivity of all the models
Comp=0; % 0 or 1

% Mean grain diameter of sand
input.dia_sand=0.2; % in cm
% Mean grain diameter of clay
input.dia_clay=0.0002; % in cm
% Porosity of sand
input.por_sand=0.30;
% Porosity of clay
input.por_clay=0.40; % ?
% Cementation exponent of sand
input.cem_exp_sand=1.5; % % clean sand (Revil and Cathles, 1999)
% Cementation exponent of clay
input.cem_exp_clay=1.5; % % pure clay (Berg, 1996)
% A mixed cementation exponent
input.cem_exp_mix=1.5;
% Saturation exponent
input.sat_exp=1.5;
% Grain density
input.dens_gr=2.650; % typically grain density for sedimentary rocks in g/cm3
% Water resistivity
input.res_wat=2; % Ohm-m

% Fraction of counterions (or amount of cations per mass of grains) located n the Stern layer
% input.fraction=0.50;
% Surface charge density of sand
input.Qs_sand=0.64; % C/m2
input.Qs_sand=input.Qs_sand*(1/(96.32*10000)); % convert to meq/cm2
% Surface charge density of clay minerals
input.Qs_clay=0.32; % C/m2
input.Qs_clay=input.Qs_clay*(1/(96.32*10000)); % convert to meq/cm2

% Inputs for Berg (2007)
input.clay_cond=150./1000;
input.sand_cond=0.1./1000;
input.dnapl_cond=5.6e-09./1000;

max_Snw=0.9;

kk=1;

for ii=0:0.15:0.3
    k=1;
    for i=0:0.1:max_Snw % DNAPL saturation
%     for i=0:1:10 % % #AA#
        % Clay content
        input.vol_cl=ii;
        input.vol_cl_matrix(kk)=input.vol_cl;
        % Mobility of the counterions in the Stern layer in m2/sV
        if input.vol_cl==0
            input.fraction=0.5;
            input.beta_sl=5.2*10E-8; 
            input.beta_dl=5.2*10E-8; % m2/sV
        else
            input.fraction=0.90;
            input.beta_sl=1.5*10E-10; 
            input.beta_dl=5.2*10E-8; % m2/sV
        end

        % DNAPL saturation
        input.Snw(k)=i;
        % Water saturation
        input.Sw(k)=1-input.Snw(k); 
%         input.Sw(k)=1; % #AA#
%         input.res_wat(k)=i+0.00001; % #AA#
%         input.cond_wat_plot(k)=input.res_wat(k); % #AA#

%         input.cond_wat=1/(input.res_wat(k)); % S/m #AA#
        input.cond_wat=0.5; % #AA#

        input.vol_cl_matrix(kk)=input.vol_cl; % For plots
        if input.vol_cl==0 % Then there is no clay
            input.por_total=input.por_sand;
        else
            % Total porosity
            input.por_total=(input.por_sand*(1-input.vol_cl))+(input.por_clay*input.vol_cl);
        end
        % Formation factor
        input.Ff(k,kk)=input.por_total^(-input.cem_exp_mix);
        % Equivalent conductance of the counterions as a function of solution conductivity
%         input.Beta_ws=(1-0.6*exp(-input.cond_wat/0.013))*0.0461*100; % Ugbo et al., (2009)
%         input.Beta_ws=(1-0.6*exp(-77*input.cond_wat))*0.046*100; % Tenchov (2004)
%         input.Beta_ws=(1-(0.83*exp(-0.5.*input.cond_wat)))*3.83; % Butler and Knight (1998) (S*cm3)/(m*meq)
        input.Beta_ws=(1-0.6*exp(-input.cond_wat/0.013))*4.78*(1E-08)*((1E+06)*(96.32)); % Revil et al., 1998 Eq., 19 convert from m2/sV to S cm3/m meq
%         input.Beta_ws=3.83*(1-(0.83*exp(-0.5.*input.cond_wat)))*(1E-08)*((1E+06)*(96.32));
        % Apparent mobility of the counterions for the polarization associated with the quadrature conductivity
        input.lamda=input.beta_sl*input.fraction*((1E+06)*(96.32)); % This is for low conductivity M and Mn
        input.Beta_revil_low=(input.beta_dl*(1-input.fraction))*((1E+06)*(96.32)); % Revil (2013a) convert from m2/sV to S cm3/m meq
        % Apparent mobility of the counterions for surface conduction
        input.Beta_revil=input.Beta_revil_low+input.lamda; % This is for M and high conductivity both DL and SL
        % Specific surface area of sand
        input.Ssp_sand=6/(2.65*input.dia_sand);
        % Specific surface area of clay
%         input.Ssp_clay=6/(1.2*input.dia_clay);
        input.Ssp_clay=80000; % 
        % Calculate CEC based on Ssp and Qs for sand and clay minerals
        input.CEC_clay=input.Qs_clay*input.Ssp_clay;
%         input.CEC_clay=0.03;
        input.CEC_sand=input.Qs_sand*input.Ssp_sand;
        % Mixed CEC, The conversion has been also confirmed with Revil et al., (2013) page 10.
        input.CECmix(k,kk)=(input.CEC_clay*input.vol_cl)+(input.CEC_sand*(1-input.vol_cl)); % 0.03 is in meq/g which I convert to C/kg multiplying by 963.20 based on Abdulsamad et al., 2020
        % Charge density
        input.Qv(k,kk)=(input.dens_gr*input.CECmix(k,kk)*((1-input.por_total)/input.por_total)); % C/m3, If I divide by 96.32*10E+6 then I convert to eq/L
        
        
%***************************************************************************************************
% Call function for each model
%***************************************************************************************************
        % Waxman and Smits (1968)
        [input] = WaxmanSmits1968(input,k,kk);
        
        % Berg (2007)
        [input] = Berg2007(input,k,kk);
        
        % Revil (2013a) for Res, M and Mn
        [input] = Revil2013(input,k,kk);
        
        
        k=k+1;
        
    end % End of Snw
    kk=kk+1;
end % End of Cl 



% Dimensionless number that characterizes polarization
R=input.lamda/(input.Beta_revil);

% Converting conductivity (S/m) to resistivity (Ohmm)
input.res_ws=1./input.cond_ws; % Ohmm
input.res_low_revil=1./input.cond_low_revil; % Ohmm
% input.res_high_revil=1./input.cond_high_revil; % Ohmm
input.res_Berg_2007=1./input.cond_Berg_2007; % Ohmm
% input.res_Berg_2007=0;
% Calculate water content
input.wat_content=input.Sw*input.por_total;

%***************************************************************************************************
%******************************************* P L O T S *********************************************
%***************************************************************************************************

%---------------------------------------------------------------------------------------------------
% RESISTIVITY VS DNAPL SATURATION
%---------------------------------------------------------------------------------------------------

        for jj=1:3 % 1 for WSM, 2 for BM and 3 for RM
            for ii=1:length(input.vol_cl_matrix)
                if jj==1
                    p1=semilogy(input.Snw,input.res_ws(:,ii),'b*--');
                    xt=max(input.Snw)+0.004;
                    yt=input.res_ws(length(input.res_ws(:,ii)),ii);
                    str=strcat('\leftarrow Cl=',num2str(input.vol_cl_matrix(ii)));
                    text(xt,yt,str)
                elseif jj==2
                    p2=semilogy(input.Snw,input.res_Berg_2007(:,ii),'k-');
%                     if ii~=1
%                         xt=max(input.Snw)+0.004;
%                         yt=input.res_Berg_2007(length(input.res_Berg_2007(:,ii)),ii);
%                         str=strcat('\leftarrow Cl=',num2str(input.vol_cl_matrix(ii)));
%                         text(xt,yt,str)
%                     end
                else 
                    p3=semilogy(input.Snw,input.res_low_revil(:,ii),'r--');
                end
                hold on
            end 
        end
%         ylim([10 1000])
    ylabel('Bulk Resistivity (Ohm.m)')
%     title('Comparative Analysis of Surface Resistivity','fontweight','bold','fontsize',14)
    xlabel('DNAPL Saturation')
    legend([p1 p2 p3],'Waxman and Smits (1968)','Berg (2007)','Revil (2013a)')
    grid on
    grid minor
    
%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
% CHARGEABILITY AND NORMALIZED CHARGEABILITY VS DNAPL SATURATION
%---------------------------------------------------------------------------------------------------

%     subplot(2,1,1)
%         for jj=1:1 % 1 for Revil and 2 for WS
%             for ii=1:length(input.vol_cl_matrix)
%                     p1=semilogy(input.Snw,input.M(:,ii),'*-');
%                     xt=max(input.Snw)+0.004;
%                     yt=input.M(length(input.M(:,ii)),ii);
%                     str=strcat('\leftarrow Cl=',num2str(input.vol_cl_matrix(ii)));
%                     text(xt,yt,str)
%                 hold on
%             end 
%         end
%     ylabel('Chargeability')
%     xlabel('DNAPL Saturation')
%     grid on
%     grid minor 
%     % NORMALIZED CHARGEABILITY
%     subplot(2,1,2)
%         for jj=1:1 % 1 for Revil and 2 for WS
%             for ii=1:length(input.vol_cl_matrix)
%                     p1=semilogy(input.Snw,input.Mn(:,ii),'*-');
%                     xt=max(input.Snw)+0.004;
%                     yt=input.Mn(length(input.Mn(:,ii)),ii);
%                     str=strcat('\leftarrow Cl=',num2str(input.vol_cl_matrix(ii)));
%                     text(xt,yt,str)
%                 hold on
%             end 
%         end
%     ylabel('Normalized Chargeability (S/m)')
%     xlabel('DNAPL Saturation')
%     grid on
%     grid minor 

%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------


%---------------------------------------------------------------------------------------------------
% RESISTIVITY VS WATER CONDUCTIVITY
%---------------------------------------------------------------------------------------------------

%         for jj=1:3 % 1 for Revil and 2 for WS
%             for ii=1:length(input.vol_cl_matrix)
%                 if jj==1
%                     p1=plot(input.cond_wat_plot,input.res_ws(:,ii),'b*');
%                     xt=max(input.Snw)+0.004;
%                     yt=input.res_ws(length(input.res_ws(:,ii)),ii);
%                     str=strcat('\leftarrow Cl=',num2str(input.vol_cl_matrix(ii)));
%                     text(xt,yt,str)
%                 elseif jj==2
%                     p2=plot(input.cond_wat_plot,input.res_Berg_2007(:,ii),'k-');
%                     if ii~=1
%                         xt=max(input.Snw)+0.004;
%                         yt=input.res_Berg_2007(length(input.res_Berg_2007(:,ii)),ii);
%                         str=strcat('\leftarrow Cl=',num2str(input.vol_cl_matrix(ii)));
%                         text(xt,yt,str)
%                     end
%                 else 
%                     p3=plot(input.cond_wat_plot,input.res_low_revil(:,ii),'r--');
%                 end
%                 hold on
%             end 
%         end
%     ylabel('Bulk Resistivity (Ohm.m)')
% %     title('Comparative Analysis of Surface Resistivity','fontweight','bold','fontsize',14)
%     xlabel('Water resistivity (Ohm.m)')
%     xlim([0 10])
%     legend([p1 p2 p3],'Waxman and Smits (1968)','Berg (2007)','Revil (2013a)')
%     grid on
%     grid minor 
    
%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------




