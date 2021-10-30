clc
clear

% Mobility of the counterions in the Diffuse layer
beta_dl=5.2*10E-8; % m2/sV
% Fraction of counterions (or amount of cations per mass of grains) located n the Stern layer
fraction=0.90;

k=1;
for i=0:0.1:1
    
    % Water conductivity (S/m)
    cond_wat(k)=i; % S/m

    
    Beta1(k)=(1-(0.83*exp(-0.5*cond_wat(k))))*3.83; % Waxman and Smits (1968)
    
    Beta2(k)=(1-(0.83*exp(-0.5.*cond_wat(k))))*3.83; % Butler and Knight (1998)
    
    Beta3(k)=(1-0.6*exp(-cond_wat(k)/0.013))*4.78*1E-08*((1E+06)*(96.32)); % Revil et al., (1998) - convert from m2/sV to S cm3/m meq
    
    Beta4(k)=(1-0.6*exp(-cond_wat(k)/0.013))*0.0461*100; % Ugbo et al., (2009) - convert from S ml/cm meq to S cm3/m meq

    Beta5(k)=(beta_dl*(1-fraction))*(((1E+06)*(96.32))); % Revil et al., (2013a) - convert from m2/sV to S cm3/m meq
    
    Beta6(k)=(1-0.6*exp(-77*cond_wat(k)))*0.046*100; % Tenchov (2014) - convert from S ml/cm meq to S cm3/m meq
    

    k=k+1;
end

plot(cond_wat,Beta1,'bo',cond_wat,Beta2,'y-',cond_wat,Beta3,'k-',cond_wat,Beta4,'g+',cond_wat,Beta5,'r*-',cond_wat,Beta6,'bo');
% title('B - Different versions')
xlabel('Water conductivity (S/m)')
ylabel('B (S.cm3/m.meq)')
grid on
grid minor
legend('Waxman and Smits (1968)','Butler and Knight (1998)','Revil et al., (1998)','Ugbo et al., (2009)','Revil (2013a)','Tenchov (2004)' )

% semilogy(cond_wat,Beta4,'r*-',cond_wat,Beta5,'k*-');
% title('Equivalent conductance of the counterions B')
% xlabel('Water conductivity (S/m)')
% ylabel('B')
% grid on
% grid minor
% legend('Revil (2013a)','Butler and Knight (1998)' )