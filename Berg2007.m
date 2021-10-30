function [input] = Berg2007(input,k,kk)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% --------------- BERG (2007) INCREMENTAL METHOD --------------

% Number of iterations for Berg (2007):
nSteps=50;
nElements=3;

% Converting all parameters to (nSteps,nElements) array size:
% Initialize a (nSteps,nElements) size array for each parameter
array3=ones(nSteps,nElements);

vol_clgrain=(input.vol_cl.*(1-input.por_clay))./(1-input.por_total);

for i=1:nSteps
    for j=1:nElements
        por_tot(i,j)=array3(i,j).*input.por_total;
        vol_clgr(i,j)=array3(i,j).*vol_clgrain;
        wat_sat(i,j)=array3(i,j).*input.Sw(k);
    end
end



for i=1:nSteps
   for j=1:nElements

       %wat_sat(i,j)=1-(k-1).*0.25;

       vol_elem1(i,j)=vol_clgr(i,j).*(1-por_tot(i,j));
       vol_elem2(i,j)=1-vol_elem1(i,j)-por_tot(i,j);
       vol_elem3(i,j)=(1-wat_sat(i,j)).*por_tot(i,j);

       if j==1
          disp_cond(i,j)=input.clay_cond;
          cem_exp(i,j)=input.cem_exp_clay;
          vol_elem_s(i,j)=vol_elem1(i,j)./nSteps;
       elseif j==2
          disp_cond(i,j)=input.sand_cond;
          cem_exp(i,j)=input.cem_exp_sand;
          vol_elem_s(i,j)=vol_elem2(i,j)./nSteps;
       elseif j==3
          disp_cond(i,j)=input.dnapl_cond;
          cem_exp(i,j)=input.sat_exp;
          vol_elem_s(i,j)=vol_elem3(i,j)./nSteps;
       end

       % Setting total water volume based on either 2 or 3 elements
       if nElements<=2
          phiT(i,j)=1-vol_elem1(i,j)-vol_elem2(i,j);
       elseif nElements==3
          phiT(i,j)=1-vol_elem1(i,j)-vol_elem2(i,j)-vol_elem3(i,j);
       end

       if ((i==1) && (j==1))
          xx=0;
          vol_sum=0.0;
       end

       vol_incr(i,j)=vol_sum+vol_elem_s(i,j);
       vol_sum=vol_incr(i,j);
       phi(i,j)=1-(vol_elem_s(i,j)./(phiT(i,j)+vol_incr(i,j)));

       % Hydrocarbon-first method
       if ((nElements<=2) && (i==1) && (j==1))
          wat_cond(i,j)=input.cond_wat.*(wat_sat(i,j).^input.sat_exp);
       else
          wat_cond(:,:)=input.cond_wat;
       end

       % For the very first iteration, the mixture conductivity is water conductivity
       % Setting mixture cond to bulk cond for each new step (set equal to water cond at start)   
       mix_cond(i,j)=xx;
       if ((i==1) && (j==1))
        mix_cond(:,:)=wat_cond(:,:);
       end

       % ------ CALCULATING BULK CONDUCTIVITY WITH NEWTON-RAPHSON ITERATION LOOP ------

       % Bulk (dispersed-water mixture) conductivity
       % Newton-Raphson method for finding the zero root of the non-linear HB equation:
       % f(x)=(disp_cond.*x.^((1./cem_exp)-1))-x.^(1./cem_exp)+(wat_cond.^(1./cem_exp).*por)-(wat_cond.^((1./cem_exp)-1).*por.*disp_cond);

       % Loop to calculate bec using Newton-Raphson
       % 'Tolerance' and 'maximum iterations' parameters for Newton-Raphson
       epc = 1; tol = 2.220446049250313e-016; maxiter = 10000;
       %epc = 1; tol = 0.000001; maxiter = 100000;
       format long;

       it=0;  % Iteration
       epc=1;   
       x=0.00001; % Initial approximation of the zero bulk electrical conductivity (bec)
       while ((epc > tol) & (it < 10000)) % Two termination criteria
             f=(disp_cond(i,j).*(x.^((1./cem_exp(i,j))-1)))-(x.^(1./cem_exp(i,j)))+((mix_cond(i,j).^(1./cem_exp(i,j))).*phi(i,j))-((mix_cond(i,j).^((1./cem_exp(i,j))-1)).*phi(i,j).*disp_cond(i,j));  % The function at x = x_k
             f1=(((1./cem_exp(i,j))-1).*disp_cond(i,j).*(x.^((1./cem_exp(i,j))-2)))-((1./cem_exp(i,j)).*(x.^((1./cem_exp(i,j))-1)));  % The derivative at x = x_k
             xx=abs(x-(f./f1));  % A new approximation for the root (added 'abs' to remove complex numbers - obtain magnitude value)
             epc=abs(xx-x); % Error
             x=xx;  % Root estimate for next round (within 'while' loop)
             it=it+1;  % Number of performed iterations
      end
  end
end

input.cond_Berg_2007(k,kk)=xx;

end





