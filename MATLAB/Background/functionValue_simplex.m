function f = functionValue_simplex(param)
    control_params = [5 6; 4 9;5 10; 4.5 5.5; 4 11];
    competing_pathway1 = [5 6; 5.1 5.5;5 10; 4.5 5.5; 4 11];
    competing_pathway2 = [5 6; 4 9;4.8 5.2; 4.5 5.5; 4 11];
    prod_brkdwn = [5 6; 4 9;5 10; 4.5 5.5; 5.4 6];
    rate_consts = competing_pathway1;
    f = functionValueProdRate(param, rate_consts);
    %f = functionValueYield(param, rate_consts);
end

function f = functionValueProdRate(param,rate_const)
  %reactor_lumped_model
  T_1 = param(1)+ 273.15;
  T_2 = param(2)+ 273.15;
  t = param(3:4);
%C1 = [ A B C D E]
%C2 = [ A B C D E F G H]
  C = [1 1 0 0 0 2 0 0];
  a = rate_const(1,1); b = rate_const(1,2);
  a1 = rate_const(2,1); b1 = rate_const(2,2); a2 = rate_const(3,1); b2 = rate_const(3,2);
  a3 = rate_const(4,1); b3 = rate_const(4,2);a4 = rate_const(5,1); b4 = rate_const(5,2);
  theta = [a b a1 b1 a2 b2 a3 b3 a4 b4];
  A1 = 10^(theta(1));
  E1 = theta(2)*1000;
  A2 = 10^(theta(3));
  E2 = theta(4)*1000;
  A3 = 10^(theta(5));
  E3 = theta(6)*1000;
  A4 = 10^(theta(7));
  E4 = theta(8)*1000;
  A5 = 10^(theta(9));
  E5 = theta(10)*1000;
  k1 = A1*exp(-E1/T_1);
  k2 = A2*exp(-E2/T_1);
  k3 = A3*exp(-E3/T_1);
  options = odeset('RelTol',1e-10,'NonNegative', [1 2 4 5]); 
  
  [t_out,C_out] = ode15s(@reactor1_odefun,[0 t(1)],C(1:5),options,[k1 k2 k3]);
%   fprintf('%d %d %d\n',C_out(end,:)');
  k1 = A1*exp(-E1/T_2);
  k2 = A2*exp(-E2/T_2);
  k3 = A3*exp(-E3/T_2);
  k4 = A4*exp(-E4/T_2);
  k5 = A5*exp(-E5/T_2);
  options = odeset('RelTol',1e-10,'NonNegative', [1 2 4 5 6 7 8]); 
  F = 1/t(2);
  C_1 = C_out(end,1)*(1/t(1))/F;
  C_2 = C_out(end,2)*(1/t(1))/F;
  C_3 = C_out(end,3)*(1/t(1))/F;
  C_4 = C_out(end,4)*(1/t(1))/F;
  C_5 = C_out(end,5)*(1/t(1))/F;
  C_6 = C(6)*(1/t(2) -1/ t(1))/F;
%   figure(11)
%   plot(t_out, C_out(:,5));
  C = [C_1 C_2 C_3 C_4 C_5 C_6 0 0];
  if(C_6 > 0) 
      [t_out,C_out] = ode45(@reactor2_odefun,[0 t(2)],C(1:8),options,[k1 k2 k3 k4 k5]);
  end
  
  C_min = min([1 C_6]);
  %fprintf('%d\n',C_out(end,7));
  if C_min <=0
      C_out1=-1;
  else
      if C_6 >0
          C_temp = C_out(end,7);
          if C_temp < 0
             C_temp = 0;
          end
          C_out1 = C_temp*F; %- max(0,.70-(C_temp/C_min)^2)+.70^2;
      else
          C_out1 = -1;
      end
     % if C_out(end,7)/C_min < .2
     %     C_out1 = 0;
     % end
  end
  
  f = -C_out1;

end

function dC = reactor1_odefun(t,Y,theta)
  dC1 = -(theta(1)+theta(2))*Y(1)*Y(2);
  dC2 = -(theta(1)+theta(2))*Y(1)*Y(2);
  dC3 = theta(1)*Y(1)*Y(2)-theta(3)*Y(3)*Y(2);
  dC4 = theta(2)*Y(1)*Y(2);
  dC5 = theta(3)*Y(3)*Y(2);
  dC = [dC1 dC2 dC3 dC4 dC5]';
end

function dC = reactor2_odefun(t,Y,theta)
  dC1 = -(theta(1)+theta(2))*Y(1)*Y(2);
  dC2 = -(theta(1)+theta(2))*Y(1)*Y(2);
  dC3 = theta(1)*Y(1)*Y(2)-theta(3)*Y(3)*Y(2)-theta(4)*Y(3)*Y(6);
  dC4 = theta(2)*Y(1)*Y(2);
  dC5 = theta(3)*Y(3)*Y(2);
  dC6 = -theta(4)*Y(3)*Y(6);
  dC7 = theta(4)*Y(3)*Y(6)-theta(5)*Y(7);
  dC8 = theta(5)*Y(7);
  dC = [dC1 dC2 dC3 dC4 dC5 dC6 dC7 dC8]';
end
