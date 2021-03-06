function simRamps 
global G_temp
Temps = [50, 70, 90, 100, 130];
numTemp = length(Temps);
tau_0 = .5;
tau_f = 20;
%C3 =1;
S = 2/3;
alpha = -log((1-S));
%Vd = 0;
%Vr = 
cond = [];
control_params = [5 5; 4 9;5 10; 5.5 5.2; 4 11];
competing_pathway1 = [6 5; 6.1 5.5;6 10; 6.5 5.2; 5 11];
competing_pathway2 = [6 5; 5 9;5.8 5.2; 6.5 5.2; 5 11];
prod_brkdwn = [6 5; 5 9;6 10; 6.5 5.2; 6.4 6];
rate_const = control_params;

std_dev = .03;

for ExpNum = 1:numTemp
    t = 0;
    rxn_time = tau_0;
    while(rxn_time < tau_f)
        if t < 60
          rxn_time = tau_0;
        else
          rxn_time = tau_0 + alpha*(t-60)/60;
        end

        %tau_integral=max(tau_0,Slope*exp(-Vd/Vr*alpha)*((t-60)/60 + tau_0/alpha));

        T = Temps(ExpNum);
        C = sing_react_oracle([rxn_time T], rate_const);
        C = C + C*std_dev.*randn(1,1);
        cond(end+1,:) = [T rxn_time C];
        %Simulated actual time
        t = t+30;
    end
end
numCond = length(cond(:,1));
theta_0 = [1 1; 1 1; 1 1; 1 1; 1 1];

%Maximum likelihood estimate parameter estimates of kinetic parameters
%
options = optimoptions('fmincon','MaxFunEvals',100,'FunValCheck','on','TolX',1e-100,'TolFun',1e-100,'TolCon',1e-100,'Algorithm','sqp', 'MaxIter', 200);
%   min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%      X                     C(X) <= 0, Ceq(X) = 0   (nonlinear constraints)
%                               LB <= X <= UB        (bounds)
% mycon returns vectors [c, ceq]
% syntax: x = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
% syntax: x = fmincon(@(x) myfun(x,a1),[1;2],[],[],[],[],[],[],@(x) mycon(x,a2),options);

%Use half the values as a training set and the other for validation
V=eye(4);
trainConds = cond(1:2:end,:);
C = trainConds(:,[3,4,8,9]);
validConds = cond(2:2:end,:);
num_train_cond = length(trainConds(:,1));
LB = [1,1;1,1;1,1;1,1;1,1];
UB = [12,12;12,12;12,12;12,12;12,12];
weight_inv = zeros(4*num_train_cond,4*num_train_cond);    
for i = 1:num_train_cond
    weight_inv(4*i-3:4*i,4*i-3:4*i) = V;
end
weight = weight_inv^(-1);
[theta, fval] = fmincon(@(theta) maxLikelihoodEstimates2(theta,C,weight),theta_0,[],[],[],[],LB,UB,[],options);
fprintf('%d\n',theta);fprintf('.....<>%d\n',fval);
%------
num_valid_cond = length(validConds(:,1));
weight_inv = zeros(4*num_valid_cond,4*num_valid_cond);    
for i = 1:num_valid_cond
    weight_inv(4*i-3:4*i,4*i-3:4*i) = V;
end
weight = weight_inv^(-1);
sum = maxLikelihoodEstimates2(theta, validConds(:,[3,4,8,9]), weight);
fprintf('Sum1<>%d\n',sum);
%-----
Z = C-G_temp(:,[1,2,6,7]);
V = Z'*Z/(num_train_cond-1);
weight_inv = zeros(4*num_train_cond,4*num_train_cond);
for i = 1:num_train_cond
    weight_inv(4*i-3:4*i,4*i-3:4*i) = V;
end
weight = weight_inv^(-1);
[theta, fval] = fmincon(@(theta) maxLikelihoodEstimates2(theta,C,weight),theta,[],[],[],[],LB,UB,[],options);
fprintf('%d\n',theta);fprintf('.....<>%d\n',fval);

C = validConds(:,[3,4,8,9]);
num_valid_cond = length(validConds(:,1));
V=eye(4);
weight_inv = zeros(4*num_valid_cond,4*num_valid_cond);    
for i = 1:num_valid_cond
    weight_inv(4*i-3:4*i,4*i-3:4*i) = std_dev*V;
end
sum = maxLikelihoodEstimates2(theta, C, weight);
fprintf('Sum2<>%d\n',sum);
cond(20:25,:)
end

function Sum = maxLikelihoodEstimates2(theta, C, weight)
global G_temp
numCond = length(C(:,1));
for j = 1:numCond
    param = [C(j,2) C(j,1)];
    G_temp(j,:) = sing_react_oracle(param,theta);
end

for i = 1:numCond
    y_hat1(4*i-3:4*i,1) = G_temp(i,[1,2,6,7])';
    y1(4*i-3:4*i,1) = C(i,:);
end

Z = y1 - y_hat1;

Sum = Z'*weight*Z;
Sum = abs(Sum);
end

function C = sing_react_oracle(param,rate_const)
    a = rate_const(1,1); b = rate_const(1,2);
  a1 = rate_const(2,1); b1 = rate_const(2,2); a2 = rate_const(3,1); b2 = rate_const(3,2);
  a3 = rate_const(4,1); b3 = rate_const(4,2);a4 = rate_const(5,1); b4 = rate_const(5,2);
  theta_true = [a b a1 b1 a2 b2 a3 b3 a4 b4];
    C = single_reactor_function(param, theta_true);
end

function C = single_reactor_function(param, theta)
    %reactor_lumped_model
  T_1 = param(1)+ 273.15;
  t = param(2)*60;
  C = [1 1 0 0 0 1 0 0];
   
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
  k4 = A4*exp(-E4/T_1);
  k5 = A5*exp(-E5/T_1);
  options = odeset('RelTol',1e-10,'NonNegative', [1 2 4 5 6 7 8]); 
 
  [t_out,C_out] = ode15s(@reactor2_odefun,[0 t],C(1:8),options,[k1 k2 k3 k4 k5]);
  
  C = C_out(end,:);
end

function C = multi_reactor_function(param, rate_const)
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
  F = 1/t(1) + 1/t(2);
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
      [t_out,C_out] = ode15s(@reactor2_odefun,[0 t(2)],C(1:8),options,[k1 k2 k3 k4 k5]);
  end
  
  C = C_out(end,:);
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