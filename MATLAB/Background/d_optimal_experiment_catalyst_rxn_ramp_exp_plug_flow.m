function d_optimal_experiment_catalyst_rxn_ramp_exp_plug_flow
   % [log(A1*exp(-E1/RT*)) log(E1/R) log(A2*exp(-E2/RT*)) log(E2/R) log(A3*exp(-E3/RT*)) log(E3/R) log(A4*exp(-E4/RT*)) log(E4/R)]
   % T* = 293 K; R = 8.314 J/mol K
   % Original log10(A) and E
%    T_star = 343;
%    b0_true = [0.67 9.6 0.38  13.2  5.05  59.8 2.99  42.0];
   %theta1_0_true = [log(10^b0_true(1)) b0_true(2)/(0.008314*T_star) log(10^b0_true(3)) b0_true(4)/(0.008314*T_star) log(10^b0_true(5)) b0_true(6)/(0.008314*T_star) log(10^b0_true(7)) b0_true(8)/(0.008314*T_star)];
   %theta1_0_true = [ 11.1563 40 ];%3.783;3.783;];
   theta1_0_true = [log(2*10^4), log(10000), log(10000), log(1.42*10^2), 40, 80];
  theta_0 = [1 1 1 1 50 80];%3 3];
   
  %exp_space = [type_of_tau upper_tau_bound(min) lower_tau_bound(min) frequency temp_type temp_start];
  %exp_time(min) = frequency/100
  % Define initial set of experiments 
  initial_ramp_design = [1 43 36 1000 1 110];
  
  [tau T_fn] = expSineConds(initial_ramp_design);
  
%   Param_bound - [Lower_bound_tau Upper_bound_tau dt_tau...
%                  Lower_bound_freq Upper_bound_freq dt_freq...
%                  Lower_temp_start Upper_temp_start dt_temp_start]
  Param_bound = [1 64 7; ...
                 1000 10000 3000; ...
                 50 170 60;];
  %% Initialize
  start1 = tic;
  Param_space=enumeration_grid(Param_bound);
  %Param_space(end+1,:) = initial_ramp_design;
  %Param_space{end+1,1} = {tau
  Param_space_a = {};
  for i =1:size(Param_space(:,1))
    [tau1, T_fn1] = expSineConds(Param_space(i,:));
    Param_space_a{i,1} = tau1; Param_space_a{i,2} = T_fn1;
  end
  ind = [];
%   a = 1000;
%   for i = 1:5
%       ind(end+1) =findsubmat(Param_space,[1 21 21 a 1 100]);
%       a = a+2000;
%   end
%   index = true(1, size(Param_space, 1));
%   index(ind(:,1)) = false;
%   Param_space = Param_space(index, :);
%   Param_space_a = Param_space_a(index,:);
  
  stop1 = toc(start1);
  
  % Simulate concentration profiles for initial experiments
  B = {};
  B{1} = getConc(theta1_0_true,tau,T_fn);
  
  [rows_B,cols_B]=size(B{1});
   
  % Introduce experimental error as normal distribution
  sig1 = 8.6E-4; % sigma^2 for response 1
  sig2 = 8.6E-4; % sigma^2 for response 2
  sig3 = 8.6E-4; % sigma^2 for response 3
  sqrt_sig1 = sqrt(sig1);
  sqrt_sig2 = sqrt(sig2);
  sqrt_sig3 = sqrt(sig3);
  r = randn(rows_B,cols_B);
%   S = [sqrt_sig1 0 0; 0 sqrt_sig2 0; 0 0 sqrt_sig3];
   B{1} = B{1} + 0.03*r.*B{1};
   
  % Least squares determination of model fit
  Aineq = [];
  Bineq = [];
  Aeq = [];
  Beq = [];
  nonlcon = [];
   
  % Unscaled bounds
  %theta1_0_true = [log(1.25), log(7.29e-3), log(8.80e-4), log(217), -18.3, 29.44, 71.014, 23.02];
  %theta_0 = [1 -1 -4 1 -10 15 35 15];
   LB = [-15 -15 -15 -15 10 10];
   UB = [ 10  10  10  10 120 200];
   
  V = eye(1);
  weight_inv = zeros(1,1);
  for i = 1:1
     weight_inv(i,i) = V;
  end
 
  weight = weight_inv^(-1);
   
  % Number of models
  n = 1;
   
  % Current Parameter estimates after initial experiments
  clear theta
  options = optimset('MaxIter',10000,'MaxFunEvals',100,'FunValCheck','on','TolX',1e-40,'TolFun',1e-40);
  
  Exp_param = {};
  Exp_param{1,1} = tau ;
  Exp_param{1,2} = T_fn;
  theta_arr = [];
  start = tic;
  theta= fmincon(@LSS_SN2,theta_0,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,options,Exp_param,B,weight);
  theta_arr(end+1,:) = theta;
  %LSS_SN2(k_param,Exp_param,B,weight) 
  stop2 = toc(start);
  % Four initial experiments from factorial design
  expt = 1;
  prevExp = initial_ramp_design;
  numInit =1;
  
  % delete orignial experiments from param space
  ind = zeros(numInit);
  for i = 1:numInit
     ind(i) = findsubmat(Param_space,prevExp(i,:));
     %removeReplicates();
  end
%   a = 1000;
%   for i = 1:5
%       ind(end+1) =findsubmat(Param_space,[1 21 21 a 1 100]);
%       a = a+2000;
%   end
  index = true(1, size(Param_space, 1));
  index(ind(:,1)) = false;
  Param_space = Param_space(index, :);
  Param_space_a = Param_space_a(index,:);
  
  clear results;
  results=zeros(expt);   

  %% Parameter Estimation
  numParm = length(theta);
  clear CI_error;
  while true
      fprintf('%d\n',expt);
      fprintf('%d %d %d %d\n',theta);
      V = calc_cov(theta,Exp_param,B);
      V_theta_inv = zeros(numParm,numParm);
      V_inv = V^(-1);
      
      for j=1:expt
          X1 = [];
          temp = getConc_with_sens(theta,Exp_param{j,1},Exp_param{j,2});
          for kl = 1:size(temp(:,1,1))
           X1(1:3,1:6) = temp(kl,:,:);
           V_theta_inv = V_theta_inv + X1'*V_inv*X1;
          end
      end
      
      % Find confidence interval on each parameter
      F_value = chi2inv(.99,numParm); % 99% confidence interval, fitting 4 parameters to N experiments
      
      cov = V_theta_inv^(-1);
      stop3 = toc(start1);
      CI = zeros(2,2);
%       for i = 1:numParm
%           CI(i,1) = theta(i) - sqrt(F_value*cov(i,i)); % 8 params
%           CI(i,2) = theta(i) + sqrt(F_value*cov(i,i)); % 8 params
%           CI_error(i) = abs((CI(i,2) - theta(i))/theta(i))*100; % Relative error in percent
%       end
%       
%       max_error = max(CI_error(:));
%       if max_error < 5 || expt > 33
%           break;
%       end
%e
      [M_row N_col]=size(Param_space);
      %Size of remaining parameter space
      start2 = tic;
      clear Arg
      for k=1:M_row
          X2 = [];
          temp = getConc_with_sens(theta,Param_space_a{k,1},Param_space_a{k,2});
          Vtemp = zeros(numParm,numParm);
          for kl = 1:size(temp(:,1,1))
           X2(1:3,1:6) = temp(kl,:,:);
           Vtemp = Vtemp + X2'*V_inv*X2;
          end
          num_exp = size(Param_space_a{k,1});
          M1_temp = V_theta_inv + Vtemp;
          Arg(k)=real(log(det(M1_temp)/num_exp(1)));
          stop2 = toc(start2);
          fprintf('%d %d\n',k,stop2);
      end
      
      
      N_1_exp = find(Arg==max(Arg));
      if length(N_1_exp) > 1
          N_1_exp = N_1_exp(1);
      end
      fprintf('Arg %d', Arg(N_1_exp));
      %Simulate next experiment
      expt = expt + 1;
      %results(expt)=min(Arg);
      prevExp(expt,:) = Param_space(N_1_exp,:);
      Exp_param(expt,:) = Param_space_a(N_1_exp,:);
      tau = Exp_param{expt,1}; T_fn = Exp_param{expt,2};
      B{expt} = getConc(theta1_0_true,tau,T_fn);
      
      [rows_B,cols_B]=size(B{expt});
      r = randn(rows_B,cols_B);
      B{expt} = B{expt} + .03*r.*B{expt}; %Add normally distributed error to 'experimental' data
      
      % currently weight is left as the identity ****(COULD UPDATE)****
      weight_inv = ones(expt,expt);
%       for i = 1:expt
%           weight_inv(i,i) = V;
%       end
      
      LB_min = [0 0 ];
      UB_max = [9 9];
      %LB = max([CI(:,1)';LB_min]);
      %UB = min([CI(:,2)';UB_max]);
      LB = [-15 -15 -15 -15 10 10];
      UB = [ 10  10  10  10 120 200];
      weight = weight_inv^(-1);
      start1 = tic;
      % Calculate new optimized theta parameters
      theta_0 = theta;
      theta= fmincon(@LSS_SN2,theta_0,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon,options,Exp_param,B,weight);
      theta_arr(end+1,:) = theta;
      % remove experiment from param space
      index = true(1, size(Param_space, 1));
      index(N_1_exp) = false;
      Param_space = Param_space(index, :);
      Param_space_a = Param_space_a(index,:);
  end
  fprintf('\n');
  fprintf('%d %d %d %d\n',theta);
  fprintf('\n');
end

function [tau_arr T_fun] = expSineConds(ramp_design)
tau_arr = [];
index = 1;
L = ramp_design(4);
tau_max = ramp_design(2);
tau_min = ramp_design(3);
tau_type = ramp_design(1);
rnd_const = 800000;
Temp_type = ramp_design(5);
T_start = ramp_design(6);
exp_time = 2*pi*rnd_const/L;
for t = 0:30:exp_time
    if tau_type == 1
        tau = sin(t*L/rnd_const);
    elseif tau_type ==2
        tau = exp(t./L).*sin(2*pi/(30*60)*L*(exp(t./L)-1))/exp(1/L);
    else
        tau = exp(-t./L).*sin(2*pi/(30*60)*L*(exp(t./L)-1))/exp(-1/L);
    end
    
    tau_arr(index) = tau;
    index = index + 1;
end

tau_fun_min = min(tau_arr);
if tau_fun_min < 0
    tau_arr = tau_arr+ abs(tau_fun_min);
end

tau_fun_max = max(tau_arr);
tau_arr = tau_arr./tau_fun_max*(tau_max-tau_min);

tau_fun_min = min(tau_arr);
if tau_fun_min < tau_min
    tau_arr = tau_arr+ abs(tau_min - tau_fun_min);
end

index = tau_arr > tau_max;
tau_arr(index) = tau_max;

index = tau_arr < tau_min;
tau_arr(index) = tau_min;

%function handle of linear temperature ramp with t in seconds
if Temp_type == 1
    T_fun = @(t) T_start + t./60;
elseif Temp_type == 2
    T_fun = @(t) T_start - t./60;
end
end

function [V] = calc_cov(k_param,Exp_param,B)
[N_expts,N_param_cond] = size(B{1}); %determine number of experiments
[M,N]=size(Exp_param);
V = [];
% Divide by the degrees of freedom to return the response covariance matrix
for i = 1:M
    tau_arr = Exp_param{i,1}; T_fn = Exp_param{i,2};
    mod_pred = getConc(k_param,tau_arr,T_fn);
    b = B{i};
    for j = 1:3
        for k = 1:3
            Z1 = b(:,k) - mod_pred(:,k);
            Z2 = b(:,j) - mod_pred(:,j);
            Z = Z1.*Z2;
            V(j,k) = sum(Z);
        end
    end
end

%V = LSS_SN2(k_param,Exp_param,B,[]);%/(N_expts - 2);
end

function sse = LSS_SN2(k_param,Exp_param,B,weight)

[M,N]=size(Exp_param); %determine number of experiments

%Determine model predWictions for each experiment
mod_pred = {};
sse = 0;
for i = 1:M
    tau_arr = Exp_param{i,1}; T_fn = Exp_param{i,2};
    mod_pred = getConc(k_param,tau_arr,T_fn);
    b = B{i};
    Z = mod_pred - b;
    sse = sse + sumsqr(Z);
end
end

function Dis_pts=enumeration_grid(Param_bound)
% Generates enumerated matrix of experimental conditions
% Input
%   Param_bound - [Lower_bound_tau Upper_bound_tau dt_tau...
%                  Lower_bound_exp_time upper_bound_exp_time dt_exp_time...
%                  Lower_bound_freq Upper_bound_freq dt_freq...
%                  Lower_temp_start Upper_temp_start dt_temp_start]

% Output
%   Dis_pts is the matrix of possible experimental conditions

%exp_space = [type_of_tau upper_tau_bound lower_tau_bound exp_time frequency temp_type temp_start];
Dis_pts = [];

q=1;

for i=1:1
    %type of tau
    for j=1:1
        %type of temperature
        for k=Param_bound(1,1):Param_bound(1,3):Param_bound(1,2)
            %upper bound of tau
            for l = Param_bound(1,1):Param_bound(1,3):k
                %lower bound of tau
                %for m = Param_bound(2,1):Param_bound(2,3):Param_bound(2,2)
                    %exp time
                    if l ~= k
                        for n = Param_bound(2,1):Param_bound(2,3):Param_bound(2,2)
                            %frequency
                            for o = Param_bound(3,1):Param_bound(3,3):Param_bound(3,2)
                                %temperature start
                                Dis_pts(q,:)=[i k l n j o];
                                q=q+1;
                            end
                        end
                    end
                %end
            end
        end
    end
end
end

function C = getConc(params,tau_arr,T_fn)
 C  = [];
 Cequiv_all = [];
 for j = 1:length(tau_arr)
     tau = tau_arr(j);
     options = odeset('RelTol',1e-10,'NonNegative', [1 2 3 4 5 6]); 
     Ca = .6; Cb = .6;
     [t_out,C_out] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params,T_fn,j);
     C(end+1,:) = [C_out(end,1) C_out(end,2) C_out(end,3)];  
 end
 
end

function C = getConc_with_sens(params,tau_arr,T_fn)
 C  = [];
 Cequiv_all = [];
 for j = 1:length(tau_arr)
    tau = tau_arr(j);
    options = odeset('RelTol',1e-10,'NonNegative', [1 2 3 4 5 6]); 
    Ca = .6; Cb = .6;
    [t_out,C_out] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params,T_fn,j);
    delta = .03;
    params1 = params; params1(1) = params(1)*(1 + delta/100);
    params2 = params; params2(2) = params(2)*(1 + delta/100);
    params3 = params; params3(3) = params(3)*(1 + delta/100);
    params4 = params; params4(4) = params(4)*(1 + delta/100);
    params5 = params; params5(5) = params(5)*(1 + delta/100);
    params6 = params; params6(6) = params(6)*(1 + delta/100);
%    params7 = params; params7(7) = params(7)*(1 + delta/100);
%     params8 = params; params8(8) = params(8)*(1 + delta/100);
    [t_out,C_out1] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params1,T_fn,j);
    [t_out,C_out2] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params2,T_fn,j);
    [t_out,C_out3] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params3,T_fn,j);
    [t_out,C_out4] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params4,T_fn,j);
    [t_out,C_out5] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params5,T_fn,j);
    [t_out,C_out6] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params6,T_fn,j);
%     [t_out,C_out7] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params7,T_fn,j);
%     [t_out,C_out8] = ode15s(@reactor_temp_down1,[0 tau],[Ca Cb 0 .06 0 0],options,params8,T_fn,j);
    dCdTheta1 = (C_out1(end,1:3) - C_out(end,1:3))/(params(1)*delta/100);
    dCdTheta2 = (C_out2(end,1:3) - C_out(end,1:3))/(params(2)*delta/100);
    dCdTheta3 = (C_out3(end,1:3) - C_out(end,1:3))/(params(3)*delta/100);
    dCdTheta4 = (C_out4(end,1:3) - C_out(end,1:3))/(params(4)*delta/100);
    dCdTheta5 = (C_out5(end,1:3) - C_out(end,1:3))/(params(5)*delta/100);
    dCdTheta6 = (C_out6(end,1:3) - C_out(end,1:3))/(params(6)*delta/100);
%     dCdTheta7 = (C_out7(end,3:5) - C_out(end,3:5))/(params(7)*delta/100);
%     dCdTheta8 = (C_out8(end,3:5) - C_out(end,3:5))/(params(8)*delta/100);
    C(end+1,:,:) = [dCdTheta1' dCdTheta2' dCdTheta3' dCdTheta4' dCdTheta5' dCdTheta6'];% dCdTheta7' dCdTheta8']; 
 end
end

function dC = reactor_temp_down1_temp(t,Y,theta,T_fn,j)
  T = T_fn(t+j*30);
  %T = T_fn;
  if T< 40
      T = 40;
  elseif T> 200
      T = 200;
  end

  T = T + 273.15;
  A1 = exp(theta(1));
  E1 = theta(2)*1000; 
%   A2 = exp(theta(2));
%   E2 = theta(4)*1000;
  k1 = A1*exp(-E1/T/8.314);
%   k2 = A2*exp(-E2/T/8.314);
  k = [k1];
  dC1 = -k(1)*Y(1)*Y(2);
  dC2 = -k(1)*Y(2)*Y(1);
%  dC3 = k(1)*Y(2)^.77*Y(1)^.77;
%  dC3 = k(1)*Y(2)*Y(1)-k(2)*Y(3);
%  dC4 = k(2)*Y(3);
%  dC = [dC1 dC2 dC3 dC4]';
  dC = [dC1 dC2]';
end

function dC = reactor_temp_down1(t,Y,theta,T_fn,j)
T = T_fn(t+j*30);
  %T = T_fn;
  if T< 50
      T = 50;
  elseif T> 200
      T = 200;
  end

  T = T + 273.15;

  k1 = exp(theta(1))*exp(-theta(5)*1000/8.314/T);
  k2 = exp(theta(2));
  k3 = exp(theta(3));
  k4 = exp(theta(4))*exp(-theta(6)*1000/8.314/T);
  %k3 = k2/kc2;
  
  %[ [1] [2] [3] [NiL] [SMNiL] {[SM*][NiL]{Menthol]} ]
  
  r1 = k1*Y(1)*Y(4);
  r2 = k2*Y(2)*Y(5);
  r3 = k3*Y(6);
  r4 = k4*Y(4);
  
  dC1 = -r1;
  dC2 = -r2;
  dC3 = r3;
  dCNiL = -r4+r3-r1;
  dCSMNiL = -r2+r1;
  dCSMNiLMen = -r3+r2;
  dC = [dC1 dC2 dC3 dCNiL dCSMNiL dCSMNiLMen]';
end