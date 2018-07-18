function test_nelder_mead_simplex
n = 4; % number of dimensions
alpha = 1; beta = 1 + 2/n; gamma = .75 - 1/(2*n); delta = 1 - 1/n;

funValMap = containers.Map('KeyType','char','ValueType', 'double');
funValArr = zeros(1,n + 1);
%***Initialization

% First point is centroid of param space
T1_span =[25, 100];
T2_span =[25, 100];
t1_span =[60, 1000];
t2_span =[30,1000];

Box = [T1_span; T2_span; t1_span; t2_span];
x0 = [mean(T1_span),mean(T2_span), mean(t1_span), mean(t2_span)]; % initial pt

%DOE
x1 = [25, 25, 100, 70]; x2 =[100,100,1000,900];
x3 = [25, 25, 1000, 900]; x4 = [100, 25, 100, 70];
x5 = [25, 100, 100, 70]; x6 = [25, 100, 1000, 900];
x7 = [100, 25, 1000, 900]; 
doeXArr = [x1;x2;x3;x4;x5;x6;x7;x0*.8];
doeFunArr = [];
for i =1:8
   x = squeeze(doeXArr(i,:));
   
   doeFunArr(i) = functionValue(x);
end
x
doeFunArr
[x, temp] = sortVertex1(doeXArr,doeFunArr);
x0 = squeeze(x(1,:));
%x0 = x0';
a = .40; % initial simplex size
p = a/n/sqrt(2)*(sqrt(n+1)+n-1);
q = a/n/sqrt(2)*(sqrt(n+1)-1) ;
x = zeros(1,n+1,n);
tol =1e-3;
numIter =20;
% Initialize simplex
x(1,1,:) = x0;
for i = 2:n+1
    x(1,i,:) = x0;
    x(1,i,i-1) = x0(i-1)*.90;
end

for i = 1:n+1
   if x(1,i,4) - x(1,i,3) > -30
       x(1,i,4) = x(1,i,3) - 60;
   end
end


%1----get experimental data here
for i = 1:n+1
    fprintf('%d %d %d %d\n',squeeze(x(1,i,:)));
    temp = functionValue(x(1,i,:));
    funValMap(vertex_to_string(x(1,i,:))) = temp;
    funValArr(i) = temp;
end
%----

ind = 1;
%*** Sort
[x(ind,:,:), funValArr] = sortVertex(x(ind,:,:),funValArr);
x
funValArr
while(true)
    x_bar = squeeze(mean(x(ind,1:end-1,:)));
    x_n = x(ind,n,:);
    f_n = funValArr(n);
    x_n_plus_one = squeeze(x(ind,n+1,:));
    f_n_plus_one = funValArr(n+1);
    
    %*** Reflection
    x_r = x_bar + alpha*(x_bar - x_n_plus_one);
    x_r = isFeasible(x_r,x_n_plus_one,Box);
    %2---------- get experimental data here
    f_r = functionValue(x_r);
    %----------
    
    x_1 = squeeze(x(ind,1,:));
    f_1 = funValArr(1);
    x(ind+1,:,:) = x(ind,:,:);
    if (f_1 <= f_r && f_r < f_n)
        x(ind+1,n+1,:) = x_r;
        funValArr(n+1) = f_r;
        [x(ind+1,:,:),funValArr] = insertionsortVertex(squeeze(x(ind+1,:,:)), funValArr);
%         fprintf('Reflection - \n');
%         disp(x);
%         fprintf('Reflection FVal - %d %d %d %d\n',funValArr');
    elseif (f_r < f_1)
    %*** Expansion
        x_e = x_bar + beta*(x_r - x_bar);
        x_new = isFeasible(x_e,x_n_plus_one,Box);
    %3 ----- get experimental data here
        f_e = functionValue(x_new);
    %------
        if f_e < f_r
           x(ind+1,n+1,:) = x_new; 
           funValArr(n+1) = f_e;
        else
           x(ind+1,n+1,:)= x_r; 
           funValArr(n+1) = f_r;
        end
        %fprintf('Expansion FVal - %d %d %d %d %d\n',funValArr');disp(x);
        [x(ind+1,:,:),funValArr] = insertionsortVertex(squeeze(x(ind+1,:,:)), funValArr);
%         fprintf('Expansion -\n');
%         disp(x);
%         fprintf('Expansion FVal - %d %d %d %d %d\n',funValArr');
    elseif (f_n <= f_r && f_r <= f_n_plus_one)    
    %*** Outside Contraction
        x_oc = x_bar - gamma*(x_r - x_bar);
        x_new = isFeasible(x_oc,x_n_plus_one,Box);
        %4 ------ get experimental data here
        f_oc = functionValue(x_new);
        %-------
        if(f_oc <= f_r)
            x(ind+1,n+1,:) = x_new; 
            funValArr(n+1) = f_oc;
            [x(ind+1,:,:),funValArr] = insertionsortVertex(squeeze(x(ind+1,:,:)), funValArr);
%             fprintf('Contraction1 -\n');
%             disp(x);
%             fprintf('Contraction1 FVal - %d %d %d %d %d\n',funValArr');
        else
            %**** Shrink
            for i = 2:n
               x(ind+1,i,:) =  x_1 + delta*(squeeze(x(ind,i,:)) - x_1);
               funValArr(i) = functionValue(x(ind,i,:));
            end
            %*** Sort
            [x(ind+1,:,:), funValArr] = sortVertex(x(ind+1,:,:),funValArr);
%             fprintf('Shrink1 -\n');
%             disp(x);
%             fprintf('Shrink2 FVal - %d %d %d %d %d\n',funValArr');
        end
    elseif (f_r >= f_n_plus_one)
    %*** Inside Contraction
       x_ic = x_bar - gamma*(x_r - x_bar);
       x_new = isFeasible(x_ic,x_n_plus_one,Box);
       %5 ------ get experimental data here
        f_ic = functionValue(x_new);
        %-------
       if f_ic < f_n_plus_one
         x(ind+1,n+1,:) = x_new;
         funValArr(n+1) = f_ic;
         [x(ind+1,:,:),funValArr] = insertionsortVertex(squeeze(x(ind+1,:,:)), funValArr);
%          fprintf('Contraction2 -\n');
%          disp(x);
%          fprintf('Contraction2 FVal - %d %d %d %d %d\n',funValArr');
       else
       %*** Shrink
            for i = 2:n
               x(ind+1,i,:) =  x_1 + delta*(squeeze(x(ind,i,:)) - x_1);
               funValArr(i) = functionValue(x(ind,i,:));
            end
            %*** Sort
            [x(ind+1,:,:), funValArr] = sortVertex(x(ind+1,:,:),funValArr);
%             fprintf('Shrink2 -\n');
%             disp(x);
%             fprintf('Shrink2 FVal - %d %d %d %d %d\n',funValArr');
       end
    end
    ind = ind +1;
    temp = std(funValArr)/mean(funValArr);
    if  (abs(temp) < tol || ind > numIter)
        break;
    end
end

length(x(:,1,1))
x
temp
fprintf('%d %d %d %d %d\n',funValArr');
a =squeeze(x(end,:,:));a';
control_params = [5 5; 4 9;5 10; 5.5 5.2; 4 11];
competing_pathway1 = [5 5; 5.1 5.5;5 10; 5.5 5.2; 4 11];
competing_pathway2 = [5 5; 4 9;4.8 5.2; 5.5 5.2; 4 11];
prod_brkdwn = [5 5; 4 9;5 10; 5.5 5.2; 5.4 6];
rate_consts = competing_pathway2;
fprintf('best yield %d\n',functionValueYield(squeeze(a(1,:)),rate_consts));
%fprintf('best prod rate %d\n',functionValueProdRate(squeeze(a(1,:)),rate_consts));
end

function x_new = isFeasible(x_new,x_f,Box)
   feasible = true;
   num_iterations = 10;
   j = 0;
   alpha = .6;
   while(true)
       if(x_new(1) < Box(1,1) || x_new(1) > Box(1,2))
            feasible = false;
       elseif (x_new(2) < Box(2,1) || x_new(2) > Box(2,2))
            feasible = false;  
       elseif (x_new(3) < Box(3,1) || x_new(3) > Box(3,2))
            feasible = false; 
       elseif (x_new(4) < Box(4,1) || x_new(4) > Box(4,2))
            feasible = false;  
       elseif (x_new(4) >= x_new(3))
           x_new(4) = x_new(3) - 30;
           feasible = false;
       end
       if ~feasible
           if j >= num_iterations
               x_new = x_f;
               break;
           else
               x_new = x_f + alpha*(x_new - x_f);
               j = j+1;
           end
           feasible = true;
       else
           break
       end
   end
end

function f = functionValueYield(param, rate_const)
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
  options = odeset('RelTol',1e-10);%'NonNegative', [1 2 4 5 6 7 8]); 
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
      [t_out,C_out] = ode15s(@reactor2_odefun,[0 t(2)],C(1:8),options,[k1 k2 k3 k4 k5]);
  end
  
  C_min = min([1 C_6])
  if C_min < 1e-4
      C_min = 0;
  end
  if C_min <0
      C_out1=-1;
  else
      if C_6 >0
          C_temp = C_out(end,7)
          if C_temp <0
             C_temp = 0;
          end
          if C_min ==0
              C_out1 = 0;
          else
            C_out1 = C_temp/C_min;
          end
      else
          C_out1 = -1;
      end
       if C_temp*F < .001
           C_out1 = -1;
       end
  end
  
  f = -C_out1

end

function f = functionValue(param)
    control_params = [5 5; 4 9;5 10; 5.5 5.2; 4 11];
    competing_pathway1 = [5 5; 5.1 5.5;5 10; 5.5 5.2; 4 11];
    competing_pathway2 = [5 5; 4 9;4.8 5.2; 5.5 5.2; 4 11];
    prod_brkdwn = [5 5; 4 9;5 10; 5.5 5.2; 5.4 6];
    rate_consts = competing_pathway2;
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
  fprintf('%d\n',C_out(end,7));
  if C_min <=0
      C_out1=-1;
  else
      if C_6 >0
          C_temp = C_out(end,7);
          if C_temp < 0
             C_temp = 0;
          end
          C_out1 = C_temp*F;
      else
          C_out1 = -1;
      end
      if C_out(end,7)/C_min < .75
          C_out1 = 0;
      end
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

function str = vertex_to_string(x)
   str = '';
   len = length(x);
   if len > 1
       for i = 1:len-1
           str = strcat(str,int2str(x(i)),',');
       end
   end
   str = strcat(str,int2str(x(end)));
end

function [x, funVal] = insertionsortVertex(x, funVal)
    % Perform one iteration of insertion sort where the last element
    % of funVal is inserted in its proper place.
    size = length(funVal);
    f_new = funVal(end);
    x_new = x(end,:);
    for i = 1: size - 1
       if f_new < funVal(i)
          funVal(i:end) = circshift(funVal(i:end)',1);
          for j =1:size-1
              x(i:end,j) = circshift(x(i:end,j),1);
          end
          %x(i,:) = x_new;
          %funVal(i) = f_new;
          break;
       end
    end
end

function [x, fun] = sortVertex(x, fun)
% bubbleSort
    n = length(fun);
    while(true)
       swapped = false;
       for i=2: n
          if fun(i-1) > fun(i) 
              temp = fun(i-1);
              fun(i-1) = fun(i);
              fun(i) = temp;
              temp = x(1,i-1,:);
              x(1,i-1,:) = x(1,i,:);
              x(1,i,:) = temp;
              swapped = true;
          end
       end
       if ~swapped
            break;
       end
    end
end

function [x, fun] = sortVertex1(x, fun)
% bubbleSort
    n = length(fun);
    while(true)
       swapped = false;
       for i=2: n
          if fun(i-1) > fun(i) 
              temp = fun(i-1);
              fun(i-1) = fun(i);
              fun(i) = temp;
              temp = x(i-1,:);
              x(i-1,:) = x(i,:);
              x(i,:) = temp;
              swapped = true;
          end
       end
       if ~swapped
            break;
       end
    end
end