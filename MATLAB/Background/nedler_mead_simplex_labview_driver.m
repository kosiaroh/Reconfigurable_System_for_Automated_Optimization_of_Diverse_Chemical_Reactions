T1_span =[25, 100];
T2_span =[25, 100];
T3_span =[25, 100];
t1_span =[1, 16.667];
t2_span =[.5,16.667];
t3_span =[.25,16.667];
n=6;
Box = [T1_span; T2_span; T3_span; t1_span; t2_span; t3_span];
maxIterations = 40;
all_data = zeros(maxIterations,n+1,n);
curr_step_arr = zeros(1,maxIterations);
funValArr_backup = zeros(maxIterations,n+1);
indData = 1;
Vr1 = 300;
Vr2 = 300;
newFunVal =0;
new_x =0;

%%--------Initialization
 Box (4:6,:) = Box(4:6,:).*60; %convert minutes to seconds
  n = 6;
  funValArr = zeros(1,n + 1);
  x0 = [mean(Box(1,:)),mean(Box(2,:)),mean(Box(3,:)),mean(Box(4,:)),mean(Box(5,:)),mean(Box(6,:))];
  tol = 1e-3;
  all_verticies = zeros(n+1,n);
  all_verticies (1,:) = x0;
  for i = 2:n+1
      all_verticies (i,:) = x0;
      all_verticies (i,i-1) = x0(i-1)*.80;
  end

  for i = 1:n+1
     if all_verticies (i,5) - all_verticies (i,4) > -30
        all_verticies (i,5) = all_verticies (i,4) - 60;
     end
    if all_verticies (i,6) - all_verticies (i,5) > -30
        all_verticies (i,6) = all_verticies (i,5) - 60;
    end
  end
  initial_simplex = true;
  doe_ind = 0;
  curr_step = 0;
  itrNum = 0;
  shrink = false;
  shrinkInd = 1;
%-------------------


while(itrNum <= maxIterations)
    
%%------First Part of loop
if initial_simplex
    doe_ind  = doe_ind + 1;
    new_x =  squeeze(all_verticies(doe_ind ,:));
elseif  (shrink && shrinkInd < n + 1)
   shrinkInd = shrinkInd + 1; 
   new_x =  squeeze(all_verticies(shrinkInd ,:));
else
    if curr_step == 0
      x_r = 0;
      f_r = 0;
    elseif curr_step == 1
       curr_step = 11;
    elseif curr_step == 2
        curr_step = 22;
    elseif curr_step == 3
        curr_step = 33;
    elseif curr_step == 4
        curr_step = 44;
    elseif curr_step == 5
        curr_step = 55;
    end
    %----Debugging purposes
    all_data(indData,:,:) = all_verticies;
    funValArr_backup(indData,:) = funValArr;
    %---
    %fprintf('%d %d %d %d %d\n',funValArr');
    [all_verticies, new_x, curr_step,f_r,x_r,funValArr] = nedler_mead_simplex(all_verticies,curr_step, funValArr,newFunVal,new_x,x_r,f_r, Box);
    curr_step_arr(indData) = curr_step;
    indData = indData + 1;
    fprintf('%d %d %d %d %d\n',funValArr');
    if curr_step == 6
       shrink = true;
    end

    itrNum = itrNum + 1;
end

T1 = new_x(1); %'C
T2 = new_x(2);  %'C
tau1 = new_x(3);
tau2 = new_x(4);
P1 = 2/3*Vr1/tau1*60; %ul/min
P2 = 1/3*Vr1/tau1*60; %ul/min
Ftot = Vr2/tau2;
P3 = Ftot - Vr1/tau1; 
%---------  
averagePeaks = functionValue_simplex3(new_x);

%%--------Second part of loop


% lastrowStr = num2str(lastrow2);
% temp_val1 = [itrNum tau1 tau2 T1 T2  Temp1 Temp2 P1 P2 P3 averagePeaks curr_step];
% xlSheets.Item(idx_con).Range(strcat('A',lastrowStr,':',lastcol2,lastrowStr)).Value2 = temp_val1;
% lastrow2 = lastrow2+ 1;

if  initial_simplex
   %Only tracking main product production rate
   funValArr(doe_ind) = averagePeaks;
   if doe_ind == n+1
       initial_simplex = false;
       curr_step =0;
       [all_verticies,funValArr] = sortVertex_simplex(all_verticies,funValArr);
    end
elseif  (shrink && shrinkInd <= n + 1)
  %Only tracking main product production rate
   funValArr(shrinkInd ) = averagePeaks;
   if shrinkInd == n+1
       shrink = false;
       shrinkInd = 1;
       curr_step =0;
    end
else
   newFunVal = averagePeaks;
end
end