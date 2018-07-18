%clear all;
T1_span =[25, 100];
T2_span =[25, 100];
T3_span =[25, 100];
t1_span =[1, 16.667];
t2_span =[.5,16.667];
t3_span =[.25,16.667];
Grad =[.001; .001; .001; .001; .001; .001];
GradOld = [.001; .001; .001; .001; .001; .001];
d_k = [.001; .001; .001; .001; .001; .001];
DoEStartPts = [];
DoEStartConv = [];
newTrajectStart = [];
n=6;
Box = [T1_span; T2_span; T3_span; t1_span; t2_span; t3_span];
Box (4:6,:) = Box(4:6,:).*60; %convert minutes to seconds
all_conv = [];
Step_size =[];
indData = 1;
Vr1 = 300;
Vr2 = 300;
Vr3 = 300;
newFunVal = 0;
new_x =0;
StepSize = [25 25 25 240 240 240];
%%--------Initialization
 weight = zeros(n,1);
 %funValArr = zeros(1,n + 1);
 x0 = [mean(Box(1,:)),mean(Box(2,:)),mean(Box(3,:)),mean(Box(4,:)),mean(Box(5,:)),mean(Box(6,:))];
 dT1 = round(.08*mean(Box(1,:)));
 dT2 = round(.08*mean(Box(2,:)));
 dT3 = round(.08*mean(Box(3,:)));
 dt1 = round(.08*mean(Box(4,:)));
 dt2 = round(.08*mean(Box(5,:)));
 dt3 = round(.08*mean(Box(6,:)));
 GradStep = 0;
 Contract = 0;
 ContractAccept = 0;
 ArmijoContract = 0;
 ArmijoAccept = 0;
 ArmijoStep = 1/2;
%  DoET = [InitT-TSpan(3),InitT-TSpan(3),InitT+TSpan(3),InitT+TSpan(3)];
%  DoEtau = [Inittau-tauSpan(3),Inittau+tauSpan(3),Inittau-tauSpan(3),Inittau+tauSpan(3)];
%  DoEratio = [1,1,1,1];
%  [YRand,IRand] = sort(rand(1,DoEExpNum));

%-------------------


old_conv = [];
all_pts = [];
all_d_k = [];
Step = [25 25 25 240 240 240];
dtparams = [dT1 dT2 dT3 dt1 dt2 dt3];
for j =1:5
 ExpNum = 1;
 GradStep = 0;
 alpha = 1;
 armjio_contract = false;
 Conv = [];
 pt = [];
 if j > 1
    dT1 = round(.05*x0(1));
    dT2 = round(.05*x0(2));
    dT3 = round(.05*x0(3));
    dt1 = round(.05*x0(4));
    dt2 = round(.05*x0(5));
    dt3 = round(.05*x0(6));
    dtparams = [dT1 dT2 dT3 dt1 dt2 dt3];
    for lm = 1:length(dtparams)
      if x0(lm) - dtparams(lm) <   Box(lm,1)
        x0(lm) =  Box(lm,1) + dtparams(lm);
      elseif  x0(lm) + dtparams(lm) >   Box(lm,2)
        x0(lm) =  Box(lm,2) - dtparams(lm);
      end
    end
    ExpNum = 2;
    x = returnDOE(n,dtparams, x0);
    DoEExpNum = length(x(:,1));
    DoEMatrix = [2:1:DoEExpNum+1]';
    Step = [25 25 25 240 240 240];
    pt = [1 x0(1) x0(2) x0(3) x0(4) x0(5) x0(6)];
    Conv(1) = old_conv(end);
 else
   x = returnDOE(n,dtparams, x0);
   DoEExpNum = length(x(:,1));
   DoEMatrix = [1:1:DoEExpNum]';
 end
 
 
 

 MaxExp = DoEExpNum + 30;
 MaxIndex = DoEExpNum + 2;
 Index = [DoEExpNum+1,DoEExpNum+2];
 %Box = [T1_span; T2_span; t1_span; t2_span];
 
 pt = [pt;DoEMatrix,x(:,1), x(:,2),x(:,3), x(:,4),x(:,5), x(:,6)];
 if j == 1
    pt(end+1,1:7)= [DoEExpNum+1 x0(1) x0(2) x0(3) x0(4) x0(5) x0(6)]; 
 end
 
 tol = 1e-7;
 Termination = 0;
 MinStep = 1;
 tConstraint = zeros(1,length(dtparams));


while(not(Termination))
    if ExpNum > DoEExpNum+5

        % Check for accepting Armijo step
        currStepSize = abs(d_k(2)*StepSize(2)*ArmijoStep^(GradStep+1));
        
        % Termination criteria
        if ExpNum > MaxExp || ...
                (Conv(end) <= 0.95*max(Conv(DoEExpNum+1:end))) || ...
                (checkConstraintsTerminate(tConstraint, d_k, StepSize, MinStep)) || ...
                abs(Conv(end) - max(Conv(DoEExpNum+1:end-1)) < tol)
            Termination=1;
        end
    end
    
    if (j > 1 && ExpNum == 1 + DoEExpNum) || (j == 1 && ExpNum == 1+DoEExpNum)
         [y, temp2] = max(Conv);
         x0 = pt(temp2,2:length(dtparams)+1);
         DoEStartPts = [DoEStartPts;x0];
         DoEStartConv = [DoEStartConv;y];
    end
    
    Optimization_vars(1) = GradStep;
    Optimization_vars(2) = ArmijoStep;
    Optimization_vars(3) = alpha;
    Optimization_vars(4) = armjio_contract ;
    Step_size(end+1:end+n) = StepSize;
    [pt,Conv,Grad,GradOld,d_k, Termination, ExpNum, Optimization_vars, tol,Options] = conjugateGradient(pt,Conv,GradOld, d_k, Grad, Termination, ExpNum, Box, x0, dtparams,j, Optimization_vars, [DoEExpNum  tol Index MaxIndex StepSize Step tConstraint ]);
    GradStep = Optimization_vars(1);
    ArmijoStep = Optimization_vars(2);
    alpha = Optimization_vars(3);
    armjio_contract = Optimization_vars(4);
    StepSize = Options(1:end/2);
    tConstraint = Options(end/2+1:end);
    indData = indData + 1;
    
    if armjio_contract
        break;
    end
    
    new_x = [pt(ExpNum, 2), pt(ExpNum, 3), pt(ExpNum, 4), pt(ExpNum, 5), pt(ExpNum, 6), pt(ExpNum, 7)];
    T1 = new_x(1); %'C
    T2 = new_x(2);  %'C
    T3 = new_x(3); 
    tau1 = new_x(4);
    tau2 = new_x(5);
    tau3 = new_x(6);
    P1 = 2/3*Vr1/(tau1*60); %ul/min
    P2 = 1/3*Vr1/(tau1*60); %ul/min
    Ftot = Vr2/(60*tau2);
    P3 = Ftot - Vr1/(tau1*60); 

    %---------  

   if validResidenceTimes(new_x)
     averagePeaks = functionValue_simplex3(new_x);
   else
     averagePeaks = 0; 
   end

%%--------Second part of loop

   %Only tracking main product production rate
   Conv(ExpNum) = -averagePeaks;
   ExpNum = ExpNum + 1;
   if norm(d_k) == 0
       break;
   end
end
  [q temp_i] = max(Conv);
  x0 = pt(temp_i,2:length(dtparams)+1);
  newTrajectStart = [newTrajectStart;x0];
  old_conv(end+1) = q;
  num_pts = length(pt(:,1));
  all_pts = [all_pts;pt];
  all_d_k = [all_d_k;d_k];
  if ~isempty(all_conv)
      len = length(all_conv(1,:));
      num = len - length(Conv);
      b = padarray(Conv',num,0,'post');
      Conv = b';
  else
      b = padarray(Conv',10,0,'post');
      Conv = b';
  end
  all_conv = [all_conv; Conv];
  if (j >= 2 && (abs(old_conv(end) - old_conv(end-1)) < tol)) || norm(d_k) == 0
    break;
  end
end