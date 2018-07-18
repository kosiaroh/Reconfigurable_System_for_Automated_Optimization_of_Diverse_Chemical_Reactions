function [pt,Conv,Grad,GradOld,d_k, Termination, ExpNum, Optimization_vars,tol, Options] = conjugateGradient(pt,Conv,GradOld, d_k, Grad, Termination, ExpNum, Box, initParams, dtParams,numTraject, Optimization_vars, Options)
DoEExpNum = Options(1);
numParams = length(dtParams);
tol = Options(2);
Index = Options(3:4);
MaxIndex = Options(5);
StepSize = Options(6:5+numParams);
Step = Options(6+numParams:5+2*numParams);
tConstraint = Options(6+2*numParams:end);

% NextExpFlag = Exp_Flags(1);
% FlowRateFlag = Exp_Flags(2);

GradStep = Optimization_vars(1);
ArmijoStep = Optimization_vars(2);
alpha = Optimization_vars(3);
armjio_contract = Optimization_vars(4);

if (numTraject > 1 && ExpNum == 1 + DoEExpNum) || (numTraject == 1 && ExpNum ==1+ DoEExpNum)
    % Calculate gradient direction after DOE
    GradOld = Grad;
%     options = optimoptions('fminunc','Algorithm','quasi-newton');
%     f  = @(weight) maxLikelihoodEstimate(weight, Conv, numParams, pt);
%     weight = zeros(numParams,1);
%     weight = fminunc(f,weight,options) ;
%     
%     [temp, Grad] =  predictedFunctionValve(weight, initParams, numParams);
    if numTraject > 1
        Grad = returnDOEGrad(Conv,pt(2:end,:), length(pt(1,:))-1, dtParams);
    else
        Grad = returnDOEGrad(Conv,pt, length(pt(1,:))-1, dtParams);
    end
    
    beta = norm(Grad)^2/norm(GradOld)^2;
    if numTraject == 1
       d_k =  Grad'; 
    else
      d_k =  Grad' + beta*d_k;
    end
end

if ExpNum > length(pt(:,1))
    [StepSize, tConstraint, params] = getGradStep(Step, d_k, 1, initParams, Box, GradStep);
    feasible = isFeasible(params,Box);
    
    % if initial proposed step is not feasible according to box constraints
    if ~feasible
        for i = 1:10
            alpha = alpha*.5;
            [StepSize, tConstraint, params] = getGradStep(Step, d_k, alpha, initParams, Box, GradStep);
            feasible = isFeasible(params,Box);
            if feasible
                pt = [pt; ExpNum, params];
                GradStep = GradStep + 1;
                break;
            end
        end 
        %Need to tell driver to calculate new DOE or terminate since we
        %reached the border of constraint(s)
        if ~feasible
            armjio_contract = true;
        end
    else 
        pt = [pt; ExpNum, params];
        GradStep = GradStep + 1;
    end
    
end


Optimization_vars(1) = GradStep;
Optimization_vars(2) = ArmijoStep;
Optimization_vars(3) = alpha;
Optimization_vars(4) = armjio_contract ;
Options = [StepSize, tConstraint];

end

function params  = checkResidenceTimeCondition(params,Box)
  % Assume Params are ordered T1 T2 ... t1 t2 ...
  numResTime = length(params)/2;
  for i = numResTime+2:length(params)
      if params(i-1) - params(i) < 30
          if params(i) - 30 > Box(i,1)
              params(i) = params(i) - 30;
          elseif params(i-1) + 30 < Box(i-1,2)
              params(i-1) = params(i-1) + 30;
          end
      end
  end
end

function [StepSize, tConstraint, params] = getGradStep(Step, d_k, alpha, initParams, Box, GradStep)
   StepSize = Step./norm(d_k);
   numParam = length(initParams);
   tConstraint = zeros(1, numParam);
   params = zeros(1, numParam);
   for i = 1:numParam
       if Box(i,1) == initParams(i) || Box(i,2) == initParams(i) 
           StepSize(i) = Step(i)/sqrt(d_k(i)^2);
           break;
       end
   end
   
   GradStep = GradStep+1;
   for i = 1:length(initParams)
       if initParams(i)+ alpha*GradStep*d_k(i)*StepSize(i) < Box(i,1) 
           params(i) = Box(i,1);
           tConstraint(i) = 1;
       elseif initParams(i)+ alpha*GradStep*d_k(i)*StepSize(i) > Box(i,2)
           params(i) = Box(i,2);
           tConstraint(i) = 1;
       else
           params(i) = round(initParams(i)+ alpha*GradStep*d_k(i)*StepSize(i));
       end
   end
   params = checkResidenceTimeCondition(params,Box);
end

function feasible = isFeasible(x,Box)
   feasible = true;
   
   for i = 1:length(x)
       if(x(i) < Box(i,1) || x(i) > Box(i,2))
           feasible = false;
           break;
       end
   end
       
   if length(x) > 2
       for i =length(x)/2+2:length(x)
           if (x(i) >= x(i-1))
               x(i) = x(i-1) - 30;
               feasible = false;
               break;
           end
       end
   end
end

% function [f, g] =  predictedFunctionValve(weight, params,numParams)
%     f = 0;
%     g = 0;
%     %
%     weightIndex = 1;
%     for i=1:numParams
%        f = f + weight(weightIndex)*params(i);
%        g(i) = weight(weightIndex);
%        weightIndex = weightIndex + 1;
% %        for j=1:numParams
% %            if i ~= j
% %               f = f + weight(weightIndex)*params(i)*params(j); 
% %               g(i) = g(i) + weight(weightIndex)*params(j);
% %               weightIndex = weightIndex + 1;
% %            end
% %        end
%        %^(-1);% + weight(2*(i-1)+3)*(params(i)^(-1))^2;
%        %*params(i);%^(-2);% -2*weight(2*(i-1)+3)*params(i)^(-3);
%     end
% end

% function Sum = maxLikelihoodEstimate(weight, Conv, numParams, params)
%     numCond = length(Conv(:));
%     Conv_pred = [];
%     for j = 1:numCond
%         [Conv_pred(j), temp] = predictedFunctionValve(weight, params(j,2:end),numParams);
%     end
% 
%     Z = Conv_pred' - Conv';
% 
%     Sum = Z'*Z;
% end
% 
% function params = getNewContractedParams(pt, Index)
%    paramSize = length(pt(1,2:end));
%    params = zeros(1,paramSize);
%    for i=1:length(pt(1,2:end))
%        params(i) = round(mean([pt(Index(1),i+1),pt(Index(2),i+1)]));
%    end
% end