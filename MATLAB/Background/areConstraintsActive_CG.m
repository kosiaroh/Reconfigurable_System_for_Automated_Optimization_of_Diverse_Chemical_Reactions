function trueOrFalse = areConstraintsActive_CG(param,Box, cg,GradStep, StepSize, ArmijoStep, MinStep)
   trueOrFalse = false;
   sum = true;
   for i = 1:length(param)
       onBorder = (param(i) == Box(i,1) || param(i) == Box(i,2));
       paramStepSize = abs(cg(i)*StepSize*ArmijoStep^(GradStep+1));
       sum = sum & onBorder;
       if  onBorder && paramStepSize < MinStep
           trueOrFalse = true;
           break;
       end
   end
   
   if  not(trueOrFalse)
       trueOrFalse = sum;
   end
end