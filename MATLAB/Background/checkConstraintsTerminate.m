function terminate = checkConstraintsTerminate(tConstraints, cg, StepSize, MinStep)
    terminate = true;
    for i = 1: length(tConstraints)
       terminate = terminate && (tConstraints(i) || abs(cg(i)*StepSize(i)) < MinStep);
    end
end