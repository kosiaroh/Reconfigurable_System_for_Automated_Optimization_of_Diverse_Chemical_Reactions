function trueOrFalse = validResidenceTimes(params)
  trueOrFalse = true;
  numResTime = length(params)/2;
  for i = numResTime+2:length(params)
      if params(i-1) - params(i) < 10
          trueOrFalse = false;
          break;
      end
  end
end