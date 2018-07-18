function grad = returnDOEGrad(Conv,pt, numparams, dtParams)
  grad = zeros(1,numparams);
  switch numparams
      case 1
          str = 'a';
      case 2
        str = 'a b';
      case 3
        str = 'a b c';
      case 4
          str = 'a b c d';
      case 5
          str = 'a b c d ab';
      case 6
          str = 'a b c d ab ac';
      case 7
          str = 'a b c d e ab ac';
      case 8
          str = 'a b c d e ab ac ad';
      case 9
          str = 'a b c d e f ab ac ad';
      case 10
          str = 'a b c d e f ab ac ad ae';
  end
  dF = fracfact(str);
  numExp = length(dF(:,1));
  for i = 1:numparams
      sum = 0;
      for j =1:numExp
         sum = sum + Conv(end - numExp+j)*dF(j,i); 
      end
      grad(i) = sum/dtParams(i)/2/(numExp/2);
  end
end