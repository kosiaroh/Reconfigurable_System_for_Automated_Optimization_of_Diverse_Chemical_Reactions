function x = returnDOE(numparams,paramDelta, centerPt)
  % assumes parameter comes in pairs
  str = '';
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
  %disp(dF)
  temp = repmat(paramDelta,[length(dF(:,1)) 1]);
  %disp(temp)
  x = repmat(centerPt, [length(dF(:,1)) 1]) + temp.* dF;
end