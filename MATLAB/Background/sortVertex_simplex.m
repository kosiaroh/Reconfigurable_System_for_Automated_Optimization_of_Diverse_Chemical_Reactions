function [x, fun] = sortVertex_simplex(x, fun)
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