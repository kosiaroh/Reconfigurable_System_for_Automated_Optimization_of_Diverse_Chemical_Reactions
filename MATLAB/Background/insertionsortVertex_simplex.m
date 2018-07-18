function [x, funVal] = insertionsortVertex_simplex(x, funVal)
    % Perform one iteration of insertion sort where the last element
    % of funVal is inserted in its proper place.
    size = length(funVal);
    f_new = funVal(end);
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