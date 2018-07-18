function x_new = isFeasible_simplex(x_new,x_f,Box)
   feasible = true;
   num_iterations = 10;
   j = 0;
   alpha = .6;
   n = length(x_new);
   while(true)
       for i = 1:n
         if(x_new(i) < Box(i,1) || x_new(i) > Box(i,2))
           feasible = false;
         end
       end
       for i=n/2+1:n-1
         if (x_new(i+1) >= x_new(i))
           x_new(i+1) = x_new(i) - 5;
           feasible = false;
         end
       end
       if ~feasible
           if j >= num_iterations
               x_new = x_f;
               break;
           else
               x_new = x_f + alpha*(x_new - x_f);
               j = j+1;
           end
           feasible = true;
       else
           break
       end
   end
end