function [all_x, new_x, curr_step,f_r,x_r,funValArr] = nedler_mead_simplex(old_x,curr_step, funValArr,newFunVal,new_x, x_r,f_r,Box)

%% Function performs a step of the simplex method to determine the next operating point
%  Parameters:
%  old_x - array of n+1 vertices of current simplex
%  curr_step - integer tells where in the algorithm we are
%              0 => start with reflection
%              1 => get refection data
%              11 => got reflection data
%              2 => get expansion data
%              22 => got expansion data
%              3 => get ouside contraction data
%              33 => got outside contraction data
%              4 => get inside contraction data
%              44 => got inside contraction data
%              6 => get data for shrink points
%  funValArr - array of function value corresponding to verticies
%  newFunVal - 
%  new_x - 
%  Box - bounds on parameters
     n = length(old_x(1,:)); % number of dimensions
     alpha = 1; beta = 1 + 2/n; gamma = .75 - 1/(2*n); delta = 1 - 1/n;
     x_bar = mean(old_x(1:end-1,:));
    x_n = old_x(n,:);
    f_n = funValArr(n);
    x_n_plus_one = old_x(n+1,:);
    f_n_plus_one = funValArr(n+1);
    while(true)% FIX
        if curr_step ==0
            %*** Reflection
            x_r = x_bar + alpha*(x_bar - x_n_plus_one);
            x_r = isFeasible_simplex(x_r,x_n_plus_one,Box);
            curr_step = 1;
        elseif curr_step == 11 
            %2---------- get experimental data here
            x_r = new_x;
            f_r = newFunVal;
            %----------
            curr_step = 0;
        end
        x_1 = squeeze(old_x(1,:));
        f_1 = funValArr(1);
        all_x = old_x;
        if (curr_step ==0 && (f_1 <= f_r && f_r < f_n))
            all_x(n+1,:) = x_r;
            funValArr(n+1) = f_r;
            [all_x,funValArr] = insertionsortVertex_simplex(all_x, funValArr);
    %         fprintf('Reflection - \n');
    %         disp(x);
    %         fprintf('Reflection FVal - %d %d %d %d\n',funValArr');
        elseif (curr_step ==22 || (curr_step == 0 && f_r < f_1))
            if curr_step ==22  
                %3 ----- get experimental data here
                f_e = newFunVal;
                %------
                if f_e < f_r
                   all_x(n+1,:) = new_x; 
                   funValArr(n+1) = f_e;
                else
                   all_x(n+1,:)= x_r; 
                   funValArr(n+1) = f_r;
                end
                %fprintf('Expansion FVal - %d %d %d %d %d\n',funValArr');disp(x);
                [all_x,funValArr] = insertionsortVertex_simplex(all_x, funValArr);
                curr_step = 0;
            else
                  %*** Expansion
                x_e = x_bar + beta*(x_r - x_bar);
                new_x = isFeasible_simplex(x_e,x_n_plus_one,Box);  
                curr_step = 2;
            end
    %         fprintf('Expansion -\n');
    %         disp(x);
    %         fprintf('Expansion FVal - %d %d %d %d %d\n',funValArr');
        elseif (curr_step ==33 || (curr_step == 0 && f_n <= f_r && f_r <= f_n_plus_one))   
            if curr_step ==33
                %-------
                %4 ------ get experimental data here
                f_oc = newFunVal;
                if(f_oc <= f_r)
                    all_x(n+1,:) = new_x; 
                    funValArr(n+1) = f_oc;
                    [all_x,funValArr] = insertionsortVertex_simplex(all_x, funValArr);
                    curr_step =0;
        %             fprintf('Contraction1 -\n');
        %             disp(x);
        %             fprintf('Contraction1 FVal - %d %d %d %d %d\n',funValArr');
                else
                    %**** Shrink
                    for i = 2:n
                       all_x(i,:) =  x_1 + delta*(all_x(i,:) - x_1);
                       
                       %funValArr(i) = functionValue(x(ind,i,:));
                    end
                    curr_step = 6;
                    %*** Sort
                    %[x(ind+1,:,:), funValArr] = sortVertex(x(ind+1,:,:),funValArr);
        %             fprintf('Shrink1 -\n');
        %             disp(x);
        %             fprintf('Shrink2 FVal - %d %d %d %d %d\n',funValArr');
                end
            else
                %*** Outside Contraction
                x_oc = x_bar - gamma*(x_r - x_bar);
                new_x = isFeasible_simplex(x_oc,x_n_plus_one,Box);
                curr_step =3;
            end
            
        elseif (curr_step ==44 || (curr_step == 0 && f_r >= f_n_plus_one))
            if curr_step ==44
                %5 ------ get experimental data here
                f_ic = newFunVal;
                %-------
                if f_ic < f_n_plus_one
                    all_x(n+1,:) = new_x;
                    funValArr(n+1) = f_ic;
                    [all_x,funValArr] = insertionsortVertex_simplex(all_x, funValArr);
                    curr_step =0;
        %          fprintf('Contraction2 -\n');
        %          disp(x);
        %          fprintf('Contraction2 FVal - %d %d %d %d %d\n',funValArr');
                else
               %*** Shrink
                    for i = 2:n
                       all_x(i,:) =  x_1 + delta*(all_x - x_1);
                       %funValArr(i) = functionValue(x(ind,i,:));
                    end
                    curr_step = 6;
                    %*** Sort
                    %[x(ind+1,:,:), funValArr] = sortVertex(x(ind+1,:,:),funValArr);
        %             fprintf('Shrink2 -\n');
        %             disp(x);
        %             fprintf('Shrink2 FVal - %d %d %d %d %d\n',funValArr');
                end
            else
                %*** Inside Contraction
                x_ic = x_bar - gamma*(x_r - x_bar);
                new_x = isFeasible_simplex(x_ic,x_n_plus_one,Box);
                curr_step =3;
            end
        end
        
        if curr_step ~=0
            break;
        end
    end
end