fbest = inf;
matlabsavefile1 = matlabsavefile;
xl_file1 =  xl_file;
if Continue
   load(matlabsavefile);
   matlabsavefile = matlabsavefile1;
   xl_file =  xl_file1;
else
   x = [];
   IRSteadystate = 0;
   ProductVec = [];
   TwopointbaselineProductVec=[];
end

all_MSData = [];
P = [0 0 0 0 0 0];
pump_chem_name = [0 0 0 0 0 0];
if not(isempty(pump1_chem)), pump_chem_name(1) = 1; end;
if not(isempty(pump2_chem)), pump_chem_name(2) = 1; end;
if not(isempty(pump3_chem)), pump_chem_name(3) = 1; end;
if not(isempty(pump4_chem)), pump_chem_name(4) = 1; end;
if not(isempty(pump5_chem)), pump_chem_name(5) = 1; end;
if not(isempty(pump6_chem)), pump_chem_name(6) = 1; end;
% ****************************************************************************************************
if Run_at_T_SP(1)
    T1_SP = T_SP(1);
end
if Run_at_T_SP(2)
    T2_SP = T_SP(2);
end
if Run_at_T_SP(3)
    T3_SP = T_SP(3);
end
if Run_at_T_SP(4)
    T4_SP = T_SP(4);
end
if Run_at_T_SP(5)
    T5_SP = T_SP(5);
end
% ****************************************************************************************************
for i = 1:5
    Bay = BayConfig(i);
    if Bay >= 1 && Bay <=4
        if  Optimize_params(i*2-1)
            % if choice is made to optimize temperature
            if T_params(i,1) < T_params(i,2)
                T_span(i,:) = [T_params(i,1) ,  T_params(i,2)];
                Opt_params(i) = 1;
            else
                no_errors = 0;
                error_message = 'Lower bound of temperature must be lower than upper bound';
            end
        end
    else
        if Optimize_params(i*2-1)
            no_errors = 0;
            if Bay == 0
                error_message = 'Can not optimize temperature on bypass ';
            elseif Bay == 6
                error_message = 'Can not optimize temperature on a separator';
            else
                error_message = 'Can not optimize temperature on a photo reactor';
            end
        end
    end 
    
    if Optimize_params(i*2)
        % if user wants to optimize flowrates
        if i == 1
            if pump_params(i,1) < pump_params(i,2)
                p_span(i,:) = [pump_params(i,1),  pump_params(i,2)];
                p_span(i+1,:) = [pump_params(i+1,1),  pump_params(i+1,2)];
                Opt_params(5+i) = 1;Opt_params(5+i+1) = 1;
            else
                no_errors = 0;
                error_message = strcat('Lower bound of flowrate ',num2str(i),' must be lower than upper bound');
            end
            
        else
            if pump_params(i+1,1) < pump_params(i+1,2)
                p_span(i+1,:) = [pump_params(i+1,1) ,  pump_params(i+1,2)];
                Opt_params(5+i+1) = 1;
            else
                no_errors = 0;
                error_message = strcat('Lower bound of flowrate ',num2str(i+1),' must be lower than upper bound');
            end
        end
    else
        if use_flowrate(i)
            if i == 1
                if P_setrate(i) < 0
                    no_errors = 0;
                    error_message = 'Pump set flowrates can not be negative';
                else
                    P(i) = P_setrate(i);
                end
            end
            
            if P_setrate(i+1) < 0
                no_errors = 0;
                error_message = 'Pump set flowrates can not be negative';
            else
                P(i+1) = P_setrate(i+1);
            end
       elseif Bay_ratio(i)
            P(i+1) = Bay_pump_ratio(i)*(sum(P(1:i)));
        else
            if i==1
                if pump_chem_name(i)
                    no_errors = 0;
                    error_message = strcat('Bay',num2str(i),' Pump ', num2str(i),' has a chemical but no settings');
                end
            end
            if pump_chem_name(i+1)
                no_errors = 0;
                error_message = strcat('Bay',num2str(i),' Pump ', num2str(i+1),' has a chemical but no settings');
            end
        end
    end
end
% ****************************************************************************************************

