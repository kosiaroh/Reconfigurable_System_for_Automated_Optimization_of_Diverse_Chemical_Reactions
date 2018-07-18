if not(system_error)

    if not(redo)
        % snobfit
        P1 = P(1);P2 = P(2);P3 = P(3);P4 = P(4);P5 = P(5);P6 = P(6);
        snobfit_labview_driver;
        
        if opt_ctr == 2
            new_x = x(req_ctr,:);
            x_ind = 1;
            
            if Opt_params(1)
                T1_SP = new_x(x_ind ); %'C
                params_1(1) = new_x(x_ind);
                x_ind  = x_ind  + 1;
            end
            if Opt_params(2)
                T2_SP = new_x(x_ind ); %'C
                params_1(3) = new_x(x_ind);
                x_ind  = x_ind  + 1;
            end
            if Opt_params(3)
                T3_SP = new_x(x_ind ); %'C
                x_ind  = x_ind  + 1;
            end
            if Opt_params(4)
                T4_SP = new_x(x_ind ); %'C
                x_ind  = x_ind  + 1;
            end
            if Opt_params(5)
                T5_SP = new_x(x_ind ); %'C
                x_ind  = x_ind  + 1;
            end
            if Opt_params(6)
                %tau1 = Vr1/(new_x( x_ind) +new_x( x_ind+1) );
                %new_x( x_ind) = tau1*60;
                params_1(2) = new_x(x_ind);
                P1=new_x( x_ind); %ul/min
                P2 =new_x( x_ind+1); %ul/min
                x_ind  = x_ind  + 2;
            end
            if Opt_params(8)
                tau2= Vr2/(new_x( x_ind)+P1+P2);
                %new_x( x_ind) = tau2*60;
                params_1(4) = new_x(x_ind);
                %P3=new_x( x_ind); %ul/minP2 = Vr4/(tau4)*60/2;
                P3 = new_x( x_ind);
                x_ind  = x_ind  + 1;
                % P2 = Vr4/(tau4)*60/2;
            elseif Bay2_ratio && not(use_flowrate2)
                P3= Bay2_pump_ratio*(P1+P2); %chnaged from  P3= Bay2_pump_ratio*(P1 + P2);
            end
            if Opt_params(9)
                tau3 =Vr3/(new_x( x_ind)+P1+P2+P3);
                P4=new_x( x_ind); %ul/min
                %P3= Bay2_pump_ratio*(P4);
                x_ind  = x_ind  + 1;
            elseif Bay3_ratio && not(use_flowrate3)
                P4= Bay3_pump_ratio*(P1+P2 +P3);
            end
            if Opt_params(10)
                tau4 = Vr4/(new_x( x_ind)+P1+P2+P3+P4);
                %P5=new_x( x_ind) ; %ul/min
                P5=new_x( x_ind);
                x_ind  = x_ind  + 1;
            elseif Bay4_ratio && not(use_flowrate4)
                P5= Bay4_pump_ratio*(P1+P2 +P3 +P4);
            end
            if Opt_params(11)
                tau5 = Vr5/(new_x( x_ind)+P1+P2+P3+P4+P5);
                %new_x( x_ind) = tau5*60;
                %P6=new_x( x_ind) ; %ul/min
                P6=new_x( x_ind);
                x_ind  = x_ind  + 1;
            elseif Bay5_ratio && not(use_flowrate5)
                P6= Bay5_pump_ratio*(P1+P2 +P3 +P4 +P4 +P5);
            end
            %---------
        end
       if syringe_pump
        Conditions_Array = [P1*2 P2 P3*2 P4*2 P5*2 P6*2 T1_SP T2_SP T3_SP T4_SP T5_SP];
       else
        Conditions_Array = [P1*2 P2*2 P3*2 P4*2 P5*2 P6*2 T1_SP T2_SP T3_SP T4_SP T5_SP];   
       end
        if opt_ctr == 2
            opt_ctr = 1;
        else
            opt_ctr = 2;
        end
    else
        opt_ctr =1;
        redo = false;
        is_redo_exp = 1;
        f_redo = f(end,1);
        if mod(CurrExp,req_tot) ~= 1
            req_ctr = req_ctr - 1;
        end
        CurrExp = CurrExp - 1;
        ncall0 = ncall0 - 1;
    end

end
P_tot = P1 + P2 + P3 + P4 + P5 + P6 ;

tau1 = 0; tau2 = 0; tau3 = 0; tau4 = 0; tau5 = 0;
P_temp = P1+P2;
if P_temp >0
    tau1 = Vr1/P_temp;
end
P_temp = P1+P2+P3;
if P_temp >0
    tau2 = Vr2/P_temp;
end
P_temp = P1+P2+P3+P4;
if P_temp >0
    tau3 = Vr3/P_temp;
end
P_temp = P1+P2+P3+P4+P5;
if P_temp >0
    tau4 = Vr4/P_temp;
end
P_temp = P1+P2+P3+P4+P5+P6;
if P_temp >0
    tau5 = Vr5/P_temp;
end
experiment_history(CurrExp,:) = [P1 P2 P3 P4 P5 P6 T1_SP T2_SP T3_SP T4_SP T5_SP 0 0];
tau_array = [tau1, tau2, tau3, tau4, tau5];