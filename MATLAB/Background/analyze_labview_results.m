if not(system_error) && not(redo)
  if opt_ctr == 1
      if throughput_yield
         f(req_ctr,:) = [-(C_measured)/total_wait_time*60 sqrt(eps)];        
      else
         f(req_ctr,:) = [-(C_measured) sqrt(eps)];   
      end
  end
  temp_write = [CurrExp, f(req_ctr,1), abs(C_measured)];
  lastrowStr4 = num2str(lastrow4);
  xlSheets.Item(idx4).Range(strcat('A',lastrowStr4,':',lastcol4,lastrowStr4)).Value2 =  temp_write;
  lastrow4 = lastrow4 +1;
  xlWorkbook.SaveAs(xl_file);
  allHPLCData(end+1,:) =  [CurrExp new_x(:)' f(req_ctr,1) abs(C_measured)];
  experiment_history(CurrExp,:) = [P1 P2 P3 P4 P5 P6 T1_SP T2_SP T3_SP T4_SP T5_SP f(req_ctr,1) abs(C_measured)];
  %allHPLCFile(end+1,:) = [num2str(CurrExp) hplc_file];

  opt_ctr = 1;
  CurrExp = CurrExp + 1;
  terminate = 0;

    if is_first_set_of_exp
       ncall0 = size(f,1);
       is_first_set_of_exp = false;
    else
       ncall0 = ncall0 +1;% size(f,1); % update function call counter
    end
    [fbestn,jbest] = min(f(:,1)); % best function value
    
    if is_redo_exp
        if f_redo == fbest
            xbest =  xbest_prev;
            fbest = fbest_prev;
        end
        is_redo_exp = 0;
    end
    
    if f(req_ctr,1) <= fbest
        fbest_prev = fbest;
        fbest = fbestn;
        xbest_prev = xbest;
        xbest = x(jbest,:);
        Tbest = [T1_SP T2_SP T3_SP T4_SP T5_SP];
        Pbest = [P1 P2 P3 P4 P5 P6];
        best_exp_num = CurrExp - 1;
        ncall0,xbest;fbest % display current number of function values,
        % best point and function value if fbest has changed
        nstop0 = 0;
    end
  
      % check stopping criterion
    if nstop0 >= nstop  || ncall0 >= ncall
       terminate = 1;
       lastrowStr_temp = num2str(10);
       lastcol_temp =  xlcolumnletter(length(xbest)+1);
       xlSheets.Item(idx3).Range(strcat('B',lastrowStr_temp,':',lastcol_temp,lastrowStr_temp)).Value2 =  xbest;
       xlWorkbook.SaveAs(xl_file);
    end

    if req_ctr == req_tot
        if ncall >= minfcall && fbestn > fbest
           nstop0 = nstop0 + 1;
        end
       opt_ctr = 3;
    else
       req_ctr = req_ctr +1;
    end

   save(matlabsavefile);
end