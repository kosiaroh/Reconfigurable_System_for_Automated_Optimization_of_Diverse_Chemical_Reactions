xbest_arr = [];
fbest_arr = [];
ncall_arr = [];
xbest_arr2 = [];
fbest_arr2 = [];
ncall_arr2 = [];
for j = 1:50
x = [];
%Initialize variables
opt_ctr = 0;ncall = 150;ncall0 = 0;nstop0 =0;
Opt_params = [1;1;1;1;1;1;1;1;1;1;1;];
T1_span = [30, 100];C1_span = [20, 200];C2_span = [20, 200];
T2_span = [30, 100];%C3_span = [1, 15];
T3_span = [30, 100];C3_span = [0.5, 3.49];
T4_span = [30, 100];C5_span = [20, 200];
T5_span = [30, 100];C6_span = [20, 200];
fcn = 'reactor_obj_fun2';
terminate = false;
is_first_set_of_exp = true;
allData = [];
ncallbest2 = 0;
ncallbest = 0;
xbest2 = [];
fbest2 = [];
while not(terminate)
    snobfit_labview_driver_local
    %x_fun = [x(req_ctr,1) x(req_ctr,2) x(req_ctr,3)];
    x(req_ctr, 4) = round(x(req_ctr, 4));
    x_fun = [x(req_ctr,1) x(req_ctr,3) x(req_ctr,2) x(req_ctr,4)];
    %x_fun = [x(req_ctr,1) x(req_ctr,4) x(req_ctr,2) x(req_ctr,5) x(req_ctr,3) x(req_ctr,6) x(req_ctr,7)];
    %x_fun = [x(req_ctr,1) x(req_ctr,5) x(req_ctr,2) x(req_ctr,6) x(req_ctr,3) x(req_ctr,7) x(req_ctr,4) x(req_ctr,8) x(req_ctr,9)];
    %x_fun = [x(req_ctr,1) x(req_ctr,6) x(req_ctr,2) x(req_ctr,7) x(req_ctr,3) x(req_ctr,8) x(req_ctr,4) x(req_ctr,9) x(req_ctr,5) x(req_ctr,10) x(req_ctr,11)];
    f(req_ctr,:) = [feval(fcn,x_fun) sqrt(eps)];
    allData = [allData;x_fun f(req_ctr,1)];
    ncall0 = ncall0 + 1;
    opt_ctr = 1;
    [fbestn,jbest] = min(f(:,1));
    if fbestn < fbest && abs((fbestn - fbest)/fbestn) > 1e-4
        ncallbest2 = ncallbest;
        xbest2 = xbest;
        fbest2 = fbest;
        fbest = fbestn;
        xbest = x(jbest,:);
        ncallbest = ncall0;
    end
    if req_ctr == req_tot
        if is_first_set_of_exp
           % ncall0 = size(f,1);
            is_first_set_of_exp = false;
        else
           % ncall0 = ncall0 + size(f,1); % update function call counter
        end
        [fbestn,jbest] = min(f(:,1)); % best function value
        if fbestn < fbest
            
             % display current number of function values,
            % best point and function value if fbest has
            % changed
            nstop0 = 0;
        elseif ncall >= minfcall,
            nstop0 = nstop0 + 1;
        end
        
        opt_ctr = 3;
    else
        req_ctr = req_ctr +1;
    end
    % check stopping criterion
    if nstop0 >= nstop || ncall0 >= ncall, break, end
end
xbest_arr(end+1,:) = xbest;
fbest_arr(end+1,:) = fbest;
ncall_arr(end+1,:) = ncallbest;
xbest_arr2(end+1,:) = xbest2;
fbest_arr2(end+1,:) = fbest2;
ncall_arr2(end+1,:) = ncallbest2;
size(xbest_arr)
%pause(3)
end
%f(j,:) = [fval sqrt(eps)];


