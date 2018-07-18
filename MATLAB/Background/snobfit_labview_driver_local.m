% opt_ctr keeps track of what step in the optimization the system is in
%        0 : initialization of the optimization
%        1 : returning from an experimental run
%        2 : need to collect experimental data
%        3 : need new set of experimental points from snobfit algorithm
if opt_ctr == 0
    % Initial system initialization
    tempdir = 'C:\Users\Admin\Dropbox (MIT)\multistep_optimization';
    file = 'C:\Users\Admin\Dropbox (MIT)\multistep_optimization\reductive_amination_snobfit2'; % filename for storing intermediate data;
    % after each call to SNOBFIT) are stored
    % in <file>.mat
    fcn = 'reactor_obj_fun2';         % objective function
    indices = find(Opt_params);
    %Box_pos = [T1_span; T2_span; T3_span; T4_span; T5_span;...
    %    C1_span; C2_span; C3_span; C4_span; C5_span ; C6_span];
    
    Box_pos = [T1_span; C1_span; C2_span; C3_span;];
    Box = Box_pos;
    u = Box(:,1)';  % lower box bounds
    v = Box(:,2)';    % upper box bounds
    n = length(u);        % problem dimension (number of variables)
    
    % meaningful default values for SNOBFIT parameters
    nreq = n+6;           % number of points to be generated in each call
    % to SNOBFIT
    dx = (v-u)'*1.e-5;    % resolution vector
    p = 0.5;              % probability of generating a point of class 4
    prt = 1;              % print level
    % prt = 0 prints ncall, xbest and fbest
    % if xbest has changed
    % prt = 1 in addition prints the points suggested
    % by SNOBFIT, their model function values and
    % classes after each call to SNOBFIT
    ctr = 0;
    
    % (The default for minfcall should be chosen such that the the local
    % quadratic fit around the best point uses at least once the maximal
    % number of points, which requires minfcall >= n*(n+3)+1+nreq.
    % In cases where function evaluations are not very expensive,
    % minfcall should be chosen much higher.)
    
    minfcall = n*(n+3)+1+nreq;
    nstop = 100;
    
    initialize_opt = 1;
    num_iteration = 1;
    f = [];
    % First call to snobfit
    params = struct('bounds',{u,v},'nreq',nreq,'p',p); % input structure
    [request,xbest,fbest] = snobfit(file,[],[],params,dx);
    
    req_ctr = 1;
    req_tot = size(request,1);
   
    skip_exp_ctr = 0;
    initial_set = 1;
elseif opt_ctr == 3
    % get next set of points
    if initial_set
        initial_set = 0;
        ncall0 = size(f,1) -  skip_exp_ctr; % function call counter
        [fbestn,jbest] = min(f(:,1));
        xbest = x(jbest,:);
        %fprintf(' ncall0    xbest(1:n)                 fbest')
        %disp([ncall0,xbest,fbestn]);
        ncall0;xbest;fbestn;
        nstop0 = 0;
        [request,xbest,fbest] = snobfit(file,x,f,params);
        
        req_ctr = 1;
        req_tot = size(request,1);
    else
        [request,xbest,fbest] = snobfit(file,x,f,params);
        
        req_ctr = 1;
        req_tot = size(request,1);
    end
end

if opt_ctr ~= 2
    flag = true;
    % verify that point is feasible, if not assign a value of zero
    % if so run experiment to find value
    while (opt_ctr == 3)|| flag
%         for j = req_ctr: req_tot
%             isFeasible = true;
%             for i = n/2+1:n-1
%                 if request( req_ctr ,i) < request( req_ctr ,i+1)
%                     temp = request( req_ctr ,i) - 30;
%                     if temp >=  u(i)
%                         request( req_ctr ,i+1) = temp;
%                     else
%                         isFeasible = false;
%                         break;
%                     end
%                 end
%             end
%             
%             x(req_ctr,:) = request(req_ctr,1:n);
%             
%             if isFeasible
%                 opt_ctr = 2; % need to perform experiment
%                 break;
%             else
%                 f(req_ctr,:) = [0 sqrt(eps)];
%                 skip_exp_ctr = skip_exp_ctr+ 1;
%                 req_ctr = req_ctr + 1;
%                 opt_ctr =  3; % if it doesn't find an experiment to perform then it will exit looking for a new set
%             end
%             
%         end
        x(req_ctr,:) = request(req_ctr,1:n);
        opt_ctr = 2; % need to perform experiment
        if req_ctr > req_tot
            opt_ctr = 1;
            req_ctr = req_ctr - 1;
        end
        
        flag = false;
    end
end