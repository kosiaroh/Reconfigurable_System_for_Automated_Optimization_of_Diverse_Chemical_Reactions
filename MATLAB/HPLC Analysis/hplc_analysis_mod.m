% % Build file name
% file_name = [data_path '\REPORT01.xls'];

% Read signals from file_name
[values,signals] = xlsread(hplc_file,'Signal','E2:E100');

% Read signal numbers fron file_name
% sig_nums = xlsread(hplc_file,'Signal','N2:N100');

% Read peak data from file_name
peak_data = xlsread(hplc_file,'Peak','D2:P1000');

%Reductive Amination
% Conc ~ areaA(254)/areaISTD(254) areaA(210)/areaISTD(210)
% ISTD_rt_min_210 = 4.1;
% ISTD_rt_max_210 = 5;

sig_num = 1;
ret_time_ISTD_254 = 0;
area_ISTD_254 = 0;
if not(isempty(peak_data)) && use_internal_std
    peak_range = find(peak_data(:,8) - ISTD_start >= 0 & ISTD_end - peak_data(:,8) >= 0);
    
    for peak_num = peak_range'
        % Check for right signal
        if peak_data(peak_num,1) == sig_num
            area_ISTD_254 = peak_data(peak_num,11)+area_ISTD_254;
            % If greater than previous area, accept
            %         if peak_data(peak_num,11) > area_ISTD_254;
            %             ret_time_ISTD_254 = peak_data(peak_num,8);
            %             area_ISTD_254 = peak_data(peak_num,11);
            %         end
        end
        
    end
end

sig_num = 1;
ret_time_react_prod = 0;
area_react_prod = 0;
% Find all peaks between minimum and maximum residence time
if not(isempty(peak_data))
peak_range = find(peak_data(:,8) - hplc_peak_start >= 0 & hplc_peak_end - peak_data(:,8) >= 0);

for peak_num = peak_range'
    % Check for right signal
    if peak_data(peak_num,1) == sig_num
       % area_react_prod = peak_data(peak_num,11)+area_react_prod;
        % If greater than previous area, accept
         if (peak_data(peak_num,11) > area_react_prod) 
             ret_time_react_prod = peak_data(peak_num,8);
             area_react_prod = peak_data(peak_num,11);
         end
    end
end
end

Conc = 0;
if area_ISTD_254 == 0
   Conc = area_react_prod; 
   %%Remove Next Line if not a simulation
   %Conc = P1*P2/T4/P3*P4*T2;
else
   Conc = (area_react_prod)/area_ISTD_254; 
    %%Remove Next Line if not a simulation
   %Conc = P1*P2/P5;
end


if Conc == 0
    if not(dcr_hplc_peak)
        Conc = -1;
    else
        Conc = 1;
    end
else
    if dcr_hplc_peak
        Conc = -Conc;
    end
end

HPLC_Area(end +1,:) = [area_react_prod area_ISTD_254];
