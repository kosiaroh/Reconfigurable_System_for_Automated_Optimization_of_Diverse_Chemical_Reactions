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

ISTD_rt_min_254 = 5;
ISTD_rt_max_254 = 9;
sig_num = 1;
ret_time_ISTD_254 = 0;
area_ISTD_254 = 0;
if not(isempty(peak_data))
peak_range = find(peak_data(:,8) - ISTD_rt_min_254 >= 0 & ISTD_rt_max_254 - peak_data(:,8) >= 0);

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
React_rt_min_aldehyde = 1.0;
React_rt_max_aldehyde = 2;
ret_time_react_aldehyde = 0;
area_react_aldehyde = 0;
% Find all peaks between minimum and maximum residence time
if not(isempty(peak_data))
peak_range = find(peak_data(:,8) - React_rt_min_aldehyde >= 0 & React_rt_max_aldehyde - peak_data(:,8) >= 0);

for peak_num = peak_range'
    % Check for right signal
    if peak_data(peak_num,1) == sig_num
        % If greater than previous area, accept
        if peak_data(peak_num,11) > area_react_aldehyde;
            ret_time_react_aldehyde = peak_data(peak_num,8);
            area_react_aldehyde = peak_data(peak_num,11);
        end
    end
end
end

sig_num = 1;
React_rt_min_prod = 1.2;
React_rt_max_prod = 2.5;
ret_time_react_prod = 0;
area_react_prod = 0;
% Find all peaks between minimum and maximum residence time
if not(isempty(peak_data))
peak_range = find(peak_data(:,8) - React_rt_min_prod >= 0 & React_rt_max_prod - peak_data(:,8) >= 0);

for peak_num = peak_range'
    % Check for right signal
    if peak_data(peak_num,1) == sig_num
        area_react_prod = peak_data(peak_num,11)+area_react_prod;
        % If greater than previous area, accept
%         if (peak_data(peak_num,11) > area_react_prod) 
%             ret_time_react_prod = peak_data(peak_num,8);
%             area_react_prod = peak_data(peak_num,11);
%         end
    end
end
end

Conc = 0;
if area_ISTD_254 == 0
   Conc = area_react_prod; 
   %Conc = area_react_prod; 
else
   Conc = (area_react_prod)/area_ISTD_254; 
   %Conc = (area_react_prod)/area_ISTD_254; 
end

HPLC_Area(end +1,:) = [area_react_aldehyde area_react_prod area_ISTD_254];
