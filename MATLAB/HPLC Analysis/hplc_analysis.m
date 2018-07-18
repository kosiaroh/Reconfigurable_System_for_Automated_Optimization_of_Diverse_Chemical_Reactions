% % Build file name
% file_name = [data_path '\REPORT01.xls'];

% Read signals from file_name
[values,signals] = xlsread(hplc_file,'Signal','E2:E100');

% Read signal numbers fron file_name
% sig_nums = xlsread(hplc_file,'Signal','N2:N100');

% Read peak data from file_name
peak_data = xlsread(hplc_file,'Peak','D2:P1000');

%Reductive Amination
% Conc ~ areaA/areaISTD
ISTD_rt_min = 3.3;
ISTD_rt_max = 3.8;

% Find all peaks between minimum and maximum residence time
% peak_range = find(peak_data(:,8) - ISTD_rt_min >= 0 & ISTD_rt_max - peak_data(:,8) >= 0);
sig_num = 1;
ret_time_ISTD = 0;
area_ISTD = 0;
% for peak_num = peak_range'
%     % Check for right signal
%     if peak_data(peak_num,1) == sig_num
%         % If greater than previous area, accept
%         if peak_data(peak_num,11) > area_ISTD;
%             ret_time_ISTD = peak_data(peak_num,8);
%             area_ISTD = peak_data(peak_num,11);
%         end
%     end
%     
% end

React_rt_min = 2;
React_rt_max = 3;
% Find all peaks between minimum and maximum residence time
peak_range = find(peak_data(:,8) - React_rt_min >= 0 & React_rt_max - peak_data(:,8) >= 0);
ret_time_react = 0;
area_react = 0;
for peak_num = peak_range'
    % Check for right signal
    if peak_data(peak_num,1) == sig_num
        % If greater than previous area, accept
        if peak_data(peak_num,11) > area_react;
            ret_time_react = peak_data(peak_num,8);
            area_react = peak_data(peak_num,11);
        end
    end
    
end

if area_ISTD == 0
   Conc = area_react; 
else
   Conc = area_react/area_ISTD; 
end
