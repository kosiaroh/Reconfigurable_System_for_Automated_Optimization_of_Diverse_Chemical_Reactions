%Setup output excel file
% Opening Excel file and inputing first peaks headers if needed
try
  xlApp = actxGetRunningServer('Excel.Application');
catch ME
  xlApp = actxserver('Excel.Application');
end

xlApp.DisplayAlerts = false; 
xlApp.Visible = false;

try
   xlWorkbook = xlApp.workbooks.Open(xl_file, 0,false);
   xlSheets = xlWorkbook.Sheets;
   xlSheetNamesArray = cell(xlSheets.Count, 1);
catch
    copyfile('C:\Users\HPLC\Dropbox (MIT)\MoD Files\Data\Optimization_template_file.xlsx', xl_file);
    %copyfile('C:\Users\Admin\Dropbox (MIT)\MoD Files\Data\Optimization_template_file.xlsx', xl_file);
    xlWorkbook = xlApp.workbooks.Open(xl_file, 0,false);
   xlSheets = xlWorkbook.Sheets;
   xlSheetNamesArray = cell(xlSheets.Count, 1);
end

  for i = 1:xlSheets.Count
     xlSheetNamesArray{i} = xlSheets.Item(i).Name;
  end

  [~, idx1] = ismember('Time_Series_Conditions', xlSheetNamesArray);
  lastcol1 = xlcolumnletter(22);

  [~, idx2] = ismember('MS Results', xlSheetNamesArray);
  lastcol2 = xlcolumnletter(2);

  [~, idx3] = ismember('System Setup', xlSheetNamesArray);
  lastcol3 = xlcolumnletter(11);

  [~, idx4] = ismember('HPLC Results', xlSheetNamesArray);
  lastcol4 = xlcolumnletter(3);

  [~, idx_ir_exp] = ismember('IR Results', xlSheetNamesArray);
  IR_i= xlSheets.Item(idx_ir_exp).Range('A1').End('xlDown').Row ; 
  if IR_i == 1048576
     IR_i = 2;
  end
 
  lastrow1 =  xlSheets.Item(idx1).Range('A1').End('xlDown').Row;
  if lastrow1 == 1048576
      lastrow1 = 2;
  end

  lastrow2 =  xlSheets.Item(idx2).Range('A1').End('xlDown').Row;
  if lastrow2 == 1048576
      lastrow2 = 2;
  end

   lastrow3= 2;

    lastrow4=  xlSheets.Item(idx4).Range('A1').End('xlDown').Row;
    if lastrow4 == 1048576
       lastrow4 =2;
    end

%catch ME
 %   xlWorkbook = xlApp.Workbooks;
 %  xlWorkbook = invoke(xlWorkbook, 'Add');
 %   WS = xlWorkbook.Worksheets;
 %   row = num2str(1);
 %   temp_write = [Bay1, Bay2, Bay3, Bay4, Bay5, 1, Vr1, Vr2, Vr3, Vr4, Vr5];
 %   xlSheets.Item(idx3).Range(strcat('A',row,':',lastcol3,row)).Value2 =  temp_write;
 %   WS.Add([], WS.Item(WS.Count));
 %   WS.Item(WS.Count).Name = 'Time_Series_Conditions';
 %    WS.Add([], WS.Item(WS.Count));
 %    WS.Item(WS.Count).Name = 'MS Results';
 %   WS.Add([], WS.Item(WS.Count));
 %    WS.Item(WS.Count).Name = 'System Setup';
 %   WS.Add([], WS.Item(WS.Count));
 %    WS.Item(WS.Count).Name = 'HPLC Results';
%end

load('C:\Users\HPLC\Dropbox (MIT)\MoD Files\MATLAB\IR_Pks.mat');
%load('C:\Users\Admin\Dropbox (MIT)\MoD Files\MATLAB\IR_Pks.mat');
temp_write = [Bay1, Bay2, Bay3, Bay4, Bay5, 1, Vr1, Vr2, Vr3, Vr4, Vr5];
lastrowStr3 = num2str(lastrow3);
xlSheets.Item(idx3).Range(strcat('A',lastrowStr3,':',lastcol3,lastrowStr3)).Value2 =  temp_write;
lastrow3 = lastrow3 +1;

temp_write = {pump1_chem, pump2_chem, pump3_chem, pump4_chem, pump5_chem, pump6_chem};
lastrowStr3 = num2str(6);
xlSheets.Item(idx3).Range(strcat('B',lastrowStr3,':','G',lastrowStr3)).Value2 =  temp_write;

lastrowStr3 = num2str(8);
xlSheets.Item(idx3).Range(strcat('B',lastrowStr3,':','G',lastrowStr3)).Value2 = Pump_concen;

experiment_history = [];