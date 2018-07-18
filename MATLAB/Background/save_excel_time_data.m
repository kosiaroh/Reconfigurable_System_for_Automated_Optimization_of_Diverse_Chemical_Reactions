   temp_time_cond= [CurrExp, 1, P1, P2, P3, P4, P5,P6, T1, T2, T3, T4, T5, T1_SP, T2_SP, T3_SP, T4_SP, T5_SP, F1, F2, Pres1, Pres2];
   lastrowStr1 = num2str(lastrow1);
   xlSheets.Item(idx1).Range(strcat('A',lastrowStr1,':',lastcol1,lastrowStr1)).Value2 =   temp_time_cond;
   lastrow1 = lastrow1+1;
  xlWorkbook.SaveAs(xl_file);