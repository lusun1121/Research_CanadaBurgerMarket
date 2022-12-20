clear;
tic;
dat=365*439;

%%revenue Wind Wave Diesel weekend Diesel2 stock stock2 len len2 len_wspd len_wvht len_disel vess_fixed vess_fixed2
%%error har2,har,har_len,har2_len,year,opening;

final2000=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_2\com2000\final.csv');
final2001=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_2\com2001\final.csv');
final_decison1998=[final2000(:,11:end);final2001(:,11:end)];
clear  final2000  final2001;

final2002=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_2\com2002\final.csv');
final2003=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_2\com2003\final.csv');
final2004=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_2\com2004\final.csv');
final2005=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_2\com2005\final.csv');
final_decison2002=[ final2002(:,11:end); final2003(:,11:end); ...
     final2004(:,11:end); final2005(:,11:end)];
clear  final2002 final2003  final2004 final2005;


final_all=[ final_decison1998; final_decison2002];
clear final_decison1998 final_decison2002;
dlmwrite('final_all2.txt', final_all);
toc;