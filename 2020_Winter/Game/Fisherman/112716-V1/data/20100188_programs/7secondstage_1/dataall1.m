clear;
tic;
dat=365*439;

final2000=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2000\final.csv');
deci2000=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2000\indiv_daily.txt');
%year, y365, vess, fishing
final2001=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2001\final.csv');
deci2001=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2001\indiv_daily.txt');
final_decison1998=[deci2000 final2000(:,1:10);deci2001 final2001(:,1:10)];
clear deci2000 final2000 deci2001 final2001;

final2002=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2002\final.csv');
deci2002=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2002\indiv_daily.txt');
final2003=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2003\final.csv');
deci2003=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2003\indiv_daily.txt');
final2004=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2004\final.csv');
deci2004=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2004\indiv_daily.txt');
final2005=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2005\final.csv');
deci2005=load('C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2005\indiv_daily.txt');
final_decison2002=[deci2002 final2002(:,1:10);deci2003 final2003(:,1:10); ...
    deci2004 final2004(:,1:10);deci2005 final2005(:,1:10)];
clear deci2002 final2002 deci2003 final2003 deci2004 final2004 deci2005 final2005;

year1998=[2000*ones(dat,1);2001*ones(dat,1)];
year2002=[2002*ones(dat,1);2003*ones(dat,1);2004*ones(dat,1);2005*ones(dat,1)];

final_all=[year1998 final_decison1998;year2002 final_decison2002];
clear final_decison1998 final_decison2002;
dlmwrite('final_all1.txt', final_all);
toc;
