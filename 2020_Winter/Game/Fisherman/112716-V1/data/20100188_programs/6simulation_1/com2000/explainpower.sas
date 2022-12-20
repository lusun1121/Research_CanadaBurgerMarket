libname pros 'C:\Documents and Settings\LHuang\Desktop\code\pros';
proc import datafile='C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\6simulation_1\com2000\reg.txt'
out=regg  dbms=dlm
replace;
delimiter=",";getnames=no;
run;
quit;


data regg1;
set regg;
mape1=abs(var2-var1)/var2;
mape2=abs(var4-var3)/var4;
run;


proc sql;
create table regg2 as
select distinct *, sum(var4) as has
from regg1;
quit;

proc sql;
create table regg3 as
select distinct sum(mape1*var4/has) as mape_d1, sum(mape2*var4/has) as mape_d2
from regg2;
quit;



proc reg data = regg1;
  model  var2=var1;
run;
quit;

proc reg data = regg1;
  model  var4=var3;
run;
quit;

/*proc corr data=regg1;
var var1 var2;
run;
