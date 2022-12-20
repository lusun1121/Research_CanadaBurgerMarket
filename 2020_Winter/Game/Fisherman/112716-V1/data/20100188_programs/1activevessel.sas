libname shrimp 'C:\Documents and Settings\LHuang\Desktop\code\shrimp';
libname pros 'C:\Documents and Settings\LHuang\Desktop\code\pros';

/***************************************************************************************/
/***************active vess only********************************************************/
/***************************************************************************************/
Proc sql;
create table active0 as
select distinct *,mean(length) as len
from shrimp.shrimp_new
where vesselid~=. and started~=. and landed~=. and pounds~=. and year>=2000
group by vesselid;
quit;
data active1;
set active0;
if started=. then started=landed;
active=landed-started+1;
if active>10 then delete;
run;

Proc sql;
create table active2 as
select distinct vesselid,started,len,round(mean(active)) as active1,sum(pounds) as poundss
from active1
group by vesselid, started;
quit;


Proc sql;
create table active3 as
select distinct *,sum(active1) as valid
from active2
group by vesselid;
quit;

data pros.active;
set active3;
where valid>100 & vesselid~=.;
run;
