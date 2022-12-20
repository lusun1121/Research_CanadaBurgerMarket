libname pros 'C:\Documents and Settings\LHuang\Desktop\code\pros';



proc sql;
create table simu1 as
select distinct year,y365,stock, prize_a, wspd,
wvht, holiday,week,opening1,sum(pounds_a) as har,sum(fishing) as Tvessel
from pros.firststage
group by year,y365;
quit;




proc sql;
create table simu2 as
select distinct y365, mean(prize_a) as prize_a1,
mean(wspd) as wspd1, mean(wvht) as wvht1
from simu1
group by y365;
quit;
proc sort data=simu1;
by  y365;
run;
proc sort data=simu2;
by  y365;
run;
data simu3;
merge simu1 simu2;
by  y365;
run;
data daily;
set simu3;
if y365=. then delete;
if prize_a=. then prize_a=prize_a1;
if prize_a=. then prize_a=2;
if prize_a>5 then prize_a=3;
if wspd=. then wspd=wspd1;
if wvht=. then wvht=wvht1; 
if wspd=0 then wspd=7;
if wvht =0 then wvht=1;
if har=. then har=0;
drop prize_a1  wspd1 wvht1;
run;

proc sort data=daily;
by year y365;
run;

PROC EXPORT DATA=daily
OUTFILE='C:\Documents and Settings\LHuang\Desktop\code\SNPL2\6simulation_1\com2000\daily0.txt'
replace;
 putnames=no;
quit;

PROC EXPORT DATA=pros.vess_fixed 
OUTFILE='C:\Documents and Settings\LHuang\Desktop\code\SNPL2\6simulation_1\com2000\indiv.txt'
replace;
 putnames=no;
quit;

proc sql;
create table ze as
select distinct year,y365,zeta
from pros.firststage;
quit;

data ze1;
set ze;
if y365=. then delete;
run;
proc sort data=ze1;
by year y365;
run;


PROC EXPORT DATA=ze1
OUTFILE='C:\Documents and Settings\LHuang\Desktop\code\SNPL2\6simulation_1\com2000\zeta.txt'
replace;
 putnames=no;
quit;

/*Do the simnulation in folder from com2004 to com2004 seperately*/
/*merge all the data in folder secondstage and do the MLE in 6MLE.sas*/



data firststage;
set pros.firststage;
if pounds_a=. then pounds_a=0;
run;

proc sql;
create table indiv_daily as
select distinct year, y365, vess, fishing
from pros.firststage
where y365~=. & year=2005; 
quit;
proc sort data=indiv_daily;
by year  y365 vess;
run;

PROC EXPORT DATA=indiv_daily 
OUTFILE='C:\Documents and Settings\LHuang\Desktop\code\SNPL2\6simulation_1\com2005\indiv_daily.txt'
replace;
 putnames=no;
quit;




