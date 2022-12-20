libname shrimp 'C:\Documents and Settings\LHuang\Desktop\code\shrimp';

libname pros 'C:\Documents and Settings\LHuang\Desktop\code\pros';
/***************************************************************************************/
/***************frame the time and  vessel**********************************************/
/***************************************************************************************/
data frame0;
do vess=1 to 439;
    do year=2000 to 2005;
	   do month=1 to 12;
	       do day=1 to 31;
		   output;
		   end;
		end;
	end;
end;
run;
data frame;
set frame0;
   	time=mdy(month,day,year);
	format time mmddyy10.;
	if time=. then delete;
	week=week(time);
	if weekday(time)=6 or weekday(time)=7 then holiday=1;
    else holiday=0;
if year=2000 then year2000=1;else year2000=0;
if year=2001 then year2001=1;else year2001=0;
if year=2002 then year2002=1;else year2002=0;
if year=2003 then year2003=1;else year2003=0;
if year=2004 then year2004=1;else year2004=0;
if year=2005 then year2005=1;else year2005=0;
run;
proc sort data=frame;
by month day;
run;
proc sort data=pros.y365;
by month day;
run;
data frame1;
merge pros.y365 frame;
by month day;
run;
/***************************************************************************************/
/***************merge with decisions****************************************************/
/***************************************************************************************/
proc sort data=frame1;
by vess year month day;
run;
proc sql;
create table deci1 as
select distinct vess,year,month,day,fishing,pounds_a
from pros.stock;
quit;

proc sort data=deci1;
by vess year month day;
run;

data deci2;
merge deci1 frame1;
by vess year month day;
run;

data frame_vess;
set deci2;
if year2000=. then delete;
if fishing=. then fishing=0;
fishing0=1-fishing;
if pounds_a=. then pounds_a=0;
yt=log(y365);
yt2=yt*yt;
yt3=yt2*yt;
yt4=yt3*yt;
run;

/***************************************************************************************/
/***************merge with len and catchability*****************************************/
/***************************************************************************************/

proc sort data=frame_vess;
by vess;
run;
proc sort data=pros.vess_fixed;
by vess;
run;
data frame1;
merge frame_vess pros.vess_fixed;
by vess;
run;
data frame1;
set frame1;
vess_fixede=exp(vess_fixed)*fishing;
run;

/***************************************************************************************/
/***************accumulated harvest*****************************************************/
/***************************************************************************************/
proc sql;
create table accun as
select distinct year, month, day, sum(pounds_a) as poun
from frame1
group by year, month, day;
quit;

proc sort data=accun ;
by year month day;
run;
DATA accun;
set accun;
if poun =. then poun=0;
run;


DATA accumulated_harn;
    set accun;
	by year;
 if first.year then poundss_sum=0;
    poundss_sum+poun;
run;


proc sort data=accumulated_harn;
by year month day;
run;

proc sort data=frame1;
by year month day;
run;

data frame_vess1;
merge frame1 accumulated_harn;
by year month day;
run;

/***************************************************************************************/
/*************************merge with zeta***********************************************/
/***************************************************************************************/

proc sort data=pros.zeta;
by year y365;
run;
proc sort data=frame_vess1;
by year y365;
run;
data frame_vess2;
merge frame_vess1 pros.zeta;
by year y365;
run;

/***************************************************************************************/
/******************predict stock Index**************************************************/
/***************************************************************************************/
data frame_stock;
set frame_vess2;
stock=exp(0.577319617*year2000+ 0.066830682*year2001+0.689766442*year2002+0.381456412*year2003
+ 0.080095902*year2004+9.854830539*yt -7.394585900*yt2+ 1.796771225*yt3
   -0.137313118*yt4 -0.000000127*poundss_sum);
if y365<11 then stock=10;
zeta1=exp(zeta);
run;

 /*data frame_stock;
set frame_vess1;
stock=exp(1176.40083*year2000 -7155.61407*year2001+743.65283*year2002+ -8868.67070*year2003
 -11167.85934*year2004 -18793.66365*year2005+ 5912.09123*yt2 -3138.06848*yt3
    +409.15885*yt4 -1*effort_suml);
run;*/


/***************************************************************************************/
/***************merge with fuel price***************************************************/
/***************************************************************************************/
PROC EXPORT 
data=pros.diesel
OUTFILE='C:\Documents and Settings\LHuang\Desktop\diesel.csv'
replace;
quit;

data diesel;
set pros.diesel;
dieselp=diesel_price;
drop diesel_price;
run;
data pros.diesel1;
set diesel;
where year>=2000;
if year=2000 then diesel_price=dieselp*0.86/0.86;
if year=2001 then diesel_price=dieselp*0.84/0.86;
if year=2002 then diesel_price=dieselp*0.82/0.86;
if year=2003 then diesel_price=dieselp*0.81/0.86;
if year=2004 then diesel_price=dieselp*0.78/0.86;
if year=2005 then diesel_price=dieselp*0.76/0.86;
run;
proc sort data=pros.diesel1;
by year week;
run;
proc sort data=frame_stock;
by year week;
run;
data frame_stock1;
merge pros.diesel1 frame_stock;
by year week;
run;
/***************************************************************************************/
/***************merge with shrimp price*************************************************/
/***************************************************************************************/

proc sql;
create table stock00 as
select distinct started,sum(value) as values,sum(pounds) as pounds_su
from shrimp.shrimp_new
group by started;
quit;



proc sql;
create table stock01 as
select distinct started,values/pounds_su as prize_a_non
from stock00;
quit;

data pros.price;
set stock01;
year=year(started);
month=month(started);
day=day(started);
if year<2000 then delete;
if year=2000 then prize_a=prize_a_non*0.86/0.86;
if year=2001 then prize_a=prize_a_non*0.84/0.86;
if year=2002 then prize_a=prize_a_non*0.82/0.86;
if year=2003 then prize_a=prize_a_non*0.81/0.86;
if year=2004 then prize_a=prize_a_non*0.78/0.86;
if year=2005 then prize_a=prize_a_non*0.76/0.86;
where started~=.;
run;




proc sort data=pros.price;
by year month day;
run;
proc sort data=frame_stock1;
by year month day;
run;
data frame_price;
merge frame_stock1 pros.price;
by year month day;
run;
/***************************************************************************************/
/***************merge with wave height and wave speed***********************************/
/***************************************************************************************/
data weather;
set pros.weather;
where year>=2000;
run;

proc sort data=weather;
by year month day;
run;
proc sort data=frame_price;
by year month day;
run;
data pros.firststage;
merge frame_price weather;
by year month day;
run;

data pros.firststage;
set pros.firststage;
where year>=2000;
if year2000=. then delete;
if year=2000 and y365<170 and y365>121 then opening1=0;
else if year=2000 and (y365>169 or y365<122) then  opening1=1;
if year=2001 and y365<153 and y365>139 then opening1=0;
else if year=2001 and (y365>152 or y365<140) then  opening1=1; 
if year=2002 and y365<153 and y365>139 then opening1=0;
else if year=2002 and (y365>152 or y365<140) then  opening1=1;
if year=2003 and y365<155 and y365>121 then opening1=0;
else if year=2003 and (y365>154 or y365<122) then  opening1=1;
if year=2004 and y365<152 and y365>121 then opening1=0;
else if year=2004 and (y365>151 or y365<122) then  opening1=1;
if year=2005 and y365<163 and y365>121 then opening1=0;
else if year=2005 and (y365>162 or y365<122) then  opening1=1; 
if month=1 & day<10 then stock=10;
run;

Proc sql;
create table gh as
select distinct year,sum(pounds_a*prize_a) as revenue,sum(pounds_a) as harvest,mean(prize_a) as price
from pros.firststage
group by year;
quit;
