libname pros 'C:\Documents and Settings\LHuang\Desktop\code\pros';
/***************************************************************************************/
/***************single day trip only****************************************************/
/***************************************************************************************/
proc sql;
create table stockn as
select distinct *, started as started_back
from pros.active
where active1>=1;
quit;

data fish;
set stockn;
run;

%macro dup;
%do i=2 %to 10;
   data stock&i;
   set stockn;
   where active1>=&i;
   started=started-1+&i;
   run;
   proc append base=fish  data=stock&i;
   run;
%end;
%mend;
%dup;

data fish_n;
set fish;
fishing=1;
run;
/***************************************************************************************/
/***************no overlap trip ********************************************************/
/***************************************************************************************/
Proc sql;
create table fish_overlap1 as
select distinct *,sum(fishing) as overlap
from fish_n
group by vesselid, started;
quit;

proc sort data=fish_overlap1;
by vesselid started;
run;

DATA fish_overlap2;
    set fish_overlap1;
	by vesselid started;
 overlap1=0; 
 if overlap=1 then overlap1=1;
 if first.started & overlap>1 then overlap1=1; 
run;


Proc sql;
create table fish_overlap3 as
select distinct *,sum(overlap1) as active2
from fish_overlap2
group by vesselid, started_back;
quit;


data pros.fish;
set fish_overlap3;
if  overlap1=0 then delete;
year=year(started);
month=month(started);
day=day(started);
pounds_a=poundss/active2;
run;


/***************************************************************************************/
/***************Tvessel*****************************************************************/
/***************************************************************************************/

proc sql;
create table fish1 as
select distinct *, sum(fishing) as Tvessel
from pros.fish
group by year, month, day;
quit;




/***************************************************************************************/
/***************id vessel***************************************************************/
/***************************************************************************************/
proc sql;
create table vessel as
select distinct vesselid, len
from fish1;
quit;

data vessel1;
 do vess=1 to 439;
   output; 
 end;
 run;
data pros.vessel;
merge vessel vessel1;
run;
proc sort data=fish1;
by vesselid;
run;
proc sort data=pros.vessel;
by vesselid;
run;
data fishh2;
merge pros.vessel fish1;
by vesselid;
run;
/***************************************************************************************/
/***************id time*****************************************************************/
/***************************************************************************************/
proc sort data=fishh2;
by month day;
run;
proc sort data=pros.y365;
by month day;
run;
data fishh3;
merge pros.y365 fishh2;
by month day;
run;



data fishh4;
set fishh3;
where year>=2000;
pounds_l=log(pounds_a);
if weekday(started)=6 or weekday(started)=7 then holiday=1;
else holiday=0;
if month=. then delete;
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
if pounds_a=. then poundss_a=0;
run;

/***************************************************************************************/
/***************accumulated harvest*****************************************************/
/***************************************************************************************/
proc sql;
create table accu as
select distinct year, month, day, sum(pounds_a) as poun
from fishh4
group by year, month, day;
quit;

proc sort data=accu;
by year month day;
run;

DATA pros.accumulated_har;
    set accu;
	by year;
 if first.year then poundss_sum=0;
    poundss_sum+poun;
run;


proc sort data=pros.accumulated_har;
by year month day;
run;

proc sort data=fishh4;
by year month day;
run;

data fishh5;
merge fishh4 pros.accumulated_har;
by year month day;
run;


/***************************************************************************************/
/***************merge with wave height and wave speed***********************************/
/***************************************************************************************/

proc sort data=pros.weather;
by year month day;
run;

PROC EXPORT 
data=pros.weather
OUTFILE='C:\Documents and Settings\LHuang\Desktop\weather.csv'
replace;
quit;


/*proc sql;
create table weather1 as
select distinct *, mean(WSPD) as WSPD_m, mean(WVHT) as WVHT_m
from pros.weather
group by year, month;
quit;*/


proc sort data=fishh5;
by year month day;
run;
data stock;
merge fishh5 pros.weather;
by year month day;
run;



data pros.stock;
set stock_n;
where year>=2000;
yt=log(y365);
yt2=yt*yt;
yt3=yt2*yt;
yt4=yt3*yt;
run;

/***************************************************************************************/
/***************regression of schafer model*********************************************/
/***************************************************************************************/

PROC GLM DATA=pros.stock;
CLASS year vess;
  model pounds_l= vess year yt yt2 yt3 yt4 Tvessel poundss_sum opening1 WSPD WVHT/noint SOLUTION;
run;
quit;


PROC GLM DATA=pros.stock;
CLASS year vess y365;
  model pounds_l= vess year yt y365 Tvessel poundss_sum opening1 WSPD WVHT/noint SOLUTION;
run;
quit;



/***************************************************************************************/
/***************IV regression of schafer model******************************************/
/***************************************************************************************/

PROC GLM DATA=pros.stock;
CLASS year  vess;
MODEL Tvessel= vess year yt yt2 yt3 yt4 holiday poundss_sum opening1 WSPD WVHT/noint SOLUTION;
  output out=stock_iv p=IVTvessel;
run;
quit;

PROC GLM DATA=stock_iv;
CLASS  year y365 vess;
  model pounds_l= vess year yt yt2 yt3 yt4 IVTvessel poundss_sum opening1 WSPD WVHT/noint SOLUTION;
  output out=stock_iv1 r=residual;
run;
quit;

data stock_iv2;
set stock_iv1;
resi=residual-IVTvessel*(-0.001270944)+Tvessel*(-0.001270944);
run;

PROC GLM DATA=stock_iv2 ;
CLASS year y365;
MODEL resi= year*y365/noint SOLUTION;
  output out=stock_iv3 r=residual1;
run;
quit;

PROC UNIVARIATE Data=stock_iv3  NOPRINT;
VAR residual1;
HISTOGRAM / NORMAL (COLOR=RED);
RUN; 

PROC UNIVARIATE Data=stock_iv3 NORMALTEST;
VAR residual1;
RUN; 

proc corr Data=stock_iv3 sscp cov plots=matrix;
var residual1 Tvessel;
run;

PROC reg DATA=stock_iv3;
MODEL residual1= Tvessel/noint;
run;
quit;

proc import datafile="C:\Documents and Settings\LHuang\Desktop\code\SNPL2\zeta.csv"
out=zeta
replace;
run;

data zeta_frame;
    do year=2000 to 2005;
	   do y365=1 to 365;
		   output;
	   end;
	end;
run;

proc sort data=zeta;
by year y365;
run;

data zeta1;
merge zeta_frame zeta;
by year y365;
run;

proc sort data=zeta1;
by year y365;
run;

data zeta2;
set zeta1;
zeta_l=lag(zeta);
run;
/*result*/

proc reg data=zeta2;
model zeta=zeta_l/noint;
	  where y365>150 and y365<300;
run;
quit;

proc reg data=zeta2;
model zeta=zeta_l;
run;
quit;


data pros.zeta;
set zeta2;
if zeta=. then zeta=0;
run;


/***************************************************************************************/
/*************Manually create csv of fixed vess*****************************************/
/***************************************************************************************/
PROC IMPORT 
DATAFILE='C:\Documents and Settings\LHuang\Desktop\code\SNPL2\vess_fixed.csv'
OUT=pros.vess_fixe
replace;
quit;
proc sort data=pros.vessel;
by vess;
run;
proc univariate data=pros.vessel;
      var len;
      histogram;
   run;

proc sort data=pros.vess_fixe;
by vess;
run;
data pros.vess_fixed;
merge pros.vessel pros.vess_fixe;
by vess;
run;
