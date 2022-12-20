libname pros 'C:\Documents and Settings\LHuang\Desktop\code\pros';


/***************************************************************************************/
/***************logit mdoel*************************************************************/
/***************************************************************************************/

proc logistic data=pros.firststage;
class vess;
model fishing0=
vess  prize_a   stock  
diesel_price holiday prize_a*prize_a prize_a*prize_a*prize_a prize_a*stock 
stock*stock stock*stock*stock
prize_a*diesel_price stock*diesel_price diesel_price*diesel_price diesel_price*diesel_price*diesel_price
prize_a*vess_fixed stock*len wspd wspd*wspd wspd*wspd*wspd wvht wvht*wvht wvht*wvht*wvht
 opening1 
vess_fixed*wspd  vess_fixed*wvht zeta1
 / link=logit;
 OUTPUT OUT=hh Pred=p;
run;



/*****produce data total catch and trips for submission*******/


proc sql;
create table catch_trip_aggregate as
select distinct year, month, day, sum(pounds_a) as catch_total, sum(fishing) as trips_total
from pros.firststage
group by year, month, day;
quit;

data catch_trip_aggregate1;
set catch_trip_aggregate;
if catch_total<200 & catch_total>0 then catch_scrubbed="< 200";
else catch_scrubbed=catch_total;
if trips_total<4 & trips_total>0 then trips_scrubbed="< 4";
else trips_scrubbed=trips_total;
run;



PROC EXPORT 
data=catch_trip_aggregate1
OUTFILE='C:\Documents and Settings\LHuang\Desktop\catch_trip_aggregate.csv'
replace;
quit;


goptions reset=all gunit=pct border cback=white;
	 symbol1 v=dot w=1;
	 symbol2 interpol=join w=5;
proc Gplot data=hh1;
plot ab*y365;
run;
quit;

/*proc sql;
create table hh1 as
select distinct year, month, day, sum(p) as prd1
from hh
group by year, month, day;
quit;*/




/***************************************************************************************/
/***************shimp price evolution***************************************************/
/***************************************************************************************/
proc sql;
create table pros.price1 as
select distinct year,month,day,time,prize_a,week
from pros.firststage;
quit;

proc sort data=pros.price1;
by year month day;
run;
data price1;
set pros.price1;
if prize_a>5 then prize_a=5;
run;

proc transreg data=price1;
    model identity(prize_a) = bspline(time/ nknots=60);
    output out=price2 p;
run;

proc reg data=price2;
 model prize_a= time_0-time_62;
 run;


proc sort data=price2;
by time;
run;

data daily_price;
set price2;
format time date9.;
if time='29feb2000'd then delete;
if time='29feb2004'd then delete;
if prize_a=. then prize_a=1.954;
keep prize_a Pprize_a;
run;

PROC EXPORT DATA=daily_price 
OUTFILE='I:\Ling\Prospectus\NPL01232012\6simulation\shrimp_price0.txt'
replace;
quit;




/*data price3;
   set price2;
   format time date9.;
run;

goptions reset=all gunit=pct border cback=white;
	 symbol1 v=dot w=1;
	 symbol2 interpol=join w=5;
axis1 order=('01Jan2000'd to '31Dec2005'd by year);
Proc GPLOT data=price3;
     Plot prize_a*time Pprize_a*time/overlay haxis=axis1;
	 run;
quit;*/

/***************************************************************************************/
/***************diesel price evolution**************************************************/
/***************************************************************************************/
data diesel0;
set pros.diesel1;
week1=week+52*(year-2000);
run;

proc sort data=diesel0;
by week1;
run;

proc transreg ss2 data=diesel0;
     model identity(diesel_price)=bspline(week1/nknots=30);
     output out=diesel1 p;
run;




proc sql;
create table d1 as
select distinct year,month,day,week
from pros.firststage;
quit;
data d2;
set d1;
week1=week+52*(year-2000);
run;
proc sort data=d2;
by week1;
run;
proc sort data=diesel1;
by week1;
run;
data d3;
merge d2 diesel1;
by week1;
run;
proc sort data=d3;
by year month day;
run;
data d4;
set d3;
if month=2 & day=29 then delete;
keep diesel_price Pdiesel_price;
run;
PROC EXPORT DATA=d4 
OUTFILE='I:\Ling\Prospectus\NPL01232012\6simulation\diesel_price0.txt'
replace;
quit;

/*goptions reset=all gunit=pct border cback=white;
SYMBOL1 v=dot i=none w=1 C=blue;      
SYMBOL2 i=join w=2 c=red;
axis1 order=(1 to 321 by 20);
Proc GPLOT data=diesel1;
     Plot diesel_price*week1=1 Pdiesel_price*week1=2/overlay haxis=axis1; 
	 run;
quit;*/

/***************************************************************************************/
/***************weather evolution*******************************************************/
/***************************************************************************************/
proc sort data=pros.weather;
by year month day;
run;

data station_lag;
set pros.weather;
wspd_l=log(wspd);
wspd_l1=lag(wspd_l);
wspd_l2=lag(wspd_l1);
wvht_l=log(wvht);
wvht_l1=lag(wvht_l);
wvht_l2=lag(wvht_l1);
run;

proc arima data=station_lag;
identify var=wspd_l nlag=1 stationarity=(adf=(1));
identify var=wvht_l nlag=1 stationarity=(adf=(1));
run;
quit;

/*wind speed evolution*/
/*wave height evolution*/
proc syslin data =station_lag sur ;
  wspd_l: model wspd_l=wspd_l1;
  wvht_l: model wvht_l=wvht_l1;
run;
