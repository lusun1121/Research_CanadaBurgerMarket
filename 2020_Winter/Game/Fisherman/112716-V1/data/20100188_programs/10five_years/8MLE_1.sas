libname pros 'C:\Documents and Settings\LHuang\Desktop\code\pros';



/*import the data*/
proc import datafile='C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\five_years\7secondstage_1\final_all1.txt'
out=fiinal_all1  dbms=dlm
replace;
delimiter=',';
getnames=no;
run;
quit;
proc import datafile='C:\Documents and Settings\LHuang\Desktop\code\SNPL2_1\five_years\7secondstage_1\final_all2.txt'
out=fiinal_all2  dbms=dlm
replace;
delimiter=',';
getnames=no;
run;
quit;

/*1year 2year 3t 4vess 5fishing 6revenue 7Wind 8Wave 9Diesel 10weekend Diesel2 stock stock2 14len len2*/

data final_all3;
set fiinal_all1;
fishing0=1-var5;
fishing1=var5;
revnue=var6/1000;
wspd=var7;
wvht=var8;
diesel=var9;
weekend=var10;
diesel2=var11;
stock1=var12;
stock2=var13;
len=var14/10;
len2=var15/100;

drop var5-var14;
if var1=2000 then year2000=1;else year2000=0;
if var1=2001 then year2001=1;else year2001=0;
if var1=2002 then year2002=1;else year2002=0;
if var1=2003 then year2003=1;else year2003=0;
if var1=2004 then year2004=1;else year2004=0;
if var1=2005 then year2005=1;else year2005=0;
run;
/*1len_wspd 2len_wvht 3len_disel 4vess_fixed 5vess_fixed2
  61error 7har2,8har,9har_len,10har2_len,11year,12opening;*/
data final_all4;
set fiinal_all2;
zeta1=var14;
len_wspd=var1/100;
len_wvht=var2/100;
len_diesel=var3/100;
vess_fixed=var4;

vess_fixed2=var5;
error=var6;
har2=var7/1000000;
har=var8/1000;
har_len=var9/10000;
har2_len=var10/10000;
opening=var12;
drop var1-var12;
run;
data pros.final_all;
merge final_all3 final_all4;
run;
data final_all1;
set pros.final_all;
where var1<=2004;
run;


/*year y365 vess fishing revenue Wind Wave Diesel weekend Diesel2 stock stock2 len len2 
 error har2,vess,vess2,har,har_len,har2_len,har3;*/
/*proc logistic data=pros.final_all;
model fishing0=
revnue wspd wvht diesel  weekend
stock1 stock2 len len2  len_wspd  len_wvht len_diesel har har_len har2 error / link=probit;
run;

proc logistic data=pros.final_all;
model fishing0=
 revnue wspd wvht diesel  weekend
stock1 stock2 len  len_wspd  len_wvht len_diesel har har_len har2 / link=logit;
run;*/

/***************************************************************************************/
/****************************************Result*****************************************/
/***************************************************************************************/


proc qlim data=final_all1;
model fishing1= revnue wspd wvht diesel  weekend
stock1  len len2  len_wspd  len_wvht len_diesel har2 har har_len  error
/NOINT Discrete(dist=logit);
restrict error=1;
OUTPUT OUT=pros.prd PROBALL;
run;

/***************************************************************************************/
/***************Iterations**************************************************************/
/***************************************************************************************/

data prd;
set pros.prd;
year=var1;
y365=var3;
vess=var4;
logi=log(Prob2_fishing1/(1-Prob2_fishing1));
keep year y365 vess logi Prob2_fishing1;
run;

proc sort data=prd;
by year y365 vess;
run;

proc sort data=pros.firststage;
by year y365 vess;
run;

data prd1;
merge pros.firststage prd;
by year y365 vess;
run;

proc glm data=prd1;
class vess;
model logi=
  prize_a   stock  
diesel_price holiday prize_a*prize_a prize_a*prize_a*prize_a prize_a*stock 
stock*stock stock*stock*stock
prize_a*diesel_price stock*diesel_price diesel_price*diesel_price diesel_price*diesel_price*diesel_price
prize_a*vess_fixed stock*len wspd wspd*wspd wspd*wspd*wspd wvht wvht*wvht wvht*wvht*wvht
 opening1 
vess_fixed*wspd  vess_fixed*wvht zeta1 vess /solution e;
run;
quit;
