libname pros 'C:\Documents and Settings\LHuang\Desktop\code\pros';

/***************************************************************************************/
/***************logit mdoel*************************************************************/
/***************************************************************************************/
data pros.firststage1;
set pros.firststage;
where year<=2004;
run;

/**result**/

proc logistic data=pros.firststage1;
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
