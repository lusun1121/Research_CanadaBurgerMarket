clear;
tic;
load indiv.txt;%vesselid len vess vess_fixed
load daily0.txt;%year,y365,stock, prize_a, wspd,wvht, holiday,8week,opening,totalhar
load shrimp_price0.txt; %price predicted
load diesel_price0.txt; %price predicted
load zeta.txt;

    year=1;
    timess=100;
    dayt=365;
    
daily=daily0(365*(year-1)+1:365*year,:);
shrimp_price=shrimp_price0(365*(year-1)+1:365*year,:);
diesel_price=diesel_price0(365*(year-1)+1:365*year,:); 
zetan=zeta(365*(year-1)+1:365*year,3);
indiv(:,5)=exp(indiv(:,4));
    mu = [0 0];
    sigma = [0.129154 0.071702; 0.071702 0.119709];
    dis=0.9998;
    vvn=439;

    vv=1:vvn;
    v_shrimp=var(shrimp_price0(:,1)-shrimp_price0(:,2));
    v_diesel=var(diesel_price0(:,1)-diesel_price0(:,2));
    EU=0.577215665;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%%%%%%%%%%%%stock dynamics coefficients%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stock_year=[0.577319625	0.066830688	0.689766452	0.381456417	0.080095907	0]; %%stock dynamics yearly dummy
    stock_parm=[9.854831206	-7.394586282	1.796771306	-0.137313124	-0.000000127];%%stock dynamics coefficients
    catch_parm=[-0.001270945	0.717147271	-0.00407091	0.022396933];%%Tvessel opening wspd wvht        
    gamma=catch_parm(1);

load ind_parm2.csv;%%1prize_a*vess  stock*vess	diesel_price*vess holiday*vess	

ind_parm=ind_parm2;
%%%probility parameter recovered from the first stage
% yearch=-2.987617892;
% probit_par=[0.110795653	0.068141184	-0.008592882	-0.829442837	-0.034409959	...
%     0.001155369	0.007109313	-0.000629725	0.000003588	-0.000222315	-0.000015618	...
%     -0.000001251	0.000000003	0.096645103	0.000441607	-0.079830256	-0.00048004	...
%     0.00001562	-0.190420152	-0.018729955	0.001475651	0.235449612	...
%     0.016413731	0.137615939];
%%revenue Wind Wave Diesel weekend Diesel2  stock stock2 len len2
%%len_wspd len_wvht len_disel vess_fixed vess_fixed2 
%%har2,har,har_len,har2_len,year*****calibrated*******,opening

% struc=[1.468587/1000	-0.122402	-0.301329	-0.020102	...
%   -0.839072	0 0.062601	0 -0.317129/10	-0.038309/100	...
%     0.123840/100	0.352540/100	0.024288/100 0 0 	-5.374121/1000000 	...
%    -1.432597/1000	0.664760/10000 0 0 0];

% yearch=-19.7419;%%CCP constant
% probit_par=[3.0774	0.111	0.1534	-0.9464	-0.6834	0.0306	-0.0151	0.00223	-0.00005	...
%     -0.00277	-9.24E-06	-0.00077	1.20E-06	0.6656	0.000837	-0.0958	0.0108...
%     -0.00037	0.0839	-0.1148	0.00913	0.6886	-0.0169	-0.0179	0.1655];
% 
% 
% struc=[1.438065/1000	-0.130760	-0.285467	-0.019618	...
%   -0.840191	0 0.060645	0 -0.304766/10	-0.038391/100	...
%     0.130252/100	0.322302/100	0.02384/100 0 0 		-5.185154/1000000 	...
%    -1.613786/1000	0.67781/10000 0 0 0];


yearch=-2.987941963;
probit_par=[0.117414027	0.063966915	-0.007352238	-0.830043794	-0.033051607	...
    0.001088292	0.00685054	-0.000565685	0.000002878	-0.000262922	...
    -0.000010657	-0.000005409	0.000000009	0.095935324	0.000443973	...
    -0.094368799	0.00044001	-0.000005967	-0.183506342	-0.021570058	...
    0.001685203	0.222659708	0.017507882	0.133714599	0.002475202];


struc=[1.433641/1000	-0.130195	-0.291469	-0.019565	...
  -0.83758	0 0.060517	0 -0.307893/10	-0.038134/100	...
    0.130261/100	0.329719/100	0.023889/100 0 0 		-5.163624/1000000 	...
   -1.61147/1000	0.67486/10000 0 0 0];





ffix=exp(indiv(:,4));%vess_fixed
valueo1=zeros(vvn,1);
revo1=zeros(vvn,1);

state_daily(:,1)=shrimp_price(:,1);%prize_a
state_daily(:,2:3)=daily(:,5:6);%wspd,wvht
state_daily(:,4)=diesel_price(:,1);%%diesel price


% for time=2:dayt;
%     state_daily(time,1)= shrimp_price(time-1,2)+sqrt(v_shrimp)*randn;   %%%shrimp price
%     r = mvnrnd(mu,sigma,2);
%     state_daily(time,2)=exp(1.088944+0.406951*log(state_daily(time-1,2))+r(1)); %%%wind speed
%     state_daily(time,3)=exp(0.039297+0.693213*log(state_daily(time-1,3))+r(2)); %%%wave height
%     state_daily(time,4)=diesel_price(time-1,2)+sqrt(v_diesel)*randn;       %%%diesel price
% end;
   
for rep=1:timess; 
state_daily0(:,1)=ones(vvn,1)*daily(1,3);%stock
state_daily0(:,2)=ones(vvn,1)*daily(1,10);%stock



accum=0;

stock_fin(1)=10;
stock_fin(26)=daily(26,3);
valueo=zeros(vvn,1);  
revo=zeros(vvn,1);





 for tt=1:dayt
     probitt=yearch+ind_parm(:,1)+probit_par(1)*state_daily(tt,1)+probit_par(2)*stock_fin(tt)+probit_par(3)*state_daily(tt,4) ...
     +probit_par(4)*daily(tt,7)+probit_par(5)*state_daily(tt,1)^2+probit_par(6)*state_daily(tt,1)^3  ...
     +probit_par(7)*state_daily(tt,1)*stock_fin(tt)+probit_par(8)*(stock_fin(tt)).^2+probit_par(9)*(stock_fin(tt)).^3 ...
     +probit_par(10)*state_daily(tt,1)*state_daily(tt,4)+probit_par(11)*state_daily(tt,4)*stock_fin(tt)+probit_par(12)*state_daily(tt,4)^2 ...
     +probit_par(13)*state_daily(tt,4)^3+probit_par(14)*state_daily(tt,1)*indiv(:,4)+probit_par(15)*stock_fin(tt).*indiv(:,2) ...
     +probit_par(16)*state_daily(tt,2)+probit_par(17)*state_daily(tt,2)^2+probit_par(18)*state_daily(tt,2)^3 ...
     +probit_par(19)*state_daily(tt,3)+probit_par(20)*state_daily(tt,3)^2+probit_par(21)*state_daily(tt,3)^3 ...
     +probit_par(22)*daily(tt,9)+probit_par(23)*state_daily(tt,2)*indiv(:,4)+probit_par(24)*state_daily(tt,3)*indiv(:,4)+probit_par(25)*exp(zetan(tt,1));

    
    prob=exp(probitt)./(1+exp(probitt));
    Tvessel(tt)=sum(prob);    
    
    share1=sum(prob)-prob;
    share_both=(sum(prob)^2-sum(prob.^2))/2-prob.*(share1);
    share_ex=1+gamma+gamma^2/2+(gamma+3*gamma^2/2)*share1+gamma^2/2*share_both;
    har=exp(catch_parm(3)*state_daily(tt,2)+catch_parm(4)*state_daily(tt,3)+catch_parm(2)*daily(tt,9))*share_ex.*ffix.*stock_fin(tt);

    har_real=har;    
    har2_real=har_real.^2;
    
    accum=accum+sum(har_real.*prob);    
    days=log(tt);
    days_m=[days days^2 days^3 days^4];
    stock_fin(tt+1)=exp(stock_year(year)+stock_parm(1:4)*days_m'+stock_parm(5)*accum);  
     if tt<10
        stock_fin(tt+1)=10;
     end
    zetan(tt+1,1)=zetan(tt,1)*(0.72201)+mvnrnd(0,0.19225,1); %%%%zeta
    harvest(tt)=sum(har_real.*prob);
   % state_daily0(tt,2)=sum(har_real);  
    
        %%%get accumulated sum;
    value_day(:,1)=state_daily(tt,1)*har_real;%%revenue
    value_day(:,2:6)= ones(vvn,1)*[state_daily(tt,2:4),daily(tt,7),state_daily(tt,4)^2];
    %%Wind Wave Diesel weekend Diesel2
    value_day(:,7:15)=(ones(vvn,9)).*[stock_fin(tt)*ones(vvn,1), stock_fin(tt)^2*ones(vvn,1), ...
        indiv(:,2), indiv(:,2).^2, indiv(:,2)*state_daily(tt,2:4), indiv(:,4), indiv(:,4).^2];
    %%stock stock2 len len2 len_wspd len_wvht len_disel vess_fixed vess_fixed2
    value_day(:,16:21)=[har2_real,har_real,har_real.*indiv(:,2), har2_real.*indiv(:,2),ones(vvn,1), ones(vvn,1)*daily(tt,9)];
    %%har2,har,har_len,har2_len,year,opening
   
    
    err_r= (-log(prob));
     
    valueo=valueo+dis^(tt-1)*prob.*(value_day*struc'+err_r);
    revo=revo+dis^(tt-1)*prob.*(value_day(:,1)*struc(1)+err_r);
    

    end
    valueo1=valueo1+1/timess*valueo;
    revo1=revo1+1/timess*revo;
end; 

adj_p_total=sum(valueo1)/struc(1)
harvest_dis=0;
for ti=1:dayt;
harvest_dis=harvest_dis+daily(ti,10)*shrimp_price(ti,2)*dis^(ti-1);
end;
actural=harvest_dis
har_totrl=sum(harvest)
toc;
figure(4)
plot(1:dayt,Tvessel,1:dayt,daily(1:dayt,11))
legend('Simulated','Actual')
xlabel('Year 2000')
ylabel('Fishing vessel number')


figure(5)
plot(1:dayt,daily(1:dayt,end-1))
xlabel('Year 2000')
ylabel('Daily harvest (pound)')


dlmwrite('Tvessel.txt', Tvessel');
dlmwrite('harvest.txt', harvest');
dlmwrite('stock_fin.csv', stock_fin);

t=1:dayt;
reg=[Tvessel',daily(t,end),harvest',daily(t,end-1)];
dlmwrite('reg.txt', reg);


hh0=[];
hh=[daily(:,end) Tvessel' daily(:,end-1)];
for h1=1:dayt
    if hh(h1,1)~=0;
        hh0=[hh0; hh(h1,:)];
    end
end

[o1 o2]=size(hh0);
weight=hh0(1:o1,3)/sum(hh0(1:o1,3));
abdis=abs(hh0(1:o1,1)-hh0(1:o1,2))./hh0(1:o1,1);
percentage_error1=sum(abdis.*weight)
t=1:dayt;
percentage_error2=sum(abs((daily(t,end-1)-harvest'))/sum(daily(t,end-1)))



