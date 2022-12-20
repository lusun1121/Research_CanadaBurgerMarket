% Copyright Ling Huang, Jan. 2009
clear all;

load indiv.txt;%vesselid len vess vess_fixed
load daily0.txt;%year,y365,stock, prize_a, wspd,wvht, holiday,9week,opening,totalhar
load shrimp_price0.txt; %price predicted
load diesel_price0.txt; %price predicted
load ind_parm2.csv;%%1prize_a*vess  stock*vess	diesel_price*vess holiday*vess	
load zeta.txt;
ind_parm=ind_parm2;
tic;

year=1; %year number
maxx=[80 50 92 70 50 50];
scc=[1/20 1/40 1/20 1/30 1/40 1/40];
sto_max=maxx(year);

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
    EU=0.577215665;
    v_shrimp=var(shrimp_price0(:,1)-shrimp_price0(:,2));
    v_diesel=var(diesel_price0(:,1)-diesel_price0(:,2));
    
    stock_year=[0.577319625	0.066830688	0.689766452	0.381456417	0.080095907	0]; %%stock dynamics yearly dummy
    stock_parm=[9.854831206	-7.394586282	1.796771306	-0.137313124	-0.000000127];%%stock dynamics coefficients
    catch_parm=[-0.001270945	0.717147271	-0.00407091	0.022396933];%%Tvessel opening wspd wvht           
    gamma=catch_parm(1);
%%%revenue Wind Wave Diesel weekend Diesel2  stock stock2 len len2
%%%len_wspd len_wvht len_disel vess_fixed vess_fixed2 
%%%har2,har,har_len,har2_len,year*****calibrated*******,opening error
% yearch=-2.987617892;
% probit_par=[0.110795653	0.068141184	-0.008592882	-0.829442837	-0.034409959	...
%     0.001155369	0.007109313	-0.000629725	0.000003588	-0.000222315	-0.000015618	...
%     -0.000001251	0.000000003	0.096645103	0.000441607	-0.079830256	-0.00048004	...
%     0.00001562	-0.190420152	-0.018729955	0.001475651	0.235449612	...
%     0.016413731	0.137615939];

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

yearch=-2.825799276;
probit_par=[0.093694825	0.066400077	-0.009232206	-0.832516647	-0.028709727	0.000980912	...
    0.005530288	-0.00055852	0.0000031	-0.000133724	-0.000006491	-0.000002561	...
    0.000000003	0.069997019	0.000473012	-0.111566316	0.00018671	-0.000004758	...
    -0.131226086	-0.012175392	0.001007597	0.240396609	0.018542828	...
    0.136400365	0.002096766];

struc=[1.140367/1000	-0.151846	-0.235786	-0.019224	...
  -0.842689	0 0.061773	0 -0.299459/10	-0.034591/100	...
    0.13124/100	0.398136/100	0.020482/100 0 0 	-4.649996/1000000 	...
  -0.411332/1000	0.552502/10000 0 0 0];

ffix=exp(indiv(:,4));%vess_fixed

%%revenue Wind Wave Diesel weekend Diesel2 stock stock2 len len2 len_wspd 
%%%len_wvht len_disel vess_fixed vess_fixed2 har2 har har*len har2_len one opening

tic;
%scale1=1/10;
%scale2=1;
scale1=scc(year); %stock
scale2=5;  %vessel number
scale3=20;  %days


%  scale1=0.5*10^5;
%  scale2=5;
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
mapping=zeros(round(sto_max/scale1)+1,floor(dayt/scale3)*3+3);

mt=floor(dayt/scale3);

for j=1:mt-1; 
    for stock1=1:scale1:sto_max; 
        for I=0:scale2:vvn;
            stock=stock1;
            ppro=0;
                zetan1=zetan(dayt-j*scale3+1,1);
            for ttt=1:scale3

    time=dayt-j*scale3+ttt;

% %%1prize_a*prize_a stock*stock 3prize_a*stock prize_a*diesel_price 5prize_a*vess_fixed
% %%stock*len 7WSPD 8WVHT 9opening1 10diesel_pr*diesel_pri 11len 12len*len 13vess_fixed
% %%14vess_fixed*WSPD 15vess_fixed*WVHT 16vess_fixed*vess_fixed
      
     probitt=yearch+ind_parm(:,1)+probit_par(1)*state_daily(time,1)+probit_par(3)*state_daily(time,4) ...
     +probit_par(4)*daily(time,7)+probit_par(5)*state_daily(time,1)^2+probit_par(6)*state_daily(time,1)^3  ...
     +probit_par(10)*state_daily(time,1)*state_daily(time,4)+probit_par(12)*state_daily(time,4)^2 ...
     +probit_par(13)*state_daily(time,4)^3+probit_par(14)*state_daily(time,1)*indiv(:,4)+ ...
     +probit_par(16)*state_daily(time,2)+probit_par(17)*state_daily(time,2)^2+probit_par(18)*state_daily(time,2)^3 ...
     +probit_par(19)*state_daily(time,3)+probit_par(20)*state_daily(time,3)^2+probit_par(21)*state_daily(time,3)^3 ...
     +probit_par(22)*daily(time,9)+probit_par(23)*state_daily(time,2)*indiv(:,4) ...
     +probit_par(24)*state_daily(time,3)*indiv(:,4)+probit_par(25)*exp(zetan1)+probit_par(2)*stock  ...
     +probit_par(7)*state_daily(time,1)*stock+probit_par(8)*(stock).^2+probit_par(9)*(stock).^3 ...
     +probit_par(11)*state_daily(time,4)*stock+probit_par(15)*stock.*indiv(:,2);
      
     prob=exp(probitt)./(1+exp(probitt));
     share1=sum(prob)-prob;
    share_both=(sum(prob)^2-sum(prob.^2))/2-prob.*(share1);
    share_ex=1+gamma+gamma^2/2+(gamma+3*gamma^2/2)*share1+gamma^2/2*share_both;
    har_new_1=exp(catch_parm(3)*state_daily(time,2)+catch_parm(4)*state_daily(time,3) ...
        +catch_parm(2)*daily(time,9))*share_ex.*ffix.*stock;

      har2_new_1=har_new_1.^2;
      har_real=har_new_1.*prob;
      err= (-log(prob));
      %%%get value;
      value0(:,1)=state_daily(time,1)*har_new_1;%%revenue
      value0(:,2:6)=ones(vvn,1)*[state_daily(time,2:4),daily(time,7),state_daily(time,4)^2];
      %%Wind Wave Diesel weekend Diesel2
      value0(:,7:15)=[stock*ones(vvn,1), stock^2*ones(vvn,1),indiv(:,2), indiv(:,2).^2, ...
         indiv(:,2)*state_daily(time,2:4),indiv(:,4), indiv(:,4).^2];
      %%stock stock2 len len2 len_wspd len_wvht len_disel vess_fixed vess_fixed2
      value0(:,16:21)=[har2_new_1,har_new_1,har_new_1.*indiv(:,2), har2_new_1.*indiv(:,2),  ...
      ones(vvn,1),ones(vvn,1)*daily(time,9)];
    
      pred_pro=prob.*(value0*struc'+err);
      pred_m=[(1:vvn)' pred_pro];
      pred_s=sortrows(pred_m,2);
          
      active=[zeros(vvn-I,1);ones(I,1)];                     
      pred_so=[active pred_s];
      pred_sor=sortrows(pred_so,2);
      active_id=pred_sor(:,1);   
        
      pred_pro1=pred_pro.*active_id;
      har_real_1=active_id.*har_real;
      
                    tt=time+1;     
                    days=log(tt);
                    days_m=[days days^2 days^3 days^4];
                    daysp=log(tt-1);
                    days_mp=[daysp daysp^2 daysp^3 daysp^4];
      stock2=stock*exp(stock_parm(1:4)*(days_m'-days_mp')+stock_parm(5)*sum(har_real_1));
      stock=stock2;
      zetan1=zetan1*(0.72201)+mvnrnd(0,0.19225,1); %%%%zeta
      ppro=ppro+dis^(time-1)*sum(pred_pro1);

      end;
stock_trim=round(stock/scale1)+1;
stock_trimn=stock_trim*(stock_trim<=size(mapping,1))+size(mapping,1)*(stock_trim>size(mapping,1));
final1(I/scale2+1,1)=ppro+mapping(stock_trimn,(mt-j+1)*3+3);      
final1(I/scale2+1,2)=I; 
    end;
    final2=sortrows(final1,1); 
    mapping(round(stock1/scale1)+1,(mt-j)*3+1)=stock1;
    mapping(round(stock1/scale1)+1,(mt-j)*3+2)=final2(end,2);
    mapping(round(stock1/scale1)+1,(mt-j)*3+3)=final2(end,1);
    end;
end;

for j=mt; 
    for stock1=1:scale1:sto_max; 

        for I=0:scale2:vvn;

            stock=stock1;
            ppro=0;
            
             zetan1=zetan(10,1);
            for ttt=10-dayt+j*scale3:scale3
    time=dayt-j*scale3+ttt;
% %%1prize_a*prize_a stock*stock 3prize_a*stock prize_a*diesel_price 5prize_a*vess_fixed
% %%stock*len 7WSPD 8WVHT 9opening1 10diesel_pr*diesel_pri 11len 12len*len 13vess_fixed
% %%14vess_fixed*WSPD 15vess_fixed*WVHT 16vess_fixed*vess_fixed
      
     probitt=yearch+ind_parm(:,1)+probit_par(1)*state_daily(time,1)+probit_par(3)*state_daily(time,4) ...
     +probit_par(4)*daily(time,7)+probit_par(5)*state_daily(time,1)^2+probit_par(6)*state_daily(time,1)^3  ...
     +probit_par(10)*state_daily(time,1)*state_daily(time,4)+probit_par(12)*state_daily(time,4)^2 ...
     +probit_par(13)*state_daily(time,4)^3+probit_par(14)*state_daily(time,1)*indiv(:,4)+ ...
     +probit_par(16)*state_daily(time,2)+probit_par(17)*state_daily(time,2)^2+probit_par(18)*state_daily(time,2)^3 ...
     +probit_par(19)*state_daily(time,3)+probit_par(20)*state_daily(time,3)^2+probit_par(21)*state_daily(time,3)^3 ...
     +probit_par(22)*daily(time,9)+probit_par(23)*state_daily(time,2)*indiv(:,4) ...
     +probit_par(24)*state_daily(time,3)*indiv(:,4)+probit_par(25)*exp(zetan1)+probit_par(2)*stock  ...
     +probit_par(7)*state_daily(time,1)*stock+probit_par(8)*(stock).^2+probit_par(9)*(stock).^3 ...
     +probit_par(11)*state_daily(time,4)*stock+probit_par(15)*stock.*indiv(:,2);
      
     prob=exp(probitt)./(1+exp(probitt));
     share1=sum(prob)-prob;
    share_both=(sum(prob)^2-sum(prob.^2))/2-prob.*(share1);
    share_ex=1+gamma+gamma^2/2+(gamma+3*gamma^2/2)*share1+gamma^2/2*share_both;
    har_new_1=exp(catch_parm(3)*state_daily(time,2)+catch_parm(4)*state_daily(time,3) ...
        +catch_parm(2)*daily(time,9))*share_ex.*ffix.*stock;

      har2_new_1=har_new_1.^2;
      har_real=har_new_1.*prob;
      err= (-log(prob));
      %%%get value;
      value0(:,1)=state_daily(time,1)*har_new_1;%%revenue
      value0(:,2:6)=ones(vvn,1)*[state_daily(time,2:4),daily(time,7),state_daily(time,4)^2];
      %%Wind Wave Diesel weekend Diesel2
      value0(:,7:15)=[stock*ones(vvn,1), stock^2*ones(vvn,1),indiv(:,2), indiv(:,2).^2, ...
         indiv(:,2)*state_daily(time,2:4),indiv(:,4), indiv(:,4).^2];
      %%stock stock2 len len2 len_wspd len_wvht len_disel vess_fixed vess_fixed2
      value0(:,16:21)=[har2_new_1,har_new_1,har_new_1.*indiv(:,2), har2_new_1.*indiv(:,2),  ...
      ones(vvn,1),ones(vvn,1)*daily(time,9)];
      %%har2 har har*len har2_len one opening
    
      pred_pro=prob.*(value0*struc'+err);
      pred_m=[(1:vvn)' pred_pro];
      pred_s=sortrows(pred_m,2);
          
      active=[zeros(vvn-I,1);ones(I,1)];                     
      pred_so=[active pred_s];
      pred_sor=sortrows(pred_so,2);
      active_id=pred_sor(:,1);   
        
      pred_pro1=pred_pro.*active_id;
      har_real_1=active_id.*har_real;
      
                    tt=time+1;     
                    days=log(tt);
                    days_m=[days days^2 days^3 days^4];
                    daysp=log(tt-1);
                    days_mp=[daysp daysp^2 daysp^3 daysp^4];
      stock2=stock*exp(stock_parm(1:4)*(days_m'-days_mp')+stock_parm(5)*sum(har_real_1));
      stock=stock2;
      zetan1=zetan1*(0.72201)+mvnrnd(0,0.19225,1); %%%%zeta
      ppro=ppro+dis^(time-1)*sum(pred_pro1);

      end;
stock_trim=round(stock/scale1)+1;
stock_trimn=stock_trim*(stock_trim<=size(mapping,1))+size(mapping,1)*(stock_trim>size(mapping,1));
final1(I/scale2+1,1)=ppro+mapping(stock_trimn,(mt-j+1)*3+3);      
final1(I/scale2+1,2)=I; 
    end;
    final2=sortrows(final1,1); 
    mapping(round(stock1/scale1)+1,(mt-j)*3+1)=stock1;
    mapping(round(stock1/scale1)+1,(mt-j)*3+2)=final2(end,2);
    mapping(round(stock1/scale1)+1,(mt-j)*3+3)=final2(end,1);
    end;
end;
toc;
dlmwrite('mapping.csv', mapping);



