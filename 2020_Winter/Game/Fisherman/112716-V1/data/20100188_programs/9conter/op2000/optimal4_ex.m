% Copyright Ling Huang, Jan. 2009
clear all;

load indiv.txt;%vesselid len vess vess_fixed
load daily0.txt;%year,y365,stock, prize_a, wspd,wvht, holiday,9week,opening,totalhar
load shrimp_price0.txt; %price predicted
load diesel_price0.txt; %price predicted
load ind_parm2.csv;%%1prize_a*vess  stock*vess	diesel_price*vess holiday*vess	
ind_parm=ind_parm2;
load mapping.csv;
load max_rev.txt;
load zeta.txt;

tic;
year=1; %year number
timess=100;
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
% 
% struc=[1.468587/1000	-0.122402	-0.301329	-0.020102	...
%   -0.839072	0 0.062601	0 -0.317129/10	-0.038309/100	...
%     0.123840/100	0.352540/100	0.024288/100 0 0 	-5.374121/1000000 	...
%    -1.432597/1000	0.664760/10000 0 0 0];


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
tt=10;
days=log(tt);
days_m=[days days^2 days^3 days^4];
stock_uni=exp(stock_year(year)+stock_parm(1:4)*(days_m'));

profit=0;
for rep=1:timess; 


%for restrict1=max_rev(end,2)
  %  for restrict2=max_rev(end,3)
   for restrict1=0
    for restrict2=0

mt=floor(dayt/scale3);
ppro=0;

stock=stock_uni;
kk(1:9)=stock_uni*ones(1,9);

for jj=1;
    j=mt-jj+1;
stock_trim=round(stock/scale1)+1;
stock_trimn=stock_trim*(stock_trim<=size(mapping,1))+size(mapping,1)*(stock_trim>size(mapping,1));
I=mapping(stock_trimn,(mt-j)*3+2);
%I =vvn;
zetan1=zetan(10,1);
    for ttt=10-dayt+j*scale3:scale3
    time=dayt-j*scale3+ttt;   
    venum(time)=I;

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
      
                    tt=time;     
                    days=log(tt);
                    days_m=[days days^2 days^3 days^4];
                    daysp=log(tt-1);
                    days_mp=[daysp daysp^2 daysp^3 daysp^4];
      stock2=stock*exp(stock_parm(1:4)*(days_m'-days_mp')+stock_parm(5)*sum(har_real_1));
      stock=stock2;
      if time<11
         stock=stock_uni;
      end
       zetan1=zetan1*(0.72201)+mvnrnd(0,0.19225,1); %%%%zeta     
      kk(time)=stock;
      harvest_op(time)=sum(har_real_1);
      ppro=ppro+dis^(time-1)*sum(pred_pro1);
      venumnn(time)=sum(prob.*active_id);
      end;
end;

%stock=daily(26,3);

for jj=2:mt;
    j=mt-jj+1;
stock_trim=round(stock/scale1)+1;
stock_trimn=stock_trim*(stock_trim<=size(mapping,1))+size(mapping,1)*(stock_trim>size(mapping,1));
I=mapping(stock_trimn,(mt-j)*3+2);
%I=vvn;
zetan1=zetan(dayt-j*scale3+1,1);
             for ttt=1:scale3
    time=dayt-j*scale3+ttt;   
   venum(time)=I;
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

      %% har2,har,har_len,har2_len,year,opening error

     pred_pro=prob.*(value0*struc'+err);
      pred_m=[(1:vvn)' pred_pro];
      pred_s=sortrows(pred_m,2);
          
      active=[zeros(vvn-I,1);ones(I,1)];                     
      pred_so=[active pred_s];
      pred_sor=sortrows(pred_so,2);
      active_id=pred_sor(:,1); 
      
      
     if time>=restrict1 &  time<=restrict2
           active_id=zeros(vvn,1);
      end
      
        
      pred_pro1=pred_pro.*active_id;
      har_real_1=active_id.*har_real;
      
                    tt=time;     
                    days=log(tt);
                    days_m=[days days^2 days^3 days^4];
                    daysp=log(tt-1);
                    days_mp=[daysp daysp^2 daysp^3 daysp^4];
      stock2=stock*exp(stock_parm(1:4)*(days_m'-days_mp')+stock_parm(5)*sum(har_real_1));
      stock=stock2;
      zetan1=zetan1*(0.72201)+mvnrnd(0,0.19225,1); %%%%zeta
      kk(time)=stock;
      harvest_op(time)=sum(har_real_1);
      ppro=ppro+dis^(time-1)*sum(pred_pro1);
      venumnn(time)=sum(prob.*active_id);

      end;
end;
hh=ppro/struc(1) ;

     end
   end
 
   profit=profit+1/timess*hh;
end;
profit
   
   
toc;


load Tvessel.txt;
load stock_fin.csv;
load harvest.txt;
stock_in(1:10)=stock_uni*ones(1,10);
  
%%stock without fishing
for tt=11:dayt         
                    days=log(tt);
                    days_m=[days days^2 days^3 days^4];
                    daysp=log(tt-1);
                    days_mp=[daysp daysp^2 daysp^3 daysp^4];
     stock_in(tt)=stock_in(tt-1)*exp(stock_parm(1:4)*(days_m'-days_mp'));
end

figure(5)
plot(1:dayt,venumnn'-Tvessel(1:dayt),'linewidth',1)
hold on;
plot(1:dayt, zeros(dayt,1), 'r','linewidth',2)
hold off
xlabel('Year 2000')
ylabel('Differnce of optimal and predicted vessel number')

figure(6)
subplot(1,2,1)
plot(1:dayt,harvest_op'-harvest(1:dayt),'linewidth',1)
hold on;
plot(1:dayt, zeros(dayt,1), 'r','linewidth',2)
hold off
xlabel('Year 2000')
ylabel('Differnce of optimal and predicted harvest')
subplot(1,2,2)
ttt=1:dayt;
plot(ttt,stock_fin(1:dayt),'b',ttt,kk(1:dayt),'r',ttt,stock_in(1:dayt),'g','linewidth',2)
%legend('Predicted','Optimal','Stock w/o harvest')
%axis([0 400 0 90])
xlabel('Year 2000')
ylabel('Stock index')


% Tvessel_trunk=[];
% Tvessel_trunk1=[];
% opTvessel_trunk=[];
% opTvessel_trunk1=[];
% har_tr=[];
% har_tr1=[];
% ophar_tr=[];
% ophar_tr1=[];
% sto_tr=[];
% sto_tr1=[];
% 
% for t=1:dayt;
%     if daily(t,7)==0;
%         Tvessel_trunk=[Tvessel_trunk Tvessel(t)];
%         opTvessel_trunk=[opTvessel_trunk venumnn(t)];
%         har_tr=[har_tr harvest(t)];
%         ophar_tr=[ophar_tr harvest_op(t)];
%     else 
%         Tvessel_trunk1=[Tvessel_trunk1 Tvessel(t)];
%         opTvessel_trunk1=[opTvessel_trunk1 venumnn(t)];
%         har_tr1=[har_tr1 harvest(t)];
%         ophar_tr1=[ophar_tr1 harvest_op(t)];
%     end
% end
% n=sum(1-daily(:,7));
% m=sum(daily(:,7));
% 
% figure(7)
% subplot(1,2,1)
% t=1:n;
% plot(t,opTvessel_trunk,'r', t, Tvessel_trunk)
% xlabel('Year 2000 Weekdays')
% ylabel('Fishing vessel number')
% 
% subplot(1,2,2)
% t=1:m;
% plot(t,opTvessel_trunk1, 'r', t, Tvessel_trunk1)
% legend('Predicted','Optimal')
% xlabel('Year 2000 Weekends')
% 
% 
% 
% stock_in(1:10)=stock_uni*ones(1,10);
%   
% %%stock without fishing
% for tt=11:dayt         
%                     days=log(tt);
%                     days_m=[days days^2 days^3 days^4];
%                     daysp=log(tt-1);
%                     days_mp=[daysp daysp^2 daysp^3 daysp^4];
%      stock_in(tt)=stock_in(tt-1)*exp(stock_parm(1:4)*(days_m'-days_mp'));
% end
% 
% figure(8)
% ttt=1:n;
% subplot(1,3,1)
% plot(ttt,ophar_tr,'r', ttt, har_tr)
% xlabel('Year 2000 weekdays')
% ylabel('Logged harvest')
% 
% ttt=1:m;
% subplot(1,3,2)
% plot(ttt,ophar_tr1,'r', ttt, har_tr1)
% %legend('Predicted','Optimal')
% xlabel('Year 2000 weekends')
% ylabel('Logged harvest')
% subplot(1,3,3)
% ttt=1:dayt;
% plot(ttt,stock_fin(1:dayt),'b',ttt,kk(1:dayt),'r',ttt,stock_in(1:dayt),'g','linewidth',2)
% %legend('Predicted','Optimal','Stock w/o harvest')
% xlabel('Year 2000')
% ylabel('Stock index')
% axis([0 400 0 95])
