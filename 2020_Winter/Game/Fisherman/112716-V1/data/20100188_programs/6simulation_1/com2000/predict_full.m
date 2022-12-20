clear;
tic;
load indiv.txt;%vesselid len vess vess_fixed
load daily0.txt;%year,y365,stock,prize_a, wspd,wvht, holiday,8week,9opening,totalhar
load shrimp_price0.txt; %price predicted
load diesel_price0.txt; %price predicted
load zeta.txt;

    year=1;
    timess=1;
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
    final=zeros(vvn*dayt,24);
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
    

load ind_parm1.csv;%%1prize_a*vess  stock*vess	diesel_price*vess holiday*vess	
ind_parm=ind_parm1;

%yearch=-19.0923;%%CCP constant
% probit_par=[2.8513	0.0824	0.1502	-0.9473	-0.6162	0.0268	-0.0175	0.00353	-0.00006	...
%     -0.00213	-0.00005	-0.00075	1.15E-06	0.6659	0.000837	-0.093	0.0102	...
%     -0.00033	0.1096	-0.1245	0.00965	0.6752	-0.0168	-0.0201];
yearch=-19.7419;%%CCP constant
probit_par=[3.0774	0.111	0.1534	-0.9464	-0.6834	0.0306	-0.0151	0.00223	-0.00005	...
    -0.00277	-9.24E-06	-0.00077	1.20E-06	0.6656	0.000837	-0.0958	0.0108...
    -0.00037	0.0839	-0.1148	0.00913	0.6886	-0.0169	-0.0179	0.1655];


%prize_a	stock	diesel_price	holiday	prize_a2 prize_a3 prize_a*stock	stock2 stock3 prize_a*diesel_price	
%stock*diesel_price	diesel_pr*diesel_pri	diesel*diesel*diesel	prize_a*vess_fixed	
%stock*len	WSPD	WSPD*WSPD	WSPD*WSPD*WSPD	WVHT	WVHT*WVHT	WVHT*WVHT*WVHT	opening1	
%	vess_fixed*WSPD	vess_fixed*WVHT zeta


ffix=exp(indiv(:,4));%vess_fixed


    m=1;
    state_daily0(:,1)=ones(vvn,1)*daily(m,3);%stock, har
    state_daily(1,1)=shrimp_price(m,1);%prize_a
    state_daily(1,2:3)=daily(m,5:6);% wspd,wvht
    state_daily(1,4)=diesel_price(m,1);%%%diesel price

%1stock, 2har,
%1prize_a, 2wspd,3wvht,4diesel
%vesselid len vess vess_fixed
    probitt=yearch+ind_parm(:,1)+probit_par(1)*state_daily(m,1)+probit_par(2)*state_daily0(:,1)+probit_par(3)*state_daily(m,4) ...
     +probit_par(4)*daily(m,7)+probit_par(5)*state_daily(m,1)^2+probit_par(6)*state_daily(m,1)^3  ...
     +probit_par(7)*state_daily(m,1)*state_daily0(:,1)+probit_par(8)*(state_daily0(:,1)).^2+probit_par(9)*(state_daily0(:,1)).^3 ...
     +probit_par(10)*state_daily(m,1)*state_daily(m,4)+probit_par(11)*state_daily(m,4)*state_daily0(:,1)+probit_par(12)*state_daily(m,4)^2 ...
     +probit_par(13)*state_daily(m,4)^3+probit_par(14)*state_daily(m,1)*indiv(:,4)+probit_par(15)*state_daily0(:,1).*indiv(:,2) ...
     +probit_par(16)*state_daily(m,2)+probit_par(17)*state_daily(m,2)^2+probit_par(18)*state_daily(m,2)^3 ...
     +probit_par(19)*state_daily(m,3)+probit_par(20)*state_daily(m,3)^2+probit_par(21)*state_daily(m,3)^3 ...
     +probit_par(22)*daily(m,9)+probit_par(23)*state_daily(m,2)*indiv(:,4)+probit_par(24)*state_daily(m,3)*indiv(:,4)+probit_par(25)*exp(zetan(m,1));
 
%%vess 1prize_a	2stock	3diesel_price	4holiday	5prize_a*prize_a	6prize_*prize_*prize_	
%%7prize_a*stock	8stock*stock	9stock*stock*stock	
%%10prize_a*diesel_price	11stock*diesel_price	12diesel_pr*diesel_pri	
%%13diesel*diesel*diesel	14prize_a*vess_fixed	15stock*len	
%%16WSPD	17WSPD*WSPD	18WSPD*WSPD*WSPD	
%%19WVHT	20WVHT*WVHT	21WVHT*WVHT*WVHT	
%%22opening1	23vess_fixed*WSPD	24vess_fixed*WVHT

   prob=exp(probitt)./(1+exp(probitt));
        share1=sum(prob);
    share_both=(sum(prob)^2-sum(prob.^2))/2-prob.*(share1);
    share_ex=1+gamma+gamma^2/2+(gamma+3*gamma^2/2)*share1+gamma^2/2*share_both;
    har=exp(catch_parm(3)*state_daily(2)+catch_parm(4)*state_daily(3)+catch_parm(2)*daily(m,9))*share_ex.*ffix.*state_daily0(:,1);
    har2=har.^2;
 %%Tvessel opening wspd wvht    
    
    drawnum = rand(vvn,1);
    fishing_new=(drawnum <= prob);
    
    har_real=fishing_new.*har;
    
    Total_vessel(m)=share1;    
    harvest(m)=sum(har_real);

    
%     effort_1=sum(indiv(:,5).*prob)+(1-prob).*indiv(:,5);%%%total effort
%     acc1=sum(daily(1:m-1,4))+effort_1;


for rep=1:timess;
    accum1=harvest(m);
    for tt=m:dayt-1
    days=log(tt);
    days_m=[days days^2 days^3 days^4];
    stock_1=exp(stock_year(year)+stock_parm(1:4)*days_m'+stock_parm(5)*accum1);
    if tt<10
        stock_1=10*ones(vvn,1);
    end

%     state_daily(tt-m+2,1)= shrimp_price(tt-m+1,2)+sqrt(v_shrimp)*randn;   %%%shrimp price
%     r = mvnrnd(mu,sigma,2);
%     state_daily(tt-m+2,2)=exp(1.088944+0.406951*log(state_daily(tt-m+1,2))+r(1)); %%%wind speed
%     state_daily(tt-m+2,3)=exp(0.039297+0.693213*log(state_daily(tt-m+1,3))+r(2)); %%%wave height
%     state_daily(tt-m+2,4)=diesel_price(tt-m+1,2)+sqrt(v_diesel)*randn;       %%%diesel price

    state_daily(tt-m+2,1)= shrimp_price(tt-m+1,2);%prize_a
    state_daily(tt-m+2,2:3)= daily(tt-m+1,5:6);% wspd,wvht
    state_daily(tt-m+2,4)=diesel_price(tt-m+1,1);%%%diesel price
    zetan(tt-m+2,1)=zetan(tt-m+1,1)*(0.72201)+mvnrnd(0,0.19225,1); %%%%zeta
  
%     probitt_new_1=yearch+ind_parm(:,1)*state_daily(tt-m+2,1)+ind_parm(:,2).*stock_1+ind_parm(:,3)*state_daily(tt-m+2,4)+ind_parm(:,4)*daily(tt,8) ...
%      +probit_par(1)*state_daily(tt-m+2,1)^2+probit_par(2)*stock_1.^2+probit_par(3)*stock_1*state_daily(tt-m+2,1) ...
%      +probit_par(4)*state_daily(tt-m+2,4)*state_daily(tt-m+2,1)+probit_par(5).*indiv(:,4)*state_daily(tt-m+2,1) ...
%      +probit_par(6)*stock_1.*indiv(:,2)+probit_par(7)*state_daily(tt-m+2,2)+probit_par(8)*state_daily(tt-m+2,3) ...
%      +probit_par(9)*daily(tt,10)+probit_par(10)*state_daily(tt-m+2,4)^2+probit_par(14)*state_daily(tt-m+2,2)*indiv(:,4) ...
%      +probit_par(15)*state_daily(tt-m+2,3)*indiv(:,4)+ffix1;
 
     probitt_new_1=yearch+ind_parm(:,1)+probit_par(1)*state_daily(tt-m+2,1)+probit_par(2)*stock_1+probit_par(3)*state_daily(tt-m+2,4) ...
     +probit_par(4)*daily(tt-m+2,7)+probit_par(5)*state_daily(tt-m+2,1)^2+probit_par(6)*state_daily(tt-m+2,1)^3  ...
     +probit_par(7)*state_daily(tt-m+2,1)*stock_1+probit_par(8)*(stock_1).^2+probit_par(9)*(stock_1).^3 ...
     +probit_par(10)*state_daily(tt-m+2,1)*state_daily(tt-m+2,4)+probit_par(11)*state_daily(tt-m+2,4)*stock_1+probit_par(12)*state_daily(tt-m+2,4)^2 ...
     +probit_par(13)*state_daily(tt-m+2,4)^3+probit_par(14)*state_daily(tt-m+2,1)*indiv(:,4)+probit_par(15)*stock_1.*indiv(:,2) ...
     +probit_par(16)*state_daily(tt-m+2,2)+probit_par(17)*state_daily(tt-m+2,2)^2+probit_par(18)*state_daily(tt-m+2,2)^3 ...
     +probit_par(19)*state_daily(tt-m+2,3)+probit_par(20)*state_daily(tt-m+2,3)^2+probit_par(21)*state_daily(tt-m+2,3)^3 ...
     +probit_par(22)*daily(tt-m+2,9)+probit_par(23)*state_daily(tt-m+2,2)*indiv(:,4)+probit_par(24)*state_daily(tt-m+2,3)*indiv(:,4) ...
     +probit_par(25)*exp(zetan(tt-m+2,1));
 
%%prize_a	stock	diesel_price	holiday	prize_a*prize_a	prize_*prize_*prize_	
%%prize_a*stock	stock*stock	stock*stock*stock	prize_a*diesel_price	
%%stock*diesel_price	diesel_pr*diesel_pri	diesel*diesel*diesel	%%
%%prize_a*vess_fixed	stock*len	WSPD	WSPD*WSPD	WSPD*WSPD*WSPD	WVHT	
%%WVHT*WVHT	WVHT*WVHT*WVHT	opening1	vess_fixed*WSPD	vess_fixed*WVHT


    prob_1=exp(probitt_new_1)./(1+exp(probitt_new_1));
    drawnum = rand(vvn,1);
    fishing_new_1=(drawnum <= prob_1);

    share_new_1=sum(prob_1);
    share_new_both_1=(sum(prob_1)^2-sum(prob_1.^2))/2-prob_1.*(share_new_1);
    share_new_ex_1=1+gamma+gamma^2/2+(gamma+3*gamma^2/2)*share_new_1+gamma^2/2*share_new_both_1;
    har_new_1=exp(catch_parm(3)*state_daily(tt-m+2,2)+catch_parm(4)*state_daily(tt-m+2,3)+catch_parm(2)*daily(tt-m+2,9))*share_new_ex_1.*ffix.*stock_1;
  
    har2_new_1=har_new_1.^2;
    har_real_1=fishing_new_1.*har_new_1;
   % effort_1=sum(indiv(:,5).*fishing_new_1);%%%new  harvest
    %accum1=accum1+effort_1;
    
    Total_vessel(tt-m+2)=share_new_1;    
    harvest(tt-m+2)=sum(har_real_1);
    
     accum1=accum1+harvest(tt-m+2);

    end; 
end;
clear state_daily0 state_daily;

toc;  
t=1:dayt;

reg=[Total_vessel',daily(t,end),harvest',daily(t,end-1), daily(t,8)];
dlmwrite('reg.txt', reg);

figure (1)
plot(t,Total_vessel, t,daily(t,end))
legend('Predicted','Observed')
xlabel('Year 2000')
ylabel('Fishing vessel number')

 figure (2)
 plot(t,harvest, t,daily(t,end-1))     

hh0=[];
hh=[daily(:,end) Total_vessel' daily(:,end-1)];
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



