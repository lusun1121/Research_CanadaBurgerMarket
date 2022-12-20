clear;
tic;
load indiv.txt;%vesselid len vess vess_fixed
load daily0.txt;%year,y365,stock,prize_a, wspd,wvht, holiday,8week,9opening,totalhar
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
    

load ind_parm2.csv;%%1prize_a*vess  stock*vess	diesel_price*vess holiday*vess	
ind_parm=ind_parm2;

%yearch=-19.0923;%%CCP constant
% probit_par=[2.8513	0.0824	0.1502	-0.9473	-0.6162	0.0268	-0.0175	0.00353	-0.00006	...
%     -0.00213	-0.00005	-0.00075	1.15E-06	0.6659	0.000837	-0.093	0.0102	...
%     -0.00033	0.1096	-0.1245	0.00965	0.6752	-0.0168	-0.0201];
yearch=-2.987941960;
probit_par=[0.117414027	0.063966915	-0.007352238	-0.830043794	-0.033051607	0.001088292	...
    0.00685054	-0.000565685	0.000002878	-0.000262922	-0.000010657	-0.000005409	...
    0.000000009	0.095935324	0.000443973	-0.094368799	0.00044001	-0.000005967	...
    -0.183506342	-0.021570058	0.001685203	0.222659708	0.017507882	...
    0.133714599	0.002475202];


%prize_a	stock	diesel_price	holiday	prize_a2 prize_a3 prize_a*stock	stock2 stock3 prize_a*diesel_price	
%stock*diesel_price	diesel_pr*diesel_pri	diesel*diesel*diesel	prize_a*vess_fixed	
%stock*len	WSPD	WSPD*WSPD	WSPD*WSPD*WSPD	WVHT	WVHT*WVHT	WVHT*WVHT*WVHT	opening1	
%	vess_fixed*WSPD	vess_fixed*WVHT zeta


ffix=exp(indiv(:,4));%vess_fixed


 
for m=1:dayt;
    state_daily0(:,1)=ones(vvn,1)*daily(m,3);%stock, har
    state_daily(m,1:3)=daily(m,4:6);%prize_a, wspd,wvht
    state_daily(m,4)=diesel_price(m,1);%%%diesel price

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
 
 
%vess prize_a	stock	diesel_price	holiday	prize_a2 prize_a3 prize_a*stock	stock2 stock3 prize_a*diesel_price	
%stock*diesel_price	diesel_pr*diesel_pri	diesel*diesel*diesel	prize_a*vess_fixed	
%stock*len	WSPD	WSPD*WSPD	WSPD*WSPD*WSPD	WVHT	WVHT*WVHT	WVHT*WVHT*WVHT	opening1	
%	vess_fixed*WSPD	vess_fixed*WVHT
    prob=exp(probitt)./(1+exp(probitt));
    share1=sum(prob)-prob;
    share_both=(sum(prob)^2-sum(prob.^2))/2-prob.*(share1);
    share_ex=1+gamma+gamma^2/2+(gamma+3*gamma^2/2)*share1+gamma^2/2*share_both;
    har=exp(catch_parm(3)*state_daily(m,2)+catch_parm(4)*state_daily(m,3)+catch_parm(2)*daily(m,9))*share_ex.*ffix.*state_daily0(:,1);
    har2=har.^2;
    
%    effort_1=sum(indiv(:,5).*prob)+(1-prob).*indiv(:,5);%%%total effort
    totalharvest_1=sum(har.*prob)+(1-prob).*har;
    acc1=sum(daily(1:m-1,10))+totalharvest_1;
    %%everybody has a different total harvest for all vessels if they deveate from the equilibrium
    %effort_0=sum(indiv(:,5).*prob)-prob.*indiv(:,5);
    %acc0=sum(daily(1:m-1,10))+effort_0;
    totalharvest_0=sum(har.*prob)-prob.*har;
    acc0=sum(daily(1:m-1,10))+totalharvest_0;

for rep=1:timess;
    value0=zeros(vvn*dayt,24);
    value1=zeros(vvn*dayt,24);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%if you choose to fish given everybody else adopts optimal strategy%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%get accumulated sum;
    value1(vvn*(m-1)+1:vvn*m,1)=state_daily(m,1)*har; %%revenue
    value1(vvn*(m-1)+1:vvn*m,2:6)=ones(vvn,1)*[state_daily(m,2:4),daily(m,7),state_daily(m,4)^2];
    %%Wind Wave Diesel weekend Diesel2  
    value1(vvn*(m-1)+1:vvn*m,7:15)= [state_daily0(:,1),state_daily0(:,1).^2,indiv(:,2),indiv(:,2).^2, ...
        indiv(:,2)*state_daily(m,2:4), indiv(:,4),indiv(:,4).^2];
    %%stock stock2 len len2 len_wspd len_wvht len_disel vess_fixed vess_fixed2
    value1(vvn*(m-1)+1:vvn*m,16:23)=[zeros(vvn,1),har2,har,har.*indiv(:,2),har2.*indiv(:,2),ones(vvn,1), ones(vvn,1)*daily(m,9),zeros(vvn,1)];
    %% 0 har2,har,har_len,har2_len,year,opening
    value1(vvn*(m-1)+1:vvn*m,24)=ones(vvn,1)*exp(zetan(m,1));

    accum1=acc1;
    accum0=acc0;

    for tt=m+1:dayt
    days=log(tt);
    days_m=[days days^2 days^3 days^4];
    stock_1=exp(stock_year(year)+stock_parm(1:4)*days_m'+stock_parm(5)*accum1);
    if tt<10
        stock_1=10*ones(vvn,1);
    end

    state_daily(tt,1)= shrimp_price(tt-1,2)+sqrt(v_shrimp)*randn;   %%%shrimp price
    r = mvnrnd(mu,sigma,2);
    state_daily(tt,2)=exp(1.088944+0.406951*log(state_daily(tt-1,2))+r(1)); %%%wind speed
    state_daily(tt,3)=exp(0.039297+0.693213*log(state_daily(tt-1,3))+r(2)); %%%wave height
    state_daily(tt,4)=diesel_price(tt-1,2)+sqrt(v_diesel)*randn;       %%%diesel price
    zetan(tt,1)=zetan(tt-1,1)*(0.72201)+mvnrnd(0,0.19225,1); %%%%zeta
    
     probitt_new_1=yearch+ind_parm(:,1)+probit_par(1)*state_daily(tt,1)+probit_par(2)*stock_1+probit_par(3)*state_daily(tt,4) ...
     +probit_par(4)*daily(tt,7)+probit_par(5)*state_daily(tt,1)^2+probit_par(6)*state_daily(tt,1)^3  ...
     +probit_par(7)*state_daily(tt,1)*stock_1+probit_par(8)*(stock_1).^2+probit_par(9)*(stock_1).^3 ...
     +probit_par(10)*state_daily(tt,1)*state_daily(tt,4)+probit_par(11)*state_daily(tt,4)*stock_1+probit_par(12)*state_daily(tt,4)^2 ...
     +probit_par(13)*state_daily(tt,4)^3+probit_par(14)*state_daily(tt,1)*indiv(:,4)+probit_par(15)*stock_1.*indiv(:,2) ...
     +probit_par(16)*state_daily(tt,2)+probit_par(17)*state_daily(tt,2)^2+probit_par(18)*state_daily(tt,2)^3 ...
     +probit_par(19)*state_daily(tt,3)+probit_par(20)*state_daily(tt,3)^2+probit_par(21)*state_daily(tt,3)^3 ...
     +probit_par(22)*daily(tt,9)+probit_par(23)*state_daily(tt,2)*indiv(:,4)+probit_par(24)*state_daily(tt,3)*indiv(:,4)+probit_par(25)*exp(zetan(tt,1));
 

    %raner=mvnrnd(0,1,vvn);
    %fishing_new_1=(probitt_new_1+raner(:,1))>0;
    prob_1=exp(probitt_new_1)./(1+exp(probitt_new_1));
    drawnum = rand(vvn,1);
    fishing_new_1=(drawnum <= prob_1);

    share_new_1=sum(prob_1)-prob_1;
    share_new_both_1=(sum(prob_1)^2-sum(prob_1.^2))/2-prob_1.*(share_new_1);
    share_new_ex_1=1+gamma+gamma^2/2+(gamma+3*gamma^2/2)*share_new_1+gamma^2/2*share_new_both_1;
    har_new_1=exp(catch_parm(3)*state_daily(tt,2)+catch_parm(4)*state_daily(tt,3)+catch_parm(2)*daily(tt,9))*share_new_ex_1.*ffix.*stock_1;
  
    har2_new_1=har_new_1.^2;
    har_real_1=fishing_new_1.*har_new_1;
    %effort_1=sum(indiv(:,5).*fishing_new_1);%%%new  harvest
    %accum1=accum1+effort_1;
    totalharvest_1=sum(har_real_1);%%%new  harvest
    accum1=accum1+totalharvest_1;

  err_r_1= fishing_new_1.*(-log(prob_1));
  
    %%%get accumulated sum;
    value1(vvn*(m-1)+1:vvn*m,1)=value1(vvn*(m-1)+1:vvn*m,1)+dis^(tt-m+1)*state_daily(tt,1)*har_real_1;%%revenue
    value1(vvn*(m-1)+1:vvn*m,2:6)=value1(vvn*(m-1)+1:vvn*m,2:6)+dis^(tt-m+1)*fishing_new_1*[state_daily(tt,2:4),daily(tt,7),state_daily(tt,4)^2];
    %%Wind Wave Diesel weekend Diesel2
    value1(vvn*(m-1)+1:vvn*m,7:15)=value1(vvn*(m-1)+1:vvn*m,7:15) ...
        +dis^(tt-m+1)*(fishing_new_1*ones(1,9)).*[stock_1, stock_1.^2, indiv(:,2), indiv(:,2).^2, indiv(:,2)*state_daily(tt,2:4), indiv(:,4), indiv(:,4).^2];
    %%stock stock2 len len2 len_wspd len_wvht len_disel vess_fixed vess_fixed2
    value1(vvn*(m-1)+1:vvn*m,16:23)=value1(vvn*(m-1)+1:vvn*m,16:23) ...
        +dis^(tt-m+1)*(fishing_new_1*ones(1,8)).*[zeros(vvn,1),har2_new_1,har_new_1,har_new_1.*indiv(:,2), har2_new_1.*indiv(:,2), ...
        ones(vvn,1),ones(vvn,1)*daily(tt,9),zeros(vvn,1)];
        %%error har2,har,har_len,har2_len,year,opening
    value1(vvn*(m-1)+1:vvn*m,24)=value1(vvn*(m-1)+1:vvn*m,24)+dis^(tt-m+1)*fishing_new_1*exp(zetan(tt,1));

    value1(vvn*(m-1)+1:vvn*m,16)=value1(vvn*(m-1)+1:vvn*m,16)+dis^(tt-m+1)*err_r_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%if you choose not to fish given everybody else adopts optimal strategy%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%get accumulated sum;

    %%everybody has a different total harvest 
    %%for all vessels if they deveate from the equilibrium

    stock_0=exp(stock_year(year)+stock_parm(1:4)*days_m'+stock_parm(5)*accum0);
    if tt<10
        stock_0=10*ones(vvn,1);
    end
    
     probitt_new_0=yearch+ind_parm(:,1)+probit_par(1)*state_daily(tt,1)+probit_par(2)*stock_0+probit_par(3)*state_daily(tt,4) ...
     +probit_par(4)*daily(tt,7)+probit_par(5)*state_daily(tt,1)^2+probit_par(6)*state_daily(tt,1)^3  ...
     +probit_par(7)*state_daily(tt,1)*stock_0+probit_par(8)*(stock_0).^2+probit_par(9)*(stock_0).^3 ...
     +probit_par(10)*state_daily(tt,1)*state_daily(tt,4)+probit_par(11)*state_daily(tt,4)*stock_0+probit_par(12)*state_daily(tt,4)^2 ...
     +probit_par(13)*state_daily(tt,4)^3+probit_par(14)*state_daily(tt,1)*indiv(:,4)+probit_par(15)*stock_0.*indiv(:,2) ...
     +probit_par(16)*state_daily(tt,2)+probit_par(17)*state_daily(tt,2)^2+probit_par(18)*state_daily(tt,2)^3 ...
     +probit_par(19)*state_daily(tt,3)+probit_par(20)*state_daily(tt,3)^2+probit_par(21)*state_daily(tt,3)^3 ...
     +probit_par(22)*daily(tt,9)+probit_par(23)*state_daily(tt,2)*indiv(:,4)+probit_par(24)*state_daily(tt,3)*indiv(:,4)+probit_par(25)*exp(zetan(tt,1));
 

%%1prize_a*vess    stock*vess	diesel*vess holiday*vess	
%%1prize_a*prize_a stock*stock 3prize_a*stock prize_a*diesel_price 5prize_a*vess_fixed
%%stock*len 7WSPD 8WVHT 9opening1 10diesel_pr*diesel_pri 11len 12len*len 13vess_fixed
%%14vess_fixed*WSPD 15vess_fixed*WVHT 16vess_fixed*vess_fixed

 
    %raner=mvnrnd(0,1,vvn);
    %fishing_new_0=(probitt_new_0+raner(:,1))>0;
    %prob_0=normcdf(probitt_new_0);
    
    prob_0=exp(probitt_new_0)./(1+exp(probitt_new_0));
    drawnum = rand(vvn,1);
    fishing_new_0=(drawnum <= prob_0);

    share_new_0=sum(prob_0)-prob_0;
    share_new_both_0=(sum(prob_0)^2-sum(prob_0.^2))/2-prob_0.*(share_new_0);
    share_new_ex_0=1+gamma+gamma^2/2+(gamma+3*gamma^2/2)*share_new_0+gamma^2/2*share_new_both_0;
    har_new_0=exp(catch_parm(3)*state_daily(tt,2)+catch_parm(4)*state_daily(tt,3)+catch_parm(2)*daily(tt,9))*share_new_ex_0.*ffix.*stock_0;
    har2_new_0=har_new_0.^2;
    har_real_0=fishing_new_0.*har_new_0;
    %effort_0=sum(indiv(:,5).*fishing_new_0);%%%new  harvest
    %accum0=accum0+effort_0;
    totalharvest_0= sum(har_real_0);%%%new total harvest
    accum0=accum0+totalharvest_0;
    %%%get accumulated sum;

     err_r_0= fishing_new_0.*(-log(prob_0));
      
    value0(vvn*(m-1)+1:vvn*m,1)=value0(vvn*(m-1)+1:vvn*m,1)+dis^(tt-m+1)*state_daily(tt,1)*har_real_0;%%revenue
    value0(vvn*(m-1)+1:vvn*m,2:6)=value0(vvn*(m-1)+1:vvn*m,2:6)+dis^(tt-m+1)*fishing_new_0*[state_daily(tt,2:4),daily(tt,7),state_daily(tt,4)^2];
    %%Wind Wave Diesel weekend Diesel2
    value0(vvn*(m-1)+1:vvn*m,7:15)=value0(vvn*(m-1)+1:vvn*m,7:15) ...
        +dis^(tt-m+1)*(fishing_new_0*ones(1,9)).*[stock_0, stock_0.^2, indiv(:,2), indiv(:,2).^2, indiv(:,2)*state_daily(tt,2:4), ...
        indiv(:,4), indiv(:,4).^2];
    %%stock stock2 len len2 len_wspd len_wvht len_disel vess_fixed vess_fixed2
    value0(vvn*(m-1)+1:vvn*m,16:23)=value0(vvn*(m-1)+1:vvn*m,16:23) ...
        +dis^(tt-m+1)*(fishing_new_0*ones(1,8)).*[zeros(vvn,1),har2_new_0,har_new_0,har_new_0.*indiv(:,2), har2_new_0.*indiv(:,2), ...
        ones(vvn,1),ones(vvn,1)*daily(tt,9),zeros(vvn,1)];
    value0(vvn*(m-1)+1:vvn*m,24)=value0(vvn*(m-1)+1:vvn*m,24) ...
        +dis^(tt-m+1)*fishing_new_0*exp(zetan(tt,1));
    %%error har2,har,har_len,har2_len,year,opening
    value0(vvn*(m-1)+1:vvn*m,16)=value0(vvn*(m-1)+1:vvn*m,16)+dis^(tt-m+1)*err_r_0;
   end; 
final=final+(value1-value0)/timess;
end;
clear  state_daily;
end;
dlmwrite('final.csv', final);
toc;
    
        