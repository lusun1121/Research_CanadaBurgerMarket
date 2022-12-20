
%clear;
load daily0.txt;   %year,y365,stock, prize_a, wspd,wvht, holiday,9week,opening,totalhar
load zeta.txt; %%year y365 zeta
    stock_year=[0.577319625	0.066830688	0.689766452	0.381456417	0.080095907	0]; %%stock dynamics yearly dummy
    stock_parm=[9.854831206	-7.394586282	1.796771306	-0.137313124	-0.000000127];%%stock dynamics coefficients
    catch_parm=[-0.001270945	0.717147271	-0.00407091	0.022396933];%%Tvessel opening wspd wvht        
 
dayt=365;
scale1=1/10;

  year=1;
  zeta1=zeta(365*(year-1)+1:365*year,3);
  daily=daily0(365*(year-1)+1:365*year,:);
  tt=10;
    days=log(tt);
    days_m=[days days^2 days^3 days^4];
stock_i(1:tt)=exp(stock_year(year)+stock_parm(1:4)*(days_m'))*ones(1,tt);
stock_i1(1:tt)=exp(stock_year(year)+stock_parm(1:4)*(days_m'))*ones(1,tt);
stock_in(1:tt)=exp(stock_year(year)+stock_parm(1:4)*(days_m'))*ones(1,tt);
stock_i2(1:tt)=stock_i1(1:tt);   


%stock without fishing
   
   effort=daily(:,4);
  % effort(1:150,1)=0;
   
for tt=11:dayt
    days=log(tt);
    days_m=[days days^2 days^3 days^4];
     daysp=log(tt-1);
     days_mp=[daysp daysp^2 daysp^3 daysp^4];
     
     stock_i(tt)=(round(stock_i(tt-1)/scale1))*scale1*exp(stock_parm(1:4)*(days_m'-days_mp')); 
      
     stock_i(tt)=stock_i(tt-1)*exp(stock_parm(1:4)*(days_m'-days_mp')); 
     stock_i1(tt)=((round(stock_i1(tt-1)/scale1)))*scale1*exp(stock_parm(1:4)*(days_m'-days_mp')+stock_parm(5)*daily(tt-1,10));
stock_i2(tt)=stock_i2(tt-1)*exp(stock_parm(1:4)*(days_m'-days_mp')+stock_parm(5)*daily(tt-1,10));
stock_i2(tt)=stock_i2(tt-1)*exp(stock_parm(1:4)*(days_m'-days_mp')+stock_parm(5)*daily(tt-1,10));
end
zeta2=zeta1(10);
zetaa=[];
for tt=11:dayt
    days=log(tt);
    days_m=[days days^2 days^3 days^4];
    %zeta2=-0.1039+(-0.8008)*zeta2+mvnrnd(0,0.43758,1);
        zeta2=-0.622+(-0.6991)*zeta2+mvnrnd(0,0.25686,1);
    zetaa=[zetaa zeta2];
    stock_in(tt)=exp(stock_year(year)+stock_parm(1:4)*days_m'+stock_parm(5)*sum(daily(1:tt-1,10))+zeta2);  
end
figure(1)
plot(1:dayt,stock_i1,'--',1:dayt,stock_in,1:dayt,stock_i,'linewidth',2)
legend('Stock','Stock accum', 'no harvest')
 xlabel('Year 2000')
 ylabel('Stock index')
 
 figure(2)
plot(1:dayt,stock_i1,'--',1:dayt,stock_i,1:dayt,stock_in,'linewidth',2)