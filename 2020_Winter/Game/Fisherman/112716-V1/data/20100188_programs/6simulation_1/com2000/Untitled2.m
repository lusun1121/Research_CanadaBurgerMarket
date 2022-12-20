x= (1:0.25:10);    y= x.^2;
E3= [ones(1, length(x))*30; 0.7*y];  % 30 above and 0.3*y below
E3(1,7)= Inf;      E3(2,10)= Inf;       E3(2,25)= NaN;
errorfill(x, y, [0.1 0.02], 0.3, E3, 'b-+');
legend('E3', 'E2', 'E1', 'x^2','Inf','-Inf','NaN');