function [X,Y,x,y] = gen_grid(stats,mle,s,n)
% n is an estimate of the number of points in a vector

fac = (n-1)/5;
x_unit = round(2*s(1),1,'significant');
x = floor((mle(1)-4*s(1))/x_unit)*x_unit:...
    x_unit/fac:...
    ceil((mle(1)+4*s(1))/x_unit)*x_unit;
y_unit = round(2*s(2),1,'significant');
y = floor((mle(2)-4*s(2))/y_unit)*y_unit:...
    y_unit/fac:...
    ceil((mle(2)+4*s(2))/y_unit)*y_unit;
[X,Y] = meshgrid(x,y);

end

