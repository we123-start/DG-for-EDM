function [P]=Lpoly(k,x)
%   the function of Legendre polynomial
l=length(x);
P=ones(k+1,l);
for i=1:k+1
if i==1
    P(1,:) = 1;
elseif i==2
    P(2,:) = x;
elseif i==3
    P(3,:) = (3*x.^2-1)/2;
elseif i==4
    P(4,:) = (5*x.^3-3*x)/2;
elseif i==5
    P(5,:) = (35*x.^4-30*x.^2+3)/8;
elseif i==6
    P(6,:) = (63*x.^5-70*x.^3+15*x)/8;
elseif i==7
    P(7,:) = (231*x.^6-315*x.^4+105*x.^2-5)/16;
else
    fprintf('There appears an error in Lpoly.m!\n')
end
end
return
end

