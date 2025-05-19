function [D_P]=D_Lpoly(k,x)
%   the function of the derivative of Legendre polynomial
l=length(x);
D_P=ones(k+1,l);
for i=1:k+1
if i==1
    D_P(1,:) = 0;
elseif i==2
    D_P(2,:) = 1;
elseif i==3
    D_P(3,:) = 3*x;
elseif i==4
    D_P(4,:) = (15*x.^2-3)/2;
elseif i==5
    D_P(5,:) = (35*x.^3-15)/2;
elseif i==6
    D_P(6,:) = (315*x.^4-210*x.^2+15)/8;
elseif i==7
    D_P(7,:) = (693*x.^5-630*x.^3+210*x)/8;
else
    fprintf('There appears an error in D_Lpoly.m!\n')
end
end
return
end

