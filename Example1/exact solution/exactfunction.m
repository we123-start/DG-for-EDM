% function [cexact] = exactfunction(Dapp, x, t)
% p = -0.5*Dapp*t*(pi/0.2)^2+i*pi/0.2*(0.2-x+0.5*t);
% ialpha = pi*sqrt(0.5*Dapp*t)/0.2;
% ibelta = pi*sqrt(0.5*Dapp*t)/0.2;
% alpha = (-0.2+x-0.5*t)/(2*sqrt(0.5*Dapp*t));
% belta = (-0.4+x-0.5*t)/(2*sqrt(0.5*Dapp*t));
% %cexact = 0.5.*real(i.*p.*(erf(sym(alpha-i*ialpha))-erf(sym(belta-i*ibelta))));
% %cexact = 0.5.*real(i.*p.*(erf(alpha)+erfi(-ialpha)-erf(belta)-erfi(-ibelta)));
% cexact = 0.5.*real(i.*exp(p).*(erfz(alpha-i*ialpha)-erfz(belta-i*ibelta)));
% end

function [cexact] = exactfunction(Dapp, x, t,a)
u=1;
Deef=Dapp/(1+a);ueef=u/(1+a);
p = -Deef*t*(pi/0.2)^2+i*pi/0.2*(0.2-x+ueef*t);
ialpha = pi*sqrt(Deef*t)/0.2;
ibelta = pi*sqrt(Deef*t)/0.2;
alpha = (-0.2+x-ueef*t)/(2*sqrt(Deef*t));
belta = (-0.4+x-ueef*t)/(2*sqrt(Deef*t));
%cexact = 0.5.*real(i.*p.*(erf(sym(alpha-i*ialpha))-erf(sym(belta-i*ibelta))));
%cexact = 0.5.*real(i.*p.*(erf(alpha)+erfi(-ialpha)-erf(belta)-erfi(-ibelta)));
cexact = 0.5.*real(i.*exp(p).*(erfz(alpha-i*ialpha)-erfz(belta-i*ibelta)));
end