function [y]=funn6_sec(x1,x2,x3,x4,x5,x6)
y=(x2+0.5).*exp(sin((0.9.*(x1+0.48)).^(10)))+x3.*x4+0*x5+0*x6;
%y=(x2+0.8).*exp(sin(((x3+0.4).*(x1+0.48)).^(10)))+x4+0*x5+0*x6;
%y=(x2+0.5).*exp(sin((x3.*(x1+0.48)).^(10)))+x4+0*x5+0*x6;
%y=exp(sin((0.9.*(x2+0.48)).^(10)))+x1.*x3+x4+0*x5+0*x6;
end
