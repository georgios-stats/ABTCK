function [y]=funnLeGtariet(x,funname)
switch funname
    case 'level1'
      y=0.5*(((6*x-2).^2).*(sin(12*x-4)))+10*(x-0.5)-0.5;
    case 'level2'
      y=2*(0.5.*(((6.*x-2).^2).*(sin(12.*x-4)))+10.*(x-0.5)-0.5)-20*x+20+sin(10.*cos(5.*x));
end