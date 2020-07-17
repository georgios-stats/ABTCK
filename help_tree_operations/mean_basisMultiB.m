

function [Gt]=mean_basisMultiB(site,Zt,S,funname2)
if (S==1)
  nn{1}=size(site{1},1);
    Gt{1}=[ones(nn{1},1)];
else    
  nn{1}=size(site{1},1);
    Gt{1}=[ones(nn{1},1)];    
for t=2:S
    %Constant mean only 
switch funname2
case 'constant'
    nn{t}=size(site{t},1);
     X1 = [ones(nn{t},1)]; 
    ISM=ismember(site{t-1},site{t});
     X2=[Zt{t-1}(ISM(:,1))];
     Gt{t}=[X1  X2];
case 'linear'
     nn{t}=size(site{t},1);
     X1 = [ones(nn{t},1) site{t}] ;
    ISM=ismember(site{t-1},site{t});
     X2=[Zt{t-1}(ISM(:,1)) site{t-1}(ISM(:,1)).*Zt{t-1}(ISM(:,1))];
     Gt{t}=[X1  X2];
case 'quadratic'
      nn{t}=size(site{t},1);
     X1 = [ones(nn{t},1)] ;
%     XXbspline_basismatrix(3,[(min(site{t})+(max(site{t})-min(x))/7):(max(x)-min(x))/7:(max(site{t})-(max(site{t})-min(site{t}))/7)],site{t});
     %X1 = [ones(nn{t},1) site{t}  site{t}.^2  site{t}.^3] ;
    ISM=ismember(site{t-1},site{t});
     X2=[Zt{t-1}(ISM(:,1)) site{t-1}(ISM(:,1)).*Zt{t-1}(ISM(:,1)) (site{t-1}(ISM(:,1)).^2).*Zt{t-1}(ISM(:,1))];
     Gt{t}=[X1  X2];
     
%      ISM=ismember(site{t-1},site{t});
%      X1 = [ones(nn{t},1)];
%      XX=bspline_basismatrix(3,[(min(site{t})+(max(site{t})-min(site{t}))/7):(max(site{t})-min(site{t}))/7:(max(site{t})-(max(site{t})-min(site{t}))/7)],site{t});
%      X2=[XX(:,1).*Zt{t-1}(ISM(:,1)) XX(:,2).*Zt{t-1}(ISM(:,1)) (XX(:,3)).*Zt{t-1}(ISM(:,1))];
%      Gt{t}=[X1  X2];
%      
%      
%      
%      ISM=ismember(site{t-1},site{t});
%      X1 = [ones(nn{t},1)];
%      XX=bspline_basismatrix(3,[(min(site{t})+(max(site{t})-min(site{t}))/8):(max(site{t})-min(site{t}))/8:(max(site{t})-(max(site{t})-min(site{t}))/8)],site{t});
%      X2=[XX(:,1).*Zt{t-1}(ISM(:,1)) XX(:,2).*Zt{t-1}(ISM(:,1)) (XX(:,3)).*Zt{t-1}(ISM(:,1)) (XX(:,4)).*Zt{t-1}(ISM(:,1))];
%      Gt{t}=[X1  X2];
     
%      
%      ISM=ismember(site{t-1},site{t});
%      X1 = [ones(nn{t},1)];
%      XX=bspline_basismatrix(2,[(min(site{t})+(max(site{t})-min(site{t}))/7):(max(site{t})-min(site{t}))/7:(max(site{t})-(max(site{t})-min(site{t}))/7)],site{t});
%      X2=[XX(:,1).*Zt{t-1}(ISM(:,1)) XX(:,2).*Zt{t-1}(ISM(:,1)) (XX(:,3)).*Zt{t-1}(ISM(:,1)) (XX(:,4)).*Zt{t-1}(ISM(:,1))];
%      Gt{t}=[X1  X2];
%      
     
case 'polythree'
     nn{t}=size(site{t},1);
     X1 = [ones(nn{t},1)] ;
    ISM=ismember(site{t-1},site{t});
     X2=[Zt{t-1}(ISM(:,1)) site{t-1}(ISM(:,1)).*Zt{t-1}(ISM(:,1))  (site{t-1}(ISM(:,1)).^2).*Zt{t-1}(ISM(:,1))  (site{t-1}(ISM(:,1)).^3).*Zt{t-1}(ISM(:,1))];
     Gt{t}=[X1  X2];
 case 'polyfour'
     nn{t}=size(site{t},1);
     X1 = [ones(nn{t},1)] ;
    ISM=ismember(site{t-1},site{t});
     X2=[Zt{t-1}(ISM(:,1)) site{t-1}(ISM(:,1)).*Zt{t-1}(ISM(:,1))  (site{t-1}(ISM(:,1)).^2).*Zt{t-1}(ISM(:,1))  (site{t-1}(ISM(:,1)).^3).*Zt{t-1}(ISM(:,1)) (site{t-1}(ISM(:,1)).^4).*Zt{t-1}(ISM(:,1))];
     Gt{t}=[X1  X2];
case 'polypente'
     nn{t}=size(site{t},1);
     X1 = [ones(nn{t},1)] ;
    ISM=ismember(site{t-1},site{t});
     X2=[Zt{t-1}(ISM(:,1)) site{t-1}(ISM(:,1)).*Zt{t-1}(ISM(:,1))  (site{t-1}(ISM(:,1)).^2).*Zt{t-1}(ISM(:,1))  (site{t-1}(ISM(:,1)).^3).*Zt{t-1}(ISM(:,1)) (site{t-1}(ISM(:,1)).^4).*Zt{t-1}(ISM(:,1)) (site{t-1}(ISM(:,1)).^5).*Zt{t-1}(ISM(:,1))];
     Gt{t}=[X1  X2];
end
end 
end
end
