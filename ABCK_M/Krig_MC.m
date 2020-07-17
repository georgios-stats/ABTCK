function [Kriging_W, Var_W]=Krig_MC(Zt,Beta_X,sigma2,AA,site,site_e,S,tau2)
% Kriging with nuget effect.
% [site1; site2] real and computer code data.  
% site = the acctual sites where we have data.
% site_e is the grided locations where we do prediction.
% ***** Usualy l_1 and l_2 equal with a*****%
% nn is a cell with the number of points in the Element (region)


if (size(site_e,1)>0)
    pseudo_sigma=cell(S,1);
    
Kriging_W=[];Var_W=[];

Zt_prev{1}=[]; Zt_prevE{1}=[];

for k=1:S
% pseudo_sigma{k}=1; 
% end

% && size(site1,1)>0 && size(site2,1)>0

if k==1
    nn=size(site{k},1); 
    X1=[ones(nn,1)];%kron(X_in,X_s);
    XX{k}=X1;
    
    nn_e=size(site_e,1); 
    Xe=[ones(nn_e,1)];
    XE{k}=Xe; % mean for all the GP's   
 else
 % Regression coeff   
    nn=size(site{k},1);
    X1=[ones(nn,1)];
    ISM11=ismember(site{k-1},site{k});
    ISM=logical(prod(ISM11,2));
    %ISM=logical(kron(ISM_in(:,1),ISM_s(:,1)));
    X2=[Zt{k-1}(ISM(:,1))];
    XX{k}=[X1  X2];
    
    nn_e=size(site_e,1); 
    Xe=[ones(nn_e,1)];
    XE{k}=[Xe Zt_prevE{k}]; % mean for all the GP's      
    
 end 
    
    Kriging_W{k}=[];
       % input covariance 
    C_In{k} =Gauss_Cov(1, AA{k}, site{k},site{k})+tau2{k}.*eye(nn);

% input covariance 
    c_In{k} =  Gauss_Cov(1, AA{k}, site{k}, site_e); 
% Kriging method    
    lambda_In{k} = (C_In{k}\c_In{k})'; 
Yt{k}=(Zt{k}-XX{k}*Beta_X{k});%reshape((Zt{k}-XX{k}*Beta_X{k}),nn_s,nn_in);
mulS=(lambda_In{k}*Yt{k});
%size(mulS1(:))
%mulS=mulS1(:);%kron(lambda_In{k},mulS1(:));
%size(mulS)
%toc

%    Kriging_W1{k}=XE{k}*Beta_X{k}+lambda{k}*(Zt{k}-XX{k}*Beta_X{k});
    Kriging_W{k}=XE{k}*Beta_X{k}+mulS;

    var_ST2=(diag(lambda_In{k}*c_In{k}));
    Var_W{k}= (ones(nn_e,1)-var_ST2).*sigma2{k};
     %Var_W{k}= (ones(n_in_e,1)-var_ST2).*sigma2{k};
%     Kriging_WLB{k}=Kriging_W{k}-1.96*sqrt(Var_W{k});
%     Kriging_WUB{k}=Kriging_W{k}+1.96*sqrt(Var_W{k});
%   Zt_prev{k+1}=[Zt_prev{k} Zt{k}];
    Zt_prevE{k+1}=[Zt_prevE{k} Kriging_W{k}];  
end
         
end
end

         
         