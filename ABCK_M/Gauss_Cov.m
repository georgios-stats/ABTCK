function [Cmat]=Gauss_Cov(sigma_2, Lambda, site1,site2)
n = size(site1,1); m = size(site2,1);
Dim=size(site1,2);
SIGMA=Lambda'*Lambda;%SIGMA=Lambda;%
inv_S=diag(1./diag(SIGMA));%inv_S=inv(SIGMA);
DD=0;
for k=1:Dim
SITE_M1n{k}=kron(site1(:,k),ones(1,m));
SITE_M1m{k}=kron(site2(:,k),ones(1,n));
DSITE_1{k}=SITE_M1n{k}-SITE_M1m{k}';
D2{k} = inv_S(k,k)*DSITE_1{k}.^2;
DD=(D2{k}+DD);
end
Q=(DD);
m3=exp(-Q); 
if n==m
  result =m3;%sigma_2.*(m3+eye(n)*tau2); %+eye(n)*0.0005
else
  result =m3;%sigma_2.*(m3); %+eye(n)*0.05
    
end
  Cmat=result;  
end