clear all; close all;
% load covariance function
addpath ./cov_functions
% add functions
addpath ./function_fold
addpath ./function_fold/conditionalLHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load functions to generate random variables 
addpath ./Random_var; 
% load tree operations 
addpath ./Multi_GP_operations;
addpath ./Multi_tree_NON_nested;

addpath ./Multi_tree_nested;
addpath ./help_tree_operations; 
addpath ./Multi_prediction
addpath ./ABCK_M
% Load separable model MCMC_operations
%addpath ./GP_operations; 
%%
% warning('off', 'MSGID')
warning('off')
 funname0='first';
% % Generate data from a two dimentional input and two dimentional output.
 x1=[0:(1/100):1];

% % First example
 y1=funnPerdikaris(x1,'level1');%+0.5*(Xxx(:,1)/3).*sin(10*cos(2*(Xxx(:,1)/12)*pi))+1+norm_r(0,0.001,size(Xxx,1));%+0.1*sin(10*cos(2*(Xxx(:,1)/12)*pi));%+0.1*cos(2*(Xxx(:,2))*pi)+0.5+0.1*Xxx(:,2);
  %figure; scaterplot(x1,x2,y1); set(gca,'YDir','normal') 
  
  
  %%

   y2=funnPerdikaris(x1,'level2'); %+0.5*(Xxx(:,1)/3).*sin(10*cos(2*(Xxx(:,1)/12)*pi))+1+norm_r(0,0.001,size(Xxx,1));%+0.1*sin(10*cos(2*(Xxx(:,1)/12)*pi));%+0.1*cos(2*(Xxx(:,2))*pi)+0.5+0.1*Xxx(:,2);
  %figure; scaterplot(x1,x2,y1); set(gca,'YDir','normal') 
  
  Dim_x=1;

    
  
 %%
for ii=1:10
%ii=1;
close all; 
lim_all=[0 1]; 
rng(11+ii); % rng(111); rng(112); rng(111)
n=60; m=20;
LH_S11=lhsdesign(n-m,1)';
%LH_S12=lhsdesign(m,1)';
LH_S12=[0:(1-0)/(m-1):1];

Dim_x=1; %Dim_y=2; M=Dim_y; 

%LH_S=lhsdesign(m,Dim_x)';%[0,0.4,0.6,1];%
SITE1=[LH_S11'; LH_S12'];


%%
mm=200; X1=[0:1/mm:1]'; 
site_e=[X1(:)];  n_e =(mm+1);
x_new=LH_S12';
% A second equation
data_y=funnPerdikaris(x_new,'level2');
data_eta=funnPerdikaris(SITE1,'level1');
SITE_new{1}{1}= SITE1; SITE_new{1}{2}=x_new;
ZT_new{1}{1}=data_eta; ZT_new{1}{2}=data_y; 
p=1; a_old=1; BB=50000;


hold on; plot(SITE1, data_eta,'b.');
plot(x_new, data_y,'r*');
% xlabel(['x'], 'interpreter','latex');
% ylabel(['$t_1$'], 'interpreter','latex');
% zlabel(['z'], 'interpreter','latex');
% hold off; set(fighdl6_R, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%        MCMC algorith    %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funname='gauss'; %funname2='anisotropic';

funname2='polythree';%funname2='polythree';%funname2='polyfour';%funname2='quadratic'; %%funname2 = 'constant';%  funname2 ='linear'; %  funname2='polyfour';%
Delta=0.4; 
Dim_x=1;
a_NODES=1;
 th_bound.min=0;
 th_bound.max=1;
  mu_Z1=mean(data_eta); 
 Mu_Z1= kron(mu_Z1,ones(m,1));
 mu_Z2=mean(data_y); 
 Mu_Z2= kron(mu_Z2,ones(n,1)); 
 AA_E{1}{1}=0.1*eye(p); AA_E{1}{2}=0.1*eye(p);
LIM_ALL{1}=[-0.0001 1.001]; 
 mu_Z1=mean(data_eta); 
 NODES{1}.number =1; %  
 NODES{1}.inter =0; % exterior 
 NODES{1}.parent_child=0; % Left child 
 NODES{1}.direction=1;
 NODES{1}.rul=LIM_ALL{1}(1,2);
 
aalpha=0.8; bbeta=2; %aalpha=0.05;  bbeta=20; %
mm_e=130; X11=[0:1/mm_e:1]';
site_e_eta=X11(:);
 PR_CURR{1}{1}=10^(-20); PR_CURR{1}{2}=10^(-20);        
prior{1}.a=2.5; prior{1}.b=2.5;
prior{2}.a=2.5; prior{2}.b=2.5;
tau1=0;%0.005;
%tau2=10^(-3);
tau2=10^(-6);
%Different choices here an impact
%pprop.a1=8; pprop.b1=4; pprop.a2=15; pprop.b2=2; pprop.probm =0.5;
pprop.a1=6; pprop.b1=4; pprop.a2=10; pprop.b2=2; pprop.probm =0.5;
%pprop.a1=1; pprop.b1=10; pprop.a2=10; pprop.b2=5; pprop.probm =0.5;% a1=2; b1=2; a2=20; b2=2; probm =0.5;
%pprop.a1=2; pprop.b1=0.5; pprop.a2=10; pprop.b2=3; pprop.probm=0.5;
%   pprop.a1=2;
 %  pprop.b1=1;%pprop.b1=0.5;%
 %  pprop.a2=20;
 %  pprop.b2=2;
 %  pprop.probm=0.5;
bound{1}.low=10^(-3); bound{1}.up=70;  bound{2}.low=0.1;bound{2}.up=60;
 a1=pprop.a1; b1=pprop.b1; a2=pprop.a2; b2=pprop.b2; probm =pprop.probm;
% 
for kk=1:1000
yxy(kk)=rand_mix_gamma(a1,b1,a2,b2,probm);
end
figure; hist(yxy)




S=2; low_l=6; a_new=1;
TAU2{1}{1}=10^(-4);
TAU2{1}{2}=10^(-7);
nn{1}{1}=size(ZT_new{1}{1},1);
nn{1}{2}=size(ZT_new{1}{2},1);


%%
[ZtBN_Q,siteBN_Q,prCurrBN_Q,AA_BQ,a_newBNQ,lim_allBNQ, KrigSBNQ,Var_krigSBNQ,Krig_rhoBNQ, Krig_W_etaMATBNQ, Krig_site_etaBNQ]=...
    BTMultiB(PR_CURR,SITE_new,ZT_new,nn,AA_E,LIM_ALL,a_new,a_NODES,NODES,pprop, aalpha, bbeta,S,low_l,site_e,n_e,BB,tau1,TAU2,Delta,bound,funname,funname2);

Pred_MeanGPNQ=mean(Krig_W_etaMATBNQ(1:BB,:))';
Real_val =  funnPerdikaris(Krig_site_etaBNQ,'level2');
MSPE_GPNQ=mean(((Pred_MeanGPNQ-Real_val).^2));
bound{1}.low=10^(-5); bound{1}.up=40; bound{2}.low=0.000001;%
bound{2}.up=20;

[ZtBN,siteBN,prCurrBN,AA_B,a_newBN,lim_allBN, KrigSBN,Var_krigSBN,Krig_rhoBN, Krig_W_etaMATBN, Krig_site_etaBN]=...
    BMulti(PR_CURR,SITE_new,ZT_new,nn,AA_E,LIM_ALL,a_new,a_NODES,NODES,pprop, aalpha, bbeta,S,low_l,site_e,n_e,BB,tau1,tau2,Delta,bound,funname,funname2);

Pred_MeanGPN=mean(Krig_W_etaMATBN(1:BB,:))';
Real_val =  funnPerdikaris(Krig_site_etaBN,'level2');
MSPE_GPN=mean(((Pred_MeanGPN-Real_val).^2));


Real_val =  funnPerdikaris(site_e,'level2');

fighdl32=figure; 
     set(findall(gcf,'type','text'),'fontSize',25)
     set(gcf,'DefaultAxesFontSize', 25)

% xlabel(['x'], 'interpreter','latex')
% ylabel(['t_1'], 'interpreter','latex')
% zlabel(['$\eta$'], 'interpreter','latex')
set(fighdl32, 'PaperUnits', 'inches', 'PaperSize', [10 6], 'PaperPosition', [0 0 10 6]);
hold on;
plot(site_e,Real_val,'r','Linewidth',1.5);
plot(SITE1(1:(n)), data_eta(1:(n)),'b.','Linewidth',1.5);
plot(x_new(1:(m)), data_y(1:(m)),'r*','MarkerSize',10);
plot(site_e,Pred_MeanGPN,'black','Linewidth',1.5);
plot(Krig_site_etaBNQ,Pred_MeanGPNQ,'green', 'Linewidth',1.5);
%saveas(fighdl32, ['Results\SimulationS1\ExS3_nstart.pdf']); 
%saveas(fighdl32, ['Results\SimulationS1\ ExS3_nstart.pdf']); 
hold off;
end