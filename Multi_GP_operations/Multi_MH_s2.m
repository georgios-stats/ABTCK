function [AA, SIG, tau2, prCurr, mu, beta_X]= Multi_MH_s2(Zt,site,AA,S,Delta,tau1,tau2)
% MH within Gibbs sampling.   
up_bound=30; ub1TAU2=0.1;ub2TAU2=0.25;
Dim_x=size(site{1},2); %nn=cell(S,1); beta_X=cell(S,1); prCurr=cell(S,1); 
nn{1}=size(site{1},1);
    X{1}=[ones(nn{1},1)];
    [mu{1}, SIG{1}, beta_X{1}, prCurr{1}]=Multi_mu_s3(Zt{1},AA{1}, site{1},X{1},tau1,tau2{1});
for k1=1:Dim_x     
    [AA{1}, prCurr{1}]=Mult_Diag_s(Zt{1}, prCurr{1}, SIG{1},AA{1}, site{1},X{1},k1,Delta,tau1,tau2{1},up_bound) ;            
end
[tau2{1}, prCurr{1}]=Mult_tau2_s(Zt{1}, prCurr{1}, SIG{1},AA{1}, site{1},X{1},Delta,tau1,tau2{1},ub1TAU2);
for t=2:S
    %Constant mean only 
    nn{t}=size(site{t},1);
     X1 = [ones(nn{t},1)] ;
    ISM=ismember(site{t-1},site{t});
     X2=[Zt{t-1}(ISM(:,1))];

     Gt{t}=[X1  X2];
    [mu{t}, SIG{t}, beta_X{t}, prCurr{t}]=Multi_mu_s3(Zt{t},AA{t}, site{t},Gt{t},tau1,tau2{t});
for k1=1:Dim_x   
    up_bound=24;
    [AA{t}, prCurr{t}]=Mult_Diag_s(Zt{t}, prCurr{t}, SIG{t},AA{t}, site{t},Gt{t},k1,Delta,tau1,tau2{t},up_bound) ;  
end
[tau2{t}, prCurr{t}]=Mult_tau2_s(Zt{t}, prCurr{t}, SIG{t},AA{t}, site{t},Gt{t},Delta,tau1,tau2{t},ub2TAU2);
end

end
