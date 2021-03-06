function [Zt_real,site_real,prCurr,AA,tau2, a_new,lim_all, KrigS,Var_krigS,Krig_site, Krig_W_etaMAT, Krig_site_eta]=...
    ABTCK_B(PR_CURR,SITE_real,ZT_real,n,AA_E,LIM_ALL,a_new,a_NODES,NODES,pprop, aalpha, bbeta,S,low_l,site_e,n_e,BB,Delta,tau1,TAU2,bound,funname)
KrigS.mean= zeros(BB,n_e); KrigS.Low= zeros(BB,n_e);
KrigS.Up= zeros(BB,n_e); Var_krigS= zeros(BB,n_e); 
funname2='constant';
%funname2='linear';
%funname2='quadratic';
%funname2='polyfour';%
%funname2='polythree';
SITE{1}{1}=[SITE_real{1}{1}; SITE_real{1}{2}];
ZT{1}{1}=[ZT_real{1}{1}; ZT_real{1}{2}];%+norm_r(0,0.02,size(ZT_real{1}{2},1))];
SITE{1}{2}=[SITE_real{1}{2}];
ZT{1}{2}=[ZT_real{1}{2}];

Krig_site_eta=[]; Krig_W_etaMAT=zeros(BB,size(site_e,1));  
Dim_x =size(SITE{1}{1},2); Krig_site=zeros(n_e,Dim_x); mm=sqrt(size(site_e,1))-1;%29;%

for i=1:BB
    % ZT_NEW, SITE_NEW,  LIM_ALL, PR_CURR FOR ALL THE NODES. 
    % Zt_new AND site_new FOR EXTERNAL NODES.
    operation = randi(5);
    %operation
    if (operation ==1)
    % Split 
    [SITE,ZT, SITE_real,ZT_real,  n, PR_CURR, AA_E,TAU2,LIM_ALL, a_new,NODES,a_NODES]= ...
    Split_Nonnest(PR_CURR,SITE,ZT,SITE_real,ZT_real,n,AA_E,S,low_l,LIM_ALL,a_new,a_NODES,NODES, pprop, aalpha, bbeta,tau1,TAU2,funname,funname2);
    elseif  (operation ==2)
            % Merge 
    if a_NODES>2 
        [SITE,ZT,SITE_real,ZT_real, n, PR_CURR, AA_E,TAU2,a_new, LIM_ALL,NODES,a_NODES]=...
    Merge_NON_nest(PR_CURR,SITE,ZT,SITE_real,ZT_real,n,AA_E,S,LIM_ALL,a_new,a_NODES,NODES, aalpha, bbeta,tau1,TAU2,funname,funname2);
    end
    elseif  (operation >2 && operation < 5)
    % Change 
    if a_NODES>2     
          [SITE,ZT,SITE_real,ZT_real, n, LIM_ALL, PR_CURR,NODES]= ...
        Change_Non_nest(PR_CURR,SITE,ZT,SITE_real,ZT_real,n,AA_E,S,low_l,LIM_ALL, a_NODES,NODES,tau1,TAU2,funname,funname2);                  
    end
    elseif  (operation ==5)
        %Swap-Rotate
    if a_NODES>2
        [SITE,ZT,SITE_real,ZT_real, n, lim_ALL, prCurr,NODES]=Swap_Nonnest(PR_CURR,SITE,ZT,SITE_real,ZT_real,n,AA_E,S,low_l,LIM_ALL, a_NODES,NODES,tau1,TAU2,funname);
    end
    end 
    %a_NODES
   for inter=1:a_NODES
   for t=1:S
    [YY,II]=sort(SITE{inter}{t}(:,1));
    SITE{inter}{t}=SITE{inter}{t}(II,:);     
    ZT{inter}{t}=ZT{inter}{t}(II,:); 
    %operation
    %SITE_real{inter}
    [YYy,IIi]=sort(SITE_real{inter}{t}(:,1));
    SITE_real{inter}{t}=SITE_real{inter}{t}(IIi,:);     
    ZT_real{inter}{t}=ZT_real{inter}{t}(IIi,:); 
   end 
   end
    
   
    % Find the external NODES so we can make the GP updates.  
   clear Zt; clear AA; clear beta_X; clear site; clear lim_all; clear tau2;
   clear site_real; clear Zt_real; clear site_add; clear New_Z;
    hh2=0; 
    
    for ii=1:a_NODES
        if (NODES{ii}.inter==0)
            hh2=hh2+1; fin_nodes(hh2)=NODES{ii};
            depth_each(hh2)=floor(NODES{ii}.number/2);
            prCurr{hh2}= PR_CURR{ii}; 
            Zt{hh2}= ZT{ii}; site{hh2}=SITE{ii}; 
            %size(Zt{hh2}{1})
            Zt_real{hh2}= ZT_real{ii}; site_real{hh2}=SITE_real{ii}; 
            AA{hh2}= AA_E{ii}; tau2{hh2}= TAU2{ii}; 
            lim_all{hh2}=LIM_ALL{ii}; 
            %nn{hh2}=n{ii};
            % GP update  
           for jj=1:1 % to inprouve convergency we may use 3 or 4;     
               [AA{hh2}, sigma2{hh2}, tau2{hh2}, prCurr{hh2}, mu{hh2}, beta_X{hh2}]= Multi_MH_s(Zt{hh2},site{hh2},AA{hh2},S,Delta,tau1,tau2{hh2},bound,funname2);
                %[AA{hh2}, sigma2{hh2}, tau2{hh2}, prCurr{hh2}, mu{hh2}, beta_X{hh2}]= Multi_MH_s2(Zt{hh2},site{hh2},AA{hh2},S,Delta,tau1,tau2{hh2});     
           end
            PR_CURR{ii}=prCurr{hh2};  AA_E{ii}=AA{hh2}; TAU2{ii}=tau2{hh2};  LIM_ALL{ii}=lim_all{hh2}; 
            beta_X{hh2}{1}(2)=1;
            for t=1:S 
            beta{hh2}{t}=beta_X{hh2}{t}(1);
            rho{hh2}{t}=beta_X{hh2}{t}(2:size(beta_X{hh2}{t}));
            %beta{hh2}{t}=beta_X{hh2}{t}(1);
            %rho{hh2}{t}=beta_X{hh2}{t}(2);
            %[hh2 beta_X{hh2}{t}(1) beta_X{hh2}{t}(2)]
            end
             %Beta_X{hh2}{1}=beta_X{hh2}{1}(1); Beta_X{hh2}{2}=beta_X{hh2}{2};
            %rho{hh2}
            % Prediction of the unobserved data
            for tt=1:(S-1)
            site_add{hh2}=site_real{hh2}{tt+1};
            site{hh2}{tt}=[site_real{hh2}{tt}; site_add{hh2}];
            site{hh2}{tt+1}=[site_real{hh2}{tt+1}];
               [New_Z{hh2}, Var_Z{hh2}]=Krig_MultiB(Zt_real{hh2},beta{hh2},rho{hh2},AA{hh2},site_real{hh2},site_add{hh2},tt,tau2{hh2},funname2);
                %Zt{hh2}{tt}=[Zt_real{hh2}{tt}; New_Z{hh2}{tt}];
                Zt{hh2}{tt}=[Zt_real{hh2}{tt}; New_Z{hh2}];
                Zt{hh2}{tt+1}=Zt_real{hh2}{tt+1};
                %size(Zt{hh2}{tt})
            end
            SITE{ii}= site{hh2};
            ZT{ii}=Zt{hh2};
            %size(ZT{ii}{1})
        end
    end     
       

M=1;
[site_e_t]= RE_site(site_e,lim_all,a_new);   
%size(site_e_t)
[Kriging_W,LB_W, UB_W, site_W]=KM_Augmented_TreeB(Zt,beta,rho,sigma2,AA,site,site_e_t,a_new,S,tau2,funname2);
XY_eta= [site_W, Kriging_W,LB_W, UB_W]; 
XY_11_eta=sortrows(XY_eta,2);
XY_12_eta=sortrows(XY_11_eta,1);
Kriging_W_etaS=XY_12_eta(:,(Dim_x+1):(Dim_x+M));
KrigS.mean(i,:)=Kriging_W_etaS;
KrigS.Low(i,:)=XY_12_eta(:,(Dim_x+2));
KrigS.Up(i,:)=XY_12_eta(:,(Dim_x+3));
Krig_W_etaMAT(i,:)=Kriging_W_etaS;%
Krig_site_eta=XY_12_eta(:,1:Dim_x);
    
predim1=2;
%figure
XX1=[]; XX2=[]; KW_U=[];
for jj=1:(mm+1)
     KW_U=[KW_U Kriging_W_etaS((1+(jj-1)*(mm+1)):(jj*(mm+1)))];
     XX1=[XX1 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),1)];
     XX2=[XX2 Krig_site_eta((1+(jj-1)*(mm+1)):(jj*(mm+1)),predim1)];    
end
if mod(i,4000)<=0  
    a_new
%  for kk=1:a_NODES
%      NODES{kk}.direction
%       NODES{kk}.rul
%  end     
imagesc(XX1(:),XX2(:),KW_U); colorbar; set(gca,'YDir','normal'); %caxis([-1.5 3]);     
 hold on;
for ii=1:a_new
    l1=lim_all{ii}(1,1); l2=lim_all{ii}(1,2);        
    %l3=lim_all{ii}((Dim_x-p),1); l4=lim_all{ii}((Dim_x-p),2);
    l3=lim_all{ii}((Dim_x),1); l4=lim_all{ii}((Dim_x),2);
    %text(xcurr(1,i), xcurr(2,i),num2str(i), )
    text((l1+l2)/2,(l3+l4)/2,num2str(ii), 'Color','r');
    plot([l1 l1],[l3 l4]);hold on; plot([l2 l2],[l3 l4]);    
    plot([l1 l2], [l3 l3]); plot([l1 l2],[l4 l4]);
end    
hold off; pause(.0001); % <- a time consuming OP
fprintf(1,'iter %d\n',i);
end    
end
end 
% hold off; 