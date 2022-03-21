addpath('./lib/BCT');
%% group-level
load('offtarget.mat');
load('sc_fa_roi246_hc.mat'); con_m=con; con_bm=con_m>0;
disc=load('sc_length_roi246_hc.mat'); disc=disc.con;
dd=load('sc_edist_roi246_hc.mat'); dd=dd.con;

% seed=1:246;
seed=[51,171]; % propagation hub

exc=union(find(offtarget),find(sum(con_m)==0));
seed=setdiff(seed,exc);

clear c_score* K y*
c_score=NaN*ones(size(con_m));
c_score1=NaN*ones(size(con_m));
c_score2=NaN*ones(size(con_m));
for r=1:length(seed)
    disp(r)
    s=seed(r);
    G=graph(con_bm);
    
    L=weight_conversion(con_m,'lengths');
    dist_m=distance_wei(L);
          
    for t=1:size(con_m,1)
        if t==s
            continue;
        end
        
        [mf,gf]=maxflow(G,s,t);
        tmp=find(gf.Edges.EndNodes(:,2)==t);
        
        A=sparse(adjacency(gf)); [l,m]=find(A);
        for n=1:length(l)
            A(l(n),m(n))=disc(l(n),m(n));
        end
        
        wei=cell(length(tmp),1);
        di=cell(length(tmp),1);
        for i=1:length(tmp)
            if i==1
                c_score(s,t)=0;
            end
            
            [~,path{i},~]=graphshortestpath(A,s,t);
            
            if length(path{i})>1
                for k=1:length(path{i})-1
                    wei{i}(k)=con_m(path{i}(k),path{i}(k+1));
                    di{i}(k)=disc(path{i}(k),path{i}(k+1));
                    A(path{i}(k),path{i}(k+1))=0;
                    A(path{i}(k+1),path{i}(k))=0;
                end
            else
                wei{i}=[];
                di{i}=[];
            end
            c_score(s,t)=c_score(s,t)+mean(wei{i})/sum(di{i});
        end
        c_score1(s,t)=dist_m(s,t); % Shortest path length
        c_score2(s,t)=dd(s,t); % Euclidean distance
        
        clear path
    end
end
clear wei di

%% correlation between network flow-based connectivity and tau W contrast
load('w_contrast_acc.mat');

clear r1 p1
for r=1:length(seed)
    s=seed(r);
    if s<=size(con_m,2)/2
        idd=1:size(con_m,2)/2;
    else
        idd=size(con_m,2)/2+1:size(con_m,2);
    end
    idd=setdiff(idd,[exc,s]);
    
    [K.r(1,s),K.p(1,s)]=corr(c_score(s,idd)',Wdelta(idd),'rows','complete');
    %     [K.r(2,s),K.p(2,s)]=corr(-c_score1(s,idd)',Wdelta(:,idd)','rows','complete');
    %     [K.r(3,s),K.p(3,s)]=corr(-c_score2(s,idd)',Wdelta(:,idd)','rows','complete');
end
