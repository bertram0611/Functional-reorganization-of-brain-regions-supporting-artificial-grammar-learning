%% Analyses of neonates' functional connectivities 
%% Fig3A: One-sample t-tests were applied to find the significant FCs for each phase (Pre-Rest, Learning, and Post-Rest phases) against a zero baseline
%% Pearson correlation For the Pre-Rest phase
clear all
clc
load TC_hbo_pre.mat
TC_hbo(2:3)=[];
N_sub = size(TC_hbo,2);
N_region = size(TC_hbo{1},2);
%
rest = cell(N_sub,1);
p_rest = cell(N_sub,1);
for i = 1:N_sub
    [rest{i},p_rest{i}]=corr(TC_hbo{i});%CORR(X) returns a P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
end
%% Fisher z transformation
clear i
for i = 1:N_sub
    for j=1:N_region
        rest{i,1}(j,j) = 0;
    end
end
clear i j
r2z = cell(N_sub,1);
for i = 1:N_sub
    r2z{i,1} = atanh(rest{i,1});
end
% save rest_z_hbo r2z
%% zero t-test to find the significant FC for the Pre-Rest phase
clear i
for i = 1:N_sub
    r2z_upper(i,:)=r2z{i,1}(triu(true(N_region),1));
end
clear j
stats=cell(1035,1);
for j = 1:1035
    [h(j),p(j),ci,stats{j}] = ttest(r2z_upper(:,j));
    tvalue(j) = stats{j}.tstat;
end
find(tvalue<=0)% to judge whether there are negative correlations
% I found no negative correlations

%&&&&& to plot redundant FC during the Pre-Rest phase
index4 = find((p<0.05/1035));
h4 = zeros(1,1035);
h4(index4) = 1;

aa = zeros(N_region);
aa(triu(true(N_region),1)) = h4;
% copy aa to a txt file named "Edge46AGneonatePreRestManyFC.edge"

%&&&&&
index0 = find(p<0.000005/1035);
h0 = zeros(1,1035);
h0(index0) = 1;

index1 = find(p<0.0000005/1035);
h1 = zeros(1,1035);
h1(index1) = 1;

index2 = find(p<0.00000005/1035);
h2 = zeros(1,1035);
h2(index2) = 1;

index3 = find(p<0.000000005/1035);
h3 = zeros(1,1035);
h3(index3) = 1;
hh = h0+h1+h2+h3;
a = zeros(N_region);
a(triu(true(N_region),1)) = hh;
% copy a to a txt file named "Edge46AGneonatePreRest.edge"
% use the edge files to plot Fig.3A using BrainNet Viewer 


%% Pearson correlation For the Learning phase
clear all
clc
load TC_hbo_learn.mat
TC_hbo([5 12 17])=[]; % exlude 3 babies 
N_sub = size(TC_hbo,2);
N_region = size(TC_hbo{1},2);
%
rest = cell(N_sub,1);
p_rest = cell(N_sub,1);
for i = 1:N_sub
    [rest{i},p_rest{i}]=corr(TC_hbo{i});%CORR(X) returns a P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
end
%% Fisher z transformation
clear i
for i = 1:N_sub
    for j=1:N_region
        rest{i,1}(j,j) = 0;
    end
end
clear i j
r2z = cell(N_sub,1);
for i = 1:N_sub
    r2z{i,1} = atanh(rest{i,1});
end
% save rest_z_hbo r2z
%% zero t-test to find the significant FC for the Learning phase
clear i
for i = 1:N_sub
    r2z_upper(i,:)=r2z{i,1}(triu(true(N_region),1));
end
clear j
stats=cell(1035,1);
for j = 1:1035
    [h(j),p(j),ci,stats{j}] = ttest(r2z_upper(:,j));
    tvalue(j) = stats{j}.tstat;
end
find(tvalue<=0)% to judge whether there are negative correlations
% I found no negative correlations

%&&&&& to plot redundant FC during the Learning phase
index4 = find((p<0.05/1035));
h4 = zeros(1,1035);
h4(index4) = 1;

aa = zeros(N_region);
aa(triu(true(N_region),1)) = h4;
% copy aa to a txt file named "Edge46AGneonateLearnManyFC.edge"

%&&&&&
index0 = find(p<0.000005/1035);
h0 = zeros(1,1035);
h0(index0) = 1;

index1 = find(p<0.0000005/1035);
h1 = zeros(1,1035);
h1(index1) = 1;

index2 = find(p<0.00000005/1035);
h2 = zeros(1,1035);
h2(index2) = 1;

index3 = find(p<0.000000005/1035);
h3 = zeros(1,1035);
h3(index3) = 1;
hh = h0+h1+h2+h3;
a = zeros(N_region);
a(triu(true(N_region),1)) = hh;
% copy a to a txt file named "Edge46AGneonateLearn.edge"
% use the edge files to plot Fig.3A using BrainNet Viewer 


%% Pearson correlation For the Post-Rest phase
clear all
clc
load TC_hbo_post.mat
% due to nan values in subject No4, so exclude it
TC_hbo([1 4 15 16 19])=[];
N_sub = size(TC_hbo,2);
N_region = size(TC_hbo{1},2);
%
rest = cell(N_sub,1);
p_rest = cell(N_sub,1);
for i = 1:N_sub
    [rest{i},p_rest{i}]=corr(TC_hbo{i});%CORR(X) returns a P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
end
%% Fisher z transformation and average FCs across subjects
clear i
for i = 1:N_sub
    for j=1:N_region
        rest{i,1}(j,j) = 0;
    end
end
clear i j
r2z = cell(N_sub,1);
for i = 1:N_sub
    r2z{i,1} = atanh(rest{i,1});
end
% save rest_z_hbo r2z
%% zero t-test to find the significant FC for the Post-Rest phase
clear i
for i = 1:N_sub
    r2z_upper(i,:)=r2z{i,1}(triu(true(N_region),1));
end
clear j
stats=cell(1035,1);
for j = 1:1035
    [h(j),p(j),ci,stats{j}] = ttest(r2z_upper(:,j));
    tvalue(j) = stats{j}.tstat;
end
find(tvalue<=0)% to judge whether there are negative correlations
% I found no negative correlations

%&&&&& to plot redundant FC during the Post-Rest phase
index4 = find((p<0.05/1035));
h4 = zeros(1,1035);
h4(index4) = 1;

aa = zeros(N_region);
aa(triu(true(N_region),1)) = h4;
% copy aa to a txt file named "Edge46AGneonatePostRestManyFC.edge"

%&&&&&
index0 = find(p<0.000005/1035);
h0 = zeros(1,1035);
h0(index0) = 1;

index1 = find(p<0.0000005/1035);
h1 = zeros(1,1035);
h1(index1) = 1;

index2 = find(p<0.00000005/1035);
h2 = zeros(1,1035);
h2(index2) = 1;

index3 = find(p<0.000000005/1035);
h3 = zeros(1,1035);
h3(index3) = 1;
hh = h0+h1+h2+h3;
a = zeros(N_region);
a(triu(true(N_region),1)) = hh;
% copy a to a txt file named "Edge46AGneonatePostRest.edge"
% use the edge files to plot Fig.3A using BrainNet Viewer 


%% Fig.3B: Compare FC differences for the contrasts of Learning vs. Pre-Rest
% Data Preparation
% Pre-Rest Pearson correlation
clear all
clc
load TC_hbo_pre.mat
TC_hbo([2 3 5 12 17])=[];
N_sub = size(TC_hbo,2);
N_region = size(TC_hbo{1},2);
%
rest = cell(N_sub,1);
p_rest = cell(N_sub,1);
for i = 1:N_sub
    [rest{i},p_rest{i}]=corr(TC_hbo{i});%CORR(X) returns a P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
end
% Fisher z transformation and average FCs across subjects
clear i
for i = 1:N_sub
    for j=1:N_region
        rest{i,1}(j,j) = 0;
    end
end
clear i j
r2z = cell(N_sub,1);
for i = 1:N_sub
    r2z{i,1} = atanh(rest{i,1});
end
save prerest_z_hbo r2z

% Learning Pearson correlation
clear all
clc
load TC_hbo_learn.mat
TC_hbo([2 3 5 12 17])=[];
N_sub = size(TC_hbo,2);
N_region = size(TC_hbo{1},2);
%
rest = cell(N_sub,1);
p_rest = cell(N_sub,1);
for i = 1:N_sub
    [rest{i},p_rest{i}]=corr(TC_hbo{i});%CORR(X) returns a P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
end
% Fisher z transformation and average FCs across subjects
clear i
for i = 1:N_sub
    for j=1:N_region
        rest{i,1}(j,j) = 0;
    end
end
clear i j
r2z = cell(N_sub,1);
for i = 1:N_sub
    r2z{i,1} = atanh(rest{i,1});
end
save learn_z_hbo r2z


%% Compute FC differences between Pre-Rest and Learning
clear
clc
load('prerest_z_hbo.mat')

 
N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

a = N_region*(N_region-1)/2;
prerest = zeros(a,N_sub);
for i = 1:N_sub
    prerest(:,i) = r2z{i}(triu(true(N_region),1));%
end

clear i r2z
load('learn_z_hbo.mat')
learn = zeros(a,N_sub);
for i = 1:N_sub
    learn(:,i) = r2z{i}(triu(true(N_region),1));%
end

for j=1:a
    [h,p,ci,stats] = ttest(prerest(j,:),learn(j,:));
    allh(j) = h;
    allp(j) = p;
    allt(j) = stats.tstat;
end

% FDR correction
[h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(allp,0.05,'pdep','yes');
%% This part is for BrainNet, only show edges with p<0.05, which show the T values 
[r,rr] = find(allp>=0.05);
% [r,rr] = find(allp>=0.01);
% [r,rr] = find(allp>=0.001);
allt(rr)=0;
a2 = zeros(N_region);
a2(triu(true(N_region),1)) = allt;
% copy a2 to a txt file named "Edge46AGneonateLearnVSPreTvalue.edge"
%%
[r,rr] = find(allh==1);
pp = allp(rr);
tt = allt(rr);

r_pos = rr(tt>0);% positive represents prerest larger than learn 
r_neg = rr(tt<0);% negative represents prerest smaller than learn

%
a_pos = zeros(1,a);
a_neg = zeros(1,a);
a_pos(r_pos) = -1;
a_neg(r_neg) = 1;
a2 = zeros(N_region);
a2(triu(true(N_region),1)) = a_pos;
% pos_matrix = (a2+a2') - eye(size(a2,1)).*diag(a2);% to make a symmetric matrix

a3 = zeros(N_region);
a3(triu(true(N_region),1)) = a_neg;
neg_matrix = (a3+a3') - eye(size(a3,1)).*diag(a3);% to make a symmetric matrix
lowertri = neg_matrix(tril(true(N_region),-1));
a2(tril(true(N_region),-1)) = lowertri;
% the lower part of a2 indicates negative edges while upper part positive
% edges


%% Fig.3B: Compare FC differences for the contrasts of Post-Rest vs. Pre-Rest
% Data Preparation
% PreRest Pearson correlation
clear all
clc
load TC_hbo_pre.mat
TC_hbo([1 2 3 4 15 16 19])=[];% Don't forget sub04, due to its nan values in some channels 
N_sub = size(TC_hbo,2);
N_region = size(TC_hbo{1},2);
%
rest = cell(N_sub,1);
p_rest = cell(N_sub,1);
for i = 1:N_sub
    [rest{i},p_rest{i}]=corr(TC_hbo{i});%CORR(X) returns a P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
end

% Fisher z transformation and average FCs across subjects
clear i
for i = 1:N_sub
    for j=1:N_region
        rest{i,1}(j,j) = 0;
    end
end
clear i j
r2z = cell(N_sub,1);
for i = 1:N_sub
    r2z{i,1} = atanh(rest{i,1});
end
save prerest_z_hbo r2z

% PostRest Pearson correlation
clear all
clc
load TC_hbo_post.mat
TC_hbo([1 2 3 4 15 16 19])=[];% Don't forget sub04, due to its nan values in some channels 
N_sub = size(TC_hbo,2);
N_region = size(TC_hbo{1},2);
%
rest = cell(N_sub,1);
p_rest = cell(N_sub,1);
for i = 1:N_sub
    [rest{i},p_rest{i}]=corr(TC_hbo{i});%CORR(X) returns a P-by-P matrix containing the pairwise linear correlation coefficient between each pair of columns in the N-by-P matrix X.
end

% Fisher z transformation and average FCs across subjects
clear i
for i = 1:N_sub
    for j=1:N_region
        rest{i,1}(j,j) = 0;
    end
end
clear i j
r2z = cell(N_sub,1);
for i = 1:N_sub
    r2z{i,1} = atanh(rest{i,1});
end
save postrest_z_hbo r2z

%% Compute FC differences between Pre-Rest and Post-Rest
clear
clc
load('prerest_z_hbo.mat')

 
N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

a = N_region*(N_region-1)/2;
prerest = zeros(a,N_sub);
for i = 1:N_sub
    prerest(:,i) = r2z{i}(triu(true(N_region),1));%
end

clear i r2z
load('postrest_z_hbo.mat')
postrest = zeros(a,N_sub);
for i = 1:N_sub
    postrest(:,i) = r2z{i}(triu(true(N_region),1));%
end

for j=1:a
    [h,p,ci,stats] = ttest(prerest(j,:),postrest(j,:));
    allh(j) = h;
    allp(j) = p;
    allt(j) = stats.tstat;
end
% FDR correction
[h1, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(allp,.05,'pdep','yes');
%% This part is for BrainNet, only show edges with p<0.005, which show the T values 
[r,rr] = find(allp>=0.05);
% [r,rr] = find(allp>=0.01);
% [r,rr] = find(allp>=0.001);
allt(rr)=0;
a2 = zeros(N_region);
a2(triu(true(N_region),1)) = allt;
% copy a2 to a text file named "Edge46AGneonatePostVSPreTvalue.edge"
%%
% [p1,p2] = sort(allp);%
% p3 =0.05/a:0.05/a:0.05;
% d =find((p1-p3)<0,1,'last');%
% dd = p2(1:d);

%
[r,rr] = find(allh==1);
pp = allp(rr);
tt = allt(rr);

r_pos = rr(tt>0);% positive represents prerest larger than postrest 
r_neg = rr(tt<0);% negative represents prerest smaller than postrest

%
a_pos = zeros(1,a);
a_neg = zeros(1,a);
a_pos(r_pos) = -1;
a_neg(r_neg) = 1;
a2 = zeros(N_region);
a2(triu(true(N_region),1)) = a_pos;
% pos_matrix = (a2+a2') - eye(size(a2,1)).*diag(a2);% to make a symmetric matrix

a3 = zeros(N_region);
a3(triu(true(N_region),1)) = a_neg;
neg_matrix = (a3+a3') - eye(size(a3,1)).*diag(a3);% to make a symmetric matrix
lowertri = neg_matrix(tril(true(N_region),-1));
a2(tril(true(N_region),-1)) = lowertri;
% the lower part of a2 indicates negative edges while upper part positive
% edges

