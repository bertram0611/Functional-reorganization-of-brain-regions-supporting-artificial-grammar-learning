%% Channel30
clear
clc
load('prerest_z_hbo.mat')
N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

%%%%%%%%%%%%%%%%%%%%%%
a = N_region*(N_region-1)/2;
prerest = zeros(a,N_sub);
for i = 1:N_sub
    prerest(:,i) = r2z{i}(triu(true(N_region),1));%
end

clear i r2z
load('postrest_z_hbo.mat')
post = zeros(a,N_sub);
for i = 1:N_sub
    post(:,i) = r2z{i}(triu(true(N_region),1));%
end

for j=1:a
    [h,p,ci,stats] = ttest(prerest(j,:),post(j,:));
    allh(j) = h;
    allp(j) = p;
    allt(j) = stats.tstat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch30prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch30prerest(i,:) = r2z{i}(30,:);% chan30
end

clear i r2z
load('postrest_z_hbo.mat')
ch30post = zeros(N_sub,N_region);
for i = 1:N_sub
    ch30post(i,:) = r2z{i}(30,:);% chan30
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch30post-ch30prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[5 6 7 8 9 11 12 13 14 15 18 19 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(1,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

[x3,x4] = find(pvalue<0.05);
sig_rvalue = rvalue(1,x4);
edge_corr = zeros(N_region);
edge_corr(30,x4) = sig_rvalue;% chan30
edge_corr_ch = tril(edge_corr,-1)'+triu(edge_corr,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rvalue1 = edge_corr_ch(triu(true(N_region),1));
[x5, x6] = find(rvalue1~=0);%
C = intersect(x2,x5');
rvalue_sig = rvalue1(C);

hh = zeros(1,1035);
hh(C) = rvalue_sig;
edge_corr1 = zeros(N_region);
edge_corr1(triu(true(N_region),1)) = hh;

%% Channel34
clear
clc
load('prerest_z_hbo.mat')
N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

%%%%%%%%%%%%%%%%%%%%%%
a = N_region*(N_region-1)/2;
prerest = zeros(a,N_sub);
for i = 1:N_sub
    prerest(:,i) = r2z{i}(triu(true(N_region),1));%
end

clear i r2z
load('postrest_z_hbo.mat')
post = zeros(a,N_sub);
for i = 1:N_sub
    post(:,i) = r2z{i}(triu(true(N_region),1));%
end

for j=1:a
    [h,p,ci,stats] = ttest(prerest(j,:),post(j,:));
    allh(j) = h;
    allp(j) = p;
    allt(j) = stats.tstat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch34prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch34prerest(i,:) = r2z{i}(34,:);% chan34
end

clear i r2z
load('postrest_z_hbo.mat')
ch34post = zeros(N_sub,N_region);
for i = 1:N_sub
    ch34post(i,:) = r2z{i}(34,:);% chan34
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch34post-ch34prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[5 6 7 8 9 11 12 13 14 15 18 19 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(2,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

% S6 Fig: FC34－17 and activation in Ch34
figure('color','w');  
scatter(dif_FC(:,17),dif_sigchan_incon_con_remain(2,:)',100,'b','.')% chan34
hold on
lsline
axis([-1 2 -0.04 0.08])%[xmin xmax ymin ymax]
set(gcf,'units','centimeters');
set(gcf,'position',[10 10 5 5])
% xlabel('channel'); ylabel('channel');
set(gca,'yTick',[-0.04 0 0.04 0.08]);
set(gca,'yTickLabel',{'-0.04', '0', '0.04', '0.08'});
% set(gca,'yTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
set(gca,'xTick',[-1 0 1 2]);
set(gca,'xTickLabel',{'-1', '0', '1', '2'});
% set(gca,'xTickLabel',[]);
set(gca,'FontSize',8, 'FontName','Arial')
% set(gca,'xaxislocation','top')
% colorbar('FontSize',8)
axis square


[x3,x4] = find(pvalue<0.05);
sig_rvalue = rvalue(1,x4);
edge_corr = zeros(N_region);
edge_corr(34,x4) = sig_rvalue;% chan34
edge_corr_ch = tril(edge_corr,-1)'+triu(edge_corr,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rvalue1 = edge_corr_ch(triu(true(N_region),1));
[x5, x6] = find(rvalue1~=0);%
C = intersect(x2,x5');
rvalue_sig = rvalue1(C);

hh = zeros(1,1035);
hh(C) = rvalue_sig;
edge_corr1 = zeros(N_region);
edge_corr1(triu(true(N_region),1)) = hh;

%% Channel35
clear
clc
load('prerest_z_hbo.mat')
N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

%%%%%%%%%%%%%%%%%%%%%%
a = N_region*(N_region-1)/2;
prerest = zeros(a,N_sub);
for i = 1:N_sub
    prerest(:,i) = r2z{i}(triu(true(N_region),1));%
end

clear i r2z
load('postrest_z_hbo.mat')
post = zeros(a,N_sub);
for i = 1:N_sub
    post(:,i) = r2z{i}(triu(true(N_region),1));%
end

for j=1:a
    [h,p,ci,stats] = ttest(prerest(j,:),post(j,:));
    allh(j) = h;
    allp(j) = p;
    allt(j) = stats.tstat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch35prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch35prerest(i,:) = r2z{i}(35,:);% chan35
end

clear i r2z
load('postrest_z_hbo.mat')
ch35post = zeros(N_sub,N_region);
for i = 1:N_sub
    ch35post(i,:) = r2z{i}(35,:);% chan35
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch35post-ch35prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[5 6 7 8 9 11 12 13 14 15 18 19 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(3,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

[x3,x4] = find(pvalue<0.05);
sig_rvalue = rvalue(1,x4);
edge_corr = zeros(N_region);
edge_corr(35,x4) = sig_rvalue;% chan35
edge_corr_ch = tril(edge_corr,-1)'+triu(edge_corr,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rvalue1 = edge_corr_ch(triu(true(N_region),1));
[x5, x6] = find(rvalue1~=0);%
C = intersect(x2,x5');
rvalue_sig = rvalue1(C);

hh = zeros(1,1035);
hh(C) = rvalue_sig;
edge_corr1 = zeros(N_region);
edge_corr1(triu(true(N_region),1)) = hh;

%% Channel38
clear
clc
load('prerest_z_hbo.mat')
N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

%%%%%%%%%%%%%%%%%%%%%%
a = N_region*(N_region-1)/2;
prerest = zeros(a,N_sub);
for i = 1:N_sub
    prerest(:,i) = r2z{i}(triu(true(N_region),1));%
end

clear i r2z
load('postrest_z_hbo.mat')
post = zeros(a,N_sub);
for i = 1:N_sub
    post(:,i) = r2z{i}(triu(true(N_region),1));%
end

for j=1:a
    [h,p,ci,stats] = ttest(prerest(j,:),post(j,:));
    allh(j) = h;
    allp(j) = p;
    allt(j) = stats.tstat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch38prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch38prerest(i,:) = r2z{i}(38,:);% chan38
end

clear i r2z
load('postrest_z_hbo.mat')
ch38post = zeros(N_sub,N_region);
for i = 1:N_sub
    ch38post(i,:) = r2z{i}(38,:);% chan38
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch38post-ch38prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[5 6 7 8 9 11 12 13 14 15 18 19 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(4,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

[x3,x4] = find(pvalue<0.05);
sig_rvalue = rvalue(1,x4);
edge_corr = zeros(N_region);
edge_corr(38,x4) = sig_rvalue;% chan38
edge_corr_ch = tril(edge_corr,-1)'+triu(edge_corr,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rvalue1 = edge_corr_ch(triu(true(N_region),1));
[x5, x6] = find(rvalue1~=0);%
C = intersect(x2,x5');
rvalue_sig = rvalue1(C);

hh = zeros(1,1035);
hh(C) = rvalue_sig;
edge_corr1 = zeros(N_region);
edge_corr1(triu(true(N_region),1)) = hh;

%% Channel39
clear
clc
load('prerest_z_hbo.mat')
N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

%%%%%%%%%%%%%%%%%%%%%%
a = N_region*(N_region-1)/2;
prerest = zeros(a,N_sub);
for i = 1:N_sub
    prerest(:,i) = r2z{i}(triu(true(N_region),1));%
end

clear i r2z
load('postrest_z_hbo.mat')
post = zeros(a,N_sub);
for i = 1:N_sub
    post(:,i) = r2z{i}(triu(true(N_region),1));%
end

for j=1:a
    [h,p,ci,stats] = ttest(prerest(j,:),post(j,:));
    allh(j) = h;
    allp(j) = p;
    allt(j) = stats.tstat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch39prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch39prerest(i,:) = r2z{i}(39,:);% chan39
end

clear i r2z
load('postrest_z_hbo.mat')
ch39post = zeros(N_sub,N_region);
for i = 1:N_sub
    ch39post(i,:) = r2z{i}(39,:);% chan39
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch39post-ch39prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[5 6 7 8 9 11 12 13 14 15 18 19 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(5,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

[x3,x4] = find(pvalue<0.05);
sig_rvalue = rvalue(1,x4);
edge_corr = zeros(N_region);
edge_corr(39,x4) = sig_rvalue;% chan39
edge_corr_ch = tril(edge_corr,-1)'+triu(edge_corr,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rvalue1 = edge_corr_ch(triu(true(N_region),1));
[x5, x6] = find(rvalue1~=0);%
C = intersect(x2,x5');
rvalue_sig = rvalue1(C);

hh = zeros(1,1035);
hh(C) = rvalue_sig;
edge_corr1 = zeros(N_region);
edge_corr1(triu(true(N_region),1)) = hh;

%% Channel43
clear
clc
load('prerest_z_hbo.mat')
N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

%%%%%%%%%%%%%%%%%%%%%%
a = N_region*(N_region-1)/2;
prerest = zeros(a,N_sub);
for i = 1:N_sub
    prerest(:,i) = r2z{i}(triu(true(N_region),1));%
end

clear i r2z
load('postrest_z_hbo.mat')
post = zeros(a,N_sub);
for i = 1:N_sub
    post(:,i) = r2z{i}(triu(true(N_region),1));%
end

for j=1:a
    [h,p,ci,stats] = ttest(prerest(j,:),post(j,:));
    allh(j) = h;
    allp(j) = p;
    allt(j) = stats.tstat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch43prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch43prerest(i,:) = r2z{i}(43,:);% chan43
end

clear i r2z
load('postrest_z_hbo.mat')
ch43post = zeros(N_sub,N_region);
for i = 1:N_sub
    ch43post(i,:) = r2z{i}(43,:);% chan43
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch43post-ch43prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[5 6 7 8 9 11 12 13 14 15 18 19 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(6,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

[x3,x4] = find(pvalue<0.05);
sig_rvalue = rvalue(1,x4);
edge_corr = zeros(N_region);
edge_corr(43,x4) = sig_rvalue;% chan43
edge_corr_ch = tril(edge_corr,-1)'+triu(edge_corr,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rvalue1 = edge_corr_ch(triu(true(N_region),1));
[x5, x6] = find(rvalue1~=0);%
C = intersect(x2,x5');
rvalue_sig = rvalue1(C);

hh = zeros(1,1035);
hh(C) = rvalue_sig;
edge_corr1 = zeros(N_region);
edge_corr1(triu(true(N_region),1)) = hh;


%% S6 Fig: Relation between Pre-Rest and Post-Rest minus Pre-Rest for FC strength of Ch17-Ch34
clear
clc
load('prerest_z_hbo.mat')

 
N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

prerest = zeros(1,N_sub);% total 10 FCs
for i = 1:N_sub
    prerest(1,i) = r2z{i}(17,34);%
end

clear i r2z
load('postrest_z_hbo.mat')
learn = zeros(1,N_sub);
for i = 1:N_sub
    postrest(1,i) = r2z{i}(17,34);%  
end

mat = zeros(N_sub,2);

mat(:,1) = prerest(1,:)';
mat(:,2) = postrest(1,:)';

%%
FC_post_pre = mat(:,2)-mat(:,1);
[rvalue, pvalue] = corr(mat(:,1),FC_post_pre(:,1));


figure('color','w');  
scatter(mat(:,1),FC_post_pre(:,1),100,'b','.')
hold on
lsline
axis([-2 2 -1 2])%[xmin xmax ymin ymax]
set(gcf,'units','centimeters');
set(gcf,'position',[10 10 5 5])
% xlabel('channel'); ylabel('channel');
set(gca,'yTick',[-1 0 1 2]);
set(gca,'yTickLabel',{'-1', '0', '1', '2'});
% set(gca,'yTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
set(gca,'xTick',[-2 -1 0 1 2]);
set(gca,'xTickLabel',{'-2','-1', '0', '1', '2'});
% set(gca,'xTickLabel',[]);
set(gca,'FontSize',8, 'FontName','Arial')
% set(gca,'xaxislocation','top')
% colorbar('FontSize',8)
axis square













