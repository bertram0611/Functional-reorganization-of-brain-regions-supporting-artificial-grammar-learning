%% Task data preparation
% 6 prefrontal channels that were significantly activated during the Test
% phase (Ch 30, Ch 34, Ch 35, Ch 38, Ch 39, Ch 43, termed seed channels from here onwards)
clear all
clc
load incongruent_allsub_2probe_final.mat
load congruent_allsub_2probe_final.mat

% for congruent condition
N_sub = length(congruent_allsub_2probe_final);
N_chan = size(congruent_allsub_2probe_final{1},1);
reprevalue_allchan_con = zeros(N_chan,1);
reprevalue_allsub_con = zeros(N_chan,N_sub);
for i = 1:N_sub
    for j = 1:N_chan
        reprevalue_allchan_con(j,:) = mean(congruent_allsub_2probe_final{i}(j,125:175),2);% early window
    end
    reprevalue_allsub_con(:,i) = reprevalue_allchan_con;
end

% for incongruent condition
clear i j
reprevalue_allchan_incon = zeros(N_chan,1);
reprevalue_allsub_incon = zeros(N_chan,N_sub);
for i = 1:N_sub
    for j = 1:N_chan
        reprevalue_allchan_incon(j,:) = mean(incongruent_allsub_2probe_final{i}(j,125:175),2);% early window
    end
    reprevalue_allsub_incon(:,i) = reprevalue_allchan_incon;
end

% if we conduct permutation test channel by channel, after removing rows with nan values
pval_all = zeros(N_chan,1);
t_orig_all = zeros(N_chan,1);
crit_t_all = cell(N_chan,1);
clear k array
for k = 1:N_chan
    array = [reprevalue_allsub_incon(k,:)' reprevalue_allsub_con(k,:)'];
    array(any(isnan(array),2),:) = [];% to delete the rows in which either column 1 or 2 is NaN
    dif = array(:,1)-array(:,2);
    [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(dif,50000);
    pval_all(k,1) = pval;
    t_orig_all(k,1) = t_orig;
    crit_t_all{k,1} = crit_t;
end

% 6 seed channels
[x1, x2] = find(pval_all<0.05);
sigchan_con = reprevalue_allsub_con(x1,:);
sigchan_incon = reprevalue_allsub_incon(x1,:);
dif_sigchan_incon_con = sigchan_incon-sigchan_con;
% supplymentary analysis for bilateral SMG(Ch2,5,13,15)
smg=[2 5 13 15];
sigchan_con_smg = reprevalue_allsub_con(smg,:);
sigchan_incon_smg = reprevalue_allsub_incon(smg,:);
dif_sigchan_incon_con_smg = sigchan_incon_smg-sigchan_con_smg;

save dif_sigchan_incon_con dif_sigchan_incon_con
save dif_sigchan_incon_con_smg dif_sigchan_incon_con_smg

%% Analyze the correlations between the strength of FCs with significant changes during the Learning phase (i.e., Learning minus Pre-Rest)
% and the degree of activation in those 6 prefrontal channels that were significantly activated during the Test phase 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch30prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch30prerest(i,:) = r2z{i}(30,:);% chan30
end

clear i r2z
load('learn_z_hbo.mat')
ch30learn = zeros(N_sub,N_region);
for i = 1:N_sub
    ch30learn(i,:) = r2z{i}(30,:);% chan30
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch30learn-ch30prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[1 4 6 7 8 9 11 12 14 15 16 17 19 20 21]);% number is correct


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(1,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

% % FC30－9 and activation in Ch30
% figure('color','w');  
% scatter(dif_FC(:,9),dif_sigchan_incon_con_remain(1,:)',100,'b','.')
% hold on
% lsline
% % caxis([-1 1]);%title('Static FC')
% set(gcf,'units','centimeters');
% set(gcf,'position',[10 10 5 5])
% % xlabel('channel'); ylabel('channel');
% % set(gca,'yTick',[6 16 26 36 46]);
% % set(gca,'yTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'yTickLabel',[]);
% % set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xTick',[6 16 26 36 46]);
% % set(gca,'xTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'xTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xaxislocation','top')
% % colorbar('FontSize',8)
% axis square


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
% save tmp1 edge_corr1 % tmp1 is matrix with r value
tt = zeros(1,1035);
tt(C) = allt(C);
edge_corrT = zeros(N_region);
edge_corrT(triu(true(N_region),1)) = tt;
% save tt1 edge_corrT 
% tt1 is matrix with t value from paired t-tests
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch34prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch34prerest(i,:) = r2z{i}(34,:);% chan34
end

clear i r2z
load('learn_z_hbo.mat')
ch34learn = zeros(N_sub,N_region);
for i = 1:N_sub
    ch34learn(i,:) = r2z{i}(34,:);% chan34
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch34learn-ch34prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[1 4 6 7 8 9 11 12 14 15 16 17 19 20 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(2,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

% Fig.4B: FC34－9 and activation in Ch34
figure('color','w');  
scatter(dif_FC(:,9),dif_sigchan_incon_con_remain(2,:)',100,'b','.')% chan34
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
% save tmp2 edge_corr1
tt = zeros(1,1035);
tt(C) = allt(C);
edge_corrT = zeros(N_region);
edge_corrT(triu(true(N_region),1)) = tt;
% save tt2 edge_corrT
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch35prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch35prerest(i,:) = r2z{i}(35,:);% chan35
end

clear i r2z
load('learn_z_hbo.mat')
ch35learn = zeros(N_sub,N_region);
for i = 1:N_sub
    ch35learn(i,:) = r2z{i}(35,:);% chan35
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch35learn-ch35prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[1 4 6 7 8 9 11 12 14 15 16 17 19 20 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(3,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

% % FC35－9 and activation in Ch35
% figure('color','w');  
% scatter(dif_FC(:,9),dif_sigchan_incon_con_remain(3,:)',100,'b','.')%Ch35
% hold on
% lsline
% % caxis([-1 1]);%title('Static FC')
% set(gcf,'units','centimeters');
% set(gcf,'position',[10 10 5 5])
% % xlabel('channel'); ylabel('channel');
% % set(gca,'yTick',[6 16 26 36 46]);
% % set(gca,'yTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'yTickLabel',[]);
% % set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xTick',[6 16 26 36 46]);
% % set(gca,'xTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'xTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xaxislocation','top')
% % colorbar('FontSize',8)
% axis square

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
% save tmp3 edge_corr1
tt = zeros(1,1035);
tt(C) = allt(C);
edge_corrT = zeros(N_region);
edge_corrT(triu(true(N_region),1)) = tt;
% save tt3 edge_corrT
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch38prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch38prerest(i,:) = r2z{i}(38,:);% chan38
end

clear i r2z
load('learn_z_hbo.mat')
ch38learn = zeros(N_sub,N_region);
for i = 1:N_sub
    ch38learn(i,:) = r2z{i}(38,:);% chan38
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch38learn-ch38prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[1 4 6 7 8 9 11 12 14 15 16 17 19 20 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(4,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

% % FC38－17 and activation in Ch38
% figure('color','w');  
% scatter(dif_FC(:,17),dif_sigchan_incon_con_remain(4,:)',100,'b','.')%Ch38
% hold on
% lsline
% % caxis([-1 1]);%title('Static FC')
% set(gcf,'units','centimeters');
% set(gcf,'position',[10 10 5 5])
% % xlabel('channel'); ylabel('channel');
% % set(gca,'yTick',[6 16 26 36 46]);
% % set(gca,'yTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'yTickLabel',[]);
% % set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xTick',[6 16 26 36 46]);
% % set(gca,'xTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'xTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xaxislocation','top')
% % colorbar('FontSize',8)
% axis square

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
% save tmp4 edge_corr1
tt = zeros(1,1035);
tt(C) = allt(C);
edge_corrT = zeros(N_region);
edge_corrT(triu(true(N_region),1)) = tt;
% save tt4 edge_corrT
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch39prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch39prerest(i,:) = r2z{i}(39,:);% chan39
end

clear i r2z
load('learn_z_hbo.mat')
ch39learn = zeros(N_sub,N_region);
for i = 1:N_sub
    ch39learn(i,:) = r2z{i}(39,:);% chan39
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch39learn-ch39prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[1 4 6 7 8 9 11 12 14 15 16 17 19 20 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(5,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

% % FC39－40 and activation in Ch39
% figure('color','w');  
% scatter(dif_FC(:,40),dif_sigchan_incon_con_remain(5,:)',100,'b','.')
% hold on
% lsline
% % caxis([-1 1]);%title('Static FC')
% set(gcf,'units','centimeters');
% set(gcf,'position',[10 10 5 5])
% % xlabel('channel'); ylabel('channel');
% % set(gca,'yTick',[6 16 26 36 46]);
% % set(gca,'yTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'yTickLabel',[]);
% % set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xTick',[6 16 26 36 46]);
% % set(gca,'xTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'xTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xaxislocation','top')
% % colorbar('FontSize',8)
% axis square

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
% save tmp5 edge_corr1
tt = zeros(1,1035);
tt(C) = allt(C);
edge_corrT = zeros(N_region);
edge_corrT(triu(true(N_region),1)) = tt;
% save tt5 edge_corrT
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = find(allh==1);
clear i j  r2z


load('prerest_z_hbo.mat')
ch43prerest = zeros(N_sub,N_region);
for i = 1:N_sub
    ch43prerest(i,:) = r2z{i}(43,:);% chan43
end

clear i r2z
load('learn_z_hbo.mat')
ch43learn = zeros(N_sub,N_region);
for i = 1:N_sub
    ch43learn(i,:) = r2z{i}(43,:);% chan43
end

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
dif_FC = ch43learn-ch43prerest;
load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[1 4 6 7 8 9 11 12 14 15 16 17 19 20 21]);


for j =1:N_region
    [rvalue(j), pvalue(j)] = corr(dif_sigchan_incon_con_remain(6,:)',dif_FC(:,j));% 1-6:channel 30,34,35,38,39,43
end

% % FC43－22 and activation in Ch43
% figure('color','w');  
% scatter(dif_FC(:,22),dif_sigchan_incon_con_remain(6,:)',100,'b','.')% dif_FC(:,2);dif_FC(:,17);dif_FC(:,19);dif_FC(:,22)
% hold on
% lsline
% % caxis([-1 1]);%title('Static FC')
% set(gcf,'units','centimeters');
% set(gcf,'position',[10 10 5 5])
% % xlabel('channel'); ylabel('channel');
% % set(gca,'yTick',[6 16 26 36 46]);
% % set(gca,'yTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'yTickLabel',[]);
% % set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xTick',[6 16 26 36 46]);
% % set(gca,'xTickLabel',{'6','16', '26', '36', '46'});
% % set(gca,'xTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xaxislocation','top')
% % colorbar('FontSize',8)
% axis square

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
% save tmp6 edge_corr1
tt = zeros(1,1035);
tt(C) = allt(C);
edge_corrT = zeros(N_region);
edge_corrT(triu(true(N_region),1)) = tt;
% save tt6 edge_corrT
%% to output all seed-based FC related with all 6 activated channels in one matrix.
clear
clc
% For r value to indicate the thickness of line
load('tmp1.mat')
a1 = edge_corr1;
load('tmp2.mat')
a2 = edge_corr1;
load('tmp3.mat')
a3 = edge_corr1;
load('tmp4.mat')
a4 = edge_corr1;
load('tmp5.mat')
a5 = edge_corr1;
load('tmp6.mat')
a6 = edge_corr1;

a=a1+a2+a3+a4+a5+a6;% copy a to a txt file named "Edge46AGneonateCorrelationActivation_Learn_Pre_forAll6channels.edge"

% For t value to indicate the thickness of line
clear
clc

load('tt1.mat')
a1 = edge_corrT;
load('tt2.mat')
a2 = edge_corrT;
load('tt3.mat')
a3 = edge_corrT;
load('tt4.mat')
a4 = edge_corrT;
load('tt5.mat')
a5 = edge_corrT;
load('tt6.mat')
a6 = edge_corrT;

a=a1+a2+a3+a4+a5+a6;

%% Fig4B FC among other non-seed channels (Ch2 9 17 19 22 40)
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
    allh(j) = h;%
    allp(j) = p;
    allt(j) = stats.tstat;
end
%%
dif_FC = learn-prerest;
[x1,x2] = find(allh==1);

% sigFC_learnpre = dif_FC(x2,:);

load('dif_sigchan_incon_con.mat')
%
dif_sigchan_incon_con_remain = dif_sigchan_incon_con(:,[1 4 6 7 8 9 11 12 14 15 16 17 19 20 21]);

clear i
for i =1:length(dif_FC)
%     [rvalue(i), pvalue(i)] = corr(dif_sigchan_incon_con_remain(1,:)',sigFC_learnpre(i,:)');
    [rvalue(i), pvalue(i)] = corr(dif_sigchan_incon_con_remain(2,:)',dif_FC(i,:)');% 1-6:channel 30,34,35,38,39,43. To change it
end


[x3,x4] = find(pvalue<0.05);
 
 
C = intersect(x2,x4);
rvalue_sig = rvalue(C);

hh = zeros(1,1035);
hh(C) = rvalue_sig;
edge_corr = zeros(N_region);
edge_corr(triu(true(N_region),1)) = hh;% edge_corr represents the correlation between activation and FC of Learn-PreRest

%% Fig.4B: Correlation between activation in Ch 34 and Ch 2-9 FC strength
figure('color','w');  
scatter(dif_FC(30,:),dif_sigchan_incon_con_remain(2,:)',100,'b','.')% 1-6:channel 30,34,35,38,39,43. To change it
hold on
lsline
% caxis([-1 1]);%title('Static FC')
set(gcf,'units','centimeters');
set(gcf,'position',[10 10 6 6])
axis([-1 2 -0.04 0.08])%[xmin xmax ymin ymax]
% xlabel('channel'); ylabel('channel');
set(gca,'yTick',[-0.04 0 0.04 0.08]);
set(gca,'yTickLabel',{'-0.04','0', '0.04', '0.08'});
% set(gca,'yTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
set(gca,'xTick',[-1 0 1 2]);
set(gca,'xTickLabel',{'-1','0', '1', '2'});
% set(gca,'xTickLabel',[]);
set(gca,'FontSize',8, 'FontName','Arial')
% set(gca,'xaxislocation','top')
% colorbar('FontSize',8)
axis square

%% Fig.4E: Correlation between Pre-Rest FCs and learning-related FCs
%　Data preparation
clear
load('prerest_z_hbo.mat')

N_sub = size(r2z,1);
N_region = size(r2z{1},2); 

prerest = zeros(10,N_sub);% total 10 FCs
for i = 1:N_sub
    prerest(1,i) = r2z{i}(2,9);%
    prerest(2,i) = r2z{i}(9,30);%
    prerest(3,i) = r2z{i}(9,34);%
    prerest(4,i) = r2z{i}(9,35);%
    prerest(5,i) = r2z{i}(17,38);%
    prerest(6,i) = r2z{i}(39,40);%
    prerest(7,i) = r2z{i}(2,43);%
    prerest(8,i) = r2z{i}(17,43);%
    prerest(9,i) = r2z{i}(19,43);%
    prerest(10,i) = r2z{i}(22,43);%   
end

clear i r2z
load('learn_z_hbo.mat')
learn = zeros(10,N_sub);
for i = 1:N_sub
    learn(1,i) = r2z{i}(2,9);%
    learn(2,i) = r2z{i}(9,30);%
    learn(3,i) = r2z{i}(9,34);%
    learn(4,i) = r2z{i}(9,35);%
    learn(5,i) = r2z{i}(17,38);%
    learn(6,i) = r2z{i}(39,40);%
    learn(7,i) = r2z{i}(2,43);%
    learn(8,i) = r2z{i}(17,43);%
    learn(9,i) = r2z{i}(19,43);%
    learn(10,i) = r2z{i}(22,43);%   
end

mat = zeros(N_sub,20);

mat(:,1) = prerest(1,:)';
mat(:,2) = learn(1,:)';

mat(:,3) = prerest(2,:)';
mat(:,4) = learn(2,:)';

mat(:,5) = prerest(3,:)';
mat(:,6) = learn(3,:)';

mat(:,7) = prerest(4,:)';
mat(:,8) = learn(4,:)';

mat(:,9) = prerest(5,:)';
mat(:,10) = learn(5,:)';

mat(:,11) = prerest(6,:)';
mat(:,12) = learn(6,:)';

mat(:,13) = prerest(7,:)';
mat(:,14) = learn(7,:)';

mat(:,15) = prerest(8,:)';
mat(:,16) = learn(8,:)';

mat(:,17) = prerest(9,:)';
mat(:,18) = learn(9,:)';

mat(:,19) = prerest(10,:)';
mat(:,20) = learn(10,:)';

csvwrite('FC_PreRest_Learn.csv', mat);% This csv file was used to produce S5 Fig.


FC_pre = mat(:,[1:2:20]);
FC_learn = mat(:,[2:2:20]);

%
FC_learn_pre = FC_learn-FC_pre;
clear i rvalue pvalue
for i =1:10 % 10 tagerted FC
    [rvalue(i), pvalue(i)] = corr(FC_pre(:,i),FC_learn_pre(:,i));
end

% Fig.4E: Ch2-9
for j = 1
figure('color','w');  
scatter(FC_pre(:,j),FC_learn_pre(:,j),100,'b','.')
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
end

% % Ch9-30
% for j = 2
% figure('color','w');  
% scatter(FC_pre(:,j),FC_learn_pre(:,j),100,'b','.')
% hold on
% lsline
% axis([-1 2 -1 2])%[xmin xmax ymin ymax]
% set(gcf,'units','centimeters');
% set(gcf,'position',[10 10 5 5])
% % xlabel('channel'); ylabel('channel');
% set(gca,'yTick',[-1 0 1 2]);
% set(gca,'yTickLabel',{'-1', '0', '1', '2'});
% % set(gca,'yTickLabel',[]);
% % set(gca,'FontSize',8, 'FontName','Arial')
% set(gca,'xTick',[-1 0 1 2]);
% set(gca,'xTickLabel',{'-1', '0', '1', '2'});
% % set(gca,'xTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xaxislocation','top')
% % colorbar('FontSize',8)
% axis square
% end

% Fig.4E: Ch9-34
for j = 3
figure('color','w');  
scatter(FC_pre(:,j),FC_learn_pre(:,j),100,'b','.')
hold on
lsline
axis([-1 2 -1 2])%[xmin xmax ymin ymax]
set(gcf,'units','centimeters');
set(gcf,'position',[10 10 5 5])
% xlabel('channel'); ylabel('channel');
set(gca,'yTick',[-1 0 1 2]);
set(gca,'yTickLabel',{'-1', '0', '1', '2'});
% set(gca,'yTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
set(gca,'xTick',[-1 0 1 2]);
set(gca,'xTickLabel',{'-1', '0', '1', '2'});
% set(gca,'xTickLabel',[]);
set(gca,'FontSize',8, 'FontName','Arial')
% set(gca,'xaxislocation','top')
% colorbar('FontSize',8)
axis square
end

% % Ch9-35
% for j = 4
% figure('color','w');  
% scatter(FC_pre(:,j),FC_learn_pre(:,j),100,'b','.')
% hold on
% lsline
% axis([-1 2 -1 2])%[xmin xmax ymin ymax]
% set(gcf,'units','centimeters');
% set(gcf,'position',[10 10 5 5])
% % xlabel('channel'); ylabel('channel');
% set(gca,'yTick',[-1 0 1 2]);
% set(gca,'yTickLabel',{'-1', '0', '1', '2'});
% % set(gca,'yTickLabel',[]);
% % set(gca,'FontSize',8, 'FontName','Arial')
% set(gca,'xTick',[-1 0 1 2]);
% set(gca,'xTickLabel',{'-1', '0', '1', '2'});
% % set(gca,'xTickLabel',[]);
% set(gca,'FontSize',8, 'FontName','Arial')
% % set(gca,'xaxislocation','top')
% % colorbar('FontSize',8)
% axis square
% end