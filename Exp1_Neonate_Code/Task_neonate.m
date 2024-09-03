%% Statistical analyses for HbO (main findings)
%% to compare correct condition and baseline
% define representative value
clear all
load incongruent_allsub_2probe_final.mat
load congruent_allsub_2probe_final.mat

% for congruent condition
N_sub = length(congruent_allsub_2probe_final);
N_chan = size(congruent_allsub_2probe_final{1},1);
reprevalue_allchan_con = zeros(N_chan,1);
reprevalue_allchan_con_baseline = zeros(N_chan,1);
reprevalue_allsub_con = zeros(N_chan,N_sub);
reprevalue_allsub_con_baseline = zeros(N_chan,N_sub);
for i = 1:N_sub
    for j = 1:N_chan
        reprevalue_allchan_con(j,:) = mean(congruent_allsub_2probe_final{i}(j,125:175),2);% early window
        reprevalue_allchan_con_baseline(j,:) = mean(congruent_allsub_2probe_final{i}(j,1:50),2);
    end
    reprevalue_allsub_con(:,i) = reprevalue_allchan_con;
    reprevalue_allsub_con_baseline(:,i) = reprevalue_allchan_con_baseline;
end

% for incongruent condition
reprevalue_allchan_incon = zeros(N_chan,1);
reprevalue_allchan_incon_baseline = zeros(N_chan,1);
reprevalue_allsub_incon = zeros(N_chan,N_sub);
reprevalue_allsub_incon_baseline = zeros(N_chan,N_sub);

clear i j
for i = 1:N_sub
    for j = 1:N_chan
        reprevalue_allchan_incon(j,:) = mean(incongruent_allsub_2probe_final{i}(j,125:175),2);% early window
        reprevalue_allchan_incon_baseline(j,:) = mean(incongruent_allsub_2probe_final{i}(j,1:50),2);
    end
    reprevalue_allsub_incon(:,i) = reprevalue_allchan_incon;
    reprevalue_allsub_incon_baseline(:,i) = reprevalue_allchan_incon_baseline;
end

%% Permutation test to compare correct condition and average value in baseline period
% if we conduct permutation test channel by channel, after removing rows with nan values
dif_allchannel = zeros(N_sub,N_chan); 
pval_all = zeros(N_chan,1);
t_orig_all = zeros(N_chan,1);
crit_t_all = cell(N_chan,1);
clear k array
for k = 1:N_chan
    array = [reprevalue_allsub_con(k,:)' reprevalue_allsub_con_baseline(k,:)'];
    array(any(isnan(array),2),:) = [];% to delete the rows in which either column 1 or 2 is NaN
    dif = array(:,1)-array(:,2);
    [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(dif,50000);
    pval_all(k,1) = pval;
    t_orig_all(k,1) = t_orig;
    crit_t_all{k,1} = crit_t;
end
index = find(pval_all<0.05); % 18 channels are significant
table_p = pval_all(index);
table_t = t_orig_all(index);
% numbers of S4 Table are from t_orig_all, pval_all 

%% Permutation test to compare incorrect condition and average value in baseline period
% if we conduct permutation test channel by channel, after removing rows with nan values
dif_allchannel = zeros(N_sub,N_chan); 
pval_all = zeros(N_chan,1);
t_orig_all = zeros(N_chan,1);
crit_t_all = cell(N_chan,1);
clear k array
for k = 1:N_chan
    array = [reprevalue_allsub_incon(k,:)' reprevalue_allsub_incon_baseline(k,:)'];
    array(any(isnan(array),2),:) = [];% to delete the rows in which either column 1 or 2 is NaN
    dif = array(:,1)-array(:,2);
    [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(dif,50000);
    pval_all(k,1) = pval;
    t_orig_all(k,1) = t_orig;
    crit_t_all{k,1} = crit_t;
end
index = find(pval_all<0.05);
table_p = pval_all(index); % no significant results
table_t = t_orig_all(index);

%% Permutation test to compare correct and incorrect conditons
% if we conduct permutation test channel by channel, after remove rows with nan values
dif_allchannel = zeros(N_sub,N_chan); 
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
index = find(pval_all<0.05);
table_p = pval_all(index); % Ch 30, 34, 35, 38, 39, 43 are significant
table_t = t_orig_all(index);
% numbers of S4 Table are from t_orig_all, pval_all 


%% Statistical analyses for HbR
%% to compare correct condition and baseline
% define representative value
clear all
load incongruent_allsub_2probe_final_hbr.mat
load congruent_allsub_2probe_final_hbr.mat

% for congruent condition
N_sub = length(congruent_allsub_2probe_final_hbr);
N_chan = size(congruent_allsub_2probe_final_hbr{1},1);
reprevalue_allchan_con = zeros(N_chan,1);
reprevalue_allchan_con_baseline = zeros(N_chan,1);
reprevalue_allsub_con = zeros(N_chan,N_sub);
reprevalue_allsub_con_baseline = zeros(N_chan,N_sub);
for i = 1:N_sub
    for j = 1:N_chan
        reprevalue_allchan_con(j,:) = mean(congruent_allsub_2probe_final_hbr{i}(j,125:175),2);% early window
        reprevalue_allchan_con_baseline(j,:) = mean(congruent_allsub_2probe_final_hbr{i}(j,1:50),2);
    end
    reprevalue_allsub_con(:,i) = reprevalue_allchan_con;
    reprevalue_allsub_con_baseline(:,i) = reprevalue_allchan_con_baseline;
end

% for incongruent condition
reprevalue_allchan_incon = zeros(N_chan,1);
reprevalue_allchan_incon_baseline = zeros(N_chan,1);
reprevalue_allsub_incon = zeros(N_chan,N_sub);
reprevalue_allsub_incon_baseline = zeros(N_chan,N_sub);

clear i j
for i = 1:N_sub
    for j = 1:N_chan
        reprevalue_allchan_incon(j,:) = mean(incongruent_allsub_2probe_final_hbr{i}(j,125:175),2);% early window
        reprevalue_allchan_incon_baseline(j,:) = mean(incongruent_allsub_2probe_final_hbr{i}(j,1:50),2);
    end
    reprevalue_allsub_incon(:,i) = reprevalue_allchan_incon;
    reprevalue_allsub_incon_baseline(:,i) = reprevalue_allchan_incon_baseline;
end

%% Permutation test to compare correct condition and average value in baseline period
% if we conduct permutation test channel by channel, after removing rows with nan values
dif_allchannel = zeros(N_sub,N_chan); 
pval_all = zeros(N_chan,1);
t_orig_all = zeros(N_chan,1);
crit_t_all = cell(N_chan,1);
clear k array
for k = 1:N_chan
    array = [reprevalue_allsub_con(k,:)' reprevalue_allsub_con_baseline(k,:)'];
    array(any(isnan(array),2),:) = [];% to delete the rows in which either column 1 or 2 is NaN
    dif = array(:,1)-array(:,2);
    [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(dif,50000);
    pval_all(k,1) = pval;
    t_orig_all(k,1) = t_orig;
    crit_t_all{k,1} = crit_t;
end
index = find(pval_all<0.05);
table_p = pval_all(index);
table_t = t_orig_all(index);

%% Permutation test to compare incorrect condition and average value in baseline period
% if we conduct permutation test channel by channel, after removing rows with nan values
dif_allchannel = zeros(N_sub,N_chan); 
pval_all = zeros(N_chan,1);
t_orig_all = zeros(N_chan,1);
crit_t_all = cell(N_chan,1);
clear k array
for k = 1:N_chan
    array = [reprevalue_allsub_incon(k,:)' reprevalue_allsub_incon_baseline(k,:)'];
    array(any(isnan(array),2),:) = [];% to delete the rows in which either column 1 or 2 is NaN
    dif = array(:,1)-array(:,2);
    [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t1(dif,50000);
    pval_all(k,1) = pval;
    t_orig_all(k,1) = t_orig;
    crit_t_all{k,1} = crit_t;
end
index = find(pval_all<0.05);
table_p = pval_all(index);
table_t = t_orig_all(index);

%% Permutation test to compare correct and incorrect conditons
% if we conduct permutation test channel by channel, after removing rows with nan values
dif_allchannel = zeros(N_sub,N_chan); 
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
index = find(pval_all<0.05);
table_p = pval_all(index);
table_t = t_orig_all(index);




















