%% S2 Fig: To plot time course after block averaging for each condition over all subjects
clear all
clc

load incongruent_allsub_2probe_final.mat
load congruent_allsub_2probe_final.mat
load incongruent_allsub_2probe_final_hbr.mat
load congruent_allsub_2probe_final_hbr.mat

D_con = cat(3, congruent_allsub_2probe_final{:});
DD_con = nanmean(D_con,3);
A_con = mean(DD_con,1);

D_con_hbr = cat(3, congruent_allsub_2probe_final_hbr{:});
DD_con_hbr = nanmean(D_con_hbr,3);
A_con_hbr = mean(DD_con_hbr,1);

D_incon = cat(3, incongruent_allsub_2probe_final{:});
DD_incon = nanmean(D_incon,3);
A_incon = mean(DD_incon,1);

D_incon_hbr = cat(3, incongruent_allsub_2probe_final_hbr{:});
DD_incon_hbr = nanmean(D_incon_hbr,3);
A_incon_hbr = mean(DD_incon_hbr,1);


% Plot
for i = 1:size(DD_con,1)
    figure('name','avgHbR','color','w');
    set(gcf,'units','centimeters');
    set(gcf,'position',[10 10 3 3]);
    
    plot(1:size(DD_con,2),DD_con(i,:),'color',[0 128/255 0], 'linewidth', 1)
    hold on
    plot(1:size(DD_incon,2),DD_incon(i,:),'color',[236/255 112/255  22/255], 'linewidth', 1)

    plot(1:size(DD_con_hbr,2),DD_con_hbr(i,:),'--','color',[0 128/255 0], 'linewidth', 1)
    plot(1:size(DD_incon_hbr,2),DD_incon_hbr(i,:),'--','color',[236/255 112/255  22/255], 'linewidth', 1)
             
%     xlabel('Time(s)'); ylabel('ΔHb(mM mm)', 'FontSize', 6);
    set(gca,'XLim',[0 300]);
    set(gca,'YLim',[-0.05 0.05]);
%   line([50 50],[-0.1,0.1]);
%   patch([50 165 165 50],[-0.1 -0.1 0.1 0.1],[0 .5 .5],'EdgeColor',[0 .5 .5],'FaceAlpha',0.3)%'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5],'facealpha',0.5);
%   axis([-50 200 -5 5])
    set(gca,'ytick',[-0.05 -0.025 0 0.025 0.05],'LineWidth',1);
    set(gca,'YTickLabel',[]);
%   set(gca,'YTickLabel',[-0.05 0 0.05],'fontsize',8, 'FontName','Arial');

    set(gca,'xTick',[0 50 100 150 200 250 300],'LineWidth',1);
    set(gca,'XTickLabel',[]);
%   set(gca,'xTickLabel',{'-5', '0', '5', '10', '15', '20', '25'},'fontsize',8, 'FontName','Arial');
%     title([num2str(i)],'fontsize',12)
    box off
    print([num2str(i)], '-dtiff', '-r600');
end

close all