function[]=Plot_FRET_tau_DandA(Data,burn_in_perc,bin_size,mmk)


% Select the burn-in part of the data
burn_in      = max([1 floor((burn_in_perc/100)*length(Data.F_all(1,:)))]);


% Set the figure type and dimensions
figure('Name','Continuous lifetime estimation','Units','inch','Position',[0 0 8 7])
set(gcf,'color','w')
clf
movegui(gcf,'center')

sigg_n1 = [];
sigg_n2 = [];
sigg_n3 = [];
max_1(1)=0;
for mme=1:Data.num_sub_sigs
    sigg_n1 =[sigg_n1,Data.signals{mme}(1,:)];
    if  mme==1
        sigg_n2 =Data.signals{mme}(2,:);
        max_1(2) = max(unique(sigg_n2));
    else
        sigg_n2 =[sigg_n2, max_1(mme)+Data.signals{mme}(2,:)];
        max_1(mme+1) = max(sigg_n2);
    end
    sigg_n3 =[sigg_n3,Data.signals{mme}(3,:)];
end

dbnn=linspace(min(sigg_n1),max(sigg_n1),40);



sig1=sigg_n2(sigg_n3==1);
sig2=sigg_n2(sigg_n3==2);

edge = 0:bin_size:max(sigg_n2)+2*bin_size;
sig_bin1 = histcounts(sig1,edge);
sig_bin2 = histcounts(sig2,edge);





h13r=subplot(7,4,[1,2,3]);
% h13r.Position = h13r.Position+[0 0.025 0 0];
plot(sig1,sigg_n1(sigg_n3==1),'.','color',[0 .7 0],'MarkerSize',1)
ylabel({'Donor','Delay times','(ns)'})
set(gca,'XTick',[])
box off
xlim([0 max(sigg_n2)])
ylim([0,max(sigg_n1)])


h137=subplot(7,4,[5,6,7]);
% h137.Position = h137.Position+[0 0.025 0 0];
plot(sig2,sigg_n1(sigg_n3==2),'.','color','r','MarkerSize',1)
ylabel({'Acceptor','Delay times','(ns)'})
set(gca,'XTick',[])
box off
xlim([0 max(sigg_n2)])
ylim([0,max(sigg_n1)])



h131=subplot(7,4,[9,10,11]);
% h131.Position = h131.Position+[0 0.04 0 0];
stairs(edge(1:end-1),sig_bin1,'color',[0 .7 0])
hold on
stairs(edge(1:end-1),sig_bin2,'color','r')
stairs(edge(1:end-1),sig_bin1+sig_bin2,'color','k')
hold off
box off
xlim([0 max(sigg_n2)])
% xlabel('Time (s)')
set(gca,'XTick',[])
ylim([0 max([max(sig_bin1+sig_bin2)])*1.1])
ylabel({'photons/',[num2str(bin_size*1000),' ms']})






subplot(7,4,[13,14,15]) % Pulse macro times
if isfield(Data,'KD_real')
   stairs(Data.signal(2,:),Data.F_fluc_real,'g')
end

ress=[];
F_Ind_points=[];
map_F=[];
hold on
for mmk=1:Data.num_sub_sigs
    F_Ind_points    = [F_Ind_points,max_1(mmk)+Data.signals{mmk}(2,:)];
    ttf=[];
for ii=1:size(Data.F_all(1,burn_in:end),2)
    ttf(:,ii)=Data.KF_val_learned{mmk}(burn_in+ii-1,Data.F_res{mmk}(burn_in+ii-1,:));
end
    ress            = [ress,ttf']; 
    map_F     = [map_F,Data.map_KF_val(1,Data.map_F_res{mmk})];
end
    F_mean_down1 = quantile(ress , 0.025);
    F_mean_up1   = quantile(ress , 0.975);
    F_mean1      = quantile(ress , 0.5);

    stairs(F_Ind_points,F_mean1,'b')
    xx = [F_Ind_points , fliplr(F_Ind_points)];
    fill(xx,[F_mean_down1 , fliplr(F_mean_up1)], 'c','facealpha', 0.2,'edgecolor', 'none')
   
    stairs(F_Ind_points,map_F,'-.','Color','m')

    set(gca,'XTick',[])
    if isfield(Data,'KD_real')
    ylim([0 max([max(F_mean_up1), max(Data.F_fluc_real),max(map_F) ])*1.2])
    else
    ylim([0 max([max(F_mean_up1),max(F_Ind_points,map_F) ])*1.2])    
    end
    ylabel({'FRET rate','(ns^{-1})'})
    box off
    xlim([0 max(sigg_n2)])



subplot(7,4,16)
dbnd=linspace(0,max(F_mean1)*1.2,50);
maxxy=0;

    for mm=1:size(Data.KF_val_learned{mmk},2)
        hh=histogram(Data.KF_val_learned{mmk}(burn_in:end,mm),dbnd,'Normalization','pdf','Orientation','Horizontal','FaceColor',lines(1));
        hold on
        loo=quantile(Data.KF_val_learned{mmk}(burn_in:end,mm),0.025);
        hii=quantile(Data.KF_val_learned{mmk}(burn_in:end,mm),0.975);
        maxxy=max([maxxy,1.2*max(unique(hh.Values))]);
        h1e=fill([0 1 1 0].*maxxy,[loo loo hii hii], 'c','facealpha', 0.2,'edgecolor', 'none');
        
        histogram(Data.KF_val_learned{mmk}(burn_in:end,mm),dbnd,'Normalization','pdf','Orientation','Horizontal','FaceColor',lines(1));
    end


for mm=1:size(Data.KF_val_learned{mmk},2)
    hee3=line(maxxy.*[0 1],[1 1]*Data.map_KF_val(mm),'LineStyle','-.','color','m','LineWidth',1);
end

if isfield(Data,'KD_real')
    for kll=1:length(Data.KF_val_real)
        h3e=line(maxxy.*[0 1],[1 1]*Data.KF_val_real(kll),'LineStyle','-.','color','g','LineWidth',1);
    end
end
if isfield(Data,'KD_real')
    ylim([0 max([max(F_mean_up1), max(Data.F_fluc_real),max(map_F) ])*1.2])
    else
    ylim([0 max([max(F_mean_up1),max(F_Ind_points,map_F) ])*1.2])    
end
% xlabel('FRET rate')
% ylabel('Pdf')
set(gca,'YTick',[])
box off
if  isfield(Data,'KD_val_real')
    h23=legend([h3e,hh,h1e,hee3],'Ground truth','Postr. prob. distr.','95% confidence interval','MAP','EdgeColor',[1 1 1],'location','Northoutside');
    h23.Position = h23.Position+[0.06 0.1 0 0];
else
    h23=legend([hh,h1e,hee3],'Postr. prob. distr.','95% confidence interval','MAP','EdgeColor',[1 1 1],'location','Northoutside');
    h23.Position = h23.Position+[0.06 0.1 0 0];
end


subplot(7,4,[17,18,19])

if isfield(Data,'KD_real')
   stairs(Data.signal(2,:),Data.tau_D_fluc_real,'g')
end

ress=[];
tau_D_Ind_points=[];
map_p =[];

hold on

for mmk=1:Data.num_sub_sigs
    
    
    tau_D_Ind_points=[tau_D_Ind_points,max_1(mmk)+Data.signals{mmk}(2,:)];
    tt=[];
    for ii=1:size(Data.F_all(1,burn_in:end),2)
        tt(:,ii)=Data.KD_val_learned{mmk}(burn_in+ii-1,Data.tau_D_res{mmk}(burn_in+ii-1,:));
    end
    ress=[ress,tt'];
    map_p = [map_p,Data.map_KD_val(1,Data.map_tau_D{mmk})];
end
    tau_D_mean_down1 = quantile(ress , 0.025);
    tau_D_mean_up1   = quantile(ress , 0.975);
    tau_D_mean1      = quantile(ress , 0.5);

    stairs(tau_D_Ind_points,tau_D_mean1,'b')
    xx = [tau_D_Ind_points , fliplr(tau_D_Ind_points)];
    fill(xx,[tau_D_mean_down1 , fliplr(tau_D_mean_up1)], 'c','facealpha', 0.2,'edgecolor', 'none')
%     set(gca,'XTick',[])
    stairs(tau_D_Ind_points,map_p,'-.','Color','m')
    
    if  isfield(Data,'KD_real')
        ylim([0 max([max(tau_D_mean_up1), max(Data.tau_D_fluc_real),max(map_p) ])*1.2])
    else
        ylim([0 max([max(tau_D_mean_up1),max(tau_D_Ind_points,map_p) ])*1.2])
    end

ylabel({'\tau_D (ns)'})
box off
set(gca,'XTick',[])

xlim([0 max(sigg_n2)])




subplot(7,4,20)
dbnd=linspace(0,max(tau_D_mean_up1)*1.2,50);
maxxy=0;
    for mm=1:size(Data.KD_val_learned{mmk},2)
        hh=histogram(Data.KD_val_learned{mmk}(burn_in:end,mm),dbnd,'Normalization','pdf','Orientation','Horizontal','FaceColor',lines(1));
        hold on
        loo=quantile(Data.KD_val_learned{mmk}(burn_in:end,mm),0.025);
        hii=quantile(Data.KD_val_learned{mmk}(burn_in:end,mm),0.975);
        maxxy=max([maxxy,1.2*max(unique(hh.Values))]);
        fill([0 1 1 0].*maxxy,[loo loo hii hii], 'c','facealpha', 0.2,'edgecolor', 'none')
        
        histogram(Data.KD_val_learned{mmk}(burn_in:end,mm),dbnd,'Normalization','pdf','Orientation','Horizontal','FaceColor',lines(1));
    end


for mm=1:size(Data.KD_val_learned{mmk},2)
    line(maxxy.*[0 1],[1 1]*Data.map_KD_val(mm),'LineStyle','-.','color','m','LineWidth',1)
end


if isfield(Data,'KD_real')
    for kll=1:length(Data.KD_val_real)
        line(maxxy.*[0 1],[1 1]*Data.KD_val_real(kll),'LineStyle','-.','color','g','LineWidth',1)
    end
end
if  isfield(Data,'KD_real')
        ylim([0 max([max(tau_D_mean_up1), max(Data.tau_D_fluc_real),max(map_p) ])*1.2])
    else
        ylim([0 max([max(tau_D_mean_up1),max(tau_D_Ind_points,map_p) ])*1.2])
    end
% xlabel('\tau_D (ns)')
% ylabel('Pdf')
set(gca,'YTick',[])
box off





subplot(7,4,[21,22,23])

if isfield(Data,'KD_real')
   stairs(Data.signal(2,:),Data.tau_A_fluc_real,'g')
end

ress=[];
tau_A_Ind_points=[];
map_A =[];

hold on

for mmk=1:Data.num_sub_sigs
    tau_A_Ind_points=[tau_A_Ind_points,max_1(mmk)+Data.signals{mmk}(2,:)];
    tt=[];
    for ii=1:size(Data.A_all(1,burn_in:end),2)
        tt(:,ii)=Data.KA_val_learned{mmk}(burn_in+ii-1,Data.tau_A_res{mmk}(burn_in+ii-1,:));
    end
    ress=[ress,tt'];
    map_A = [map_A,Data.map_KA_val(1,Data.map_tau_A{mmk})];
end
    tau_A_mean_down1 = quantile(ress , 0.025);
    tau_A_mean_up1   = quantile(ress , 0.975);
    tau_A_mean1      = quantile(ress , 0.5);

    stairs(tau_A_Ind_points,tau_A_mean1,'b')
    xx = [tau_A_Ind_points , fliplr(tau_A_Ind_points)];
    fill(xx,[tau_A_mean_down1 , fliplr(tau_A_mean_up1)], 'c','facealpha', 0.2,'edgecolor', 'none')
    set(gca,'XTick',[])
    stairs(tau_A_Ind_points,map_A,'-.','Color','m')
    
    if  isfield(Data,'KD_real')
        ylim([0 max([max(tau_A_mean_up1), max(Data.tau_A_fluc_real),max(map_A) ])*1.2])
    else
        ylim([0 max([max(tau_A_mean_up1),max(tau_A_Ind_points,map_A) ])*1.2])
    end

ylabel({'\tau_A (ns)'})
box off

xlim([0 max(sigg_n2)])



subplot(7,4,24)
dbnd=linspace(0,max(tau_A_mean_up1)*1.2,50);
maxxy=0;
    for mm=1:size(Data.KA_val_learned{mmk},2)
        hh=histogram(Data.KA_val_learned{mmk}(burn_in:end,mm),dbnd,'Normalization','pdf','Orientation','Horizontal','FaceColor',lines(1));
        hold on
        loo=quantile(Data.KA_val_learned{mmk}(burn_in:end,mm),0.025);
        hii=quantile(Data.KA_val_learned{mmk}(burn_in:end,mm),0.975);
        maxxy=max([maxxy,1.2*max(unique(hh.Values))]);
        fill([0 1 1 0].*maxxy,[loo loo hii hii], 'c','facealpha', 0.2,'edgecolor', 'none')
        
        histogram(Data.KA_val_learned{mmk}(burn_in:end,mm),dbnd,'Normalization','pdf','Orientation','Horizontal','FaceColor',lines(1));
    end


for mm=1:size(Data.KA_val_learned{mmk},2)
    line(maxxy.*[0 1],[1 1]*Data.map_KA_val(mm),'LineStyle','-.','color','m','LineWidth',1)
end


if isfield(Data,'KD_real')
    for kll=1:length(Data.KA_val_real)
        line(maxxy.*[0 1],[1 1]*Data.KA_val_real(kll),'LineStyle','-.','color','g','LineWidth',1)
    end
end
if  isfield(Data,'KD_real')
        ylim([0 max([max(tau_A_mean_up1), max(Data.tau_A_fluc_real),max(map_A) ])*1.2])
    else
        ylim([0 max([max(tau_A_mean_up1),max(tau_A_Ind_points,map_A) ])*1.2])
    end
xlabel('Pdf')
set(gca,'YTick',[])
box off



subplot(7,4,[25 26 27])
tt=0;
for mm=1:length(Data.signals)
if  isfield(Data,'KD_real')
stairs(Data.signals{mm}(2,:),Data.F_fluc_real./(Data.F_fluc_real+(1./min(Data.KD_val_real))) , 'g')
end
hold on
stairs(Data.signals{mm}(2,:)+tt,Data.map_KF_val(Data.map_F_res{mm})'./(Data.map_KF_val(Data.map_F_res{mm})'+(1./min(Data.map_KD_val))) , 'b')
tt=max(Data.signals{mm}(2,:))+tt;
end
ylabel('FRET efficiency')
xlabel('Time (s)')
xlim([0 max(sigg_n2)])
ylim([0 1])
end
