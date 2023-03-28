function[] = Plot_transition_rates(Data,burn_in_perc,mmk)

figure('Name','Continuous lifetime estimation','Units','inch','Position',[0 0 7 4])
set(gcf,'color','w')
clf
movegui(gcf,'center')

burn_in      = max([1 floor((burn_in_perc/100)*length(Data.KF_val_learned{mmk}(:,1)))]);


bb0=subplot(3,4,1);
for nm=1:Data.Num_states_F
    for mm=1:Data.Num_states_F
        if nm~=mm
        plot(permute(Data.KF_rate_learned(nm,mm,:),[1,3,2]),'b')
        hold on
        plot(permute(Data.KF_rate_learned(nm,mm,1:burn_in),[1,3,2]),'r')
        line([1 size(Data.KF_rate_learned,3)] ,Data.map_KF_rate(nm,mm)*[1 1],'Color','m','LineStyle','--')
        end
    end
end
if isfield(Data,'F_fluc_real')
   for nm=1:length(Data.KF_real)
       for mm=1:length(Data.KF_real)
        if nm~=mm
           line([1 size(Data.KF_rate_learned,3)] ,Data.KF_real(nm,mm)*[1 1],'Color','g','LineStyle','--')
        end
       end
   end
end
xlabel('Iteration')
ylabel('K_F (s^{-1})')
ylim([0 max(max(max(Data.KF_rate_learned(:,:,burn_in:end))))+eps])
bb0.Position=bb0.Position.*[1 1 1.1 1];



bb1=subplot(3,4,2);
TT=linspace(0,max(max(max(Data.KF_rate_learned(:,:,burn_in:end)))),100);
maxx=0;
for nm=1:Data.Num_states_F
     for mm=1:Data.Num_states_F
        if nm~=mm
    hh1=histogram(Data.KF_rate_learned(nm,mm,burn_in:end),TT,...
        'Normalization','pdf','Orientation','Horizontal','EdgeColor','none','FaceColor','b');
    hold on
    if  max(unique(hh1.Values))>maxx
        maxx = max(unique(hh1.Values));
    end
        end
     end
end

for m=1:Data.Num_states_F
     for mm=1:Data.Num_states_F
        if nm~=mm
    line([0 1.1]*maxx ,Data.map_KF_rate(m,mm)*[1 1],'Color','m','LineStyle','--')
    hold on
        end
     end
end

 if  isfield(Data,'F_fluc_real')
     for nm=1:length(Data.KF_real)
          for mm=1:length(Data.KF_real)
        if nm~=mm
        line([0 1.1]*maxx,Data.KF_real(nm,mm)*[1 1],'Color','g','LineStyle','--')
        end
          end
     end
 end
    
plot(gampdf(TT,Data.KF_rate_alpha ,Data.KF_rate_beta),TT,'k')

set(gca,'YTick',[])
xlabel('pdf')
ylim([0 max(max(max(Data.KF_rate_learned(:,:,burn_in:end))))+eps])
bb1.Position=bb1.Position.*[1 1 0.6 1];



bb2=subplot(3,4,3);
for nm=1:Data.Num_states_F
plot(Data.KF_val_learned{mmk}(:,nm),'b')
hold on
plot(Data.KF_val_learned{mmk}(1:burn_in,nm),'r')
line([1 size(Data.KF_rate_learned,3)] ,Data.map_KF_val(nm)*[1 1],'Color','m','LineStyle','--')
end
if  isfield(Data,'F_fluc_real')
    for nm=1:length(Data.KF_val_real)
       line([1 size(Data.KF_rate_learned,3)] ,Data.KF_val_real(nm)*[1 1],'Color','g','LineStyle','--')
    end
end

xlabel('Iteration')
ylabel('\mu_F (ns^{-1})')
ylim([0 max(max(Data.KF_val_learned{mmk}(burn_in:end,:)))])
bb2.Position=bb2.Position.*[1 1 1.1 1];


bb3=subplot(3,4,4);
TT=linspace(0,max(max(Data.KF_val_learned{mmk}(burn_in:end,:))),100);
maxx=0;
for nm=1:length(Data.KF_val_learned(1,:))
hh1=histogram(Data.KF_val_learned{mmk}(burn_in:end,nm),TT,...
    'Normalization','pdf','Orientation','Horizontal','EdgeColor','none','FaceColor','b');
hold on
if maxx<max(unique(hh1.Values))
    maxx=max(unique(hh1.Values));
end
end

plot(gampdf(TT,Data.KF_val_alpha ,Data.KF_val_beta),TT,'k')
if  isfield(Data,'F_fluc_real')
    for mm=1:length(Data.KF_val_real)
    line([0 1.1]*maxx,Data.KF_val_real(mm)*[1 1],'Color','g','LineStyle','--')
    end
end

for m=1:length(Data.KF_val_learned{mmk}(1,:))
line([0 1.1]*maxx ,Data.map_KF_val(m)*[1 1],'Color','m','LineStyle','--')
end

set(gca,'YTick',[])
xlabel('pdf')
ylim([0 max(max(Data.KF_val_learned{mmk}(burn_in:end,:)))])
bb3.Position=bb3.Position.*[1 1 0.6 1];






bb4=subplot(3,4,5);
for nm=1:Data.Num_states_D
    for mm=1:Data.Num_states_D
           if mm~=nm
plot(permute(Data.KD_rate_learned(nm,mm,:),[1,3,2]),'b')
hold on
plot(permute(Data.KD_rate_learned(nm,mm,1:burn_in),[1,3,2]),'r')
line([1 size(Data.KF_rate_learned,3)] ,Data.map_KD_rate(nm,mm)*[1 1],'Color','m','LineStyle','--')
           end
    end
end

if isfield(Data,'F_fluc_real')
   for mn=1:length(Data.KD_real)
       for mm=1:length(Data.KD_real)
           if mm~=mn
       line([1 size(Data.KF_rate_learned,3)],Data.KD_real(mn,mm)*[1 1],'Color','g','LineStyle','--')
           end
       end
   end
end
xlabel('Iteration')
ylabel('K_D (s^{-1})')
if Data.Num_states_D*(Data.Num_states_D-1)>0
ylim([0 max(max(max(Data.KD_rate_learned)))+0.001])
end
bb4.Position=bb4.Position.*[1 1 1.1 1];


bb5=subplot(3,4,6);

if Data.Num_states_D*(Data.Num_states_D-1)>0
TT=linspace(0,max(max(max(Data.KD_rate_learned(:,:,burn_in:end)))),100);
end
maxx=0;
for nm=1:Data.Num_states_D
    for mm=1:Data.Num_states_D
           if mm~=nm
hh1=histogram(Data.KD_rate_learned(nm,mm,burn_in:end),TT,...
    'Normalization','pdf','Orientation','Horizontal','EdgeColor','none','FaceColor','b');
hold on
if maxx<max(unique(hh1.Values))
    maxx=max(unique(hh1.Values));
end
           end
    end
end
plot(gampdf(TT,Data.KD_rate_alpha ,Data.KD_rate_beta),TT,'k')

for nm=1:Data.Num_states_D
    for mm=1:Data.Num_states_D
        if mm~=nm
           line([0 1.1]*maxx ,Data.map_KD_rate(nm,mm)*[1 1],'Color','m','LineStyle','--')
        end
    end
end
    
if isfield(Data,'F_fluc_real')
   for mn=1:length(Data.KD_real)
       for mm=1:length(Data.KD_real)
           if mm~=mn
       line([0 1.1]*maxx,Data.KD_real(mn,mm)*[1 1],'Color','g','LineStyle','--')
           end
       end
   end
end


set(gca,'YTick',[])
xlabel('pdf')
if Data.Num_states_D>0
ylim([0 max(max(max(Data.KD_rate_learned(:,:,burn_in:end))))+0.01])
end
bb5.Position=bb5.Position.*[1 1 0.6 1];



bb6=subplot(3,4,7);
for nm=1:Data.Num_states_D
plot(Data.KD_val_learned{mmk}(:,nm),'b')
hold on
plot(Data.KD_val_learned{mmk}(1:burn_in,nm),'r')
line([1 size(Data.KF_rate_learned,3)] ,Data.map_KD_val(nm)*[1 1],'Color','m','LineStyle','--')
end
if  isfield(Data,'F_fluc_real')
    for nm=1:length(Data.KD_val_real)
    line([1 size(Data.KF_rate_learned,3)] ,Data.KD_val_real(nm)*[1 1],'Color','g','LineStyle','--')
    end
end


xlabel('Iteration')
ylabel('\tau_D (ns)')
ylim([0 max(max(Data.KD_val_learned{mmk}(burn_in:end,:)))])
bb6.Position=bb6.Position.*[1 1 1.1 1];



bb7=subplot(3,4,8);
maxx=0;
TT=linspace(0,max(max(Data.KD_val_learned{mmk})),100);
for nm=1:Data.Num_states_D 
hh1=histogram(Data.KD_val_learned{mmk}(burn_in:end,nm),TT,...
    'Normalization','pdf','Orientation','Horizontal','EdgeColor','none','FaceColor','b');
hold on
if maxx<max(unique(hh1.Values))
   maxx=max(unique(hh1.Values));
end
end
plot(gampdf(TT,Data.KD_val_alpha ,Data.KD_val_beta),TT,'k')
if  isfield(Data,'F_fluc_real')
    for nm=1:length(Data.KD_val_real)
    line([0 1.1]*max(maxx),Data.KD_val_real(nm)*[1 1],'Color','g','LineStyle','--')
    end
end

for m=1:length(Data.KD_val_learned{mmk}(1,:))
line([0 1.1]*maxx ,Data.map_KD_val(m)*[1 1],'Color','m','LineStyle','--')
end

set(gca,'YTick',[])
xlabel('pdf')
ylim([0 max(max(Data.KD_val_learned{mmk}(burn_in:end,:)))])
bb7.Position=bb7.Position.*[1 1 0.6 1];












bb4=subplot(3,4,9);
for nm=1:Data.Num_states_A
    for mm=1:Data.Num_states_A
           if mm~=nm
plot(permute(Data.KA_rate_learned(nm,mm,:),[1,3,2]),'b')
hold on
plot(permute(Data.KA_rate_learned(nm,mm,1:burn_in),[1,3,2]),'r')
line([1 size(Data.KA_rate_learned,3)] ,Data.map_KA_rate(nm,mm)*[1 1],'Color','m','LineStyle','--')
           end
    end
end

if isfield(Data,'F_fluc_real')
   for mn=1:length(Data.KA_real)
       for mm=1:length(Data.KA_real)
           if mm~=mn
       line([1 size(Data.KA_rate_learned,3)],Data.KA_real(mn,mm)*[1 1],'Color','g','LineStyle','--')
           end
       end
   end
end
xlabel('Iteration')
ylabel('K_A (s^{-1})')
if Data.Num_states_A*(Data.Num_states_A-1)>0
ylim([0 max(max(max(Data.KA_rate_learned(:,:,burn_in:end))))+0.001])
end
bb4.Position=bb4.Position.*[1 1 1.1 1];


bb5=subplot(3,4,10);

if Data.Num_states_A*(Data.Num_states_A-1)>0
TT=linspace(0,max(max(max(Data.KA_rate_learned))),100);
end
maxx=0;
for nm=1:Data.Num_states_A
    for mm=1:Data.Num_states_A
           if mm~=nm
hh1=histogram(Data.KA_rate_learned(nm,mm,burn_in:end),TT,...
    'Normalization','pdf','Orientation','Horizontal','EdgeColor','none','FaceColor','b');
hold on
if maxx<max(unique(hh1.Values))
    maxx=max(unique(hh1.Values));
end
           end
    end
end
plot(gampdf(TT,Data.KA_rate_alpha ,Data.KA_rate_beta),TT,'k')

for nm=1:Data.Num_states_A
    for mm=1:Data.Num_states_A
        if mm~=nm
          line([0 1.1]*maxx ,Data.map_KA_rate(nm,mm)*[1 1],'Color','m','LineStyle','--')
        end
    end
end
    
if isfield(Data,'F_fluc_real')
   for mn=1:length(Data.KA_real)
       for mm=1:length(Data.KA_real)
           if mm~=mn
       line([0 1.1]*maxx,Data.KA_real(mn,mm)*[1 1],'Color','g','LineStyle','--')
           end
       end
   end
end


set(gca,'YTick',[])
xlabel('pdf')
if Data.Num_states_A>0
ylim([0 max(max(max(Data.KA_rate_learned(:,:,burn_in:end))))+0.001])
end
bb5.Position=bb5.Position.*[1 1 0.6 1];









bb6=subplot(3,4,11);
for nm=1:Data.Num_states_A
plot(Data.KA_val_learned{mmk}(:,nm),'color','b')
hold on
plot(Data.KA_val_learned{mmk}(1:burn_in,nm),'color','r')
line([1 size(Data.KA_rate_learned,3)] ,Data.map_KA_val(nm)*[1 1],'Color','m','LineStyle','--')
end
if  isfield(Data,'F_fluc_real')
    for nm=1:length(Data.KA_val_real)
    line([1 size(Data.KA_rate_learned,3)] ,Data.KA_val_real(nm)*[1 1],'Color','g','LineStyle','--')
    end
end


xlabel('Iteration')
ylabel('\tau_A (ns)')
ylim([0 max(max(Data.KA_val_learned{mmk}(burn_in:end,:)))])
bb6.Position=bb6.Position.*[1 1 1.1 1];





bb7=subplot(3,4,12);
maxx=0;
TT=linspace(0,max(max(Data.KA_val_learned{mmk})),100);
for nm=1:Data.Num_states_A 
hh1=histogram(Data.KA_val_learned{mmk}(burn_in:end,nm),TT,...
    'Normalization','pdf','Orientation','Horizontal','EdgeColor','none','FaceColor','b');
hold on
if maxx<max(unique(hh1.Values))
   maxx=max(unique(hh1.Values));
end
end
plot(gampdf(TT,Data.KA_val_alpha ,Data.KA_val_beta),TT,'k')
if  isfield(Data,'F_fluc_real')
    for nm=1:length(Data.KA_val_real)
    line([0 1.1]*max(maxx),Data.KA_val_real(nm)*[1 1],'Color','g','LineStyle','--')
    end
end

for m=1:length(Data.KA_val_learned{mmk}(1,:))
    line([0 1.1]*maxx ,Data.map_KA_val(m)*[1 1],'Color','m','LineStyle','--')
end

set(gca,'YTick',[])
xlabel('pdf')
ylim([0 max(max(Data.KA_val_learned{mmk}(burn_in:end,:)))])
bb7.Position=bb7.Position.*[1 1 0.6 1];
