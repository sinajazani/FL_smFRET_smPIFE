function[] = Plot_mu_ex(Data,burn_in_perc)



burn_in       = max([1 floor((burn_in_perc/100)*length(Data.KA_val_learned(1,:)))]);


subplot(4,1,1)
for mj=1:Data.Num_states_D
    for mmk=1:length(Data.mu_D_ex_learned)
    plot(Data.mu_D_ex_learned{mmk}(:,mj),'b')
    hold on
    plot(Data.mu_D_ex_learned{mmk}(1:burn_in,mj),'r')
    
    line([1,length(Data.mu_D_ex_learned{mmk}(:,mj))],Data.map_mu_D_ex{mmk}(mj).*[1 1],'Color','m')
    end
end
if  isfield(Data,'ex_D_real')
    for mj=1:length(Data.ex_D_real)
        line([1 length(Data.mu_D_ex_learned{1}(:,1))] ,Data.ex_D_real(mj)*[1 1],'Color','g','LineStyle','--')
    end
end
ylabel('\mu_{D,ex}')


subplot(4,1,2)
for mj=1:Data.Num_states_A
    for mmk=1:length(Data.mu_D_ex_learned)
    plot(Data.mu_A_ex_learned{mmk}(:,mj),'b')
    hold on
    plot(Data.mu_A_ex_learned{mmk}(1:burn_in,mj),'r')
    line([1,length(Data.mu_A_ex_learned{mmk}(:,mj))],Data.map_mu_A_ex{mmk}(mj).*[1 1],'Color','m')
    end
end
if isfield(Data,'ex_A_real')
    for mj=1:length(Data.ex_A_real)
        line([1 length(Data.mu_D_ex_learned{1}(:,1))] ,Data.ex_A_real(mj)*[1 1],'Color','g','LineStyle','--')
    end
end
ylabel('\mu_{A,ex}')


subplot(4,1,3)
for mmk=1:length(Data.mu_D_ex_learned)
plot(Data.mu_BD_learned(mmk,:),'b')
hold on
plot(Data.mu_BD_learned(mmk,1:burn_in),'r')
line([1,length(Data.mu_BD_learned(1,:))],Data.map_mu_BD(mmk).*[1 1],'Color','m')
end
if isfield(Data,'Back_D_real')
line([1 length(Data.mu_D_ex_learned{1}(:,1))],[1 1]*Data.pulse_priod*Data.Back_D_real*10^-9,'Color','g','LineStyle','--')
end
ylabel('\mu_{BD}')

subplot(4,1,4)
for mmk=1:length(Data.mu_D_ex_learned)
plot(Data.mu_BA_learned(mmk,:),'b')
hold on
plot(Data.mu_BA_learned(mmk,1:burn_in),'r')
line([1,length(Data.mu_BA_learned(1,:))],Data.map_mu_BA(mmk).*[1 1],'Color','m')
end
if isfield(Data,'Back_A_real')
line([1 length(Data.mu_D_ex_learned{1}(:,1))],[1 1]*Data.pulse_priod*Data.Back_A_real*10^-9,'Color','g','LineStyle','--')
end
ylabel('\mu_{BA}')




end


