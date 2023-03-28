function[]=Plot_synthetic_data(Data)
subplot(5,1,1)

stairs(Data.signal(2,:),Data.F_fluc_real)
ylim([0 max([1 max(Data.F_fluc_real)*1.1])])
xlim([0 max(Data.signal(2,:))])
ylabel('K_F (ns^{-1})')
xlabel('Time (s)')

subplot(5,1,2)
stairs(Data.signal(2,:),Data.tau_D_fluc_real)
ylim([0 max([1 max(Data.tau_D_fluc_real)*1.1])])
xlim([0 max(Data.signal(2,:))])
ylabel('\tau_D (ns)')
xlabel('Time (s)')

subplot(5,1,3)
stairs(Data.signal(2,:),Data.tau_A_fluc_real)
ylim([0 max([1 max(Data.tau_A_fluc_real)*1.1])])
xlim([0 max(Data.signal(2,:))])
ylabel('\tau_A (ns)')
xlabel('Time (s)')

subplot(5,1,4)
plot(Data.signal(2,Data.signal(3,:)==1),Data.signal(1,Data.signal(3,:)==1),'.','color','g','MarkerSize',2)
ylabel({'Donor','Photon arriving','time (ns)'})
xlabel('Time (s)')
xlim([0 max(Data.signal(2,:))])

subplot(5,1,5)
plot(Data.signal(2,Data.signal(3,:)==2),Data.signal(1,Data.signal(3,:)==2),'.','color','r','MarkerSize',2)
ylabel({'Acceptor','Photon arriving','time (ns)'})
xlabel('Time (s)')
xlim([0 max(Data.signal(2,:))])