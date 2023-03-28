function[Data] = Read_exp_files( Exp_type , Sigs , Acept_ratio )



Data.pulse_priod    =      50 ;
Data.max_bin        =      23.5 ;


% Add the folder contains of the experimental data sets into the directory path
addpath('exp_data_sets')

% Pre_calculations
mm           = 0 ;
Data.signals = [];


% Choose the data set you are interested to analyze
switch Exp_type
    case 1
         load('a3Ddata.mat')      % a3D data sets
         n=1;

         % n=1
         Data.IRF_D_mean  = 2.1+donT0(n)*0.04 ; % Offset of the IRF in nano second
         Data.IRF_A_mean  = 2.35+accT0(n)*0.04; % Offset of the IRF in nano second
    
         % n=2
%          Data.IRF_D_mean  = 2.2+donT0(n)*0.04 ; % Offset of the IRF in nano second
%          Data.IRF_A_mean  = 2.55+accT0(n)*0.042; % Offset of the IRF in nano second   

         % Standard deviation of the IRF
         Data.IRF_D_str   = 0.26   ;
         Data.IRF_A_str   = 0.26   ;

         for m=1:length(Sigs)
             Data.signal = [];
             signn = photontrajectories{n}(cumindexall{n}(Sigs(m))+1:cumindexall{n}(Sigs(m)+1),3)';
             deett = photontrajectories{n}(cumindexall{n}(Sigs(m))+1:cumindexall{n}(Sigs(m)+1),4)';
             ddg   = photontrajectories{n}(cumindexall{n}(Sigs(m))+1:cumindexall{n}(Sigs(m)+1),2)'.*10^-3 ;
             ddg   = ddg - min(ddg);

             binaa            = (rand(1,size(deett,2))<Acept_ratio);
             deett=(deett-1);  deett(deett==0) = 2;  Data.signal(3,:) = deett;
    
             signn(Data.signal(3,:)==1) = (signn(Data.signal(3,:)==1)-4000)*2./1000+(donT0(n)*0.04)+0.2;
             signn(Data.signal(3,:)==2) = (signn(Data.signal(3,:)==2)-3600)*2./1000+(accT0(n)*0.04)+0.8;

             Data.signal      = Data.signal(:,binaa);
             Data.signal(1,:) = signn(binaa);Data.signal(2,:) = ddg(binaa);mm=mm+1;Data.signals{mm} = Data.signal;
         end
    case 2
         load('gpWdata.mat')     % gpW data sets

         %%% gpW IRF
         Data.IRF_D_mean  = 2 ; % Offset of the IRF in nano second
         Data.IRF_A_mean  = 1.5; % Offset of the IRF in nano second

         % Standard deviation of the IRF
         Data.IRF_D_str  =  0.2      ; 
         Data.IRF_A_str   = 0.2   ;

         for m=1:length(Sigs)
             Data.signal = [];
    
             signn = photontrajectories(cumindexall(Sigs(m))+1:cumindexall(Sigs(m)+1),3)';
             deett = photontrajectories(cumindexall(Sigs(m))+1:cumindexall(Sigs(m)+1),4)';
             ddg   = photontrajectories(cumindexall(Sigs(m))+1:cumindexall(Sigs(m)+1),2)'.*10^-3 ;
             ddg   = ddg - min(ddg);
    
             binaa            = (rand(1,size(deett,2))<Acept_ratio);
             deett=(deett-1);  deett(deett==0) = 2;  Data.signal(3,:) = deett;
    
             signn(Data.signal(3,:)==1) = (signn(Data.signal(3,:)==1)-4000)*2./1000;
             signn(Data.signal(3,:)==2) = (signn(Data.signal(3,:)==2)-3750)*2./1000;

             Data.signal      = Data.signal(:,binaa);
             Data.signal(1,:) = signn(binaa);Data.signal(2,:) = ddg(binaa);mm=mm+1;Data.signals{mm} = Data.signal;
         end

    case 3
         load('WWdata.mat')     % WWdomain data sets 

         %%% WW IRF
         Data.IRF_D_mean  = 2 ; % Offset of the IRF in nano second
         Data.IRF_A_mean  = 1.5 ; % Offset of the IRF in nano second

         % Standard deviation of the IRF
         Data.IRF_D_str   = 0.2   ; 
         Data.IRF_A_str   = 0.2   ;

         for m=1:length(Sigs)
             Data.signal = [];
    
             signn = photontrajectories(cumindexall(Sigs(m))+1:cumindexall(Sigs(m)+1),3)';
             deett = photontrajectories(cumindexall(Sigs(m))+1:cumindexall(Sigs(m)+1),4)';
             ddg   = photontrajectories(cumindexall(Sigs(m))+1:cumindexall(Sigs(m)+1),2)'.*10^-3 ;
             ddg   = ddg - min(ddg);
    
             binaa            = (rand(1,size(deett,2))<Acept_ratio);
             deett=(deett-1);  deett(deett==0) = 2;  Data.signal(3,:) = deett;
    
             signn(Data.signal(3,:)==1) = (signn(Data.signal(3,:)==1)-3910)*2./1000;
             signn(Data.signal(3,:)==2) = (signn(Data.signal(3,:)==2)-3660)*2./1000;

             Data.signal      = Data.signal(:,binaa);
             Data.signal(1,:) = signn(binaa);Data.signal(2,:) = ddg(binaa);mm=mm+1;Data.signals{mm} = Data.signal;
         end

end







% dde = load('DonorIRFCH2.dat')    ;
% aae = load('AcceptorIRFCH1.dat') ;


subplot(4,1,1)
dbnd = 0:.1:Data.pulse_priod ;
dbnd1 = 0:.01:Data.pulse_priod;
sig=[]; tt=0;
for mf=1:length(Data.signals)
    fig=Data.signals{mf}; fig(2,:)=fig(2,:)+tt(mf); sig=[sig,fig]; tt(mf+1) = sig(2,end);
end
histogram(sig(1,sig(3,:)==1)  ,dbnd,'Normalization','pdf','FaceColor','g','FaceAlpha',0.5); hold on
%hh=plot((0:length(dde)-4000)*0.002 ,150*dde(4000:end)./sum(dde(4000:end)),'k');
plot(dbnd1, 0.5*exp(-0.5*((dbnd1-Data.IRF_D_mean)./Data.IRF_D_str).^2),'k' )
%legend(hh,'IRF','EdgeColor',[1 1 1])
xlim([ 0 20]); ylabel('Donor'); set(gca,'YTick',[]); box off


subplot(4,1,2)
histogram(sig(1,sig(3,:)==2)  ,dbnd,'Normalization','pdf','FaceColor','r','FaceAlpha',0.5); hold on
%plot((0:length(aae)-3600).*0.002 ,60*aae(3600:end)./sum(aae(3600:end)),'k')
plot(dbnd1, 0.23*exp(-0.5*((dbnd1-Data.IRF_A_mean)./Data.IRF_A_str).^2) ,'k')
xlim([ 0 20]); box off; ylabel('Acceptor'); set(gca,'YTick',[]); xlabel('Delay time (ns)')

% Plot the Single photon arrivals collected by the Donor detector
subplot(4,1,3)
plot(sig(2,sig(3,:)==1),sig(1,sig(3,:)==1),'.','color','b')
hold on
for mf=1:length(Data.signals)+1
    line(tt(mf)*[1 1],[0 Data.pulse_priod],'LineStyle','-','Color','k')
end
ylabel({'Donor','Photon arriving time (ns)'})
xlabel('Time (s)'); xlim([0 max(sig(2,:))]); ylim([ 0 20])

% Plot the Single photon arrivals collected by the Acceptor detector
subplot(4,1,4)
plot(sig(2,sig(3,:)==2),sig(1,sig(3,:)==2),'.','color','r')
hold on
for mf=1:length(Data.signals)+1
    line(tt(mf)*[1 1],[0 Data.pulse_priod],'LineStyle','-','Color','k')
end
ylabel({'Acceptor','Photon arriving time (ns)'}); xlabel('Time (s)'); xlim([0 max(sig(2,:))]);ylim([ 0 20])


disp(['Total number of traces=' num2str(length(Data.signals))])
disp(['Total number of photons=' num2str(length(sig(2,:)))])


end