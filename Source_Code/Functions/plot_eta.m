function[]=plot_eta(Data, burn_in_perc)

% Select the burn-in part of the data
burn_in      = max([1 floor((burn_in_perc/100)*length(Data.eta_D_learned))]);

dbnd=linspace(0,1,50);

subplot(2,2,1)
plot(Data.eta_D_learned(1:end))
hold on
plot(Data.eta_D_learned(1:burn_in),'r')
if isfield(Data,'eta_D_real')
   line([1 length(Data.eta_D_learned)],Data.eta_D_real.*[1 1],'color','g')
end
xlabel('Iteration')
ylabel('Donor channel efficiency')
ylim([0 1])

subplot(2,2,2)
hhD=histogram(Data.eta_D_learned(burn_in:end),dbnd,'Normalization','pdf','Orientation','Horizontal');
xlabel('Post. prob. distr.')
hold on
xx=0:0.01:1;
plot((max(unique(hhD.Values))).*(xx.^(Data.alpha_eta_D-1)).*((1-xx).^(Data.beta_eta_D-1)),xx)
ylim([0 1])



subplot(2,2,3)
plot(Data.eta_A_learned(1:end))
hold on
plot(Data.eta_A_learned(1:burn_in),'r')
if isfield(Data,'eta_A_real')
   line([1 length(Data.eta_D_learned)],Data.eta_A_real.*[1 1],'color','g')
end
xlabel('Iteration')
ylabel('Acceptor channel efficiency')
ylim([0 1])



subplot(2,2,4)
hhA=histogram(Data.eta_A_learned(burn_in:end),dbnd,'Normalization','pdf','Orientation','Horizontal');
xlabel('Post. prob. distr.')
hold on
xx=0:0.01:1;
plot((max(unique(hhA.Values))).*(xx.^(Data.alpha_eta_A-1)).*((1-xx).^(Data.beta_eta_A-1)),xx)
ylim([0 1])