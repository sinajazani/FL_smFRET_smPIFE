clear all

% IRF fit ----------
chmergenum=20;
irfch1=cumsum(load('AcceptorIRFCH1.dat'));
irfch2=cumsum(load('DonorIRFCH2.dat'));

irf(:,1)=irfch1(chmergenum*(1:floor(length(irfch1))/chmergenum));
irf(:,2)=irfch2(chmergenum*(1:floor(length(irfch2))/chmergenum));
irf=[irf(1,:); diff(irf,1,1)];
figure;plot(1:length(irf(:,1)),irf,'o-');

irfcut=[irf(236:295,1) irf(261:320,2)];

% % double gaussian irf
% irffitparam=fitirf(1:length(irfcut(:,1)),irfcut,[[10; 5] [10; 5]]);
% figure;plot(1:length(irfcut),irfcut,'o');hold;
% xrightid1=find(1:length(irfcut) > irffitparam(1,1),1);
% plot(1:length(irfcut),irffitparam(4,1)*[exp(-((1:xrightid1-1)'-irffitparam(1,1)).^2/(2*irffitparam(2,1)^2)); exp(-((xrightid1:length(irfcut(:,1)))'-irffitparam(1,1)).^2/2/irffitparam(3,1)^2)] + irffitparam(5,1),'r-');
% xrightid1=find(1:length(irfcut) > irffitparam(1,2),1);
% plot(1:length(irfcut),irffitparam(4,2)*[exp(-((1:xrightid1-1)'-irffitparam(1,2)).^2/(2*irffitparam(2,2)^2)); exp(-((xrightid1:length(irfcut(:,1)))'-irffitparam(1,2)).^2/2/irffitparam(3,2)^2)] + irffitparam(5,2),'r-');
% hold;
% irfoff=diff(irffitparam([2 3],:),1,1)/sqrt(2*pi);

% gamma irf
irffitparam=fitirfgamma(1:length(irfcut(:,1)),irfcut,[[10; sqrt(6); sqrt(10); max(irfcut(:,1))] [10; sqrt(5); sqrt(10); max(irfcut(:,2))]]);
irffitparam(2:3,:)=irffitparam(2:3,:).^2;
%%
figure;
plot(1:length(irfcut),irfcut(:,1),'o','color','r')
hold on
plot(1:length(irfcut),irfcut(:,2),'o','color','g')
xrightid1=find(1:length(irfcut) > irffitparam(1,1),1);

plot(1:length(irfcut),irffitparam(4,1)/gamma(irffitparam(2,1))*[zeros(xrightid1-1,1); ...
    exp(-((xrightid1:length(irfcut(:,1)))'-irffitparam(1,1))*irffitparam(3,1)).*(((xrightid1:length(irfcut(:,1)))'-irffitparam(1,1))*irffitparam(3,1)).^(irffitparam(2,1)-1)],'r-','Color','r');


xrightid1=find(1:length(irfcut) > irffitparam(1,2),1);
plot(1:length(irfcut),irffitparam(4,2)/gamma(irffitparam(2,2))*[zeros(xrightid1-1,1); ...
    exp(-((xrightid1:length(irfcut(:,1)))'-irffitparam(1,2))*irffitparam(3,2)).*(((xrightid1:length(irfcut(:,1)))'-irffitparam(1,2))*irffitparam(3,2)).^(irffitparam(2,2)-1)],'r-','Color','g');
hold;
irfoff=irffitparam(2,:)./irffitparam(3,:);