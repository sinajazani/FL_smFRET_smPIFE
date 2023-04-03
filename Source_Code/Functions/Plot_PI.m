function[]=Plot_PI(Data,burn_in_perc)

burn_in       = max([1 floor((burn_in_perc/100)*length(Data.tau_A_res{1}(:,1)))]);

dd1           = [] ; dd2          = [] ; sigg_n2      = [] ; max_1(1)     =  0 ;
for mmk=1:Data.num_sub_sigs
    if  mmk==1
       
        sigg_n2 =Data.signals{mmk}(2,:);
       
        max_1(2) = max(unique(sigg_n2));
    else
        sigg_n2 =[sigg_n2, max_1(mmk)+Data.signals{mmk}(2,:)];
        max_1(mmk+1) = max(sigg_n2);
    end
    ttp=[];
    for ii=burn_in:size(Data.tau_D_res{1}(:,:),1)
        ttp(ii-burn_in+1,:)=Data.mu_D_ex_learned{mmk}(ii-1,Data.tau_D_res{mmk}(ii-1,:));
    end
    dd1=[dd1;exp(-mean(ttp))];
    
    ttp1=[];
    for ii=burn_in:size(Data.tau_A_res{1}(:,:),1)
        ttp1(ii-burn_in+1,:)=Data.mu_A_ex_learned{mmk}(ii-1,Data.tau_A_res{mmk}(ii-1,:));
    end
    dd2=[dd2;exp(-mean(ttp1))];
    
end
dd3 = exp(-mean(Data.mu_BD_learned(burn_in:end)));
dd4 = exp(-mean(Data.mu_BA_learned(burn_in:end)));

DD = [(1-dd1).*dd2.*dd3.*dd4 ; (1-dd2).*dd1.*dd3.*dd4 ; (1-dd3).*dd1.*dd2.*dd4 ; (1-dd4).*dd1.*dd2.*dd3];
DD = DD./sum(DD,1);

for m=1:4
    subplot(4,1,m)
    line(sigg_n2,DD(m,:)' ,'Color','b'  )
    hold on
    if isfield(Data,'PI_real')
       line(sigg_n2,Data.PI_real(m,:)' ,'Color','g','LineStyle','--'  )
    end
    ylim([0,1])
end
end