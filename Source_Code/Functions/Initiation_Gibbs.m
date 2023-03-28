%% Pre-calculation and preparing the parameter to be used inside the Gibbs Sampler
%
%  This code is written by Sina Jazani (09/14/2021)
%  Contact: sjazani1@jhu.edu

function[tester            , map              , KF_rate_learned  , KD_rate_learned   , ...
         IRF_D_str2        , IRF_A_str2       , eta_A_learned    , signal_time       , ...
         sigs_size         , eta_D_learned    , mu_D_ex_learned  , mu_A_ex_learned   , ...
         mu_BD_learned     , mu_BA_learned    , pi_accept        , prop_mu_D_ex      , ...
         alpha_mu_D_ex     , beta_mu_D_ex     , prop_mu_A_ex     , alpha_mu_A_ex     , ...
         beta_mu_A_ex      , max_bin          , prop_mu_BD       , alpha_mu_BD       , ...
         beta_mu_BD        , prop_mu_BA       , alpha_mu_BA      , beta_mu_BA        , ...
         KF_val_prop       , KF_val_alpha     , KF_val_beta      , AII               , ...
         A_all             , sigg11           , KD_val_prop      , KD_val_alpha      , ...
         KD_val_beta       , KD_val_learned   , KF_val_learned   , sigg22            , ...
         FII               , DII              , F_res            , tau_D_res         , ...
         indx1             , indx2            , sigg1            , sigg2             , ...
         pi_learned        , KF_rate_alpha    , KF_rate_beta     , Num_states_F      , ...
         KD_rate_alpha     , KD_rate_beta     , Num_states_D     , F_all             , ...
         D_all             , F_all_alpha      , F_all_base       , D_all_alpha       , ...
         D_all_base        , Num_states_A     , KA_rate_alpha    , KA_rate_beta      , ...
         KA_val_prop       , KA_val_alpha     , KA_val_beta      , KA_rate_learned   , ...
         A_all_alpha       , A_all_base       , tau_A_res        , KA_val_learned    , ...
         IRFDDstr2         , IRFAAstr2        ]= Initiation_Gibbs(Data)



   tester            = false;
   
   Num_states_F      = Data.Num_states_F ;
   Num_states_D      = Data.Num_states_D ;
   Num_states_A      = Data.Num_states_A ;
   
   KF_rate_alpha     = Data.KF_rate_alpha ;
   KF_rate_beta      = Data.KF_rate_beta;
   KD_rate_alpha     = Data.KD_rate_alpha ;
   KD_rate_beta      = Data.KD_rate_beta;
   KA_rate_alpha     = Data.KA_rate_alpha ;
   KA_rate_beta      = Data.KA_rate_beta;
   
   map               = Data.map(end);

   IRF_D_str         = Data.IRF_str(1);
   IRF_A_str         = Data.IRF_str(2);
   
   IRF_D_str_sqrt2   = IRF_D_str*sqrt(2) ;
   IRF_D_str2        = (IRF_D_str^2)       ;
   IRF_A_str_sqrt2   = IRF_A_str*sqrt(2) ;
   IRF_A_str2        = (IRF_A_str^2)       ;
   
   IRFDDstr2 = IRF_D_str./sqrt(2);
   IRFAAstr2 = IRF_A_str./sqrt(2);

   eta_A_learned     = Data.eta_A_learned(end);
   eta_D_learned     = Data.eta_D_learned(end);
   
   mu_BD_learned     = Data.mu_BD_learned(:,end);
   mu_BA_learned     = Data.mu_BA_learned(:,end);
   
   pi_accept         = Data.pi_accept           ;
   prop_mu_D_ex      = Data.prop_mu_D_ex          ;
   alpha_mu_D_ex     = Data.alpha_mu_D_ex         ; 
   beta_mu_D_ex      = Data.beta_mu_D_ex          ;
   
   prop_mu_A_ex      = Data.prop_mu_A_ex       ;
   alpha_mu_A_ex     = Data.alpha_mu_A_ex      ;
   beta_mu_A_ex      = Data.beta_mu_A_ex       ;
   
   prop_mu_BD        = Data.prop_mu_BD          ;
   alpha_mu_BD       = Data.alpha_mu_BD         ;
   beta_mu_BD        = Data.beta_mu_BD          ;
   
   prop_mu_BA        = Data.prop_mu_BA          ;
   alpha_mu_BA       = Data.alpha_mu_BA         ;
   beta_mu_BA        = Data.beta_mu_BA          ; 
   
   KF_val_prop       = Data.KF_val_prop         ;
   KF_val_alpha      = Data.KF_val_alpha        ;
   KF_val_beta       = Data.KF_val_beta         ;
   KD_val_prop       = Data.KD_val_prop         ;
   KD_val_alpha      = Data.KD_val_alpha        ;
   KD_val_beta       = Data.KD_val_beta         ;
   KA_val_prop       = Data.KA_val_prop         ;
   KA_val_alpha      = Data.KA_val_alpha        ;
   KA_val_beta       = Data.KA_val_beta         ;
if isfield(Data,'max_bin')
    if  Data.max_bin<Data.pulse_priod 
        max_bin      = Data.max_bin./2;
    else
        max_bin      = Data.pulse_priod./2      ;
    end
else
    max_bin          = Data.pulse_priod./2      ;
end


   
        
   signal_time           = [];
   indx1                 = [];
   indx2                 = [];
   pi_learned            = []; 

   
   KD_rate_learned       = [];
   KF_rate_learned       = [];
   KA_rate_learned       = [];
   
   KD_rate_learned(:,:,1)  = Data.KD_rate_learned(:,:,end);
   KF_rate_learned(:,:,1)  = Data.KF_rate_learned(:,:,end);
   KA_rate_learned(:,:,1)  = Data.KA_rate_learned(:,:,end);
   
   
   F_all                 = Data.F_all(:,end) ;
   D_all                 = Data.D_all(:,end) ;
   A_all                 = Data.A_all(:,end) ;
   
   
   
   F_all_alpha           = Data.F_all_alpha  ;
   F_all_base            = Data.F_all_base ;

   D_all_alpha           = Data.D_all_alpha  ;
   D_all_base            = Data.D_all_base ;
   
   A_all_alpha           = Data.A_all_alpha  ;
   A_all_base            = Data.A_all_base ;
   
   for  mmk=1:Data.num_sub_sigs  
       
        mu_D_ex_learned{mmk}   = Data.mu_D_ex_learned{mmk}(end,:);
        mu_A_ex_learned{mmk}   = Data.mu_A_ex_learned{mmk}(end,:);
        
        KD_val_learned(mmk,:)   = Data.KD_val_learned{mmk}(end,:);
        KF_val_learned(mmk,:)   = Data.KF_val_learned{mmk}(end,:);
        KA_val_learned(mmk,:)   = Data.KA_val_learned{mmk}(end,:);

        TF                     = Data.KF_rate_learned(:,:,end)-diag(sum(Data.KF_rate_learned(:,:,end),1));
        TD                     = Data.KD_rate_learned(:,:,end)-diag(sum(Data.KD_rate_learned(:,:,end),1));
        TA                     = Data.KA_rate_learned(:,:,end)-diag(sum(Data.KA_rate_learned(:,:,end),1));
        
        sigs_size(mmk,1)       = Data.sigs_size(mmk)           ;
        signal_time{mmk}       = diff(Data.signals{mmk}(2,:) ) ;
        

        [vecF , valF]=eig(TF); vallF = diag(valF); vec1F=vecF^-1;
        [vecD , valD]=eig(TD); vallD = diag(valD); vec1D=vecD^-1;
        [vecA , valA]=eig(TA); vallA = diag(valA); vec1A=vecA^-1;

        for n = 1:sigs_size(mmk)-1
            FII{mmk}(:,:,n)    = real(vecF*diag(exp(vallF*signal_time{mmk}(n)))*(vec1F));
            DII{mmk}(:,:,n)    = real(vecD*diag(exp(vallD*signal_time{mmk}(n)))*(vec1D));
            AII{mmk}(:,:,n)    = real(vecA*diag(exp(vallA*signal_time{mmk}(n)))*(vec1A));
        end
            
        F_res{mmk}             = Data.F_res{mmk}(end,:)    ;
        tau_D_res{mmk}         = Data.tau_D_res{mmk}(end,:);
        tau_A_res{mmk}         = Data.tau_A_res{mmk}(end,:);
      
        
        indx1{mmk}             = find(Data.signals{mmk}(3,:)==1); % index of Donor channel detected photons
        indx2{mmk}             = find(Data.signals{mmk}(3,:)==2); % index of Acceptor channel detected  photons
        
        sigg1{mmk}             = Data.IRF_mean(1)-Data.signals{mmk}(1,indx1{mmk});
        sigg2{mmk}             = Data.IRF_mean(2)-Data.signals{mmk}(1,indx2{mmk});
        
        sigg11{mmk}            = sigg1{mmk}./IRF_D_str_sqrt2;
        sigg22{mmk}            = sigg2{mmk}./IRF_A_str_sqrt2;


        A1                     = exp(-mu_D_ex_learned{mmk}(end,tau_D_res{mmk}));
        A2                     = exp(-mu_A_ex_learned{mmk}(end,tau_A_res{mmk}));
        A3                     = exp(-mu_BD_learned(mmk,1));
        A4                     = exp(-mu_BA_learned(mmk,1));
        
        B                      = [(1-A1).*A2.*A3.*A4 ; (1-A2).*A1.*A3.*A4 ; (1-A3).*A1.*A2.*A4 ; (1-A4).*A1.*A2.*A3 ];
        pi_learned{mmk}        = (B./sum(B,1));
    end     

end