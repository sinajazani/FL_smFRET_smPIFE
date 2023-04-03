function[ Data ] = Initiation_param(Data)


if  isfield(Data,'KD_real')
    Data.signals           = []          ; 
    Data.signals{1}        = Data.signal ;
end

Data.F_all_base            = ones(Data.Num_states_F,1)./Data.Num_states_F;
Data.F_all                 = dirichletRnd(Data.F_all_alpha*Data.F_all_base) ;


Data.D_all_base            = ones(Data.Num_states_D,1)./Data.Num_states_D;
Data.D_all                 = dirichletRnd(Data.D_all_alpha*Data.D_all_base) ;

Data.A_all_base            = ones(Data.Num_states_A,1)./Data.Num_states_A;
Data.A_all                 = dirichletRnd(Data.A_all_alpha*Data.A_all_base) ;

Data.sigs_size             = []; 

Data.F_res                 = []; 
Data.tau_D_res             = []; 
Data.tau_A_res             = []; 

Data.KF_val_learned        = []; 
Data.KD_val_learned        = []; 
Data.KA_val_learned        = []; 

Data.KF_rate_learned       = [];
Data.KD_rate_learned       = []; 
Data.KA_rate_learned       = [];

Data.mu_D_ex_learned       = []; 
Data.mu_A_ex_learned       = [];

Data.mu_BD_learned         = []; 
Data.mu_BA_learned         = []; 


% Number of individual time traces
Data.num_sub_sigs = length(Data.signals);

% Constructing of the transition matrix
Data.KF_rate_learned     = gamrnd(Data.KF_rate_alpha ,Data.KF_rate_beta ,Data.Num_states_F,Data.Num_states_F).*(ones(Data.Num_states_F)-eye(Data.Num_states_F)) ;
Data.KD_rate_learned     = gamrnd(Data.KD_rate_alpha ,Data.KD_rate_beta ,Data.Num_states_D,Data.Num_states_D).*(ones(Data.Num_states_D)-eye(Data.Num_states_D)) ;
Data.KA_rate_learned     = gamrnd(Data.KA_rate_alpha ,Data.KA_rate_beta ,Data.Num_states_A,Data.Num_states_A).*(ones(Data.Num_states_A)-eye(Data.Num_states_A)) ;

mmtF                     = Data.KF_rate_learned-diag(sum(Data.KF_rate_learned,1));
mmtD                     = Data.KD_rate_learned-diag(sum(Data.KD_rate_learned,1));
mmtA                     = Data.KA_rate_learned-diag(sum(Data.KA_rate_learned,1));


for mmk=1:Data.num_sub_sigs
    
    Data.KA_val_learned{mmk}      = gamrnd(Data.KA_val_alpha   , Data.KA_val_beta   , 1,Data.Num_states_A) ;
    Data.KD_val_learned{mmk}      = gamrnd(Data.KD_val_alpha   , Data.KD_val_beta   , 1,Data.Num_states_D) ;
    Data.KF_val_learned{mmk}      = gamrnd(Data.KF_val_alpha   , Data.KF_val_beta   , 1,Data.Num_states_F) ;

    
    Data.mu_D_ex_learned{mmk}(1,:)  = gamrnd(Data.alpha_mu_D_ex  , Data.beta_mu_D_ex  , 1 , Data.Num_states_D ) ;
    Data.mu_A_ex_learned{mmk}(1,:)  = gamrnd(Data.alpha_mu_A_ex  , Data.beta_mu_A_ex  , 1 , Data.Num_states_A ) ;
    
    Data.mu_BD_learned(mmk,1)       = gamrnd(Data.alpha_mu_BD    , Data.beta_mu_BD    )                ;
    Data.mu_BA_learned(mmk,1)       = gamrnd(Data.alpha_mu_BA    , Data.beta_mu_BA    )                ;
    
    Data.F_res{mmk}(1)              = Discrete_sampler(Data.F_all)                                ;
    Data.tau_D_res{mmk}(1)          = Discrete_sampler(Data.D_all)                                ;
    Data.tau_A_res{mmk}(1)          = Discrete_sampler(Data.A_all)                                ;
    
    Data.sigs_size(mmk)             = size(Data.signals{mmk},2)                                ;
    for jl=1:Data.sigs_size(mmk)-1
        TA                          = expm(mmtA*(Data.signals{mmk}(2,jl+1)-Data.signals{mmk}(2,jl))) ;
        Data.tau_A_res{mmk}(1,jl+1) = Discrete_sampler(TA(:,Data.tau_A_res{mmk}(jl)))                ;
        TF                          = expm(mmtF*(Data.signals{mmk}(2,jl+1)-Data.signals{mmk}(2,jl))) ;
        Data.F_res{mmk}(1,jl+1)     = Discrete_sampler(TF(:,Data.F_res{mmk}(jl)))                    ;
        TD                          = expm(mmtD*(Data.signals{mmk}(2,jl+1)-Data.signals{mmk}(2,jl))) ;
        Data.tau_D_res{mmk}(1,jl+1) = Discrete_sampler(TD(:,Data.tau_D_res{mmk}(jl)))                ;
    end
end

Data.tau_accept                     = zeros(2,Data.num_sub_sigs)  ;
Data.pi_accept                      = zeros(2,Data.num_sub_sigs)  ;
Data.eta_accept                     = [0 ; 0]  ;


% Set the values of MAP estimates
Data.map             = -Inf                            ;
Data.RAW_map         = Data.map                        ;
Data.map_KF_rate     = Data.KF_rate_learned(:,end)     ;
Data.map_KD_rate     = Data.KD_rate_learned(:,end)     ;
Data.map_KA_rate     = Data.KA_rate_learned(:,end)     ;

Data.map_F_all       = Data.F_all(:,end)               ; 
Data.map_D_all       = Data.D_all(:,end)               ;
Data.map_A_all       = Data.A_all(:,end)               ;

Data.map_KA_val      = Data.KA_val_learned             ;
Data.map_KD_val      = Data.KD_val_learned             ;
Data.map_KF_val      = Data.KF_val_learned             ;

Data.map_F_res       = Data.F_res                      ; 
Data.map_tau_D       = Data.tau_D_res                  ;
Data.map_tau_A       = Data.tau_A_res                  ;

Data.map_mu_D_ex     = Data.mu_D_ex_learned            ;
Data.map_mu_A_ex     = Data.mu_A_ex_learned            ;

Data.map_mu_BD       = Data.mu_BD_learned              ; 
Data.map_mu_BA       = Data.mu_BA_learned              ;

Data.map_eta_D       = Data.eta_D_learned(end)         ;
Data.map_eta_A       = Data.eta_A_learned(end)         ; 

end
