%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0
%%%%%   Clear the Wordspace & add the function folder to the path    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

addpath('Functions')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1
%%%%%%%%%%%%%%%%%%%      Generate Synthetic Data      %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data.Num_pulses     =   1000000                                           ;
Data.pulse_priod    =        50                                           ; % The time scale is in nano second

% Backgroand photon emission rates for donor and acceptor channels
Data.Back_D_real    =       500                                           ; % Average number of photons per second collected in donor channel
Data.Back_A_real    =       500                                           ; % Average number of photons per second collected in acceptor channel

% IRFS (Gausian)
Data.IRF_mean       =       [ 1.914 , 1.914 ]                             ; % Offset of the IRF in nano second for [Donor , Acceptor]
Data.IRF_str        =       [ 0.22  , 0.24  ]                             ; % Standard deviation of the IRF for [Donor , Acceptor]

% Detection cross-talks of bouth channels
Data.eta_D_real     =         0.05                                        ;
Data.eta_A_real     =         0.05                                        ;

Data.ex_A_real      =        10^-6*[ 1  1 ]                               ; % Absorbsion cross-talk of the acceptor rate times pulse width
Data.ex_D_real      =        10^-3*[ 2  2 ]                               ; % Absorbsion cross-talk of the donor rate times pulse width

Data.KF_val_real    =       [  0.1  , 0.5  ]                              ; % FRET rates of each state
Data.KD_val_real    =       [  4    , 4    ]                              ; % Donor lifetime of each PIFE state 
Data.KA_val_real    =       [  4    , 4    ]                              ; % Acceptor lifetime of each PIFE state of acceptor

% Transition rates
Data.KF_real        =       [ -1    , 2    ;...
                               1    ,-2    ]*50                           ;

Data.KD_real        =       [ -1    , 2    ; ...
                               1    ,-2    ]*50                           ;

Data.KA_real        =       [ -4    , 6    ; ...
                               4    ,-6    ]*50                           ;
                   
% Genrate the data
Data                =   Generative_model(Data)                            ;

% Plot the generated data set
Plot_synthetic_data(Data)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2
%%%%%%%%%%%%%%%%%       Import Experimental Data       %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

Exp_type        = 1        ;  % Choose the data set you are interested to analyze (1=a3D , 2=gpw , 3=WW domain )
Acept_ratio     = 1        ;  % You can reduce the number of photons, randomly selected based on the defined acceptance ratio. Value between [0 1]
Sigs            = [ 3,4  ] ; % Select the traces you want to analyze

Data            = Read_exp_files( Exp_type , Sigs , Acept_ratio );


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3
%%%%%%    Initilaize the log of background photon emission rate    %%%%%%%%
%%%%%%         per each pixel through inducing points              %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% Number of states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data.Num_states_F          =      2                                       ; % Number of FRET states
Data.Num_states_D          =      1                                       ; % Number of Donor_PIFE states
Data.Num_states_A          =      1                                       ; % Number of Acceptor_PIFE states


%%%%%%%%%%%%%%%%     Parameters of transition rates    %%%%%%%%%%%%%%%%%%%%
Data.KF_rate_alpha         =      1                                       ; 
Data.KF_rate_beta          =  10000                                       ;

Data.F_all_alpha           =      1                                       ;

Data.KD_rate_alpha         =      1                                       ;
Data.KD_rate_beta          =  10000                                       ;

Data.D_all_alpha           =      1                                       ;

Data.KA_rate_alpha         =      1                                       ;
Data.KA_rate_beta          =  10000                                       ;

Data.A_all_alpha           =      1                                       ;


%%%%%%%%%%%%%%%%     Parameters of excitation rates    %%%%%%%%%%%%%%%%%%%%
Data.prop_mu_D_ex          =   1000                                       ;
Data.alpha_mu_D_ex         =      1                                       ;
Data.beta_mu_D_ex          =      0.001                                   ;

Data.prop_mu_A_ex          =   1000                                       ;
Data.alpha_mu_A_ex         =      1                                       ;
Data.beta_mu_A_ex          =      0.001                                   ;


%%%%%%%%%%%%%%%%     Parameters of background rates    %%%%%%%%%%%%%%%%%%%%
Data.prop_mu_BD            =  10000                                       ;
Data.alpha_mu_BD           =      1                                       ;
Data.beta_mu_BD            =      0.001                                   ;

Data.prop_mu_BA            =  10000                                       ;
Data.alpha_mu_BA           =      1                                       ;
Data.beta_mu_BA            =      0.001                                   ;


%%%%%%%%%%%%%%%%     Parameters of cross-talk ratios     %%%%%%%%%%%%%%%%%%
Data.sample_eta            =   true   ; % Please set to 'false' if the cross-talk ratio is known

Data.eta_prop_D            =      0.5                                     ;
Data.alpha_eta_D           =      0.5                                     ;
Data.beta_eta_D            =     10                                       ;

Data.eta_prop_A            =      0.5                                     ;
Data.alpha_eta_A           =      0.5                                     ;
Data.beta_eta_A            =     10                                       ;

if  Data.sample_eta 
    Data.eta_D_learned     = betarnd(Data.alpha_eta_D , Data.beta_eta_D)  ;
    Data.eta_A_learned     = betarnd(Data.alpha_eta_A , Data.beta_eta_A)  ;
else % pre calibrated values
    Data.eta_D_learned     =      0.05                                    ;
    Data.eta_A_learned     =      0.05                                    ;
end

%%%%%%%%%%%%%%%%     Parameters of FRET rate values    %%%%%%%%%%%%%%%%%%%%
Data.KF_val_prop           =   1000                                       ;
Data.KF_val_alpha          =      1                                       ;
Data.KF_val_beta           =     10                                       ;

%%%%%%%%%%%%%%%%      Parameters of donor lifetimes    %%%%%%%%%%%%%%%%%%%%
Data.KD_val_prop           =   1000                                       ;
Data.KD_val_alpha          =      2                                       ;
Data.KD_val_beta           =      4                                       ;

%%%%%%%%%%%%%%%%     Parameters of acceptor lifetimes    %%%%%%%%%%%%%%%%%%
Data.KA_val_prop           =   1000                                       ;
Data.KA_val_alpha          =      2                                       ;
Data.KA_val_beta           =      4                                       ;

%%%%%%%%%%%%%%%%           Data saving control           %%%%%%%%%%%%%%%%%%
Data.save_size             =     10                                       ;
Data.save_on_off           = true                                         ;

%%%%%%%%%%%%%%%%        Parameter pre-calculations       %%%%%%%%%%%%%%%%%%
Data                       = Initiation_param(Data)                       ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4
%%%%%%%%%%%%%%%%%             Gibbs sampler            %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Max_iter  = 500 ;
tic; Data = Gibbs_sampler(Data ,  Max_iter); toc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5
%%%%%%%%%%%%%%%%%             PLOT RESULTS             %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAP estimatation of log joint posterior
plot(Data.map); hold on ; plot(Data.RAW_map) ; ylim([min([Data.map,Data.RAW_map]) max([Data.map,Data.RAW_map])+100]) ; set(gca,'YScale','log')

%% Transition rates
mmk           =   1     ; % Trace number (Generative model only creates one trace)
burn_in_perc  =  50     ;  Plot_transition_rates( Data , burn_in_perc ,mmk )

%% FRET+PIFE
mmk           =   1     ; % Trace number (Generative model only creates one trace)
bin_size      =   0.001 ; % Bin size in unit of second
burn_in_perc  =  90     ; Plot_FRET_tau_DandA( Data , burn_in_perc , bin_size ,mmk)
%% ETA
burn_in_perc  =  60     ; plot_eta(Data, burn_in_perc)

%%
burn_in_perc  =  60     ; Plot_mu_ex(Data,burn_in_perc)

%%
burn_in_perc  =  90     ; Plot_PI(Data,burn_in_perc)
