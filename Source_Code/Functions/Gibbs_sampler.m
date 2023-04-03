%% The MCMC chain (Gibbs sampler) for single photon FRET+PIFE analysis
%  This function is main sampler
%
%  This code is written by Sina Jazani (09/14/2021)
%  Contact: sjazani1@jhu.edu

function[Data] = Gibbs_sampler(Data , Max_iter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Data                  :  A structure array contains of all parameters and variables
   % Max_iter              :  Number of iteration per run
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Data                  :  The updatetd structure array contains of all parameters and variables
   
   % Shuffle the random generator
   rng('shuffle');

   % Initiation of the parameters for the Parallel loop and Gibbs sampler
   [  tester            , map               , KF_rate_learned   , KD_rate_learned   , ...
      IRF_D_str2        , IRF_A_str2        , eta_A_learned     , signal_time       , ...
      sigs_size         , eta_D_learned     , mu_D_ex_learned   , mu_A_ex_learned   , ...
      mu_BD_learned     , mu_BA_learned     , pi_accept         , prop_mu_D_ex      , ...
      alpha_mu_D_ex     , beta_mu_D_ex      , prop_mu_A_ex      , alpha_mu_A_ex     , ...
      beta_mu_A_ex      , max_bin           , prop_mu_BD        , alpha_mu_BD       , ...
      beta_mu_BD        , prop_mu_BA        , alpha_mu_BA       , beta_mu_BA        , ...
      KF_val_prop       , KF_val_alpha      , KF_val_beta       , AII               , ...
      A_all             , sigg11            , KD_val_prop       , KD_val_alpha      , ...
      KD_val_beta       , KD_val_learned    , KF_val_learned    , sigg22            , ...
      FII               , DII               , F_res             , tau_D_res         , ...
      indx1             , indx2             , sigg1             , sigg2             , ...
      pi_learned        , KF_rate_alpha     , KF_rate_beta      , Num_states_F      , ...
      KD_rate_alpha     , KD_rate_beta      , Num_states_D      , F_all             , ...
      D_all             , F_all_alpha       , F_all_base        , D_all_alpha       , ...
      D_all_base        , Num_states_A      , KA_rate_alpha     , KA_rate_beta      , ...
      KA_val_prop       , KA_val_alpha      , KA_val_beta       , KA_rate_learned   , ...
      A_all_alpha       , A_all_base        , tau_A_res         , KA_val_learned    , ...
      IRFDDstr2         , IRFAAstr2         ]= Initiation_Gibbs(Data);

   
   % Start of the MCMC chain, Gibbs sampling iteration FOR loop
   for i=1:Max_iter
       
       for mmk=1:Data.num_sub_sigs
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%           Sample the FRET state trajectory            %%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if Num_states_F>1
              [ F_res{mmk}  ] = sampling_FRET_HMM( ...
                 ... 
              sigs_size(mmk)       , IRF_D_str2       , IRF_A_str2     , FII{mmk}                , ...
              KA_val_learned(mmk,:), tau_D_res{mmk}   , eta_D_learned  , KF_val_learned(mmk,:)'  , ...
              eta_A_learned        , pi_learned{mmk}  , max_bin        , KD_val_learned(mmk,:)   , ...
              indx1{mmk}           , indx2{mmk}       , sigg1{mmk}     , sigg2{mmk}              , ...
              F_all                , Num_states_F     , tau_A_res{mmk} , IRFDDstr2               , ...
              IRFAAstr2            , sigg11{mmk}      , sigg22{mmk}   );
           end
        

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%      Sample the acceptor PIFE state trajectory       %%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if Num_states_A>1
              [ tau_A_res{mmk}   ] = sampling_tau_A_HMM( ...
                  ... 
              sigs_size(mmk)       , IRF_D_str2           , IRF_A_str2           , KA_val_learned(mmk,:)' , ...
              F_res{mmk}           , eta_D_learned        , eta_A_learned        , max_bin                , ...
              mu_D_ex_learned{mmk} , mu_A_ex_learned{mmk} , mu_BD_learned(mmk,1) , tau_D_res{mmk}         , ...
              mu_BA_learned(mmk,1) , KD_val_learned(mmk,:), KF_val_learned(mmk,:), AII{mmk}               , ...
              indx1{mmk}           , indx2{mmk}           , sigg1{mmk}           , sigg2{mmk}             , ...
              A_all                , Num_states_A         , IRFDDstr2            , IRFAAstr2              , ...
              sigg11{mmk}          , sigg22{mmk}          );
           end


           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%        Sample the donor PIFE state trajectory         %%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if Num_states_D>1
              [ tau_D_res{mmk}     , pi_learned{mmk} ] = sampling_tau_D_HMM( ...
                 ... 
              sigs_size(mmk)       , IRF_D_str2            , IRF_A_str2           , KA_val_learned(mmk,:) , ...
              F_res{mmk}           , eta_D_learned         , eta_A_learned        , max_bin               , ...
              mu_D_ex_learned{mmk} , mu_A_ex_learned{mmk}  , mu_BD_learned(mmk,1) , tau_A_res{mmk}        , ...
              mu_BA_learned(mmk,1) , KD_val_learned(mmk,:)', KF_val_learned(mmk,:), DII{mmk}              , ...
              indx1{mmk}           , indx2{mmk}            , sigg1{mmk}           , sigg2{mmk}            , ...
              D_all                , Num_states_D          , IRFDDstr2            , IRFAAstr2             , ...
              sigg11{mmk}          , sigg22{mmk}         );
           end
     
           
      
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%       Sample the donor and ccceptor excitation rates         %%%%%%%
           %%%%%%                            and                               %%%%%%%  
           %%%%%%  donor and acceptor channel background photon emission rates %%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           [mu_D_ex_learned{mmk} , mu_A_ex_learned{mmk} , mu_BD_learned(mmk,1) , mu_BA_learned(mmk,1) , ...
            pi_learned{mmk}      , pi_accept(:,mmk)     , log_like(1,mmk)                             ] = sampling_pi_tau( ...
                 ... 
           pi_accept(:,mmk)      , IRF_D_str2           , IRF_A_str2           , F_res{mmk}           , ...
           tau_D_res{mmk}        , KA_val_learned(mmk,:), eta_D_learned        , eta_A_learned        , ...
           max_bin               , mu_D_ex_learned{mmk} , mu_A_ex_learned{mmk} , mu_BD_learned(mmk,1) , ...
           mu_BA_learned(mmk,1)  , pi_learned{mmk}      , prop_mu_D_ex         , alpha_mu_D_ex        , ...
           beta_mu_D_ex          , prop_mu_A_ex         , alpha_mu_A_ex        , beta_mu_A_ex         , ...
           prop_mu_BD            , alpha_mu_BD          , beta_mu_BD           , prop_mu_BA           , ...
           alpha_mu_BA           , beta_mu_BA           , KD_val_learned(mmk,:), KF_val_learned(mmk,:), ...
           indx1{mmk}            , indx2{mmk}           , sigg1{mmk}           , sigg2{mmk}           , ...
           Num_states_D          , Num_states_A         , tau_A_res{mmk}       , IRFDDstr2            , ...
           IRFAAstr2             , sigg11{mmk}          , sigg22{mmk}         );
       end

     
        [log_like , KA_val_learned  , KD_val_learned  , KF_val_learned  ] = sampling_pi_tau1( ...
                 ... 
        IRF_D_str2          , IRF_A_str2       , F_res               , tau_D_res       , ...
        KA_val_learned      , eta_D_learned    , eta_A_learned       , max_bin         , ...
        KD_val_learned      , KF_val_learned   , pi_learned          , indx1           , ...
        indx2               , sigg1            , sigg2               , Num_states_D    , ...
        Num_states_A        , tau_A_res        , KA_val_prop         , KA_val_alpha    , ...
        KA_val_beta         , KF_val_prop      , KF_val_alpha        , KF_val_beta     , ...
        KD_val_prop         , KD_val_alpha     , KD_val_beta         , Num_states_F    , ...
        log_like            , IRFDDstr2        , IRFAAstr2           , sigg11          , ...
        sigg22              );



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%         Sample FRET state transition rates         %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  Num_states_F==1
            KF_log_like = 0;
        else
            [ KF_rate_learned , FII     , KF_log_like ] = sampling_Transit_rates( ...
             ...
            F_res           , KF_rate_alpha  , KF_rate_beta  , sigs_size          , signal_time  , ...
            KF_rate_learned , FII            , Num_states_F  , Data.num_sub_sigs );
         
            F_all  = sampling_rate_all( F_res   , Num_states_F   , F_all_alpha    , F_all_base  , Data.num_sub_sigs );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%        Sample donor PIFE state transition rates        %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  Num_states_D==1
            KD_log_like = 0;
        else
            [ KD_rate_learned , DII      , KD_log_like  ] = sampling_Transit_rates( ...
             ...
            tau_D_res       , KD_rate_alpha  , KD_rate_beta  , sigs_size          , signal_time  , ...
            KD_rate_learned , DII            , Num_states_D  , Data.num_sub_sigs  );
        
            D_all  = sampling_rate_all( tau_D_res , Num_states_D  , D_all_alpha    , D_all_base    , Data.num_sub_sigs );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%       Sample acceptor PIFE state transition rates       %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if   Num_states_A==1
             KA_log_like = 0;
        else
            [ KA_rate_learned , AII        , KA_log_like  ] = sampling_Transit_rates( ...
             ...
             tau_A_res       , KA_rate_alpha  , KA_rate_beta  , sigs_size           , signal_time  , ...
             KA_rate_learned , AII            , Num_states_A  , Data.num_sub_sigs  );
         
             A_all  = sampling_rate_all( tau_A_res  , Num_states_A   , A_all_alpha    , A_all_base     , Data.num_sub_sigs );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%    Sample the donor and ccceptor channel efficiencies     %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Data.sample_eta
          [ eta_D_learned   , eta_A_learned       , Data.eta_accept    , log_like           ] = sampling_eta( ...
                 ... 
           Data.eta_accept  , Data.num_sub_sigs   , IRF_D_str2         , IRF_A_str2         , ...
           F_res            , tau_D_res           , KA_val_learned     , Data.alpha_eta_D   , ...
           Data.beta_eta_D  , Data.eta_prop_D     , Data.eta_prop_A    , Data.alpha_eta_A   , ...
           Data.beta_eta_A  , eta_D_learned       , eta_A_learned      , pi_learned         , ...
           max_bin          , KD_val_learned      , KF_val_learned     , indx1              , ...
           indx2            , sigg1               , sigg2              , log_like           , ...
           tau_A_res        , IRFDDstr2           , IRFAAstr2          , sigg11             , ...
           sigg22          );
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%                     MAP estimate                      %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ map                , Data.map_KF_rate   , Data.map_KD_rate   , Data.map_KA_rate   , ...
          Data.map_KF_val    , Data.map_KD_val    , Data.map_KA_val    , Data.map_F_res     , ...
          Data.map_tau_D     , Data.map_tau_A     , Data.map_eta_D     , Data.map_eta_A     , ...
          Data.map_mu_D_ex   , Data.map_mu_A_ex   , Data.map_mu_BD     , Data.map_mu_BA     , ...
          Data.map_F_all     , Data.map_D_all     , Data.map_A_all     , RAW_map            ] = MAP_estimate( ...
          ...
          map                , KF_rate_learned    , KD_rate_learned    , KA_rate_learned    , ...
          KF_val_learned     , KD_val_learned     , KA_val_learned     , eta_D_learned      , ...
          eta_A_learned      , mu_D_ex_learned    , mu_A_ex_learned    , mu_BD_learned      , ...
          mu_BA_learned      , F_res              , tau_D_res          , tau_A_res          , ...
          log_like           , Data.alpha_mu_D_ex , Data.beta_mu_D_ex  , Data.alpha_mu_A_ex , ...
          Data.beta_mu_A_ex  , Data.alpha_mu_BD   , Data.beta_mu_BD    , Data.alpha_mu_BA   , ...
          Data.beta_mu_BA    , Data.KD_rate_alpha , Data.KD_rate_beta  , Data.KF_rate_alpha , ...
          Data.KF_rate_beta  , KD_log_like        , KF_log_like        , Data.map_KF_rate   , ...
          Data.map_KD_rate   , Data.map_KA_rate   , Data.map_KF_val    , Data.map_KD_val    , ...
          Data.map_KA_val    , Data.map_F_res     , Data.map_tau_D     , Data.map_tau_A     , ...
          Data.map_eta_D     , Data.map_eta_A     , Data.map_mu_D_ex   , Data.map_mu_A_ex   , ...
          F_all_base         , D_all_base         , A_all_base         , Data.map_mu_BD     , ...
          Data.map_mu_BA     , Data.num_sub_sigs  , Data.alpha_eta_D   , Data.beta_eta_D    , ...
          Data.alpha_eta_A   , Data.beta_eta_A    , Data.KA_val_alpha  , Data.KA_val_beta   , ...
          Data.KF_val_alpha  , Data.KF_val_beta   , Data.KD_val_alpha  , Data.KD_val_beta   , ...
          Data.map_F_all     , Data.map_D_all     , Data.map_A_all     , F_all              , ...
          D_all              , A_all              , Num_states_F       , Num_states_D       , ...
          Num_states_A       , F_all_alpha        , D_all_alpha        , A_all_alpha        , ...
          KA_log_like        , Data.KA_rate_alpha , Data.KA_rate_beta );
           
          Data.map     = [Data.map     , map     ];
          Data.RAW_map = [Data.RAW_map , RAW_map ];
          
          
      if (floor(i/Data.save_size)*Data.save_size-i)==0
          for mmk=1:Data.num_sub_sigs
              Data.F_res{mmk}           = [ Data.F_res{mmk}           ; F_res{mmk}             ];
              Data.tau_D_res{mmk}       = [ Data.tau_D_res{mmk}       ; tau_D_res{mmk}         ];
              Data.tau_A_res{mmk}       = [ Data.tau_A_res{mmk}       ; tau_A_res{mmk}         ];
              Data.mu_D_ex_learned{mmk} = [ Data.mu_D_ex_learned{mmk} ; mu_D_ex_learned{mmk}   ];
              Data.mu_A_ex_learned{mmk} = [ Data.mu_A_ex_learned{mmk} ; mu_A_ex_learned{mmk}   ];
              Data.KF_val_learned{mmk}  = [ Data.KF_val_learned{mmk}  ; KF_val_learned(mmk,:)  ];
              Data.KD_val_learned{mmk}  = [ Data.KD_val_learned{mmk}  ; KD_val_learned(mmk,:)  ];
              Data.KA_val_learned{mmk}  = [ Data.KA_val_learned{mmk}  ; KA_val_learned(mmk,:)  ];
          end
          Data.KF_rate_learned   = cat(3,Data.KF_rate_learned , KF_rate_learned  );
          Data.KD_rate_learned   = cat(3,Data.KD_rate_learned , KD_rate_learned  );
          Data.KA_rate_learned   = cat(3,Data.KA_rate_learned , KA_rate_learned  );
          Data.F_all             = [ Data.F_all               , F_all            ];
          Data.D_all             = [ Data.D_all               , D_all            ];
          Data.A_all             = [ Data.A_all               , A_all            ];
          if Data.sample_eta
             Data.eta_A_learned  = [ Data.eta_A_learned       , eta_A_learned    ];
             Data.eta_D_learned  = [ Data.eta_D_learned       , eta_D_learned    ];
          end
          Data.mu_BD_learned     = [ Data.mu_BD_learned       , mu_BD_learned    ];
          Data.mu_BA_learned     = [ Data.mu_BA_learned       , mu_BA_learned    ];
          Data.pi_accept         = pi_accept;
%           Data.tau_accept        = tau_accept;
          
          if  (floor(i/(10*Data.save_size))*(10*Data.save_size)-i)==0
              if  Data.save_on_off
               % Save the data in .mat format
%                save('param','Data','-v7.3','-nocompression')
               save('param','Data','-v7.3')
              end
          end
           
          tester = false;
      else
          tester =true;
      end
          
      
   end
  
   if tester
       for mmk=1:Data.num_sub_sigs
           Data.F_res{mmk}           = [ Data.F_res{mmk}           ; F_res{mmk}             ];
           Data.tau_D_res{mmk}       = [ Data.tau_D_res{mmk}       ; tau_D_res{mmk}         ];
           Data.tau_A_res{mmk}       = [ Data.tau_A_res{mmk}       ; tau_A_res{mmk}         ];
           Data.mu_D_ex_learned{mmk} = [ Data.mu_D_ex_learned{mmk} ; mu_D_ex_learned{mmk}   ];
           Data.mu_A_ex_learned{mmk} = [ Data.mu_A_ex_learned{mmk} ; mu_A_ex_learned{mmk}   ];
           Data.KF_val_learned{mmk}  = [ Data.KF_val_learned{mmk}  ; KF_val_learned(mmk,:)  ];
           Data.KD_val_learned{mmk}  = [ Data.KD_val_learned{mmk}  ; KD_val_learned(mmk,:)  ];
           Data.KA_val_learned{mmk}  = [ Data.KA_val_learned{mmk}  ; KA_val_learned(mmk,:)  ];
       end
       Data.KF_rate_learned   = cat(3,Data.KF_rate_learned , KF_rate_learned  );
       Data.KD_rate_learned   = cat(3,Data.KD_rate_learned , KD_rate_learned  );
       Data.KA_rate_learned   = cat(3,Data.KA_rate_learned , KA_rate_learned  );
       Data.F_all             = [ Data.F_all               , F_all            ];
       Data.D_all             = [ Data.D_all               , D_all            ];
       Data.A_all             = [ Data.A_all               , A_all            ];
       if Data.sample_eta
          Data.eta_A_learned  = [ Data.eta_A_learned       , eta_A_learned    ];
          Data.eta_D_learned  = [ Data.eta_D_learned       , eta_D_learned    ];
       end
       Data.mu_BD_learned     = [ Data.mu_BD_learned       , mu_BD_learned    ];
       Data.mu_BA_learned     = [ Data.mu_BA_learned       , mu_BA_learned    ];
       Data.pi_accept         = pi_accept  ;
%        Data.tau_accept        = tau_accept ;
   end

   if  Data.save_on_off
       % Save the data in .mat format
%        save('param','Data','-v7.3','-nocompression')
       save('param','Data','-v7.3')
   end
end
