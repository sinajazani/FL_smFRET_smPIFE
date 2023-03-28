%% Calculate the Maximum of posteriori (MAP) and corresponding values of the parameters of interests
%
%  This code is written by Sina Jazani (09/14/2021)
%  Contact: sjazani1@jhu.edu
function[ map             , map_KF_rate       , map_KD_rate       , map_KA_rate      , ...
          map_KF_val      , map_KD_val        , map_KA_val        , map_F_res        , ...
          map_tau_D_res   , map_tau_A_res     , map_eta_D         , map_eta_A        , ...
          map_mu_D_ex     , map_mu_A_ex       , map_mu_BD         , map_mu_BA        , ...
          map_F_all       , map_D_all         , map_A_all         , RAW_map          ] = MAP_estimate( ...
          ...
          map             , KF_rate_learned   , KD_rate_learned   , KA_rate_learned  , ...
          KF_val_learned  , KD_val_learned    , KA_val_learned    , eta_D_learned    , ...
          eta_A_learned   , mu_D_ex_learned   , mu_A_ex_learned   , mu_BD_learned    , ...
          mu_BA_learned   , F_res             , tau_D_res         , tau_A_res        , ...
          log_like        , alpha_mu_D_ex     , beta_mu_D_ex      , alpha_mu_A_ex    , ...
          beta_mu_A_ex    , alpha_mu_BD       , beta_mu_BD        , alpha_mu_BA      , ...
          beta_mu_BA      , KD_rate_alpha     , KD_rate_beta      , KF_rate_alpha    , ...
          KF_rate_beta    , KD_log_like       , KF_log_like       , map_KF_rate      , ...
          map_KD_rate     , map_KA_rate       , map_KF_val        , map_KD_val       , ...
          map_KA_val      , map_F_res         , map_tau_D_res     , map_tau_A_res    , ...
          map_eta_D       , map_eta_A         , map_mu_D_ex       , map_mu_A_ex      , ...
          F_all_base      , D_all_base        , A_all_base        , map_mu_BD        , ...
          map_mu_BA       , num_sub_sigs      , alpha_eta_D       , beta_eta_D       , ...
          alpha_eta_A     , beta_eta_A        , KA_val_alpha      , KA_val_beta      , ...
          KF_val_alpha    , KF_val_beta       , KD_val_alpha      , KD_val_beta      , ...
          map_F_all       , map_D_all         , map_A_all         , F_all            , ...
          D_all           , A_all             , Num_states_F      , Num_states_D     , ...
          Num_states_A    , F_all_alpha       , D_all_alpha       , A_all_alpha      , ...
          KA_log_like     , KA_rate_alpha     , KA_rate_beta       )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % sigs_size             :  The length of the trajectory, equivalent to number of detetcted photons
   % IRF_D_str_sqrt2       :  A pre-calculated value for the likelihood. (sqrt(2*IRF_D_str^2))
   % IRF_D_str2            :  A pre-calculated value for the likelihood. IRF_D_str^2
   % IRF_A_str_sqrt2       :  A pre-calculated value for the likelihood. (sqrt(2*IRF_A_str^2))
   % IRF_A_str2            :  A pre-calculated value for the likelihood. IRF_A_str^2
   % tau_A                 :  Acceptor lifetimes for any detected photon
   % tau_D_res             :  Inverse of donor lifetimes for any detected photon
   % F_res                 :  FRET rates for each state
   % eta_D                 :  Cross-talk ratio from donor channel
   % eta_A                 :  Cross-talk ratio from acceptor channel
   % PI                    :  Weights on the source of photon detections (Donor , Acceptor , Donor backgound , Acceptor background)
   % T_pulse               :  Inter-pulse time
   % indx1                 :  Index over photons collected in the donor channel
   % indx2                 :  Index over photons collected in the acceptor channel
   % sigg1                 :  Pre-calculated values for the likelihood
   % sigg2                 :  Pre-calculated values for the likelihood
   % Num_states            :  Nember of FRET states
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % map                   :  Value of the MAP estimate
   % KF_rate_learned       :  The likelihood for any detected photon
   % KD_rate_learned       :  The likelihood for any detected photon
   % KF_val_learned        :  The likelihood for any detected photon
   % F_res                 :  The likelihood for any detected photon
   % tau_D_res             :  The likelihood for any detected photon      
   % KD_val_learned        :  The likelihood for any detected photon
   % KA_val_learned        :  The likelihood for any detected photon  
   % eta_D_learned         :  The likelihood for any detected photon  
   % eta_A_learned         :  The likelihood for any detected photon  
   % mu_D_ex_learned       :  The likelihood for any detected photon
   % mu_A_ex_learned       :  The likelihood for any detected photon
   % mu_BD_learned         :  The likelihood for any detected photon 
   % mu_BA_learned         :  The likelihood for any detected photon 
   % RAW_map               :  The likelihood for any detected photon  
   % map_F_all             :  The likelihood for any detected photon 
   % map_D_all             :  The likelihood for any detected photon 
   % map_pi_F_rate         :  The likelihood for any detected photon
   % map_pi_D_rate         :  The likelihood for any detected photon
   % map_pi_A_rate         :  The likelihood for any detected photon
   % map_tau_A             :  The likelihood for any detected photon
   % map_A_all             :  The likelihood for any detected photon
   % map_KA_rate           :  The likelihood for any detected photon
   % tau_A_res             :  The likelihood for any detected photon
   % KA_rate_learned       :  The likelihood for any detected photon
      
      
      
% Calculating the joint posteior probability basssed on the sampled values of parameters of interest


map_new=0;ssA=zeros(Num_states_A,1);ssP=zeros(Num_states_D,1);ssF=zeros(Num_states_F,1);


% Prior on excitation rates for all sub traces
for mmk=1:num_sub_sigs
    map_new =map_new + sum(-alpha_mu_D_ex*log(beta_mu_D_ex)+gammaln(alpha_mu_D_ex) +(alpha_mu_D_ex-1)*log(mu_D_ex_learned{mmk})-mu_D_ex_learned{mmk}/beta_mu_D_ex);
    map_new =map_new + sum(-alpha_mu_A_ex*log(beta_mu_A_ex)+gammaln(alpha_mu_A_ex) +(alpha_mu_A_ex-1)*log(mu_A_ex_learned{mmk})-mu_A_ex_learned{mmk}/beta_mu_A_ex);
    map_new =map_new + sum(-KF_val_alpha*log(KF_val_beta) +gammaln(KF_val_alpha) +(KF_val_alpha-1)*log(KF_val_learned(mmk,:))-KF_val_learned(mmk,:)/KF_val_beta,'all');
    map_new =map_new + sum(-KD_val_alpha*log(KD_val_beta) +gammaln(KD_val_alpha) +(KD_val_alpha-1)*log(KD_val_learned(mmk,:))-KD_val_learned(mmk,:)/KD_val_beta,'all');
    map_new =map_new + sum(-KA_val_alpha*log(KA_val_beta) +gammaln(KA_val_alpha) +(KA_val_alpha-1)*log(KA_val_learned(mmk,:))-KA_val_learned(mmk,:)/KA_val_beta,'all');
    ssA = ssA + histcounts(tau_A_res{mmk}(1),1:Num_states_A+1)';
    ssP = ssP + histcounts(tau_D_res{mmk}(1),1:Num_states_D+1)';
    ssF = ssF + histcounts(F_res{mmk}(1),1:Num_states_F+1)';
end

% Prior on beta_{hyper,D}
map_new = map_new+gammaln(sum(D_all_alpha*D_all_base+ssP))+sum((D_all_alpha*D_all_base+ssP-1).*log(D_all+eps)-gammaln(D_all_alpha*D_all_base+ssP+eps));

% Prior on beta_{hyper,F}
map_new = map_new+gammaln(sum(F_all_alpha*F_all_base+ssF))+sum((F_all_alpha*F_all_base+ssF-1).*log(F_all+eps)-gammaln(F_all_alpha*F_all_base+ssF+eps)); 

% Prior on beta_{hyper,A}
map_new = map_new+gammaln(sum(A_all_alpha*A_all_base+ssA))+sum((A_all_alpha*A_all_base+ssA-1).*log(A_all+eps)-gammaln(A_all_alpha*A_all_base+ssA+eps)); 


map_new =map_new+ sum(log_like) + KD_log_like  + KF_log_like + KA_log_like+...
      sum( -KF_rate_alpha*log(KF_rate_beta) +gammaln(KF_rate_alpha) +(KF_rate_alpha-1)*log(KF_rate_learned)-KF_rate_learned/KF_rate_beta,'all','omitnan')+...
      sum( -KD_rate_alpha*log(KD_rate_beta) +gammaln(KD_rate_alpha) +(KD_rate_alpha-1)*log(KD_rate_learned)-KD_rate_learned/KD_rate_beta,'all','omitnan')+...
      sum( -KA_rate_alpha*log(KA_rate_beta) +gammaln(KA_rate_alpha) +(KA_rate_alpha-1)*log(KA_rate_learned)-KA_rate_learned/KA_rate_beta,'all','omitnan')+...
     +sum(-alpha_mu_BD*log(beta_mu_BD)+gammaln(alpha_mu_BD) +(alpha_mu_BD-1)*log(mu_BD_learned)-mu_BD_learned/beta_mu_BD...
           -alpha_mu_BA*log(beta_mu_BA)+gammaln(alpha_mu_BA) +(alpha_mu_BA-1)*log(mu_BA_learned)-mu_BA_learned/beta_mu_BA)...
          +gammaln(alpha_eta_D+beta_eta_D)-gammaln(alpha_eta_D)-gammaln(beta_eta_D)+(alpha_eta_D-1)*log(eta_D_learned)+(beta_eta_D-1)*log(1-eta_D_learned)...
          +gammaln(alpha_eta_A+beta_eta_A)-gammaln(alpha_eta_A)-gammaln(beta_eta_A)+(alpha_eta_A-1)*log(eta_A_learned)+(beta_eta_A-1)*log(1-eta_A_learned);
       

% In the case of numercal under fellow, model might return imaginary number, and we weill stop the algorithm for further inverstigation
if imag(map_new)>0
   keyboard
end
      
% Save the value of the joint posterior 
RAW_map = map_new;
      
% Update the MAP estimations based on the values of the joint posterior
if  map_new>map
    map             = map_new           ; map_KF_val      = KF_val_learned    ; map_KD_val      = KD_val_learned    ;
    map_KA_val      = KA_val_learned    ; map_F_res       = F_res             ; map_tau_D_res   = tau_D_res         ;
    map_tau_A_res   = tau_A_res         ; map_eta_D       = eta_D_learned     ; map_eta_A       = eta_A_learned     ;
    map_mu_D_ex     = mu_D_ex_learned   ; map_mu_A_ex     = mu_A_ex_learned   ; map_mu_BD       = mu_BD_learned     ;
    map_mu_BA       = mu_BA_learned     ; map_KF_rate     = KF_rate_learned   ; map_KD_rate     = KD_rate_learned   ;
    map_KA_rate     = KA_rate_learned   ; map_F_all       = F_all             ; map_D_all       = D_all             ;
    map_A_all       = A_all             ;
      
end

end