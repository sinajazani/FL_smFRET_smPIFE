%% Sampling acceptor lifetime state trajectory through time
%  In this function we apply a forward filtering and backward sampling approach
%
%  This code is written by Sina Jazani (09/14/2021)
%  Contact: sjazani1@jhu.edu
function[ A_K  ] = sampling_tau_A_HMM(...
    ...
sigs_size  , IRF_D_str2    , IRF_A_str2 , tau_A      , ...
F_K        , eta_D         , eta_A      , T_pulse    , ...
mu_D_ex    , mu_A_ex       , mu_BD      , D_K        , ...
mu_BA      , KD_val        , KF_val     , AII        , ...
indx1      , indx2         , sigg1      , sigg2      , ...
A11        , Num_states    , IRFDDstr2  , IRFAAstr2  , ...
sigg11     , sigg22        )
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % sigs_size             :  The length of the trajectory, equivalent to number of detetcted photons
   % IRF_D_str_sqrt2       :  A pre-calculated value for the likelihood. (sqrt(2*IRF_D_str^2))
   % IRF_D_str2            :  A pre-calculated value for the likelihood. IRF_D_str^2
   % IRF_A_str_sqrt2       :  A pre-calculated value for the likelihood. (sqrt(2*IRF_A_str^2))
   % IRF_A_str2            :  A pre-calculated value for the likelihood. IRF_A_str^2
   % AII                   :  Transition probability matrixes
   % tau_A                 :  Acceptor lifetime for each donor state
   % F_K                   :  FRET state trajectory
   % D_K                   :  Donor-PIFE state trajectory
   % eta_D                 :  Cross-talk ratio from donor channel
   % eta_A                 :  Cross-talk ratio from acceptor channel
   % KF_val                :  FRET rate values for each FRET state
   % T_pulse               :  Inter-pulse time
   % KP_val                :  Donor lifetime for each donor state
   % indx1                 :  Index over photons collected in the donor channel
   % indx2                 :  Index over photons collected in the acceptor channel
   % sigg1                 :  Pre-calculated values for the likelihood
   % sigg2                 :  Pre-calculated values for the likelihood
   % A11                   :  Weights on the initial acceptor lifetime states
   % Num_states            :  Nember of acceptor lifetime states
   % A_K                   :  Acceptor lifetime state trajectory
   % mu_D_ex               :  Donor excitation rates for each state
   % mu_A_ex               :  Acceptor excitation rates for each state
   % mu_BD                 :  Donor background photon emission rate
   % mu_BA                 :  Acceptor background photon emission rate
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % A_K                   :  The sampled acceptor lifetime state trajectory

   
% Assign the space for the acceptor lifetime state trajectory, filter, weights on photons and the likelihood
A_K                 = nan(1          , sigs_size )                        ;
A                   = nan(Num_states , sigs_size )                        ;
PI                  = nan(Num_states , sigs_size , 4)                     ;
likelihood          = nan(Num_states , sigs_size )                        ;

% Calculate the weights on the source of photon detections (Donor , Acceptor , Donor backgound , Acceptor background)
A1                  = exp(-mu_D_ex(1,D_K))                                ;
A2                  = exp(-mu_A_ex)'                                      ;
A3                  = exp(-mu_BD)                                         ;
A4                  = exp(-mu_BA)                                         ;

PI(:,:,1)           = (1-A1).*(A2.*A3.*A4)                                ;
PI(:,:,2)           = ((1-A2).*A3.*A4).*A1                                ;
PI(:,:,3)           = ((1-A3).*A4.*A2).*A1                                ;
PI(:,:,4)           = ((1-A4).*A3.*A2).*A1                                ;
PI                  = PI./sum(PI,3)                                       ;

% Pre-calculations
KT1                 = (1./KD_val(1,D_K(indx1)))+KF_val(1,F_K(indx1))      ;
KT2                 = (1./KD_val(1,D_K(indx2)))+KF_val(1,F_K(indx2))      ;

% AA1                 = exp(KT1.*(0.5*IRF_D_str2*KT1+sigg1)+log(erfc(IRFDDstr2*KT1+sigg11)));
% AA2                 = exp(KT2.*(0.5*IRF_A_str2*KT2+sigg2)+log(erfc(IRFAAstr2*KT2+sigg22)));

AA1                 = exp(KT1.*(0.5*IRF_D_str2*KT1+sigg1)+log_erfc(IRFDDstr2*KT1+sigg11));
AA2                 = exp(KT2.*(0.5*IRF_A_str2*KT2+sigg2)+log_erfc(IRFAAstr2*KT2+sigg22));


DD1                 = exp((0.5*IRF_D_str2./tau_A+sigg1)./tau_A+log(erfc(IRFDDstr2./tau_A+sigg11)));
DD2                 = exp((0.5*IRF_A_str2./tau_A+sigg2)./tau_A+log(erfc(IRFAAstr2./tau_A+sigg22)));


% Calculate the likelihood for all acceptor lifetime states at any detectted photons in both channels 
likelihood(:,indx1) = (((PI(:,indx1,1).*AA1./KD_val(1,D_K(indx1)))+(PI(:,indx1,3)/T_pulse)).*(1-eta_D))+...
                      (((PI(:,indx1,1).*KF_val(1,F_K(indx1)).*((DD1-AA1)./(KT1.*tau_A-1)))+(PI(:,indx1,2).*DD1./tau_A)+(PI(:,indx1,4)/T_pulse)).*eta_A) ;
                
likelihood(:,indx2) = (((PI(:,indx2,1).*AA2./KD_val(1,D_K(indx2)))+(PI(:,indx2,3)/T_pulse)).*eta_D)+...
                      (((PI(:,indx2,1).*KF_val(1,F_K(indx2)).*((DD2-AA2)./(KT2.*tau_A-1)))+(PI(:,indx2,2).*DD2./tau_A)+(PI(:,indx2,4)/T_pulse)).*(1-eta_A)) ;

% Calculate the first filter effected by the prior on the first acceptor lifetime state
A(:,1)              = A11.*likelihood(:,1)                                ;
A(:,1)              = A(:,1)./sum(A(:,1))                                 ;

% Construct the computational stable forward filters
for n = 1:(sigs_size-1)
    A(:,n+1)        = likelihood(:,n+1).*(AII(:,:,n)*A(:,n))              ;
    A(:,n+1)        = A(:,n+1)./sum(A(:,n+1))                             ;
end

% Sample the last acceptor lifetime state directly from the last filter
P                   = cumsum(A(:,sigs_size))                              ;
A_K(1,sigs_size)    = find( P(end)*rand(1)<=P, 1 )                        ;

% Sample the rest of the acceptor lifetime states by marching backward
for n = (sigs_size-1):-1:1
    P               = cumsum(AII(A_K(1,n+1),:,n)'.*A(:,n))                ;
    A_K(1,n)        = find( P(end)*rand(1) <= P, 1 )                      ;
end

end