%% Sampling FRET state trajectory through time
%  In this function we apply a forward filtering and backward sampling approach
%
%  This code is written by Sina Jazani (03/28/2023)
%  Contact: sjazani1@jhu.edu
function[ F_K ] = sampling_FRET_HMM(...
    ...
sigs_size   , IRF_D_str2    , IRF_A_str2   , FII       , ...
tau_A       , tau_D         , eta_D        , KF_val    , ...
eta_A       , PI            , T_pulse      , KP_val    , ...
indx1       , indx2         , sigg1        , sigg2     , ...
F11         , Num_states    , A_res        , IRFDDstr2 , ...
IRFAAstr2   , sigg11        , sigg22       )
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % sigs_size             :  The length of the trajectory, equivalent to number of detetcted photons
   % IRF_D_str_sqrt2       :  A pre-calculated value for the likelihood. (sqrt(2*IRF_D_str^2))
   % IRF_D_str2            :  A pre-calculated value for the likelihood. IRF_D_str^2
   % IRF_A_str_sqrt2       :  A pre-calculated value for the likelihood. (sqrt(2*IRF_A_str^2))
   % IRF_A_str2            :  A pre-calculated value for the likelihood. IRF_A_str^2
   % FII                   :  Transition probability matrixes
   % tau_A                 :  Acceptor lifetime for each donor state
   % tau_D                 :  Donor lifetime state trajectory
   % eta_D                 :  Cross-talk ratio from donor channel
   % eta_A                 :  Cross-talk ratio from acceptor channel
   % KF_val                :  FRET rate values for each FRET state
   % PI                    :  Weights on the source of photon detections (Donor , Acceptor , Donor backgound , Acceptor background)
   % T_pulse               :  Inter-pulse time
   % KD_val                :  Donor lifetime for each donor state
   % indx1                 :  Index over photons collected in the donor channel
   % indx2                 :  Index over photons collected in the acceptor channel
   % sigg1                 :  Pre-calculated values for the likelihood
   % sigg2                 :  Pre-calculated values for the likelihood
   % F11                   :  Weights on the initial FRET states
   % Num_states            :  Nember of FRET states
   % A_res                 :  Acceptor lifetime state trajectory
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % F_K                   :  The sampled FRET state trajectory
   
   
% Assign the space for the FRET state trajectory and the filter
F_K                 = nan(1         ,sigs_size);
A                   = nan(Num_states,sigs_size);
likelihood          = nan(Num_states,sigs_size);


% Pre-calculations
tau_A1              = tau_A(1,A_res(indx1))                               ;
tau_A2              = tau_A(1,A_res(indx2))                               ;

KT1                 = 1./KP_val(1,tau_D(indx1))+KF_val                    ;
KT2                 = 1./KP_val(1,tau_D(indx2))+KF_val                    ;

% AA1                 = exp(KT1.*(0.5*IRF_D_str2*KT1+sigg1)+log(erfc(IRFDDstr2*KT1+sigg11)));
% AA2                 = exp(KT2.*(0.5*IRF_A_str2*KT2+sigg2)+log(erfc(IRFAAstr2*KT2+sigg22)));

AA1                 = exp(KT1.*(0.5*IRF_D_str2*KT1+sigg1)+log_erfc(IRFDDstr2*KT1+sigg11));
AA2                 = exp(KT2.*(0.5*IRF_A_str2*KT2+sigg2)+log_erfc(IRFAAstr2*KT2+sigg22));

DD1                 = exp((0.5*IRF_D_str2./tau_A1+sigg1)./tau_A1).*erfc(IRFDDstr2./tau_A1+sigg11);
DD2                 = exp((0.5*IRF_A_str2./tau_A2+sigg2)./tau_A2).*erfc(IRFAAstr2./tau_A2+sigg22);

% Calculate the likelihood for all FRET states at any detectted photons in both channels 
likelihood(:,indx1) = (((PI(1,indx1).*AA1./KP_val(1,tau_D(indx1)))+(PI(3,indx1)./T_pulse)).*(1-eta_D))+...
                      (((PI(1,indx1).*KF_val.*((DD1-AA1)./(KT1.*tau_A1-1)))+(PI(2,indx1).*DD1./tau_A1)+(PI(4,indx1)./T_pulse)).*eta_A) ;
               
likelihood(:,indx2) = (((PI(1,indx2).*AA2./KP_val(1,tau_D(indx2)))+(PI(3,indx2)./T_pulse)).*eta_D) +...
                      (((PI(1,indx2).*KF_val.*((DD2-AA2)./(KT2.*tau_A2-1)))+(PI(2,indx2).*DD2./tau_A2)+(PI(4,indx2)./T_pulse)).*(1-eta_A)) ;

% Calculate the first filter effected by the prior on the first FRET state
A(:,1)              = F11.*likelihood(:,1)                                ;
A(:,1)              = A(:,1)./sum(A(:,1))                                 ;

% Construct the computational stable forward filters
for n = 1:(sigs_size-1)
    A(:,n+1)        = likelihood(:,n+1).*(FII(:,:,n)*A(:,n))              ;
    A(:,n+1)        = A(:,n+1)./sum(A(:,n+1))                             ;
end

% Sample the last FRET state directly from the last filter
P                   = cumsum(A(:,sigs_size))                              ;
F_K(1,sigs_size)    = find( P(end)*rand(1)<=P, 1 )                        ;

% Sample the rest of the FRET states by marching backward
for n = (sigs_size-1):-1:1
    P               = cumsum(FII(F_K(1,n+1),:,n)'.*A(:,n))                ;
    F_K(1,n)        = find( P(end)*rand(1)<=P, 1 )                        ;
end


end
% End of the function