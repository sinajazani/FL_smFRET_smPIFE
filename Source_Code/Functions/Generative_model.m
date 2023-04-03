function[ Data ]   = Generative_model(Data)


% Erase the Data vectors just in case! 
Data.signal             = []; Data.PI_real            = [];
Data.F_fluc_real        = []; Data.tau_D_fluc_real    = []; Data.tau_A_fluc_real    = [];     

% Pre localize the space
Num_states_F            = length(Data.KF_real);
Num_states_D            = length(Data.KD_real);
Num_states_A            = length(Data.KA_real);

% Transition probabilities of FRET and PIFE
if  sum(diag(Data.KF_real)+sum(Data.KF_real-diag(diag(Data.KF_real)),1)')==0
    TF                  = expm((Data.KF_real)*Data.pulse_priod*10^-9);
else
    disp('The input value is not correct. We normalize the FRET transition matrix!')
    TF                  = expm((Data.KF_real-diag(diag(Data.KF_real))-diag(sum(Data.KF_real-diag(diag(Data.KF_real)),1)))*Data.pulse_priod*10^-9);
end

if  sum(diag(Data.KD_real)+sum(Data.KD_real-diag(diag(Data.KD_real)),1)')==0
    TD                  = expm((Data.KD_real)*Data.pulse_priod*10^-9);
else
    disp('The input value is not correct. We normalize the donor-PIFE transition matrix!')
    TD                  = expm((Data.KD_real-diag(diag(Data.KD_real))-diag(sum(Data.KD_real-diag(diag(Data.KD_real)),1)))*Data.pulse_priod*10^-9);
end

if  sum(diag(Data.KA_real)+sum(Data.KA_real-diag(diag(Data.KA_real)),1)')==0
    TA                  = expm((Data.KA_real)*Data.pulse_priod*10^-9);
else
    disp('The input value is not correct. We normalize the acceptor-PIFE transition matrix!')
    TA                  = expm((Data.KA_real-diag(diag(Data.KA_real))-diag(sum(Data.KA_real-diag(diag(Data.KA_real)),1)))*Data.pulse_priod*10^-9);
end

F_state(1)              = Discrete_sampler(ones(1,Num_states_F)./Num_states_F);
D_state(1)              = Discrete_sampler(ones(1,Num_states_D)./Num_states_D);
A_state(1)              = Discrete_sampler(ones(1,Num_states_A)./Num_states_A);

for ii=2:Data.Num_pulses
    F_state(ii,1)       = find( rand(1)<=cumsum(TF(:,F_state(ii-1,1))), 1 )                  ; 
    D_state(ii,1)       = find( rand(1)<=cumsum(TD(:,D_state(ii-1,1))), 1 )                  ;
    A_state(ii,1)       = find( rand(1)<=cumsum(TA(:,A_state(ii-1,1))), 1 )                  ;
end

F_fluc_real             = Data.KF_val_real(1,F_state)                     ;
tau_D_fluc_real         = Data.KD_val_real(1,D_state)                     ;
tau_A_fluc_real         = Data.KA_val_real(1,A_state)                     ;

% Calculate the weights over each pulse. This weight is not normalized and sum of it is equal to the photon detection probability for each pulse.
pii1(1,:)               = exp(-Data.ex_D_real(1,D_state));
pii2(1,:)               = exp(-Data.ex_A_real(1,A_state));
pii3                    = exp(-Data.pulse_priod*Data.Back_D_real.*(10^-9));
pii4                    = exp(-Data.pulse_priod*Data.Back_A_real.*(10^-9));
PII                     = [ (1-pii1).*pii2.*pii3*pii4 ; ...
                            (1-pii2).*pii1.*pii3*pii4 ; ...
                            (1-pii3).*pii1.*pii2*pii4 ; ...
                            (1-pii4).*pii1.*pii2*pii3 ];  

% Sample the events of photon detetction for each excitation pulse
Active_pulses           = sum(PII,1)>rand(1,Data.Num_pulses);

% Calculate the weight on the source of the detected photons
Data.PI_real            = PII(:,Active_pulses)./sum(PII(:,Active_pulses),1);

% Number of detcted photons
num_photons              = size(Data.PI_real,2);

% Sample the source of the photons
PPo                     = cumsum(Data.PI_real,1);
Events                  = 1+sum((repmat(rand(1,num_photons),4,1) <= PPo)==0,1);

% Save the FRET and PIFE states of the molecule at the time we detect the photons
Data.F_fluc_real        = F_fluc_real(1,Active_pulses)     ;
Data.tau_D_fluc_real    = tau_D_fluc_real(1,Active_pulses) ;
Data.tau_A_fluc_real    = tau_A_fluc_real(1,Active_pulses) ;

% Pre caculation of the KT
K_T                     = (Data.F_fluc_real)+(1./Data.tau_D_fluc_real);

% Order of pulses
pulses                  = 1:Data.Num_pulses;

% The time photons are detected
ii                      = pulses(Active_pulses)*Data.pulse_priod*10^-9;

for i=1:num_photons
    
    switch Events(i)
        case 1 % Donor excitation
             % FRET events 
             FRET = Discrete_sampler([1./Data.tau_D_fluc_real(i),Data.F_fluc_real(i)]);  
             
             if  FRET==1  % Donor excitation-> Donor emission
                 if   Data.eta_D_real>rand() % Acceptor detection
                      Data.signal(:,i)  = [(exprnd(1/K_T(i))+Data.IRF_mean(2)+Data.IRF_str(2)*randn());ii(i);2];
                 else                        % Donor detection
                      Data.signal(:,i)  = [(exprnd(1/K_T(i))+Data.IRF_mean(1)+Data.IRF_str(1)*randn());ii(i);1];
                 end
             else     % Donor excitaiton -> FRET -> Acceptor emission
                 if   Data.eta_A_real>rand() % Donor detection
                      Data.signal(:,i)  = [(exprnd(1/K_T(i))+exprnd(Data.tau_A_fluc_real(i))+Data.IRF_mean(1)+Data.IRF_str(1)*randn());ii(i);1];
                 else                        % Acceptor detection
                      Data.signal(:,i)  = [(exprnd(1/K_T(i))+exprnd(Data.tau_A_fluc_real(i))+Data.IRF_mean(2)+Data.IRF_str(2)*randn());ii(i);2];
                 end
             end
        case 2  % Acceptor excitation
            if   Data.eta_A_real>rand() % Donor detection
                 Data.signal(:,i)  = [(exprnd(Data.tau_A_fluc_real(i))+Data.IRF_mean(1)+Data.IRF_str(1)*randn());ii(i);1];
            else                        % Acceptor detection
                 Data.signal(:,i)  = [(exprnd(Data.tau_A_fluc_real(i))+Data.IRF_mean(2)+Data.IRF_str(2)*randn());ii(i);2];
            end
        case 3   % Donor background
            if   Data.eta_D_real>rand() % Acceptor detection
                 Data.signal(:,i)  = [Data.pulse_priod*rand();ii(i);2];
            else                        % Donor detection
                 Data.signal(:,i)  = [Data.pulse_priod*rand();ii(i);1];
            end
        case 4   % Acceptor background
            if   Data.eta_A_real>rand() % Donor detection
                 Data.signal(:,i)  = [Data.pulse_priod*rand();ii(i);1];
            else                        % Acceptor detection
                 Data.signal(:,i)  = [Data.pulse_priod*rand();ii(i);2];
            end
    end
end

end