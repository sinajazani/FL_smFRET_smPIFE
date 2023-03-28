%% Sampling the state transition rates
%
%  This code is written by Sina Jazani (03/28/2023)
%  Contact: sjazani1@jhu.edu
function[ KFDA_rate_le , FDAKK , KFDA_log_like ] = sampling_Transit_rates( ...
...
FDA_res        , KFDA_rate_alpha    , KFDA_rate_beta     , sigs_size      , signal_time    , ...
KFDA_rate_le     , FDAKK            , Num_states       , num_sub_sigs    )
      

% Gamma +  Dirichlet direct sampling

% Parameters of hyper-prior
KFDA_PI_alp = 2;
KFDA_PI_bet = ones(Num_states-1,1)./(Num_states-1);

matt      = nan(Num_states);
Q         = zeros(Num_states);
TT        = zeros(Num_states,1);
for mmk=1:num_sub_sigs
    for k=1:sigs_size(mmk,1)-1
        Q(FDA_res{mmk}(k+1),FDA_res{mmk}(k))=Q(FDA_res{mmk}(k+1),FDA_res{mmk}(k))+1;
        TT(FDA_res{mmk}(k))=TT(FDA_res{mmk}(k))+signal_time{mmk}(k);
    end
end

indx=1:Num_states;
for mm=1:Num_states
    dd           = (indx~=mm);
    matt(mm,mm)  = -gamrnd(KFDA_rate_alpha+sum(Q(dd,mm)) , 1/(1./KFDA_rate_beta+TT(mm)));
    matt(dd,mm)  = -matt(mm,mm).*dirichletRnd(KFDA_PI_alp*KFDA_PI_bet+Q(dd,mm));
end
matt=matt-diag(diag(matt)); matt=matt-diag(sum(matt,1));

KFDA_rate_le   = matt.*(eye(Num_states)==0);
KFDA_log_like  = 0;

% Calculate the eigenvector and eigenvalue of the transition rate materix
% in order to calculate the matrix expoential in a faster way
[vecP , valP]=eig(matt); vallP = diag(valP); vec1P=vecP^-1;

for mmk=1:num_sub_sigs
    for k=1:sigs_size(mmk,1)-1
            FDAKK{mmk}(:,:,k)   = real(vecP*diag(exp(vallP*signal_time{mmk}(k)))*(vec1P));
            KFDA_log_like       = KFDA_log_like +log(FDAKK{mmk}(FDA_res{mmk}(k+1) , FDA_res{mmk}(k) , k) );
    end
end

end