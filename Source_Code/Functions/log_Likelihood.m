%% Pre-calculation of the log likelihood
%
%  This code is written by Sina Jazani (03/28/2023)
%  Contact: sjazani1@jhu.edu
function[ log_lik ]= log_Likelihood(...
    ...
F_res       , tau_D_res    , tau_A        , eta_D       , ...
eta_A       , IRF_D_str2   , IRF_A_str2   , PI          , ...
T_pulse     , Indx1        , Indx2        , sigg1       , ...
sigg2       , IRFDDstr2    , IRFAAstr2    , sigg11      , ...
sigg22      )


KT1     = tau_D_res(1,Indx1)+F_res(1,Indx1) ;
tau_A1  = tau_A(1,Indx1);

KT2     = tau_D_res(1,Indx2)+F_res(1,Indx2) ;
tau_A2  = tau_A(1,Indx2);


% Original calculation
% AA1     = exp(KT1.*(0.5*IRF_D_str2.*KT1+sigg1)).*erfc(IRFDDstr2*KT1+sigg11);
% AA2     = exp(KT2.*(0.5*IRF_A_str2.*KT2+sigg2)).*erfc(IRFAAstr2*KT2+sigg22);


% Approximation of log(erfc) for large values>27, it is 6 times slower
AA1     = exp(KT1.*(0.5*IRF_D_str2.*KT1+sigg1)+log_erfc(IRFDDstr2*KT1+sigg11));
AA2     = exp(KT2.*(0.5*IRF_A_str2.*KT2+sigg2)+log_erfc(IRFAAstr2*KT2+sigg22));



% AA11     = exp(KT1.*(0.5*IRF_D_str2.*KT1+sigg1)+double(log(erfc(vpa(IRFDDstr2*KT1+sigg11)))));
% AA22     = exp(KT2.*(0.5*IRF_A_str2.*KT2+sigg2)+double(log(erfc(vpa(IRFAAstr2*KT2+sigg22)))));

% In the case very large KF (>500), the scaled erfc might be helpful to avoid numerical error. 
% AA1     = exp(KT1.*(0.5*IRF_D_str2.*KT1+sigg1)-(IRFDDstr2*KT1+sigg11).^2+log(erfcx(IRFDDstr2*KT1+sigg11)));
% AA2     = exp(KT2.*(0.5*IRF_A_str2.*KT2+sigg2)-(IRFAAstr2*KT2+sigg22).^2+log(erfcx(IRFAAstr2*KT2+sigg22)));
% 
DD1     = exp((0.5*IRF_D_str2./tau_A1+sigg1)./tau_A1).*erfc(IRFDDstr2./tau_A1+sigg11);
DD2     = exp((0.5*IRF_A_str2./tau_A2+sigg2)./tau_A2).*erfc(IRFAAstr2./tau_A2+sigg22);

log_lik = sum(log( (((PI(1,Indx1).*AA1.*tau_D_res(1,Indx1))+(PI(3,Indx1)/T_pulse)).*(1-eta_D))+...
                   (((PI(1,Indx1).*F_res(1,Indx1).*(DD1-AA1)./(KT1.*tau_A1-1))+(PI(2,Indx1).*DD1./tau_A1)+(PI(4,Indx1)/T_pulse)).*eta_A) +eps) )+...
...
          sum(log( (((PI(1,Indx2).*AA2.*tau_D_res(1,Indx2))+(PI(3,Indx2)/T_pulse)).*eta_D)+...
                   (((PI(1,Indx2).*F_res(1,Indx2).*(DD2-AA2)./(KT2.*tau_A2-1))+(PI(2,Indx2).*DD2./tau_A2)+(PI(4,Indx2)/T_pulse)).*(1-eta_A)) +eps) );

end