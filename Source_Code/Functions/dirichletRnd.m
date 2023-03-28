% Sampling from a dirichlet distribution
% written by Sina Jazani (09/27/2021)
%
% email: (sina Jazani) sjazani1@jh.edu

function x = dirichletRnd(a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % a                             :  An array of weights
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % x                             :  Sampled value
     
     
% take a sample from a dirichlet distribution
p = length(a);
x = gamrnd(a,1,p,1);
x = x ./ sum(x);

% x = randg(a);
% x = x./sum(x,2);
% 
% % check for underflows and turn to stick-breakings if so
% ind = isnan(x(:,1));
% if any(ind)
%     disp('Dir underflows')
%      keyboard
%     K = size(a,1);
%     x(:,1) = zeros(K,1);
%     for k=1:K
%         x(k,1) = betarnd(a(k,1),sum(a(k+1:K,1),1)+realmin).*(1-sum(x(:,1),1));
%     end
%      
% end

end

