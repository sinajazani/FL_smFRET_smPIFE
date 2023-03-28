%% Direct sampling of the transition rate ratio to other states
%
%  This code is written by Sina Jazani (03/28/2023)
%  Contact: sjazani1@jhu.edu
function[ FDA_all ] = sampling_rate_all( FDA_res      , Num_states     , all_alpha    , all_base      ,num_sub_sigs )
      
ss=zeros(Num_states,1);
for ii=1:num_sub_sigs
    ss = ss + histcounts(FDA_res{ii}(1),1:Num_states+1)';
end

FDA_all = dirichletRnd(ss+all_alpha*all_base);

end