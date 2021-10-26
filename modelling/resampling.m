function p_value = resampling(data)

%% ERP modelling code relative to:
% Mancini F, Pepe A, Bernacchia, A, Di Stefano G, Mouraux A, Iannetti GD. (2017)
% Characterising the short-term habituation of event-related evoked potentials
% eNeuro

% written in Matlab R2016b by F Mancini, fm456@cam.ac.uk

nperms = 1000; % number of resampling iterations
for perms = 1:nperms
    
    tmp = shuffle(data);
    [fit_tmp, gof_tmp, bic_tmp(:,:,perms)] = postsvd_fitting(tmp);
    
    for rank_ord = 1:60 
        winning_mod_bic_perms(perms,rank_ord) = find(bic_tmp(:,rank_ord,perms)==min(squeeze((bic_tmp(:,rank_ord,perms)))));
        hab_indx(perms,rank_ord) = winning_mod_bic_perms(perms,rank_ord)~=4;
    end
    
end

for rank_ord = 1:5
    p_value(rank_ord) = sum(hab_indx(:,rank_ord))./nperms ;
end