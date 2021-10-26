%% ERP habituation modelling code and data relative to:
% 
% Mancini F, Pepe A, Bernacchia A, Di Stefano G, Mouraux A, Iannetti GD. (2018)
% Characterising the short-term habituation of event-related evoked
% potentials. E-neuro.
% 
% Written in Matlab R2016b by F Mancini, fm456@cam.ac.uk
% It requires the Curve Fitting Toolbox in Matlab


%% INPUT
% wave.abeta.sub (subject,frame,trial,channel)
% wave.abeta.avgsub (frame,trial,channel)
% wave.abeta.sub_rrCc (subject,frame,trial)
% wave.abeta.sub_rrCc (frame,trial)
% wave.abeta.chan
% wave.adelta.sub (subject,frame,trial,channel)
% wave.adelta.avgsub (frame,trial,channel)
% wave.adelta.sub_rrCc (subject,frame,trial)
% wave.adelta.sub_rrCc (frame,trial)
% wave.adelta.chan

close all
clear all
load('wave_data.mat');


%% Singular value decomposition of the EEG matrix
% The EEG group average (sample, trial) is used here as the input for the
% SVD.
% The 'svd' function returns U (left-singular matrix, i.e. ERP component), S (diagonal matrix
% of singular values), V (right-singular matrix, i.e. habituation component).

% a-beta, vertex wave
[U.b.cz, S.b.cz, V.b.cz] = svd(squeeze(wave.abeta.avgsub(:,:,wave.abeta.chan.Cz)));

% a-beta, lateralized wave
[U.b.cc, S.b.cc, V.b.cc] = svd(wave.abeta.avgsub_rrCc);

% a-delta, vertex wave
[U.d.cz, S.d.cz, V.d.cz] = svd(squeeze(wave.adelta.avgsub(:,:,wave.adelta.chan.Cz)));

% a-delta, lateralized wave
[U.d.cc, S.d.cc, V.d.cc] = svd(wave.adelta.avgsub_rrCc);

%% Estimate SVD significance by performing SVD on single-subject noise traces..
% and calculating confidence intervals for statistical comparison

svd_noise_estimation;

%% Modelling of the right-singular vectors (habituation components)

[fitresult_b_cz, gof_b_cz, bic_b_cz] = postsvd_fitting(V.b.cz);
[fitresult_b_cc, gof_b_cc, bic_b_cc] = postsvd_fitting(V.b.cc);
[fitresult_d_cz, gof_d_cz, bic_d_cz] = postsvd_fitting(V.d.cz);
[fitresult_d_cc, gof_d_cc, bic_d_cc] = postsvd_fitting(V.d.cc);

%% Compare decay models according to BIC

for rank_ord = 1:60
    winning_mod_bic.b.cz(rank_ord) = find(bic_b_cz(:,rank_ord)==min((bic_b_cz(:,rank_ord))));
    winning_mod_bic.b.cc(rank_ord) = find(bic_b_cc(:,rank_ord)==min((bic_b_cc(:,rank_ord))));
    winning_mod_bic.d.cz(rank_ord) = find(bic_d_cz(:,rank_ord)==min((bic_d_cz(:,rank_ord))));
    winning_mod_bic.d.cc(rank_ord) = find(bic_d_cc(:,rank_ord)==min((bic_d_cc(:,rank_ord))));
end


%% Estimate p-values of winning models using resampling statistics

p_value_b_cz = resampling(V.b.cz);
p_value_b_cc = resampling(V.b.cc);
p_value_d_cz = resampling(V.d.cz);
p_value_d_cc = resampling(V.d.cc);

%% Make figure of winning models for each rank and condition

plot_winningmod(winning_mod_bic);

%% Make figures of SVD results

% a-beta
svd_plots(S.b.cz,U.b.cz,V.b.cz,fitresult_b_cz,winning_mod_bic.b.cz,'A-beta, Cz',S_noise_mean.b.cz, S_noise_ste.b.cz, U_noise_mean.b.cz, U_noise_ste.b.cz, V_noise_mean.b.cz, V_noise_ste.b.cz)
svd_plots(S.b.cc,U.b.cc,V.b.cc,fitresult_b_cc,winning_mod_bic.b.cc,'A-beta, Cc',S_noise_mean.b.cc, S_noise_ste.b.cc, U_noise_mean.b.cc, U_noise_ste.b.cc, V_noise_mean.b.cc, V_noise_ste.b.cc)

% a-delta
svd_plots(S.d.cz,U.d.cz,V.d.cz,fitresult_d_cz,winning_mod_bic.d.cz,'A-delta, Cz',S_noise_mean.d.cz, S_noise_ste.d.cz, U_noise_mean.d.cz, U_noise_ste.d.cz, V_noise_mean.d.cz, V_noise_ste.d.cz)
svd_plots(S.d.cc,U.d.cc,V.d.cc,fitresult_d_cc,winning_mod_bic.d.cc,'A-delta, Cc',S_noise_mean.d.cc, S_noise_ste.d.cc, U_noise_mean.d.cc, U_noise_ste.d.cc, V_noise_mean.d.cc, V_noise_ste.d.cc)

%% Save results

save('results.mat')
