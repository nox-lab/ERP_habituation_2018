%% ERP habituation modelling code and data relative to:
% 
% Mancini F, Pepe A, Bernacchia A, Di Stefano G, Mouraux A, Iannetti GD. (2018)
% Characterising the short-term habituation of event-related evoked
% potentials. E-neuro.
% 
% Written in Matlab R2016b by F Mancini, fm456@cam.ac.uk
% It requires the Curve Fitting Toolbox in Matlab



%% calculate single-subject noise traces and perform SVD on these traces, in each condition

for sub = 1:size(wave.abeta.sub,1)
    
    data_sub = squeeze(wave.abeta.sub(sub,:,:,wave.abeta.chan.Cz));
    data_subset = squeeze(wave.abeta.sub(:,:,:,wave.abeta.chan.Cz));
    data_subset(sub,:,:) = [];
    data_subset_mean=squeeze(mean(squeeze(data_subset),1)); %calculate group mean after excluding subject
    data_sub_res = data_sub - data_subset_mean; %substract the group mean from single-subject data
    [U_res.b.cz(:,:,sub), S_res, V_res.b.cz(:,:,sub)] = svd(data_sub_res); %perform SVD
    S_noise.b.cz(sub,:)=diag(S_res)'; % save singular values
    
end

for sub = 1:size(wave.adelta.sub,1)
    
    data_sub = squeeze(wave.adelta.sub(sub,:,:,wave.adelta.chan.Cz));
    data_subset = squeeze(wave.adelta.sub(:,:,:,wave.adelta.chan.Cz));
    data_subset(sub,:,:) = [];
    data_subset_mean=squeeze(mean(squeeze(data_subset),1));
    data_sub_res = data_sub - data_subset_mean;
    [U_res.d.cz(:,:,sub), S_res, V_res.d.cz(:,:,sub)] = svd(data_sub_res);
    S_noise.d.cz(sub,:)=diag(S_res)';
    
end

for sub = 1:size(wave.abeta.sub_rrCc,1)
    
    data_sub = squeeze(wave.abeta.sub_rrCc(sub,:,:));
    data_subset = squeeze(wave.abeta.sub_rrCc(:,:,:));
    data_subset(sub,:,:) = [];
    data_subset_mean=squeeze(mean(squeeze(data_subset),1));
    data_sub_res = data_sub - data_subset_mean;
    [U_res.b.cc(:,:,sub), S_res, V_res.b.cc(:,:,sub)] = svd(data_sub_res);
    S_noise.b.cc(sub,:)=diag(S_res)';
    
end

for sub = 1:size(wave.adelta.sub_rrCc,1)
    
    data_sub = squeeze(wave.adelta.sub_rrCc(sub,:,:));
    data_subset = squeeze(wave.adelta.sub_rrCc(:,:,:));
    data_subset(sub,:,:) = [];
    data_subset_mean=squeeze(mean(squeeze(data_subset),1));
    data_sub_res = data_sub - data_subset_mean;
    [U_res.d.cc(:,:,sub), S_res, V_res.d.cc(:,:,sub)] = svd(data_sub_res);
    S_noise.d.cc(sub,:)=diag(S_res)';
    
end


%% Calculate mean and standard error of the singular, left-singular and right-singular vectors of noise traces  

N = size(wave.abeta.sub,1);
S_noise_mean.b.cz = mean(S_noise.b.cz,1)./sqrt(N);
S_noise_mean.b.cc = mean(S_noise.b.cc,1)./sqrt(N);
S_noise_mean.d.cz = mean(S_noise.d.cz,1)./sqrt(N);
S_noise_mean.d.cc = mean(S_noise.d.cc,1)./sqrt(N);

S_noise_ste.b.cz = ste(S_noise.b.cz,1);
S_noise_ste.b.cc = ste(S_noise.b.cc,1);
S_noise_ste.d.cz = ste(S_noise.d.cz,1);
S_noise_ste.d.cc = ste(S_noise.d.cc,1);

U_noise_mean.b.cz = mean(U_res.b.cz,3)./sqrt(N);
U_noise_mean.b.cc = mean(U_res.b.cc,3)./sqrt(N);
U_noise_mean.d.cz = mean(U_res.d.cz,3)./sqrt(N);
U_noise_mean.d.cc = mean(U_res.d.cc,3)./sqrt(N);

V_noise_mean.b.cz = mean(V_res.b.cz,3)./sqrt(N);
V_noise_mean.b.cc = mean(V_res.b.cc,3)./sqrt(N);
V_noise_mean.d.cz = mean(V_res.d.cz,3)./sqrt(N);
V_noise_mean.d.cc = mean(V_res.d.cc,3)./sqrt(N);

for rank = 1:60
    tmp= squeeze(U_res.b.cz(:,rank,:));
    U_noise_ste.b.cz(:,rank) = ste(tmp,2);
    
    tmp= squeeze(U_res.b.cc(:,rank,:));
    U_noise_ste.b.cc(:,rank) = ste(tmp,2);
    
    tmp= squeeze(U_res.d.cz(:,rank,:));
    U_noise_ste.d.cz(:,rank) = ste(tmp,2);
  
    tmp= squeeze(U_res.d.cc(:,rank,:));
    U_noise_ste.d.cc(:,rank) = ste(tmp,2);
  
    
    tmp= squeeze(V_res.b.cz(:,rank,:));
    V_noise_ste.b.cz(:,rank) = ste(tmp,2);
    
    tmp= squeeze(V_res.b.cc(:,rank,:));
    V_noise_ste.b.cc(:,rank) = ste(tmp,2);
    
    tmp= squeeze(V_res.d.cz(:,rank,:));
    V_noise_ste.d.cz(:,rank) = ste(tmp,2);
  
    tmp= squeeze(V_res.d.cc(:,rank,:));
    V_noise_ste.d.cc(:,rank) = ste(tmp,2); 
    
end

%% Test null Hypothesis that singular values of 'signal+noise' = singular values of 'noise'

nst = 2.33; % p = 0.01, one tail

for rank_ord = 1:length(S_noise_ste.b.cz)
    
    if S.b.cz(rank_ord,rank_ord) > (S_noise_mean.b.cz(rank_ord) + S_noise_ste.b.cz(rank_ord).*nst)
        S_sig.b.cz(rank_ord) = 1;
    else
        S_sig.b.cz(rank_ord) = 0;
    end
    
    if S.b.cc(rank_ord,rank_ord) > (S_noise_mean.b.cc(rank_ord) + S_noise_ste.b.cc(rank_ord).*nst)
        S_sig.b.cc(rank_ord) = 1;
      else
        S_sig.b.cc(rank_ord) = 0;
    end
    
    if S.d.cz(rank_ord,rank_ord) > (S_noise_mean.d.cz(rank_ord) + S_noise_ste.d.cz(rank_ord).*nst)
        S_sig.d.cz(rank_ord) = 1;
    else
        S_sig.d.cz(rank_ord) = 0;
    end
    
    if S.d.cc(rank_ord,rank_ord) > (S_noise_mean.d.cc(rank_ord) + S_noise_ste.d.cc(rank_ord).*nst)
        S_sig.d.cc(rank_ord) = 1;
     else
        S_sig.d.cc(rank_ord) = 0;
    end
       
end


clear('data_sub','data_sub_res','data_subset','data_subset_mean','tmp');

