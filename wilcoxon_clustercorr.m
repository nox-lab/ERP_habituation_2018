function wilcoxon_clustercorr(wave)

% Mancini F, Pepe A, Bernacchia A, Di Stefano G, Mouraux A, Iannetti GD. (2018)
% Characterising the short-term habituation of event-related evoked
% potentials. E-neuro.
%
% Written for Matlab R2016b

% The nonparametric cluster-based permutation approach for statistical testing 
% assumes that true neural activity will tend to generate signal changes over contiguous 
% time points. 
% 
% The input required is a 2D matrix of EEG data (rows=subjects or epochs, columns=samples).
% 
% First, the EEG waveform is compared by means of a point-by-point, one-sample 
% Wilcoxon signed-rank test against zero. 
% 
% Then, clusters of contiguous time points above the critical z-score for a 
% non-parametric two-sided test are identified, and an estimate of the magnitude 
% of each cluster is obtained by computing the sum of the z-scores constituting 
% each cluster (cluster-level statistic). Random permutation testing of the subject-specific 
% EEG waveforms of the different conditions (performed independently for every subject)
% is then used to obtain a reference distribution of mean cluster magnitude. 
% 
% Finally, the proportion of random partitions that results in a larger cluster-level statistic 
% than the observed one is calculated. Clusters in the observed data are regarded 
% as significant if they have a magnitude exceeding the threshold of the percentiles 
% defined by variable 'cluster_threshold' (specify 1 value for 1-sided test, 
% 2 values for a 2-sided test). 

% Output: 
% 'actual_tres_pvalue': uncorrected p-value for the point-by-point signed rank test
% 'actual_tres_Zvalue': uncorrected Z-score for the point-by-point signed rank test
% 'outdata_Zvalue': cluster-corrected Z-score for the point-by-point signed rank test
% 'outdata_pvalue': cluster-corrected p-value for the point-by-point signed rank test
% 'cluster_distribution': distribution of clusters

%% SET PARAMETERS & INITIALIZE

outfile='wilcoxon_output.mat';
figname='wilcoxon_output.eps';

alpha=0.05;
permutation=1; %1= perform permutation testing; 0=no permutation testing
num_permutations=1000;
cluster_statistic='perc_mean'; 
% cluster_threshold=95; this is for a one-sided test, p < 0.05
cluster_threshold=[2.5 97.5]; % two-sided test, p < 0.05


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

% load('wave_data.mat');

%% PREPARE DATA

% stim1.abeta.cz = squeeze(wave.abeta.sub(:,:,1,wave.abeta.chan.Cz));
% stim1.abeta.cc = squeeze(wave.abeta.sub_rrCc(:,:,1));
% stim1.adelta.cz = squeeze(wave.adelta.sub(:,:,1,wave.adelta.chan.Cz));
% stim1.adelta.cc = squeeze(wave.adelta.sub_rrCc(:,:,1));
% 
% tmp = squeeze(wave.abeta.sub(:,:,[6:60],wave.abeta.chan.Cz));
% stim660.abeta.cz = squeeze(nanmean(tmp,3)); %(sub,frame)
% 
% tmp=[];
% tmp = squeeze(wave.abeta.sub_rrCc(:,:,[6:60]));
% stim660.abeta.cc = squeeze(nanmean(tmp,3)); %(sub,frame)
% 
% tmp=[];
% tmp = squeeze(wave.adelta.sub(:,:,[6:60],wave.adelta.chan.Cz));
% stim660.adelta.cz = squeeze(nanmean(tmp,3)); %(sub,frame)
% 
% tmp=[];
% tmp = squeeze(wave.adelta.sub_rrCc(:,:,[6:60]));
% stim660.adelta.cc = squeeze(nanmean(tmp,3)); %(sub,frame)
% 
% 
% data=stim660.adelta.cz;  %% CHANGE THIS


%init cluster_distribution
cluster_distribution.mean_statistic=[];
cluster_distribution.max_statistic=[];


%% ONE_SAMPLE WILCOXON TEST ON THE OBSERVED DATA

for d = 1:size(data,2)
    [actual_tres_pvalue(d),H,STATS]=signrank(data(:,d));
    actual_tres_Zvalue(d)=STATS.zval;
end

% figure;plot(actual_tres_Zvalue);

%% CLUSTER-BASED PERMUTATION TESTING

%cluster thresholding?
if permutation==1
    
    %figure to draw evolution of criticals
    hf=figure;
    disp(['Performing cluster-based thresholding']);
    
    %loop
    blobsizes=[];
    for iter=1:num_permutations
        disp(['Permutation : ' num2str(iter)]);
        %loop through d
        for d=1:size(data,2)
            %permutation
            tmp = reshape(data(randperm(size(data,1)*size(data,2))),size(data));
            
            %wilcoxon (output is tres with p values and Z values tres_pvalue / tres_Zvalue)
            [p,H,STATS]=signrank(tmp(:,d));
            tres_pvalue(d)=p*H;
            tres_Zvalue(d)=STATS.zval*H;
            
        end
        
        %blobology
        RLL=bwlabel(tres_Zvalue);
        RLL_size=[];
        blobpos=1;
        for i=1:max(max(RLL))
            ff=find(RLL==i);
            v=sum(sum(abs(RLL(ff))));
            if v>0;
                RLL_size(blobpos)=v;
                blobpos=blobpos+1;
            end
        end
        if isempty(RLL_size);
            RLL_size=0;
        end
        
        %blob summary
        blob_size.size(iter)=mean(abs(RLL_size));
        blob_size_max.size(iter)=max(abs(RLL_size));
        
        %critical
        switch cluster_statistic
            case 'perc_mean'
                criticals=prctile(blob_size.size,cluster_threshold);
            case 'perc_max'
                criticals=prctile(blob_size_max.size,cluster_threshold);
        end
        
        tp_plot(iter,:)=squeeze(criticals(:,1));
        plot(tp_plot);
        drawnow;
    end
    
    %process actual data
    outdata_pvalue=actual_tres_pvalue;
    outdata_Zvalue=actual_tres_Zvalue;
    tres=zeros(size(outdata_pvalue));
    tp=find(outdata_pvalue<alpha);
    tres(tp)=1;
    blob_size=[];
    
    
    %loop through blobs
    tp2=bwlabel(tres);
    toutput_Zvalues=zeros(size(outdata_Zvalue));
    toutput_pvalues=ones(size(outdata_pvalue));
    for i=1:max(max(tp2))
        %sum Zvalues
        idx=find(tp2==i);
        blob_size=sum(sum(abs(outdata_Zvalue(idx))));
        disp(['B' num2str(i) ': ' num2str(blob_size)]);
        if sum(sum(tres(find(tp2==i))))>0
            if abs(blob_size)>criticals
                disp('FOUND a significant cluster!');
                toutput_Zvalues(idx)=outdata_Zvalue(idx);
                toutput_pvalues(idx)=outdata_pvalue(idx);
            end
        end
    end
    outdata_Zvalue=toutput_Zvalues;
    outdata_pvalue=toutput_pvalues;
    fig=figure;plot(outdata_Zvalue);
    axis([0 1024 -5 5])
    axis ij
    saveas(fig,figname);
    
    %cluster_distribution
    cluster_distribution.mean_statistic=blob_size;
    cluster_distribution.max_statistic=blob_size_max;
    save(outfile,'actual_tres_pvalue','actual_tres_Zvalue','outdata_Zvalue','outdata_pvalue','cluster_distribution');
    
end