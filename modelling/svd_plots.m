function svd_plots(S,U,V,fitresult,winmod,title_txt,Snoise_mean,Snoise_ste,Unoise_mean,Unoise_ste,Vnoise_mean,Vnoise_ste)
%% ERP modelling code relative to:
% Mancini F, Pepe A, Bernacchia, A, Di Stefano G, Mouraux A, Iannetti GD. (2017)
% Characterising the short-term habituation of event-related evoked potentials
% eNeuro

% written in Matlab R2016b by F Mancini, fm456@cam.ac.uk

% number of standard errors for statistical comparison
nst = 2.58; % 2-tails, p 0.01
nst1 = 2.33; % 1-tail, p 0.01

% to convert from samples to seconds use function sample2sec.m


%% plot singular values for signal+noise and noise
figure('Name',title_txt); 
hold on
shadedErrorBar([],Snoise_mean,Snoise_ste.*nst1,'lineprops','-r','patchSaturation',0.1); %noise displayed as red line ± CIs
plot(diag(S),'ok','markerfacecolor','k'); %singla+noise
xlabel('Rank order')
ylabel('Singular value')
title(title_txt)
hold off

%% plot left- and right-singular values for signal+noise and noise at the first 3 ranks
% superimpose decay models to right-singular vectors only if a decay model
% wins over the model of non-habituation

figure('Name',title_txt);
hold on
subplot(2,3,1),shadedErrorBar([],Unoise_mean(:,1),Unoise_ste(:,1).*nst,'lineprops','-r','patchSaturation',0.1);
hold on
subplot(2,3,1),plot(U(:,1),'b','linewidth',1.5);
xlabel('Sample')
ylabel('Left-singular vector')
axis([0 1025 -0.12 0.1])
% axis ij
subplot(2,3,2),shadedErrorBar([],Unoise_mean(:,2),Unoise_ste(:,2).*nst,'lineprops','-r','patchSaturation',0.1);
hold on
subplot(2,3,2),plot(U(:,2),'m','linewidth',1.5);
xlabel('Sample')
axis([0 1025 -0.1 0.1])
% axis ij
subplot(2,3,3),shadedErrorBar([],Unoise_mean(:,3),Unoise_ste(:,3).*nst,'lineprops','-r','patchSaturation',0.1);
hold on
subplot(2,3,3),plot(U(:,3),'g','linewidth',1.5);
xlabel('Sample')
axis([0 1025 -0.1 0.1])
% axis ij


subplot(2,3,4),shadedErrorBar([],Vnoise_mean(:,1),Vnoise_ste(:,1).*nst,'lineprops','-r','patchSaturation',0.1);
hold on
subplot(2,3,4),plot( [1:60], V(:,1),'ob','markerfacecolor','b');
hold on
if winmod(1) ~= 4
    subplot(2,3,4),plot(fitresult{winmod(1),1},'k-', [1:60], V(:,1),'ob');
else
    subplot(2,3,4),plot( [1:60], V(:,1),'ob');
end
xlabel('Trial number')
ylabel('Right-singular vector')
axis([0 60 -0.5 0.5])
legend off

subplot(2,3,5),shadedErrorBar([],Vnoise_mean(:,2),Vnoise_ste(:,2).*nst,'lineprops','-r','patchSaturation',0.1);
hold on
subplot(2,3,5),plot(V(:,2),'om','markerfacecolor','m');
hold on
if winmod(2) ~= 4
    subplot(2,3,5),plot(fitresult{winmod(2),2},'m-', [1:60], V(:,2),'om' );
end
xlabel('Trial number')
ylabel(' ')
axis([0 60 -0.5 1])
legend off

subplot(2,3,6),shadedErrorBar([],Vnoise_mean(:,3),Vnoise_ste(:,3).*nst,'lineprops','-r','patchSaturation',0.1);
hold on
subplot(2,3,6),plot(V(:,3),'og','markerfacecolor','g');
hold on
if winmod(3) ~= 4
    subplot(2,3,6),plot(fitresult{winmod(3),3},'k-', [1:60], V(:,3),'og' );
end
xlabel('Trial number')
ylabel(' ')
axis([0 60 -0.5 1])
legend off


