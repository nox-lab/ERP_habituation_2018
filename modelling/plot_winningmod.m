function plot_winningmod(winning_mod_bic)

%% ERP modelling code relative to:
% Mancini F, Pepe A, Bernacchia, A, Di Stefano G, Mouraux A, Iannetti GD. (2017)
% Characterising the short-term habituation of event-related evoked potentials
% eNeuro

% written in Matlab R2016b by F Mancini, fm456@cam.ac.uk

% % winning model (1-4) is color-coded in the figure:
% % pink = equation 1 wins
% % green = equation 2 wins
% % blue = = equation 3 wins
% % white = = equation 4 wins (no decay)

load('cmap.mat');
figure('Name','Winning models according to BIC'); hold on

subplot(4,1,1),imagesc(winning_mod_bic.b.cz);
colormap(cmap)
title('A-beta, Cz')
axis([0.5 60 0.5 1])

subplot(4,1,2),imagesc(winning_mod_bic.d.cz);
colormap(cmap)
title('A-delta, Cz')
axis([0.5 60 0.5 1])

subplot(4,1,3),imagesc(winning_mod_bic.b.cc);
colormap(cmap)
title('A-beta, Cc')
axis([0.5 60 0.5 1])

subplot(4,1,4),imagesc(winning_mod_bic.d.cc);
colormap(cmap)
title('A-delta, Cc')
axis([0.5 60 0.5 1])
xlabel('Rank number')

hold off