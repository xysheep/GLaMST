addpath src
addpath lib
clear
prefix = 'data\real';
[b, c, d, bselected, cselected, dselected] = prepareadj2(prefix);
figtr = cell2mat(struct2cell(treestats(b)))';
fpeng = cell2mat(struct2cell(treestats(c)))';
comp = full([figtr;fpeng]);
%% Visualize the trees
close all;
figure('pos',[300 300 800 600]);
h=draw_hierarchy_tree(c, cselected);
title('Reconstructed tree from GLaMST')
figure('pos',[1100 300 800 600]);
draw_hierarchy_tree(b, bselected);
title('Reconstructed tree from Igtree')
figure('pos',[1900 300 800 600]);
draw_hierarchy_tree(d, dselected);
title('Reconstructed tree from dnapars')