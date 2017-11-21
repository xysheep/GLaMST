addpath src
addpath lib
prefix = 'data\real';
[b, c, bselected, cselected] = prepareadj2(prefix);
figtr = cell2mat(struct2cell(treestats(b)))';
fpeng = cell2mat(struct2cell(treestats(c)))';
comp = full([figtr;fpeng]);
%% Visualize the trees
figure;
draw_hierarchy_tree(c, cselected);
title('Reconstructed tree from Peng')
figure;
draw_hierarchy_tree(b, bselected);
title('Reconstructed tree from Igtree')