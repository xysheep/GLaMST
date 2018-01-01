addpath src
addpath lib
prefix = 'data\real';
[b, c, d, bselected, cselected, dselected] = prepareadj2(prefix);
figtr = cell2mat(struct2cell(treestats(b)))';
fpeng = cell2mat(struct2cell(treestats(c)))';
fdnapars = cell2mat(struct2cell(treestats(d)))';
comp = full([figtr;fpeng;fdnapars]);
%% Visualize the trees
figure('pos',[300 300 800 600]);
draw_hierarchy_tree(c, cselected);
title('Reconstructed tree from GLaMST')
figure('pos',[300 300 800 600]);
draw_hierarchy_tree(b, bselected);
title('Reconstructed tree from Igtree')
figure('pos',[300 300 800 600]);
draw_hierarchy_tree(d, dselected);
title('Reconstructed tree from dnapars')

save('trees.mat','b','c','d','bselected','cselected','dselected');