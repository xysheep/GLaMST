function treeplotA(A)
[~, p] = max(A);
p(1) = 0;
treeplot(p);