function horizonatal_hierarchy(adj, nodeposition)
adj = adj + adj';
adj(adj>1) = 1;
