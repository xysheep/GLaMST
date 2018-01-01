function node_positions = spring_layout(graph_adj, playout)

A = sparse(graph_adj);
[X,~,~]=kamada_kawai_spring_layout_mex(...
    A, 1e-20, 20000, 1, ...  % adj, tolerance, max iteration, spring_constant
    playout, 1, 0, 'matrix');     % progressive_opt, options.edge_length, edge_weights, edge_weight_opt;
node_positions = X';         % use spring embedding to determine node positoin


node_positions = node_positions - repmat((max(node_positions,[],2)+min(node_positions,[],2))/2,1,size(node_positions,2));
node_positions = node_positions/max(abs(node_positions(:)))*50;
