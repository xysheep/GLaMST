function view_tree_here(mst_adj, node_color, is_selected)

N = size(mst_adj,1);
if ~exist('node_color')
    node_color = zeros(1,N);
end
if ~exist('is_selected')
    is_selected = zeros(1,N);
end
node_positions = spring_layout(mst_adj);
node_size = 7*ones(1,N); node_size(1) = 15;
draw_SPADE_tree_annotation(mst_adj, node_positions, node_size, node_color, [-1,1], 1, 0, is_selected, 'jet', [], []);
for k=1:size(node_positions,2) 
    text(node_positions(1,k)+2, node_positions(2,k), num2str(k), 'FontSize', 7); 
end

