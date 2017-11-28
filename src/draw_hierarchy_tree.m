function draw_hierarchy_tree(reconstructed_directed_adj, reconstructed_is_selected)
reconstructed_parent_vector = (1:size(reconstructed_directed_adj,1))*reconstructed_directed_adj;
x = treelayout(reconstructed_parent_vector);
depth_level = zeros(size(reconstructed_parent_vector)); depth_level(1)=1; level=1;
tmp = depth_level;
while ismember(0,depth_level)
    tmp = tmp*reconstructed_directed_adj;
    depth_level(tmp~=0) = level + 1;
    level = level + 1;
end
% plot(x,(max(depth_level)-depth_level+1)/max(depth_level),'o')

node_positions = [x;(max(depth_level)-depth_level+1)/max(depth_level)];
node_positions(1,:) = ((node_positions(1,:) - min(node_positions(1,:)))/(max(node_positions(1,:))-min(node_positions(1,:))) - 0.5)*2*50;
node_positions(2,:) = ((node_positions(2,:) - min(node_positions(2,:)))/(max(node_positions(2,:))-min(node_positions(2,:))) - 0.5)*2*50;
node_size = 7*ones(1,size(node_positions,2)); node_size(1) = 15;
node_color_data = reconstructed_is_selected;
node_color_data(node_color_data==0)=NaN;
draw_SPADE_tree_annotation(reconstructed_directed_adj+reconstructed_directed_adj', node_positions, node_size, node_color_data, [-1,1], 1, 0, reconstructed_is_selected, 'jet', [], []);
