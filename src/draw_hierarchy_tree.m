function h=draw_hierarchy_tree(reconstructed_directed_adj, reconstructed_is_selected)
reconstructed_parent_vector = (1:size(reconstructed_directed_adj,1))*(reconstructed_directed_adj>0);
x = treelayout(reconstructed_parent_vector);
d2leave = leavedist(reconstructed_directed_adj);
depth_level = zeros(size(reconstructed_parent_vector)); depth_level(1)=1; level=1;
tmp = depth_level;
while ismember(0,depth_level)
    tmp = tmp*reconstructed_directed_adj;
    depth_level(tmp~=0) = level + 1;
    level = level + 1;
end
for d = unique(depth_level)
    idx = find(depth_level==d);
    % pos = pos - median(pos);
    if d == 1
        [~, ind] = sort(d2leave(idx));
        layout = 1;
    else
        pidx = reconstructed_parent_vector(idx);
        % [~, ind] = sortrows([d2leave(pidx), pidx', d2leave(idx)]);
        [~, ind] = sortrows([x(pidx)', pidx', d2leave(idx)],[1, 2, -3]);
        layout = sort(x(pidx));
    end
    %x(idx(ind)) = (1:length(idx)) - length(idx)/2;%pos(ind);
    for i = 2:length(layout)
        if layout(i) == layout(i-1)
            layout(i:end) = layout(i:end) + 1;
        end
    end
    x(idx(ind)) = layout;%(1:length(idx))/(length(idx));
    %x(idx(ind)) = 1:length(idx);%
end
% plot(x,-depth_level,'.')
% hold on 
% [i, j] = find(reconstructed_directed_adj >= 1 | reconstructed_directed_adj' >=1);
% for k = 1:length(i)
%     line(x([i(k) j(k)]), -depth_level([i(k) j(k)]));
% end

node_positions = [x;(max(depth_level)-depth_level+1)/max(depth_level)];
% g = graph(reconstructed_directed_adj+reconstructed_directed_adj');
% % h = plot(g, 'MarkerSize',1,'Layout','force','XStart',node_positions(1,:),...
% %     'YStart',node_positions(2,:),'Iterations',1000);
% % % layout(h,'force');
% % h.YData = node_positions(2,:);
% highlight(h, find(reconstructed_is_selected), 'NodeColor','r')
% highlight(h, 1, 'NodeColor','r', 'MarkerSize',10);
node_positions = spring_layout(reconstructed_directed_adj + reconstructed_directed_adj', node_positions');
node_positions(2,:) = (max(depth_level)-depth_level+1)/max(depth_level);

node_positions(1,:) = ((node_positions(1,:) - min(node_positions(1,:)))/(max(node_positions(1,:))-min(node_positions(1,:))) - 0.5)*2*50;
node_positions(2,:) = ((node_positions(2,:) - min(node_positions(2,:)))/(max(node_positions(2,:))-min(node_positions(2,:))) - 0.5)*2*50;

node_size = ones(1,size(node_positions,2)); 
node_size(reconstructed_is_selected==1) = 2;
node_size(1) = 10;
node_color_data = reconstructed_is_selected * [1 1 1];
%node_color_data(node_color_data==0)=NaN;
[~, h]=draw_SPADE_tree_annotation(reconstructed_directed_adj+reconstructed_directed_adj', node_positions, node_size, node_color_data, [-1,1], 1, 0, reconstructed_is_selected, 'jet', [], []);
