function node_positions = spring_layout_v2(adj, ini_pos)
G = graph(adj + adj');

fh = figure(1);
set(fh, 'Visible','off')
h = plot(G);
layout(h, 'force', 'XStart',ini_pos(:,1),'YStart',ini_pos(:,1),'Iterations',500);
node_positions = [h.XData;h.YData];
close(fh);
node_positions = node_positions - repmat((max(node_positions,[],2)+min(node_positions,[],2))/2,1,size(node_positions,2));
node_positions = node_positions/max(abs(node_positions(:)))*50;
