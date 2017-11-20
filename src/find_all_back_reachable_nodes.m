%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reachable_nodes_ind = find_all_back_reachable_nodes(adj,end_node)
% for a directed adjacency matrix, given a set of end_node, find all nodes
% that are back reachable from the end_node's
reachable_nodes_ind = [];
reachable_queue = end_node(:);
while ~isempty(reachable_queue)
    reachable_queue = [reachable_queue; find(adj(:,reachable_queue(1))==1)];
    reachable_nodes_ind = [reachable_nodes_ind; reachable_queue(1)];
    reachable_queue(1) = [];
end
reachable_nodes_ind = unique(reachable_nodes_ind);
