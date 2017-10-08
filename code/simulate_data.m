function simulate_data(filename,sequence_length,num_tree_nodes,sample_size,operation_probability)

if ~exist('sequence_length')
    sequence_length = 80;
end
if ~exist('num_tree_nodes') || isempty(num_tree_nodes)
    num_tree_nodes = randsample(20:80,1);
end
if ~exist('sample_size') || isempty(sample_size)
    sample_size = randsample(2:min(24, num_tree_nodes-1),1);
end
if ~exist('operation_probability')
%     operation_probability = [0.98, 0.01, 0.01]; % 
%     operation_probability = [0.4,0.3,0.3]; % mutation, insertion, deletion
    operation_probability = [1, 0, 0]; % 
end



fprintf('Simulation:\n');
fprintf('    simulate a tree with %d nodes, and %d observed sequences (include root) .... %5d', num_tree_nodes, sample_size+1, 0);

% generate root sequence" 
base_pairs = 'ATCG';
root_sequence = base_pairs(randsample(1:4,sequence_length,true));

% generate mutated cells and ground truth tree structure
nodes = {root_sequence};
directed_adj = zeros(num_tree_nodes);
i=1;
longb = 0;
while 1
%     if length(longb) > 1
%         w = longb.^2./sum(longb.^2);
%     else
%         w = 1;
%     end
    % randomly pick a node
    node_ind = randsample(1:length(nodes),1);
    % randomly pick an operation
    operation_ind = randsample(1:3, 1, true,  operation_probability);
    % randomly pick a position in this seq
    if operation_ind==1 || operation_ind==3
        position_ind = randsample(1:length(nodes{node_ind}),1); % mutation and deletion much happen at an existing position
    else
        position_ind = randsample(0:length(nodes{node_ind}),1); % insertion can happen before the first position (after the 0th position), therefore, this operartion allows position_ind to be 0
    end
    if operation_ind==1
        new_sequence = nodes{node_ind};
        tmp = setdiff(base_pairs,new_sequence(position_ind));
        new_sequence(position_ind) = tmp(randsample(1:3,1));
    elseif operation_ind==2 % insertion  
        new_sequence = nodes{node_ind};
        new_sequence = [new_sequence(1:position_ind),base_pairs(randsample(1:4,1)),new_sequence(position_ind+1:end)];
    else % operation_ind must be 3 deletion
        new_sequence = nodes{node_ind};
        new_sequence(position_ind)=[];
    end
    % add the node to the ground truth tree
    if ~ismember({new_sequence},nodes)
        nodes = [nodes;{new_sequence}];
        directed_adj(node_ind,length(nodes))=1;
        directed_adj(length(nodes),node_ind)=0;
        longb(length(nodes)) = longb(node_ind) + 1;
        i=i+1;
    end
    fprintf('\b\b\b\b\b%5d', i);
    if i==num_tree_nodes
        fprintf('\n', i);
        break;
    end
end
adj = double((directed_adj + directed_adj')~=0); 
% stats = treestats(directed_adj);
% randomly select the observable sequences
observable_ind = randsample(2:num_tree_nodes, sample_size);
is_selected = zeros(1,num_tree_nodes);
is_selected(1) = 1;
is_selected(observable_ind)=1;


% % visualize ground truth
% N = size(adj,1);
% node_positions = spring_layout(adj);
% node_size = 7*ones(1,N); node_size(1) = 15;
% node_color = is_selected; 
% subplot(1,2,1);
% draw_SPADE_tree_annotation(adj, node_positions, node_size, node_color, [-1,1], 1, 0, is_selected, 'jet', [], []);
% % for k=1:size(node_positions,2), text(node_positions(1,k)+2, node_positions(2,k), num2str(k), 'FontSize', 7); end


% remove unnecessary nodes and edges
reachable_nodes_ind = find_all_back_reachable_nodes(directed_adj,find(is_selected==1));
to_remove = setdiff(1:size(directed_adj,1),reachable_nodes_ind);

% to_remove = [];
% for i=1:size(adj,1)
%     if is_selected(i)==1
%         continue;
%     end
%     e = zeros(size(adj,1),1); e(i)=1;
%     for k=1:size(adj,1)
%         e = ((eye(size(directed_adj))+directed_adj)'*e)~=0;
%     end
%     if sum(e==1 & is_selected(:)==1)==0
%         to_remove = [to_remove, i];
%     end
% end

adj(to_remove,:)=[];
adj(:,to_remove)=[];
directed_adj(to_remove,:)=[];
directed_adj(:,to_remove)=[];
nodes(to_remove)=[];
% node_positions(:,to_remove)=[];
% node_size(:,to_remove)=[];
% node_color(:,to_remove)=[];
is_selected(to_remove)=[];



% % visualize ground truth
% N = size(adj,1);
% subplot(1,2,2);
% draw_SPADE_tree_annotation(adj, node_positions, node_size, node_color, [-1,1], 1, 0, is_selected, 'jet', [], []);
% for k=1:size(node_positions,2), text(node_positions(1,k)+2, node_positions(2,k), num2str(k), 'FontSize', 7); end


observed_sequences = nodes(is_selected==1);
save(filename)

fprintf('    obtained tree size: %d nodes, %d observed sequences.\n', length(nodes), length(observed_sequences));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reachable_nodes_ind = find_all_back_reachable_nodes(adj,end_node)
reachable_nodes_ind = [];
reachable_queue = end_node(:);
while ~isempty(reachable_queue)
    reachable_queue = [reachable_queue; find(adj(:,reachable_queue(1))==1)];
    reachable_nodes_ind = [reachable_nodes_ind; reachable_queue(1)];
    reachable_queue(1) = [];
end
reachable_nodes_ind = unique(reachable_nodes_ind);