function [all_sequences,mst_adj,reconstructed_observed_indicator, reconstructed_directed_adj] = reconstruct_tree_minimun_tree_size(observed_sequences)


fprintf('Reconstruction:\n');
tic;
fprintf('    Initialization compute pairwise edit distance: %5d / %5d',0,0);
pairwise_dist = zeros(length(observed_sequences),length(observed_sequences));
for i=1:length(observed_sequences)
    for j=i+1:length(observed_sequences)
        % [v,V,edge_list,edit_operations, edit_dist] = EditDistance_all_faster(observed_sequences{i},observed_sequences{j});
        % pairwise_dist(i,j) = min(edit_dist);
        [v,~] = EditDistance_only(observed_sequences{i},observed_sequences{j});
        pairwise_dist(i,j) = v;
        pairwise_dist(j,i) = pairwise_dist(i,j);
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%5d / %5d',i,j);
    end
end
toc
seed_adj = zeros(size(pairwise_dist));
[mst_adj,adj2, ~] = mst_from_dist_matrix(pairwise_dist, seed_adj);
seed_adj = sparse(double(adj2==1));

tic
all_sequences = observed_sequences;
reconstructed_indicator = zeros(size(all_sequences)); reconstructed_indicator(1) = 1;
reconstructed_observed_indicator = zeros(size(all_sequences)); reconstructed_observed_indicator(1) = 1;
reconstructed_directed_adj = zeros(size(all_sequences));

fprintf('    Iteratived build the tree: %5d / %5d',sum(reconstructed_indicator),length(observed_sequences)-sum(reconstructed_observed_indicator));

while sum(reconstructed_indicator==0)~=0
    tic
    
    tmp_adj = mst_adj;
    tmp_adj(reconstructed_indicator==1,:)=0;
    tmp_adj(:,reconstructed_indicator==1)=0;
    
    [~, ~, ~,components] = extract_connected_component(tmp_adj);
    components(:,sum(components(reconstructed_indicator==1,:),1)==1)=[];  % remove components corresponding to already reconstructed portion

    merged_parent_ind = [];
    merged_components = [];
    for i=1:size(components,2)
        tmp_parent_ind = find(sum(mst_adj(:,components(:,i)==1),2)~=0 & reconstructed_indicator==1);
        for j=1:length(tmp_parent_ind)
            if ~ismember(tmp_parent_ind(j), merged_parent_ind)
                merged_parent_ind = [merged_parent_ind, tmp_parent_ind(j)];
                merged_components = [merged_components, components(:,i)];
            else
                tmp_I = find(ismember(merged_parent_ind, tmp_parent_ind(j))==1);
                merged_components(:,tmp_I) = merged_components(:,tmp_I) + components(:,i);
            end
        end    
    end
    parent_ind = merged_parent_ind;
    components = merged_components;
    

    % for each components, find the best operations with the best score
    tmp_operations = [];
    tmp_counts = [];
    tmp_parents = [];
    for i=1:size(components,2)
        children_ind = find(components(:,i)==1);
        all_unique_operations = cell(0);
        all_unique_operations_weights = [];
        for j=1:length(children_ind)
            %fprintf('\n i = %d\tj = %d\n',parent_ind(i),children_ind(j));
            [~,~,~,~, ~, unique_operations, unique_operations_weights] = EditDistance_all_fastest(all_sequences{parent_ind(i)},all_sequences{children_ind(j)});
            all_unique_operations = [all_unique_operations; unique_operations];
            all_unique_operations_weights = [all_unique_operations_weights; unique_operations_weights];
        end
        [a,b] = grpstats(all_unique_operations_weights, all_unique_operations, {'numel','gname'});
        tmp_counts = [tmp_counts; a];
        tmp_operations = [tmp_operations;b];
        tmp_parents = [tmp_parents; parent_ind(i)*ones(size(a))];
    end
    
    % pick the highest count
    [~,ind] = sort(tmp_counts,'descend');
    
    while 1
        best_parent = tmp_parents(ind(1));
        best_operation = tmp_operations{ind(1)};
    
        % create a new node from the best_parent and best_operation
        new_sequence = all_sequences{best_parent};
        tmp = regexp(best_operation,' ','split');
        if isequal(tmp{1},'mutate')
            new_sequence(str2num(tmp{3})) = tmp{5};
        elseif isequal(tmp{1},'delete')
            new_sequence(str2num(tmp{3})) = [];
        elseif isequal(tmp{1},'insert')
            new_sequence = [new_sequence(1:str2num(tmp{5})), tmp{2}, new_sequence(str2num(tmp{5})+1:length(new_sequence))];
        end
    
        if ismember({new_sequence},all_sequences(length(observed_sequences)+1:end)) || ismember({new_sequence},all_sequences(reconstructed_indicator==1)) 
            ind(1)=[];
        else
            break;
        end
    end
    
    
    % add the new sequence to the all_sequences
    ind = find(ismember(all_sequences, new_sequence));
    if ~isempty(ind)
        reconstructed_indicator(ind) = 1;
        reconstructed_observed_indicator(ind) = 1; 
        seed_adj(ind, best_parent) = 1;
        seed_adj(best_parent, ind) = 1;
        reconstructed_directed_adj(best_parent, ind) = 1;
        reconstructed_directed_adj(ind, best_parent) = 0;
        [mst_adj,adj2, cost_value] = mst_from_dist_matrix(pairwise_dist, seed_adj);
        1;
    else
        all_sequences = [all_sequences; {new_sequence}];
        reconstructed_indicator(end+1) = 1;
        reconstructed_observed_indicator(end+1) = 0;
        pairwise_dist(length(all_sequences),length(all_sequences)) = 0;
        for i=1:length(all_sequences)-1
            % [v,V,edge_list,edit_operations, edit_dist] = EditDistance_all_faster(all_sequences{end},all_sequences{i});
            % pairwise_dist(i,length(all_sequences)) = min(edit_dist);
            [v,~] = EditDistance_only(all_sequences{end},all_sequences{i});
            pairwise_dist(i,length(all_sequences)) = v;
            pairwise_dist(length(all_sequences),i) = pairwise_dist(i,length(all_sequences));
        end
        seed_adj(length(all_sequences), best_parent) = 1;
        seed_adj(best_parent, length(all_sequences)) = 1;
        reconstructed_directed_adj(best_parent, length(all_sequences)) = 1;
        reconstructed_directed_adj(length(all_sequences), best_parent) = 0;
        [mst_adj,adj2, cost_value] = mst_from_dist_matrix(pairwise_dist, seed_adj);
    end
        
    % display progress
%     fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%5d / %5d',sum(reconstructed_indicator),length(observed_sequences)-sum(reconstructed_observed_indicator));

    fprintf('    Iteratived build the tree: %5d / %5d',sum(reconstructed_indicator),length(observed_sequences)-sum(reconstructed_observed_indicator));
    toc
    save reconstruct_tree_minimun_tree_size_tmp.mat

end



% remove unnecessary nodes and edges
fprintf('    Trim the tree for un-necessary nodes ... \n'); 
reachable_nodes_ind = find_all_back_reachable_nodes(reconstructed_directed_adj,find(reconstructed_observed_indicator==1));
to_remove = setdiff(1:size(reconstructed_directed_adj,1),reachable_nodes_ind);
all_sequences(to_remove) = [];
mst_adj(to_remove,:) = [];
mst_adj(:,to_remove) = [];
reconstructed_observed_indicator(to_remove) = [];
reconstructed_directed_adj(to_remove,:) = [];
reconstructed_directed_adj(:,to_remove) = [];
pairwise_dist(to_remove,:) = [];
pairwise_dist(:,to_remove) = [];



% remove unnecessary nodes and edges
fprintf('    Rewire the tree to further reduce size ... %d \n', length(all_sequences)); 
while 1
    tic
    nodes_to_rewire = union(find(reconstructed_observed_indicator==1), find(sum(reconstructed_directed_adj,2)>1));  % observed node, and nodes at branching points
    tmp_directed_adj = reconstructed_directed_adj;
    tmp_directed_adj(nodes_to_rewire,:) = 0;
    cost_specific_to_nodes_to_rewire=[];
    cost_specific_to_rewire_the_node = [];
    for i=2:length(nodes_to_rewire)
        back_reachable_nodes = setdiff(find_all_back_reachable_nodes(tmp_directed_adj,nodes_to_rewire(i)), nodes_to_rewire(i));
        forward_reachable_nodes = setdiff(find_all_back_reachable_nodes(reconstructed_directed_adj',nodes_to_rewire(i)), nodes_to_rewire(i));
        if length(back_reachable_nodes)~=0
            cost_specific_to_nodes_to_rewire(i) = length(back_reachable_nodes);
        else
            cost_specific_to_nodes_to_rewire(i) = (1-reconstructed_observed_indicator(find(reconstructed_directed_adj(:,nodes_to_rewire(i))==1)))*0.1;
        end
        unqualified_rewire_destination = unique([back_reachable_nodes(:); nodes_to_rewire(i); find(reconstructed_directed_adj(:,nodes_to_rewire(i))==1); forward_reachable_nodes(:)]);
        tmp = pairwise_dist(:,nodes_to_rewire(i));  
        tmp(unqualified_rewire_destination) = size(reconstructed_directed_adj,1)+10;
        if (min(tmp)-1)~=0
            cost_specific_to_rewire_the_node(i) = min(tmp)-1; 
        else
            cost_specific_to_rewire_the_node(i) = (1-(sum(reconstructed_observed_indicator(tmp==min(tmp)))~=0))*0.1;
        end
    end
    % [cost_specific_to_nodes_to_rewire;cost_specific_to_rewire_the_node]
    if sum(cost_specific_to_nodes_to_rewire>cost_specific_to_rewire_the_node)==0  % if no new wiring is good, break
        break;
    end
    
    % find the new wiring destination
    rewire_scores = cost_specific_to_nodes_to_rewire-cost_specific_to_rewire_the_node;
    [~,i] = max(rewire_scores);
    node_to_rewire = nodes_to_rewire(i);
    back_reachable_nodes = setdiff(find_all_back_reachable_nodes(tmp_directed_adj,nodes_to_rewire(i)), nodes_to_rewire(i));
    forward_reachable_nodes = setdiff(find_all_back_reachable_nodes(reconstructed_directed_adj',nodes_to_rewire(i)), nodes_to_rewire(i));
    unqualified_rewire_destination = unique([back_reachable_nodes(:); nodes_to_rewire(i); find(reconstructed_directed_adj(:,nodes_to_rewire(i))==1); forward_reachable_nodes(:)]);
    tmp = pairwise_dist(:,nodes_to_rewire(i));  
    tmp(unqualified_rewire_destination) = size(reconstructed_directed_adj,1)+10;
    rewire_destinations = find(tmp==min(tmp));
    if sum(reconstructed_observed_indicator(rewire_destinations))~=0
        rewire_destinations = rewire_destinations(reconstructed_observed_indicator(rewire_destinations)==1);
    end
    rewire_destination = rewire_destinations(1);
    
    % perform the rewiring
    parent_of_node_to_rewire = find(reconstructed_directed_adj(:,node_to_rewire)==1);
    reconstructed_directed_adj(parent_of_node_to_rewire, node_to_rewire)=0;
    mst_adj(parent_of_node_to_rewire, node_to_rewire)=0;
    mst_adj(node_to_rewire, parent_of_node_to_rewire)=0;

    reachable_nodes_ind = find_all_back_reachable_nodes(reconstructed_directed_adj,find(reconstructed_observed_indicator==1));
    to_remove = setdiff(1:size(reconstructed_directed_adj,1),reachable_nodes_ind);
    
    if pairwise_dist(rewire_destination, node_to_rewire)==1
        reconstructed_directed_adj(rewire_destination, node_to_rewire) = 1;
        reconstructed_directed_adj(node_to_rewire, rewire_destination) = 0;
        mst_adj(rewire_destination, node_to_rewire)=1;
        mst_adj(node_to_rewire, rewire_destination)=1;
    else
        % [~,~,~,edit_operations, edit_dist, ~, ~] = EditDistance_all_faster(all_sequences{rewire_destination},all_sequences{node_to_rewire});
        % [~,I] = min(edit_dist); edit_operations = edit_operations{I};
        [~,~,~,~, edit_dist, ~, ~] = EditDistance_all_fastest(all_sequences{rewire_destination},all_sequences{node_to_rewire});
        best_parent = rewire_destination;
        for i=1:edit_dist 
            % [~,~,~,edit_operations, edit_dist, ~, ~] = EditDistance_all_faster(all_sequences{best_parent},all_sequences{node_to_rewire});
            % [~,I] = min(edit_dist); edit_operations = edit_operations{I};
            % best_operation = edit_operations{1};
            [~,~,~,~, ~, unique_operations, unique_operations_weights] = EditDistance_all_fastest(all_sequences{best_parent},all_sequences{node_to_rewire});
            unique_operations_position = zeros(size(unique_operations));
            for k=1:length(unique_operations)
                tmp = [unique_operations{k}(strfind(unique_operations{k},'positioin')+10 : end), ' '];
                tmp = tmp(1:find(tmp==' ',1)-1);
                unique_operations_position(k) = str2num(tmp);
            end
            tmp = find(unique_operations_position==min(unique_operations_position(unique_operations_weights == max(unique_operations_weights))) & unique_operations_weights == max(unique_operations_weights)); 
            tmp = tmp(1);
            best_operation = unique_operations{tmp};
            
            new_sequence = all_sequences{best_parent};
            % create a new node from the best_parent and best_operation
            tmp = regexp(best_operation,' ','split');
            if isequal(tmp{1},'mutate')
                new_sequence(str2num(tmp{3})) = tmp{5};
            elseif isequal(tmp{1},'delete')
                new_sequence(str2num(tmp{3})) = [];
            elseif isequal(tmp{1},'insert')
                new_sequence = [new_sequence(1:str2num(tmp{5})), tmp{2}, new_sequence(str2num(tmp{5})+1:length(new_sequence))];
            end
            
            if isequal(new_sequence, all_sequences{node_to_rewire})
                reconstructed_directed_adj(best_parent, node_to_rewire) = 1;
                reconstructed_directed_adj(node_to_rewire, best_parent) = 0;
                mst_adj(best_parent, node_to_rewire)=1;
                mst_adj(node_to_rewire, best_parent)=1;
                break;
            end
            
            if ismember({new_sequence},all_sequences(setdiff(1:end,to_remove)))
                new_sequence_exist = setdiff(find(ismember(all_sequences, {new_sequence})), to_remove);
                if ismember(new_sequence_exist, find_all_back_reachable_nodes(reconstructed_directed_adj',node_to_rewire));
                    reconstructed_directed_adj(best_parent, new_sequence_exist) = 1;
                    reconstructed_directed_adj(new_sequence_exist, best_parent) = 0;
                    mst_adj(best_parent, new_sequence_exist)=1;
                    mst_adj(new_sequence_exist, best_parent)=1;
                end
                best_parent = find(ismember(all_sequences, {new_sequence}));
            else
                all_sequences = [all_sequences; {new_sequence}];
                reconstructed_indicator(end+1) = 1;
                reconstructed_observed_indicator(end+1) = 0;
                pairwise_dist(length(all_sequences),length(all_sequences)) = 0;
                for i=1:length(all_sequences)-1
                    % [~,~,~,~, edit_dist] = EditDistance_all_faster(all_sequences{end},all_sequences{i});
                    % pairwise_dist(i,length(all_sequences)) = min(edit_dist);
                    [v,V] = EditDistance_only(all_sequences{end},all_sequences{i});
                    pairwise_dist(i,length(all_sequences)) = v;
                    pairwise_dist(length(all_sequences),i) = pairwise_dist(i,length(all_sequences));
                end
                mst_adj(length(all_sequences), best_parent) = 1;
                mst_adj(best_parent, length(all_sequences)) = 1;
                reconstructed_directed_adj(best_parent, length(all_sequences)) = 1;
                reconstructed_directed_adj(length(all_sequences), best_parent) = 0;
                best_parent = length(all_sequences);
            end
            
        end
    end
    
    reachable_nodes_ind = find_all_back_reachable_nodes(reconstructed_directed_adj,find(reconstructed_observed_indicator==1));
    to_remove = setdiff(1:size(reconstructed_directed_adj,1),reachable_nodes_ind);
    all_sequences(to_remove) = [];
    mst_adj(to_remove,:) = [];
    mst_adj(:,to_remove) = [];
    reconstructed_observed_indicator(to_remove) = [];
    reconstructed_directed_adj(to_remove,:) = [];
    reconstructed_directed_adj(:,to_remove) = [];
    pairwise_dist(to_remove,:) = [];
    pairwise_dist(:,to_remove) = [];
    
    fprintf('    Rewire the tree to further reduce size ... %d   ', length(all_sequences)); 
    toc
end




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

