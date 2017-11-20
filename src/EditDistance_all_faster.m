function [v,V, edge_list, all_edit_operations, edit_dist, unique_operations, unique_operations_weights] = EditDistance_all_faster(string1,string2)
% Given two strings s1 and s2 (DNA or protein), the edit distance between 
% them is the minimum number of operations required to convert s1 to s2. 
% Allowed operations are:
%   Replacing one character of string by another character.
%   Deleting a character from string
%   Adding a character to string
% Cost of each type of operation are all the same, 1


% compute edit distance matrix


m=length(string1);
n=length(string2);
V=zeros(m+1,n+1);
edge_list = [];
for i=1:1:m
    V(i+1,1)=i;
end
for j=1:1:n
    V(1,j+1)=j;
end
for i=1:m
    for j=1:n
        if (string1(i) == string2(j))
            V(i+1,j+1)=V(i,j);
        else
            V(i+1,j+1)=1+min(min(V(i+1,j),V(i,j+1)),V(i,j));
        end
    end
end
v=V(m+1,n+1);


% % construct node names for subsequent convenience
% tmp = num2str([reshape(repmat((0:m)',1,n+1),(m+1)*(n+1),1),reshape(repmat(0:n,m+1,1),(m+1)*(n+1),1)],'%d,');
% tmp = [tmp,repmat('+',size(tmp,1),1)]';
% tmp = tmp(:)';
% tmp(tmp==' ')=[];
% tmp(end)=[];
% nodes = regexp(tmp,'+','split')';
% 
% 
% % delete nodes that are not back-reachable
% start_node = 1;             % find(ismember(nodes, {[num2str(0),',',num2str(0),',']}));
% end_node = length(nodes);   % find(ismember(nodes, {[num2str(m),',',num2str(n),',']}));
% reachable_nodes_ind = sort(find_all_back_reachable_nodes(adj,end_node));
% adj = adj(reachable_nodes_ind,reachable_nodes_ind);
% nodes = nodes(reachable_nodes_ind);

tmp = [m,n]; 
nodes = zeros((m+n)*5,2);
adj = zeros((m+n)*5,(m+n)*5);
nodes(1,:) = tmp;
num_nodes = 1;
current_queue = 1;
while ~isempty(current_queue)
    current_node_position = nodes(current_queue(1),:);
    i = current_node_position(1);
    j = current_node_position(2);
    if i~=0 && j~=0 && V(i+1,j+1) == V(i,j) && (string1(i) == string2(j))
        ind = find(nodes(1:num_nodes,1)==i-1 & nodes(1:num_nodes,2)==j-1);
        if isempty(ind)
            nodes(num_nodes+1,:) = [i-1,j-1];
            adj(num_nodes+1, current_queue(1)) = 1;
            adj(current_queue(1), num_nodes+1) = 0;
            num_nodes = num_nodes + 1;
            current_queue = [current_queue, num_nodes];
        else
            adj(ind, current_queue(1)) = 1;
            adj(current_queue(1), ind) = 0;
        end
    else
        if j~=0 && V(i+1,j+1) == V(i+1,j)+1
            ind = find(nodes(1:num_nodes,1)==i & nodes(1:num_nodes,2)==j-1);
            if isempty(ind)
                nodes(num_nodes+1,:) = [i,j-1];
                adj(num_nodes+1, current_queue(1)) = 1;
                adj(current_queue(1), num_nodes+1) = 0;
                num_nodes = num_nodes + 1;
                current_queue = [current_queue, num_nodes];
            else
                adj(ind, current_queue(1)) = 1;
                adj(current_queue(1), ind) = 0;
            end
        end
        if i~=0 && V(i+1,j+1) == V(i,j+1)+1
            ind = find(nodes(1:num_nodes,1)==i-1 & nodes(1:num_nodes,2)==j);
            if isempty(ind)
                nodes(num_nodes+1,:) = [i-1,j];
                adj(num_nodes+1, current_queue(1)) = 1;
                adj(current_queue(1), num_nodes+1) = 0;
                num_nodes = num_nodes + 1;
                current_queue = [current_queue, num_nodes];
            else
                adj(ind, current_queue(1)) = 1;
                adj(current_queue(1), ind) = 0;
            end
        end
        if i~=0 && j~=0 && V(i+1,j+1) == V(i,j)+1
            ind = find(nodes(1:num_nodes,1)==i-1 & nodes(1:num_nodes,2)==j-1);
            if isempty(ind)
                nodes(num_nodes+1,:) = [i-1,j-1];
                adj(num_nodes+1, current_queue(1)) = 1;
                adj(current_queue(1), num_nodes+1) = 0;
                num_nodes = num_nodes + 1;
                current_queue = [current_queue, num_nodes];
            else
                adj(ind, current_queue(1)) = 1;
                adj(current_queue(1), ind) = 0;
            end
        end
    end
    current_queue(1) = [];
end

if num_nodes<size(nodes,1)
    nodes(num_nodes+1:end,:)=[];
    adj(num_nodes+1:end,:)=[];
    adj(:,num_nodes+1:end)=[];
end
[~,I] = sort(nodes(:,1) + nodes(:,2)*max(max(nodes)));
adj = adj(I,I);
nodes = nodes(I,:);


% find all paths from start to end node
start_node = 1;             % find(ismember(nodes, {[num2str(0),',',num2str(0),',']}));
end_node = length(nodes);   % find(ismember(nodes, {[num2str(m),',',num2str(n),',']}));
[all_paths] = find_all_paths(adj, start_node, end_node);


all_edit_operations = [];
all_edit_operations_length = [];
for k=1:length(all_paths)
    edit_path = nodes(all_paths{k},:);
    edit_operations = [];
    previous_state = edit_path(1,:);
    for i=2:length(edit_path)
        current_state = edit_path(i,:);
        if current_state(1) == previous_state(1)+1 && current_state(2) == previous_state(2)+1  && string1(current_state(1))~=string2(current_state(2))
            edit_operations = [edit_operations; {['mutate positioin ', num2str(current_state(1)), ' to ', string2(current_state(2))]}];
        end
        if current_state(1) == previous_state(1)+1 && current_state(2) == previous_state(2)  
            edit_operations = [edit_operations; {['delete positioin ', num2str(current_state(1))]}];
        end
        if current_state(1) == previous_state(1) && current_state(2) == previous_state(2)+1  
            edit_operations = [edit_operations; {['insert ', string2(current_state(2)), ' after positioin ', num2str(previous_state(1))]}];
        end
        previous_state = current_state;
    end
    
    all_edit_operations{k} = edit_operations;
    all_edit_operations_length(k) = length(edit_operations);
end


% remove duplicates
paths_to_remove = [];
for i=1:length(all_edit_operations)
    for j=i+1:length(all_edit_operations)
        if isequal(all_edit_operations(i),all_edit_operations(j))
            paths_to_remove = [paths_to_remove,j];
        end
    end
end
all_edit_operations(paths_to_remove)=[];
all_edit_operations_length(paths_to_remove)=[];


% compute summary: min dist, unique operations and their weights
all_individual_operations = [];
for k=1:length(all_edit_operations)
    edit_dist(k) = length(all_edit_operations{k});
    all_individual_operations = [all_individual_operations;all_edit_operations{k}];
end
unique_operations = unique(all_individual_operations);
unique_operations_weights = [];
for k=1:length(unique_operations)
    unique_operations_weights(k,1) = sum(ismember(all_individual_operations, unique_operations(k)))/length(all_edit_operations);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reachable_nodes_ind = find_all_back_reachable_nodes(adj,end_node)
reachable_nodes_ind = [];
reachable_queue = end_node;
while ~isempty(reachable_queue)
    reachable_queue = [reachable_queue; find(adj(:,reachable_queue(1))==1)];
    reachable_nodes_ind = [reachable_nodes_ind; reachable_queue(1)];
    reachable_queue(1) = [];
end
reachable_nodes_ind = unique(reachable_nodes_ind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all_paths] = find_all_paths(adj, start_node, end_node)

all_paths = [];
path_queue{1} = start_node;
while ~isempty(path_queue)
    if path_queue{1}(end)==end_node
        all_paths = [all_paths;path_queue(1)];
    else
        possible_next_node = find(adj(path_queue{1}(end),:)==1);
        if ~isempty(possible_next_node)
            for i=1:length(possible_next_node)
                path_queue(end+1) = {[path_queue{1},possible_next_node(i)]};
            end
        end        
    end
    path_queue(1)=[];
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist, previous_vertex, shortest_paths] = Dijkstra(graph, source)
% [dist, previous_vertex] = Dijkstra(graph, source)
%     This function finds the shortest graph distance from a source node to every other node in the graph      
%     
%     graph is an adjacency matrix, doesn't matter whether weighted or not
%     source is the starting node
%     dist is the approximated geodesics, length of the shortest path on graph 
%     previous_vertex is: to reach this point via shortest path, what is the previous point before arriving here
% 
% http://en.wikipedia.org/wiki/Dijkstra's_algorithm

graph = graph - diag(diag(graph));

% initialize algorithm
N = size(graph,1);
dist = zeros(1,N) + Inf;
previous_vertex = zeros(1,N);

% handle the source point
dist(source)=0;
previous_vertex(source) = source;

Q = 1:N; % All nodes in the graph are unoptimized - thus are in Q

while ~isempty(Q)  % The main loop
    [min_dist, u] = min(dist(Q)); u = Q(u); % u is the vertex in Q with smallest distance in dist[] ;
    if isinf(min_dist) 
        break;  % all remaining vertices are inaccessible from source
    end
    Q = setdiff(Q,u);
    for v = intersect(Q(:)', find(graph(u,:)~=0))
        alt = dist(u) + graph(u,v);
        if alt<dist(v)
            dist(v) = alt;
            previous_vertex(v) = u;
        end
    end
end

for i=1:N
    shortest_paths{i,1} = i;
end

for target=1:N
    while shortest_paths{target}(1)~=source && previous_vertex(shortest_paths{target}(1))~=0
        shortest_paths{target} = [shortest_paths{previous_vertex(shortest_paths{target}(1))},shortest_paths{target}];
        if previous_vertex(shortest_paths{target}(1))==0
            break;
        end
    end
end
