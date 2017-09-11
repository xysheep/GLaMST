function [v,V, edge_list, all_edit_operations, edit_dist, unique_operations, unique_operations_weights] = EditDistance_all_fastest(string1,string2)
% Given two strings s1 and s2 (DNA or protein), the edit distance between 
% them is the minimum number of operations required to convert s1 to s2. 
% Allowed operations are:
%   Replacing one character of string by another character.
%   Deleting a character from string
%   Adding a character to string
% Cost of each type of operation are all the same, 1

%unused variables
all_edit_operations = [];
edge_list = [];
edit_dist =[];
% compute edit distance matrix
m=length(string1);
n=length(string2);
V=zeros(m+1,n+1);

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
% [v,V] = editDist_only(string1,string2,m,n); 


% [i,j] = size(V);
% step = v;
% unique_operations = cell(v,1);
% unique_operations_weights = ones(v,1);
% while (i>1 || j >1)
%     if i==1
%         j = j - 1;
%         unique_operations{step} = sprintf('insert %s after position %d',string2(j),i - 1);
%         step = step - 1;
%     elseif j==1
%         i = i - 1;
%         unique_operations{step} = sprintf('delete position %d',i);
%         step = step - 1;
%     else
%         if V(i,j) == V(i-1,j-1)
%             i = i - 1;
%             j = j - 1;
%         elseif V(i,j) == V(i-1,j-1) + 1
%             i = i - 1;
%             j = j - 1;
%             unique_operations{step} = sprintf('mutate position %d to %s',i, string2(j));
%             step = step - 1;
%         elseif V(i,j) == V(i,j-1) + 1
%             j = j - 1;
%             unique_operations{step} = sprintf('insert %s after position %d',string2(j),i - 1);
%             step = step - 1;
%         elseif V(i,j) == V(i-1,j) + 1
%             i = i - 1;
%             unique_operations{step} = sprintf('delete position %d',i);
%             step = step - 1;
%         end
%     end
% %     fprintf('i=%d\tj=%d\n',i,j);
% end

% build a graph with only the paths that lead to the end node
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


% legacy variables from the pervious version called faster
all_edit_operations = [];
edit_dist = v; 

% find all unique edit operations and their weights
start_node = 1;             
end_node = length(nodes);   
total_num_path = find_num_path(adj, start_node, end_node);

[i,j] = find(adj==1);
edit_operations = [];
edit_operations_weights = [];
for k=1:length(i)
    previous_state = nodes(i(k),:);
    current_state = nodes(j(k),:);
    if current_state(1) == previous_state(1)+1 && current_state(2) == previous_state(2)+1  && string1(current_state(1))==string2(current_state(2))
        continue;
    end
    if current_state(1) == previous_state(1)+1 && current_state(2) == previous_state(2)+1  && string1(current_state(1))~=string2(current_state(2))
        edit_operations = [edit_operations; {['mutate positioin ', num2str(current_state(1)), ' to ', string2(current_state(2))]}];
    end
    if current_state(1) == previous_state(1)+1 && current_state(2) == previous_state(2)  
        edit_operations = [edit_operations; {['delete positioin ', num2str(current_state(1))]}];
    end
    if current_state(1) == previous_state(1) && current_state(2) == previous_state(2)+1  
        edit_operations = [edit_operations; {['insert ', string2(current_state(2)), ' after positioin ', num2str(previous_state(1))]}];
    end
    if length(edit_operations)~=length(unique(edit_operations))
        1;
    end
    
    tmp_adj = adj;
    tmp_adj(i(k), j(k)) = 0; 
    tmp_num_path = find_num_path(tmp_adj, start_node, end_node);
    edit_operations_weights = [edit_operations_weights, 1 - tmp_num_path/total_num_path];    
end

unique_operations = unique(edit_operations);
unique_operations_weights = [];
for k=1:length(unique_operations)
    unique_operations_weights(k,1) = sum(edit_operations_weights(ismember(edit_operations, unique_operations(k))));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function num_path = find_num_path(adj, start_node, end_node)
% adj needs to be directed, and no cycles

% num_path = zeros(size(adj,1),1);
% num_path(start_node)=1;
% for i=1:size(adj,1)
%     tmp = max(num_path,adj'*num_path);
%     if isequal(tmp,num_path)
%         break;
%     else
%         num_path = tmp;
%     end
% end
% num_path = num_path(end_node);

num_path = zeros(size(adj,1),1);
num_path(start_node)=1;
nodes_updated = 1;
while ~isempty(nodes_updated)
    nodes_to_update = find(sum(adj(nodes_updated,:),1)~=0);
    nodes_updated = [];
    for k=1:length(nodes_to_update)
        new_num = sum(num_path(adj(:,nodes_to_update(k))==1));
        if num_path(nodes_to_update(k)) ~= new_num
            num_path(nodes_to_update(k)) = new_num;
            nodes_updated = [nodes_updated, nodes_to_update(k)];
        end
    end
end
num_path = num_path(end_node);

% num_path = 0;
% current_k_step_adj = eye(size(adj)); 
% for k=1:size(adj,1)
%     current_k_step_adj = current_k_step_adj * adj;
%     num_path = num_path + current_k_step_adj(start_node, end_node);
% end

