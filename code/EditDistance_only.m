function [v,V] = EditDistance_only(string1,string2)
% Given two strings s1 and s2 (DNA or protein), the edit distance between 
% them is the minimum number of operations required to convert s1 to s2. 
% Allowed operations are:
%   Replacing one character of string by another character.
%   Deleting a character from string
%   Adding a character to string
% Cost of each type of operation are all the same, 1


% compute edit distance matrix

[v,V] = EditDistance_only_cpp(string1,string2);
% m=length(string1);
% n=length(string2);
% V=zeros(m+1,n+1);
% edge_list = [];
% for i=1:1:m
%     V(i+1,1)=i;
% end
% for j=1:1:n
%     V(1,j+1)=j;
% end
% for i=1:m
%     for j=1:n
%         if (string1(i) == string2(j))
%             V(i+1,j+1)=V(i,j);
%         else
%             V(i+1,j+1)=1+min(min(V(i+1,j),V(i,j+1)),V(i,j));
%         end
%     end
% end
% v=V(m+1,n+1);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function reachable_nodes_ind = find_all_back_reachable_nodes(adj,end_node)
% reachable_nodes_ind = [];
% reachable_queue = end_node;
% while ~isempty(reachable_queue)
%     reachable_queue = [reachable_queue; find(adj(:,reachable_queue(1))==1)];
%     reachable_nodes_ind = [reachable_nodes_ind; reachable_queue(1)];
%     reachable_queue(1) = [];
% end
% reachable_nodes_ind = unique(reachable_nodes_ind);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [all_paths] = find_all_paths(adj, start_node, end_node)
% 
% all_paths = [];
% path_queue{1} = start_node;
% while ~isempty(path_queue)
%     if path_queue{1}(end)==end_node
%         all_paths = [all_paths;path_queue(1)];
%     else
%         possible_next_node = find(adj(path_queue{1}(end),:)==1);
%         if ~isempty(possible_next_node)
%             for i=1:length(possible_next_node)
%                 path_queue(end+1) = {[path_queue{1},possible_next_node(i)]};
%             end
%         end        
%     end
%     path_queue(1)=[];
% end
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dist, previous_vertex, shortest_paths] = Dijkstra(graph, source)
% % [dist, previous_vertex] = Dijkstra(graph, source)
% %     This function finds the shortest graph distance from a source node to every other node in the graph      
% %     
% %     graph is an adjacency matrix, doesn't matter whether weighted or not
% %     source is the starting node
% %     dist is the approximated geodesics, length of the shortest path on graph 
% %     previous_vertex is: to reach this point via shortest path, what is the previous point before arriving here
% % 
% % http://en.wikipedia.org/wiki/Dijkstra's_algorithm
% 
% graph = graph - diag(diag(graph));
% 
% % initialize algorithm
% N = size(graph,1);
% dist = zeros(1,N) + Inf;
% previous_vertex = zeros(1,N);
% 
% % handle the source point
% dist(source)=0;
% previous_vertex(source) = source;
% 
% Q = 1:N; % All nodes in the graph are unoptimized - thus are in Q
% 
% while ~isempty(Q)  % The main loop
%     [min_dist, u] = min(dist(Q)); u = Q(u); % u is the vertex in Q with smallest distance in dist[] ;
%     if isinf(min_dist) 
%         break;  % all remaining vertices are inaccessible from source
%     end
%     Q = setdiff(Q,u);
%     for v = intersect(Q(:)', find(graph(u,:)~=0))
%         alt = dist(u) + graph(u,v);
%         if alt<dist(v)
%             dist(v) = alt;
%             previous_vertex(v) = u;
%         end
%     end
% end
% 
% for i=1:N
%     shortest_paths{i,1} = i;
% end
% 
% for target=1:N
%     while shortest_paths{target}(1)~=source && previous_vertex(shortest_paths{target}(1))~=0
%         shortest_paths{target} = [shortest_paths{previous_vertex(shortest_paths{target}(1))},shortest_paths{target}];
%         if previous_vertex(shortest_paths{target}(1))==0
%             break;
%         end
%     end
% end
