function [s] = treestats(adj)
% adj is a N-by-N matrix where N is number of nodes. adj(i,j)=1 means j is
% one child of i. 
s = struct;
s.SIZE = size(adj,1);
s.OD_Root = sum(adj(1,:));
s.OD_Avg = mean(sum(adj,2));
s.OD_Ratio = s.OD_Root./mean(sum(adj(2:end,:),2));
%% Find distance from all nodes to root
n = size(adj, 1);
depth = zeros(n, 1);
queue = 1;
while ~isempty(queue)
    idx = queue(1);
    queue = [queue(2:end) find(adj(idx,:))];
    if idx ~= 1
        depth(idx) = depth(adj(:,idx)==1) + 1;
    end
end
max_depth = max(depth);
if sum(depth==max_depth)>1
    s.T = max_depth;
else
    s.T = max_depth;
end
%%  Find all leaves 
leaf_idx = find(sum(adj,2)==0)';
PLs = depth(leaf_idx);
s.PL_Min = min(PLs);
s.PL_Avg = mean(PLs);
%% Find the first split nodes
split_idx = find(sum(adj, 2) > 1)';
if isempty(split_idx)
    s.DRSN_Min = nan;
    s.DLFSN_Min = nan;
    s.DLFSN_Avg = nan;
    s.DASN_Avg = nan;
    s.DASN_Min = nan;
    return
end
[~, i ] = min(depth(split_idx));
s.DRSN_Min = depth(split_idx(i));
s.DLFSN_Min = s.PL_Min - s.DRSN_Min;
s.DLFSN_Avg = s.PL_Avg - s.DRSN_Min;

%% Find all split nodes
if length(split_idx) == 1
    s.DASN_Avg = nan;
    s.DASN_Min = nan;
    return
end
pair_SN_dists = pdist2(depth(split_idx), depth(split_idx), 'cityblock');
s.DASN_Min = min(min(pair_SN_dists(pair_SN_dists>0)));
queue = 1;
SNqueue = 0;
dist2SN = zeros(n,1);
while ~isempty(queue)
    idx = queue(1);
    queue = [queue(2:end) find(adj(idx,:))];
    dist2SN(idx) = SNqueue(1);
    if ismember(idx, split_idx)
        SNqueue = [SNqueue(2:end) ones(size(find(adj(idx,:))))];
    elseif SNqueue(1) ~= 0
        SNqueue = [SNqueue(2:end) SNqueue(1)+1];
    else
        SNqueue = [SNqueue(2:end) SNqueue(1)];
    end
end
SN_dist = dist2SN(split_idx);
s.DASN_Avg = mean(SN_dist);




