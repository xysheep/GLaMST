function str_newick = newick(adj, nodeid, rootid,dist2p)
if length(nodeid) ~= size(adj,1)
    nodeid = cellfun(@(x)num2str(x), num2cell(1:size(adj,1)),'UniformOutput',false);
end
if sum(adj(rootid,:))==0
    str_newick = [nodeid{rootid},':',num2str(dist2p+1)];
elseif sum(adj(rootid,:))==1
    str_newick = newick(adj, nodeid, find(adj(rootid,:)==1), dist2p+1);
else
    children = num2cell(find(adj(rootid,:)==1));
    subnodes = cellfun(@(ri)newick(adj, nodeid, ri, 0),children, 'UniformOutput',false); 
    str_newick = ['(',strjoin(subnodes,','),')'];
end
