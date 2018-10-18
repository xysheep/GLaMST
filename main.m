function main(fastafilename, varargin)
pnames = {'rewire'};
dflts = {1};
[rewire] = internal.stats.parseArgs(pnames,dflts,varargin{:});
[h, s] = fastaread(fastafilename);
[reconstructed_nodes,~,observed, reconstructed_directed_adj] = reconstruct_tree_minimun_tree_size([s(1);unique(s(2:end))'], rewire);

fileID = fopen([fastafilename '.out.tree'],'w');
fastaID = fopen([fastafilename '.out.fa'],'w');
k_obs = 1;
k_ins = 1;

findcellarray = @(query) cellfun(@(x) strcmp(x,query), s);
relation = cell2mat(cellfun(@(q) findcellarray(q), ...
    reconstructed_nodes(observed==1), 'UniformOutput',false));
observed_idx = find(observed==1);

node_id = cell(size(reconstructed_nodes));
for i = 1:length(reconstructed_nodes)
    output = [i find(reconstructed_directed_adj(i,:))];
    if length(output)>1
        fprintf(fileID, '%s\n',num2str(output));
    end
    
    if i==1
        originnames = ['"', strjoin(h(relation(observed_idx==i,:)),'"\t"'), '"'];
        fprintf(fastaID, '>G.L.\t%s\n',originnames);
        node_id{i} = 'G.L.';
    elseif observed(i)
        originnames = ['"', strjoin(h(relation(observed_idx==i,:)),'"\t"'), '"'];
        fprintf(fastaID, '>Obs%d\t%s\n',k_obs,originnames);
        k_obs = k_obs + 1;
        node_id{i} = sprintf('Obs%d',k_obs);
    elseif ~observed(i)
        fprintf(fastaID, '>index_%d\n',k_ins);
        k_ins = k_ins + 1;
        node_id{i} = sprintf('index_%d',k_ins);
    end
    fprintf(fastaID, '%s\n', reconstructed_nodes{i});
end
fclose(fileID);
fclose(fastaID);
    
fnewick = fopen([fastafilename '.out.newick'],'w');
fprintf(fnewick, newick(reconstructed_directed_adj, node_id, 1, 0));
fclose(fnewick);