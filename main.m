function main(fastafilename, varargin)
pnames = {'rewire'};
dflts = {1};
[rewire] = internal.stats.parseArgs(pnames,dflts,varargin{:});
[~, s] = fastaread(fastafilename);
[reconstructed_nodes,~,observed, reconstructed_directed_adj] = reconstruct_tree_minimun_tree_size([s(1);unique(s(2:end))'], rewire);
fileID = fopen([fastafilename '.out.tree'],'w');
fastaID = fopen([fastafilename '.out.fa'],'w');
k_obs = 1;
k_ins = 1;
for i = 1:length(reconstructed_nodes)
    output = [i find(reconstructed_directed_adj(i,:))];
    if length(output)>1
        fprintf(fileID, '%s\n',num2str(output));
    end
    if i==1
        fprintf(fastaID, '>G.L.\n');
    elseif observed(i)
        fprintf(fastaID, '>Obs%d\n',k_obs);
        k_obs = k_obs + 1;
    elseif ~observed(i)
        fprintf(fastaID, '>index_%d\n',k_ins);
        k_ins = k_ins + 1;
    end
    fprintf(fastaID, '%s\n', reconstructed_nodes{i});
end
fclose(fileID);
fclose(fastaID);
    