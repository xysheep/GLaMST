function Write_tree_toFile(reconstructed_nodes,reconstructed_directed_adj, reconstructed_is_selected, outfilename)
fileID = fopen([outfilename '.out.tree'],'w');
fastaID = fopen([outfilename '.out.fa'],'w');
k_obs = 1;
k_ins = 1;
for i = 1:length(reconstructed_nodes)
    output = [i find(reconstructed_directed_adj(i,:))];
    if length(output)>1
        fprintf(fileID, '%s\n',num2str(output));
    end
    if i==1
        fprintf(fastaID, '>G.L.\n');
    elseif reconstructed_is_selected(i)
        fprintf(fastaID, '>Obs%d\n',k_obs);
        k_obs = k_obs + 1;
    elseif ~reconstructed_is_selected(i)
        fprintf(fastaID, '>index_%d\n',k_ins);
        k_ins = k_ins + 1;
    end
    fprintf(fastaID, '%s\n', reconstructed_nodes{i});
end
fclose(fileID);
fclose(fastaID);