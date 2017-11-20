function [realadj, Igtreeadj, Pengadj] = prepareadj(prefix)
%% Igtree results
filename = [prefix,'.pir.out.tree'];
fid = fopen(filename);
Igtreeadj = sparse([]);
while ~feof(fid)
    newline = fgetl(fid);
    cells = split(newline);
    numbers = cell2mat(cellfun(@str2num, cells(1:end-1), 'UniformOutput',false));
    Igtreeadj(numbers(1)+1,numbers(2:end)+1) = 1;
end
fclose(fid);
Igtreeadj(max(size(Igtreeadj)),max(size(Igtreeadj))) = 0;
%% Real tree 
filename = [prefix,'.mat'];
load(filename,'directed_adj'); 
realadj = sparse(directed_adj);
%% Peng reconstructed tree
filename = [prefix,'.out.mat'];
load(filename,'reconstructed_directed_adj'); 
Pengadj = sparse(reconstructed_directed_adj);

% figure;
% subplot(1,3,1)
% treeplotA(realadj);
% title('Real')
% subplot(1,3,2)
% treeplotA(Pengadj);
% title('Peng')
% subplot(1,3,3)
% treeplotA(Igtreeadj); 
% title('igtree')
