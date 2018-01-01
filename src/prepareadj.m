function [realadj, Igtreeadj, Pengadj, PengNWadj, phyadj, physize] = prepareadj(prefix)
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
%% Peng reconstructed tree without rewire
filename = [prefix,'.norewire.out.mat'];
load(filename,'reconstructed_directed_adj'); 
PengNWadj = sparse(reconstructed_directed_adj);
%% Phylip reconstructed tree
filename = [prefix,'.phy.out.edges'];
fid = fopen(filename);
phyadj = sparse([]);
while ~feof(fid)
    newline = fgetl(fid);
    cells = split(newline);
    numbers = cell2mat(cellfun(@str2num, cells(1:end), 'UniformOutput',false));
    phyadj(numbers(1)+1,numbers(2:end)+1) = 1;
end
fclose(fid);
phyadj(max(size(phyadj)),max(size(phyadj))) = 0;



% 
% 
% figure;
% subplot(2,2,1)
% treeplotA(realadj);
% title('Real')
% subplot(2,2,2)
% treeplotA(Pengadj);
% title('Peng')
% subplot(2,2,3)
% treeplotA(Igtreeadj); 
% title('igtree')
% subplot(2,2,4)
% treeplotA(phyadj); 
% title('phylip')
