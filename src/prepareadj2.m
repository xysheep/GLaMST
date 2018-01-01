function [Igtreeadj, Pengadj, phyadj, Igtreeselected, Pengselected, phyadjselected] = prepareadj2(prefix)
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
Igtreeselected = zeros(size(Igtreeadj,1), 1);
Igtreeselected(csvread('data/real.out.id')+1) = 1;
%% Peng reconstructed tree
filename = [prefix,'.out.mat'];
load(filename,'reconstructed_directed_adj', 'reconstructed_is_selected'); 
Pengadj = sparse(reconstructed_directed_adj);
Pengselected = reconstructed_is_selected;
%% dnapars
filename = [prefix,'.dnapars.tree'];
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
phyadjselected = zeros(size(phyadj,1), 1);
phyadjselected(csvread('data/real.dnapars.id')+1) = 1;
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
