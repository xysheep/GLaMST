%% Generate the data
addpath code
addpath read_write
op = {[1, 0, 0],[0.98, 0.01, 0.01],[0.9, 0.05, 0.05]};
op = [op op op];
foldername = strrep('sim#','#',cellstr(num2str((1:length(op))'))');
rng(9464);
len = kron([300 80 20],[1 1 1]);
nnodes = kron([150 80 40],[1 1 1]);
nruns = kron([500 500 500],[1 1 1]);
for i = 1:length(op)
    newfolder = ['data/Simulated/',foldername{i}];
    if exist(newfolder,'dir')~=7
        mkdir(newfolder);
    end
    for k=1:nruns(i)
        file_ind = k;
        filename = ['data/Simulated/',foldername{i},'/',num2str(file_ind,'%05d'),'.mat'];
        %if exist(filename,'file')~=2
            simulate_data(filename,len(i),[],[],op{i});
        %end
        fastafilename = [filename(1:end-3),'fasta'];
        if exist(fastafilename,'file')~=2
            data = [];
            load(filename,'observed_sequences')
            for t = 1:length(observed_sequences)
                if t==1
                    data(t).Sequence = observed_sequences{t};
                    data(t).Header = 'G.L.';
                else
                    data(t).Sequence = observed_sequences{t};
                    data(t).Header = num2str(t, 'Obs%d');
                end
            end
            fastawrite(fastafilename, data);
        end  
    end
end

%% Run the analysis
for i = 1:length(op)
    for k=1:nruns(i)
        file_ind = k;
        filename = ['data/Simulated/',foldername{i},'/',num2str(file_ind,'%05d'),'.mat'];
        load(filename,'nodes','observed_sequences');
        outfilename = [filename(1:end-3) 'out.mat'];
%         if exist(outfilename,'file') ~=2
            [reconstructed_nodes,reconstructed_directed_adj,reconstructed_is_selected] = reconstruct_tree_minimun_tree_size(observed_sequences);
            save(outfilename,'reconstructed_nodes','reconstructed_directed_adj','reconstructed_is_selected')
%         else
%             load(outfilename,'reconstructed_nodes');
%         end
        
    end
end


% %% Merge Igtree results
% for i = 1:length(op)
%     igtreefilename = ['sim',num2str(i),'.table.csv'];
%     igdata = csvread(igtreefilename,1,3);  
%     tree_size{i}(:,3) = igdata;
% end
% 
% stats = zeros(length(op),3);
% stats1 = zeros(length(op),3);
% for i = 1:length(op)
%     truth = tree_size{i}(:,1);
%     igtree = tree_size{i}(:,3);
%     stats(i,:) = [sum(truth<igtree) sum(truth==igtree) sum(truth>igtree)];
%     peng = tree_size{i}(:,2);
%     stats1(i,:) = [sum(truth<peng) sum(truth==peng) sum(truth>peng)];
% end

% 
% csvwrite('pengtree.csv',treesize);
% 
% 
% 
% 
% d= csvread('summary.csv',1,1);
% metric = zeros(5,3);
% metric(1,1) = sum(d(:,2)> d(:,1));
% metric(2,1) = sum(d(:,2)< d(:,1));
% metric(3,1) = sum(d(:,2)==d(:,1));
% metric(1,2) = sum(d(:,3)> d(:,1));
% metric(2,2) = sum(d(:,3)< d(:,1));
% metric(3,2) = sum(d(:,3)==d(:,1));
% metric(1,3) = sum(d(:,3)> d(:,1) & d(:,2)> d(:,1));
% metric(2,3) = sum(d(:,3)< d(:,1) & d(:,2)< d(:,1));
% metric(3,3) = sum(d(:,3)==d(:,1) & d(:,2)==d(:,1));
% metric(4,1) = sum(d(:,2)< d(:,1) & d(:,3)>=d(:,1));
% metric(4,2) = sum(d(:,3)< d(:,1) & d(:,2)>=d(:,1));
% metric(5,1) = sum(d(:,2)==d(:,1) & d(:,3)> d(:,1));
% metric(5,2) = sum(d(:,3)==d(:,1) & d(:,2)> d(:,1));