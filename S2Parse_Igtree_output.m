addpath code
addpath read_write
foldername = strrep('data/Simulated/sim#','#',cellstr(num2str((1:9)'))');
nruns = 500;
directed_adj = [];
stats = cell(9,3);
for ind_fd = 1:length(foldername)
    db = 0;
    dc = 0;
    lsa = zeros(500,12);
    lsb = zeros(500,12);
    lsc = zeros(500,12);
    win = [0 0 0 0 0 0];
    for i = 1:nruns
        prefix = [foldername{ind_fd}, '/',num2str(i,'%05d')];
        [a, b, c] = prepareadj(prefix);
        db = db + abs(size(b,1)-size(a,1));
        dc = dc + abs(size(c,1)-size(a,1));
        if size(b,1) > size(a,1)
            win(1) = win(1) + 1;
        elseif size(b,1) < size(a,1)
            win(2) = win(2) + 1;
        else
            win(3) = win(3) + 1;
        end
        if size(c,1) > size(a,1)
            win(4) = win(4) + 1;
        elseif size(c,1) < size(a,1)
            win(5) = win(5) + 1;
        else
            win(6) = win(6) + 1;
        end
        lsa(i, :) = cell2mat(struct2cell(treestats(a)))';
        lsb(i, :) = cell2mat(struct2cell(treestats(b)))';
        lsc(i, :) = cell2mat(struct2cell(treestats(c)))';
        %fprintf('Prefix: %s\t%d\t%d\t%d\n',prefix, size(a, 1), size(b, 1), size(c, 1));
    end
    fprintf('Sim%d: DIgtree %5d\tDPeng %5d  %5d\t%5d\t%5d%5d\t%5d\t%5d\n',ind_fd, db, dc, win); 
    stats{ind_fd, 1} = lsa;
    stats{ind_fd, 2} = lsb;
    stats{ind_fd, 3} = lsc;
    break
end
save('data/simulation.mat','stats');
%% Generate nrmse for all features on each simulation setting
load('data/simulation.mat','stats');
output = zeros(9, 24);
for i = 1:9
    output(i, (0:11)*2+1) = mean((stats{i,2}-stats{i,1}).^2,'omitnan');
    output(i, (1:12)*2) = mean((stats{i,3}-stats{i,1}).^2,'omitnan');
%     output(i, (0:11)*2+1) = mean(stats{i,2}-stats{i,1},'omitnan');
%     output(i, (1:12)*2) = mean(stats{i,3}-stats{i,1},'omitnan');
end
%% Gengerate 9-by-12 figures that compare distribution of tree feature difference.
close all
titles = {'Tree size', 'OD Root', 'OD Avg', 'OD Ratio',...
    'T', 'PL Min', 'PL Avg', 'DRSN Min', 'DLFSN Min', ...
    'DLFSN Avg', 'DASN Min', 'DASN Avg'};	
figs = cell(9,1);
for i = 1:9
    figs{i} = figure('pos',[100 100 1200 800]);
    for j = 1:12
        subplot(3,4,j);
        hold on
%         ksdensity(stats{i,2}(:,j)-stats{i,1}(:,j));
%         ksdensity(stats{i,3}(:,j)-stats{i,1}(:,j));
        hist(log2(abs(stats{i,3}(:,j)-stats{i,1}(:,j))./abs(stats{i,2}(:,j)-stats{i,1}(:,j))), 20);
        hold off
        title(titles{j})
    end
    %legend('Igtree', 'Peng')
end
figure('pos', [200 200 1400 600])
for i = [3 6 9]
    subplot(2,3,3+i/3)
    boxplot((stats{i,3}-stats{i,1})./abs(prctile(stats{i,1}, 99)-prctile(stats{i,1}, 1)));
    axis([0 13 -1.5 1.5])
    xticks(1:12)
    xticklabels(titles)
    xtickangle(45)
    title(sprintf('Diff Peng and Real (sim%d)',i))
    subplot(2,3,i/3)
    boxplot((stats{i,2}-stats{i,1})./abs(prctile(stats{i,1}, 99)-prctile(stats{i,1}, 1)));
    axis([0 13 -1.5 1.5])
    title(sprintf('Diff igtree and Real (sim%d)',i))
    xticks(1:12)
    xticklabels(titles)
    xtickangle(45)
end










