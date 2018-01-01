addpath src
addpath lib
foldername = strrep('d:/GLaMST/Simulated/sim#','#',cellstr(num2str((1:9)'))');
nruns = 500;
directed_adj = [];
stats = cell(9,5);
for ind_fd = 1:length(foldername)
    db = 0;
    dc = 0;
    dd = 0;
    de = 0;
    lsa = zeros(500,12);
    lsb = zeros(500,12);
    lsc = zeros(500,12);
    lsd = zeros(500,12);
    lse = zeros(500,12);
    win = zeros(4,3);
    for i = 1:nruns
        prefix = [foldername{ind_fd}, '/',num2str(i,'%05d')];
        [a, b, c, d, e] = prepareadj(prefix);
        db = db + abs(size(b,1)-size(a,1));
        dc = dc + abs(size(c,1)-size(a,1));
        dd = dc + abs(size(d,1)-size(a,1));
        de = de + abs(size(e,1)-size(a,1));
        if size(b,1) > size(a,1)
            win(1,1) = win(1,1) + 1;
        elseif size(b,1) < size(a,1)
            win(1,2) = win(1,2) + 1;
        else
            win(1,3) = win(1,3) + 1;
        end
        if size(c,1) > size(a,1)
            win(2,1) = win(2,1) + 1;
        elseif size(c,1) < size(a,1)
            win(2,2) = win(2,2) + 1;
        else
            win(2,3) = win(2,3) + 1;
        end
        if size(d,1) > size(a,1)
            win(3,1) = win(3,1) + 1;
        elseif size(d,1) < size(a,1)
            win(3,2) = win(3,2) + 1;
        else
            win(3,3) = win(3,3) + 1;
        end
        if size(e,1) > size(a,1)
            win(4,1) = win(4,1) + 1;
        elseif size(e,1) < size(a,1)
            win(4,2) = win(4,2) + 1;
        else
            win(4,3) = win(4,3) + 1;
        end
        lsa(i, :) = cell2mat(struct2cell(treestats(a)))';
        lsb(i, :) = cell2mat(struct2cell(treestats(b)))';
        lsc(i, :) = cell2mat(struct2cell(treestats(c)))';
        lsd(i, :) = cell2mat(struct2cell(treestats(d)))';
        lse(i, :) = cell2mat(struct2cell(treestats(e)))';
        %fprintf('Prefix: %s\t%d\t%d\t%d\n',prefix, size(a, 1), size(b, 1), size(c, 1));
        %fprintf('Sim%d_%05d:%3d\t%3d\t%3d\t%3d\t%3d\n',ind_fd,i, size(a, 1), size(b, 1), size(c, 1),size(d, 1), se)
    end
    fprintf('Sim%d: DIgtree %5d\tDPeng %5d\t DNWPeng%5d\tDPhylip%5d\n',ind_fd, db, dc, dd,de);
    disp(win)
    stats{ind_fd, 1} = lsa;
    stats{ind_fd, 2} = lsb;
    stats{ind_fd, 3} = lsc;
    stats{ind_fd, 4} = lsd;
    stats{ind_fd, 5} = lse;
end
save('data/simulation.mat','stats');
%% Generate nrmse for all features on each simulation setting
load('data/simulation.mat','stats');
output = zeros(9, 36);
for i = 1:9
    output(i, (0:11)*3+1) = mean((stats{i,2}-stats{i,1}).^2,'omitnan');
    output(i, (0:11)*3+3) = mean((stats{i,3}-stats{i,1}).^2,'omitnan');
    output(i, (0:11)*3+2) = mean((stats{i,5}-stats{i,1}).^2,'omitnan');
%     output(i, (0:11)*2+1) = mean(stats{i,2}-stats{i,1},'omitnan');
%     output(i, (1:12)*2) = mean(stats{i,3}-stats{i,1},'omitnan');
end
%% Gengerate 9-by-12 figures that compare distribution of tree feature difference.
close all
titles = {'Tree size', 'OD Root', 'OD Avg', 'OD Ratio',...
    'T', 'PL Min', 'PL Avg', 'DRSN Min', 'DLFSN Min', ...
    'DLFSN Avg', 'DASN Min', 'DASN Avg'};
%% Compare tree features
close all
figure('pos',[100 50 1200 900]);
selected = 1:12;
k = [1 5 7 10 6 8 9 12 2 3 4 11];
for j=1:12%1:length(selected)
    i = selected(k(j));
    subplot(3,4,j)
    hold on;
    plot(output(:,(3*(i-1)+1)), '-s','LineWidth',3,'Color', [0.5 0.5 0.5]);
    plot(output(:,(3*(i-1)+2)), '-.^','LineWidth',3,'Color', [0.75 0.75 0.75]);
    plot(output(:,(3*(i-1)+3)), '-ok','LineWidth',2);
    title(sprintf('%s',titles{i}));
    if j == 1
        legend({'igtree','dnapars','GLaMST'},'Location','northwest')
    end
    axis tight;
    xlabel('Simulation')
    ylabel('MSE')
    xticks(1:9);
end
