addpath src
addpath lib
foldername = strrep('data/Simulated/sim#','#',cellstr(num2str((1:9)'))');
nruns = 500;
directed_adj = [];
stats = cell(9,4);
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
        [a, b, c, d, e, se] = prepareadj(prefix);
        db = db + abs(size(b,1)-size(a,1));
        dc = dc + abs(size(c,1)-size(a,1));
        dd = dc + abs(size(d,1)-size(a,1));
        de = de + abs(se-size(a,1));
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
        if se > size(a,1)
            win(4,1) = win(4,1) + 1;
        elseif se < size(a,1)
            win(4,2) = win(4,2) + 1;
        else
            win(4,3) = win(4,3) + 1;
        end
        lsa(i, :) = cell2mat(struct2cell(treestats(a)))';
        lsb(i, :) = cell2mat(struct2cell(treestats(b)))';
        lsc(i, :) = cell2mat(struct2cell(treestats(c)))';
        lsd(i, :) = cell2mat(struct2cell(treestats(d)))';
        %lse(i, :) = cell2mat(struct2cell(treestats(e)))';
        %fprintf('Prefix: %s\t%d\t%d\t%d\n',prefix, size(a, 1), size(b, 1), size(c, 1));
        %fprintf('Sim%d_%05d:%3d\t%3d\t%3d\t%3d\t%3d\n',ind_fd,i, size(a, 1), size(b, 1), size(c, 1),size(d, 1), se)
    end
    fprintf('Sim%d: DIgtree %5d\tDPeng %5d\t DNWPeng%5d\tDPhylip%5d\n',ind_fd, db, dc, dd,de);
    disp(win)
    stats{ind_fd, 1} = lsa;
    stats{ind_fd, 2} = lsb;
    stats{ind_fd, 3} = lsc;
    stats{ind_fd, 4} = lsd;
    %stats{ind_fd, 5} = lse;
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
%% hist plot, bad visualization
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
%% boxplot
close all
for i = 1:9
    figure('Name',num2str(i),'pos', [200 200 1400 600])
    subplot(1,2,1)
    boxplot((stats{i,3}-stats{i,1})./abs(prctile(stats{i,1}, 99)-prctile(stats{i,1}, 1)));
    if i <9 
        axis([0 13 -1.5 1.5])
    else 
        axis([0 13 -2.2 2.2])
    end
    xticks(1:12)
    xticklabels(titles)
    xtickangle(45)
    title(sprintf('Diff Peng and Real (sim%d)',i))
    subplot(1,2,2)
    boxplot((stats{i,2}-stats{i,1})./abs(prctile(stats{i,1}, 99)-prctile(stats{i,1}, 1)));
    if i <9 
        axis([0 13 -1.5 1.5])
    else 
        axis([0 13 -2.2 2.2])
    end
    title(sprintf('Diff igtree and Real (sim%d)',i))
    xticks(1:12)
    xticklabels(titles)
    xtickangle(45)
    print(sprintf('manuscript/figures/FeatureDiffSim%d.jpeg',i),'-djpeg')
end
%% Density curve on tree size
close all
for j = 1:12
    figure('Name',num2str(i),'pos', [200 200 1400 600]);
    for i = 1:9
        subplot(3,3,i);
        hold on
        h2 = histogram(stats{i, 3}(:,j)-stats{i, 1}(:, j));
        h2.FaceAlpha = 0.5;
        h1 = histogram(stats{i, 2}(:,j)-stats{i, 1}(:, j));
        h1.FaceAlpha = 0.5;
        bwidth = min(h1.BinWidth, h2.BinWidth);
        h1.BinWidth = bwidth;
        h2.BinWidth = bwidth;
        h2.BinLimits = h1.BinLimits;
        % plot([0 0],[0 0.26])
        % axis([-20 20 0 0.26])
        xlabel(sprintf('%s difference to real tree', titles{j}));
        ylabel('Frequency');
        hold off
        legend('Peng','Igtree');
    end
    print(sprintf('manuscript/figures/FeatureHist%d.jpeg',j),'-djpeg')
end
%% Scatter plot for each feature

figure;
for i_ft = 1:12
    subplot(3, 4, i_ft);
    hold on;

    c = kron([1 1 1], [1 2 3]);
    shape = {'d','d','d','^','^','^','o','o','o'};
    colormap jet;
    rng(9464);
    for i = 1:9
        d2 = abs(stats{i, 2}(:,i_ft)-stats{i, 1}(:, i_ft)) + rand(500, 1)/5 * i;
        d3 = abs(stats{i, 3}(:,i_ft)-stats{i, 1}(:, i_ft)) + rand(500, 1)/5 * i;
        s = scatter(d2, d3, 25 * ones(size(d2)), c(i * ones(size(d2))), 'filled');
        s.Marker = shape{i};
    end
    plot([-0.5 9999], [-0.5 9999]);
    xlabel('Igtree');
    ylabel('Peng');
    title('Tree feature difference from real tree');
    legend({'1','2','3','4','5','6','7','8','9'})
    axis([-0.5 max([d2;d3])+0.5 -0.5 max([d2;d3])+0.5])
end

%% Grouped boxplot for each feature
for i_ft = 1:12
    subplot(3,4,i_ft)
    d = [];
    sims = [];
    algs = [];
    for i = [3 6 9]
        d2 = abs(stats{i, 2}(:,i_ft)-stats{i, 1}(:, i_ft));
        d3 = abs(stats{i, 3}(:,i_ft)-stats{i, 1}(:, i_ft));
        d = [d;d2;d3];
        sims = [sims; i * ones(size([d2;d3]))];
        algs = [algs; ones(size(d2)); 2 * ones(size(d3))];
    end
    boxplot(d, {sims, algs},'factorgap',[10,0]);
    title(titles{i_ft})
end

