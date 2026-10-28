%% MAPseq: projections of superficial SC neurons (mice M1-M3)
% Analyses and figures of the MAPseq data. Run mapseq_preprocess.m first.
% Sections are meant to be run in order (Run Section): later sections use
% variables computed in earlier ones.
%
% Data (from mapseq_preprocessed.mat; columns = samples, see sample_names)
%   B1, B2, B3, BigB              molecule counts (per mouse / all mice)
%   B1norm, B2norm, B3norm        spike-normalised, per mouse
%   BigBnorm                      BigB / spike-in molecules summed over mice
%   m1, m2, m3, maxes             soma sample of each barcode
%   sSC                           ipsilateral sSC samples (column vector)
%   projTargets                   target columns analysed (see targetSet)
%
% Requires on the path: cbrewer, hex2rgb, customcolormap, magma,
% alluvialflow, mysigstar (see README) and colours from color_schemes.mat.

clearvars; close all;

P = setup_paths();
load(fullfile(P.processed, 'mapseq_preprocessed.mat'));
load(fullfile(P.resources, 'color_schemes.mat'));   % col, col_bias_med_merge, col_bias_lat_merge, div_cmap3, div_cmap_grp3

% Targets analysed: longIpsi, longContra, longIpsiContra, dSCipsi,
% dSCcontra, dSCipsiContra or longIpsi_dSCipsi (see mapseq_preprocess.m)
targetSet = 'longIpsi';
projTargets = targetSets.(targetSet);

set(0, 'defaultAxesFontSize', 10);

%% count distribution in target areas (all mice matrices in one)

id=sSC;

%make big matrix of all mice and max expressions
% (maxes, the soma sample of each barcode, comes from mapseq_preprocess.m)

BCodes = ismember(maxes,id);
data=BigB(BCodes,:);
data_proj=data(:,projTargets);
data_all=reshape(data_proj,[],1); data_all=data_all(data_all>0);
figure(1)
histogram(data_all)
set(gca,'xtick', 0.5:30.5, 'xticklabels', 0:30, 'box', 'off', 'TickDir', 'out', 'TickLength', [0 0])
ylabel('# barcode counts')
xlabel('Barcode count in targets')
xlim([0.7 25])

%% projection patterns of unique barcodes

DATA=BigBnorm;

BCodes = find(ismember(maxes,sSC));
data=BigB(BCodes,:);
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;

threshold = 0;
id_all=sSC;
id_med=[3 7 11]; %medial SC
id_lat=[2 6 10]; %lateral SC
id_tot=[id_med;id_lat];
names_medlat={'medial SCs','lateral SCs'};

scnt = 0;
Numbers = zeros(length(sSC),2); allareas = [];
for s = 1:size(id_all,2) % go through medial and lateral injection sites
    scnt = scnt +1;
    % id = sSC(s);
    % id = id_tot(s,:);
    id=id_all;
    BCodes = find(ismember(maxes,id));
    Numbers(scnt,1) = length(BCodes);
    % for each barcode find the # of areas
    areas = [];
    for b = 1:length(BCodes)
        curr = DATA(BCodes(b),:);
        areas = [areas; curr];
    end

    %     figure

    allareas = [allareas; areas];


    %filter projection targets, and remove rows of all 0
    targets_bc=areas(:, projTargets);
    zero_rows=all(targets_bc==0, 2);
    targets_bc=targets_bc(~zero_rows, :);

    %normalize values to sum of each row
    norm_vals=sum(targets_bc,2);
    targets_bc=targets_bc./norm_vals;

    %sort the array
    [~,loc]=max(targets_bc,[],2);
    [~,ix]=sort(loc);


% %     %number of barcode counts
%     targetsMatrix = zeros(length(projTargets),length(projTargets));
%     targetsMatrix2 = zeros(length(projTargets),1);
%     for pr = 1:length(projTargets)
%         tmp = find(keepData(:,projTargets(pr))>0);
%         targetsMatrix2(pr,1) = length(tmp);
%         for pr2 = pr:length(projTargets)
%             tmp2 = find(keepData(tmp,projTargets(pr2))>0);
%             targetsMatrix(pr,pr2) = 1/length(tmp)*length(tmp2);
%             targetsMatrix(pr2,pr) = 1/length(tmp)*length(tmp2);
%
%         end
%     end

    % figure(100)
    % subplot(3,3,scnt)
    % imagesc(targets_bc(ix,:));
    % cols=cbrewer('seq', 'Reds', 300);
    % cols(1,:)=[1 1 1]; %make 0 values value white
    % cols=cols(1:280,:); %make max values less dark
    % colormap(cols)
    % set(gca,'xticklabel',sample_names(projTargets), 'xtick', 1:5,'ytickLabels',[],...
    %     'ytick', [], 'FontSize', 7)
    % xtickangle(45);
    % xlabel('Projection targets', 'FontSize',8)
    % % ylabel('Area B', 'FontSize',10)
    % cb=colorbar;
    % ylabel(cb, 'Projection strength');
    % name_str= strjoin(SampleName{3}(sSC(s)), ', ');
    % title(['Projection patterns of ' name_str ' neurons'], 'FontSize', 10)


    %fraction of neurons innervating a single area out of all neurons that
    %innervate that area
    nonzero_counts=zeros(1, size(targets_bc, 2)); %number of neurons projecting to each area

    for i = 1:size(targets_bc, 2)
        nonzero_counts(i) = nnz(targets_bc(:, i));
    end
    ones_counts = zeros(1, size(targets_bc, 2)); %number of neurons that project only one particular area
    for i = 1:size(targets_bc, 2)
        ones_counts(i)=sum(targets_bc(:, i) == 1);
    end
    frac_ones=ones_counts./nonzero_counts;

    figure(s); set(gcf, 'position', [286.6000 242.6000 305.6000 272])
    bar(frac_ones);
    xlim([0.2 5.8]);
    xticks(1:5);
    set(gca,'xticklabel',sample_names(projTargets), 'box', 'off', 'TickDir', 'out')
    xtickangle(0)
    ylabel('fraction of dedicated input')
    % title(['Dedicated inputs of ', names_medlat{s}, ' neurons'], 'FontSize', 11)
    title('Dedicated inputs of SCs neurons', 'FontSize', 11);

end


%------ipsi vs contra projection patterns of all targets
figure(3)
all_sorted=allareas(:,projTargets);

%filter projection targets, and remove rows of all 0
zero_rows=all(all_sorted==0, 2);
all_sorted=all_sorted(~zero_rows, :);

%normalize values to sum of each row
norm_vals=sum(all_sorted,2);
all_sorted=all_sorted./norm_vals;

%sort the array
[~,loc]=max(all_sorted,[],2);
[~,ix]=sort(loc);

%separate ipsi and contra

%plot
imagesc(all_sorted(ix,:));
cols=cbrewer('seq', 'Reds', 300);
cols(1,:)=[1 1 1]; %make 0 values value white
cols=cols(1:280,:); %make max values less dark
colormap(cols)
% colormap('gray')
set(gca,'xticklabel', sample_names(projTargets), 'xtick', 1:5,'ytickLabels',[],...
    'ytick', [], 'FontSize', 7)
xtickangle(45);
xlabel('Projection targets', 'FontSize',10)
ylabel('Barcodes', 'FontSize',10)
cb=colorbar;
ylabel(cb, 'Projection strength', 'FontSize', 8);


%% projection patterns of unique barcodes - MED LAT

coloritos=[hex2rgb("#E8AE7D");hex2rgb("#6D9CC5");hex2rgb("#D06D74");hex2rgb("#38A68A")];
DATA = BigBnorm;
% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<=cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
id_med=[3 7 11]; %medial SC
id_lat=[2 6 10]; %lateral SC
id_tot={id_med;id_lat};
names_medlat={'upper sSC','lower sSC'};

% id_tot={id_med;id_lat;14};
% names_medlat={'upper sSC','lower sSC', 'caudal sSC'};

names=sample_names(projTargets);

scnt = 0;
Numbers = zeros(length(sSC),2); allareas = []; all_targets=[];
% for s = orderOfSC % go through injection sites
for s = 1:size(id_tot,1) % go through injection sites
    scnt = scnt +1;
    id = id_tot{s};
    BCodes= find(ismember(maxes,id));
    Numbers(scnt,1) = length(BCodes);
    % for each barcode find the # of areas
    areas = [];
    for b = 1:length(BCodes)
        curr = DATA(BCodes(b),:);
        areas = [areas; curr];
    end

    %     figure

    allareas = [allareas; areas];


    %filter projection targets, and remove rows of all 0
    targets_bc=areas(:, projTargets);

    zero_rows=all(targets_bc==0, 2);
    targets_bc=targets_bc(~zero_rows, :);


    %normalize values to sum of each row
    norm_vals=sum(targets_bc,2);
    targets_bc=targets_bc./norm_vals;
    all_targets=[all_targets ; targets_bc];

    %---SORT---

    % Preallocate a sorted matrix
    sortedMatrix = [];

    % Step 1: Handle rows with a value of 1 in any column
    for i = 1:length(projTargets)
        % Extract rows where column `col` has a value of 1
        rowsWithOnes = targets_bc(targets_bc(:, i) == 1, :);

        % Append these rows to the sorted matrix
        sortedMatrix = [sortedMatrix; rowsWithOnes];
    end

    % Step 2: Handle rows with no values of 1
    remainingRows = targets_bc(~any(targets_bc == 1, 2), :);

    % Sort the remaining rows by descending values in each column (priority order: col1, col2, ...)
    [~, sortIndices] = sortrows(remainingRows, -(1:length(projTargets)));
    sortedMatrix = [sortedMatrix; remainingRows(sortIndices, :)];

    % %group dSC regions into medial lateral anterior and posterior
    % if exist('dSCmi', 'var')
    %     sortedMatrix = [sortedMatrix(:,1:5) sum(sortedMatrix(:,dSCmi),2) ...
    %         sum(sortedMatrix(:,dSCli),2) sum(sortedMatrix(:,dSCai),2) ...
    %         sum(sortedMatrix(:,dSCpi),2)];
    %     names=[sample_names(projTargets(1:5)) sample_names(end-3:end)];
    % end

    % sorted_data = data(sort_idx, :);
    % [~,loc]=max(targets_bc,[],2);
    % [~,ix]=sort(loc);
    n_neurons=size(sortedMatrix,1);
    disp(n_neurons)

    figure(1); set(gcf, 'Position',[85.8000 282.6000 1192 379.4000]);

    if s==1
        ax(1)=subplot(1,size(id_tot,1),scnt);
        % imagesc(targets_bc(ix,:));
        imagesc(sortedMatrix);
        % cols_1=cbrewer('seq', 'Oranges', 300);
        cols_1=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#FEF7EE','#F9DFC2','#F4B985', '#E79357', '#B86C44', '#764630'},300);
        cols_1=flipud(cols_1);
        cols_1(1,:)=[1 1 1];
        cols_1=cols_1(1:250,:);
        colormap(ax(1),cols_1)
    elseif s==2
        ax(2)=subplot(1,size(id_tot,1),scnt);
        imagesc(sortedMatrix);
        % cols_2=cbrewer('seq', 'Blues', 300);
        cols_2=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#F7FAFE','#D5DFEE','#9DBFD9', '#5A8EBF', '#2B599C', '#182D62'},300);
        cols_2=flipud(cols_2);
        cols_2(1,:)=[1 1 1];
        cols_2=cols_2(1:250,:);
        colormap(ax(2),cols_2)
    else
        ax(3)=subplot(1,size(id_tot,1),scnt);
        imagesc(sortedMatrix);
        cols_3=cbrewer('seq', 'Greens', 300);
        cols_3(1,:)=[1 1 1];
        colormap(ax(3),cols_3)
    end
    % cols(1,:)=[1 1 1]; %make 0 values value white
    % cols=cols(1:280,:); %make max values less dark
    % colormap(cols)


    set(gca,'xticklabel',names, 'xtick', 1:size(sortedMatrix,2),'ytickLabels',[],...
        'ytick', [], 'FontSize', 14)
    % xtickangle(90);
    % xlabel('Projection targets', 'FontSize',12)
    ylabel('Barcodes', 'FontSize',14)
    cb=colorbar;
    ylabel(cb, 'Projection strength');
    title([names_medlat{s}, ' projections (n=', num2str(n_neurons),')'], 'FontSize', 14)

    %fraction of neurons innervating a single area out of all neurons that
    %innervate that area
    figure(2); set(gcf, 'position', [445.8000 241 816 272])
    subplot(1,size(id_tot,1),scnt)
    nonzero_counts=zeros(1, size(targets_bc, 2)); %number of neurons projecting to each area

    for i = 1:size(targets_bc, 2)
        nonzero_counts(i) = nnz(targets_bc(:, i));
    end
    ones_counts = zeros(1, size(targets_bc, 2)); %number of neurons that project only one particular area
    for i = 1:size(targets_bc, 2)
        ones_counts(i)=sum(targets_bc(:, i) == 1);
    end
    frac_ones=ones_counts./nonzero_counts;

    b=bar(frac_ones, 'EdgeColor', 'none');
    b.FaceColor = 'flat';
    if s>3
       b.CData=[0.5 0.5 0.5];
    elseif s==3
        b.CData=col(s-1,:);
    else
        b.CData=coloritos(s,:);
    end
    % xlim([0.2 5.8]);
    xticks(1:length(projTargets));
    set(gca,'xticklabel',sample_names(projTargets), 'box', 'off', 'TickDir', 'out',...
        'FontSize', 11)
    xtickangle(90)
    ylim([0 1])
    ylabel('fraction of dedicated output')
    % title(['Dedicated inputs of ', names_medlat{s}, ' neurons'], 'FontSize', 11)
    title(['dedicated outputs of ', names_medlat{s}, ' neurons'], 'FontSize', 12);

end

% %fraction of neurons innervating a single area out of all neurons that
% %innervate that area
% nonzero_counts=zeros(1, size(all_targets, 2)); %number of neurons projecting to each area
%
% for i = 1:size(all_targets, 2)
%     nonzero_counts(i) = nnz(all_targets(:, i));
% end
% ones_counts = zeros(1, size(all_targets, 2)); %number of neurons that project only one particular area
% for i = 1:size(all_targets, 2)
%     ones_counts(i)=sum(all_targets(:, i) == 1);
% end
% frac_ones=ones_counts./nonzero_counts;
%

%% projection patterns of unique barcodes and individual mice - med lat

coloritos=[hex2rgb("#E8AE7D");hex2rgb("#6D9CC5");hex2rgb("#D06D74");hex2rgb("#38A68A")];
DATA = BigBnorm;
% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<=cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
id_med=[3 7 11]; %medial SC
id_lat=[2 6 10]; %lateral SC
id_tot={id_med;id_lat};
names_medlat={'upper sSC','lower sSC'};

%apply cutoff to single mice
DATA1=B1norm;
cutoff_tmp=B1(:,:)<=cutoff;
DATA1(cutoff_tmp)=0;

DATA2=B2norm;
cutoff_tmp=B2(:,:)<=cutoff;
DATA2(cutoff_tmp)=0;

DATA3=B3norm;
cutoff_tmp=B3(:,:)<=cutoff;
DATA3(cutoff_tmp)=0;


scnt = 0;
Numbers = {}; allareas = []; all_targets=[];
% for s = orderOfSC % go through injection sites
for s = 1:size(id_tot,1) % go through adults and juv
    scnt = scnt +1;
    id = id_tot{s};
    BCodes1 = find(ismember(m1,id));
    BCodes2 = find(ismember(m2,id));
    BCodes3 = find(ismember(m3,id));
    BCodes = find(ismember(maxes,id));
    Numbers{1}(scnt,1) = length(BCodes1);
    Numbers{2}(scnt,1) = length(BCodes2);
    Numbers{3}(scnt,1) = length(BCodes3);
    % for each barcode find the # of areas
    areas1 = []; areas2 = []; areas3 = []; areas=[];
    for b = 1:length(BCodes1)
        curr = DATA1(BCodes1(b),:);
        areas1 = [areas1;curr];
    end
    for b = 1:length(BCodes2)
        curr = DATA2(BCodes2(b),:);
        areas2 = [areas2;curr];
    end
    for b = 1:length(BCodes3)
        curr = DATA3(BCodes3(b),:);
        areas3 = [areas3;curr];
    end
    for b = 1:length(BCodes)
        curr = DATA(BCodes(b),:);
        areas = [areas;curr];
    end


    %filter projection targets, and remove rows of all 0
    targets_bc1=areas1(:,projTargets);
    targets_bc2=areas2(:,projTargets);
    targets_bc3=areas3(:,projTargets);
    targets_bc=areas(:, projTargets);

    zero_rows1=all(targets_bc1==0, 2);
    zero_rows2=all(targets_bc2==0, 2);
    zero_rows3=all(targets_bc3==0, 2);
    zero_rows=all(targets_bc==0, 2);

    targets_bc1=targets_bc1(~zero_rows1, :);
    targets_bc2=targets_bc2(~zero_rows2, :);
    targets_bc3=targets_bc3(~zero_rows3, :);
    targets_bc=targets_bc(~zero_rows, :);

    n_targets1 = sum(targets_bc1 > 0, 2); % number of targets per neuron
    n_targets2 = sum(targets_bc2 > 0, 2);
    n_targets3 = sum(targets_bc3 > 0, 2);

    nNeurons1 = size(targets_bc1, 1);
    nNeurons2 = size(targets_bc2, 1);
    nNeurons3 = size(targets_bc3, 1);

    % frac_by_ntargets: rows = mice, cols = [1 target, 2 targets, 3+]
    frac_by_ntargets{scnt}(1,:) = [sum(n_targets1==1)/nNeurons1, sum(n_targets1==2)/nNeurons1, sum(n_targets1>=3)/nNeurons1];
    frac_by_ntargets{scnt}(2,:) = [sum(n_targets2==1)/nNeurons2, sum(n_targets2==2)/nNeurons2, sum(n_targets2>=3)/nNeurons2];
    frac_by_ntargets{scnt}(3,:) = [sum(n_targets3==1)/nNeurons3, sum(n_targets3==2)/nNeurons3, sum(n_targets3>=3)/nNeurons3];

    nNeurons_store{scnt} = [nNeurons1; nNeurons2; nNeurons3];


    %normalize values to sum of each row
    norm_vals1=sum(targets_bc1,2);
    targets_bc1=targets_bc1./norm_vals1;
    norm_vals2=sum(targets_bc2,2);
    targets_bc2=targets_bc2./norm_vals2;
    norm_vals3=sum(targets_bc3,2);
    targets_bc3=targets_bc3./norm_vals3;

    norm_vals=sum(targets_bc,2);
    targets_bc=targets_bc./norm_vals;
    all_targets=[all_targets ; targets_bc];


    %---SORT---

    % Preallocate a sorted matrix
    sortedMatrix = [];

    % Step 1: Handle rows with a value of 1 in any column
    for i = 1:length(projTargets)
        % Extract rows where column `col` has a value of 1
        rowsWithOnes = targets_bc(targets_bc(:, i) == 1, :);

        % Append these rows to the sorted matrix
        sortedMatrix = [sortedMatrix; rowsWithOnes];
    end

    % Step 2: Handle rows with no values of 1
    remainingRows = targets_bc(~any(targets_bc == 1, 2), :);

    % Sort the remaining rows by descending values in each column (priority order: col1, col2, ...)
    [~, sortIndices] = sortrows(remainingRows, -(1:length(projTargets)));
    sortedMatrix = [sortedMatrix; remainingRows(sortIndices, :)];

    % %group dSC regions into medial lateral anterior and posterior
    % if exist('dSCmi', 'var')
    %     sortedMatrix = [sortedMatrix(:,1:5) sum(sortedMatrix(:,dSCmi),2) ...
    %         sum(sortedMatrix(:,dSCli),2) sum(sortedMatrix(:,dSCai),2) ...
    %         sum(sortedMatrix(:,dSCpi),2)];
    %     names=[sample_names(projTargets(1:5)) sample_names(end-3:end)];
    % end

    % sorted_data = data(sort_idx, :);
    % [~,loc]=max(targets_bc,[],2);
    % [~,ix]=sort(loc);
    disp(size(sortedMatrix,1))

    figure(1); set(gcf, 'Position',[85.8000 282.6000 1192 379.4000]);

    if s==1
        ax(1)=subplot(1,size(id_tot,1),scnt);
        % imagesc(targets_bc(ix,:));
        imagesc(sortedMatrix);
        % cols_1=cbrewer('seq', 'Oranges', 300);
        cols_1=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#FEF7EE','#F9DFC2','#F4B985', '#E79357', '#B86C44', '#764630'},300);
        cols_1=flipud(cols_1);
        cols_1(1,:)=[1 1 1];
        cols_1=cols_1(1:250,:);
        colormap(ax(1),cols_1)
    elseif s==2
        ax(2)=subplot(1,size(id_tot,1),scnt);
        imagesc(sortedMatrix);
        % cols_2=cbrewer('seq', 'Blues', 300);
        cols_2=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#F7FAFE','#D5DFEE','#9DBFD9', '#5A8EBF', '#2B599C', '#182D62'},300);
        cols_2=flipud(cols_2);
        cols_2(1,:)=[1 1 1];
        cols_2=cols_2(1:250,:);
        colormap(ax(2),cols_2)
    else
        ax(3)=subplot(1,size(id_tot,1),scnt);
        imagesc(sortedMatrix);
        cols_3=cbrewer('seq', 'Greens', 300);
        cols_3(1,:)=[1 1 1];
        colormap(ax(3),cols_3)
    end
    % cols(1,:)=[1 1 1]; %make 0 values value white
    % cols=cols(1:280,:); %make max values less dark
    % colormap(cols)


    set(gca,'xticklabel',names, 'xtick', 1:size(sortedMatrix,2),'ytickLabels',[],...
        'ytick', [], 'FontSize', 14)
    % xtickangle(90);
    % xlabel('Projection targets', 'FontSize',12)
    ylabel('Barcodes', 'FontSize',14)
    cb=colorbar;
    ylabel(cb, 'Projection strength');
    title([names_medlat{s}, ' projections'], 'FontSize', 14)

    %fraction of neurons innervating a single area out of all neurons that
    %innervate that area
    figure(2); set(gcf, 'position', [445.8000 241 816 272])
    subplot(1,2,scnt)
    nonzero_counts=zeros(1, size(targets_bc, 2)); %number of neurons projecting to each area

    nonzero_counts1=zeros(1, size(targets_bc1, 2));
    nonzero_counts2=zeros(1, size(targets_bc2, 2));
    nonzero_counts3=zeros(1, size(targets_bc3, 2));

    for i = 1:size(targets_bc, 2)
        nonzero_counts1(i) = nnz(targets_bc1(:, i));
        nonzero_counts2(i) = nnz(targets_bc2(:, i));
        nonzero_counts3(i) = nnz(targets_bc3(:, i));
        nonzero_counts(i) = nnz(targets_bc(:, i));
    end
    ones_counts = zeros(1, size(targets_bc, 2)); %number of neurons that project only one particular area

    ones_counts1 = zeros(1, size(targets_bc1, 2));
    ones_counts2 = zeros(1, size(targets_bc2, 2));
    ones_counts3 = zeros(1, size(targets_bc3, 2));

    for i = 1:size(targets_bc, 2)
        ones_counts1(i)=sum(targets_bc1(:, i) == 1);
        ones_counts2(i)=sum(targets_bc2(:, i) == 1);
        ones_counts3(i)=sum(targets_bc3(:, i) == 1);
        ones_counts(i)=sum(targets_bc(:, i) == 1);
    end
    frac_ones=ones_counts./nonzero_counts;

    frac_ones1=ones_counts1./nonzero_counts1;
    frac_ones2=ones_counts2./nonzero_counts2;
    frac_ones3=ones_counts3./nonzero_counts3;


    b=bar(frac_ones);
    b.FaceColor = 'flat';
    if s>3
       b.CData=[0.5 0.5 0.5];
    elseif s==3
        b.CData=col(s-1,:);
    else
        b.CData=coloritos(s,:);
    end
    b.EdgeColor='none';
    hold on;
    plot(1:length(projTargets),frac_ones1,'o','Color', [.3 .3 .3], 'MarkerSize', 5, LineWidth=1.5); hold on
    plot(1:length(projTargets),frac_ones2,'o','Color', [.3 .3 .3], 'MarkerSize', 5, LineWidth=1.5); hold on
    plot(1:length(projTargets),frac_ones3,'o','Color', [.3 .3 .3], 'MarkerSize', 5, LineWidth=1.5)


    % xlim([0.2 5.8]);
    xticks(1:length(projTargets));
    set(gca,'xticklabel',sample_names(projTargets), 'box', 'off', 'TickDir', 'out',...
        'FontSize', 11)
    xtickangle(90)
    ylim([0 1])
    ylabel('fraction of dedicated output')
    % title(['Dedicated inputs of ', names_medlat{s}, ' neurons'], 'FontSize', 11)
    title(['dedicated outputs of ', names_medlat{s}, ' neurons'], 'FontSize', 12);

end

%plot stacked bar with fraction of dedicated outputs
stackData = [frac_by_ntargets{1}(1,:);  % Mouse1 medial
             frac_by_ntargets{2}(1,:);  % Mouse1 lateral
             frac_by_ntargets{1}(2,:);  % Mouse3 medial
             frac_by_ntargets{2}(2,:);  % Mouse3 lateral
             frac_by_ntargets{1}(3,:);  % Mouse2 medial
             frac_by_ntargets{2}(3,:)]; % Mouse2 lateral

nNeuronsAll = [nNeurons_store{1}(1); nNeurons_store{2}(1);
               nNeurons_store{1}(2); nNeurons_store{2}(2);
               nNeurons_store{1}(3); nNeurons_store{2}(3)];

% mouseLabels = {'Mouse 1','Mouse 3','Mouse 2'};
mouseLabels = {'Mouse 1','Mouse 2','Mouse 3'};


% Custom x positions with gaps between mice
xpos = [1 2 4 5 7 8];
barW = 0.7;

stackColors = [0.25 0.25 0.25;   % dark gray  - 1 target
               0.55 0.55 0.55;   % mid gray   - 2 targets
               0.82 0.82 0.82];  % light gray - 3+ targets

labels = {'1 target','2 targets','3+ targets'};

figure; set(gcf, 'Position', [200 200 550 380]);
hold on;

% Draw stacked bars as rectangles
for i = 1:6
    bottom = 0;
    for j = 1:3
        h = fill([xpos(i)-barW/2, xpos(i)+barW/2, xpos(i)+barW/2, xpos(i)-barW/2], ...
                 [bottom, bottom, bottom+stackData(i,j), bottom+stackData(i,j)], ...
                 stackColors(j,:), 'EdgeColor', 'none');
        % Only tag for legend on the first bar group
        if i == 1
            h.DisplayName = labels{j};
        else
            h.HandleVisibility = 'off';
        end
        bottom = bottom + stackData(i,j);
    end
    text(xpos(i), 1.03, ['n=' num2str(nNeuronsAll(i))], ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

set(gca, 'XTick', xpos, ...
    'XTickLabel', {'med','lat','med','lat','med','lat'}, ...
    'FontSize', 11, 'Box', 'off', 'TickDir', 'out');
xlim([0 9]);
ylim([0 1.12]);
ylabel('Fraction of neurons');

for m = 1:3
    text(mean(xpos((m-1)*2+1 : m*2)), -0.12, mouseLabels{m}, ...
        'HorizontalAlignment', 'center', 'FontSize', 11);
end

legend('Location', 'bestoutside', 'Box', 'off');
title('Projection target count per neuron');
hold off;

%% projection patterns and counts of unique barcodes - MED LAT


DATA = BigB;
% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<=cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
id_med=[3 7 11]; %medial SC
id_lat=[2 6 10]; %lateral SC
id_tot={id_med;id_lat};
names_medlat={'upper sSC','lower sSC'};

% id_tot={id_med;id_lat;14};
% names_medlat={'upper sSC','lower sSC', 'caudal sSC'};

names=sample_names(projTargets);

scnt = 0;
Numbers = zeros(length(sSC),2); allareas = []; all_targets=[];

for s = 1:size(id_tot,1) % go through injection sites
    scnt = scnt +1;
    id = id_tot{s};
    BCodes= find(ismember(maxes,id));
    Numbers(scnt,1) = length(BCodes);
    % for each barcode find the # of areas
    areas = [];
    for b = 1:length(BCodes)
        curr = DATA(BCodes(b),:);
        areas = [areas; curr];
    end

    %     figure

    allareas = [allareas; areas];


    %filter projection targets, and remove rows of all 0
    targets_bc=areas(:, projTargets);


    zero_rows=all(targets_bc==0, 2);
    targets_bc=targets_bc(~zero_rows, :);

    % new_areas=areas(:,id);
    % new_areas=new_areas(~zero_rows, :);
    % max_areas=max(new_areas,[],2);
    % mean_expression=mean(max_areas);


    %normalize values to the mean of max norm counts at injection site
    % targets_bc_norm=targets_bc./mean_expression;
    targets_bc_norm=targets_bc;
    all_targets=[all_targets ; targets_bc_norm];

    norm_vals=sum(targets_bc,2); %just to sort the matrix later
    targets_bc=targets_bc./norm_vals;

    % %sort the array
    % [~,loc]=max(targets_bc,[],2);
    % [~,ix]=sort(loc);
    % targets_bc=targets_bc(ix,:);


    %---SORT---

    % Preallocate a sorted matrix
    sortedMatrix = [];

    % Step 1: Handle rows with a value of 1 in any column
    for i = 1:length(projTargets)
        % Extract rows where column `col` has a value of 1
        rowsWithOnes = targets_bc_norm(targets_bc(:, i) == 1, :);
        rowsWithOnes=sort(rowsWithOnes,1, 'descend');

        % Append these rows to the sorted matrix
        sortedMatrix = [sortedMatrix; rowsWithOnes];
    end

    % Step 2: Handle rows with no values of 1
    remainingRows = targets_bc_norm(~any(targets_bc == 1, 2), :);

    % Sort the remaining rows by descending values in each column (priority order: col1, col2, ...)
    [~, sortIndices] = sortrows(remainingRows, -(1:length(projTargets)));
    sortedMatrix = [sortedMatrix; remainingRows(sortIndices, :)];


    disp(size(sortedMatrix,1))

    figure(1); set(gcf, 'Position',[85.8000 126.6000 962.4000 580.8000]);

    if s==1
        ax(1)=subplot(1,size(id_tot,1),scnt);
        % imagesc(targets_bc(ix,:));
        imagesc(sortedMatrix);
        % cols_1=cbrewer('seq', 'Oranges', 300);
        cols_1=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#FEF7EE','#F9DFC2','#F4B985', '#E79357', '#B86C44', '#764630'},300);
        cols_1=flipud(cols_1);
        cols_1(1,:)=[1 1 1];
        cols_1=cols_1(1:250,:);
        colormap(ax(1),cols_1)
    elseif s==2
        ax(2)=subplot(1,size(id_tot,1),scnt);
        imagesc(sortedMatrix);
        % cols_2=cbrewer('seq', 'Blues', 300);
        cols_2=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#F7FAFE','#D5DFEE','#9DBFD9', '#5A8EBF', '#2B599C', '#182D62'},300);
        cols_2=flipud(cols_2);
        cols_2(1,:)=[1 1 1];
        cols_2=cols_2(1:250,:);
        colormap(ax(2),cols_2)
    else
        ax(3)=subplot(1,size(id_tot,1),scnt);
        imagesc(sortedMatrix);
        cols_3=cbrewer('seq', 'Greens', 300);
        cols_3(1,:)=[1 1 1];
        colormap(ax(3),cols_3)
    end
    % cols(1,:)=[1 1 1]; %make 0 values value white
    % cols=cols(1:280,:); %make max values less dark
    % colormap(cols)


    set(gca,'xticklabel',names, 'xtick', 1:size(sortedMatrix,2),'ytickLabels',[],...
        'ytick', [], 'FontSize', 14)
    % xtickangle(90);
    % xlabel('Projection targets', 'FontSize',12)
    clim([0 30])
    ylabel('Barcodes', 'FontSize',14)
    cb=colorbar;
    ylabel(cb, 'Barcode counts');
    title([names_medlat{s}, ' projections'], 'FontSize', 14)


end


%% projection patterns of unique barcodes - sSCam sSCm sSCal sSCl sSCc

DATA = BigBnorm;
% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<=cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
id_tot={[3 7]; 11; [2 6]; 10; 14};
% id_tot={[3 7 11 2 6 10 14]};

names_medlat={'sSCam', 'sSCm', 'sSCal' ,'sSCl', 'sSCc'};
% names_medlat={'sSC'};

names=sample_names(projTargets);

scnt = 0;
Numbers = zeros(length(sSC),2); allareas = []; all_targets=[]; neuronnum=0;
% for s = orderOfSC % go through injection sites
for s = 1:size(id_tot,1) % go through injection sites
    scnt = scnt +1;
    id = id_tot{s};
    BCodes= find(ismember(maxes,id));
    Numbers(scnt,1) = length(BCodes);
    % for each barcode find the # of areas
    areas = [];
    for b = 1:length(BCodes)
        curr = DATA(BCodes(b),:);
        areas = [areas; curr];
    end

    %     figure

    allareas = [allareas; areas];


    %filter projection targets, and remove rows of all 0
    targets_bc=areas(:, projTargets);

    zero_rows=all(targets_bc==0, 2);
    targets_bc=targets_bc(~zero_rows, :);
    temp_num=targets_bc(:,1);
    temp_num(all(temp_num == 0, 2), :)=[];
    neuronnum=neuronnum+size(temp_num,1);


    %normalize values to sum of each row
    norm_vals=sum(targets_bc,2);
    targets_bc=targets_bc./norm_vals;
    all_targets=[all_targets ; targets_bc];

    %---SORT---

    % Preallocate a sorted matrix
    sortedMatrix = [];

    % Step 1: Handle rows with a value of 1 in any column
    for i = 1:length(projTargets)
        % Extract rows where column `col` has a value of 1
        rowsWithOnes = targets_bc(targets_bc(:, i) == 1, :);

        % Append these rows to the sorted matrix
        sortedMatrix = [sortedMatrix; rowsWithOnes];
    end

    % Step 2: Handle rows with no values of 1
    remainingRows = targets_bc(~any(targets_bc == 1, 2), :);

    % Sort the remaining rows by descending values in each column (priority order: col1, col2, ...)
    [~, sortIndices] = sortrows(remainingRows, -(1:length(projTargets)));
    sortedMatrix = [sortedMatrix; remainingRows(sortIndices, :)];

    % %group dSC regions into medial lateral anterior and posterior
    % if exist('dSCmi', 'var')
    %     sortedMatrix = [sortedMatrix(:,1:5) sum(sortedMatrix(:,dSCmi),2) ...
    %         sum(sortedMatrix(:,dSCli),2) sum(sortedMatrix(:,dSCai),2) ...
    %         sum(sortedMatrix(:,dSCpi),2)];
    %     names=[sample_names(projTargets(1:5)) sample_names(end-3:end)];
    % end

    % sorted_data = data(sort_idx, :);
    % [~,loc]=max(targets_bc,[],2);
    % [~,ix]=sort(loc);


    figure(1); set(gcf, 'Position',[422.6000 114.6000 1192 620.8000]);
    subplot(2,3,scnt);
    imagesc(sortedMatrix);
    cols=cbrewer('seq', 'Purples', 300);
    cols(1,:)=[1 1 1];
    % cols=cols(1:280,:); %make max values less dark
    colormap(cols)


    set(gca,'xticklabel',names, 'xtick', 1:size(sortedMatrix,2),'ytickLabels',[],...
        'ytick', [], 'FontSize', 8)
    xtickangle(90);
    % xlabel('Projection targets', 'FontSize',12)
    ylabel('Barcodes', 'FontSize',9)
    cb=colorbar;
    ylabel(cb, 'Projection strength');
    title([names_medlat{s}, ' projections', ' (',num2str(size(sortedMatrix,1)),')'], 'FontSize', 10)

    %fraction of neurons innervating a single area out of all neurons that
    %innervate that area
    figure(2); set(gcf, 'position', [321 169 816 473.6000])
    subplot(2,3,scnt)
    nonzero_counts=zeros(1, size(targets_bc, 2)); %number of neurons projecting to each area

    for i = 1:size(targets_bc, 2)
        nonzero_counts(i) = nnz(targets_bc(:, i));
    end
    ones_counts = zeros(1, size(targets_bc, 2)); %number of neurons that project only one particular area
    for i = 1:size(targets_bc, 2)
        ones_counts(i)=sum(targets_bc(:, i) == 1);
    end
    frac_ones=ones_counts./nonzero_counts;

    b=bar(frac_ones);
    b.FaceColor = 'flat';
    b.EdgeAlpha=0;
    b.CData=cols(200,:);


    % xlim([0.2 5.8]);
    xticks(1:length(projTargets));
    set(gca,'xticklabel',sample_names(projTargets), 'box', 'off', 'TickDir', 'out',...
        'FontSize', 8)
    xtickangle(90)
    ylim([0 1])
    ylabel('fraction of dedicated output', 'FontSize', 8)
    title(['dedicated outputs of ', names_medlat{s}, ' neurons'], 'FontSize', 8);

end

% %fraction of neurons innervating a single area out of all neurons that
% %innervate that area
% nonzero_counts=zeros(1, size(all_targets, 2)); %number of neurons projecting to each area
%
% for i = 1:size(all_targets, 2)
%     nonzero_counts(i) = nnz(all_targets(:, i));
% end
% ones_counts = zeros(1, size(all_targets, 2)); %number of neurons that project only one particular area
% for i = 1:size(all_targets, 2)
%     ones_counts(i)=sum(all_targets(:, i) == 1);
% end
% frac_ones=ones_counts./nonzero_counts;
%


%% with hierarchical clustering
DATA = BigBnorm;
% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
id_med=[3 7 11]; %medial SC
id_lat=[2 6 10]; %lateral SC
id_tot=[id_med,id_lat];
names_medlat={'medial sSC','lateral sSC', 'sSC'};


scnt = 1;
Numbers = zeros(length(sSC),2); allareas = []; all_targets=[];
% for s = orderOfSC % go through injection sites

id = id_tot;
BCodes= find(ismember(maxes,id));
Numbers(scnt,1) = length(BCodes);
% for each barcode find the # of areas
areas = [];
for b = 1:length(BCodes)
    curr = DATA(BCodes(b),:);
    areas = [areas; curr];
end

%     figure

allareas = [allareas; areas];


%filter projection targets, and remove rows of all 0
targets_bc=areas(:, projTargets);

zero_rows=all(targets_bc==0, 2);
targets_bc=targets_bc(~zero_rows, :);

%normalize values to sum of each row
norm_vals=sum(targets_bc,2);
targets_bc=targets_bc./norm_vals;
all_targets=[all_targets ; targets_bc];

%sort the array
[~,loc]=max(targets_bc,[],2);
[~,ix]=sort(loc);

figure(1); set(gcf, 'Position',[85.8000 282.6000 480.8000 379.4000]);
imagesc(targets_bc(ix,:));
cols=cbrewer('seq', 'Reds', 300);
cols(1,:)=[1 1 1]; %make 0 values value white
cols=cols(1:280,:); %make max values less dark

colormap(cols)
set(gca,'xticklabel',sample_names(projTargets), 'xtick', 1:5,'ytickLabels',[],...
    'ytick', [], 'FontSize', 9)
xtickangle(0);
xlabel('Projection targets', 'FontSize',10)
% ylabel('Area B', 'FontSize',10)
cb=colorbar;
ylabel(cb, 'Projection strength');
title(['Projection patterns of ', names_medlat{2}, ' neurons'], 'FontSize', 11)


% do hierarchical clustering on neurons
matrix=targets_bc(ix,:);
cols=cbrewer('seq', 'Reds', 300);
cols(1,:)=[1 1 1]; %make 0 values value white
cols=cols(1:280,:); %make max values less dark
clust=clustergram(matrix,'OptimalLeafOrder','false','RowPDist','cosine','Standardize',3,'Cluster',1,'Colormap',cols,...
    'Symmetric', 'false','LogTrans','false');
set(clust,'ColumnLabels',sample_names(projTargets), 'ColumnLabelsRotate', 0)
% addTitle(clust,['Projection patterns of ', names_medlat{1}, ' neurons'], 'FontSize', 11)


%% do k-means clustering
id=id_med;
if length(id)==length(id_med) && isequal(id,id_med)
    reg=1;
elseif length(id)==length(id_lat) && isequal(id,id_lat)
    reg=2;
else
    reg=3;
end
matrix=targets_bc(ix,:);

%remove single projecting neurons
nonzeros_per_row = sum(matrix ~= 0, 2);
matrix = matrix(nonzeros_per_row >= 2, :);
% % check cluster stability
% clusterstability_kmeans=evalclusters(matrix,'kmeans','gap','Klist',[1:20],...
%     'Distance','cosine');
% figure;plot(clusterstability_kmeans.CriterionValues);
% xlabel('number of clusters')
% ylabel('GAP criterion')
% box off
% xlim([1 20])
%
% clusterstability_kmeans2=evalclusters(matrix,'kmeans','silhouette','Klist',[1:20],...
%     'Distance','cosine');
% figure;plot(clusterstability_kmeans2.CriterionValues);
% xlabel('number of clusters')
% ylabel('silhouettee criterion')
% box off

% %cluster MAPseq data by starting from high k and decreaseing based on
% % threshold
% k=10; %chosen cluster number
%
% % Set a threshold for the minimum number of points in a cluster
% minClusterSize = 5;
%
% while true
%     % Perform K-means clustering
%     [T,C,sumd,D]=kmeans(matrix,k,'Distance','cosine', 'Replicates', 5);
%
%     % Count the number of points in each cluster
%     clusterSizes = histcounts(T, k);
%
%     % Check if all clusters meet the size threshold
%     if all(clusterSizes >= minClusterSize)
%         break;  % Exit loop if all clusters meet the size threshold
%     else
%         % If any cluster is below the threshold, adjust the number of clusters
%         % You can increase the number of clusters, or handle it as needed
%         k = k - 1;  % Decrease the number of clusters
%     end
% end

%cluster MAPseq data
k=10; %chosen cluster number

% Set a threshold for the minimum number of points in a cluster
minClusterSize = 0;
% Perform K-means clustering
[T,C,sumd,D]=kmeans(matrix,k,'Distance','cosine', 'MaxIter',1000);
% Count the number of points in each cluster
clusterSizes = histcounts(T, k);
% Identify clusters that meet the size threshold
clustersToKeep = find(clusterSizes >= minClusterSize);
% Create a new index array to store the new cluster assignments
newIdx = zeros(size(T));  % Initialize new index array
for i = 1:numel(clustersToKeep)
    newIdx(T == clustersToKeep(i)) = i;  % Reassign indices for valid clusters
end
% Filter centroids to keep only those that meet the threshold
newCentroids = C(clustersToKeep, :);  % Keep only centroids of valid clusters

% Set new indices for points that belong to removed clusters to 0 or another label
newIdx(isnan(newIdx)) = 0;  % Mark points from removed clusters with 0 (or any other label you prefer)
C=newCentroids;
T=newIdx;

%sort centroid matrix based on max value of column
sorted_C = zeros(size(C));
sorted_indices = zeros(size(C, 1), 1);  % This will store the row indices

selected_rows = false(size(C, 1), 1);

% First pass: sort based on the max value of each column for the first 5 rows
for i = 1:size(C, 2)  % Loop over each column (5 in this case)
    [~, idx] = max(C(:, i));  % Find the row index of the max value in column i
    sorted_C(i, :) = C(idx, :);  % Move the row with max value in column i to row i
    sorted_indices(i) = idx;  % Track the original row index
    selected_rows(idx) = true;  % Mark this row as selected
end

% Second pass: sort the remaining rows
remaining_rows = find(~selected_rows);  % Find the indices of the unselected rows
remaining_C = C(remaining_rows, :);  % Extract the remaining rows

%sort the remaining rows normally
if ~isempty(remaining_C)
    [sorted_C(sum(selected_rows)+1:end,:), idx]=sortrows(remaining_C, 'descend');
    sorted_indices(sum(selected_rows)+1:end)=remaining_rows(idx);
end
% % Loop again to find max values for the remaining rows
% for i = 1:size(remaining_C, 1)  % Loop over each column
%     [~, idx] = max(remaining_C(:, i));  % Find the row index of the max value
%     sorted_C(sum(selected_rows) + i, :) = remaining_C(idx, :);  % Place in rows of sorted C
%     sorted_indices(sum(selected_rows) + i) = remaining_rows(idx);  % Track the original row index
%     remaining_C(idx, :) = [];  % Remove this row from remaining_C
%     remaining_rows(idx) = [];  % Remove this row's index from remaining_rows
% end

%------plot cluster centroids matrix as bubbles
% Get the size of the matrix
[m, n] = size(sorted_C);
% Generate the x and y coordinates for each bubble
[X, Y] = meshgrid(1:n, 1:m);
% Flatten the coordinates and the data for scatter plotting
X = X(:);
Y = Y(:);
bubbleColors = sorted_C(:);  % Flatten the data values for coloring
% Plot the bubbles (same size, color representing magnitude)
figure(1); set(gcf, 'Position',[268.2000 205.8000 339.2000 456.2000]);
scatter(X, Y, 400, bubbleColors, 'filled');  % '500' sets the constant bubble size
% imagesc(sorted_C);
cols=cbrewer('seq', 'Blues', 300);
cols(1,:)=[1 1 1]; %make 0 values value white
cols=cols(20:300,:); %make min values less light
colormap(cols)
set(gca,'xticklabel',sample_names(projTargets), 'xtick', 1:5,...
    'ytick', 1:size(sorted_C,1), 'FontSize', 10, 'ticklength', [0,0],...
    'YDir', 'reverse')
ax=gca;
ax.XAxis.Color=[1 1 1];
ax.XAxis.TickLabelColor=[0 0 0];
ax.XAxis.Label.Color=[0 0 0];
ax.YAxis.Color=[1 1 1];
ax.YAxis.TickLabelColor=[0 0 0];
ax.YAxis.Label.Color=[0 0 0];
xtickangle(0);
% xlabel('Projection targets', 'FontSize',10)
ylabel('k-means cluster', 'FontSize',11)
xlim([0.5 5.5])
ylim([0.5 size(sorted_C,1)+0.5])
cb=colorbar;
ylabel(cb, 'Projection strength');
cb.Ticks = [0 0.99];
cb.TickLabels = {"0","1"};
cb.FontSize=10;
% title(['Projection patterns of ', names_medlat{2}, ' neurons'], 'FontSize', 11)

%sort matrix based on clusters
[~, index_positions] = ismember(1:size(C,1), sorted_indices);
matrix(T==0,:)=[];
T(T==0)=[];
new_T = index_positions(T);% Replace values in B using the index positions
[sorted_T, sortedIndices] = sort(new_T);
sorted_matrix=matrix(sortedIndices,:);

figure(2); set(gcf, 'Position',[654.6000 281.8000 480.8000 379.2000]);
imagesc(sorted_matrix);
cols=cbrewer('seq', 'Blues', 300);
cols(1,:)=[1 1 1]; %make 0 values value white
% cols=cols(1:280,:); %make max values less dark
colormap(cols)
set(gca,'xticklabel',sample_names(projTargets), 'xtick', 1:5,'ytickLabels',[],...
    'ytick', [], 'FontSize', 9)
xtickangle(0);
xlabel('Projection targets', 'FontSize',10)
% ylabel('Area B', 'FontSize',10)
cb=colorbar;
ylabel(cb, 'Projection strength');
title(['Projection patterns of ', names_medlat{reg}, ' neurons'], 'FontSize', 11)

%t-SNE
figure(3); set(gcf, 'position', [488 258.6000 557 403.4000])
tsne_clust=tsne(sorted_matrix, 'Distance','cosine');
% tsne_clust=tsne(matrix);
cmap = colormap(turbo); %divide colormap into k values
indices = round(linspace(1, size(cmap, 1), size(sorted_C,1)));  % Equally spaced indices
group_colors = cmap(indices, :);  % Get the 10 selected colors
single_colors = group_colors(sorted_T, :);
% sizes=repmat(15,size(sorted_C,1),2);
% gscatter(tsne_clust(:,1),tsne_clust(:,2),sorted_T, group_colors,'.', sizes);
scatter(tsne_clust(:,1),tsne_clust(:,2),10,single_colors, 'filled',...
    'MarkerFaceAlpha',0.7);
set(gca, 'FontSize', 11, 'color', 'none')
xlabel('t-SNE 1', 'FontSize', 12);
ylabel('t-SNE 2', 'FontSize', 12);
title(['t-SNE projection of ', names_medlat{reg}, ' k-means clusters'])
axis equal
hold on;
hLegend = gobjects(size(sorted_C,1), 1);
%Create hLegend(i) scatter plots for the legend (one for each of the 10 colors)
for i = 1:size(sorted_C,1)
    hLegend(i)=scatter(nan, nan, 50, group_colors(i, :), 'filled');  % Invisible scatter points for legend
    hold on;
end

%Add the legend
legend(hLegend,arrayfun(@(x) sprintf('k %d', x), 1:size(sorted_C,1),...
    'UniformOutput', false), 'Location','northeastoutside');
legend boxoff

hold off;

%% plot projection strength for each SCs region and mice
DATA1 = B1norm;
DATA2 = B2norm;
DATA3 = B3norm;
threshold = 0;
% cmp = colormap('prism');
% cmp = [0 0 0;cmp];

scnt = 0; Numbers={};
Numbers{1} = zeros(length(sSC),2);Numbers{2} = zeros(length(sSC),2);
Numbers{3} = zeros(length(sSC),2); allKeep1 =[]; allKeep2=[]; allKeep3=[];
all_targets1=[]; all_targets2=[]; all_targets3=[];
for s = orderOfSC % go through injection sites
    scnt = scnt +1;
    id = sSC(s);
    BCodes1 = find(m1==id);
    BCodes2 = find(m2==id);
    BCodes3 = find(m3==id);
    Numbers{1}(scnt,1) = length(BCodes1);
    Numbers{2}(scnt,1) = length(BCodes2);
    Numbers{3}(scnt,1) = length(BCodes3);
    % for each barcode find the # of areas
    areas1 = []; areas2 = []; areas3 = [];
    for b = 1:length(BCodes1)
        curr = find(DATA1(BCodes1(b),:)>threshold);
        areas1 = [areas1;length(curr)];
    end
    for b = 1:length(BCodes2)
        curr = find(DATA2(BCodes2(b),:)>threshold);
        areas2 = [areas2;length(curr)];
    end
    for b = 1:length(BCodes3)
        curr = find(DATA3(BCodes3(b),:)>threshold);
        areas3 = [areas3;length(curr)];
    end
    %     figure
    [~,a2] = sort(areas1,'descend');
    data1 = DATA1(BCodes1(a2),:);
    [~,a2] = sort(areas2,'descend');
    data2 = DATA2(BCodes2(a2),:);
    [~,a2] = sort(areas3,'descend');
    data3 = DATA3(BCodes3(a2),:);
    keepData1 = data1; keepData2 = data2; keepData3 = data3;
    allKeep1 = [allKeep1; keepData1];
    allKeep2 = [allKeep2; keepData2];
    allKeep3 = [allKeep3; keepData3];

    %barcodes in each target area
    targets_bc1=data1(:, projTargets);
    zero_rows=all(targets_bc1==0, 2);
    targets_bc1=targets_bc1(~zero_rows, :);

    targets_bc2=data2(:, projTargets);
    zero_rows=all(targets_bc2==0, 2);
    targets_bc2=targets_bc2(~zero_rows, :);

    targets_bc3=data3(:, projTargets);
    zero_rows=all(targets_bc3==0, 2);
    targets_bc3=targets_bc3(~zero_rows, :);

    %normalize values to sum of each row
    norm_vals=sum(targets_bc1,2);
    targets_bc1=targets_bc1./norm_vals;
    all_targets1=[all_targets1 ; targets_bc1];

    norm_vals=sum(targets_bc2,2);
    targets_bc2=targets_bc2./norm_vals;
    all_targets2=[all_targets2 ; targets_bc2];

    norm_vals=sum(targets_bc3,2);
    targets_bc3=targets_bc3./norm_vals;
    all_targets3=[all_targets3 ; targets_bc3];

    % tmp = find(keepData<threshold);
    % keepData(tmp) = 0;

    %projection strength

    targetsMatrix2_m1 = zeros(length(projTargets),1);
    targetsMatrix2_m2 = zeros(length(projTargets),1);
    targetsMatrix2_m3 = zeros(length(projTargets),1);
    for pr = 1:length(projTargets)
        targetsMatrix2_m1(pr,1) = sum(keepData1(:,projTargets(pr)));
        targetsMatrix2_m2(pr,1) = sum(keepData2(:,projTargets(pr)));
        targetsMatrix2_m3(pr,1) = sum(keepData3(:,projTargets(pr)));
    end

    %normalize
    targetsMatrix2_m1=targetsMatrix2_m1./sum(targetsMatrix2_m1);
    targetsMatrix2_m2=targetsMatrix2_m2./sum(targetsMatrix2_m2);
    targetsMatrix2_m3=targetsMatrix2_m3./sum(targetsMatrix2_m3);

    targetsMatrix1_m1 = zeros(length(projTargets),length(projTargets));
    targetsMatrix1_m2 = zeros(length(projTargets),length(projTargets));
    targetsMatrix1_m3 = zeros(length(projTargets),length(projTargets));
    for pr = 1:length(projTargets)
        tmp1 = find(keepData1(:,projTargets(pr)));
        tmp2 = find(keepData2(:,projTargets(pr)));
        tmp3 = find(keepData3(:,projTargets(pr)));
        val1_m1 = keepData1(tmp1,projTargets(pr));
        val1_m2 = keepData2(tmp2,projTargets(pr));
        val1_m3 = keepData3(tmp3,projTargets(pr));
        for pr2 = 1:length(projTargets)
            val2_m1 = keepData1(tmp1,projTargets(pr2));
            val2_m2 = keepData2(tmp2,projTargets(pr2));
            val2_m3 = keepData3(tmp3,projTargets(pr2));
            diff1 = val2_m1-val1_m1; diff1=sum(diff1(diff1<0));
            diff2 = val2_m2-val1_m2; diff2=sum(diff2(diff2<0));
            diff3 = val2_m3-val1_m3; diff3=sum(diff3(diff3<0));
            proj_fraction1=(sum(val1_m1)+diff1)/sum(val1_m1);
            proj_fraction2=(sum(val1_m2)+diff2)/sum(val1_m2);
            proj_fraction3=(sum(val1_m3)+diff3)/sum(val1_m3);
            targetsMatrix1_m1(pr2,pr) = proj_fraction1;
            targetsMatrix1_m2(pr2,pr) = proj_fraction2;
            targetsMatrix1_m3(pr2,pr) = proj_fraction3;
        end
    end


    %calculate mean
    targetsMatrix2_mean=mean(cat(3, targetsMatrix2_m1,targetsMatrix2_m2,targetsMatrix2_m3),3,'omitnan');
    targetsMatrix1_mean=mean(cat(3, targetsMatrix1_m1,targetsMatrix1_m2,targetsMatrix1_m3),3,'omitnan');
    %set diagonal values to NaN
    targetsMatrix1_mean(1:5:end)= NaN;


    figure
    %set nan values to transparent
    imAlpha=ones(size(targetsMatrix1_mean));
    imAlpha(isnan(targetsMatrix1_mean))=0;
    imagesc(targetsMatrix1_mean,'AlphaData',imAlpha);
    set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
    'color', [0.5 0.5 0.5])
    xlabel('Area A', 'FontSize',10)
    ylabel('Area B', 'FontSize',10)
    col1=magma(300);
    colormap(col1); colorbar;
    title([SampleName{3}(sSC(s)) ': subset of area A projections projecting to area B'])
%     title([SampleName{3}(sSC(s)) ': % barcodes found in >1 long-range targets'])
% %
    Numbers{1}(scnt,2) = sum(targetsMatrix2_m1);
    Numbers{2}(scnt,2) = sum(targetsMatrix2_m2);
    Numbers{3}(scnt,2) = sum(targetsMatrix2_m3);
    targetsMatrix2_all=[targetsMatrix2_m1 targetsMatrix2_m2 targetsMatrix2_m3];

    figure(100)
    subplot(3,3,scnt)
    % bar(targetsMatrix2_mean); hold on
    % scatter(1:10, targetsMatrix2_all, 10, 'filled')
    plot(1:5,targetsMatrix2_m1,'-o','Color', '#6A2C70', LineWidth=2); hold on
    plot(1:5,targetsMatrix2_m2,'-o','Color', '#B83B5E', LineWidth=2); hold on
    plot(1:5,targetsMatrix2_m3,'-o','Color', '#FFAF45', LineWidth=2)
    xlim([0.5 5.5]); xticks(1:5);
    set(gca,'xticklabel',sample_names(projTargets))
    xtickangle(90)
    ylabel('projection strength')
    title(SampleName{3}{sSC(s)})
    legend('m1', 'm2', 'm3');

%     ylabel('# of barcodes')
%     title([SampleName{3}{sSC(s)},' | total: ',int2str(Numbers(scnt,1)),', in targets: ',int2str(Numbers(scnt,2)) ])
end

%% plot projection strength for each SCs region and mice - MED LAT
DATA1 = B1norm;
DATA2 = B2norm;
DATA3 = B3norm;
DATA = BigBnorm;
threshold = 0;

data = BigB;
cutoff = max(data(:,1));
cutoff_tmp = BigB(:,:) <= cutoff;
DATA(cutoff_tmp) = 0;

cutoff_tmp = B1(:,:) <= cutoff;
DATA1(cutoff_tmp) = 0;

cutoff_tmp = B2(:,:) <= cutoff;
DATA2(cutoff_tmp) = 0;

cutoff_tmp = B3(:,:) <= cutoff;
DATA3(cutoff_tmp) = 0;

id_med = [3 7 11];
id_lat = [2 6 10];
id_ant = [2 3 6 7];
id_post = [10 11 14];
id_tot = {id_med; id_lat; id_ant; id_post};
names_medlat = {'upper sSC', 'lower sSC', 'anterior sSC', 'posterior sSC'};

scnt = 0; Numbers = {};
Numbers{1} = zeros(length(sSC),2); Numbers{2} = zeros(length(sSC),2);
Numbers{3} = zeros(length(sSC),2);
all_targets1 = []; all_targets2 = []; all_targets3 = []; all_targets = [];
targets_med = []; targets_lat = []; targets_ant = []; targets_post = [];
targets_med_tot = []; targets_lat_tot = []; targets_ant_tot = []; targets_post_tot = [];
targets_med1 = []; targets_med2 = []; targets_med3 = [];
targets_lat1 = []; targets_lat2 = []; targets_lat3 = [];
targets_med1_tot = []; targets_med2_tot = []; targets_med3_tot = [];
targets_lat1_tot = []; targets_lat2_tot = []; targets_lat3_tot = [];
targetsMatrix2_ind = {}; targetsMatrix2_all = []; targetsMatrix2_all_nonorm = [];
proj_prob = []; targetMatrices2 = {};

for s = 1:size(id_tot,1)
    scnt = scnt + 1;
    id = id_tot{s};
    BCodes1 = find(ismember(m1, id));
    BCodes2 = find(ismember(m2, id));
    BCodes3 = find(ismember(m3, id));
    BCodes = find(ismember(maxes, id));
    Numbers{1}(scnt,1) = length(BCodes1);
    Numbers{2}(scnt,1) = length(BCodes2);
    Numbers{3}(scnt,1) = length(BCodes3);

    areas1 = []; areas2 = []; areas3 = []; areas = [];
    for b = 1:length(BCodes1)
        curr = find(DATA1(BCodes1(b),:) > threshold);
        areas1 = [areas1; length(curr)];
    end
    for b = 1:length(BCodes2)
        curr = find(DATA2(BCodes2(b),:) > threshold);
        areas2 = [areas2; length(curr)];
    end
    for b = 1:length(BCodes3)
        curr = find(DATA3(BCodes3(b),:) > threshold);
        areas3 = [areas3; length(curr)];
    end
    for b = 1:length(BCodes)
        curr = find(DATA(BCodes(b),:) > threshold);
        areas = [areas; length(curr)];
    end

    % Probability of projection
    new_data = data(BCodes,:);
    data_binary_tmp = double(new_data(:,projTargets) > cutoff);
    data_binary = data_binary_tmp(sum(data_binary_tmp,2) ~= 0,:);
    sum_projections = sum(data_binary);
    neuronnum = length(BCodes);
    proj_prob(s,:) = sum_projections ./ neuronnum;

    [~,a2] = sort(areas1, 'descend');
    data1 = DATA1(BCodes1(a2),:);
    [~,a2] = sort(areas2, 'descend');
    data2 = DATA2(BCodes2(a2),:);
    [~,a2] = sort(areas3, 'descend');
    data3 = DATA3(BCodes3(a2),:);
    data_all = DATA(BCodes,:);
    data_all(data_all(:,:) < threshold) = 0;
    keepData1 = data1; keepData2 = data2; keepData3 = data3;

    % Barcodes in each target area
    targets_bc1 = data1(:, projTargets);
    targets_bc1(targets_bc1 < threshold) = 0;
    targets_bc1_all = targets_bc1;
    targets_bc1 = targets_bc1(~all(targets_bc1 == 0, 2), :);

    targets_bc2 = data2(:, projTargets);
    targets_bc2(targets_bc2 < threshold) = 0;
    targets_bc2_all = targets_bc2;
    targets_bc2 = targets_bc2(~all(targets_bc2 == 0, 2), :);

    targets_bc3 = data3(:, projTargets);
    targets_bc3(targets_bc3 < threshold) = 0;
    targets_bc3_all = targets_bc3;
    targets_bc3 = targets_bc3(~all(targets_bc3 == 0, 2), :);

    targets_bc = data_all(:, projTargets);
    targets_bc(targets_bc < threshold) = 0;
    targets_bc_all = targets_bc;
    targets_bc = targets_bc(~all(targets_bc == 0, 2), :);

    % Normalize to sum of each row
    all_targets1 = [all_targets1; targets_bc1 ./ sum(targets_bc1, 2)];
    all_targets2 = [all_targets2; targets_bc2 ./ sum(targets_bc2, 2)];
    all_targets3 = [all_targets3; targets_bc3 ./ sum(targets_bc3, 2)];

    norm_vals = sum(targets_bc, 2);

    if s == 1
        targets_med = targets_bc; targets_med_tot = targets_bc_all;
        targets_med1 = targets_bc1; targets_med2 = targets_bc2; targets_med3 = targets_bc3;
        targets_med1_tot = targets_bc1_all; targets_med2_tot = targets_bc2_all; targets_med3_tot = targets_bc3_all;
        n_per_mouse = [size(targets_bc1,1), size(targets_bc2,1), size(targets_bc3,1)];
    elseif s == 2
        targets_lat = targets_bc; targets_lat_tot = targets_bc_all;
        targets_lat1 = targets_bc1; targets_lat2 = targets_bc2; targets_lat3 = targets_bc3;
        targets_lat1_tot = targets_bc1_all; targets_lat2_tot = targets_bc2_all; targets_lat3_tot = targets_bc3_all;
        n_per_mouse = [size(targets_bc1,1), size(targets_bc2,1), size(targets_bc3,1)];
    elseif s == 3
        targets_ant = targets_bc; targets_ant_tot = targets_bc_all;
        n_per_mouse = [size(targets_bc1,1), size(targets_bc2,1), size(targets_bc3,1)];
    else
        targets_post = targets_bc; targets_post_tot = targets_bc_all;
        n_per_mouse = [size(targets_bc1,1), size(targets_bc2,1), size(targets_bc3,1)];
    end
    targets_bc = targets_bc ./ norm_vals;
    all_targets = [all_targets; targets_bc];

    % Projection strength per mouse
    targetsMatrix2_m1 = zeros(length(projTargets), 1);
    targetsMatrix2_m2 = zeros(length(projTargets), 1);
    targetsMatrix2_m3 = zeros(length(projTargets), 1);
    targetsMatrix2_comb = zeros(length(projTargets), 1);
    for pr = 1:length(projTargets)
        targetsMatrix2_m1(pr,1) = sum(keepData1(:, projTargets(pr)));
        targetsMatrix2_m2(pr,1) = sum(keepData2(:, projTargets(pr)));
        targetsMatrix2_m3(pr,1) = sum(keepData3(:, projTargets(pr)));
        targetsMatrix2_comb(pr,1) = sum(data_all(:, projTargets(pr)));
    end

    targetMatrices2{s}(:,1) = targetsMatrix2_m1;
    targetMatrices2{s}(:,2) = targetsMatrix2_m2;
    targetMatrices2{s}(:,3) = targetsMatrix2_m3;

    % Normalize
    targetsMatrix2_m1 = targetsMatrix2_m1 ./ sum(targetsMatrix2_m1);
    targetsMatrix2_m2 = targetsMatrix2_m2 ./ sum(targetsMatrix2_m2);
    targetsMatrix2_m3 = targetsMatrix2_m3 ./ sum(targetsMatrix2_m3);

    targetsMatrix2_nonorm = targetsMatrix2_comb;
    targetsMatrix2_comb = targetsMatrix2_comb ./ sum(targetsMatrix2_comb);

    targetsMatrix2_all = [targetsMatrix2_all; targetsMatrix2_comb'];
    targetsMatrix2_all_nonorm = [targetsMatrix2_all_nonorm; targetsMatrix2_nonorm'];

    % Conditional probability matrix
    targetsMatrix1_all = zeros(length(projTargets), length(projTargets));
    for pr = 1:length(projTargets)
        tmp_all = find(data_all(:, projTargets(pr)));
        val1_all = data_all(tmp_all, projTargets(pr));
        for pr2 = 1:length(projTargets)
            val2_all = data_all(tmp_all, projTargets(pr2));
            diff_all = val2_all - val1_all; diff_all = sum(diff_all(diff_all < 0));
            proj_fraction_all = (sum(val1_all) + diff_all) / sum(val1_all);
            targetsMatrix1_all(pr2, pr) = proj_fraction_all;
        end
    end

    targetsMatrix1_all(1:length(projTargets)+1:end) = NaN;

    figure
    imAlpha = ones(size(targetsMatrix1_all));
    imAlpha(isnan(targetsMatrix1_all)) = 0;
    imagesc(targetsMatrix1_all, 'AlphaData', imAlpha);
    set(gca, 'xticklabel', sample_names(projTargets), 'yticklabel', sample_names(projTargets), ...
        'color', [0.5 0.5 0.5], 'FontSize', 14)
    xlabel('Area A', 'FontSize', 14)
    ylabel('Area B', 'FontSize', 14)
    xticks(1:length(projTargets)); yticks(1:length(projTargets));

    if s == 1
        col1 = customcolormap([0 0.2 0.4 0.6 0.8 1], {'#FEF7EE','#F9DFC2','#F4B985','#E79357','#B86C44','#764630'}, 300);
        col1 = flipud(col1);
    else
        col1 = customcolormap([0 0.2 0.4 0.6 0.8 1], {'#F7FAFE','#D5DFEE','#9DBFD9','#5A8EBF','#2B599C','#182D62'}, 300);
        col1 = flipud(col1);
    end
    col1(1,:) = [1 1 1];
    col1 = col1(1:250,:);
    colormap(col1); colorbar; clim([0 0.6])
    title([names_medlat{s} ': subset of area A outputs projecting to area B'], 'FontSize', 13)

    Numbers{1}(scnt,2) = sum(targetsMatrix2_m1);
    Numbers{2}(scnt,2) = sum(targetsMatrix2_m2);
    Numbers{3}(scnt,2) = sum(targetsMatrix2_m3);
    targetsMatrix2_ind{s} = [targetsMatrix2_m1 targetsMatrix2_m2 targetsMatrix2_m3];

    % Per-mouse projection strength with neuron counts in legend
    figure(100); set(gcf, 'Position', [97.8000 337.8000 950.2000 324.2000]);
    subplot(1, size(id_tot,1), scnt)
    plot(1:length(projTargets), targetsMatrix2_m1, '-o', 'Color', '#6A2C70', 'LineWidth', 2); hold on
    plot(1:length(projTargets), targetsMatrix2_m2, '-o', 'Color', '#B83B5E', 'LineWidth', 2);
    plot(1:length(projTargets), targetsMatrix2_m3, '-o', 'Color', '#FFAF45', 'LineWidth', 2)
    xlim([0.5 length(projTargets)+0.5]); xticks(1:length(projTargets)); ylim([0 0.9])
    set(gca, 'xticklabel', sample_names(projTargets), 'TickDir', 'out', 'box', 'off', 'FontSize', 11)
    xtickangle(90)
    ylabel('projection strength')
    title(names_medlat{s})
    legend(sprintf('m1 (n=%d)', n_per_mouse(1)), ...
           sprintf('m2 (n=%d)', n_per_mouse(2)), ...
           sprintf('m3 (n=%d)', n_per_mouse(3)), ...
           'Location', 'best');
    legend('boxoff');
    legend('boxoff');
end

% Paired scatter plot: upper vs lower per mouse, color-coded by mouse
figure(5); set(gcf, 'Position', [473 285 669.6000 296])
x_pos = [(1:length(projTargets))-0.2; (1:length(projTargets))+0.2];

% Define per-mouse shades of upper (orange) and lower (blue)
upper_colors = [0.95 0.65 0.3;   % mouse 1 - light orange
                0.85 0.45 0.15;  % mouse 2 - medium orange
                0.65 0.3 0.05];  % mouse 3 - dark orange
lower_colors = [0.5 0.7 0.9;    % mouse 1 - light blue
                0.25 0.5 0.75;   % mouse 2 - medium blue
                0.1 0.3 0.55];   % mouse 3 - dark blue

hold on;

% Get per-mouse neuron counts for upper and lower
n_upper = [size(targets_med1,1), size(targets_med2,1), size(targets_med3,1)];
n_lower = [size(targets_lat1,1), size(targets_lat2,1), size(targets_lat3,1)];

% Plot lines first (behind dots)
for m = 1:3
    for i = 1:length(projTargets)
        l = plot([x_pos(1,i) x_pos(2,i)], [targetsMatrix2_ind{1}(i,m) targetsMatrix2_ind{2}(i,m)], ...
            'LineWidth', 1.2, 'color', [0.75 0.75 0.75]);
        uistack(l, 'bottom');
    end
end

% Plot scatter dots per mouse
h_upper = gobjects(3,1);
h_lower = gobjects(3,1);
for m = 1:3
    h_upper(m) = scatter(x_pos(1,:), targetsMatrix2_ind{1}(:,m)', 40, 'filled', ...
        'MarkerFaceColor', upper_colors(m,:), 'MarkerEdgeColor', 'none');
    h_lower(m) = scatter(x_pos(2,:), targetsMatrix2_ind{2}(:,m)', 40, 'filled', ...
        'MarkerFaceColor', lower_colors(m,:), 'MarkerEdgeColor', 'none');
end

xlim([0.5 length(projTargets)+0.5]); xticks(1:length(projTargets)); ylim([0 0.9])
set(gca, 'xticklabel', sample_names(projTargets), 'FontSize', 14, 'TickDir', 'out', 'box', 'off')
ylabel('projection strength')

legend([h_upper; h_lower], ...
    {sprintf('upper m1 (n=%d)', n_upper(1)), ...
     sprintf('upper m2 (n=%d)', n_upper(2)), ...
     sprintf('upper m3 (n=%d)', n_upper(3)), ...
     sprintf('lower m1 (n=%d)', n_lower(1)), ...
     sprintf('lower m2 (n=%d)', n_lower(2)), ...
     sprintf('lower m3 (n=%d)', n_lower(3))}, ...
    'Location', 'best', 'FontSize', 9);
legend('boxoff');


% Per-mouse log2 fold change breakdown
nTargets = length(projTargets);
target_names = sample_names(projTargets);

obs_med = mean(targets_med, 1);
obs_lat = mean(targets_lat, 1);
obs_log2fc = log2(obs_med ./ obs_lat);

mouse_log2fc = zeros(3, nTargets);
mouse_n_med = [size(targets_med1,1), size(targets_med2,1), size(targets_med3,1)];
mouse_n_lat = [size(targets_lat1,1), size(targets_lat2,1), size(targets_lat3,1)];

med_means = {mean(targets_med1,1), mean(targets_med2,1), mean(targets_med3,1)};
lat_means = {mean(targets_lat1,1), mean(targets_lat2,1), mean(targets_lat3,1)};

fprintf('\n=== PER-MOUSE LOG2 FOLD CHANGE ===\n');
fprintf('%-12s %12s %12s %12s %12s\n', 'Target', 'Mouse1', 'Mouse2', 'Mouse3', 'Pooled');
fprintf('%s\n', repmat('-', 1, 60));
for i = 1:nTargets
    for m = 1:3
        if med_means{m}(i) == 0 && lat_means{m}(i) == 0
            mouse_log2fc(m,i) = NaN;
        elseif lat_means{m}(i) == 0
            mouse_log2fc(m,i) = Inf;
        elseif med_means{m}(i) == 0
            mouse_log2fc(m,i) = -Inf;
        else
            mouse_log2fc(m,i) = log2(med_means{m}(i) / lat_means{m}(i));
        end
    end
    fprintf('%-12s %11.3f %11.3f %11.3f %11.3f\n', ...
        target_names{i}, mouse_log2fc(1,i), mouse_log2fc(2,i), mouse_log2fc(3,i), obs_log2fc(i));
end

fprintf('\nNeuron counts: Med = [%d, %d, %d], Lat = [%d, %d, %d]\n', ...
    mouse_n_med, mouse_n_lat);

fprintf('\nDirection consistency (all mice same sign as pooled?):\n');
for i = 1:nTargets
    valid = ~isnan(mouse_log2fc(:,i));
    if sum(valid) == 0
        fprintf('  %s: no valid data\n', target_names{i});
    else
        same_dir = all(sign(mouse_log2fc(valid,i)) == sign(obs_log2fc(i)));
        fprintf('  %s: %s (%d/%d mice valid)\n', target_names{i}, string(same_dir), sum(valid), 3);
    end
end

% %plot ipsi vs contra for medial and lateral
% ipsi_contra_mat=[sum(targetsMatrix2_all(1,1:length(projTargets)/2)), sum(targetsMatrix2_all(2,1:length(projTargets)/2));...
%     sum(targetsMatrix2_all(1,length(projTargets)/2+1:length(projTargets))), sum(targetsMatrix2_all(2,length(projTargets)/2+1:length(projTargets)))];
%
% figure(11)
% coloritos=[col(3,:);col(4,:)];
% h=bar(1:2,ipsi_contra_mat, 'FaceAlpha', 1, 'EdgeColor', 'none');
% for k = 1:length(h)
%     h(k).FaceColor = 'flat';
%     h(k).CData = repmat(coloritos(k, :), size(ipsi_contra_mat, 2), 1);
% end
% set(gca,'xticklabel',{"ipsi", "contra"},...
%      'TickDir', 'out', 'box', 'off', 'FontSize',14)
% ylabel('Projection strength', 'FontSize', 14)
% xticks(1:2); ylim([0 1])
% legend('upper', 'lower')

% Among projecting neurons, what fraction projects ipsi vs contra?
n_ipsi_targets = floor(length(projTargets)/2);

% targets_med / targets_lat are projecting neurons only (zero rows removed)
n_ipsi_upper = sum(any(targets_med(:, 1:n_ipsi_targets) > 0, 2));
n_contra_upper = sum(any(targets_med(:, n_ipsi_targets+1:end) > 0, 2));
n_ipsi_lower = sum(any(targets_lat(:, 1:n_ipsi_targets) > 0, 2));
n_contra_lower = sum(any(targets_lat(:, n_ipsi_targets+1:end) > 0, 2));

frac_ipsi_upper = n_ipsi_upper / (n_ipsi_upper+n_contra_upper);
frac_contra_upper = n_contra_upper / (n_ipsi_upper+n_contra_upper);
frac_ipsi_lower = n_ipsi_lower / (n_ipsi_lower+n_contra_lower);
frac_contra_lower = n_contra_lower / (n_ipsi_lower+n_contra_lower);

ipsi_contra_neurons = [frac_ipsi_upper, frac_ipsi_lower; ...
                       frac_contra_upper, frac_contra_lower];

figure(11); set(gcf, 'position', [488 317 336.2000 345]);
coloritos = [col(3,:); col(4,:)];
h = bar(1:2, ipsi_contra_neurons, 'FaceAlpha', 1, 'EdgeColor', 'none');
for k = 1:length(h)
    h(k).FaceColor = 'flat';
    h(k).CData = repmat(coloritos(k, :), 2, 1);
end
set(gca, 'xticklabel', {"ipsi", "contra"}, ...
    'TickDir', 'out', 'box', 'off', 'FontSize', 14)
ylabel('fraction of projecting neurons', 'FontSize', 14)
xticks(1:2);
ylim([0 1]);
legend('upper', 'lower');
legend('boxoff');

fprintf('Upper: %d/%d ipsi (%.1f%%), %d/%d contra (%.1f%%)\n', ...
    n_ipsi_upper, size(targets_med,1), 100*frac_ipsi_upper, ...
    n_contra_upper, size(targets_med,1), 100*frac_contra_upper);
fprintf('Lower: %d/%d ipsi (%.1f%%), %d/%d contra (%.1f%%)\n', ...
    n_ipsi_lower, size(targets_lat,1), 100*frac_ipsi_lower, ...
    n_contra_lower, size(targets_lat,1), 100*frac_contra_lower);

%% projection strength and probability decomposition

% Binary matrix: whether a neuron projects (non-zero)
proj_med_binary = targets_med > 0;
proj_lat_binary = targets_lat > 0;
proj_ant_binary = targets_ant > 0;
proj_post_binary = targets_post > 0;

proj_med_bin1=targets_med1 > 0;
proj_med_bin2=targets_med2 > 0;
proj_med_bin3=targets_med3 > 0;

proj_lat_bin1=targets_lat1 > 0;
proj_lat_bin2=targets_lat2 > 0;
proj_lat_bin3=targets_lat3 > 0;

% Binary matrix: whether a neuron projects (including zero)
proj_med_binary_all = targets_med_tot > 0;
proj_lat_binary_all = targets_lat_tot > 0;
proj_ant_binary_all = targets_ant_tot > 0;
proj_post_binary_all = targets_post_tot > 0;

proj_med_bin1_all=targets_med1_tot > 0;
proj_med_bin2_all=targets_med2_tot > 0;
proj_med_bin3_all=targets_med3_tot > 0;

proj_lat_bin1_all=targets_lat1_tot > 0;
proj_lat_bin2_all=targets_lat2_tot > 0;
proj_lat_bin3_all=targets_lat3_tot > 0;

% Probability = fraction only of projecting neurons
prob_med = mean(proj_med_binary, 1);
prob_lat = mean(proj_lat_binary, 1);
prob_ant = mean(proj_ant_binary, 1);
prob_post = mean(proj_post_binary, 1);

prob_med1 = mean(proj_med_bin1, 1);
prob_med2 = mean(proj_med_bin2, 1);
prob_med3 = mean(proj_med_bin3, 1);

prob_lat1 = mean(proj_lat_bin1, 1);
prob_lat2 = mean(proj_lat_bin2, 1);
prob_lat3 = mean(proj_lat_bin3, 1);

%Probability = fraction of projecting neurons out of all neurons
prob_med_abs = mean(proj_med_binary_all, 1);
prob_lat_abs = mean(proj_lat_binary_all, 1);
prob_ant_abs = mean(proj_ant_binary_all, 1);
prob_post_abs = mean(proj_post_binary_all, 1);

prob_med_abs1 = mean(proj_med_bin1_all, 1);
prob_med_abs2 = mean(proj_med_bin2_all, 1);
prob_med_abs3 = mean(proj_med_bin3_all, 1);

prob_lat_abs1 = mean(proj_lat_bin1_all, 1);
prob_lat_abs2 = mean(proj_lat_bin2_all, 1);
prob_lat_abs3 = mean(proj_lat_bin3_all, 1);


%mean innervation strength of projecting neurons only
mean_innerv_med = arrayfun(@(c) mean(targets_med(targets_med(:,c) ~= 0, c)), 1:size(targets_med,2));
% mean_innerv_med=mean_innerv_med./sum(mean_innerv_med);
mean_innerv_lat = arrayfun(@(c) mean(targets_lat(targets_lat(:,c) ~= 0, c)), 1:size(targets_lat,2));
% mean_innerv_lat=mean_innerv_lat./sum(mean_innerv_lat);
mean_innerv_ant = arrayfun(@(c) mean(targets_ant(targets_ant(:,c) ~= 0, c)), 1:size(targets_ant,2));
mean_innerv_post = arrayfun(@(c) mean(targets_post(targets_post(:,c) ~= 0, c)), 1:size(targets_post,2));

mean_innerv_med1 = arrayfun(@(c) mean(targets_med1(targets_med1(:,c) ~= 0, c)), 1:size(targets_med1,2));
mean_innerv_med2 = arrayfun(@(c) mean(targets_med2(targets_med2(:,c) ~= 0, c)), 1:size(targets_med2,2));
mean_innerv_med3 = arrayfun(@(c) mean(targets_med3(targets_med3(:,c) ~= 0, c)), 1:size(targets_med3,2));

mean_innerv_lat1 = arrayfun(@(c) mean(targets_lat1(targets_lat1(:,c) ~= 0, c)), 1:size(targets_lat1,2));
mean_innerv_lat2 = arrayfun(@(c) mean(targets_lat2(targets_lat2(:,c) ~= 0, c)), 1:size(targets_lat2,2));
mean_innerv_lat3 = arrayfun(@(c) mean(targets_lat3(targets_lat3(:,c) ~= 0, c)), 1:size(targets_lat3,2));

%mean overall projection strength (this is mean of neurons, and equal to prob * innervation strength)
mean_strength_med = mean(targets_med, 1);
mean_strength_med=mean_strength_med./sum(mean_strength_med);
mean_strength_lat = mean(targets_lat, 1);
mean_strength_lat=mean_strength_lat./sum(mean_strength_lat);

mean_strength_med1 = mean(targets_med1, 1);
mean_strength_med1=mean_strength_med1./sum(mean_strength_med1);
mean_strength_med2 = mean(targets_med2, 1);
mean_strength_med2=mean_strength_med2./sum(mean_strength_med2);
mean_strength_med3 = mean(targets_med3, 1);
mean_strength_med3=mean_strength_med3./sum(mean_strength_med3);

mean_strength_lat1 = mean(targets_lat1, 1);
mean_strength_lat1=mean_strength_lat1./sum(mean_strength_lat1);
mean_strength_lat2 = mean(targets_lat2, 1);
mean_strength_lat2=mean_strength_lat2./sum(mean_strength_lat2);
mean_strength_lat3 = mean(targets_lat3, 1);
mean_strength_lat3=mean_strength_lat3./sum(mean_strength_lat3);


% prob * mean innervation
mean_strength_med_v2 = prob_med_abs .* mean_innerv_med;
mean_strength_lat_v2 = prob_lat_abs .* mean_innerv_lat;
mean_strength_ant_v2 = prob_ant_abs .* mean_innerv_ant;
mean_strength_post_v2 = prob_post_abs .* mean_innerv_post;
% mean_med_both=mean_strength_med_v2./sum(mean_strength_med_v2+mean_strength_lat_v2, 'omitnan');
% mean_lat_both=mean_strength_lat_v2./sum(mean_strength_med_v2+mean_strength_lat_v2, 'omitnan');
% mean_ant_both=mean_strength_ant_v2./sum(mean_strength_ant_v2+mean_strength_post_v2, 'omitnan');
% mean_post_both=mean_strength_post_v2./sum(mean_strength_ant_v2+mean_strength_post_v2, 'omitnan');

mean_med_both=mean_strength_med_v2./sum(mean_strength_med_v2, 'omitnan');
mean_lat_both=mean_strength_lat_v2./sum(mean_strength_lat_v2, 'omitnan');
mean_ant_both=mean_strength_ant_v2./sum(mean_strength_ant_v2, 'omitnan');
mean_post_both=mean_strength_post_v2./sum(mean_strength_post_v2, 'omitnan');

mean_strength_med1_v2 = prob_med_abs1 .* mean_innerv_med1;
mean_strength_med2_v2 = prob_med_abs2 .* mean_innerv_med2;
mean_strength_med3_v2 = prob_med_abs3 .* mean_innerv_med3;

mean_strength_lat1_v2 = prob_lat_abs1 .* mean_innerv_lat1;
mean_strength_lat2_v2 = prob_lat_abs2 .* mean_innerv_lat2;
mean_strength_lat3_v2 = prob_lat_abs3 .* mean_innerv_lat3;

mean_med1_both=mean_strength_med1_v2./sum(mean_strength_med1_v2+mean_strength_lat1_v2, 'omitnan');
mean_med2_both=mean_strength_med2_v2./sum(mean_strength_med2_v2+mean_strength_lat2_v2, 'omitnan');
mean_med3_both=mean_strength_med3_v2./sum(mean_strength_med3_v2+mean_strength_lat3_v2, 'omitnan');

mean_lat1_both=mean_strength_lat1_v2./sum(mean_strength_med1_v2+mean_strength_lat1_v2, 'omitnan');
mean_lat2_both=mean_strength_lat2_v2./sum(mean_strength_med2_v2+mean_strength_lat2_v2, 'omitnan');
mean_lat3_both=mean_strength_lat3_v2./sum(mean_strength_med3_v2+mean_strength_lat3_v2, 'omitnan');
%
% mean_med_both=mean_strength_med_v2./sum(mean_strength_med_v2+mean_strength_lat_v2+...
%     mean_strength_ant_v2+mean_strength_post_v2);
% mean_lat_both=mean_strength_lat_v2./sum(mean_strength_med_v2+mean_strength_lat_v2+...
%     mean_strength_ant_v2+mean_strength_post_v2);
% mean_ant_both=mean_strength_ant_v2./sum(mean_strength_med_v2+mean_strength_lat_v2+...
%     mean_strength_ant_v2+mean_strength_post_v2);
% mean_post_both=mean_strength_post_v2./sum(mean_strength_med_v2+mean_strength_lat_v2+...
%     mean_strength_ant_v2+mean_strength_post_v2);


% Compare binomial proportions to check if probs are different
p_prob=[]; h_prob=[]; stats_prob={};
for i=1:length(projTargets)
    [p_prob(i),h_prob(i),stats_prob{i}] = fishertest([sum(proj_med_binary_all(:,i)), sum(~proj_med_binary_all(:,i));
                          sum(proj_lat_binary_all(:,i)), sum(~proj_lat_binary_all(:,i))], 'Alpha',0.01);

end

p_prob_rank=[]; h_prob_rank=[]; stats_prob_rank={};
for i=1:length(projTargets)
    [p_prob_rank(i),h_prob_rank(i),stats_prob_rank{i}] = ranksum(proj_med_binary_all(:,i), proj_lat_binary_all(:,i), 'alpha', 0.01);

end

%check innervation strengths (only projecting neurons)
p_innerv=[]; h_innerv=[]; stats_innerv={};
for i = 1:length(projTargets)
    new_meds=targets_med(:,i);
    new_meds=new_meds(new_meds~=0);
    new_lats=targets_lat(:,i);
    new_lats=new_lats(new_lats~=0);
    [p_innerv(i), h_innerv(i), stats_innerv{i}] = ranksum(new_meds, new_lats, 'alpha', 0.01,...
        'tail', 'both');
end

% %check projection strengths (projecting and non-projecting neurons)
p_strength=[]; h_strength=[]; stats_strength={};
for i = 1:length(projTargets)
    [p_strength(i), h_strength(i), stats_strength{i}] = ranksum(targets_ant_tot(:,i), targets_post_tot(:,i), 'alpha', 0.01);
end

% p_strength=[]; h_strength=[]; stats_strength={};
% h_norm1=[]; p_norm1=[];
% h_norm2=[]; p_norm2=[];
% for i = 1:length(projTargets)
%     [h_strength(i), p_strength(i), stats_strength{i}] = ttest2(targets_ant_tot(:,i), targets_post_tot(:,i), 'Alpha', 0.01, 'Vartype', 'unequal');
%     [h_norm1(i), p_norm1(i)] = swtest(targets_med_tot(:,i));
%     [h_norm2(i), p_norm2(i)] = swtest(targets_lat_tot(:,i));
% end

%% PROJECTION PROBABILITY ANALYSIS (Bootstrap)
targets = sample_names(projTargets);
nTargets = numel(targets);

% Data: successes and totals
a_med = sum(proj_med_binary_all, 1);
n_med = size(proj_med_binary_all, 1);
a_lat = sum(proj_lat_binary_all, 1);
n_lat = size(proj_lat_binary_all, 1);

% Observed probabilities
obs_p_med = a_med / n_med;
obs_p_lat = a_lat / n_lat;

% BOOTSTRAP
nBootstrap = 10000;

bootstrap_prob_med = zeros(nBootstrap, nTargets);
bootstrap_prob_lat = zeros(nBootstrap, nTargets);

for i = 1:nTargets
    bootstrap_prob_med(:,i) = bootstrp(nBootstrap, @mean, proj_med_binary_all(:,i));
    bootstrap_prob_lat(:,i) = bootstrp(nBootstrap, @mean, proj_lat_binary_all(:,i));
end

bootstrap_log2fold = log2((bootstrap_prob_med + eps) ./ (bootstrap_prob_lat + eps));
bootstrap_prob_diff = bootstrap_prob_med - bootstrap_prob_lat;

median_log2fold = median(bootstrap_log2fold, 1);
median_prob_diff = median(bootstrap_prob_diff, 1);
median_p_med = median(bootstrap_prob_med, 1);
median_p_lat = median(bootstrap_prob_lat, 1);

% MULTIPLE COMPARISON CORRECTION
% Toggle: 'none', 'fcr', or 'bonferroni'
correction_method = 'fcr';

% Always compute uncorrected CIs
ci_log2fold_uncorr = quantile(bootstrap_log2fold, [0.025, 0.975], 1);
ci_prob_diff_uncorr = quantile(bootstrap_prob_diff, [0.025, 0.975], 1);
ci_med_uncorr = quantile(bootstrap_prob_med, [0.025, 0.975])';
ci_lat_uncorr = quantile(bootstrap_prob_lat, [0.025, 0.975])';

if strcmp(correction_method, 'none')
    ci_log2fold = ci_log2fold_uncorr;
    ci_prob_diff = ci_prob_diff_uncorr;
    ci_med = ci_med_uncorr;
    ci_lat = ci_lat_uncorr;
    correction_label = 'uncorrected';

elseif strcmp(correction_method, 'bonferroni')
    alpha_corrected = 0.05 / nTargets;
    ci_log2fold = quantile(bootstrap_log2fold, [alpha_corrected/2, 1-alpha_corrected/2], 1);
    ci_prob_diff = quantile(bootstrap_prob_diff, [alpha_corrected/2, 1-alpha_corrected/2], 1);
    ci_med = quantile(bootstrap_prob_med, [alpha_corrected/2, 1-alpha_corrected/2])';
    ci_lat = quantile(bootstrap_prob_lat, [alpha_corrected/2, 1-alpha_corrected/2])';
    correction_label = sprintf('Bonferroni (alpha=%.3f)', alpha_corrected);

elseif strcmp(correction_method, 'fcr')
    % Compute uncorrected p-values for ranking
    p_values = zeros(1, nTargets);
    for i = 1:nTargets
        p_values(i) = 2 * min(mean(bootstrap_log2fold(:,i) <= 0), mean(bootstrap_log2fold(:,i) >= 0));
        p_values(i) = max(p_values(i), 1/nBootstrap);
    end

    [~, p_rank_idx] = sort(p_values);
    ci_log2fold = zeros(2, nTargets);
    ci_prob_diff = zeros(2, nTargets);
    ci_med = zeros(nTargets, 2);
    ci_lat = zeros(nTargets, 2);

    for rank = 1:nTargets
        idx = p_rank_idx(rank);
        alpha_adj = 0.05 * rank / nTargets;
        ci_log2fold(:, idx) = quantile(bootstrap_log2fold(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_prob_diff(:, idx) = quantile(bootstrap_prob_diff(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_med(idx, :) = quantile(bootstrap_prob_med(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_lat(idx, :) = quantile(bootstrap_prob_lat(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
    end
    correction_label = 'FCR-adjusted';
end

fprintf('Correction method: %s\n\n', correction_label);

% PRINT RESULTS
fprintf('=== PROJECTION PROBABILITY: BOOTSTRAP LOG2 FOLD CHANGE (%s) ===\n', correction_label);
fprintf('%-12s %12s %12s %15s %20s %8s\n', 'Target', 'P(upper)', 'P(lower)', 'Median log2FC', 'CI', 'Sig');
fprintf('%s\n', repmat('-', 1, 85));
for i = 1:nTargets
    sig = '';
    if ci_log2fold(1,i) > 0 || ci_log2fold(2,i) < 0
        sig = '*';
    end
    fprintf('%-12s %11.3f %11.3f %14.3f [%7.3f, %7.3f] %s\n', ...
        targets{i}, obs_p_med(i), obs_p_lat(i), ...
        median_log2fold(i), ci_log2fold(1,i), ci_log2fold(2,i), sig);
end
fprintf('\n* = CI excludes zero (%s)\n', correction_label);

fprintf('\n=== PROBABILITY DIFFERENCES (%s) ===\n', correction_label);
fprintf('%-12s %15s %20s %8s\n', 'Target', 'Median diff', 'CI', 'Sig');
fprintf('%s\n', repmat('-', 1, 60));
for i = 1:nTargets
    sig = '';
    if ci_prob_diff(1,i) > 0 || ci_prob_diff(2,i) < 0
        sig = '*';
    end
    fprintf('%-12s %14.3f [%7.3f, %7.3f] %s\n', ...
        targets{i}, median_prob_diff(i), ci_prob_diff(1,i), ci_prob_diff(2,i), sig);
end

% FIGURE: Projection probabilities and log2 fold change
figure('Position', [100 124 500 550]);

% Top panel: probabilities
subplot(2,1,1); hold on;
x = 1:nTargets;
bar_width = 0.25;
b1 = bar(x - bar_width/2 - 0.02, median_p_med, bar_width, 'FaceColor', col(3,:), 'EdgeColor', 'none');
b2 = bar(x + bar_width/2 + 0.02, median_p_lat, bar_width, 'FaceColor', col(4,:), 'EdgeColor', 'none');

errorbar(x - bar_width/2 - 0.02, median_p_med, median_p_med - ci_med(:,1)', ci_med(:,2)' - median_p_med, ...
    'Color', [0.35 0.35 0.35], 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 0);
errorbar(x + bar_width/2 + 0.02, median_p_lat, median_p_lat - ci_lat(:,1)', ci_lat(:,2)' - median_p_lat, ...
    'Color', [0.35 0.35 0.35], 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 0);

xticks(x); xticklabels(targets);
ylabel('Projection probability');
ylim([0 max([median_p_med median_p_lat]) * 1.3]);
xlim([0.3 nTargets + 0.7]);
legend('upper', 'lower', 'Location', 'northeast');
legend('boxoff');
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12);

% Bottom panel: log2 fold change with density curves
subplot(2,1,2); hold on;
curve_width = 0.3;

for i = 1:nTargets
    data_i = bootstrap_log2fold(:,i);
    data_trimmed = data_i(data_i >= quantile(data_i, 0.0025) & data_i <= quantile(data_i, 0.9975));

    [f, xi] = ksdensity(data_trimmed, 'NumPoints', 100, 'BoundaryCorrection', 'reflection');
    f_normalized = f / max(f) * curve_width;
    x_curve = x(i) + f_normalized;

    fill([x(i), x_curve, x(i)], [xi(1), xi, xi(end)], [0.55 0.55 0.55], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none');

    plot([x(i), x(i)], [ci_log2fold(1,i), ci_log2fold(2,i)], ...
        'k-', 'LineWidth', 2.5);

    plot(x(i), median_log2fold(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
end

yline(0, '--k', 'LineWidth', 1);
xticks(x); xticklabels(targets);
ylabel('log_2 fold change');
xlim([0.3 nTargets + 0.7]);
ylim([-4 4]);
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12);

% FIGURE: Probability differences with density curves
figure('Position', [150 200 500 350]);
hold on;
curve_width = 0.3;

for i = 1:nTargets
    data_i = bootstrap_prob_diff(:,i);
    data_trimmed = data_i(data_i >= quantile(data_i, 0.0025) & data_i <= quantile(data_i, 0.9975));

    [f, xi] = ksdensity(data_trimmed, 'NumPoints', 100, 'BoundaryCorrection', 'reflection');
    f_normalized = f / max(f) * curve_width;
    x_curve = x(i) + f_normalized;

    fill([x(i), x_curve, x(i)], [xi(1), xi, xi(end)], [0.55 0.55 0.55], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none');

    plot([x(i), x(i)], [ci_prob_diff(1,i), ci_prob_diff(2,i)], ...
        'k-', 'LineWidth', 2.5);

    plot(x(i), median_prob_diff(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
end

yline(0, '--k', 'LineWidth', 1);
xticks(x); xticklabels(targets);
ylabel('Probability difference (upper - lower)');
xlim([0.3 nTargets + 0.7]);
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12);

fprintf('\n=== ANALYSIS COMPLETE ===\n');

% FIGURE: Pie charts of projection percentages
pct_med = (a_med / n_med) * 100;
pct_lat = (a_lat / n_lat) * 100;

n_proj_med = sum(any(proj_med_binary_all, 2));
n_proj_lat = sum(any(proj_lat_binary_all, 2));

figure('Position', [200 250 900 400]);

% Upper field pie chart
subplot(1,2,1);
p1 = pie(pct_med);
title(sprintf('Upper sSC (n=%d)', n_proj_med), 'FontSize', 14);
colors_upper = flipud(cbrewer('seq', 'Oranges', nTargets+2));
colors_upper = colors_upper(2:end, :);
for i = 1:nTargets
    set(p1(2*i-1), 'FaceColor', colors_upper(i,:), 'EdgeColor', 'w', 'LineWidth', 1, 'FaceAlpha', 0.7);
end
legend(targets, 'Location', 'eastoutside', 'FontSize', 11);
legend('boxoff');

% Lower field pie chart
subplot(1,2,2);
p2 = pie(pct_lat);
title(sprintf('Lower sSC (n=%d)', n_proj_lat), 'FontSize', 14);
colors_lower = flipud(cbrewer('seq', 'Blues', nTargets+1));
for i = 1:nTargets
    set(p2(2*i-1), 'FaceColor', colors_lower(i,:), 'EdgeColor', 'w', 'LineWidth', 1, 'FaceAlpha', 0.7);
end
legend(targets, 'Location', 'eastoutside', 'FontSize', 11);
legend('boxoff');
%% plot innervation strength

names=sample_names(projTargets);
tot_targets={targets_med, targets_lat};
area=1:length(projTargets);
figure(20); set(gcf, 'position', [94.6000 381.8000 1.3136e+03 360.8000]);
til=tiledlayout(1,5, 'TileSpacing', 'loose', 'Padding','compact');
set(til, 'OuterPosition', [0 0 0.98 0.98]);
coloritos=[hex2rgb("#E8AE7D");hex2rgb("#6D9CC5");hex2rgb("#D06D74");hex2rgb("#38A68A")];
for a = 1:length(projTargets)
    nexttile;
    title(names{a}, 'FontWeight', 'normal', 'FontAngle', 'italic')
    hold on
    for t=1:2
        temp_targets=tot_targets{t}(:,a);
        new_targets = temp_targets(temp_targets ~= 0);

        x_vector=[0.5 2];
        h=boxchart(repmat(x_vector(t),size(new_targets,1),1),new_targets,'BoxFaceColor',coloritos(t,:), 'BoxWidth',0.3,...
            'BoxMedianLineColor',coloritos(t,:), 'LineWidth', 0.9,'MarkerStyle','none',...
            'WhiskerLineColor', [0.5 0.5 0.5],'WhiskerLineStyle','-',...
            'BoxFaceAlpha',0.9, 'BoxEdgeColor','none'); hold on;
        med = median(new_targets); % skip zeros if needed
        me_an= mean(new_targets);

        % Get box x-position (for ungrouped boxchart, it's 1)
        xpos = h.XData(1);  % usually 1 unless you're plotting groups

        % Manually draw thicker median line and mean dot
        hold on
        plot([xpos-0.15, xpos+0.15], [med, med], '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
        plot(xpos, me_an, 'o', 'MarkerSize',4.5, 'MarkerEdgeColor', 'none',...
            'MarkerFaceColor', [0.2 0.2 0.2]);

        swarm=swarmchart((repmat(x_vector(t),size(new_targets,1),1)+0.6),new_targets', 10,...
            'filled','MarkerFaceColor',coloritos(t,:), 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);%plot swarmchart
        swarm.XJitterWidth=0.5;
        uistack(swarm, 'bottom');
        hold on;

    end
    set(gca, 'xtick', [0.8 2.3], 'xticklabel', {'med', 'lat'}, 'FontSize', 13);
    ylabel({'Innervation strength'}); ylim([0 0.17]);
    ax=gca;ax.YAxis.FontSize=13;
end


% %check conditional strength
% p_cond=[]; h_cond=[]; stats_cond={};
% for i = 1:length(projTargets)
%     [p_cond(i), h_cond(i), stats_cond{i}] = ranksum(cond_strength_med(:,i), cond_strength_lat(:,i),'alpha', 0.01);
% end


%% INNERVATION STRENGTH ANALYSIS (Bootstrap)
% Compares mean innervation strength (among projecting neurons only)
% between upper and lower field sSC populations

targets = sample_names(projTargets);
nTargets = length(targets);

% BOOTSTRAP
nBootstrap = 10000;

bootstrap_mean_med = zeros(nBootstrap, nTargets);
bootstrap_mean_lat = zeros(nBootstrap, nTargets);

fprintf('=== INNERVATION STRENGTH BOOTSTRAP ===\n');
for i = 1:nTargets
    % Non-zero neurons only
    data_med_i = targets_med(:,i);
    data_med_i = data_med_i(data_med_i ~= 0);

    data_lat_i = targets_lat(:,i);
    data_lat_i = data_lat_i(data_lat_i ~= 0);

    fprintf('Target %s: n_med=%d, n_lat=%d\n', targets{i}, length(data_med_i), length(data_lat_i));

    if ~isempty(data_med_i)
        bootstrap_mean_med(:,i) = bootstrp(nBootstrap, @mean, data_med_i);
    else
        bootstrap_mean_med(:,i) = NaN;
    end

    if ~isempty(data_lat_i)
        bootstrap_mean_lat(:,i) = bootstrp(nBootstrap, @mean, data_lat_i);
    else
        bootstrap_mean_lat(:,i) = NaN;
    end
end

bootstrap_log2fold = log2((bootstrap_mean_med + eps) ./ (bootstrap_mean_lat + eps));
bootstrap_diff = bootstrap_mean_med - bootstrap_mean_lat;

median_log2fold = median(bootstrap_log2fold, 1);
median_diff = median(bootstrap_diff, 1);
median_mean_med = median(bootstrap_mean_med, 1);
median_mean_lat = median(bootstrap_mean_lat, 1);

% MULTIPLE COMPARISON CORRECTION
% Toggle: 'none', 'fcr', or 'bonferroni'
correction_method = 'fcr';

% Always compute uncorrected CIs
ci_log2fold_uncorr = quantile(bootstrap_log2fold, [0.025, 0.975], 1);
ci_diff_uncorr = quantile(bootstrap_diff, [0.025, 0.975], 1);
ci_med_uncorr = quantile(bootstrap_mean_med, [0.025, 0.975])';
ci_lat_uncorr = quantile(bootstrap_mean_lat, [0.025, 0.975])';

if strcmp(correction_method, 'none')
    ci_log2fold = ci_log2fold_uncorr;
    ci_diff = ci_diff_uncorr;
    ci_med = ci_med_uncorr;
    ci_lat = ci_lat_uncorr;
    correction_label = 'uncorrected';

elseif strcmp(correction_method, 'bonferroni')
    alpha_corrected = 0.05 / nTargets;
    ci_log2fold = quantile(bootstrap_log2fold, [alpha_corrected/2, 1-alpha_corrected/2], 1);
    ci_diff = quantile(bootstrap_diff, [alpha_corrected/2, 1-alpha_corrected/2], 1);
    ci_med = quantile(bootstrap_mean_med, [alpha_corrected/2, 1-alpha_corrected/2])';
    ci_lat = quantile(bootstrap_mean_lat, [alpha_corrected/2, 1-alpha_corrected/2])';
    correction_label = sprintf('Bonferroni (alpha=%.3f)', alpha_corrected);

elseif strcmp(correction_method, 'fcr')
    p_values = zeros(1, nTargets);
    for i = 1:nTargets
        p_values(i) = 2 * min(mean(bootstrap_log2fold(:,i) <= 0), mean(bootstrap_log2fold(:,i) >= 0));
        p_values(i) = max(p_values(i), 1/nBootstrap);
    end

    [~, p_rank_idx] = sort(p_values);
    ci_log2fold = zeros(2, nTargets);
    ci_diff = zeros(2, nTargets);
    ci_med = zeros(nTargets, 2);
    ci_lat = zeros(nTargets, 2);

    for rank = 1:nTargets
        idx = p_rank_idx(rank);
        alpha_adj = 0.05 * rank / nTargets;
        ci_log2fold(:, idx) = quantile(bootstrap_log2fold(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_diff(:, idx) = quantile(bootstrap_diff(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_med(idx, :) = quantile(bootstrap_mean_med(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_lat(idx, :) = quantile(bootstrap_mean_lat(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
    end
    correction_label = 'FCR-adjusted';
end

fprintf('\nCorrection method: %s\n\n', correction_label);

% PRINT RESULTS
fprintf('=== INNERVATION STRENGTH: BOOTSTRAP LOG2 FOLD CHANGE (%s) ===\n', correction_label);
fprintf('%-12s %12s %12s %15s %20s %8s\n', 'Target', 'Med mean', 'Lat mean', 'Median log2FC', 'CI', 'Sig');
fprintf('%s\n', repmat('-', 1, 85));
for i = 1:nTargets
    sig = '';
    if ci_log2fold(1,i) > 0 || ci_log2fold(2,i) < 0
        sig = '*';
    end
    fprintf('%-12s %11.4f %11.4f %14.3f [%7.3f, %7.3f] %s\n', ...
        targets{i}, median_mean_med(i), median_mean_lat(i), ...
        median_log2fold(i), ci_log2fold(1,i), ci_log2fold(2,i), sig);
end
fprintf('\n* = CI excludes zero (%s)\n', correction_label);

% FIGURE: Innervation strength with swarm + log2 fold change
figure('Position', [165 374.6000 669.6000 378.4000], 'Color', 'w');
x = 1:nTargets;

% Top panel: innervation strengths with swarm plots
subplot(2,1,1); hold on;

errorbar(x - 0.22, median_mean_med, ...
    median_mean_med - ci_med(:,1)', ci_med(:,2)' - median_mean_med, 'o', ...
    'Color', col(3,:)-0.2, 'MarkerFaceColor', col(3,:)-0.2, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 5, 'LineWidth', 2, 'CapSize', 0);
errorbar(x + 0.08, median_mean_lat, ...
    median_mean_lat - ci_lat(:,1)', ci_lat(:,2)' - median_mean_lat, 'o', ...
    'Color', col(4,:)-0.2, 'MarkerFaceColor', col(4,:)-0.2, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 5, 'LineWidth', 2, 'CapSize', 0);

% Swarm plots for medial
for a = 1:nTargets
    temp_targets = targets_med(:,a);
    new_targets = temp_targets(temp_targets ~= 0);
    if ~isempty(new_targets)
        xvector = x(a) - 0.1;
        swarm = swarmchart(repmat(xvector, length(new_targets), 1), new_targets, 10, ...
            'filled', 'MarkerFaceColor', col(3,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        swarm.XJitterWidth = 0.15;
        uistack(swarm, 'bottom');
    end
end

% Swarm plots for lateral
for a = 1:nTargets
    temp_targets = targets_lat(:,a);
    new_targets = temp_targets(temp_targets ~= 0);
    if ~isempty(new_targets)
        xvector = x(a) + 0.2;
        swarm = swarmchart(repmat(xvector, length(new_targets), 1), new_targets, 10, ...
            'filled', 'MarkerFaceColor', col(4,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
        swarm.XJitterWidth = 0.15;
        uistack(swarm, 'bottom');
    end
end

xticks(x); xticklabels(targets);
ylim([0 0.025]);
xlim([0.5 nTargets+0.5]);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('Innervation strength', 'FontSize', 12);

% Bottom panel: log2 fold change with density curves
subplot(2,1,2); hold on;
curve_width = 0.3;

for i = 1:nTargets
    data_i = bootstrap_log2fold(:,i);
    data_trimmed = data_i(data_i >= quantile(data_i, 0.0025) & data_i <= quantile(data_i, 0.9975));

    [f, xi] = ksdensity(data_trimmed, 'NumPoints', 100, 'BoundaryCorrection', 'reflection');
    f_normalized = f / max(f) * curve_width;
    x_curve = x(i) + f_normalized;

    fill([x(i), x_curve, x(i)], [xi(1), xi, xi(end)], [0.55 0.55 0.55], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none');

    plot([x(i), x(i)], [ci_log2fold(1,i), ci_log2fold(2,i)], ...
        'k-', 'LineWidth', 2.5);

    plot(x(i), median_log2fold(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
end

yline(0, '--k', 'LineWidth', 1);
xticks(x); xticklabels(targets);
ylim([-1.5 1.5]);
xlim([0.5 nTargets+0.5]);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('log_2 fold change', 'FontSize', 12);

fprintf('\n=== INNERVATION STRENGTH ANALYSIS COMPLETE ===\n');

%% PROJECTION STRENGTH ANALYSIS (Bootstrap)
% Overall projection strength = probability * innervation strength
% Compares upper vs lower field sSC neurons across targets

coloritos=[hex2rgb("#E8AE7D");hex2rgb("#6D9CC5");hex2rgb("#D06D74");hex2rgb("#38A68A")];

nTargets = length(projTargets);
target_names = sample_names(projTargets);


data_med = targets_med_tot;
data_lat = targets_lat_tot;

[n_med, ~] = size(data_med);
[n_lat, ~] = size(data_lat);

fprintf('=== DATA SUMMARY ===\n');
fprintf('Upper field neurons: %d\n', n_med);
fprintf('Lower field neurons: %d\n', n_lat);
fprintf('Targets: %d\n\n', nTargets);

for i = 1:nTargets
    fprintf('%s: Med mean=%.3f (%.1f%% non-zero), Lat mean=%.3f (%.1f%% non-zero)\n', ...
        target_names{i}, mean(data_med(:,i)), 100*mean(data_med(:,i)>0), ...
        mean(data_lat(:,i)), 100*mean(data_lat(:,i)>0));
end

% BOOTSTRAP
nBootstrap = 10000;

bootstrap_mean_med = zeros(nBootstrap, nTargets);
bootstrap_mean_lat = zeros(nBootstrap, nTargets);

for i = 1:nTargets
    bootstrap_mean_med(:,i) = bootstrp(nBootstrap, @mean, data_med(:,i));
    bootstrap_mean_lat(:,i) = bootstrp(nBootstrap, @mean, data_lat(:,i));
end

bootstrap_log2fold = log2((bootstrap_mean_med + eps) ./ (bootstrap_mean_lat + eps));

% Fractional allocation per bootstrap sample
bootstrap_frac_med = bootstrap_mean_med ./ sum(bootstrap_mean_med, 2);
bootstrap_frac_lat = bootstrap_mean_lat ./ sum(bootstrap_mean_lat, 2);
bootstrap_frac_diff = bootstrap_frac_med - bootstrap_frac_lat;

median_log2fold = median(bootstrap_log2fold, 1);
median_mean_med = median(bootstrap_mean_med, 1);
median_mean_lat = median(bootstrap_mean_lat, 1);
median_frac_med = median(bootstrap_frac_med, 1);
median_frac_lat = median(bootstrap_frac_lat, 1);
median_frac_diff = median(bootstrap_frac_diff, 1);

% MULTIPLE COMPARISON CORRECTION
% Toggle: 'none', 'fcr', or 'bonferroni'
correction_method = 'fcr';

% Uncorrected CIs
ci_log2fold_uncorr = quantile(bootstrap_log2fold, [0.025, 0.975], 1);
ci_med_uncorr = quantile(bootstrap_mean_med, [0.025, 0.975])';
ci_lat_uncorr = quantile(bootstrap_mean_lat, [0.025, 0.975])';
ci_frac_med_uncorr = quantile(bootstrap_frac_med, [0.025, 0.975])';
ci_frac_lat_uncorr = quantile(bootstrap_frac_lat, [0.025, 0.975])';
ci_frac_diff_uncorr = quantile(bootstrap_frac_diff, [0.025, 0.975], 1);

if strcmp(correction_method, 'none')
    ci_log2fold = ci_log2fold_uncorr;
    ci_med = ci_med_uncorr;
    ci_lat = ci_lat_uncorr;
    ci_frac_med = ci_frac_med_uncorr;
    ci_frac_lat = ci_frac_lat_uncorr;
    ci_frac_diff = ci_frac_diff_uncorr;
    correction_label = 'uncorrected';

elseif strcmp(correction_method, 'bonferroni')
    alpha_corrected = 0.05 / nTargets;
    ci_log2fold = quantile(bootstrap_log2fold, [alpha_corrected/2, 1-alpha_corrected/2], 1);
    ci_med = quantile(bootstrap_mean_med, [alpha_corrected/2, 1-alpha_corrected/2])';
    ci_lat = quantile(bootstrap_mean_lat, [alpha_corrected/2, 1-alpha_corrected/2])';
    ci_frac_med = quantile(bootstrap_frac_med, [alpha_corrected/2, 1-alpha_corrected/2])';
    ci_frac_lat = quantile(bootstrap_frac_lat, [alpha_corrected/2, 1-alpha_corrected/2])';
    ci_frac_diff = quantile(bootstrap_frac_diff, [alpha_corrected/2, 1-alpha_corrected/2], 1);
    correction_label = sprintf('Bonferroni (alpha=%.3f)', alpha_corrected);

elseif strcmp(correction_method, 'fcr')
    p_values = zeros(1, nTargets);
    for i = 1:nTargets
        p_values(i) = 2 * min(mean(bootstrap_log2fold(:,i) <= 0), mean(bootstrap_log2fold(:,i) >= 0));
        p_values(i) = max(p_values(i), 1/nBootstrap);
    end

    [~, p_rank_idx] = sort(p_values);
    ci_log2fold = zeros(2, nTargets);
    ci_med = zeros(nTargets, 2);
    ci_lat = zeros(nTargets, 2);
    ci_frac_med = zeros(nTargets, 2);
    ci_frac_lat = zeros(nTargets, 2);
    ci_frac_diff = zeros(2, nTargets);

    for rank = 1:nTargets
        idx = p_rank_idx(rank);
        alpha_adj = 0.05 * rank / nTargets;
        ci_log2fold(:, idx) = quantile(bootstrap_log2fold(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_med(idx, :) = quantile(bootstrap_mean_med(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_lat(idx, :) = quantile(bootstrap_mean_lat(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_frac_med(idx, :) = quantile(bootstrap_frac_med(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_frac_lat(idx, :) = quantile(bootstrap_frac_lat(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
        ci_frac_diff(:, idx) = quantile(bootstrap_frac_diff(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
    end
    correction_label = 'FCR-adjusted';
end

fprintf('\nCorrection method: %s\n\n', correction_label);

% PRINT RESULTS - ABSOLUTE STRENGTH
fprintf('=== PROJECTION STRENGTH: LOG2 FOLD CHANGE (%s) ===\n', correction_label);
fprintf('%-12s %12s %12s %15s %20s %8s\n', 'Target', 'Upper mean', 'Lower mean', 'Median log2FC', 'CI', 'Sig');
fprintf('%s\n', repmat('-', 1, 85));
for i = 1:nTargets
    sig = '';
    if ci_log2fold(1,i) > 0 || ci_log2fold(2,i) < 0
        sig = '*';
    end
    fprintf('%-12s %11.4f %11.4f %14.3f [%7.3f, %7.3f] %s\n', ...
        target_names{i}, median_mean_med(i), median_mean_lat(i), ...
        median_log2fold(i), ci_log2fold(1,i), ci_log2fold(2,i), sig);
end
fprintf('\n* = CI excludes zero (%s)\n', correction_label);

% PRINT RESULTS - FRACTIONAL ALLOCATION
fprintf('\n=== FRACTIONAL ALLOCATION (%s) ===\n', correction_label);
fprintf('%-12s %12s %12s %15s %20s %8s\n', 'Target', 'Upper %%', 'Lower %%', 'Diff (U-L)', 'CI', 'Sig');
fprintf('%s\n', repmat('-', 1, 85));
for i = 1:nTargets
    sig = '';
    if ci_frac_diff(1,i) > 0 || ci_frac_diff(2,i) < 0
        sig = '*';
    end
    fprintf('%-12s %10.1f%% %10.1f%% %13.1f%% [%6.1f%%, %6.1f%%] %s\n', ...
        target_names{i}, 100*median_frac_med(i), 100*median_frac_lat(i), ...
        100*median_frac_diff(i), 100*ci_frac_diff(1,i), 100*ci_frac_diff(2,i), sig);
end

% FIGURE: Ribbon stacked bar + log2 fold change
figure('Position', [100 200 600 550]);
x = 1:nTargets;

% Top panel: ribbon stacked bar
subplot(2,1,1); hold on;

dataMatrix = [median_mean_med; median_mean_lat];
% coloritos = [col(1,:); col(2,:)];
coloritos = [coloritos(3,:); coloritos(4,:)];
[M, N] = size(dataMatrix);
sortedData = zeros(M, N);
colorOrder = zeros(M, N);

for i = 1:N
    [sortedCol, sortIdx] = sort(dataMatrix(:,i), 'ascend');
    sortedData(:,i) = sortedCol;
    colorOrder(:,i) = sortIdx;
end

h = bar(sortedData', 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.3, ...
    'FaceAlpha', 1, 'EdgeColor', 'none');

for k = 1:M
    h(k).FaceColor = 'flat';
    for i = 1:N
        datasetIndex = colorOrder(k,i);
        h(k).CData(i,:) = coloritos(datasetIndex, :);
    end
end

% Ribbons
barWidth = 0.15;
xEdges = [];
for i = 1:N
    xEdges = [xEdges, i-barWidth, i+barWidth];
end

bottom = zeros(M, N);
top = zeros(M, N);
for j = 1:N
    [sortedVals, sortIdx] = sort(dataMatrix(:, j), 'ascend');
    cum = 0;
    for k = 1:M
        idx = sortIdx(k);
        bottom(idx, j) = cum;
        cum = cum + sortedVals(k);
        top(idx, j) = cum;
    end
end

for j = 1:N-1
    xL = xEdges(2*j);
    xR = xEdges(2*j+1);
    for m = 1:M
        [xTop, yTop] = get_curves2(xL, top(m,j), xR, top(m,j+1));
        [xBot, yBot] = get_curves2(xL, bottom(m,j), xR, bottom(m,j+1));
        patchX = [xTop; flipud(xBot)];
        patchY = [yTop; flipud(yBot)];
        ptch = patch(patchX, patchY, coloritos(m,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        uistack(ptch, 'bottom');
    end
end

% Add value labels on bars
cumulativeHeights = zeros(N,1);
for k = 1:M
    for i = 1:N
        val = sortedData(k,i);
        if val >= 0.00005
            text(i, cumulativeHeights(i) + val/2, sprintf('%.2f', val*1000), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 9, 'Color', 'w', 'FontWeight', 'bold');
        end
        cumulativeHeights(i) = cumulativeHeights(i) + val;
    end
end

xticks(1:nTargets); xticklabels(target_names);
xlim([0.5 nTargets+0.5]);
ylim([0 max(sum(dataMatrix,1))*1.1]);
ylabel('avg proj strength per neuron', 'FontSize', 12);
legend([h(1), h(2)], 'upper', 'lower', 'Location', 'best', 'FontSize', 12);
legend('boxoff');
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12);

ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.3f');

% Bottom panel: log2 fold change with density curves
subplot(2,1,2); hold on;
curve_width = 0.3;

for i = 1:nTargets
    data_i = bootstrap_log2fold(:,i);
    data_trimmed = data_i(data_i >= quantile(data_i, 0.0025) & data_i <= quantile(data_i, 0.9975));

    [f, xi] = ksdensity(data_trimmed, 'NumPoints', 100, 'BoundaryCorrection', 'reflection');
    f_normalized = f / max(f) * curve_width;
    x_curve = x(i) + f_normalized;

    fill([x(i), x_curve, x(i)], [xi(1), xi, xi(end)], [0.55 0.55 0.55], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none');

    plot([x(i), x(i)], [ci_log2fold(1,i), ci_log2fold(2,i)], ...
        'k-', 'LineWidth', 2.5);

    plot(x(i), median_log2fold(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
end

yline(0, '--k', 'LineWidth', 1);
xticks(x); xticklabels(target_names);
ylim([-4 4]);
xlim([0.5 nTargets+0.5]);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
ylabel('log_2 fold change', 'FontSize', 12);
set(gcf, 'renderer', 'painters');

% FIGURE: Average projection strength (bar plot with corrected CIs)
figure('Position', [200 200 500 350]);
hold on;
x_pos = 1:nTargets;
bar_width = 0.25;

b1 = bar(x_pos - bar_width/2 - 0.02, median_mean_med, bar_width, 'FaceColor', col(3,:), 'EdgeColor', 'none');
b2 = bar(x_pos + bar_width/2 + 0.02, median_mean_lat, bar_width, 'FaceColor', col(4,:), 'EdgeColor', 'none');

errorbar(x_pos - bar_width/2 - 0.02, median_mean_med, ...
    median_mean_med - ci_med(:,1)', ci_med(:,2)' - median_mean_med, ...
    'Color', [0.35 0.35 0.35], 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 0);
errorbar(x_pos + bar_width/2 + 0.02, median_mean_lat, ...
    median_mean_lat - ci_lat(:,1)', ci_lat(:,2)' - median_mean_lat, ...
    'Color', [0.35 0.35 0.35], 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 0);

xticks(x_pos); xticklabels(target_names);
ylabel('avg projection strength per neuron', 'FontSize', 12);
xlim([0.3 nTargets + 0.7]);
legend([b1, b2], 'upper', 'lower', 'Location', 'best', 'FontSize', 12);
legend('boxoff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

fprintf('\n=== PROJECTION STRENGTH ANALYSIS COMPLETE ===\n');


% FIGURE: Fractional allocation with fold change
figure('Position', [250 250 500 550]);
x_pos = 1:nTargets;

% Top panel: fractional allocation bars
subplot(2,1,1); hold on;
bar_width = 0.25;

b1 = bar(x_pos - bar_width/2 - 0.02, 100*median_frac_med, bar_width, 'FaceColor', col(3,:), 'EdgeColor', 'none');
b2 = bar(x_pos + bar_width/2 + 0.02, 100*median_frac_lat, bar_width, 'FaceColor', col(4,:), 'EdgeColor', 'none');

errorbar(x_pos - bar_width/2 - 0.02, 100*median_frac_med, ...
    100*(median_frac_med - ci_frac_med(:,1)'), 100*(ci_frac_med(:,2)' - median_frac_med), ...
    'Color', [0.35 0.35 0.35], 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 0);
errorbar(x_pos + bar_width/2 + 0.02, 100*median_frac_lat, ...
    100*(median_frac_lat - ci_frac_lat(:,1)'), 100*(ci_frac_lat(:,2)' - median_frac_lat), ...
    'Color', [0.35 0.35 0.35], 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 0);

xticks(x_pos); xticklabels(target_names);
ylabel('% of total output', 'FontSize', 12);
xlim([0.3 nTargets + 0.7]);
legend([b1, b2], 'upper', 'lower', 'Location', 'best', 'FontSize', 12);
legend('boxoff');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);

% Bottom panel: fractional allocation log2 fold change
subplot(2,1,2); hold on;
curve_width = 0.3;

bootstrap_frac_log2fold = log2((bootstrap_frac_med + eps) ./ (bootstrap_frac_lat + eps));
median_frac_log2fold = median(bootstrap_frac_log2fold, 1);

% Compute corrected CIs for fractional log2fold
if strcmp(correction_method, 'none')
    ci_frac_log2fold = quantile(bootstrap_frac_log2fold, [0.025, 0.975], 1);
elseif strcmp(correction_method, 'bonferroni')
    alpha_corrected = 0.05 / nTargets;
    ci_frac_log2fold = quantile(bootstrap_frac_log2fold, [alpha_corrected/2, 1-alpha_corrected/2], 1);
elseif strcmp(correction_method, 'fcr')
    ci_frac_log2fold = zeros(2, nTargets);
    for rank = 1:nTargets
        idx = p_rank_idx(rank);
        alpha_adj = 0.05 * rank / nTargets;
        ci_frac_log2fold(:, idx) = quantile(bootstrap_frac_log2fold(:, idx), [alpha_adj/2, 1-alpha_adj/2]);
    end
end

for i = 1:nTargets
    data_i = bootstrap_frac_log2fold(:,i);
    data_trimmed = data_i(data_i >= quantile(data_i, 0.0025) & data_i <= quantile(data_i, 0.9975));

    [f, xi] = ksdensity(data_trimmed, 'NumPoints', 100, 'BoundaryCorrection', 'reflection');
    f_normalized = f / max(f) * curve_width;
    x_curve = x_pos(i) + f_normalized;

    fill([x_pos(i), x_curve, x_pos(i)], [xi(1), xi, xi(end)], [0.55 0.55 0.55], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none');

    plot([x_pos(i), x_pos(i)], [ci_frac_log2fold(1,i), ci_frac_log2fold(2,i)], ...
        'k-', 'LineWidth', 2.5);

    plot(x_pos(i), median_frac_log2fold(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
end

yline(0, '--k', 'LineWidth', 1);
xticks(x_pos); xticklabels(target_names);
ylabel('log_2 fold change', 'FontSize', 12);
xlim([0.3 nTargets + 0.7]);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12);
%% projection strength and probability - MED LAT

%overall projection strength (mix of projection strength and probability)
% coloritos=[col(3,:);col(4,:)];
coloritos=[hex2rgb("#E8AE7D");hex2rgb("#6D9CC5");hex2rgb("#D06D74");hex2rgb("#38A68A")];
figure(10)
% all_strength_both=[mean_strength_med;mean_strength_lat; mean_strength_ant; mean_strength_post];
all_strength_both=[mean_med_both;mean_lat_both;mean_ant_both;mean_post_both];
h=bar(1:length(projTargets),all_strength_both, 'FaceAlpha', 1, 'EdgeColor', 'none');
for k = 1:length(h)
    h(k).FaceColor = 'flat';
    h(k).CData = repmat(coloritos(k, :), size(all_strength_both, 2), 1);
end
set(gca,'xticklabel',sample_names(projTargets),...
     'TickDir', 'out', 'box', 'off', 'FontSize',14)
% xlabel('Projection targets', 'FontSize', 11)
ylabel('Overall projection strength', 'FontSize', 14)
xticks(1:length(projTargets)); ylim([0 0.45])
legend('upper', 'lower')


% innvervation strength of projecting neurons
figure(11)
strength_projecting_both=[mean_innerv_med;mean_innerv_lat];
h=bar(1:length(projTargets),strength_projecting_both, 'FaceAlpha', 1, 'EdgeColor', 'none');
for k = 1:length(h)
    h(k).FaceColor = 'flat';
    h(k).CData = repmat(coloritos(k, :), size(strength_projecting_both, 2), 1);
end
set(gca,'xticklabel',sample_names(projTargets),...
     'TickDir', 'out', 'box', 'off', 'FontSize',14)
% xlabel('Projection targets', 'FontSize', 11)
ylabel('Mean innervation strength', 'FontSize', 14)
xticks(1:length(projTargets)); %ylim([0 0.08])
legend('upper', 'lower')


%projection probability
figure(12); set(gcf, 'position', [359.4000 333 524.8000 323.2000]);
%normalize projection probability
prob_both=[prob_med_abs;prob_lat_abs];
h=bar(1:length(projTargets),prob_both, 'FaceAlpha', 1, 'EdgeColor', 'none'); hold on
for k = 1:length(h)
    h(k).FaceColor = 'flat';
    h(k).CData = repmat(coloritos(k, :), size(prob_both, 2), 1);
end
set(gca,'xticklabel',sample_names(projTargets),...
     'TickDir', 'out', 'box', 'off', 'FontSize',14)
% xlabel('Projection targets', 'FontSize', 11)
ylabel('Projection probability', 'FontSize', 14)
xticks(1:length(projTargets)); %ylim([0 0.7])
legend('upper', 'lower')

%plot sig star
for i=1:length(projTargets)
    xvector=i;
    % yl=max([log_exp(i) log_obs(i)])+0.3;
    yl=max(prob_both(:,i))+0.005;
    % sigstars on top
    mysigstar(gca, xvector, yl, h_prob(i),'k', 'down', 'contextual');
end
hold off;


%% ribbon plot

%first do stacked
coloritos=[col(3,:);col(4,:);col(1,:);col(2,:)];
% coloritos=[hex2rgb("#E8AE7D");hex2rgb("#6D9CC5");hex2rgb("#D06D74");hex2rgb("#38A68A")];
% coloritos=[hex2rgb("#E8AE7D");hex2rgb("#6D9CC5")];
% coloritos=[hex2rgb("#D06D74");hex2rgb("#38A68A")];
% coloritos=[col(1,:);col(2,:)];
figure(10); set(gcf, 'position', [439.4000 299.4000 619.2000 390.4000])
% dataMatrix=[mean_med_both;mean_lat_both;mean_ant_both;mean_post_both];
dataMatrix=[mean_med_both;mean_lat_both];
% dataMatrix=[mean_strength_med;mean_strength_lat];
% dataMatrix=[mean_ant_both;mean_post_both];
[M, N] = size(dataMatrix);


% Preallocate stack matrix and color index
sortedData = zeros(M, N);     % 2 rows: one for bottom, one for top
colorOrder = zeros(M, N);    % to track which color goes where

% Sort each pair and keep color order
for i = 1:N
    [sortedCol, sortIdx] = sort(dataMatrix(:,i), 'ascend');
    sortedData(:,i) = sortedCol;
    colorOrder(:,i) = sortIdx;  % Track original dataset for coloring
end

h = bar(sortedData','stacked', 'FaceColor', 'flat', 'BarWidth',0.3,...
    'FaceAlpha', 1, 'EdgeColor', 'none'); hold on


for k = 1:M
    h(k).FaceColor = 'flat';
    for i = 1:N
        datasetIndex = colorOrder(k,i);  % which original dataset
        h(k).CData(i,:) = coloritos(datasetIndex, :);
    end
end

% add values text
cumulativeHeights = zeros(N,1);  % One for each bar


for k = 1:M
    for i = 1:N
        val = sortedData(k,i);
        if val >= 0.045  % Only label if tall enough (adjust threshold here)
            xpos = i;
            ypos = cumulativeHeights(i) + val/2;  % Center of the bar segment
            text(xpos, ypos, sprintf('%.2f', val), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10, ...
                'Color', 'w', ...       % white text for contrast
                'FontWeight', 'bold');
        end
        cumulativeHeights(i) = cumulativeHeights(i) + val;  % update stack height
    end
end

% add ribbons as patches
xCenters = 1:N;
barWidth=0.15;
xEdges=[];
for i=1:N
    xEdges= [xEdges, i-barWidth, i+barWidth];
end
bottom = zeros(M, N);
top = zeros(M, N);

for j = 1:N
    [sortedVals, sortIdx] = sort(dataMatrix(:, j), 'ascend');
    cum = 0;
    for k = 1:M
        idx = sortIdx(k);
        bottom(idx, j) = cum;
        cum = cum + sortedVals(k);
        top(idx, j) = cum;
    end
end


for j = 1:N-1  % column pairs: j and j+1
    xL = xEdges(2*j);     % right edge of bar j
    xR = xEdges(2*j+1);   % left edge of bar j+1

    for m = 1:M  % loop over datasets
        % Y-values: follow this dataset's segment
        y1_bot = bottom(m, j);
        y1_top = top(m, j);
        y2_bot = bottom(m, j+1);
        y2_top = top(m, j+1);

        [xTop, yTop] = get_curves(xL, y1_top, xR, y2_top);
        [xBot, yBot] = get_curves(xL, y1_bot, xR, y2_bot);

        % Patch coordinates: from bar j to j+1
        patchX = [xTop; flipud(xBot)];
        patchY = [yTop; flipud(yBot)];

        ptch=patch(patchX, patchY, coloritos(m,:), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.4);
        uistack(ptch, 'bottom')
    end
end

%plot sig star
for i=1:length(projTargets)
    xvector=i;
    % yl=max([log_exp(i) log_obs(i)])+0.3;
    yl=cumulativeHeights(i)+0.03;
    % sigstars on top
    mysigstar(gca, xvector, yl, p_strength(i),'k', 'down', 'contextual');
end


hold off;

set(gca,'xticklabel',sample_names(projTargets),...
     'TickDir', 'out', 'box', 'off', 'FontSize',14)
% xlabel('Projection targets', 'FontSize', 11)
ylabel('Overall projection strength', 'FontSize', 14)
xticks(1:length(projTargets));
xlim([0.5 5.5])
% ylim([0 0.7])
% legend([h(2), h(3), h(4),h(1)],'upper', 'lower', 'anterior', 'posterior')
legend([h(2),h(1)], 'upper', 'lower')
% legend([h(1),h(2)], 'anterior', 'posterior')


%% alluvial plot of medial and lateral sSC projections to long-distance targets
% left_labels = {{'upper';'sSC'}, {'lower';'sSC'}};
% left_labels = {{'anterior';'sSC'}, {'posterior';'sSC'}};
% dataMatrix_fin=[dataMatrix_fin;dataMatrix];
left_labels = {{'upper';'sSC'}, {'lower';'sSC'},{'anterior';'sSC'}, {'posterior';'sSC'}};
right_labels = sample_names(projTargets);
figure; set(gcf, 'position', [425 225 623.2000 471.2000])
coloritos=[col(3,:);col(4,:);col(5,:);col(2,:)];
alluvialflow(dataMatrix, left_labels, right_labels, 'sSC projections to downstream targets',coloritos);
% alluvialflow(dataMatrix_fin, left_labels, right_labels,'',coloritos);
ax = gca;
ax.Position = [0.18 0.18 0.71 0.71];


%% alluvial plot of medial and lateral sSC projections to dSC
left_labels = {'upper sSC', 'lower sSC'};
right_labels = names;
targetsMatrix2_normtoregions=targetsMatrix2_all; %normalized to single sSC regions
targetsMatrix2_normtoall=targetsMatrix2_all_nonorm./sum(targetsMatrix2_all_nonorm,'all'); %normalized to the whole sSC
figure
alluvialflow(targetsMatrix2_normtoall, left_labels, right_labels, 'sSC spatial projections to dSC',coloritos);
%% plot bubblechart of proj strength and conditional prob for each SCs region - MED LAT

DATA=BigBnorm;
threshold=0;

% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<=cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
% cmp = colormap('prism');
% cmp = [0 0 0;cmp];
id_med=[3 7 11]; %medial SC
id_lat=[2 6 10]; %lateral SC
id_tot={id_med;id_lat};
names_medlat={'upper sSC','lower sSC'};

all_targets=[];
targetsMatrix_all={};
targetsMatrix2_all=[]; targetsMatrix2_all_nonorm=[];
for s = 1:size(id_tot,1) % go through medial and lateral injection sites
    id = id_tot{s};
    BCodes = find(ismember(maxes,id));
    % for each barcode find the # of areas
    areas=[];
    for b = 1:length(BCodes)
        curr = find(DATA(BCodes(b),:)>threshold);
        areas = [areas;length(curr)];
    end
    %     figure
    [~,ai] = sort(areas,'descend');
    % data_all = DATA(BCodes(ai),:);
    data_all = DATA(BCodes,:);
    data_all(data_all(:,:)<threshold)=0;

    %barcodes in each target area
    targets_bc=data_all(:, projTargets);
    targets_bc(targets_bc(:,:)<threshold)=0;
    zero_rows=all(targets_bc==0, 2);
    targets_bc=targets_bc(~zero_rows, :);

    %normalize values to sum of each row

    norm_vals=sum(targets_bc,2);
    targets_bc=targets_bc./norm_vals;
    all_targets=[all_targets ; targets_bc];

    %projection strength
    targetsMatrix2_comb = zeros(length(projTargets),1);
    for pr = 1:length(projTargets)
        targetsMatrix2_comb(pr,1) = sum(data_all(:,projTargets(pr)));
    end

    %normalize
    targetsMatrix2_nonorm=targetsMatrix2_comb; %keep this to save non-normalized matrix
    targetsMatrix2_comb=targetsMatrix2_comb./sum(targetsMatrix2_comb);


    targetsMatrix2_all=[targetsMatrix2_all; targetsMatrix2_comb'];
    targetsMatrix2_all_nonorm=[targetsMatrix2_all_nonorm; targetsMatrix2_nonorm'];


    targetsMatrix1_all = zeros(length(projTargets),length(projTargets));
    for pr = 1:length(projTargets)
        tmp_all=find(data_all(:,projTargets(pr)));
        val1_all= data_all(tmp_all,projTargets(pr));
        for pr2 = 1:length(projTargets)
            val2_all= data_all(tmp_all,projTargets(pr2));
            diff_all=val2_all-val1_all; diff_all=sum(diff_all(diff_all<0));
            proj_fraction_all=(sum(val1_all)+diff_all)/sum(val1_all);
            targetsMatrix1_all(pr2,pr) = proj_fraction_all;
        end
    end


    %or use the combined matrix
    targetsMatrix1_all(1:length(projTargets)+1:end)= NaN;
    targetsMatrix_all{s}=targetsMatrix1_all;

end


for s=1:length(targetsMatrix_all)
    figure
    %set nan values to transparent
    imAlpha=ones(size(targetsMatrix_all{s}));
    imAlpha(isnan(targetsMatrix_all{s}))=0;
    imagesc(targetsMatrix_all{s},'AlphaData',imAlpha);
    set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
        'color', [0.5 0.5 0.5], 'FontSize', 14)
    xlabel('Area A', 'FontSize',14)
    ylabel('Area B', 'FontSize',14)
    xticks(1:length(projTargets)); yticks(1:length(projTargets));
    if s==1
        % col1=cbrewer('seq', 'Oranges', 300);
        col1=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#FEF7EE','#F9DFC2','#F4B985', '#E79357', '#B86C44', '#764630'},300);
        col1=flipud(col1);
    else
        % col1=cbrewer('seq', 'Blues', 300);
        col1=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#F7FAFE','#D5DFEE','#9DBFD9', '#5A8EBF', '#2B599C', '#182D62'},300);
        col1=flipud(col1);
    end
    % col1=magma(300);
    col1(1,:)=[1 1 1];
    col1=col1(1:250,:);
    colormap(col1); colorbar; clim([0 0.8])
    title([names_medlat{s} ': subset of area A outputs projecting to area B'],...
        'FontSize', 13)


    figure(100+s); set(gcf, 'Position', [1 1 736.8000 507.2000]);

    bubbleSizes = zeros(length(projTargets),length(projTargets));
    for i=1:length(projTargets)
        for j=1:length(projTargets)
            bubbleSizes(i,j)=targetsMatrix_all{s}(i,j)*(targetsMatrix2_all(s,i)+targetsMatrix2_all(s,j)); %weighted conditional prob
            % bubbleSizes(i,j)=sqrt(targetsMatrix2_all(s,i)*targetsMatrix2_all(s,j)); %geometric mean
            % bubbleSizes(i,j)=2*(targetsMatrix2_all(s,i)*targetsMatrix2_all(s,j))/...
            %     (targetsMatrix2_all(s,i)+targetsMatrix2_all(s,j)); %harmonic mean
        end
    end

    % Generate coordinates for bubble positions based on array position
    [X, Y] = meshgrid(1:size(bubbleSizes, 2), 1:size(bubbleSizes, 1));
    X = X(:);
    Y = flip(Y(:));

    bubbleSizes2=bubbleSizes(:);
    bubbleSizes2(bubbleSizes2==0) = NaN;

    colormap(col1);
    bubblechart(X,Y,bubbleSizes2, targetsMatrix_all{s}(:), 'MarkerFaceAlpha',1)
    bubblesize([0.2 27])
    xlabel('Target area', 'FontSize', 15)
    ylabel('Source region', 'FontSize', 15)
    title([names_medlat{s} ': weighted contitional probability'], 'FontSize', 15)


    axis equal
    set(gca, 'Box', 'on', ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'TickDir', 'in', 'TickLength', [0 0], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', ...
        'Color', 'none','XTick',1:length(projTargets), 'YTick',1:length(projTargets),...
        'XTickLabel',sample_names(projTargets),...
        'YTickLabel',flip(sample_names(projTargets)), 'FontSize', 14)
    xlim([0.3 size(targetsMatrix_all{s}, 2) + 0.7]); % Add some padding to the x-axis
    ylim([0.3 size(targetsMatrix_all{s}, 1) + 0.7]);
    xtickangle(45)

    %
    % blgd = bubblelegend('Projection strength','Style','vertical',...
    %     'BubbleSizeOrder','descending',...
    %     'box','off',...
    %     'Location','northeastoutside',...
    %     'NumBubbles',3,...
    %     'FontSize', 8);

    % blgd.LimitLabels={'0.003', '0.18','0.44'};
    mx=max(targetsMatrix_all{s}, [], 'all');
    % clim([0 ceil(mx*10)/10
    clim([0 0.8])
    % cb=colorbar('Position',[0.790988056460368,0.203958908811558,0.024972855591744,0.302350239453425]);
    cb=colorbar('eastoutside');
    ylabel(cb, 'conditional prob');
    cb.FontSize=14;


    figure(s+2)
    %set nan values to transparent
    imAlpha=ones(size(targetsMatrix_all{s}));
    imAlpha(isnan(targetsMatrix_all{s}))=0;
    imagesc(bubbleSizes,'AlphaData',imAlpha);
    set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
        'color', [0.5 0.5 0.5], 'FontSize', 14)
    xlabel('Area A', 'FontSize',14)
    ylabel('Area B', 'FontSize',14)
    xticks(1:length(projTargets)); yticks(1:length(projTargets));
    if s==1
        % col1=cbrewer('seq', 'Oranges', 300);
        col1=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#FEF7EE','#F9DFC2','#F4B985', '#E79357', '#B86C44', '#764630'},300);
        col1=flipud(col1);
    else
        % col1=cbrewer('seq', 'Blues', 300);
        col1=customcolormap([0 0.2 0.4 0.6 0.8 1], {'#F7FAFE','#D5DFEE','#9DBFD9', '#5A8EBF', '#2B599C', '#182D62'},300);
        col1=flipud(col1);
    end
    % col1=magma(300);
    col1(1,:)=[1 1 1];
    col1=col1(1:250,:);
    colormap(col1); colorbar; clim([0 0.8])
    title([names_medlat{s} ' weighted conditional probabilities'],...
        'FontSize', 13)
    %     title([SampleName{3}(sSC(s)) ': % barcodes found in >1 long-range targets'])

end


coloritos=[col(3,:);col(4,:)];
figure(10)
h=bar(1:length(projTargets),targetsMatrix2_all, 'FaceAlpha', 1, 'EdgeColor', 'none');
for k = 1:length(h)
    h(k).FaceColor = 'flat';
    h(k).CData = repmat(coloritos(k, :), size(targetsMatrix2_all, 2), 1);
end
set(gca,'xticklabel',sample_names(projTargets),...
     'TickDir', 'out', 'box', 'off', 'FontSize',14)
% xlabel('Projection targets', 'FontSize', 11)
ylabel('Projection strength', 'FontSize', 14)
xticks(1:length(projTargets)); ylim([0 1])
legend('upper', 'lower')


%% plot projection strength for each SCs region and mice - sSCam sSCm sSCal sSCl sSCc
DATA1 = B1norm;
DATA2 = B2norm;
DATA3 = B3norm;
DATA=BigBnorm;

% projTargets=40:44;
projTargets=[40 42 41 43 44];

% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<=cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;

% id_tot={[3 7]; 11; [2 6]; 10; 14};
id_tot={[3 7]; [2 6]; 11; 10; 14};

names_medlat={'sSCam', 'sSCm', 'sSCal' ,'sSCl', 'sSCc'};
names=sample_names(projTargets);
threshold=0;


scnt = 0; Numbers={};
Numbers{1} = zeros(length(sSC),2);Numbers{2} = zeros(length(sSC),2);
Numbers{3} = zeros(length(sSC),2); allKeep1 =[]; allKeep2=[]; allKeep3=[];
all_targets1=[]; all_targets2=[]; all_targets3=[]; all_targets=[];
targetsMatrix2_ind={}; targetsMatrix2_all=[]; targetsMatrix2_all_nonorm=[];
neuronnum=0;
for s = 1:size(id_tot,1) % go through medial and lateral injection sites
    scnt = scnt +1;
    id = id_tot{s};
    BCodes1 = find(ismember(m1,id));
    BCodes2 = find(ismember(m2,id));
    BCodes3 = find(ismember(m3,id));
    BCodes = find(ismember(maxes,id));
    Numbers{1}(scnt,1) = length(BCodes1);
    Numbers{2}(scnt,1) = length(BCodes2);
    Numbers{3}(scnt,1) = length(BCodes3);
    % for each barcode find the # of areas
    areas1 = []; areas2 = []; areas3 = []; areas=[];
    for b = 1:length(BCodes1)
        curr = find(DATA1(BCodes1(b),:)>threshold);
        areas1 = [areas1;length(curr)];
    end
    for b = 1:length(BCodes2)
        curr = find(DATA2(BCodes2(b),:)>threshold);
        areas2 = [areas2;length(curr)];
    end
    for b = 1:length(BCodes3)
        curr = find(DATA3(BCodes3(b),:)>threshold);
        areas3 = [areas3;length(curr)];
    end
    for b = 1:length(BCodes)
        curr = find(DATA(BCodes(b),:)>threshold);
        areas = [areas;length(curr)];
    end
    %     figure
    [~,a2] = sort(areas1,'descend');
    data1 = DATA1(BCodes1(a2),:);
    [~,a2] = sort(areas2,'descend');
    data2 = DATA2(BCodes2(a2),:);
    [~,a2] = sort(areas3,'descend');
    data3 = DATA3(BCodes3(a2),:);
    [~,ai] = sort(areas,'descend');
    % data_all = DATA(BCodes(ai),:);
    data_all = DATA(BCodes,:);
    data_all(data_all(:,:)<threshold)=0;
    data_proj=data_all(:,projTargets);
    data_proj(sum(data_proj, 2) == 0, :) = [];
    neuronnum=neuronnum+size(data_proj,1);
    keepData1 = data1; keepData2 = data2; keepData3 = data3;

    %projection strength

    targetsMatrix2_m1 = zeros(length(projTargets),1);
    targetsMatrix2_m2 = zeros(length(projTargets),1);
    targetsMatrix2_m3 = zeros(length(projTargets),1);
    targetsMatrix2_comb = zeros(length(projTargets),1);
    for pr = 1:length(projTargets)
        targetsMatrix2_m1(pr,1) = sum(keepData1(:,projTargets(pr)));
        targetsMatrix2_m2(pr,1) = sum(keepData2(:,projTargets(pr)));
        targetsMatrix2_m3(pr,1) = sum(keepData3(:,projTargets(pr)));
        targetsMatrix2_comb(pr,1) = sum(data_all(:,projTargets(pr)));
    end

    %normalize
    targetsMatrix2_m1=targetsMatrix2_m1./sum(targetsMatrix2_m1);
    targetsMatrix2_m2=targetsMatrix2_m2./sum(targetsMatrix2_m2);
    targetsMatrix2_m3=targetsMatrix2_m3./sum(targetsMatrix2_m3);
    targetsMatrix2_nonorm=targetsMatrix2_comb; %keep this to save non-normalized matrix
    targetsMatrix2_comb=targetsMatrix2_comb./sum(targetsMatrix2_comb);

    targetsMatrix2_all=[targetsMatrix2_all; targetsMatrix2_comb'];
    targetsMatrix2_all_nonorm=[targetsMatrix2_all_nonorm; targetsMatrix2_nonorm'];

    targetsMatrix1_m1 = zeros(length(projTargets),length(projTargets));
    targetsMatrix1_m2 = zeros(length(projTargets),length(projTargets));
    targetsMatrix1_m3 = zeros(length(projTargets),length(projTargets));
    targetsMatrix1_all = zeros(length(projTargets),length(projTargets));
    for pr = 1:length(projTargets)
        tmp1 = find(keepData1(:,projTargets(pr)));
        tmp2 = find(keepData2(:,projTargets(pr)));
        tmp3 = find(keepData3(:,projTargets(pr)));
        tmp_all=find(data_all(:,projTargets(pr)));
        val1_m1 = keepData1(tmp1,projTargets(pr));
        val1_m2 = keepData2(tmp2,projTargets(pr));
        val1_m3 = keepData3(tmp3,projTargets(pr));
        val1_all= data_all(tmp_all,projTargets(pr));
        for pr2 = 1:length(projTargets)
            val2_m1 = keepData1(tmp1,projTargets(pr2));
            val2_m2 = keepData2(tmp2,projTargets(pr2));
            val2_m3 = keepData3(tmp3,projTargets(pr2));
            val2_all= data_all(tmp_all,projTargets(pr2));
            diff1 = val2_m1-val1_m1; diff1=sum(diff1(diff1<0));
            diff2 = val2_m2-val1_m2; diff2=sum(diff2(diff2<0));
            diff3 = val2_m3-val1_m3; diff3=sum(diff3(diff3<0));
            diff_all=val2_all-val1_all; diff_all=sum(diff_all(diff_all<0));
            proj_fraction1=(sum(val1_m1)+diff1)/sum(val1_m1);
            proj_fraction2=(sum(val1_m2)+diff2)/sum(val1_m2);
            proj_fraction3=(sum(val1_m3)+diff3)/sum(val1_m3);
            proj_fraction_all=(sum(val1_all)+diff_all)/sum(val1_all);
            targetsMatrix1_m1(pr2,pr) = proj_fraction1;
            targetsMatrix1_m2(pr2,pr) = proj_fraction2;
            targetsMatrix1_m3(pr2,pr) = proj_fraction3;
            targetsMatrix1_all(pr2,pr) = proj_fraction_all;
        end
    end

    %calculate mean
    % targetsMatrix2_mean=mean(cat(3, targetsMatrix2_m1,targetsMatrix2_m2,targetsMatrix2_m3),3,'omitnan');

    % %calculate mean from single matrices
    % targetsMatrix1_mean=mean(cat(3, targetsMatrix1_m1,targetsMatrix1_m2,targetsMatrix1_m3),3,'omitnan');
    % %set diagonal values to NaN
    % targetsMatrix1_mean(1:11:end)= NaN;
    % targetsMatrix1_m1(1:11:end)= NaN;
    % targetsMatrix1_m2(1:11:end)= NaN;
    % targetsMatrix1_m3(1:11:end)= NaN;

    %or use the combined matrix
    targetsMatrix1_all(1:length(projTargets)+1:end)= NaN;
    targetsMatrix1_m1(1:length(projTargets)+1:end)= NaN;
    targetsMatrix1_m2(1:length(projTargets)+1:end)= NaN;
    targetsMatrix1_m3(1:length(projTargets)+1:end)= NaN;


    figure
    %set nan values to transparent
    imAlpha=ones(size(targetsMatrix1_all));
    imAlpha(isnan(targetsMatrix1_all))=0;
    imagesc(targetsMatrix1_all,'AlphaData',imAlpha);
    set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
    'color', [0.5 0.5 0.5], 'FontSize', 12)
    xlabel('Area A', 'FontSize',12)
    ylabel('Area B', 'FontSize',12)
    xticks(1:length(projTargets)); yticks(1:length(projTargets));
    col1=cbrewer('seq', 'Purples', 300);
    % col1=magma(300);
    colormap(col1); colorbar; clim([0 0.8])
    title([names_medlat{s} ': subset of area A outputs projecting to area B'],...
        'FontSize', 13)
%     title([SampleName{3}(sSC(s)) ': % barcodes found in >1 long-range targets'])
% %
    Numbers{1}(scnt,2) = sum(targetsMatrix2_m1);
    Numbers{2}(scnt,2) = sum(targetsMatrix2_m2);
    Numbers{3}(scnt,2) = sum(targetsMatrix2_m3);
    targetsMatrix2_ind{s}=[targetsMatrix2_m1 targetsMatrix2_m2 targetsMatrix2_m3];

    figure(100); set(gcf, 'Position', [97.8000 337.8000 950.2000 324.2000]);
    subplot(2,3,scnt)
    % bar(targetsMatrix2_mean); hold on
    % scatter(1:10, targetsMatrix2_all, 10, 'filled')
    plot(1:length(projTargets),targetsMatrix2_m1,'-o','Color', '#6A2C70', LineWidth=2); hold on
    plot(1:length(projTargets),targetsMatrix2_m2,'-o','Color', '#B83B5E', LineWidth=2); hold on
    plot(1:length(projTargets),targetsMatrix2_m3,'-o','Color', '#FFAF45', LineWidth=2)
    xlim([0.5 length(projTargets)+0.5]); xticks(1:length(projTargets)); ylim([0 0.9])
    set(gca,'xticklabel',sample_names(projTargets))
    xtickangle(90)
    ylabel('projection strength')
    title(names_medlat{s})
    legend('m1', 'm2', 'm3');

%     ylabel('# of barcodes')
%     title([SampleName{3}{sSC(s)},' | total: ',int2str(Numbers(scnt,1)),', in targets: ',int2str(Numbers(scnt,2)) ])
end


% coloritos=[col(3,:);col(4,:)];
figure(10)
h=bar(1:length(projTargets),targetsMatrix2_all', 'FaceAlpha', 0.7, 'EdgeColor', 'none');
for k = 1:length(h)
    h(k).FaceColor = 'flat';
    % h(k).CData = repmat(coloritos(k, :), size(targetsMatrix2_all, 1), 1);
end
set(gca,'xticklabel',sample_names(projTargets),...
     'TickDir', 'out', 'box', 'off', 'FontSize',12)
% xlabel('Projection targets', 'FontSize', 11)
ylabel('Projection strength', 'FontSize', 12)
xticks(1:length(projTargets)); ylim([0 1])
% legend('sSCam', 'sSCm', 'sSCal', 'sSCl','sSCc')
legend('sSCam', 'sSCal','sSCm', 'sSCl','sSCc')

%% alluvial plot of sSC projections to dSC
% left_labels = {'sSCam', 'sSCm', 'sSCal' ,'sSCl', 'sSCc'};
left_labels = {'sSCam', 'sSCal' ,'sSCm', 'sSCl', 'sSCc'};
right_labels = names;
colortable=[0 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.9290 0.6940 0.1250;
            0.4940 0.1840 0.5560;
            0.4660 0.6740 0.1880;
            0.3010 0.7450 0.9330;
            0.3010 0.7450 0.9330];
targetsMatrix2_all(isnan(targetsMatrix2_all))=0;
targetsMatrix2_normtoregions=targetsMatrix2_all; %normalized to single sSC regions
targetsMatrix2_normtoall=targetsMatrix2_all_nonorm./sum(targetsMatrix2_all_nonorm,'all'); %normalized to the whole sSC
%normalize to medial-lateral and caudal sSC
norm_medlatcaud=sum(targetsMatrix2_all_nonorm,2);
norm_medlatcaud(1:2)=sum(norm_medlatcaud(1:2));
norm_medlatcaud(3:4)=sum(norm_medlatcaud(3:4));
targetsMatrix2_normtomedlatcaud=targetsMatrix2_all_nonorm./norm_medlatcaud;
alluvialflow(targetsMatrix2_normtoregions, left_labels, right_labels, ['sSC spatial projections to dSC (n=', num2str(neuronnum),')'],colortable);

%% plot bubblechart of proj strength and conditional prob for each SCs region - sSC am al m l c

DATA=BigBnorm;
threshold=0;

% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<=cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
% cmp = colormap('prism');
% cmp = [0 0 0;cmp];
id_tot={[3 7]; 11; [2 6]; 10; 14};

names_medlat={'sSCam', 'sSCm', 'sSCal' ,'sSCl', 'sSCc'};
names=sample_names(projTargets);
threshold=0;


all_targets=[];
targetsMatrix_all={};
targetsMatrix2_all=[]; targetsMatrix2_all_nonorm=[];neuronnum=0; n_neurons=[];
for s = 1:size(id_tot,1) % go through medial and lateral injection sites
    id = id_tot{s};
    BCodes = find(ismember(maxes,id));
    % for each barcode find the # of areas
    areas=[];
    for b = 1:length(BCodes)
        curr = find(DATA(BCodes(b),:)>threshold);
        areas = [areas;length(curr)];
    end
    %     figure
    [~,ai] = sort(areas,'descend');
    % data_all = DATA(BCodes(ai),:);
    data_all = DATA(BCodes,:);
    data_all(data_all(:,:)<threshold)=0;

    %barcodes in each target area
    targets_bc=data_all(:, projTargets);
    targets_bc(targets_bc(:,:)<threshold)=0;
    zero_rows=all(targets_bc==0, 2);
    targets_bc=targets_bc(~zero_rows, :);
    temp_num=targets_bc(:,1);
    temp_num(all(temp_num == 0, 2), :)=[];
    neuronnum=neuronnum+size(temp_num,1);

    %normalize values to sum of each row

    norm_vals=sum(targets_bc,2);
    targets_bc=targets_bc./norm_vals;
    all_targets=[all_targets ; targets_bc];
    n_neurons(s)=size(targets_bc,1);

    %projection strength
    targetsMatrix2_comb = zeros(length(projTargets),1);
    for pr = 1:length(projTargets)
        targetsMatrix2_comb(pr,1) = sum(data_all(:,projTargets(pr)));
    end

    %normalize
    targetsMatrix2_nonorm=targetsMatrix2_comb; %keep this to save non-normalized matrix
    targetsMatrix2_comb=targetsMatrix2_comb./sum(targetsMatrix2_comb);


    targetsMatrix2_all=[targetsMatrix2_all; targetsMatrix2_comb'];
    targetsMatrix2_all_nonorm=[targetsMatrix2_all_nonorm; targetsMatrix2_nonorm'];


    targetsMatrix1_all = zeros(length(projTargets),length(projTargets));
    for pr = 1:length(projTargets)
        tmp_all=find(data_all(:,projTargets(pr)));
        val1_all= data_all(tmp_all,projTargets(pr));
        for pr2 = 1:length(projTargets)
            val2_all= data_all(tmp_all,projTargets(pr2));
            diff_all=val2_all-val1_all; diff_all=sum(diff_all(diff_all<0));
            proj_fraction_all=(sum(val1_all)+diff_all)/sum(val1_all);
            targetsMatrix1_all(pr2,pr) = proj_fraction_all;
        end
    end


    %or use the combined matrix
    targetsMatrix1_all(1:length(projTargets)+1:end)= NaN;
    targetsMatrix_all{s}=targetsMatrix1_all;

end


for s=1:length(targetsMatrix_all)
    figure
    %set nan values to transparent
    imAlpha=ones(size(targetsMatrix_all{s}));
    imAlpha(isnan(targetsMatrix_all{s}))=0;
    imagesc(targetsMatrix_all{s},'AlphaData',imAlpha);
    set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
        'color', [0.5 0.5 0.5], 'FontSize', 14)
    xlabel('Area A', 'FontSize',14)
    ylabel('Area B', 'FontSize',14)
    xticks(1:length(projTargets)); yticks(1:length(projTargets));
    col1=cbrewer('seq', 'Purples', 300);
    col1(1,:)=[1 1 1];
    col1=col1(1:250,:);
    colormap(col1); colorbar; clim([0 0.8])
    title([names_medlat{s} ': subset of area A outputs projecting to area B'],...
        'FontSize', 13)


    figure(100+s); set(gcf, 'Position', [1 1 736.8000 507.2000]);

    bubbleSizes = zeros(length(projTargets),length(projTargets));
    for i=1:length(projTargets)
        for j=1:length(projTargets)
            bubbleSizes(i,j)=targetsMatrix_all{s}(i,j)*(targetsMatrix2_all(s,i)+targetsMatrix2_all(s,j)); %weighted conditional prob
            % bubbleSizes(i,j)=sqrt(targetsMatrix2_all(s,i)*targetsMatrix2_all(s,j)); %geometric mean
            % bubbleSizes(i,j)=2*(targetsMatrix2_all(s,i)*targetsMatrix2_all(s,j))/...
            %     (targetsMatrix2_all(s,i)+targetsMatrix2_all(s,j)); %harmonic mean
        end
    end

    % Generate coordinates for bubble positions based on array position
    [X, Y] = meshgrid(1:size(bubbleSizes, 2), 1:size(bubbleSizes, 1));
    X = X(:);
    Y = flip(Y(:));

    bubbleSizes2=bubbleSizes(:);
    bubbleSizes2(bubbleSizes2==0) = NaN;

    colormap(col1);
    bubblechart(X,Y,bubbleSizes2, targetsMatrix_all{s}(:), 'MarkerFaceAlpha',1)
    bubblesize([0.2 27])
    xlabel('Target area', 'FontSize', 15)
    ylabel('Source region', 'FontSize', 15)
    title([names_medlat{s} ': weighted cond prob (', 'n=', num2str(n_neurons(s)),')'], 'FontSize', 15)


    axis equal
    set(gca, 'Box', 'on', ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'TickDir', 'in', 'TickLength', [0 0], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', ...
        'Color', 'none','XTick',1:length(projTargets), 'YTick',1:length(projTargets),...
        'XTickLabel',sample_names(projTargets),...
        'YTickLabel',flip(sample_names(projTargets)), 'FontSize', 14)
    xlim([0.3 size(targetsMatrix_all{s}, 2) + 0.7]); % Add some padding to the x-axis
    ylim([0.3 size(targetsMatrix_all{s}, 1) + 0.7]);
    xtickangle(45)

    %
    % blgd = bubblelegend('Projection strength','Style','vertical',...
    %     'BubbleSizeOrder','descending',...
    %     'box','off',...
    %     'Location','northeastoutside',...
    %     'NumBubbles',3,...
    %     'FontSize', 8);

    % blgd.LimitLabels={'0.003', '0.18','0.44'};
    mx=max(targetsMatrix_all{s}, [], 'all');
    % clim([0 ceil(mx*10)/10
    clim([0 0.8])
    % cb=colorbar('Position',[0.790988056460368,0.203958908811558,0.024972855591744,0.302350239453425]);
    cb=colorbar('eastoutside');
    ylabel(cb, 'conditional prob');
    cb.FontSize=14;


    figure(s+2)
    %set nan values to transparent
    imAlpha=ones(size(targetsMatrix_all{s}));
    imAlpha(isnan(targetsMatrix_all{s}))=0;
    imagesc(bubbleSizes,'AlphaData',imAlpha);
    set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
        'color', [0.5 0.5 0.5], 'FontSize', 14)
    xlabel('Area A', 'FontSize',14)
    ylabel('Area B', 'FontSize',14)
    xticks(1:length(projTargets)); yticks(1:length(projTargets));
    col1=cbrewer('seq', 'Purples', 300);
    col1=col1(1:250,:);
    col1(1,:)=[1 1 1];
    colormap(col1); colorbar; clim([0 0.8])
    title([names_medlat{s} ' weighted cond prob'],...
        'FontSize', 13)
    %     title([SampleName{3}(sSC(s)) ': % barcodes found in >1 long-range targets'])

end


figure(10)
h=bar(1:length(projTargets),targetsMatrix2_all', 'FaceAlpha', 0.7, 'EdgeColor', 'none');
for k = 1:length(h)
    h(k).FaceColor = 'flat';
    % h(k).CData = repmat(coloritos(k, :), size(targetsMatrix2_all, 1), 1);
end
set(gca,'xticklabel',sample_names(projTargets),...
     'TickDir', 'out', 'box', 'off', 'FontSize',12)
% xlabel('Projection targets', 'FontSize', 11)
ylabel('Projection strength', 'FontSize', 12)
xticks(1:length(projTargets)); ylim([0 1])
legend('sSCam', 'sSCm', 'sSCal', 'sSCl','sSCc')
%% plot projection strength for each target area and mice
DATA1 = B1norm;
DATA2 = B2norm;
DATA3 = B3norm;
DATA = BigBnorm;

BCodes = find(ismember(maxes,sSC));
data=BigB(BCodes,:);
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;

threshold=0;

targetsMatrix2_m1 = zeros(length(projTargets),length(sSC));
targetsMatrix2_m2 = zeros(length(projTargets),length(sSC));
targetsMatrix2_m3 = zeros(length(projTargets),length(sSC));
hm_matrix = zeros(length(projTargets),length(sSC));
allKeep1=[]; allKeep2=[]; allKeep3=[];

for s = orderOfSC % go through injection sites
    id = sSC(s);
    BCodes1 = find(m1==id);
    BCodes2 = find(m2==id);
    BCodes3 = find(m3==id);
    BCodes = find(maxes==id);
    % for each barcode find the # of areas
    areas1 = []; areas2 = []; areas3 = []; areas=[];
    for b = 1:length(BCodes1)
        curr = find(DATA1(BCodes1(b),1:34)>threshold);
        areas1 = [areas1;length(curr)];
    end
    for b = 1:length(BCodes2)
        curr = find(DATA2(BCodes2(b),1:34)>threshold);
        areas2 = [areas2;length(curr)];
    end
    for b = 1:length(BCodes3)
        curr = find(DATA3(BCodes3(b),1:34)>threshold);
        areas3 = [areas3;length(curr)];
    end
    for b = 1:length(BCodes)
        curr = find(DATA(BCodes(b),1:34)>threshold);
        areas = [areas;length(curr)];
    end
    %     figure
    [~,a2] = sort(areas1,'descend');
    data1 = DATA1(BCodes1(a2),1:34);
    [~,a2] = sort(areas2,'descend');
    data2 = DATA2(BCodes2(a2),1:34);
    [~,a2] = sort(areas3,'descend');
    data3 = DATA3(BCodes3(a2),1:34);
    [~,ai] = sort(areas,'descend');
    data_all = DATA(BCodes(ai),1:34);
    keepData1 = data1; keepData2 = data2; keepData3 = data3;
    allKeep1=[allKeep1;keepData1];
    allKeep2=[allKeep2;keepData2];
    allKeep3=[allKeep3;keepData3];

    %projection strength

    for pr = 1:length(projTargets)
        targetsMatrix2_m1(pr,s) = sum(keepData1(:,projTargets(pr)));
        targetsMatrix2_m2(pr,s) = sum(keepData2(:,projTargets(pr)));
        targetsMatrix2_m3(pr,s) = sum(keepData3(:,projTargets(pr)));
        hm_matrix(pr,s) = sum(data_all(:,projTargets(pr)));
    end

end

%normalize
targetsMatrix2_m1=targetsMatrix2_m1./sum(targetsMatrix2_m1,2);
targetsMatrix2_m2=targetsMatrix2_m2./sum(targetsMatrix2_m2,2);
targetsMatrix2_m3=targetsMatrix2_m3./sum(targetsMatrix2_m3,2);
hm_matrix=hm_matrix./sum(hm_matrix,2);

%separate ipsi and contra and med and lat

ix=[2 4 6 1 3 5 7];
targetsMatrix2_m1=targetsMatrix2_m1(:, ix);
targetsMatrix2_m2=targetsMatrix2_m2(:, ix);
targetsMatrix2_m3=targetsMatrix2_m3(:, ix);
hm_matrix=hm_matrix(:, ix);

target_names=sample_names(projTargets);
source_names=SampleName{3}(sSC);
%
figure(100)
for i=1:length(projTargets)
    subplot(3,4,i)
    % bar(targetsMatrix2_mean); hold on
    % scatter(1:10, targetsMatrix2_all, 10, 'filled')
    plot(1:7,targetsMatrix2_m1(i,:),'-o','Color', '#6A2C70', LineWidth=2); hold on
    plot(1:7,targetsMatrix2_m2(i,:),'-o','Color', '#B83B5E', LineWidth=2); hold on
    plot(1:7,targetsMatrix2_m3(i,:),'-o','Color', '#FFAF45', LineWidth=2)
    xlim([0.5 7.5]); xticks(1:7); ylim([0 0.7]);
    set(gca,'xticklabel',source_names(ix))
    xtickangle(90)
    ylabel('projection strength')
    title(target_names(i))
    % legend('m1', 'm2', 'm3', 'FontSize', 7);
end


%make heatmap of combined matrices with target areas x SCs region
figure(2); set(gcf, 'Position', [811.4000 242.6000 415.2000 420]);
imAlpha=ones(size(hm_matrix));
imAlpha(isnan(hm_matrix))=0;
imagesc(hm_matrix,'AlphaData',imAlpha);
set(gca,'xticklabel',source_names(ix), 'yticklabel',target_names,...
    'color', [0.5 0.5 0.5], 'FontSize', 8)
xticks(1:7); xtickangle(45)
xlabel('Source region', 'FontSize',10)
ylabel('Target area', 'FontSize',10)
col1=magma(300);
colormap(col1); colorbar; clim([0 0.55])
title('Projection strength from source regions')


%% all co-labelled projections together: projection fraction - mean across mice
% tmp = allKeep<threshold;
% allKeep(tmp) = 0;
allKeep=[allKeep1;allKeep2;allKeep3];

%projection strength
targetsMatrix2_m1 = zeros(length(projTargets),1);
targetsMatrix2_m2 = zeros(length(projTargets),1);
targetsMatrix2_m3 = zeros(length(projTargets),1);
targetsMatrix2_all=zeros(length(projTargets),1);
for pr = 1:length(projTargets)
    targetsMatrix2_m1(pr,1) = sum(allKeep1(:,projTargets(pr)));
    targetsMatrix2_m2(pr,1) = sum(allKeep2(:,projTargets(pr)));
    targetsMatrix2_m3(pr,1) = sum(allKeep3(:,projTargets(pr)));
    targetsMatrix2_all(pr,1)= sum(allKeep(:,projTargets(pr)));
end
targetsMatrix2_m1=targetsMatrix2_m1./sum(targetsMatrix2_m1);
targetsMatrix2_m2=targetsMatrix2_m2./sum(targetsMatrix2_m2);
targetsMatrix2_m3=targetsMatrix2_m3./sum(targetsMatrix2_m3);
targetsMatrix2_all=targetsMatrix2_all./sum(targetsMatrix2_all);

% targetsMatrix2_m1(isnan(targetsMatrix2_m1))=0;
% targetsMatrix2_m2(isnan(targetsMatrix2_m2))=0;
% targetsMatrix2_m3(isnan(targetsMatrix2_m3))=0;


figure(1)
set(gcf, 'Position', [353 247.4000 840.8000 283.2000])
% bar(targetsMatrix2_mean); hold on
% scatter(1:10, targetsMatrix2_all, 20, 'filled')
subplot(1,2,1)
plot(1:5,targetsMatrix2_m1,'-o','Color', '#6A2C70', LineWidth=2); hold on
plot(1:5,targetsMatrix2_m2,'-o','Color', '#B83B5E', LineWidth=2); hold on
plot(1:5,targetsMatrix2_m3,'-o','Color', '#FFAF45', LineWidth=2)
xlim([0.5 5.5]); xticks(1:10);
ylim([0 0.7]);
set(gca,'xticklabel',sample_names(projTargets))
xtickangle(90)
ylabel('projection strength')
title('ipsilateral SCs')
legend('m1', 'm2', 'm3', 'FontSize', 12);

% %ipsi and contra
% subplot(1,2,2)
% iy=[1 3 5 7 9 2 4 6 8 10];
% ipsi1=sum(targetsMatrix2_m1(iy(1:5))); contra1=sum(targetsMatrix2_m1(iy(6:end)));
% ipsi2=sum(targetsMatrix2_m2(iy(1:5))); contra2=sum(targetsMatrix2_m2(iy(6:end)));
% ipsi3=sum(targetsMatrix2_m3(iy(1:5))); contra3=sum(targetsMatrix2_m3(iy(6:end)));
% plot(1:2, [ipsi1 contra1], '-o', 'Color', '#6A2C70', LineWidth=2); hold on
% plot(1:2, [ipsi2 contra2], '-o', 'Color', '#B83B5E', LineWidth=2); hold on
% plot(1:2, [ipsi3 contra3], '-o', 'Color', '#FFAF45', LineWidth=2);
% xlim([0.5 2.5]); xticks(1:2);
% set(gca,'xticklabel', ["ipsi", "contra"])


%fraction of neurons innervating a single area out of all neurons that
%innervate that area
nonzero_counts1=zeros(1, size(all_targets1, 2)); %number of neurons projecting to each area
nonzero_counts2=zeros(1, size(all_targets2, 2));
nonzero_counts3=zeros(1, size(all_targets3, 2));
for i = 1:size(all_targets1, 2)
    nonzero_counts1(i) = nnz(all_targets1(:, i));
    nonzero_counts2(i) = nnz(all_targets2(:, i));
    nonzero_counts3(i) = nnz(all_targets3(:, i));
end
ones_counts1 = zeros(1, size(all_targets1, 2)); %number of neurons that project only one particular area
ones_counts2 = zeros(1, size(all_targets2, 2));
ones_counts3 = zeros(1, size(all_targets3, 2));
for i = 1:size(all_targets1, 2)
    ones_counts1(i)=sum(all_targets1(:, i) == 1);
    ones_counts2(i)=sum(all_targets2(:, i) == 1);
    ones_counts3(i)=sum(all_targets3(:, i) == 1);
end
frac_ones1=ones_counts1./nonzero_counts1; frac_ones1(isnan(frac_ones1))=0;
frac_ones2=ones_counts2./nonzero_counts2; frac_ones2(isnan(frac_ones2))=0;
frac_ones3=ones_counts3./nonzero_counts3; frac_ones3(isnan(frac_ones3))=0;

% frac_ones_mean=mean(cat(3, frac_ones1',frac_ones2',frac_ones3'),3,'omitnan');
cat_frac=cat(3, frac_ones1',frac_ones2',frac_ones3');
cat_frac(isnan(cat_frac))=0;
frac_ones_mean=mean(cat_frac,3);
frac_ones_all=[frac_ones1' frac_ones2' frac_ones3'];

figure(2)
% bar(frac_ones_mean); hold on
% scatter(1:10, frac_ones_all, 20, 'filled')
plot(1:5,frac_ones1,'-o','Color', '#6A2C70', LineWidth=2); hold on
plot(1:5,frac_ones2,'-o','Color', '#B83B5E', LineWidth=2); hold on
plot(1:5,frac_ones3,'-o','Color', '#FFAF45', LineWidth=2)
xlim([0.5 5.5]); xticks(1:5);
set(gca,'xticklabel',sample_names(projTargets))
xtickangle(90)
ylabel('fraction of dedicated input')
title('ipsilateral SCs')
legend('m1', 'm2', 'm3', 'FontSize', 12);


%piechart and barplot of SCs fractions
figure(5)
subplot(1,2,1)
SCs_fractions1=sum(allKeep1(:,sSC))/sum(allKeep1(:,sSC), 'all');
SCs_fractions2=sum(allKeep2(:,sSC))/sum(allKeep2(:,sSC), 'all');
SCs_fractions3=sum(allKeep3(:,sSC))/sum(allKeep3(:,sSC), 'all');
SCs_fractions=mean([SCs_fractions1; SCs_fractions2; SCs_fractions3], 1);

[SCs_fractions,id]=sort(SCs_fractions, 'descend');
labels=string(SampleName{3}(sSC));
labels=labels(id');
p=piechart(SCs_fractions, labels);
p.Direction = "counterclockwise";
p.EdgeColor = "white"; p.LineWidth = 1;
colors=cbrewer('qual','Set3',7);
colormap(colors);

subplot(1,2,2)
% bar(SCs_fractions); hold on
SCs_fractions_all=[SCs_fractions1' SCs_fractions2' SCs_fractions3'];
% scatter(1:7, SCs_fractions_all, 20, 'filled'); hold on
plot(1:7,SCs_fractions1,'-o','Color', '#6A2C70', LineWidth=2); hold on
plot(1:7,SCs_fractions2,'-o','Color', '#B83B5E', LineWidth=2); hold on
plot(1:7,SCs_fractions3,'-o','Color', '#FFAF45', LineWidth=2)
xlim([0.5 7.5]); xticks(1:7);
set(gca,'xticklabel',string(SampleName{3}(sSC)), 'FontSize', 8);
xtickangle(45)
ylabel('soure region ratio', 'FontSize',8)
title('ipsilateral SCs fractions')
legend('m1', 'm2', 'm3', 'FontSize', 12);
set(gcf, 'Position', [227.4000 248.2000 956.0000 300])


targetsMatrix1_m1 = zeros(length(projTargets),length(projTargets));
targetsMatrix1_m2 = zeros(length(projTargets),length(projTargets));
targetsMatrix1_m3 = zeros(length(projTargets),length(projTargets));
targetsMatrix1_all = zeros(length(projTargets),length(projTargets));
for pr = 1:length(projTargets)
    tmp1 = find(allKeep1(:,projTargets(pr)));
    tmp2 = find(allKeep2(:,projTargets(pr)));
    tmp3 = find(allKeep3(:,projTargets(pr)));
    tmp_all=find(allKeep(:,projTargets(pr)));
    val1_m1 = allKeep1(tmp1,projTargets(pr));
    val1_m2 = allKeep2(tmp2,projTargets(pr));
    val1_m3 = allKeep3(tmp3,projTargets(pr));
    val1_all= allKeep(tmp_all,projTargets(pr));
    for pr2 = 1:length(projTargets)
        val2_m1 = allKeep1(tmp1,projTargets(pr2));
        val2_m2 = allKeep2(tmp2,projTargets(pr2));
        val2_m3 = allKeep3(tmp3,projTargets(pr2));
        val2_all= allKeep(tmp_all,projTargets(pr2));
        diff1 = val2_m1-val1_m1; diff1=sum(diff1(diff1<0));
        diff2 = val2_m2-val1_m2; diff2=sum(diff2(diff2<0));
        diff3 = val2_m3-val1_m3; diff3=sum(diff3(diff3<0));
        diff_all=val2_all-val1_all; diff_all=sum(diff_all(diff_all<0));
        proj_fraction1=(sum(val1_m1)+diff1)/sum(val1_m1);
        proj_fraction2=(sum(val1_m2)+diff2)/sum(val1_m2);
        proj_fraction3=(sum(val1_m3)+diff3)/sum(val1_m3);
        proj_fraction_all=(sum(val1_all)+diff_all)/sum(val1_all);
        targetsMatrix1_m1(pr2,pr) = proj_fraction1;
        targetsMatrix1_m2(pr2,pr) = proj_fraction2;
        targetsMatrix1_m3(pr2,pr) = proj_fraction3;
        targetsMatrix1_all(pr2,pr) = proj_fraction_all;
    end
end


% %calculate mean from single matrices
% targetsMatrix1_mean=mean(cat(3, targetsMatrix1_m1,targetsMatrix1_m2,targetsMatrix1_m3),3,'omitnan');
% %set diagonal values to NaN
% targetsMatrix1_mean(1:11:end)= NaN;
% targetsMatrix1_m1(1:11:end)= NaN;
% targetsMatrix1_m2(1:11:end)= NaN;
% targetsMatrix1_m3(1:11:end)= NaN;

%or use the combined matrix
targetsMatrix1_all(1:11:end)= NaN;
targetsMatrix1_m1(1:11:end)= NaN;
targetsMatrix1_m2(1:11:end)= NaN;
targetsMatrix1_m3(1:11:end)= NaN;

figure(3)
%make nan values transparent
imAlpha=ones(size(targetsMatrix1_all));
imAlpha(isnan(targetsMatrix1_all))=0;
imagesc(targetsMatrix1_all,'AlphaData',imAlpha, [0 0.8]);
set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
    'color', [0.5 0.5 0.5])
xlabel('Area A', 'FontSize',10)
ylabel('Area B', 'FontSize',10)
col1=magma(300);
colormap(col1);

colorbar;
title('ipislateral SCs: subset of area A projections projecting to area B')

figure(4)
subplot(1,3,1)
imAlpha=ones(size(targetsMatrix1_m1));
imAlpha(isnan(targetsMatrix1_m1))=0;
imagesc(targetsMatrix1_m1,'AlphaData',imAlpha, [0 0.8]);
set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
    'color', [0.5 0.5 0.5])
xlabel('Area A', 'FontSize',8)
ylabel('Area B', 'FontSize',8)
colormap(col1);
colorbar;
title('ipislateral SCs M1: subset of area A projections projecting to area B','FontSize',8)

subplot(1,3,2)
imAlpha=ones(size(targetsMatrix1_m2));
imAlpha(isnan(targetsMatrix1_m2))=0;
imagesc(targetsMatrix1_m2,'AlphaData',imAlpha, [0 0.8]);
set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
    'color', [0.5 0.5 0.5])
xlabel('Area A', 'FontSize',8)
ylabel('Area B', 'FontSize',8)
colormap(col1);
colorbar;
title('ipislateral SCs M2: subset of area A projections projecting to area B','FontSize',8)

subplot(1,3,3)
imAlpha=ones(size(targetsMatrix1_m3));
imAlpha(isnan(targetsMatrix1_m3))=0;
imagesc(targetsMatrix1_m3,'AlphaData',imAlpha, [0 0.8]);
set(gca,'xticklabel',sample_names(projTargets), 'yticklabel',sample_names(projTargets),...
    'color', [0.5 0.5 0.5])
xlabel('Area A', 'FontSize',8)
ylabel('Area B', 'FontSize',8)
colormap(col1);
colorbar;
title('ipislateral SCs M3: subset of area A projections projecting to area B', 'FontSize',8)
set(gcf, 'Position', [-7.8000 397 1.5392e+03 302.4000])

%% calculate ratio of contra vs ipsilateral to estimate upper bound of false positive


areasleft=[34 32 36 38 30 ];
areasright=[35 33 37 39 31];

% %only deep SC
% areasleft=[16 17 20 21 24 25 28];
% areasright=[18 19 22 23 26 27 29];


left_right=sum(BigBnorm(:,areasright));
left_left=sum(BigBnorm(:,areasleft));
tot_sums=left_right+left_left;
left_right=left_right./tot_sums;
left_left=left_left./tot_sums;


bar([left_left', left_right']);
currentLabels=sample_names(areasleft);
newLabels = cellfun(@(label) label(1:end-1), currentLabels, 'UniformOutput', false);
set(gca,'xticklabel',newLabels,...
     'TickDir', 'out', 'box', 'off', 'FontSize',12)
ylabel('Projection strength', 'FontSize', 12)
title('ipsi vs contra projections')

% left_right=sum(BigBnorm(:,areasright),2);
% left_left=sum(BigBnorm(:,areasleft),2);
% r_left=left_right./left_left;
% r_left_f=r_left(left_right<left_left);
%
% r=r_left_f;
% rmean=mean(r);
% rmedian=median(r);
% rstd=std(r);
% r_sderror=std(r)/sqrt(length(r));
%
% %make a plot of normalized barcodecounts
%
% left_right=BigData(:,areasright);
% left_left=BigData(:,areasleft);
%
%
% ipsi=[left_left];
% contra=[left_right];
%
% ipsimean=mean(log2(sum(left_left,1)));
% ipsistd=std(log2(sum(left_left,1)));
%
% contramean=mean(log2(sum(left_right,1)));
% contrastd=std(log2(sum(left_right,1)));
%
% x=log2(sum(left_left,1));
% y=log2(sum(left_right,1));
%
% % save('summedcounts.mat','x','y')


%% calculate projection efficiency
DATA1 = B1;
DATA2 = B2;
DATA3 = B3;
DATA = BigBnorm;

% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<=cutoff;
DATA(cutoff_tmp)=0; %apply filter to normalize data
threshold=0;

scnt = 0; Numbers={};
Numbers{1} = zeros(length(sSC),2);Numbers{2} = zeros(length(sSC),2);
Numbers{3} = zeros(length(sSC),2); Numbers{4} = zeros(length(sSC),2);

for s = orderOfSC % go through injection sites
    scnt = scnt +1;
    id = sSC(s);
    BCodes1 = find(m1==id);
    BCodes2 = find(m2==id);
    BCodes3 = find(m3==id);
    BCodes = find(ismember(maxes,id));
    Numbers{1}(scnt,1) = length(BCodes1);
    Numbers{2}(scnt,1) = length(BCodes2);
    Numbers{3}(scnt,1) = length(BCodes3);
    Numbers{4}(scnt,1) = length(BCodes);
    % for each barcode find the # of areas
    areas1 = []; areas2 = []; areas3 = []; areas=[];
    for b = 1:length(BCodes1)
        curr = find(DATA1(BCodes1(b),1:34)>0);
        areas1 = [areas1;length(curr)];
    end
    for b = 1:length(BCodes2)
        curr = find(DATA2(BCodes2(b),1:34)>0);
        areas2 = [areas2;length(curr)];
    end
    for b = 1:length(BCodes3)
        curr = find(DATA3(BCodes3(b),1:34)>0);
        areas3 = [areas3;length(curr)];
    end
    for b = 1:length(BCodes)
        curr = find(DATA(BCodes(b),1:34)>0);
        areas = [areas;length(curr)];
    end
    %     figure
    [~,a2] = sort(areas1,'descend');
    data1 = DATA1(BCodes1(a2),1:39);
    [~,a2] = sort(areas2,'descend');
    data2 = DATA2(BCodes2(a2),1:39);
    [~,a2] = sort(areas3,'descend');
    data3 = DATA3(BCodes3(a2),1:39);
    [~,ai] = sort(areas,'descend');
    data_all = DATA(BCodes(ai),1:39);
    keepData1 = data1; keepData2 = data2; keepData3 = data3;

    % targetsMatrix2_m1 = zeros(length(projTargets),1);
    % targetsMatrix2_m2 = zeros(length(projTargets),1);
    % targetsMatrix2_m3 = zeros(length(projTargets),1);


    %for efficiency to each area (can go over 100% if summed)
    % for pr = 1:length(projTargets)
    %     tmp1 = find(keepData1(:,projTargets(pr))>0);
    %     targetsMatrix2_m1(pr,1) = length(tmp1);
    %
    % end

    %for total efficiency, same barcode in different areas is counted once

    tmp1 = keepData1(:,projTargets)>threshold;
    tmp_zero1=sum(all(tmp1 == 0, 2));
    tmp_nnz1=length(BCodes1)-tmp_zero1;

    tmp2 = keepData2(:,projTargets)>threshold;
    tmp_zero2=sum(all(tmp2 == 0, 2));
    tmp_nnz2=length(BCodes2)-tmp_zero2;

    tmp3 = keepData3(:,projTargets)>threshold;
    tmp_zero3=sum(all(tmp3 == 0, 2));
    tmp_nnz3=length(BCodes3)-tmp_zero3;

    tmp = data_all(:,projTargets)>threshold;
    tmp_zero=sum(all(tmp == 0, 2));
    tmp_nnz=length(BCodes)-tmp_zero;


    Numbers{1}(scnt,2) = tmp_nnz1;
    Numbers{2}(scnt,2) = tmp_nnz2;
    Numbers{3}(scnt,2) = tmp_nnz3;
    Numbers{4}(scnt,2) = tmp_nnz;

end

% %plot it
figure(1)
proj_eff(:,1)=Numbers{1}(:,2)./Numbers{1}(:,1);
proj_eff(:,2)=Numbers{2}(:,2)./Numbers{2}(:,1);
proj_eff(:,3)=Numbers{3}(:,2)./Numbers{3}(:,1);
proj_eff(:,4)=Numbers{4}(:,2)./Numbers{4}(:,1);


% bar(mean(proj_eff,2)); hold on
% scatter(1:7, proj_eff, 20, 'filled'); hold on

plot(1:7,proj_eff(:,1),'-o','Color', '#6A2C70', LineWidth=2); hold on
plot(1:7,proj_eff(:,2),'-o','Color', '#B83B5E', LineWidth=2); hold on
plot(1:7,proj_eff(:,3),'-o','Color', '#FFAF45', LineWidth=2)
xlim([0.5 7.5]); xticks(1:7);
ylim([0 1]);
set(gca,'xticklabel',string(SampleName{3}(sSC)), 'FontSize', 8);
xtickangle(45)
ylabel('projection efficiency', 'FontSize',10)
title('Efficiency of barcode projection', 'FontSize', 12)
legend('m1', 'm2', 'm3', 'FontSize', 12);
set(gcf, 'Position', [421 217 580 436])

numbers_med=Numbers{4}([2 4 6],:); numbers_med=sum(numbers_med,1);
numbers_lat=Numbers{4}([1 3 5],:); numbers_lat=sum(numbers_lat,1);
proj_eff_med=numbers_med(:,2)./numbers_med(:,1);
proj_eff_lat=numbers_lat(:,2)./numbers_lat(:,1);

coloritos=[col(3,:);col(4,:)];
figure(10); set(gcf, 'position', [488 242 193.8000 335.8000]);
h=bar(1:2,[proj_eff_med proj_eff_lat]', 'FaceAlpha', 0.7, 'EdgeColor', 'none');
h.FaceColor='flat';
h.CData=coloritos;
set(gca,'xticklabel',{'SCs med', 'SCs lat'},...
     'TickDir', 'out', 'box', 'off')
ylabel('Projection efficiency', 'FontSize', 11)
xticks(1:2); xlim([0 3]); ylim([0 1])
xtickangle(45)


%% all co-labeled projections together: barcode counts

allKeep=[allKeep1;allKeep2;allKeep3];
tmp = find(allKeep<threshold);
allKeep(tmp) = 0;

targetsMatrix = zeros(length(projTargets),length(projTargets));
targetsMatrix2 = zeros(length(projTargets),1);
for pr = 1:length(projTargets)
    tmp = find(allKeep(:,projTargets(pr))>0);
    targetsMatrix2(pr,1) = length(tmp);
    for pr2 = pr:length(projTargets)
        tmp2 = find(allKeep(tmp,projTargets(pr2))>0);
        targetsMatrix(pr,pr2) = 1/length(tmp)*length(tmp2);
        targetsMatrix(pr2,pr) = 1/length(tmp)*length(tmp2);

    end
end
figure
imagesc(targetsMatrix)
set(gca,'xticklabel',SampleName{3}(projTargets))
set(gca,'yticklabel',SampleName{3}(projTargets))
col1=magma(300);
colormap(col1); colorbar;
title('co-labelled sSC projections')

% figure
% imagesc(targetsMatrix,[0 0.4])
% set(gca,'xticklabel',SampleName{3}(projTargets))
% set(gca,'yticklabel',SampleName{3}(projTargets))
% col1=magma(300);
% colormap(col1); colorbar;
% title('co-labelled sSC projections')

%% medial vs lateral for each main SCs region - projections normalized to SCs region

DATA1 = B1norm;
DATA2 = B2norm;
DATA3 = B3norm;
DATA = BigBnorm;

% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
threshold=0;

%group SCaa and SCa into SCanterior and SCcm into posterior
custom_id{1}=[3 7]; custom_id{2}=[2 6];
custom_id{3}=11; custom_id{4}=10; custom_id{5}=14;

scnt = 0;
coll1 = []; coll2 = []; coll3 = []; coll_all=[];
Numbers={};
Numbers{1} = zeros(length(sSC),2);Numbers{2} = zeros(length(sSC),2);
Numbers{3} = zeros(length(sSC),2); allKeep1 =[]; allKeep2=[]; allKeep3=[]; allKeep=[];
for s = 1:length(custom_id) % go through injection sites
    scnt = scnt +1;
    id = custom_id{s};
    BCodes1 = find(ismember(m1,id));
    BCodes2 = find(ismember(m2,id));
    BCodes3 = find(ismember(m3,id));
    BCodes = find(ismember(maxes,id));
    Numbers{1}(scnt,1) = length(BCodes1);
    Numbers{2}(scnt,1) = length(BCodes2);
    Numbers{3}(scnt,1) = length(BCodes3);
    % for each barcode find the # of areas
    areas1 = []; areas2 = []; areas3 = []; areas=[];
    for b = 1:length(BCodes1)
        curr = find(DATA1(BCodes1(b),1:34)>threshold);
        areas1 = [areas1;length(curr)];
    end
    for b = 1:length(BCodes2)
        curr = find(DATA2(BCodes2(b),1:34)>threshold);
        areas2 = [areas2;length(curr)];
    end
    for b = 1:length(BCodes3)
        curr = find(DATA3(BCodes3(b),1:34)>threshold);
        areas3 = [areas3;length(curr)];
    end
    for b = 1:length(BCodes)
        curr = find(DATA(BCodes(b),1:34)>threshold);
        areas = [areas;length(curr)];
    end
    %     figure
    [~,a2] = sort(areas1,'descend');
    data1 = DATA1(BCodes1(a2),1:34);
    [~,a2] = sort(areas2,'descend');
    data2 = DATA2(BCodes2(a2),1:34);
    [~,a2] = sort(areas3,'descend');
    data3 = DATA3(BCodes3(a2),1:34);
    [~,ai] = sort(areas,'descend');
    data_all = DATA(BCodes(ai),1:34);
    keepData1 = data1; keepData2 = data2; keepData3 = data3;
    allKeep1 = [allKeep1; keepData1];
    allKeep2 = [allKeep2; keepData2];
    allKeep3 = [allKeep3; keepData3];
    allKeep = [allKeep; data_all];


    % tmp = find(keepData3<threshold);
    % keepData3(tmp) = 0;

    targetsMatrix2_m1 = zeros(length(projTargets),1);
    targetsMatrix2_m2 = zeros(length(projTargets),1);
    targetsMatrix2_m3 = zeros(length(projTargets),1);
    targetsMatrix2_all = zeros(length(projTargets),1);

    for pr = 1:length(projTargets)
        targetsMatrix2_m1(pr,1) = sum(keepData1(:,projTargets(pr)));
        targetsMatrix2_m2(pr,1) = sum(keepData2(:,projTargets(pr)));
        targetsMatrix2_m3(pr,1) = sum(keepData3(:,projTargets(pr)));
        targetsMatrix2_all(pr,1) = sum(data_all(:,projTargets(pr)));
    end

    %normalize
    targetsMatrix2_m1=targetsMatrix2_m1./sum(targetsMatrix2_m1);
    targetsMatrix2_m2=targetsMatrix2_m2./sum(targetsMatrix2_m2);
    targetsMatrix2_m3=targetsMatrix2_m3./sum(targetsMatrix2_m3);
    targetsMatrix2_all=targetsMatrix2_all./sum(targetsMatrix2_all);

    targetsMatrix1_m1 = zeros(length(projTargets),length(projTargets));
    targetsMatrix1_m2 = zeros(length(projTargets),length(projTargets));
    targetsMatrix1_m3 = zeros(length(projTargets),length(projTargets));
    targetsMatrix1_all = zeros(length(projTargets),length(projTargets));
    for pr = 1:length(projTargets)
        tmp1 = find(allKeep1(:,projTargets(pr)));
        tmp2 = find(allKeep2(:,projTargets(pr)));
        tmp3 = find(allKeep3(:,projTargets(pr)));
        tmp_all=find(allKeep(:,projTargets(pr)));
        val1_m1 = allKeep1(tmp1,projTargets(pr));
        val1_m2 = allKeep2(tmp2,projTargets(pr));
        val1_m3 = allKeep3(tmp3,projTargets(pr));
        val1_all= allKeep(tmp_all,projTargets(pr));
        for pr2 = 1:length(projTargets)
            val2_m1 = allKeep1(tmp1,projTargets(pr2));
            val2_m2 = allKeep2(tmp2,projTargets(pr2));
            val2_m3 = allKeep3(tmp3,projTargets(pr2));
            val2_all= allKeep(tmp_all,projTargets(pr2));
            diff1 = val2_m1-val1_m1; diff1=sum(diff1(diff1<0));
            diff2 = val2_m2-val1_m2; diff2=sum(diff2(diff2<0));
            diff3 = val2_m3-val1_m3; diff3=sum(diff3(diff3<0));
            diff_all=val2_all-val1_all; diff_all=sum(diff_all(diff_all<0));
            proj_fraction1=(sum(val1_m1)+diff1)/sum(val1_m1);
            proj_fraction2=(sum(val1_m2)+diff2)/sum(val1_m2);
            proj_fraction3=(sum(val1_m3)+diff3)/sum(val1_m3);
            proj_fraction_all=(sum(val1_all)+diff_all)/sum(val1_all);
            targetsMatrix1_m1(pr2,pr) = proj_fraction1;
            targetsMatrix1_m2(pr2,pr) = proj_fraction2;
            targetsMatrix1_m3(pr2,pr) = proj_fraction3;
            targetsMatrix1_all(pr2,pr) = proj_fraction_all;
        end
    end


    Numbers{1}(scnt,2) = sum(targetsMatrix2_m1);
    Numbers{2}(scnt,2) = sum(targetsMatrix2_m2);
    Numbers{3}(scnt,2) = sum(targetsMatrix2_m3);
    coll1 = [coll1 targetsMatrix2_m1];
    coll2 = [coll2 targetsMatrix2_m2];
    coll3 = [coll3 targetsMatrix2_m3];
    coll_all = [coll_all targetsMatrix2_all];

end


% coll=(mean(cat(3,coll1,coll2,coll3),3,'omitnan'))'; %this is from the mean of the 3 mice matrices
% coll=coll(:,iy);

% coll_all=coll_all'; %this is from the big combined matrix of all mice


% coll_med = [coll(2,:); coll(4,:); coll(6,:)]; %SCc has only been injected in medial site, so leave it out
% coll_lat = [coll(1,:); coll(3,:); coll(5,:)];
% coll_med = [coll_all(2,:); coll_all(4,:); coll_all(6,:)]; %SCc has only been injected in medial site, so leave it out
% coll_lat = [coll_all(1,:); coll_all(3,:); coll_all(5,:)];

coll_med = [coll_all(:,1) coll_all(:,3)];
coll_lat = [coll_all(:,2) coll_all(:,4)];


% figure('position',[489 277 700 312]);
names=sample_names(projTargets);

%all together
figure(1); set(gcf, 'position', [785.8000 360.2000 430.4000 312])
h1=imagesc(coll_all);
set(gca,'xtick',1:5,'xticklabel',{'SCam','SCal','SCpm','SCpl', 'SCpp'}, 'TickLength', [0 0])
set(gca,'ytick',1:5,'yticklabel',names)
title('Projection strength of main SCs regions')
col1=magma(300);
colormap(h1.Parent, col1);
cb = colorbar;
ylabel(cb, 'Projection strength');

% %medial
% subplot(2,2,1)
% imagesc(coll_med)
% set(gca,'xtick',1:10,'xticklabel',SampleName{3}(projTargets))
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCam','SCcm','SCc'})
% title('medial')
% col1=magma(300);
% colormap(col1);
% cb = colorbar;
% ylabel(cb, 'Projection strength');
%
% %lateral
% subplot(2,2,3)
% imagesc(coll_lat)
% set(gca,'xtick',1:10,'xticklabel',SampleName{3}(projTargets))
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCam','SCcm','SCc'})
% title('lateral')
% cb = colorbar;
% ylabel(cb, 'Projection strength');

% %absolute difference
% coll_abs=abs(coll_med-coll_lat);
% subplot(2,1,1)
% h=imagesc(coll_abs);
% set(gca,'xtick',1:10,'xticklabel',SampleName{3}(projTargets))
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCam','SCcm','SCc'})
% title('|medial-lateral|')
% colormap(h.Parent, div_mgblackgr2);
% cb = colorbar;
% ylabel(cb, 'Absolute difference');
%
% %projection difference
% coll_diff=coll_med-coll_lat;
% subplot(3,1,1)
% h=imagesc(coll_diff);
% set(gca,'xtick',1:10,'xticklabel',names)
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCam','SCcm','SCc'})
% title('medial-lateral')
% colormap(h.Parent, div_mgblackgr2);
% cb = colorbar; clim([-0.3 0.3])
% ylabel(cb, 'Projection difference');

% %color-mapped percent difference ---- it's the same as bias index
% coll_diff=((coll_med-coll_lat)./(coll_med+coll_lat))*100;
% mx=max(coll_diff, [], 'all'); mn=min(coll_diff, [], 'all');
% max_val=max([mx abs(mn)]);
% subplot(4,1,3)
% h1=imagesc(coll_diff, [-max_val max_val]);
% set(gca,'xtick',1:5,'xticklabel',names, 'TickLength', [0 0])
% set(gca,'ytick',1:5,'yticklabel',{'SCa','SCp','SCcm','SCc'})
% title('% medial-lateral')
% colormap(h1.Parent, div_cmap2);
% cb = colorbar;
% ylabel(cb, 'Percent difference');

%color-mapped bias index
figure(2); set(gcf, 'position', [534.6000 161.8000 211.2000 360.8000])
coll_bias=log2(coll_med./coll_lat);
mx=ceil(max(coll_bias(coll_bias~= Inf)));
mn=floor(min(coll_bias(coll_bias~= -Inf)));
inf_val=max(mx, abs(mn));
mx_id = coll_bias == Inf; coll_bias(mx_id)=inf_val;
mn_id = coll_bias == -Inf; coll_bias(mn_id)=-inf_val;
h2=imagesc(coll_bias, [-inf_val inf_val]);
set(gca,'xtick',1:5,'xticklabel',{'SCa','SCp','SCcm','SCc'}, 'TickLength', [0 0],...
    'FontSize', 12)
set(gca,'ytick',1:5,'yticklabel',names, 'FontSize', 12)
title('upper/lower bias', 'FontSize', 12)
colormap(h2.Parent, div_cmap_grp3);
cb = colorbar;
ylabel(cb, 'Bias index','FontSize', 11);

%% medial vs lateral for grouped SCs regions

DATA1 = B1norm;
DATA2 = B2norm;
DATA3 = B3norm;
DATA = BigBnorm;

% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<=cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
threshold=0;

id_med=[3 7 11]; %medial SC
id_lat=[2 6 10]; %lateral SC
id_antero=[2 3 6 7];
id_postero=[10 11 14];
id_tot={id_med;id_lat;id_antero;id_postero};
names_medlat={'medial SCs','lateral SCs'};

scnt = 0;
coll1 = []; coll2 = []; coll3 = []; coll_all=[];
Numbers={};
Numbers{1} = zeros(length(sSC),2);Numbers{2} = zeros(length(sSC),2);
Numbers{3} = zeros(length(sSC),2); allKeep1 =[]; allKeep2=[]; allKeep3=[]; allKeep=[];
for s = 1:size(id_tot,1) % go through medial and lateral injection sites
    scnt = scnt +1;
    id = id_tot{s};
    BCodes1 = find(ismember(m1,id));
    BCodes2 = find(ismember(m2,id));
    BCodes3 = find(ismember(m3,id));
    BCodes = find(ismember(maxes,id));
    Numbers{1}(scnt,1) = length(BCodes1);
    Numbers{2}(scnt,1) = length(BCodes2);
    Numbers{3}(scnt,1) = length(BCodes3);
    % for each barcode find the # of areas
    areas1 = []; areas2 = []; areas3 = []; areas=[];
    for b = 1:length(BCodes1)
        curr = find(DATA1(BCodes1(b),1:size(B1norm,2))>threshold);
        areas1 = [areas1;length(curr)];
    end
    for b = 1:length(BCodes2)
        curr = find(DATA2(BCodes2(b),1:size(B2norm,2))>threshold);
        areas2 = [areas2;length(curr)];
    end
    for b = 1:length(BCodes3)
        curr = find(DATA3(BCodes3(b),1:size(B3norm,2))>threshold);
        areas3 = [areas3;length(curr)];
    end
    for b = 1:length(BCodes)
        curr = find(DATA(BCodes(b),1:size(BigBnorm,2))>threshold);
        areas = [areas;length(curr)];
    end
    %     figure
    [~,a2] = sort(areas1,'descend');
    data1 = DATA1(BCodes1(a2),1:size(B1norm,2));
    [~,a2] = sort(areas2,'descend');
    data2 = DATA2(BCodes2(a2),1:size(B2norm,2));
    [~,a2] = sort(areas3,'descend');
    data3 = DATA3(BCodes3(a2),1:size(B3norm,2));
    [~,ai] = sort(areas,'descend');
    data_all = DATA(BCodes(ai),1:size(BigBnorm,2));
    keepData1 = data1; keepData2 = data2; keepData3 = data3;
    allKeep1 = [allKeep1; keepData1];
    allKeep2 = [allKeep2; keepData2];
    allKeep3 = [allKeep3; keepData3];
    allKeep = [allKeep; data_all];


    % tmp = find(keepData3<threshold);
    % keepData3(tmp) = 0;

    targetsMatrix2_m1 = zeros(length(projTargets),1);
    targetsMatrix2_m2 = zeros(length(projTargets),1);
    targetsMatrix2_m3 = zeros(length(projTargets),1);
    targetsMatrix2_all = zeros(length(projTargets),1);

    for pr = 1:length(projTargets)
        targetsMatrix2_m1(pr,1) = sum(keepData1(:,projTargets(pr)));
        targetsMatrix2_m2(pr,1) = sum(keepData2(:,projTargets(pr)));
        targetsMatrix2_m3(pr,1) = sum(keepData3(:,projTargets(pr)));
        targetsMatrix2_all(pr,1) = sum(data_all(:,projTargets(pr)));
    end

    %normalize
    targetsMatrix2_m1=targetsMatrix2_m1./sum(targetsMatrix2_m1);
    targetsMatrix2_m2=targetsMatrix2_m2./sum(targetsMatrix2_m2);
    targetsMatrix2_m3=targetsMatrix2_m3./sum(targetsMatrix2_m3);
    targetsMatrix2_all=targetsMatrix2_all./sum(targetsMatrix2_all);


    Numbers{1}(scnt,2) = sum(targetsMatrix2_m1);
    Numbers{2}(scnt,2) = sum(targetsMatrix2_m2);
    Numbers{3}(scnt,2) = sum(targetsMatrix2_m3);
    coll1 = [coll1 targetsMatrix2_m1];
    coll2 = [coll2 targetsMatrix2_m2];
    coll3 = [coll3 targetsMatrix2_m3];
    coll_all = [coll_all targetsMatrix2_all];

end

%separate ipsi and contra;
% iy=[1 3 5 7 9 2 4 6 8 10];
% iy=[1 3 5 7 9 2 4 6 8 10 11:17 18 19 22 23 26 27 30 20 21 24 25 28 29 31];
% coll=(mean(cat(3,coll1,coll2,coll3),3,'omitnan'))'; %this is from the mean of the 3 mice matrices
% coll=coll(:,iy);

coll_all=coll_all'; %this is from the big combined matrix of all mice
% coll_all=coll_all(:,iy);


coll_med = coll_all(1,:); %SCc has only been injected in medial site, so leave it out
coll_lat = coll_all(2,:);

coll_ant= coll_all(3,:);
coll_post= coll_all(4,:);

% coll_med = coll(1,:); %SCc has only been injected in medial site, so leave it out
% coll_lat = coll(2,:);


figure('position',[422.6000 522.6000 593.6000 150.4000]);
names=sample_names(projTargets);


% %color-mapped percent difference M-L
% coll_diff=((coll_med-coll_lat)./(coll_med+coll_lat))*100;
% mx=max(coll_diff, [], 'all'); mn=min(coll_diff, [], 'all');
% max_val=max([mx abs(mn)]);
% subplot(2,1,1)
% h1=imagesc(coll_diff, [-max_val max_val]);
% set(gca,'xtick',1:31,'xticklabel',names, 'TickLength', [0 0])
% set(gca,'ytick',[])
% title('% medial-lateral')
% colormap(h1.Parent, div_cmap2);
% cb = colorbar; clim([-100 100])
% ylabel(cb, 'Percent difference');
%
% %color-mapped percent difference A-P
% coll_diff2=((coll_ant-coll_post)./(coll_ant+coll_post))*100;
% mx=max(coll_diff2, [], 'all'); mn=min(coll_diff2, [], 'all');
% max_val=max([mx abs(mn)]);
% subplot(2,1,2)
% h2=imagesc(coll_diff2, [-max_val max_val]);
% set(gca,'xtick',1:31,'xticklabel',names, 'TickLength', [0 0])
% set(gca,'ytick',[])
% title('% anterior-posterior')
% colormap(h2.Parent, div_cmap_grp2);
% cb = colorbar; clim([-100 100])
% ylabel(cb, 'Percent difference');

%color-mapped bias index M-L
subplot(2,1,1)
coll_bias=log2(coll_med./coll_lat);
mx=ceil(max(coll_bias(coll_bias~= Inf)));
mn=floor(min(coll_bias(coll_bias~= -Inf)));
inf_val=max(mx, abs(mn));
mx_id = coll_bias == Inf; coll_bias(mx_id)=inf_val;
mn_id = coll_bias == -Inf; coll_bias(mn_id)=-inf_val;
h1=imagesc(coll_bias, [-inf_val inf_val]);
set(gca,'xtick',1:31,'xticklabel',names, 'TickLength', [0 0])
set(gca,'ytick',[])
title('medial/lateral bias')
colormap(h1.Parent, div_cmap3);
cb = colorbar;
ylabel(cb, 'Bias index');

%color-mapped bias index A-P
subplot(2,1,2)
coll_bias2=log2(coll_ant./coll_post);
mx=ceil(max(coll_bias2(coll_bias2~= Inf)));
mn=floor(min(coll_bias2(coll_bias2~= -Inf)));
inf_val=max(mx, abs(mn));
mx_id = coll_bias2 == Inf; coll_bias2(mx_id)=inf_val;
mn_id = coll_bias2 == -Inf; coll_bias2(mn_id)=-inf_val;
h2=imagesc(coll_bias2, [-inf_val inf_val]);
set(gca,'xtick',1:31,'xticklabel',names, 'TickLength', [0 0])
set(gca,'ytick',[])
title('anterior/posterior bias')
colormap(h2.Parent, div_cmap_grp3);
cb = colorbar;
ylabel(cb, 'Bias index');

%% horizontal bar plot of medial-lateral and antero-caudal bias index

figure(1); set(gcf, 'position', [305 350.6000 743 311.4000]);

%MEDIAL-LATERAL
subplot(1,2,1)
h=barh(1:length(projTargets),coll_bias,'BarWidth',0.7, 'FaceColor', 'flat');%plot bars at mean values
max_val=max(abs(coll_bias));
% max_val=100;
normalized_vals = (coll_bias + max_val) / (max_val*2);
colormap(div_cmap_grp3);
cmap = colormap;  % Get the current colormap

% Map the normalized data to colormap
num_colors = size(cmap, 1);
color_indices = round(normalized_vals * (num_colors - 1)) + 1;
color_indices(isnan(color_indices))=1;
h.CData = cmap(color_indices, :);
xlim(round([-max_val, max_val]));
set(gca,'TickLength',[0.0025 0.0025],'TickDir', 'out', 'FontSize',12,'FontName', 'Arial',...
    'GridColor',[.5 .5 .5],'GridLineStyle', '--','XGrid', 'on','color', 'none',...
    'box', 'off', 'YDir', 'reverse', 'XDir', 'reverse', 'xtick', -3:3)
ax = gca; ax.GridAlpha = 0.5;YMinorGrid = 'on';
ax.YAxis.Visible = 'off';
xlabel(["\fontsize{13}Bias index";"\fontsize{12}log_{2}(%up/%low)"], 'FontName', 'Arial');
% xlabel("Percent difference");

%add texts with y tick labels
text(repmat(round(max_val)+0.1,1,length(projTargets)),1:length(projTargets),names, 'HorizontalAlignment', 'right', 'FontSize', 13);

text(-max_val/2,0,'Lower','Color',col(4,:),'fontsize',13,'FontWeight','bold', 'HorizontalAlignment', 'center');
text(max_val/2,0,'Upper','Color',col(3,:),'fontsize',13,'FontWeight','bold',  'HorizontalAlignment', 'center');


%ANTERIOR-POSTERIOR
subplot(1,2,2)
h=barh(1:length(projTargets),coll_bias2,'BarWidth',0.7, 'FaceColor', 'flat');%plot bars at mean values
max_val=max(abs(coll_bias2));
% max_val=100;
normalized_vals = (coll_bias2 + max_val) / (max_val*2);
% cols1=cbrewer('div', 'PiYG',300);
% cols1=flipud(cols1);
colormap(div_cmap3);
cmap = colormap;  % Get the current colormap
uistack(h, 'top');


% Map the normalized data to colormap
num_colors = size(cmap, 1);
color_indices = round(normalized_vals * (num_colors - 1)) + 1;
color_indices(isnan(color_indices))=1;
h.CData = cmap(color_indices, :);
xlim(round([-max_val, max_val]));
set(gca,'TickLength',[0.0025 0.0025],'TickDir', 'out', 'FontSize',12,'FontName', 'Arial',...
    'GridColor',[.5 .5 .5],'GridLineStyle', '--','XGrid', 'on','color', 'none',...
    'box', 'off', 'YDir', 'reverse', 'XDir', 'reverse', 'xtick', -2:2)
ax = gca; ax.GridAlpha = 0.5;YMinorGrid = 'on';
ax.YAxis.Visible = 'off';
xlabel(["\fontsize{13}Bias index";"\fontsize{12}log_{2}(%nas/%temp)"], 'FontName', 'Arial');
% xlabel("Percent difference");

%add texts with y tick labels
text(repmat(round(max_val)+0.1,1,length(projTargets)),1:length(projTargets),names, 'HorizontalAlignment', 'right', 'FontSize', 13);

text(-max_val/2,0,'Temporal','Color',col(2,:),'fontsize',13,'FontWeight','bold', 'HorizontalAlignment', 'center');
text(max_val/2,0,'Nasal','Color',col(1,:),'fontsize',13,'FontWeight','bold',  'HorizontalAlignment', 'center');


%% medial vs lateral for each main SCs region - projections normalized to target region

DATA1 = B1norm;
DATA2 = B2norm;
DATA3 = B3norm;
DATA = BigBnorm;

% BCodes = find(ismember(maxes,sSC));
data=BigB;
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes
cutoff_tmp=BigB(:,:)<cutoff; %cutoff based on non-norm data because of integer values
DATA(cutoff_tmp)=0;
threshold=0;


%group SCaa and SCa into SCanterior and SCcm into posterior
custom_id{1}=[3 7]; custom_id{2}=[2 6];
custom_id{3}=11; custom_id{4}=10; custom_id{5}=14;

scnt = 0;
coll1 = []; coll2 = []; coll3 = []; coll_all=[];
Numbers={};
Numbers{1} = zeros(length(sSC),2);Numbers{2} = zeros(length(sSC),2);
Numbers{3} = zeros(length(sSC),2); allKeep1 =[]; allKeep2=[]; allKeep3=[]; allKeep=[];
for s = 1:length(custom_id) % go through injection sites
    scnt = scnt +1;
    id = custom_id{s};
    BCodes1 = find(ismember(m1,id));
    BCodes2 = find(ismember(m2,id));
    BCodes3 = find(ismember(m3,id));
    BCodes = find(ismember(maxes,id));
    Numbers{1}(scnt,1) = length(BCodes1);
    Numbers{2}(scnt,1) = length(BCodes2);
    Numbers{3}(scnt,1) = length(BCodes3);
    % for each barcode find the # of areas
    areas1 = []; areas2 = []; areas3 = []; areas=[];
    for b = 1:length(BCodes1)
        curr = find(DATA1(BCodes1(b),1:34)>threshold);
        areas1 = [areas1;length(curr)];
    end
    for b = 1:length(BCodes2)
        curr = find(DATA2(BCodes2(b),1:34)>threshold);
        areas2 = [areas2;length(curr)];
    end
    for b = 1:length(BCodes3)
        curr = find(DATA3(BCodes3(b),1:34)>threshold);
        areas3 = [areas3;length(curr)];
    end
    for b = 1:length(BCodes)
        curr = find(DATA(BCodes(b),1:34)>threshold);
        areas = [areas;length(curr)];
    end
    %     figure
    [~,a2] = sort(areas1,'descend');
    data1 = DATA1(BCodes1(a2),1:34);
    [~,a2] = sort(areas2,'descend');
    data2 = DATA2(BCodes2(a2),1:34);
    [~,a2] = sort(areas3,'descend');
    data3 = DATA3(BCodes3(a2),1:34);
    [~,ai] = sort(areas,'descend');
    data_all = DATA(BCodes(ai),1:34);
    keepData1 = data1; keepData2 = data2; keepData3 = data3;
    allKeep1 = [allKeep1; keepData1];
    allKeep2 = [allKeep2; keepData2];
    allKeep3 = [allKeep3; keepData3];
    allKeep = [allKeep; data_all];


    % tmp = find(keepData3<threshold);
    % keepData3(tmp) = 0;

    targetsMatrix2_m1 = zeros(length(projTargets),1);
    targetsMatrix2_m2 = zeros(length(projTargets),1);
    targetsMatrix2_m3 = zeros(length(projTargets),1);
    targetsMatrix2_all = zeros(length(projTargets),1);

    for pr = 1:length(projTargets)
        targetsMatrix2_m1(pr,1) = sum(keepData1(:,projTargets(pr)));
        targetsMatrix2_m2(pr,1) = sum(keepData2(:,projTargets(pr)));
        targetsMatrix2_m3(pr,1) = sum(keepData3(:,projTargets(pr)));
        targetsMatrix2_all(pr,1) = sum(data_all(:,projTargets(pr)));
    end

    % %normalize
    % targetsMatrix2_m1=targetsMatrix2_m1./sum(targetsMatrix2_m1);
    % targetsMatrix2_m2=targetsMatrix2_m2./sum(targetsMatrix2_m2);
    % targetsMatrix2_m3=targetsMatrix2_m3./sum(targetsMatrix2_m3);
    % targetsMatrix2_all=targetsMatrix2_all./sum(targetsMatrix2_all);

    coll1 = [coll1 targetsMatrix2_m1];
    coll2 = [coll2 targetsMatrix2_m2];
    coll3 = [coll3 targetsMatrix2_m3];
    coll_all = [coll_all targetsMatrix2_all];

end

%normalize based on target region
coll_all=coll_all./sum(coll_all,2);


% coll=(mean(cat(3,coll1,coll2,coll3),3,'omitnan'))'; %this is from the mean of the 3 mice matrices
% coll=coll(:,iy);

coll_all=coll_all'; %this is from the big combined matrix of all mice


% coll_med = [coll(2,:); coll(4,:); coll(6,:)]; %SCc has only been injected in medial site, so leave it out
% coll_lat = [coll(1,:); coll(3,:); coll(5,:)];
% coll_med = [coll_all(2,:); coll_all(4,:); coll_all(6,:)]; %SCc has only been injected in medial site, so leave it out
% coll_lat = [coll_all(1,:); coll_all(3,:); coll_all(5,:)];

coll_med = [coll_all(1,:); coll_all(3,:)];
coll_lat = [coll_all(2,:); coll_all(4,:)];


names=sample_names(projTargets);

%all together
figure(1); set(gcf, 'position', [888.2000 248.2000 412 312])
h1=imagesc(coll_all);
set(gca,'xtick',1:10,'xticklabel',names, 'TickLength', [0 0])
set(gca,'ytick',1:5,'yticklabel',{'SCam','SCal','SCpm','SCpl', 'SCpp'})
title('Projection strength from SCs regions')
col1=magma(300);
colormap(h1.Parent, col1);
cb = colorbar;
ylabel(cb, 'Projection strength');

% %medial
% subplot(2,2,1)
% imagesc(coll_med)
% set(gca,'xtick',1:10,'xticklabel',SampleName{3}(projTargets))
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCam','SCcm','SCc'})
% title('medial')
% col1=magma(300);
% colormap(col1);
% cb = colorbar;
% ylabel(cb, 'Projection strength');
%
% %lateral
% subplot(2,2,3)
% imagesc(coll_lat)
% set(gca,'xtick',1:10,'xticklabel',SampleName{3}(projTargets))
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCam','SCcm','SCc'})
% title('lateral')
% cb = colorbar;
% ylabel(cb, 'Projection strength');

% %absolute difference
% coll_abs=abs(coll_med-coll_lat);
% subplot(2,1,1)
% h=imagesc(coll_abs);
% set(gca,'xtick',1:10,'xticklabel',SampleName{3}(projTargets))
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCam','SCcm','SCc'})
% title('|medial-lateral|')
% colormap(h.Parent, div_mgblackgr2);
% cb = colorbar;
% ylabel(cb, 'Absolute difference');
%
% %projection difference
% coll_diff=coll_med-coll_lat;
% subplot(3,1,1)
% h=imagesc(coll_diff);
% set(gca,'xtick',1:10,'xticklabel',names)
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCam','SCcm','SCc'})
% title('medial-lateral')
% colormap(h.Parent, div_mgblackgr2);
% cb = colorbar; clim([-0.3 0.3])
% ylabel(cb, 'Projection difference');

% %color-mapped percent difference ---- it's the same as bias index
% coll_diff=((coll_med-coll_lat)./(coll_med+coll_lat))*100;
% mx=max(coll_diff, [], 'all'); mn=min(coll_diff, [], 'all');
% max_val=max([mx abs(mn)]);
% subplot(4,1,3)
% h1=imagesc(coll_diff, [-max_val max_val]);
% set(gca,'xtick',1:10,'xticklabel',names, 'TickLength', [0 0])
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCp','SCcm','SCc'})
% title('% medial-lateral')
% colormap(h1.Parent, div_cmap2);
% cb = colorbar;
% ylabel(cb, 'Percent difference');
%
% %color-mapped bias index
% subplot(4,1,4)
% coll_bias=log2(coll_med./coll_lat);
% mx=ceil(max(coll_bias(coll_bias~= Inf)));
% mn=floor(min(coll_bias(coll_bias~= -Inf)));
% inf_val=max(mx, abs(mn));
% mx_id = coll_bias == Inf; coll_bias(mx_id)=inf_val;
% mn_id = coll_bias == -Inf; coll_bias(mn_id)=-inf_val;
% h2=imagesc(coll_bias, [-inf_val inf_val]);
% set(gca,'xtick',1:10,'xticklabel',names, 'TickLength', [0 0])
% set(gca,'ytick',1:4,'yticklabel',{'SCa','SCp','SCcm','SCc'})
% title('medial/lateral bias')
% colormap(h2.Parent, div_cmap2);
% cb = colorbar;
% ylabel(cb, 'Bias index');
%% calculate projection motifs

% projTargets = [34 35 32 33 36 37 38 39 30 31, ...
%     16 17 20 21 24 25 28]; %add SC target areas
motif_start=1; %0 to include motifs not projecting to any area,
% 1 to exclude those not projecting to any area
% length(projTargets)+1 to exclude those projecting to <2 areas

% id=sSC;
id_med=[3; 7; 11]; %medial SC
id_lat=[2; 6; 10]; %lateral SC
id_ant=[2 3 6 7];
id_post=[10 11 14];
id_tot=[3; 7; 11; 2; 6; 10; 14];
id=id_lat;

DATA = BigB;
% DATA = B3;

%add columns denoting which rows are from which mouse
mouse_id=[repmat(1,size(B1,1),1);repmat(2,size(B2,1),1);repmat(3,size(B3,1),1)];
% mouse_id=[repmat(1,size(B3,1),1)];
DATA=[DATA mouse_id];

%filter barcodes by max expression
BCodes = find(ismember(maxes,id));
% BCodes = find(ismember(m3,id));
data=DATA(BCodes,:);
data_norm=BigBnorm(BCodes,:);
% data_norm=B3norm(BCodes,:);
cutoff=max(data(:,1)); %calculate the cutoff value for these barcodes


%binarize input matrix
data_binary_tmp=double(data(:,projTargets)>cutoff);
data_binary=data_binary_tmp(sum(data_binary_tmp,2)~=0,:);
data_nonbinary=data(:,[projTargets size(data,2)]);
data_nonbinary_norm=data_norm(:,projTargets); %normalized data
cutoff_tmp=data_nonbinary(:,1:length(projTargets))<=cutoff; %cutoff based on non-norm data because of integer values
data_nonbinary_norm(cutoff_tmp)=0; %apply filter to normalize data
data_nonbinary_norm=data_nonbinary_norm(sum(data_binary_tmp,2)~=0,:); %deletes rows of all zeroes

% data_nonbinary(data_nonbinary(:,1:10)<cutoff)=0;
% data_nonbinary=data_nonbinary(sum(data_binary_tmp,2)~=0,:);

targetnum=size(data_binary,2);

%uniq
[data_binary_uniq,ia,ib]=unique(data_binary_tmp,'rows');
data_counts=accumarray(ib,1);

%calculate relative probability matrix of neurons projecting to A also
%projecting to B

prob_matrix = zeros(length(projTargets),length(projTargets));
prob_matrix_tot=zeros(length(projTargets),length(projTargets));
for pr = 1:length(projTargets)
    tmp_all=find(data_binary(:,pr));
    val1_all= data_binary(tmp_all,pr);
    for pr2 = 1:length(projTargets)
        val2_all= data_binary(tmp_all,pr2);
        diff_all=val2_all-val1_all; diff_all=sum(diff_all(diff_all<0));
        proj_fraction_all=(sum(val1_all)+diff_all)/sum(val1_all);
        proj_fraction_tot=(sum(val1_all)+diff_all)/length(data_binary);
        prob_matrix(pr2,pr) = proj_fraction_all;
        prob_matrix_tot(pr2,pr) = proj_fraction_tot;
    end
end


%count how many cells project to each area
sum_projections=sum(data_binary);
bifurc_p=0; %value of 0 means that all projections both to vLGN and dLGN are actually only to vLGN,
            %value of 1 means they are all true bifurcations
% sum_projections(1)=sum_projections(1)+sum_projections(3)*prob_matrix(1,3)*(1-bifurc_p)...
%     +sum_projections(4)*prob_matrix(1,4)*(1-bifurc_p);
% sum_projections(3)=sum_projections(3)+sum_projections(4)*prob_matrix(3,4)*(1-bifurc_p);

%get first order probabilities
neuronnum=length(BCodes);
p=sum_projections./neuronnum;
% update projection probability to LP to account for axons going to dLGN and vLGN
% p(1)=p(1)+p(3)*prob_matrix(1,3)*(1-bifurc_p)+p(4)*prob_matrix(1,4)*(1-bifurc_p);
% %update projection probability to dLGN to account for axons going to dvLGN
% p(3)=p(3)+p(4)*prob_matrix(3,4)*(1-bifurc_p);
% % %update projection probability to Pbg to account for axons going to MLR
% p(2)=p(2)+p(5)*prob_matrix(2,5)*(1-bifurc_p);


%construct motif matrix and associated expected probabilities from first
%order stats

%single proj
motif_tmp=zeros(targetnum,targetnum);
p_expected_tmp=zeros(targetnum,1);
counter=1;
for target1=1:targetnum
    motif_tmp(counter,target1)=1;
    p_expected_tmp(counter)=prod(p(logical(motif_tmp(counter,:))))*prod(1-p(~logical(motif_tmp(counter,:))));
    counter=counter+1;

end
%remove doublecounts
[motif_uno,ia,ib]=unique(motif_tmp(sum(motif_tmp,2)==1,:),'rows', 'stable');
p_expected_tmp2=p_expected_tmp(sum(motif_tmp,2)==1);
p_expected_tmp3_uno=p_expected_tmp2(ia);

%bifrucations
motif_tmp=zeros(targetnum^2,targetnum);
p_expected_tmp=zeros(targetnum^2,1);
counter=1;
for target1=1:targetnum
    for target2=1:targetnum
        motif_tmp(counter,target1)=1;
        motif_tmp(counter,target2)=1;
        p_expected_tmp(counter)=prod(p(logical(motif_tmp(counter,:))))*prod(1-p(~logical(motif_tmp(counter,:))));
        counter=counter+1;
    end
end
%remove doublecounts
[motif_bi,ia,ib]=unique(motif_tmp(sum(motif_tmp,2)==2,:),'rows', 'stable');
p_expected_tmp2=p_expected_tmp(sum(motif_tmp,2)==2);
p_expected_tmp3_bi=p_expected_tmp2(ia);

%trifrucations
motif_tmp=zeros(targetnum^3,targetnum);
p_expected_tmp=zeros(targetnum^3,1);
counter=1;
for target1=1:targetnum
    for target2=1:targetnum
        for target3=1:targetnum
            motif_tmp(counter,target1)=1;
            motif_tmp(counter,target2)=1;
            motif_tmp(counter,target3)=1;
            p_expected_tmp(counter)=prod(p(logical(motif_tmp(counter,:))))*prod(1-p(~logical(motif_tmp(counter,:))));
            counter=counter+1;
        end
    end
end
%remove doublecounts
[motif_tri,ia,ib]=unique(motif_tmp(sum(motif_tmp,2)==3,:),'rows', 'stable');
p_expected_tmp2=p_expected_tmp(sum(motif_tmp,2)==3);
p_expected_tmp3_tri=p_expected_tmp2(ia);

% quadfrucations
motif_tmp=zeros(targetnum^4,targetnum);
p_expected_tmp=zeros(targetnum^4,1);
counter=1;
for target1=1:targetnum
    for target2=1:targetnum
        for target3=1:targetnum
            for target4=1:targetnum
                motif_tmp(counter,target1)=1;
                motif_tmp(counter,target2)=1;
                motif_tmp(counter,target3)=1;
                motif_tmp(counter,target4)=1;
                p_expected_tmp(counter)=prod(p(logical(motif_tmp(counter,:))))*prod(1-p(~logical(motif_tmp(counter,:))));
                counter=counter+1;
            end
        end
    end
end
%remove doublecounts
[motif_quad,ia,ib]=unique(motif_tmp(sum(motif_tmp,2)==4,:),'rows', 'stable');
p_expected_tmp2=p_expected_tmp(sum(motif_tmp,2)==4);
p_expected_tmp3_quad=p_expected_tmp2(ia);

% 5 targets
motif_tmp=zeros(targetnum^5,targetnum);
p_expected_tmp=zeros(targetnum^5,1);
counter=1;
for target1=1:targetnum
    for target2=1:targetnum
        for target3=1:targetnum
            for target4=1:targetnum
                for target5=1:targetnum
                    motif_tmp(counter,target1)=1;
                    motif_tmp(counter,target2)=1;
                    motif_tmp(counter,target3)=1;
                    motif_tmp(counter,target4)=1;
                    motif_tmp(counter,target5)=1;
                    p_expected_tmp(counter)=prod(p(logical(motif_tmp(counter,:))))*prod(1-p(~logical(motif_tmp(counter,:))));
                    counter=counter+1;
                end
            end
        end
    end
end
%remove doublecounts
[motif_quint,ia,ib]=unique(motif_tmp(sum(motif_tmp,2)==5,:),'rows', 'stable');
p_expected_tmp2=p_expected_tmp(sum(motif_tmp,2)==5);
p_expected_tmp3_quint=p_expected_tmp2(ia);

% 6 targets
motif_tmp=zeros(targetnum^6,targetnum);
p_expected_tmp=zeros(targetnum^6,1);
counter=1;
for target1=1:targetnum
    for target2=1:targetnum
        for target3=1:targetnum
            for target4=1:targetnum
                for target5=1:targetnum
                    for target6=1:targetnum
                        motif_tmp(counter,target1)=1;
                        motif_tmp(counter,target2)=1;
                        motif_tmp(counter,target3)=1;
                        motif_tmp(counter,target4)=1;
                        motif_tmp(counter,target5)=1;
                        motif_tmp(counter,target6)=1;
                        p_expected_tmp(counter)=prod(p(logical(motif_tmp(counter,:))))*prod(1-p(~logical(motif_tmp(counter,:))));
                        counter=counter+1;
                    end
                end
            end
        end
    end
end
%remove doublecounts
[motif_sex,ia,ib]=unique(motif_tmp(sum(motif_tmp,2)==6,:),'rows', 'stable');
p_expected_tmp2=p_expected_tmp(sum(motif_tmp,2)==6);
p_expected_tmp3_sex=p_expected_tmp2(ia);


%combine all motifs and probabilities
motif=[motif_uno;motif_bi;motif_tri;motif_quad;motif_quint;motif_sex];
p_all=[p_expected_tmp3_uno;p_expected_tmp3_bi;p_expected_tmp3_tri;...
    p_expected_tmp3_quad;p_expected_tmp3_quint;p_expected_tmp3_sex];


%add projections to 0 areas
motif=[zeros(1,length(projTargets));motif];
p_all=[1-sum(p_all);p_all];

%count up the number of observed
observed=zeros(size(p_all));

[d,loc]=ismember(motif,data_binary_uniq,'rows');
observed(d)=data_counts(loc(loc~=0));

%calculate expected
expected=neuronnum.*p_all;

%calculate pvalues
% pvalue=binocdf(observed,bifrucnum,p_expected);
%left-tail
pvalue_left=binocdf(observed,ceil(neuronnum),p_all);
%right-tail
pvalue_right=1 - binocdf(observed - 1, ceil(neuronnum), p_all);

%apply filtering before computing significance
observed_filtered=observed(motif_start+1:end);
expected_filtered=expected(motif_start+1:end);
expected_prob_filtered=p_all(motif_start+1:end);
observed_prob=observed./sum(observed);
observed_prob_filtered=observed_prob(motif_start+1:end);
motif_filtered=motif(motif_start+1:end,:);
obs_keep=observed_filtered>=2; exp_keep=expected_filtered>=2; keep=obs_keep+exp_keep;
observed_nnz=observed_filtered(keep>0);
expected_nnz=expected_filtered(keep>0);
observed_prob_nnz=observed_prob_filtered(keep>0);
expected_prob_nnz=expected_prob_filtered(keep>0);
diff=abs(observed_nnz-expected_nnz); diff_id=diff>=2;
motifs=motif_filtered(keep>0,:);


pvalue_left=pvalue_left(motif_start+1:end);
pvalue_right=pvalue_right(motif_start+1:end);
pvalue_left=pvalue_left(keep>0);
pvalue_right=pvalue_right(keep>0);
pvalue_all=[pvalue_left;pvalue_right];

%decide on significance using bonferoni correction
% alpha=0.05;
% motifnum=size(motifs,1);
% alpha_adj=alpha/motifnum;
% overorunder=zeros(size(pvalue_left));
% overorunder(pvalue_left<(alpha_adj/2))=-1;
% overorunder(pvalue_right<(alpha_adj/2))=1;
% significance=overorunder~=0;

% %decide on significance using holm-bonferoni correction
% m = length(pvalue_all); % Total tests (both tails)
% [sorted_p, indices] = sort(pvalue_all);
% alpha_holm = 0.05 ./ (m:-1:1);
% reject = sorted_p < alpha_holm';
% rejected_hypotheses = false(1, m);
% rejected_hypotheses(indices) = reject;
% significance=(rejected_hypotheses(1:m/2)+rejected_hypotheses(m/2+1:end))';

% %decide on significance using benjamini-hochberg correction
m = length(pvalue_all);
[sorted_p, indices] = sort(pvalue_all);
alpha_bh = (1:m) / m * 0.05;
reject = sorted_p < alpha_bh';
rejected_hypotheses = false(1, m);
rejected_hypotheses(indices) = reject;
significance=(rejected_hypotheses(1:m/2)+rejected_hypotheses(m/2+1:end))';


sig_nnz=significance.*(diff_id);

%
% [padj_left,alpha_adj_left] = multicmp (pvalue_left,'fdr');% Holm-Bonferroni (type: 'down'), Hochberg-Bonferroni (type: 'up'), or Benjamini-Hochberg (type: 'fdr') correction
% [padj_right,alpha_adj_right] = multicmp (pvalue_right,'fdr');% Holm-Bonferroni (type: 'down'), Hochberg-Bonferroni (type: 'up'), or Benjamini-Hochberg (type: 'fdr') correction
% overorunder=zeros(size(p_all));
% overorunder(padj_left<0.05)=-1;
% overorunder(padj_right<0.05)=+1;
% significance=overorunder~=0;

% %----make plots
% %filter based on which motifs you want to show
% observed_filtered=observed(motif_start+1:end);
% expected_filtered=expected(motif_start+1:end);
% motif_filtered=motif(motif_start+1:end,:);
% significance_filtered=significance(motif_start+1:end);
%
% obs_keep=observed_filtered>0; exp_keep=expected_filtered>3; keep=obs_keep+exp_keep;
% observed_nnz=observed_filtered(keep>0);
% expected_nnz=expected_filtered(keep>0);
% diff=abs(observed_nnz-expected_nnz); diff_id=diff>=2;
% sig_nnz=significance_filtered(keep>0);
% sig_nnz=sig_nnz.*(diff_id);
% motifs=motif_filtered(keep>0,:);

vert_motifs=flipud(rot90(motifs,1));


%plot observed vs expected
figure(1); set(gcf, 'Position', [445.8000 145 827.2000 480])
til=tiledlayout(3,6, 'TileSpacing', 'tight', 'Padding','compact');
set(til, 'OuterPosition', [0 0 0.98 0.98]);
%set conditions for medial and lateral
if isequal(id,id_med)
    new_cols=[0.7 0.7 0.7; col(3,:)];
    til.Title.String='Upper sSC projection motifs';
    col_motifs=cbrewer('seq', 'YlOrBr', 30);
    col_motifs=flipud(col_motifs(7:end,:));
elseif isequal(id,id_lat)

    new_cols=[0.7 0.7 0.7; col(4,:)];
    til.Title.String='Lower sSC projection motifs';
    col_motifs=cbrewer('seq', 'PuBu', 30);
    col_motifs=flipud(col_motifs(15:end,:));
else
    new_cols=[col_bias_med_merge(3,:); col_bias_lat_merge(3,:)];
    til.Title.String='sSC projection motifs';
    col_motifs=cbrewer('seq', 'Greys', 30);
    col_motifs=flipud(col_motifs(7:end,:));
end


til.Title.FontSize=13;
til.Title.FontWeight='bold';
%plot horizontal density curve
% log_obs=log2(observed_nnz); log_exp=log2(expected_nnz);

ax1=nexttile(3,[2 4]);

% diffVals = max(expected_nnz - test, 0);
% dataToPlot = [test diffVals];
% xpos=1:length(observed_nnz);
% new_cols= [new_cols(1, :);
%      0.4 0.4 0.4;
%      new_cols(2, :)];

% bar([log_exp, log_obs], 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 0.9);
% bar(xpos-0.15,dataToPlot,'stacked', 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 1, 'BarWidth', 0.25); hold on
% bar(xpos+0.15,observed_nnz, 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 1, 'BarWidth',0.25);
bar([expected_nnz,observed_nnz], 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 1);
colororder(new_cols)
xlim([0.5 size(vert_motifs,2)+0.5])

% ylim([0 10]);
ylim([0 max([expected_nnz;observed_nnz])+3])
hold on

% title('observed and expected counts of polyfrucations')
% legend('expected','observed');


%plot sig star
for i=1:length(observed_nnz)
    xvector=i;
    % yl=max([log_exp(i) log_obs(i)])+0.3;
    yl=max([expected_nnz(i) observed_nnz(i)])+2;
    % sigstars on top
    mysigstar(gca, xvector, yl, sig_nnz(i),'k', 'down', 'singular');
end
set(ax1, 'XMinorTick', 'off', 'box', 'off', 'xticklabels', [], 'TickDir','out', ...
    'FontSize', 11)
ax1.XTick=[];
ylabel('motif count', 'FontSize', 12)
% ylabel('log2 motif count', 'FontSize', 12)
%number of cells projecting to each area

col_motifs=[0.9 0.9 0.9;col_motifs]; %make 0 values value light gray

ax1_5=nexttile(13, [1 2]);
barh(1:length(projTargets),sum_projections, 'FaceColor', col_motifs(4,:), 'EdgeAlpha',0);
set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'box', 'off', 'YGrid', 'off', 'TickDir', 'out',...
    'color', 'none', 'YColor', 'none','XTick', 50:50:max(sum_projections));
ax1_5.XAxis.FontSize=9;
xlabel(["No. of cells projec-";"ting to each region"], 'FontSize',10)
ylim([0.5 length(projTargets)+0.5]);
xlim([0 max(sum_projections)])

%motif patterns
ax2=nexttile(15,[1 4]);
% Get the size of the matrix
[m, n] = size(vert_motifs);
% Generate the x and y coordinates for each bubble
[X, Y] = meshgrid(1:n, 1:m);
% Flatten the coordinates and the data for scatter plotting
X = X(:);
Y = Y(:);
%replace values by sum of columns to categorize mono,bi,tri furcations
colSums = repmat(sum(vert_motifs, 1),m,1);
vert_motifs2=vert_motifs;
vert_motifs2(vert_motifs==1)=colSums(vert_motifs == 1);
bubbleColors = vert_motifs2(:);  % Flatten the data values for coloring
max_furc=max(bubbleColors);
%Plot the bubbles
sizes=25+15*vert_motifs(:);
scatter(X, Y, sizes, bubbleColors, 'filled');
% col_motifs=col_motifs([1,end-max_furc+1:end],:);
colormap(col_motifs([1,[10,15]],:))

% imagesc(vert_motifs);
% inv_gray = flipud(grays);
% inv_gray=inv_gray(1:200,:);
% colormap(inv_gray)
set(ax2, 'XMinorTick', 'off', 'YMinorTick', 'off', 'box', 'on', 'xticklabels', [], ...
    'yticklabels', sample_names(projTargets), 'ytick', 1:length(projTargets),...
    'xtick', 1:size(vert_motifs,2), 'ticklength', [0,0], 'YDir', 'reverse',...
    'FontSize', 11)
ax2.XTick=[]; ax2.YTick=1:17; ax2.TickLength=[0 0];
xlabel('Motifs', 'FontSize', 13)

ax2.XAxis.Color=[1 1 1];
ax2.XAxis.TickLabelColor=[0 0 0];
ax2.XAxis.Label.Color=[0 0 0];
% ax2.YAxis.Color=[1 1 1];
ax2.YAxis.TickLabelColor=[0 0 0];
ax2.YAxis.Label.Color=[0 0 0];
xtickangle(0);
ylim([0.5 length(projTargets)+0.5])
xlim([0.5 size(vert_motifs,2)+0.5])


%motifs names
names=sample_names(projTargets);
motifs_names={};

% Iterate through each column
for i = 1:size(vert_motifs, 2)
    % Find the indices of rows that contain 1 in the current column
    indices = find(vert_motifs(:, i) == 1);

    % Get the names corresponding to these indices
    selected_names = names(indices);

    % concatenate the names
    cat_names=strjoin(selected_names, ' - ');

    %place the motif name into the motifs_names variable
    motifs_names{i}=cat_names;
end

observed_prb=observed_nnz./neuronnum; %for probabilities instead of counts
expected_prb=expected_nnz./neuronnum; %for probabilities instead of counts

%% calculate cumulative distribution and difference of medial and lateral sSC motifs
if isequal(id,id_med)
    motifs_med=motifs;
    motifs_names_med=motifs_names;
    observed_med=observed_prb;
    observed_med_counts=observed_nnz;
    expected_med=expected_prb;
    sig_med=sig_nnz;
    neurons_med=neuronnum;
elseif isequal(id,id_lat)
    motifs_lat=motifs;
    motifs_names_lat=motifs_names;
    observed_lat=observed_prb;
    observed_lat_counts=observed_nnz;
    expected_lat=expected_prb;
    sig_lat=sig_nnz;
    neurons_lat=neuronnum;
elseif isequal(id,id_ant)
    motifs_med=motifs;
    motifs_names_med=motifs_names;
    observed_med=observed_prb;
    observed_med_counts=observed_nnz;
    expected_med=expected_prb;
    sig_med=sig_nnz;
    neurons_med=neuronnum;
elseif isequal(id,id_post)
    motifs_lat=motifs;
    motifs_names_lat=motifs_names;
    observed_lat=observed_prb;
    observed_lat_counts=observed_nnz;
    expected_lat=expected_prb;
    sig_lat=sig_nnz;
    neurons_lat=neuronnum;
end

col1=col(3,:);
col2=col(4,:);

combined_motifs=[motifs_med;motifs_lat]; %merge motifs of medial and lateral
[unique_motifs, ~, unique_indices] = unique(combined_motifs, 'rows', 'stable'); %keep only unique ones

num_motifs = size(unique_motifs, 1);
shared_observed = zeros(num_motifs, 2); % Column 1: med, Column 2: lat
shared_observed_counts = zeros(num_motifs, 2); % Column 1: med, Column 2: lat
shared_expected = NaN(num_motifs, 2); % Column 1: med, Column 2: lat
shared_sig = NaN(num_motifs,2);
shared_motifs_names = cell(num_motifs, 1); % Store motif names

for i = 1:num_motifs
    % Check if motif exists in A
    idx_A = find(ismember(motifs_med, unique_motifs(i, :), 'rows'), 1);
    if ~isempty(idx_A)
        shared_observed(i,1) = observed_med(idx_A); % Assign to A column
        shared_observed_counts(i,1) = observed_med_counts(idx_A);
        shared_expected(i,1) = expected_med(idx_A);
        shared_sig(i,1) = sig_med(idx_A);
        shared_motifs_names{i} = motifs_names_med{idx_A}; % Keep motif name
    end

    % Check if motif exists in B
    idx_B = find(ismember(motifs_lat, unique_motifs(i, :), 'rows'), 1);
    if ~isempty(idx_B)
        shared_observed(i,2) = observed_lat(idx_B); % Assign to B column
        shared_observed_counts(i,2) = observed_lat_counts(idx_B);
        shared_expected(i,2) = expected_lat(idx_B);
        shared_sig(i,2) = sig_lat(idx_B);
        if isempty(shared_motifs_names{i})
            shared_motifs_names{i} = motifs_names_lat{idx_B};
        end
    end
end

% -------compute empirical cumulative distribution function

% Compute differences
diff_medial = observed_med - expected_med;
diff_lateral = observed_lat - expected_lat;

%check difference in distribution
[h, p, ksstat] = kstest2(diff_medial, diff_lateral); %K-S test


% [h, p, ci, stats] = vartest2(diff_medial, diff_lateral); %F-test


[FA, XA] = ecdf(diff_medial);
[FB, XB] = ecdf(diff_lateral);

figure; set(gcf, 'position', [488 287.4000 374.6000 374.6000]);
hold on;
stairs(XB, FB, 'Color', col2, 'LineWidth', 2.5);
stairs(XA, FA, 'Color', col1, 'LineWidth', 2.5);

xlabel('Observed - Expected');
ylabel('Cumulative Distribution');
% legend({'medial sSC', 'lateral sSC'});
% title('Cumulative Distribution of Observed - Expected Motif Prob');
xlim([-0.015 0.015])
% xlim([-90 90])
xline(0,'--','Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
set(gca, 'TickDir', 'out', 'xtick', [-0.02 -0.01 0 0.01 0.02], 'ytick', [0:0.2:1],...
    'xticklabels',{'-0.02', '-0.01', '0', '0.01','0.02'}, 'yticklabels',...
    {'0','','','','','1'},'FontSize', 11, 'LineWidth', 1);
% set(gca, 'TickDir', 'out', 'ytick', [0:0.2:1],'yticklabels',...
%     {'0','','','','','1'},'FontSize', 11, 'LineWidth', 1);
hold off;


%-------compare deviation from binomial model for medial and lateral sSC
shared_deviation=(shared_observed-shared_expected)./shared_expected;
lim=5; %limit of value
shared_deviation(~isnan(shared_deviation)) = max(min(shared_deviation(~isnan(shared_deviation)), lim), -lim);
x_pos=1:size(shared_deviation,1);
offset=0.1;
alpha_sig=shared_sig; alpha_sig(alpha_sig==0)=0.4; alpha_sig(isnan(alpha_sig))=0;

figure; set(gcf, 'position', [488 41.8000 374.6000 740.8000]);
hold on;
for i=1:size(shared_deviation,1)
    scatter(shared_deviation(i, 1),x_pos(i)-offset, 35, 'filled',...
        'MarkerFaceColor',col1, 'MarkerFaceAlpha', alpha_sig(i,1));
    scatter(shared_deviation(i, 2),x_pos(i)+offset, 35, 'filled',...
        'MarkerFaceColor',col2, 'MarkerFaceAlpha', alpha_sig(i,2));
end

xlabel('(Observed - Expected)/Expected');
% ylabel('Motif');
% title('Cumulative Distribution of Observed - Expected Motif Prob');
xlim([-lim-1 lim+1])
ylim([0.5 size(shared_deviation,1)+0.5])
xline(0,'--','Color', [0.6 0.6 0.6], 'LineWidth', 1);
set(gca, 'TickDir', 'out', 'ytick', [1:size(shared_deviation,1)],...
    'yticklabels', shared_motifs_names,'FontSize', 11, 'LineWidth', 1, 'YDir', 'reverse');
% set(gca, 'TickDir', 'out', 'ytick', [0:0.2:1],'yticklabels',...
%     {'0','','','','','1'},'FontSize', 11, 'LineWidth', 1);
% legend({'upper sSC', 'lower sSC'}, 'Location', 'southeast');
hold off;


%-----statistical analysis on probability or proportion of motifs for med
%and lat
p_values_medlat = zeros(num_motifs,1);
for i = 1:size(shared_observed_counts,1)
    x1 = shared_observed_counts(i,1);
    x2 = shared_observed_counts(i,2);

    table = [x1, neurons_med - x1;
             x2, neurons_lat - x2];

    [~, p] = fishertest(table);  % Two-sided Fisher’s exact test
    p_values_medlat(i) = p;
end

m = length(p_values_medlat);

% Bonferroni correction
p_adj_bonf = min(p_values_medlat * m, 1);  % cap at 1

% FDR correction (Benjamini-Hochberg
[p_adj_fdr,~] = multicmp (p_values_medlat,'fdr');


%plot probabilities of projection motifs
figure; set(gcf, 'position', [873.8000 41.8000 375.2000 740.8000]);
hold on
% for i=1:size(shared_observed,1)
%     scatter(shared_observed(i, 1),x_pos(i)-offset, 35, 'filled',...
%         'MarkerFaceColor',col(3,:), 'MarkerFaceAlpha', 1);
%     scatter(shared_observed(i, 2),x_pos(i)+offset, 35, 'filled',...
%         'MarkerFaceColor',col(4,:), 'MarkerFaceAlpha', 1);
% end


% barh(x_pos-offset,shared_observed(:, 1), 'BarWidth', 0.17,...
%     'EdgeColor', 'none', 'FaceColor', col1, 'FaceAlpha', 1);
% barh(x_pos+offset,shared_observed(:, 2), 'BarWidth', 0.17,...
%         'EdgeColor', 'none', 'FaceColor', col2, 'FaceAlpha', 1);

b1=barh(x_pos-offset,shared_observed(:, 1), 'BarWidth', 0.17,...
    'EdgeColor', 'none', 'FaceColor', col1, 'FaceAlpha', 1);
b2=barh(x_pos+offset,shared_observed(:, 2), 'BarWidth', 0.17,...
        'EdgeColor', 'none', 'FaceColor', col2, 'FaceAlpha', 1);

% xlabel('Observed motif prob');
xlabel('Observed motif prob');
% ylabel('Motif');
% title('Cumulative Distribution of Observed - Expected Motif Prob');
% xlim([0 max(shared_observed,[],'all')+0.001])
% ylim([0.5 size(shared_deviation,1)+0.5])
xlim([0 max(shared_observed,[],'all')+0.005])
ylim([0.5 size(shared_deviation,1)+0.5])
% xline(0,'--','Color', [0.6 0.6 0.6], 'LineWidth', 1);
% set(gca, 'TickDir', 'out', 'ytick', [1:size(shared_deviation,1)],...
%     'yticklabels', shared_motifs_names,'FontSize', 11, 'LineWidth', 1, 'YDir', 'reverse');
set(gca, 'TickDir', 'out', 'ytick', [1:size(shared_deviation,1)],...
    'yticklabels', shared_motifs_names,'FontSize', 11, 'LineWidth', 1, 'YDir', 'reverse');
% set(gca, 'TickDir', 'out', 'ytick', [0:0.2:1],'yticklabels',...
%     {'0','','','','','1'},'FontSize', 11, 'LineWidth', 1);
% legend({'upper sSC', 'lower sSC'}, 'Location', 'southeast');


for jj=1:length(x_pos)
    xvector=x_pos(jj);
    yl=max([shared_observed(jj,1) shared_observed(jj,2)])+0.003;
    % sigstars on top
    mysigstar(gca, yl,xvector,  p_adj_fdr(jj));
end


hold off;


%% calculate dropfraction - add or subtract random projections to/from data

observed_motif_counts = observed;  % Replace with your actual observed motif counts
N = sum(observed_motif_counts);  % Total number of neurons
n_motifs = length(observed_motif_counts);  % Number of motifs (2^n_areas)
n_areas = length(p);

% Binomial model (dropfraction = 1)
expected_counts_binomial = zeros(1, n_motifs);

% Compute expected counts under binomial model
for i = 1:n_motifs
    % Convert motif index to binary vector
    motif_mask = motif(i,:);
    % Compute probability of this motif
    motif_prob = prod(p.^motif_mask .* (1 - p).^(1 - motif_mask));
    expected_counts_binomial(i) = N * motif_prob;
end

% Define a range of dropfractions
dropfractions = 0.1:0.1:2.0; % False negatives (<1), no change (=1), false positives (>1)

% Number of repetitions
n_repeats = 400;

% Initialize matrices for storing results
p_values = zeros(length(dropfractions), n_motifs, n_repeats);

% Initialize matrix to store modified motif counts
modified_counts = zeros(length(dropfractions), n_motifs,n_repeats);

% Loop over dropfractions
for df_idx = 1:length(dropfractions)
    dropfraction = dropfractions(df_idx);

    for rep = 1:n_repeats
        % Start with observed counts
        adjusted_counts = observed_motif_counts;

        % Calculate relative frequencies of motifs based on projection probabilities
        relative_probs = expected_counts_binomial / sum(expected_counts_binomial); % Normalized

        for i = 2:n_motifs  % Skip the all-zero motif (index 1)
            if dropfraction < 1
                % Simulate false negatives: Reduce projections proportionally
                lost_counts = binornd(observed_motif_counts(i), 1 - dropfraction);
                adjusted_counts(i) = adjusted_counts(i) - lost_counts;
            elseif dropfraction > 1
                % Simulate false positives: Add projections proportionally
                additional_counts = binornd(observed_motif_counts(i), dropfraction - 1);
                adjusted_counts(i) = adjusted_counts(i) + additional_counts;
            end
        end


        % Adjust the all-zero motif (index 1)
        if dropfraction < 1
            % Add neurons lost from other motifs
            lost_projections = sum(observed_motif_counts(2:end)) - sum(adjusted_counts(2:end));
            adjusted_counts(1) = observed_motif_counts(1) + lost_projections;
        elseif dropfraction > 1
            % Redistribute false positives among all motifs (including all-zero motif)
            gained_projections = sum(adjusted_counts(2:end)) - sum(observed_motif_counts(2:end));
            for i = 2:n_motifs
                gain_fraction = relative_probs(i) / sum(relative_probs(2:end)); % Proportional allocation
                adjusted_counts(i) = adjusted_counts(i) + round(gain_fraction * gained_projections);
            end
            adjusted_counts(1) = max(0, observed_motif_counts(1) - sum(adjusted_counts(2:end)) + sum(observed_motif_counts(2:end)));
        end

        % Store the adjusted motif counts
        modified_counts(df_idx, :,rep) = adjusted_counts;
        % Compute p-value for each motif using binomial test
        p_motifs = expected_counts_binomial / N;
        % Compute left and right tail probabilities
        p_left = binocdf(adjusted_counts', N, p_motifs);  % Left tail (probability of k or fewer)
        p_right = 1 - binocdf(adjusted_counts' - 1, N, p_motifs);  % Right tail (probability of k or more)
        % Two-tailed p-value: multiply by 2 and take the minimum
        p_values(df_idx, :, rep) = 2 * min(p_left, p_right);

    end
end

% Replace NaN values with 0
p_values(isnan(p_values)) = 0;
p_values(p_values>=1)=1;
% Compute mean and confidence intervals
mean_p_values = mean(p_values, 3);  % Mean across repetitions
conf_int_low = quantile(p_values, 0.025, 3);  % Lower bound (2.5th percentile)
conf_int_high = quantile(p_values, 0.975, 3);  % Upper bound (97.5th percentile)

% Plot results

for i = 1:n_motifs
    figure(i)

    % Plot mean p-value as a line
    plot(dropfractions, mean_p_values(:, i), 'b-', 'LineWidth', 1.5);
    hold on;

    % Plot shaded confidence interval
    fill([dropfractions, fliplr(dropfractions)], ...
         [conf_int_low(:, i)', fliplr(conf_int_high(:, i)')], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    xlabel('Dropfraction');
    ylabel('P-Value');
    title(sprintf('Motif %d', i));  % Motif index (1 = All-Zero)
    ylim([0 1]);  % P-values range from 0 to 1
    grid on;
end

% Overall figure title
sgtitle('P-Value Trends with Confidence Intervals for Each Motif');


modified_counts_filtered=modified_counts(:,motif_start+1:end);
expected_counts_binomial_filtered=expected_counts_binomial(motif_start+1:end);

N = neuronnum;  % Example: 1000 neurons

% Probabilities of projecting to each area
p_1 = p(1);  % Probability of projecting to Area 1
p_2 = p(2);  % Probability of projecting to Area 2
p_3 = p(3);  % Probability of projecting to Area 3
p_4 = p(4);  % Probability of projecting to Area 4
p_5 = p(5);  % Probability of projecting to Area 5

% % Dependency factor for Area 1 given Area 2
dependency_factor = 1;  % 30% dependency

% Initialize the projection probability matrix
projections = zeros(N, 5);

% Compute probabilities for Area 2
% projections(:, 1) = p_1;
projections(:, 4) = p_4;  % Independent probability for Area 2

% Compute probabilities for Area 1 with dependency on Area 2
% Base probability for neurons not projecting to Area 2
projections(:, 3) = p_3 * (1 - projections(:, 4));
% Adjust probability for neurons projecting to Area 2
projections(:, 3) = projections(:, 3) + projections(:, 4) * dependency_factor;

% Compute probabilities for other areas (independent)
projections(:, 1) = p_1;  % Probability for Area 3
projections(:, 2) = p_2;  % Probability for Area 4
projections(:, 5) = p_5;  % Probability for Area 5

% Motif analysis based on probabilities
motifs = [];
for i = 0:31  % 2^5 = 32 possible motifs
    % Convert motif index to binary vector
    motif_mask = bitget(i, 5:-1:1);  % Binary vector of length 5
    % Compute probability of this motif by multiplying probabilities
    motif_prob = prod(projections.^repmat(motif_mask, N, 1) .* (1 - projections).^(1 - repmat(motif_mask, N, 1)), 2);
    motif_count = sum(motif_prob);  % Total expected number of neurons with this motif
    % Store motif and count
    motifs = [motifs; motif_mask, motif_count];
end


%% piechart showing the proportion of mono, bi, trifurcating neurons

%how many different motifs
mtf_num=sum(motif,2);
obs_num=observed(observed>0); mtf_num=mtf_num(observed>0);

% logical indexing
is_one=(mtf_num == 1);
is_two=(mtf_num == 2);
% is_three=(mtf_num == 3);
is_greater_than_three=(mtf_num >= 3);

%sum counts
total_counts_one = sum(obs_num(is_one));
total_counts_two = sum(obs_num(is_two));
% total_counts_three = sum(obs_num(is_three));
total_counts_greater_than_three = sum(obs_num(is_greater_than_three));
% total_counts=[total_counts_one, total_counts_two, total_counts_three, total_counts_greater_than_three];
% total_counts=[total_counts_one+16, total_counts_two-16];
total_counts=[total_counts_one, total_counts_two, total_counts_greater_than_three];
% total_counts_percent

figure(1); set(gcf, 'position', [488 349.8000 416.2000 312.2000]);
% labels=string({"1 target";"2 targets";"3 targets";">3 targets"});
labels=string({"1 target";"2 targets";">=3 targets";});
p=piechart(total_counts, labels);
p.Direction = "counterclockwise";
p.EdgeColor = "white"; p.LineWidth = 1;
p.FaceAlpha=0.7;
p.StartAngle=90;
p.FontSize=13;
% colors=cbrewer('qual','Set3',7);
colors=flipud(cbrewer('seq','Oranges',4));
colors=colors(1:1:end,:);
colororder(colors)


%% volcano plot

% alpha_adj=alpha/motifnum;
% overorunder=zeros(size(p_all));
% overorunder(pvalue<(alpha_adj/2))=-1;
% overorunder(pvalue>(1-alpha_adj/2))=1;
%
% pval_min=nan*ones(size(pvalue));
% pval_min(pvalue<0.5)=pvalue(pvalue<0.5);
% pval_min(pvalue>0.5)=1-pvalue(pvalue>0.5);
% pval_min_filtered=pval_min(motif_start+1:end);
% mavolcanoplot(observed,expected,pval_min,'LogTrans',1,'PlotOnly',1,...
%      'Foldchange',2,'PCutoff',alpha_adj/2)
%-------do your own volcano plot

%datasets of ratios and pvalues
num_points=length(observed_nnz);

ratios=observed_nnz./expected_nnz;
% pvalues=pval_min_filtered(keep>0);
% pcutoff=alpha_adj/2; %p-value cutoff, comment if you use benjamini hochberg
pvalues=[pvalue_all(1:length(pvalue_all)/2) pvalue_all(length(pvalue_all)/2+1:end)];
pvalues=min(pvalues,[],2);
% pcutoff=zeros(length(observed_nnz)*2,1);
% pcutoff(indices)=alpha_bh; %p-value cutoff
% pcutoff=[pcutoff(1:length(pcutoff)/2) pcutoff(length(pcutoff)/2+1:end)];
% pcutoff=min(pcutoff,[],2);
% pcutoff=pcutoff.*sig_nnz;

%convert values to log
logpvalues= -log10(pvalues);
logratios= log2(ratios);

%replace -inf and inf value with a finite number
% logratios(logratios==-inf)=floor(min(logratios(logratios~=-inf)));
logratios(logratios==-inf)=-6;
logratios(logratios==inf)=ceil(max(logratios(logratios~=inf)));

% Determine point sizes and colors
pointSizes = ones(num_points, 1) * 15; % Default size
pointColors = repmat([0.6 0.6 0.6], num_points, 1); % Default color

% Increase size and change color for points above cutoffs
for i = 1:num_points
    if sig_nnz(i) && logratios(i) > 0 %&& sig_nnz(i) %overrepresented
        pointSizes(i) = 20; % Increase size
        % pointColors(i, :) = [0.384 0.729 0.098]; % Change color
        pointColors(i, :) = col(4,:);
    elseif sig_nnz(i) && logratios(i) < 0 %&& sig_nnz(i)%underrepresented
        pointSizes(i) = 20; % Increase size
        % pointColors(i, :) = [0.839 0.443 0.165]; % Change color
        pointColors(i, :) = col(4,:);

    else
        continue
    end
end

% Create the volcano plot
figure(2); set(gcf, 'Position', [420.2000 401.8000 347.2000 260.2000])
scatter(logratios, logpvalues, pointSizes, pointColors, 'filled');
xlabel(["Effect size";"\fontsize{11}log_{2}(observed/expected)"]);
% ylabel(["Significance";"\fontsize{9}-log_{10}({\itP} value)"]);
ylabel(["Significance";"\fontsize{11}{\itP} value (log_{10} scale)"]);
set(gca, 'Color', 'none')

% Add a vertical line at fold change cutoff
hold on;
line([0 0], ylim, 'Color', [0.7 0.7 0.7], 'LineStyle', '--', 'LineWidth', 1);

% Add a horizontal line at p-value cutoff
% xlim([-6 6]);
xlim([-3.5 3.5]);
xlims=xlim;
% line(xlims, [-log10(pcutoff) -log10(pcutoff)], 'Color', [0.7 0.7 0.7], 'LineStyle', '--',...
%     'LineWidth', 1);
% text(xlims(1)+0.4, -log10(pcutoff)-0.5, '{\itP} value cutoff', 'Fontsize', 9); %add text for p-value cutoff
%label p-values without logarithmic scale
ytickvalues=get(gca, 'Ytick');
newyticklabels=arrayfun(@(v) sprintf('-10^{%g}', v), ytickvalues, 'UniformOutput', false);
newyticklabels{1}=0; newyticklabels{2}=''; newyticklabels{4}=''; newyticklabels{6}='';
set(gca, 'YTickLabel', newyticklabels);

% Add text labels to significant points
for i = 1:num_points
    if sig_nnz(i) && logratios(i) < 0 %&& sig_nnz(i)
        text(logratios(i), logpvalues(i) +0.1, motifs_names{i}, 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    elseif sig_nnz(i) && logratios(i) > 0 %&& sig_nnz(i)
        text(logratios(i), logpvalues(i) +0.1, motifs_names{i}, 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end

%add underrepresented and overrepresented text
ylims=ylim;
% text(xlims(1)/2, ylims(2)+0.5, 'Underrepresented', 'Color', col(4,:), 'Fontsize', 12,...
%     'HorizontalAlignment', 'center');
% text(xlims(2)/2, ylims(2)+0.5, 'Overrepresented', 'Color', col(4,:), 'Fontsize', 12,...
%      'HorizontalAlignment', 'center');


set(gca, 'TickDir', 'out')

hold off;

%% projection motif analysis
id_under=sig_nnz & logratios<0; %underreppresented
id_over=sig_nnz & logratios>0; %overrepresented
exp=expected_nnz(id_under); exp=[exp;expected_nnz(id_over)];
sig_motifs=motifs(id_under,:); sig_names=motifs_names(id_under);
sig_motifs=[sig_motifs; motifs(id_over,:)]; sig_names=[sig_names motifs_names(id_over)];
id_sig=sum(sig_motifs,2)>1;
sig_motifs=sig_motifs(id_sig,:); sig_names=sig_names(id_sig); exp=exp(id_sig);

sig_motifs=motifs; sig_names=motifs_names;
exp=expected_nnz;

% plot all the motifs

figure(1); set(gcf, 'Position', [299.4000 131.4000 643.2000 457.6000])
til=tiledlayout('flow', 'TileSpacing', 'loose', 'Padding','compact');
set(til, 'OuterPosition', [0 0 0.98 0.98]);
m_cols=colormap('lines');
for t=1:size(sig_motifs,1)
    motif_num=t;
    mot=logical(sig_motifs(motif_num,:));
    comp_mot=~mot; %complement array to mask columns that should be 0
    dt=data_nonbinary_norm(:,mot); dt_comp=data_nonbinary_norm(:,comp_mot);
    row_index = all(dt ~= 0, 2) & all(dt_comp==0,2); %find rows that match motif
    filtered_dt=data_nonbinary_norm(row_index,:); %extract data that matches motif
    % mouse_col=filtered_dt(:,11); %extract mouse id column
    norm_vals=max(filtered_dt(:,1:length(projTargets)),[],2); %maximum value of each row for normalization
    %add random small value to norm vals
    epsilon = 10^(-3); % Small number to ensure changes are minimal
    random_offsets = epsilon * (2 * rand(size(norm_vals)) - 1); % Random values between -epsilon and +epsilon
    norm_vals2 = norm_vals + random_offsets; % Add the small random values to the original values

    norm_dt=filtered_dt(:,1:length(projTargets))./norm_vals2; %normalize data by maximum value of row
    if ~isempty(norm_dt)
        nexttile;
    end


    % %make colormap based on mouse
    % m_cols=[];
    % col1='#6A2C70'; col1=hex2rgb(col1);
    % col2='#B83B5E'; col2=hex2rgb(col2);
    % col3='#FFAF45'; col3=hex2rgb(col3);
    % for m=1:length(mouse_col)
    %     if mouse_col(m)==1
    %         m_cols=[m_cols; col1];
    %     elseif mouse_col(m)==2
    %         m_cols=[m_cols; col2];
    %     else
    %         m_cols=[m_cols; col3];
    %     end
    % end

    %plot
    for i=1:size(norm_dt,1)
        plot(norm_dt(i,:), 'Marker', 'o','LineWidth',1, 'Color',m_cols(i,:), 'MarkerSize',5);
        hold on
    end
    set(gca, 'box', 'off', 'color', 'none', 'xticklabel', names, 'YTick', 0:0.2:1, ...
        'TickDir', 'out', 'Xtick', 1:length(projTargets), 'XTickLabelRotation', 90);
    xlim([0.5 length(projTargets)+0.5]); ylim([-0.03 1.05]); xlims=xlim; ylims=ylim;
    title(sig_names(motif_num), 'FontWeight', 'normal', 'FontSize', 10)
    ylabel('Projection strength', 'FontSize', 9)
    text(xlims(2)-3, ylims(2)-0.2, [strcat("Observed: ", string(size(norm_dt,1)));...
        strcat("Expected: ", string(round(exp(t))))], 'Fontsize', 9);


end


%% Local functions
function [xOut, yOut] = get_curves2(x1, y1, x2, y2)
% Smooth (cosine) curves from (x1, y1) to (x2, y2), one column per curve
    Npoints = 15;
    t = linspace(0, pi, Npoints);
    c = (1-cos(t))./2;
    Ncurves = numel(y1);
    yOut = repmat(y1, Npoints, 1) + repmat(y2 - y1, Npoints, 1) .* repmat(c', 1, Ncurves);
    xOut = repmat(linspace(x1, x2, Npoints)', 1, Ncurves);
end

function [x, y] = get_curves(x1, y1, x2, y2)
% x1, x2: scalar x coordinates of line start, end
% y1, y2: vectors of y coordinates of line start/ends
    Npoints = 15;
    t = linspace(0, pi, Npoints);
    c = (1-cos(t))./2; % Normalized curve

    Ncurves = numel(y1);
	% Starting R2016b, the following line could be written simply as:
    %   y = y1 + (y2 - y1) .* c';
    y = repmat(y1, Npoints, 1) + repmat(y2 - y1, Npoints,1) .* repmat(c', 1, Ncurves);
    x = repmat(linspace(x1, x2, Npoints)', 1, Ncurves);
end  % get_curve
