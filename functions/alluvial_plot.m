%% --- --- prepare data -
clc;clear data;clear ttd;clear ttm;
%load(fullfile(spath,'OD.mat')); load(fullfile(spath,'Coord.mat'));
load(fullfile(spath,'OD.mat'));
load(fullfile(spath,'Coord.mat'));
myType = flip(unique(OD.type),2);

list_id = [];myexp_labels = {};nn = 0;n_Data=[];
for t = 1: length(myType)
    typeName = myType{t};   
    tableid = find(contains(OD.type,typeName));
    list_id = [list_id tableid];
    data(t).type = typeName;
    data(t).table = load(fullfile(spath,[typeName,'.mat']));
    thr = []; thr = data(t).table.Threshold;
    n_Data(1,t)=size(data(t).table.CoordMerge,2);
    for j = 1:length(tableid)
        nn = nn+1;
        myexp_labels{nn,1} = [typeName,OD.group{tableid(j)}];
        myexp_labels{nn,2} = OD.exp{tableid(j)};
    end
end
set(0,'defaultAxesFontSize',10);

% read Lat & Med data
for t = 1:2
    labels = BBB(:,[1,3,6]); ll = length(BBB);
    density = zeros(ll,n_Data(1,t));percent = zeros(ll,n_Data(1,t));

    dataCur = data(t).table.TT;data_tmp = [];
    threCur = data(t).table.Threshold;
    percent = zeros(ll,n_Data(1,t)); density = zeros(ll,n_Data(1,t));
    nrcount = zeros(ll,n_Data(1,t));
    for i = 1:ll
        data_tmp = [cell2mat(dataCur(i,7))];
        if ~isempty(data_tmp)
            nrcount(i,:) = data_tmp;
            percent(i,:) = [cell2mat(dataCur(i,9))];%percentage based on total inputs - SCs
            density(i,:) = [cell2mat(dataCur(i,14))];%density normalized to SC density
        end
    end

    density_err = std(density,[],2);%standard deviation
    %percent_err = std(percent,[],2)./sqrt(n_Data);%standard error of the mean (SEM)
    percent_err = std(percent,[],2);%standard deviation

    % data
    if t == 2
        countLat = nrcount;densityLat = density;
        percentLat = percent;threshlLat = threCur;
        errPerLat = percent_err;
        errDenLat = density_err;
    else
        countMed = nrcount;densityMed = density;
        percentMed = percent;threshlMed = threCur;
        errPerMed = percent_err;
        errDenMed = density_err;
    end

end

% remove SC
idSC = find(startsWith(labels(:,1),'SC'));%removes all SC layers from analysis
%filter out ctx
idnonCTX = find(~(cell2mat(labels(:,3))==1));

% filter
nanIDX1 = [];nanIDX2 = [];nanIDX3 = [];
thresT1 = 0;   nanIDX1 = find((threshlLat == thresT1) & (threshlMed == thresT1));%threshold  == 1
dmLat = mean(densityLat,2);dmMed = mean(densityMed,2);
thresT2 = 0;  nanIDX2 = find((dmLat < thresT2) & (dmMed < thresT2)); %normalized_density > thresT2
pmLat = mean(percentLat,2);pmMed = mean(percentMed,2);
thresT3 = 0.3; nanIDX3 = find((pmLat < thresT3) & (pmMed < thresT3)); %percentage > thresT3
nanIDX = unique([nanIDX1; nanIDX2; nanIDX3;idSC]);
n_Area = ll;n_include = n_Area - length(nanIDX);
clc;disp(['total areas = ',mat2str(n_Area),' - ',mat2str(length(nanIDX))]);
disp(['Inlcude ' ,mat2str(n_include),' areas for both circuits']);

countLat(nanIDX,:) = [];countMed(nanIDX,:) = []; %sum(countLat, 'all') counts the total number of lateral inputs across experiments
labels(nanIDX,:) = [];
percentLat(nanIDX,:) = [];percentMed(nanIDX,:) = [];
pmMed(nanIDX) = [];pmLat(nanIDX) = [];
densityLat(nanIDX,:) = [];densityMed(nanIDX,:) = [];
dmMed(nanIDX) = [];dmLat(nanIDX) = [];

errDenLat(nanIDX) = [];errDenMed(nanIDX) = [];
errPerLat(nanIDX) = [];errPerMed(nanIDX) = [];
disp(['Total number of included medial input cells: ' ,num2str(sum(countMed,'all'))]);
disp(['Total number of included lateral input cells: ' ,num2str(sum(countLat,'all'))]);
disp('data loaded');


%-----Med + Lat
%abbr = {};
plot_seq_perden=[];
labels_perden={};

for ff = 1
    if ff == 1 % percentage
        labelname = '\%';
        property = 'percentage';
        xx = mean(percentMed,2);
        yy = mean(percentLat,2);
        errY = errPerLat;errX = errPerMed;
        individuals_m = percentMed;
        individuals_l = percentLat;
    elseif ff == 2 % density
        labelname = 'norm_den';
        property = 'density';
        xx = mean(densityMed,2);
        yy = mean(densityLat,2);
        errY = errDenLat;errX = errDenMed;
        individuals_m = densityMed;
        individuals_l = densityLat;
    else
    end
end


%sort data depending on med+lat value
zz=xx+yy;
[ss, plot_seq] = sortrows(zz, 'descend');newlabels = labels(plot_seq);
xx=xx(plot_seq);yy=yy(plot_seq);
errY=errY(plot_seq);errX=errX(plot_seq);
individuals_m=individuals_m(plot_seq,:);individuals_l=individuals_l(plot_seq,:);
clearvars ss zz;
plot_seq_perden(:,ff)=plot_seq;
labels_perden(:,ff)=newlabels;
abbr = labels(plot_seq,1:2);
writetable(cell2table(abbr),fullfile(spath,'med+lat_abbr.csv'));

nGroups = 2;groupwidth = min(0.7, nGroups/(nGroups + 1));
lr = length(yy);
%h = bar(1:lr,[xx yy],'FaceColor','flat','BarWidth',0.98);hold on;
%group data so that you can access medial and lateral with for loop
%x positions for medial and lateral
x={};
for k = 1:nGroups
    x{1,k} = (1:lr) - groupwidth/2 + (2*k-1) * groupwidth/(2*nGroups);
end
ydata=[xx yy];%y mean positions for medial and lateral
errors=[errX errY];
individuals={individuals_m individuals_l};

%% Alluvial plot

right_labels = newlabels'; %43 areas
left_labels = {'Med', 'Lat'}; %med + lat

data=[];
data(1,:)=xx; data(2,:)=yy;

alluvialflow(data, left_labels, right_labels, 'Percentage input to medial and lateral ntsr', col);

set(gca, 'XDir', 'reverse');
set(gcf, 'renderer', 'Painters');