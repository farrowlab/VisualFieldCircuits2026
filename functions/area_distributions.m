%% --- prepare data for heatmap of distribution overlap between med and lat
av_mod=readNPY([atlasPath filesep 'annotation_volume_10um_merged_cortical_areas.npy']);
clc;clear data;clear ttd;clear ttm;
%load(fullfile(spath,'OD.mat')); load(fullfile(spath,'Coord.mat'));
celltypes = {'grp'};
list_id = [];myexp_labels = {};nn = 0;n_Data=[];OD={}; Coord={};tt=0;
n_total={};
for c=1:length(celltypes)
    
    OD{c}=load(fullfile(spath,['OD_',celltypes{c},'.mat']));
    Coord{c}=load(fullfile(spath,['Coord_',celltypes{c},'.mat']));
    
    myType = flip(unique(OD{c}.OD_cell.type),2);
        
    for t = 1: length(myType)
        typeName = myType{t};   
        tableid = find(contains(OD{c}.OD_cell.type,typeName));
        list_id = [list_id tableid];
        tt=tt+1;
        data(tt).type = typeName;
        data(tt).table = load(fullfile(spath,[celltypes{c}, '_', typeName,'.mat']));
        thr = []; thr = data(tt).table.Threshold;
        n_Data(c,t)=size(data(tt).table.CoordMerge,2); 
        for j = 1:length(tableid)
            nn = nn+1;
            myexp_labels{nn,1} = [typeName,OD{c}.OD_cell.group{tableid(j)}];
            myexp_labels{nn,2} = OD{c}.OD_cell.exp{tableid(j)};
            n_total{c,t}(j)=OD{c}.OD_cell.nTotalminusSC{tableid(j)};
        end
    end
end
set(0,'defaultAxesFontSize',10);

n_Data=reshape(n_Data.', 1, []);
totalMed={}; totalLat={};
countLat={};densityLat={};percentLat={};threshlLat={};errPerLat={};errDenLat={};
countMed={};densityMed={};percentMed={};threshlMed={};errPerMed={};errDenMed={};
labels = BBB(:,[1,3,6,8]); ll = length(BBB);
% read Lat & Med data
for t = 1:2 %2 if 1 celltype, 2 if 4
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
            density(i,:) = [cell2mat(dataCur(i,16))];%density normalized to sum of all densitiess exclusing SCs
        end
    end

    % density_err = std(density,[],2);%standard deviation
    density_err = std(density,[],2)./sqrt(n_Data);%standard error of the mean (SEM)
    percent_err = std(percent,[],2)./sqrt(n_Data);%standard error of the mean (SEM)
    % percent_err = std(percent,[],2);%standard deviation

    % data
    if mod(t, 2) == 0       
        tt=t/2;
        countLat{tt} = nrcount;densityLat{tt} = density;
        percentLat{tt} = percent;threshlLat{tt} = threCur;
        errPerLat{tt} = percent_err;
        errDenLat{tt} = density_err;
    else
        tt=1+0.5*(t-1);
        countMed{tt} = nrcount;densityMed{tt} = density;
        percentMed{tt} = percent;threshlMed{tt} = threCur;
        errPerMed{tt} = percent_err;
        errDenMed{tt} = density_err;
    end

end

% remove SC
idSC = find(startsWith(labels(:,1),'SCs'));%removes all SC layers from analysis
%remove MBsen,MBmot,MBsta
idMBdiv=find(ismember(labels(:,1),["MBsen", "MBmot", "MBsta"]));
%filter out ctx
idnonCTX = find(~(cell2mat(labels(:,3))==1));

%filter
dmLat={}; dmMed= {}; pmLat={}; pmMed={}; nanIDX={};
for c=1:length(celltypes)
    nanIDX1 = [];nanIDX2 = [];nanIDX3 = [];
    thresT1 = 0;   nanIDX1 = find((threshlLat{c} == thresT1) & (threshlMed{c} == thresT1));%threshold  == 1
    dmLat{c} = mean(densityLat{c},2);dmMed{c} = mean(densityMed{c},2);
    thresT2 = 0.3;  nanIDX2 = find((dmLat{c} < thresT2) & (dmMed{c} < thresT2)); %normalized_density > thresT2
    pmLat{c} = median(percentLat{c},2);pmMed{c} = median(percentMed{c},2);
    thresT3 = 0.35; nanIDX3 = find((pmLat{c} < thresT3) & (pmMed{c} < thresT3)); %percentage > thresT3
    nanIDX{c} = unique([nanIDX1; nanIDX2; nanIDX3;idSC; idMBdiv]);
end

comb_nanIDX=vertcat(nanIDX{:}); %combine all cells in nanIDX
unique_nanIDX=unique(comb_nanIDX);
labels(unique_nanIDX,:) = [];

for c=1:length(celltypes)
    countLat{c}(unique_nanIDX,:) = [];countMed{c}(unique_nanIDX,:) = []; %sum(countLat, 'all') counts the total number of lateral inputs across experiments
    percentLat{c}(unique_nanIDX,:) = [];percentMed{c}(unique_nanIDX,:) = [];
    pmMed{c}(unique_nanIDX) = [];pmLat{c}(unique_nanIDX) = [];
    densityLat{c}(unique_nanIDX,:) = [];densityMed{c}(unique_nanIDX,:) = [];
    dmMed{c}(unique_nanIDX) = [];dmLat{c}(unique_nanIDX) = [];
    errDenMed{c}(unique_nanIDX) = []; errDenLat{c}(unique_nanIDX) = [];
    errPerMed{c}(unique_nanIDX) = []; errPerLat{c}(unique_nanIDX) = [];
    nTot{c,1}=sum(countMed{c},1); nTot{c,2}=sum(countLat{c},1);
    disp(['Total number of included ', celltypes{c}, ' medial input cells: ' ,num2str(sum(countMed{c},'all'))]);
    disp(['Total number of included ', celltypes{c},' lateral input cells: ' ,num2str(sum(countLat{c},'all'))]);
    disp(['Cell counts for each ', celltypes{c}, ' brain:']);
    disp(['Med ', num2str(sum(countMed{c},1))]);
    disp(['Lat ', num2str(sum(countLat{c},1))]);
    disp(['Average cell counts per ', celltypes{c}, ' brain: Med ', num2str(round(mean(sum(countMed{c},1)))), ', Lat ',...
    num2str(round(mean(sum(countLat{c},1))))]);
end
% thresT1 = 0;   nanIDX1 = find(threshlMed == thresT1);%threshold  == 1
% dmMed = mean(densityMed,2);
% thresT2 = 0.2;  nanIDX2 = find(dmMed < thresT2); %normalized_density > thresT2
% pmMed = mean(percentMed,2);
% thresT3 = 0.2; nanIDX3 = find(pmMed < thresT3); %percentage > thresT3
% nanIDX = unique([nanIDX1; nanIDX2; nanIDX3;idSC; idMBdiv]);
n_Area = ll;n_include = n_Area - length(unique_nanIDX);
disp(['total areas = ',mat2str(n_Area),' - ',mat2str(length(unique_nanIDX))]);
disp(['Include ' ,mat2str(n_include),' areas for both circuits']);

abbr = {};

labelname = 'log2$\mathrm{^{(\%Med/\%Lat)}}$';
property = 'percentage';
%pmMed = mean(percentMed,2); just to know the definition of pmMed and pmLat
%pmLat = mean(percentLat,2);
individuals_m={}; individuals_l={}; lgData={};

%sort data based on bias
for c=1:length(celltypes)
    individuals_m{c} = percentMed{c};
    individuals_l{c} = percentLat{c};
    %         lg_ind=log2(percentMed./percentLat);
    lgData{c} = log2(pmMed{c}./pmLat{c});
    % ---- in the case of Inf/-Inf ---- %
    mx = ceil(max(lgData{c}(lgData{c} ~= Inf))); mn = floor(min(lgData{c}(lgData{c} ~= -Inf)));
    inf_val=max(mx, abs(mn));
    mx_id = lgData{c} == Inf; lgData{c}(mx_id,:) = inf_val;
    mn_id = lgData{c} == -Inf;lgData{c}(mn_id,:) = -inf_val;
end


figName = ['Summary_Bias_',property];

% ---- sort the data ---- %
[xx, plot_seq] = sortrows(lgData{1}, 'descend');xx_labels = labels(plot_seq,1);xx_labels1 = labels(plot_seq,3);
% lg_ind=lg_ind(plot_seq,:);
abbr{1,1} = property;
abbr{1,2} = labels(plot_seq,1:2);

for c=1:length(celltypes)
    lgData{c}=lgData{c}(plot_seq);
end



%% calculate distribution overlaps of areas along three main brain axes

%distrubition ML=3; AP=1; DV=2;
celltypes={'grp'};

areas=string(xx_labels);
ML_ksd=[];
AP_ksd=[];
DV_ksd=[];

ML_kstest=[]; AP_kstest=[]; DV_kstest=[];
ML_pval=[]; AP_pval=[]; DV_pval=[];
for a=1:length(areas)
    Coordinates_tot=[];coord_tot={};
    for t = 1:length(myType)
        coll = col(t+2,:);
        typeName = myType{t};
        dt = load(fullfile(spath,[celltypes{1}, '_',typeName,'.mat']));
        dt.TT_title;
        labels = dt.TT(:,1:3);
        idSC = find(ismember(labels(:,1),areas(a)));
        threCur= dt.Threshold;
        myCoord = dt.CoordMerge;
%         for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
%             if threCur(i,1)==0
%                 myCoord(i,:)={[]};
%             end
%         end
        myCoord=myCoord(idSC,:);

        %plot transformed coordinates
        %figure(f);
        [sa, sl] = size(myCoord);coordXYZ = [];
        %     theta=deg2rad(-28.022);
        %     phi=deg2rad(-8.29);
        Coord_type=[];
        for i = 1:sl
            coordXYZ = cell2mat(myCoord(:,i));
            if ~isempty(coordXYZ)
                coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:); %keep only ipsilateral
                coordXYZ=coordXYZ+randn(size(coordXYZ))/300;%add some small random noise to visualize scatter points that would be on the same pixel
                coordXYZ(:,1)=coordXYZ(:,1)+(-0.05 + (0.05+0.05).*rand(size(coordXYZ,1),1)); %jitter points along AP axis to cover slice spacing effect
                %         X=coordXYZ(:,1);
                %         Y=coordXYZ(:,3);
                %         Z=coordXYZ(:,2);
                %         Xnew=X*cos(phi)-Y*sin(theta)*sin(phi)-Z*cos(theta)*sin(phi);
                %         Ynew=Y*cos(theta)-Z*sin(theta);
                %         Znew=X*sin(phi)+Y*sin(theta)*cos(phi)+Z*cos(theta)*cos(phi);
                %         coordXYZ=[Xnew Znew Ynew];
                Coordinates_tot=[Coordinates_tot; coordXYZ];
                Coord_type=[Coord_type; coordXYZ];
            else
                coordXYZ=[NaN NaN NaN];
                Coord_type=[Coord_type; coordXYZ];
            end
        end
        coord_tot{1,t}=Coord_type;
    end
    %---ML---%
    [y_med,~] = ksdensity(coord_tot{1,1}(:,3),'Support', [min(coord_tot{1,1}(:,3))-0.1 ...
    max(coord_tot{1,1}(:,3))+0.1], 'NumPoints',300); %ML
    x_med=linspace(min(coord_tot{1,1}(:,3)), max(coord_tot{1,1}(:,3)),300);
    %plot(x_med_horizontal,y_med_horizontal,'LineWidth',2, 'Color',col(3,:)); hold on;
    [y_lat,~] = ksdensity(coord_tot{1,2}(:,3),'Support', [min(coord_tot{1,2}(:,3))-0.1 ...
    max(coord_tot{1,2}(:,3))+0.1], 'NumPoints',300);
    x_lat=linspace(min(coord_tot{1,2}(:,3)), max(coord_tot{1,2}(:,3)),300);
    %plot(x_lat_horizontal,y_lat_horizontal,'LineWidth',2, 'Color',col(4,:));
    newx=linspace(min([min(coord_tot{1,1}(:,3)) min(coord_tot{1,2}(:,3))]), ...
        max([max(coord_tot{1,1}(:,3)) max(coord_tot{1,2}(:,3))]),300); %shared x-axis between the two distributions
    new_y_med=interp1(x_med,y_med,newx); %interpolate y data to match shared x-axis
    new_y_lat=interp1(x_lat,y_lat,newx);
    pmin=min(new_y_med,new_y_lat);
    pmax=max(new_y_med,new_y_lat);
    overlap=sum(pmin)./sum(pmax);
    ML_ksd=[ML_ksd; overlap];

    [~, pval,stat] = kstest2(coord_tot{1,1}(:,3), coord_tot{1,2}(:,3));
    ML_kstest=[ML_kstest; stat];
    ML_pval=[ML_pval; pval];

    %---AP---%
    [y_med,~] = ksdensity(coord_tot{1,1}(:,1),'Support', [min(coord_tot{1,1}(:,1))-0.1 ...
    max(coord_tot{1,1}(:,1))+0.1], 'NumPoints',300); %AP
    x_med=linspace(min(coord_tot{1,1}(:,1)), max(coord_tot{1,1}(:,1)),300);
    %plot(x_med_horizontal,y_med_horizontal,'LineWidth',2, 'Color',col(3,:)); hold on;
    [y_lat,~] = ksdensity(coord_tot{1,2}(:,1),'Support', [min(coord_tot{1,2}(:,1))-0.1 ...
    max(coord_tot{1,2}(:,1))+0.1], 'NumPoints',300);
    x_lat=linspace(min(coord_tot{1,2}(:,1)), max(coord_tot{1,2}(:,1)),300);
    %plot(x_lat_horizontal,y_lat_horizontal,'LineWidth',2, 'Color',col(4,:));
    newx=linspace(min([min(coord_tot{1,1}(:,1)) min(coord_tot{1,2}(:,1))]), ...
        max([max(coord_tot{1,1}(:,1)) max(coord_tot{1,2}(:,1))]),300); %shared x-axis between the two distributions
    new_y_med=interp1(x_med,y_med,newx); %interpolate y data to match shared x-axis
    new_y_lat=interp1(x_lat,y_lat,newx);
    pmin=min(new_y_med,new_y_lat);
    pmax=max(new_y_med,new_y_lat);
    overlap=sum(pmin)./sum(pmax);
    AP_ksd=[AP_ksd; overlap];

    [~, pval,stat] = kstest2(coord_tot{1,1}(:,1), coord_tot{1,2}(:,1));
    AP_kstest=[AP_kstest; stat];
    AP_pval=[AP_pval; pval];

    %---DV---%
    [y_med,~] = ksdensity(coord_tot{1,1}(:,2),'Support', [min(coord_tot{1,1}(:,2))-0.1 ...
    max(coord_tot{1,1}(:,2))+0.1], 'NumPoints',300); %DV
    x_med=linspace(min(coord_tot{1,1}(:,2)), max(coord_tot{1,1}(:,2)),300);
    %plot(x_med_horizontal,y_med_horizontal,'LineWidth',2, 'Color',col(3,:)); hold on;
    [y_lat,~] = ksdensity(coord_tot{1,2}(:,2),'Support', [min(coord_tot{1,2}(:,2))-0.1 ...
    max(coord_tot{1,2}(:,2))+0.1], 'NumPoints',300);
    x_lat=linspace(min(coord_tot{1,2}(:,2)), max(coord_tot{1,2}(:,2)),300);
    %plot(x_lat_horizontal,y_lat_horizontal,'LineWidth',2, 'Color',col(4,:));
    newx=linspace(min([min(coord_tot{1,1}(:,2)) min(coord_tot{1,2}(:,2))]), ...
        max([max(coord_tot{1,1}(:,2)) max(coord_tot{1,2}(:,2))]),300); %shared x-axis between the two distributions
    new_y_med=interp1(x_med,y_med,newx); %interpolate y data to match shared x-axis
    new_y_lat=interp1(x_lat,y_lat,newx);
    pmin=min(new_y_med,new_y_lat);
    pmax=max(new_y_med,new_y_lat);
    overlap=sum(pmin)./sum(pmax);
    DV_ksd=[DV_ksd; overlap];

    [~, pval,stat] = kstest2(coord_tot{1,1}(:,2), coord_tot{1,2}(:,2));
    DV_kstest=[DV_kstest; stat];
    DV_pval=[DV_pval; pval];

end

%correct for multiple comparisons
[ML_padj,~] = multicmp(ML_pval,'down');% Holm-Bonferroni (type: 'down'), Hochberg-Bonferroni (type: 'up'), or Benjamini-Hochberg (type: 'fdr') correction
[AP_padj,~] = multicmp(AP_pval,'down');
[DV_padj,~] = multicmp(DV_pval,'down');
retinotopy_overlap=table(areas, ML_ksd, AP_ksd, DV_ksd);
retinotopy_difference=table(areas, ML_kstest, AP_kstest, DV_kstest);
%uses the maximum absolute difference between the cumulative distribution functions (cdfs) 
%of the distributions calculated in the kstest2 stat
retinotopy_significance=table(areas, ML_padj, AP_padj, DV_padj);
numbers=[ML_ksd';AP_ksd';DV_ksd'];

%col1=flipud(col1);

%% plot heatmap
col1=magma(300);
col1=flipud(col1);
f=figure(1);set(f,'Position', [73 422.6000 992 148.8000]);
% til=tiledlayout(length(par_seq),8, 'TileSpacing','none', 'Padding','tight');
% set(til, 'OuterPosition', [0 0 0.8 1]);
hm=heatmap(numbers, 'Colormap',col1,'ColorLimits',[0,1],...%[round(min(min(numbers)),2) round(max(max(numbers)),2)],...
    'GridVisible','off', 'CellLabelColor','none');
hm.Title= 'Separation of medial and lateral input distributions along three main brain axes';
hm.XDisplayLabels = xx_labels;
hm.YDisplayLabels={'ML', 'AP', 'DV'};
s = struct(hm);
s.XAxis.TickLabelRotation = 90; 
set(hm, 'OuterPosition', [0 0 1 0.9]);
tl_pos=get(s.Axes.Title, 'Position');
set(s.Axes.Title, 'FontSize', 11, 'Position', [tl_pos(1) tl_pos(2)-0.2 tl_pos(3)]);
cl=s.Colorbar;
set(cl, 'YDir', 'reverse')
cl.Ticks=[0 0.5 1];
cl.YTickLabel={'1', '0.5', '0'};
cl.Label.String=["separation"; "fraction"];
cl.Label.FontSize=9;


%hm.XDisplayLabels.TickLabelRotation=90;
% hm.YDisplayLabels = nan(size(hm.YDisplayData));
% hm.XDisplayLabels = nan(size(hm.XDisplayData));
% hm.YDisplayLabels = nan(size(hm.YDisplayData));

%% plot 2D scatter of single areas with ksdensity distributions

Coordinates_tot=[];
coord_tot={};
area=["APN"];

AX1=3; %distrubition along x-axis; ML=3; AP=1; DV=2;
AX2=2; %distribution along y-axis

for t = 1:length(myType)
    coll = col(t+2,:);
    typeName = myType{t};
    dt = load(fullfile(spath,[celltypes{1}, '_',typeName,'.mat']));
    dt.TT_title;
    labels = dt.TT(:,1:3);
    idSC = find(ismember(labels(:,1),area));
    threCur= dt.Threshold;
    myCoord = dt.CoordMerge;
%     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
%         if threCur(i,1)==0
%             myCoord(i,:)={[]};
%         end
%     end
    myCoord=myCoord(idSC,:);

    %plot transformed coordinates
    %figure(f);
    [sa, sl] = size(myCoord);coordXYZ = [];
        % theta=deg2rad(-28.022); %along M-L axis
        % phi=deg2rad(-8.29); %along D-V axis
        theta=deg2rad(0); %along M-L axis
        phi=deg2rad(0); %along D-V axis
    Coord_type=[];
    for i = 1:sl
        coordXYZ = cell2mat(myCoord(:,i));
        if ~isempty(coordXYZ)
            coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:); %keep only ipsilateral
            coordXYZ=coordXYZ+randn(size(coordXYZ))/300;%add some small random noise to visualize scatter points that would be on the same pixel
            coordXYZ(:,1)=coordXYZ(:,1)+(-0.05 + (0.05+0.05).*rand(size(coordXYZ,1),1)); %jitter points along AP axis to cover slice spacing effect
            X=coordXYZ(:,1);
            Y=coordXYZ(:,3);
            Z=coordXYZ(:,2);
            Xnew=X*cos(phi)-Y*sin(theta)*sin(phi)-Z*cos(theta)*sin(phi);
            Ynew=Y*cos(theta)-Z*sin(theta);
            Znew=X*sin(phi)+Y*sin(theta)*cos(phi)+Z*cos(theta)*cos(phi);
            coordXYZ=[Xnew Znew Ynew];
            Coordinates_tot=[Coordinates_tot; coordXYZ];
            Coord_type=[Coord_type; coordXYZ];
        else
        coordXYZ=[NaN NaN NaN];
        Coord_type=[Coord_type; coordXYZ]; 
        end
    end
    coord_tot{1,t}=Coord_type;
end

t=3; %1 for ntsr, 3 for grp
f=figure(1);
fact=2;
til=tiledlayout(3*fact,3*fact, 'TileSpacing','none', 'Padding','tight');
set(til, 'OuterPosition', [0 0 0.98 0.98]);
til.Title.String=area;
til.Title.FontSize=14;
til.Title.FontWeight='bold';
%plot horizontal density curve
nexttile([1 2*fact+1]);
% h=histfit(-coord_tot{1,1}(:,3),100,'kernel');delete(h(1));h(2).Visible='off'; %medial
% med_x=h(2).XData; med_y=h(2).YData;
% h=histfit(-coord_tot{1,2}(:,3),100,'kernel');delete(h(1));h(2).Visible='off'; %lateral
% lat_x=h(2).XData; lat_y=h(2).YData;
% plot(med_x, med_y/max(max(med_y), max(lat_y)),'LineWidth',2, 'Color',col(3,:)); hold on;
% plot(lat_x, lat_y/max(max(med_y), max(lat_y)),'LineWidth',2, 'Color',col(4,:));
[y_med,~] = ksdensity(coord_tot{1,1}(:,AX1),'Support', [min(coord_tot{1,1}(:,AX1))-0.1 ...
    max(coord_tot{1,1}(:,AX1))+0.1], 'NumPoints',300); %ML
x_med=linspace(min(coord_tot{1,1}(:,AX1)), max(coord_tot{1,1}(:,AX1)),300);
plot(x_med,y_med,'LineWidth',2, 'Color',col(t,:)); hold on;
[y_lat,~] = ksdensity(coord_tot{1,2}(:,AX1), 'Support', [min(coord_tot{1,2}(:,AX1))-0.1 ...
    max(coord_tot{1,2}(:,AX1))+0.1], 'NumPoints',300);%'BandWidth', 0.5
x_lat=linspace(min(coord_tot{1,2}(:,AX1)), max(coord_tot{1,2}(:,AX1)),300);
plot(x_lat,y_lat,'LineWidth',2, 'Color',col(t+1,:));

newx=linspace(min([min(coord_tot{1,1}(:,AX1)) min(coord_tot{1,2}(:,AX1))]), ...
    max([max(coord_tot{1,1}(:,AX1)) max(coord_tot{1,2}(:,AX1))]),300); %shared x-axis between the two distributions
new_y_med=interp1(x_med,y_med,newx); %interpolate y data to match shared x-axis
new_y_lat=interp1(x_lat,y_lat,newx);
pmin=min(new_y_med,new_y_lat);
pmax=max(new_y_med,new_y_lat);

overlap_horizontal=sum(pmin)./sum(pmax); 

[~, pval_horizontal,stat_horizontal] = kstest2(coord_tot{1,1}(:,AX1), coord_tot{1,2}(:,AX1));

xlim(gca,[(min(Coordinates_tot(:,AX1))-0.3) (max(Coordinates_tot(:,AX1))+0.3)]); 
set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');

ax1=nexttile(4*fact-1,[2*fact+1 2*fact+1]);
%plot scatter points
s_med=scatter(coord_tot{1,1}(:,AX1), coord_tot{1,1}(:,AX2), 10, 'filled', CData=col(t,:), MarkerFaceAlpha=0.6); hold on;
s_lat=scatter(coord_tot{1,2}(:,AX1), coord_tot{1,2}(:,AX2), 10, 'filled', CData=col(t+1,:), MarkerFaceAlpha=0.6);hold on;
uistack(s_med, 'top');
set(gca,'box', 'on', 'Color', 'none');
xlim(gca,[(min(Coordinates_tot(:,AX1))-0.3) (max(Coordinates_tot(:,AX1))+0.3)]); 
ylim(gca,[(min(Coordinates_tot(:,AX2))-0.3) (max(Coordinates_tot(:,AX2))+0.3)]);
if AX2==2
    set(gca, 'YDir', 'reverse')
end
ax=gca; ax.FontSize=12;
%exes=get(ax,XTick);

if AX1==1
    xaxis='A-P';
    ax.XTickLabel = split(num2str(xticks*(-1)))';
elseif AX1==2
    xaxis='D-V';
    ax.XTickLabel = split(num2str(xticks*(-1)))';
else
    xaxis='M-L';
end

if AX2==1
    yaxis='A-P';
    ax.YTickLabel = split(num2str(yticks*(-1)))';
elseif AX2==2
    yaxis='D-V';
    ax.YTickLabel = split(num2str(yticks*(-1)))';
else
    yaxis='M-L';
end

xlabel(xaxis, 'FontSize', 14); ylabel(yaxis,'FontSize', 14);

set(gca, 'FOntSize', 12);

%plot vertical density curve
nexttile(6*fact, [2*fact+1 1])
[y_med_vertical,~] = ksdensity(coord_tot{1,1}(:,AX2),'Support', [min(coord_tot{1,1}(:,AX2))-0.1 ...
    max(coord_tot{1,1}(:,AX2))+0.1], 'NumPoints',300);%, 'BandWidth', 0.2);
x_med_vertical=linspace(min(coord_tot{1,1}(:,AX2)), max(coord_tot{1,1}(:,AX2)),300);
plot(y_med_vertical,x_med_vertical,'LineWidth',2, 'Color',col(t,:)); hold on;
[y_lat_vertical,~] = ksdensity(coord_tot{1,2}(:,AX2),'Support', [min(coord_tot{1,2}(:,AX2))-0.1 ...
    max(coord_tot{1,2}(:,AX2))+0.1], 'NumPoints',300);
x_lat_vertical=linspace(min(coord_tot{1,2}(:,AX2)), max(coord_tot{1,2}(:,AX2)),300);
plot(y_lat_vertical,x_lat_vertical,'LineWidth',2, 'Color',col(t+1,:));

newx_vertical=linspace(min([min(coord_tot{1,1}(:,AX2)) min(coord_tot{1,2}(:,AX2))]), ...
    max([max(coord_tot{1,1}(:,AX2)) max(coord_tot{1,2}(:,AX2))]),300);
new_y_med_vertical=interp1(x_med_vertical,y_med_vertical,newx_vertical);
new_y_lat_vertical=interp1(x_lat_vertical,y_lat_vertical,newx_vertical);
pmin_vertical=min(new_y_med_vertical,new_y_lat_vertical);
pmax_vertical=max(new_y_med_vertical,new_y_lat_vertical);

overlap_vertical=sum(pmin_vertical)./sum(pmax_vertical);

[~, pval_vertical,stat_vertical] = kstest2(coord_tot{1,1}(:,AX2), coord_tot{1,2}(:,AX2)); %use kstest2 as a measure of overlap

ylim(gca,[(min(Coordinates_tot(:,AX2))-0.3) (max(Coordinates_tot(:,AX2))+0.3)]);
set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
if AX2==2
    set(gca, 'YDir', 'reverse')
end

%set(gca,'box','on', 'Color', 'none');%put 'color' zerocolor when plotting imagesc
%make figure square so that scatter axes are equal
x_range=(max(Coordinates_tot(:,AX1))+0.3) -(min(Coordinates_tot(:,AX1))-0.3);
y_range=(max(Coordinates_tot(:,AX2))+0.3) -(min(Coordinates_tot(:,AX2))-0.3);
xy_ratio=y_range./x_range;
fig_width=400;
set(f, "Position", [572.2000 61 fig_width fig_width*xy_ratio], 'renderer','Painters');


