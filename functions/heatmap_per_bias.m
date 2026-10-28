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
            density(i,:) = [cell2mat(dataCur(i,16))];%density normalized to sum of all densitiess exclusing SCs
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
idSC = find(startsWith(labels(:,1),'SCs'));%removes all SC layers from analysis
%remove MBsen,MBmot,MBsta
idMBdiv=find(ismember(labels(:,1),["MBsen", "MBmot", "MBsta"]));
%filter out ctx
idnonCTX = find(~(cell2mat(labels(:,3))==1));

% filter
nanIDX1 = [];nanIDX2 = [];nanIDX3 = [];
thresT1 = 0;   nanIDX1 = find((threshlLat == thresT1) & (threshlMed == thresT1));%threshold  == 1
dmLat = mean(densityLat,2);dmMed = mean(densityMed,2);
thresT2 = 0.3;  nanIDX2 = find((dmLat < thresT2) & (dmMed < thresT2)); %normalized_density > thresT2
pmLat = mean(percentLat,2);pmMed = mean(percentMed,2);
thresT3 = 0.2; nanIDX3 = find((pmLat < thresT3) & (pmMed < thresT3)); %percentage > thresT3
nanIDX = unique([nanIDX1; nanIDX2; nanIDX3;idSC; idMBdiv]);
n_Area = ll;n_include = n_Area - length(nanIDX);
clc;disp(['total areas = ',mat2str(n_Area),' - ',mat2str(length(nanIDX))]);
disp(['Include ' ,mat2str(n_include),' areas for both circuits']);

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
disp(['Average cell counts per brain: Med ', num2str(round(mean(sum(countMed,1)))), ', Lat ',...
    num2str(round(mean(sum(countLat,1))))]);

%% Med + Lat preprocess data and group by major brain region

%percentage data
labelname = '% of total inputs';
property = 'percentage';
xx = mean(percentMed,2);
yy = mean(percentLat,2);
errY = errPerLat;errX = errPerMed;
individuals_m = percentMed;
individuals_l = percentLat;

%sort data depending on med+lat value of percentage
zz=xx+yy;
[~, plot_seq] = sortrows(zz, 'descend');newlabels = labels(plot_seq,:); 
[newlabels, seq2]=sortrows(newlabels,3,'ascend');%sort again based on major brain region
plot_seq=plot_seq(seq2); %update original sorting sequence with second sorting
xx=xx(plot_seq);yy=yy(plot_seq); 
errY=errY(plot_seq);errX=errX(plot_seq); 
individuals_m=individuals_m(plot_seq,:);individuals_l=individuals_l(plot_seq,:);
clearvars zz;

%calculate x positions for medial and lateral
nGroups = 2;groupwidth = min(0.7, nGroups/(nGroups))+0.7;
gr = length(yy);
x={};
for k = 1:nGroups
    x{1,k} = (1:2:gr*2) - groupwidth/2 + (2*k-1) * groupwidth/(2*nGroups);
end

%put med and lat data in one variable
ydata=[xx yy];%y mean positions for medial and lateral
errors=[errX errY];
individuals={individuals_m individuals_l};



% ------bias index calculation-----%
bias_name = 'log2$\mathrm{^{(\%Med/\%Lat)}}$';
lgData = log2(pmMed./pmLat);
lgData=lgData(plot_seq);%sort bias data based on above sorting
CI_upper = 0.37; CI_lower = -CI_upper;
  
%in the case of Inf/-Inf
mx = ceil(max(lgData(lgData ~= Inf))); mn = floor(min(lgData(lgData ~= -Inf)));
inf_val=max(mx, abs(mn));
mx_id = lgData == Inf; lgData(mx_id,:) = inf_val;
mn_id = lgData == -Inf;lgData(mn_id,:) = -inf_val;

%categorize bias index by color
nbins=6;
lr = length(lgData);
colors=zeros(lr, 3);
%medial
xxup=lgData>=CI_upper;
med_val=lgData(xxup); med_grp=linspace(min(med_val), max(med_val), nbins);
med_grp=fliplr(med_grp);
med_col=[];
for i=1:length(med_grp)-1
    for j=1:length(med_val)
        if med_val(j)<=med_grp(i)&& med_val(j)>=med_grp(i+1)
            med_col(j,:)=col_bias_med(i,:);
        end
    end
end
colors(xxup,:)=med_col;  

%grays
xx1=(lgData<CI_upper)&(lgData>CI_lower);
colors(xx1,:)=repmat([0.7, 0.7, 0.7], length(xx1(xx1)),1);

%lateral
xxdown=lgData<=CI_lower;
lat_val=lgData(xxdown); lat_grp=linspace(min(lat_val), max(lat_val), nbins);
lat_grp=fliplr(lat_grp);
lat_col=[];
for i=1:length(lat_grp)-1
    for j=1:length(lat_val)
        if lat_val(j)<=lat_grp(i)&& lat_val(j)>=lat_grp(i+1)
            lat_col(j,:)=col_bias_lat(end-i,:);
        end
    end
end
colors(xxdown,:)=lat_col;  

%check color order
% for i=1:length(colors)
%     rectangle('Position' ,[1, i, 0.5 0.5], 'FaceColor', colors(i, :), 'EdgeColor', 'none');
%     hold on
% end


%% ------------------------------PLOT with colors-------------------------------%
par_seq=cell2mat(newlabels(:,3));
num_areas=length(unique(par_seq)); %number of parental areas
f=figure(1);set(f,'Position', [488.2000 52.2000 322.4000 709.6000]);
til=tiledlayout(length(par_seq),8, 'TileSpacing','none', 'Padding','tight');

%-------CTX------%
nr=length(find(par_seq==1));
%bias index
ax1=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==1,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax1, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax1.YLabel.String="CTX";
ax1.YLabel.FontSize=8;
ax1.XAxisLocation="top";
ax1.XLabel.String=["Bias"; "index"];
ax1.XLabel.FontSize=9;


%heatmaps
ind_med=(individuals_m(par_seq==1,:));
ind_lat=(individuals_l(par_seq==1,:));
max_val1=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax2=nexttile([nr 3]);
col1=customcolormap([0 0.5 1], [0.6588 0.3294 0;1.0000 0.6 0.2; 1 1 1]);

hm=heatmap(ind_med,'Colormap',col1,'ColorLimits',[0 ceil(max_val1)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));
hm.Title="medial";
set(struct(hm).Axes.Title,'FontWeight','normal', 'FontSize', 9);

%lateral heatmap
ax3=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col1,'ColorLimits',[0 ceil(max_val1)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));
hm.Title="lateral";
set(struct(hm).Axes.Title,'FontWeight','normal', 'FontSize', 9);

%labels
ax4=nexttile([nr 1]);
curr_lab=newlabels(par_seq==1);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',4.5);
set(ax4, 'color', 'none');CurYLim=ylim;
axis off;

linkaxes([ax1 ax4], 'y');


%-------sCTX------%
nr=length(find(par_seq==2));
%bias index
ax5=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==2,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax5, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax5.YLabel.String="sCTX";
ax5.YLabel.FontSize=8;

%heatmaps
ind_med=(individuals_m(par_seq==2,:));
ind_lat=(individuals_l(par_seq==2,:));
max_val2=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax6=nexttile([nr 3]);
col2=customcolormap([0 1], [0.4706 0.3137 0; 1 1 1]);

hm=heatmap(ind_med,'Colormap',col2,'ColorLimits',[0 ceil(max_val2)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%lateral heatmap
ax7=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col2,'ColorLimits',[0 ceil(max_val2)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%labels
ax8=nexttile([nr 1]);
curr_lab=newlabels(par_seq==2);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',4.5);
set(ax8, 'color', 'none');
axis off;

linkaxes([ax5 ax8], 'y');


%-------TH------%
nr=length(find(par_seq==3));
%bias index
ax9=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==3,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax9, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax9.YLabel.String="TH";
ax9.YLabel.FontSize=8;
ylim([0 nr]);

%heatmaps
ind_med=(individuals_m(par_seq==3,:));
ind_lat=(individuals_l(par_seq==3,:));
max_val3=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax10=nexttile([nr 3]);
col3=customcolormap([0 0.5 1], [0.5686 0 0.0745; 0.9216 0.2118 0.3059;1 1 1]);

hm=heatmap(ind_med,'Colormap',col3,'ColorLimits',[0 ceil(max_val3)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%lateral heatmap
ax11=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col3,'ColorLimits',[0 ceil(max_val3)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%labels
ax12=nexttile([nr 1]);
curr_lab=newlabels(par_seq==3);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',4.5);
set(ax12, 'color', 'none');
axis off;

linkaxes([ax9 ax12], 'y');


%-------HYP------%
nr=length(find(par_seq==4));
%bias index
ax13=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==4,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax13, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax13.YLabel.String="HYP";
ax13.YLabel.FontSize=8;
ylim([0 nr]);

%heatmaps
ind_med=(individuals_m(par_seq==4,:));
ind_lat=(individuals_l(par_seq==4,:));
max_val4=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax15=nexttile([nr 3]);
col4=customcolormap([0 0.5 1], [0.4314 0.0941 0.4314; 0.6902 0.2902 0.6902;1 1 1]);

hm=heatmap(ind_med,'Colormap',col4,'ColorLimits',[0 ceil(max_val4)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%lateral heatmap
ax16=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col4,'ColorLimits',[0 ceil(max_val4)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%labels
ax17=nexttile([nr 1]);
curr_lab=newlabels(par_seq==4);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',4.5);
set(ax17, 'color', 'none');
axis off;

linkaxes([ax13 ax17], 'y');

%-------MB------%
nr=length(find(par_seq==5));
%bias index
ax18=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==5,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax18, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax18.YLabel.String="MB";
ax18.YLabel.FontSize=8;
ylim([0 nr]);

%heatmaps
ind_med=(individuals_m(par_seq==5,:));
ind_lat=(individuals_l(par_seq==5,:));
max_val5=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax19=nexttile([nr 3]);
col5=customcolormap([0 0.5 1], [0 0.3333 0.5412; 0.3059 0.6314 0.8314;1 1 1]);

hm=heatmap(ind_med,'Colormap',col5,'ColorLimits',[0 ceil(max_val5)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%lateral heatmap
ax20=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col5,'ColorLimits',[0 ceil(max_val5)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%labels
ax21=nexttile([nr 1]);
curr_lab=newlabels(par_seq==5);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',4.5);
set(ax21, 'color', 'none');
axis off;

linkaxes([ax18 ax21], 'y');

%-------HB------%
nr=length(find(par_seq==6));
%bias index
ax22=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==6,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax22, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax22.YLabel.String="HB";
ax22.YLabel.FontSize=8;
ax22.XLabel.String="Animal";
ax22.XLabel.FontSize=6;

%heatmaps
ind_med=(individuals_m(par_seq==6,:));
ind_lat=(individuals_l(par_seq==6,:));
max_val6=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax23=nexttile([nr 3]);
col6=customcolormap([0 0.5 1], [0.0706 0.5020 0.0039; 0.3059 0.7608 0.2353;1 1 1]);

hm=heatmap(ind_med,'Colormap',col6,'ColorLimits',[0 ceil(max_val6)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
%hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));
hm.FontSize=6;

%lateral heatmap
ax24=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col6,'ColorLimits',[0 ceil(max_val6)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
%hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));
hm.FontSize=6;

%labels
ax25=nexttile([nr 1]);
curr_lab=newlabels(par_seq==6);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',4.5);
set(ax25, 'color', 'none');
axis off;

linkaxes([ax22 ax25], 'y');


set(f,'renderer','Painters');


% --------------PLOT COLORBARS---------%
f2=figure(2);set(f2, 'Position', [488 633 560 129]);
til=tiledlayout(1,6, 'TileSpacing','none', 'Padding','tight');

ax1=nexttile; %CTX
colormap(ax1,col1);
cl=colorbar('AxisLocation', 'in', 'Location', 'west');
pos=get(cl, 'Position');
cl.TickLength=pos(3);
cl.Ticks = [0 1/3 2/3 1];
cl.TickLabels = {'0', [], [], num2str(ceil(max_val1))};
cl.YAxisLocation='right';
cl.Label.String="% input";
cl.Label.FontSize=8;
set(ax1, 'color', 'none');
%cl.FontSize=4;
axis off;

ax2=nexttile; %sCTX
colormap(ax2,col2);
cl=colorbar('AxisLocation', 'in', 'Location', 'west');
pos=get(cl, 'Position');
cl.TickLength=pos(3);
cl.Ticks = [0 1/3 2/3 1];
cl.TickLabels = {'0', [], [], num2str(ceil(max_val2))};
cl.YAxisLocation='right';
cl.Label.String="% input";
cl.Label.FontSize=8;
set(ax2, 'color', 'none');
%cl.FontSize=4;
axis off;

ax3=nexttile; %TH
colormap(ax3,col3);
cl=colorbar('AxisLocation', 'in', 'Location', 'west');
pos=get(cl, 'Position');
cl.TickLength=pos(3);
cl.Ticks = [0 1/3 2/3 1];
cl.TickLabels = {'0', [], [], num2str(ceil(max_val3))};
cl.YAxisLocation='right';
cl.Label.String="% input";
cl.Label.FontSize=8;
set(ax3, 'color', 'none');
%cl.FontSize=4;
axis off;

ax4=nexttile; %HYP
colormap(ax4,col4);
cl=colorbar('AxisLocation', 'in', 'Location', 'west');
pos=get(cl, 'Position');
cl.TickLength=pos(3);
cl.Ticks = [0 1/3 2/3 1];
cl.TickLabels = {'0', [], [], num2str(ceil(max_val4))};
cl.YAxisLocation='right';
cl.Label.String="% input";
cl.Label.FontSize=8;
set(ax4, 'color', 'none');
%cl.FontSize=4;
axis off;

ax5=nexttile; %MB
colormap(ax5,col5);
cl=colorbar('AxisLocation', 'in', 'Location', 'west');
pos=get(cl, 'Position');
cl.TickLength=pos(3);
cl.Ticks = [0 1/3 2/3 1];
cl.TickLabels = {'0', [], [], num2str(ceil(max_val5))};
cl.YAxisLocation='right';
cl.Label.String="% input";
cl.Label.FontSize=8;
set(ax5, 'color', 'none');
%cl.FontSize=4;
axis off;

ax6=nexttile; %HB
colormap(ax6,col6);
cl=colorbar('AxisLocation', 'in', 'Location', 'west');
pos=get(cl, 'Position');
cl.TickLength=pos(3);
cl.Ticks = [0 1/3 2/3 1];
cl.TickLabels = {'0', [], [], num2str(ceil(max_val6))};
cl.YAxisLocation='right';
cl.Label.String="% input";
cl.Label.FontSize=8;
set(ax6, 'color', 'none');
%cl.FontSize=4;
axis off;

set(f2,'renderer','Painters');


%% ------------------------------PLOT without colors-------------------------------%
par_seq=cell2mat(newlabels(:,3));
num_areas=length(unique(par_seq)); %number of parental areas
f=figure(3);set(f,'Position', [656.2000 97 445.6000 627.2000]);
til=tiledlayout(length(par_seq),8, 'TileSpacing','none', 'Padding','tight');
set(til, 'OuterPosition', [0 0 0.8 1]);

%-------CTX------%
nr=length(find(par_seq==1));
%bias index
ax1=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==1,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax1, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax1.YLabel.String="CTX";
ax1.YLabel.Color=[0.8784 0.5020 0.1490];
ax1.YLabel.FontSize=8;
ax1.XAxisLocation="top";
ax1.XLabel.String=["Bias"; "index"];
ax1.XLabel.FontSize=9;


%heatmaps
ind_med=(individuals_m(par_seq==1,:));
ind_lat=(individuals_l(par_seq==1,:));
max_val1=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax2=nexttile([nr 3]);
col1=custom_cmap('whitemagma', 300);
%col1=customcolormap([0 1], [0 0 0; 1 1 1]);

hm=heatmap(ind_med,'Colormap',col1,'ColorLimits',[0 ceil(max_val1)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));
hm.Title="medial";
set(struct(hm).Axes.Title,'FontWeight','normal', 'FontSize', 9);

%lateral heatmap
ax3=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col1,'ColorLimits',[0 ceil(max_val1)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));
hm.Title="lateral";
set(struct(hm).Axes.Title,'FontWeight','normal', 'FontSize', 9);

%labels
ax4=nexttile([nr 1]);
curr_lab=newlabels(par_seq==1);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',7);
set(ax4, 'color', 'none');CurYLim=ylim;
axis off;

linkaxes([ax1 ax4], 'y');




% %-------sCTX------%
% nr=length(find(par_seq==2));
% %bias index
% ax5=nexttile([nr 1]);
% CurXLim=xlim;
% coll=colors(par_seq==2,:);
% for i=1:nr
%     p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
%        'EdgeColor','none');
%     hold on;
% end
% set(ax5, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% % ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
% ax5.YLabel.String="sCTX";
% ax5.YLabel.Color=[0.4314    0.2118    0.0471];
% ax5.YLabel.FontSize=8;
% 
% %heatmaps
% ind_med=(individuals_m(par_seq==2,:));
% ind_lat=(individuals_l(par_seq==2,:));
% max_val2=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
% %medial heatmap
% ax6=nexttile([nr 3]);
% %col2=customcolormap([0 1], [0.4706 0.3137 0; 1 1 1]);
% 
% hm=heatmap(ind_med,'Colormap',col1,'ColorLimits',[0 ceil(max_val2)],'ColorbarVisible','off',...
%     'GridVisible','off', 'CellLabelColor','none'); 
% hm.XDisplayLabels = nan(size(hm.XDisplayData));
% hm.YDisplayLabels = nan(size(hm.YDisplayData));
% 
% %lateral heatmap
% ax7=nexttile([nr 3]);
% 
% hm=heatmap(ind_lat,'Colormap',col1,'ColorLimits',[0 ceil(max_val2)],'ColorbarVisible','off',...
%     'GridVisible','off', 'CellLabelColor','none'); 
% hm.XDisplayLabels = nan(size(hm.XDisplayData));
% hm.YDisplayLabels = nan(size(hm.YDisplayData));
% 
% %labels
% ax8=nexttile([nr 1]);
% curr_lab=newlabels(par_seq==2);
% text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',6);
% set(ax8, 'color', 'none');
% axis off;
% 
% linkaxes([ax5 ax8], 'y');
% 

%-------TH------%
nr=length(find(par_seq==3));
%bias index
ax9=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==3,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax9, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax9.YLabel.String="TH";
ax9.YLabel.Color=[0.8784    0.0353    0.1490];
ax9.YLabel.FontSize=8;
ylim([0 nr]);

%heatmaps
ind_med=(individuals_m(par_seq==3,:));
ind_lat=(individuals_l(par_seq==3,:));
max_val3=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax10=nexttile([nr 3]);
%col3=customcolormap([0 0.5 1], [0.5686 0 0.0745; 0.9216 0.2118 0.3059;1 1 1]);

hm=heatmap(ind_med,'Colormap',col1,'ColorLimits',[0 ceil(max_val3)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%lateral heatmap
ax11=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col1,'ColorLimits',[0 ceil(max_val3)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%labels
ax12=nexttile([nr 1]);
curr_lab=newlabels(par_seq==3);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',7);
set(ax12, 'color', 'none');
axis off;

linkaxes([ax9 ax12], 'y');


%-------HYP------%
nr=length(find(par_seq==4));
%bias index
ax13=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==4,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax13, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax13.YLabel.String="HYP";
ax13.YLabel.Color=[0.6392    0.0588    0.4471];
ax13.YLabel.FontSize=8;
ylim([0 nr]);

%heatmaps
ind_med=(individuals_m(par_seq==4,:));
ind_lat=(individuals_l(par_seq==4,:));
max_val4=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax15=nexttile([nr 3]);
%col4=customcolormap([0 0.5 1], [0.4314 0.0941 0.4314; 0.6902 0.2902 0.6902;1 1 1]);

hm=heatmap(ind_med,'Colormap',col1,'ColorLimits',[0 ceil(max_val4)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%lateral heatmap
ax16=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col1,'ColorLimits',[0 ceil(max_val4)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%labels
ax17=nexttile([nr 1]);
curr_lab=newlabels(par_seq==4);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',7);
set(ax17, 'color', 'none');
axis off;

linkaxes([ax13 ax17], 'y');

%-------MB------%
nr=length(find(par_seq==5));
%bias index
ax18=nexttile([nr 1]);
CurXLim=xlim;
coll=colors(par_seq==5,:);
for i=1:nr
    p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
       'EdgeColor','none');
    hold on;
end
set(ax18, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
% ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
ax18.YLabel.String="MB";
ax18.YLabel.Color=[0.196	0.533	0.741];
ax18.YLabel.FontSize=8;
ax18.XLabel.String="Animal";
ax18.XLabel.FontSize=8;
ylim([0 nr]);

%heatmaps
ind_med=(individuals_m(par_seq==5,:));
ind_lat=(individuals_l(par_seq==5,:));
max_val5=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
%medial heatmap
ax19=nexttile([nr 3]);
%col5=customcolormap([0 0.5 1], [0 0.3333 0.5412; 0.3059 0.6314 0.8314;1 1 1]);

hm=heatmap(ind_med,'Colormap',col1,'ColorLimits',[0 ceil(max_val5)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
% hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%lateral heatmap
ax20=nexttile([nr 3]);

hm=heatmap(ind_lat,'Colormap',col1,'ColorLimits',[0 ceil(max_val5)],'ColorbarVisible','off',...
    'GridVisible','off', 'CellLabelColor','none'); 
% hm.XDisplayLabels = nan(size(hm.XDisplayData));
hm.YDisplayLabels = nan(size(hm.YDisplayData));

%labels
ax21=nexttile([nr 1]);
curr_lab=newlabels(par_seq==5);
text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',7);
set(ax21, 'color', 'none');
axis off;

linkaxes([ax18 ax21], 'y');

%-------HB------%
nr=length(find(par_seq==6));
if nr>0
    %bias index
    ax22=nexttile([nr 1]);
    CurXLim=xlim;
    coll=colors(par_seq==6,:);
    for i=1:nr
        p=patch([0 CurXLim(2) CurXLim(2) 0], [(i-1) (i-1) (i) (i)],coll(i,:),...
           'EdgeColor','none');
        hold on;
    end
    set(ax22, 'YDir', 'reverse', 'color', 'none', 'XTick',[],'YTick',[],'Layer', 'top','box','on');
    % ax1.YAxis.Label.Visible='on'; ax1.YAxis.Label.Color='k';
    ax22.YLabel.String="HB";
    ax22.YLabel.Color=[0.1451    0.7216    0.0588];
    ax22.YLabel.FontSize=8;
    ax22.XLabel.String="Animal";
    ax22.XLabel.FontSize=6;
    
    %heatmaps
    ind_med=(individuals_m(par_seq==6,:));
    ind_lat=(individuals_l(par_seq==6,:));
    max_val6=max(max(ind_med,[],'all'), max(ind_lat,[],'all'));
    %medial heatmap
    ax23=nexttile([nr 3]);
    %col6=customcolormap([0 0.5 1], [0.0706 0.5020 0.0039; 0.3059 0.7608 0.2353;1 1 1]);
    
    hm=heatmap(ind_med,'Colormap',col1,'ColorLimits',[0 ceil(max_val6)],'ColorbarVisible','off',...
        'GridVisible','off', 'CellLabelColor','none'); 
    %hm.XDisplayLabels = nan(size(hm.XDisplayData));
    hm.YDisplayLabels = nan(size(hm.YDisplayData));
    hm.FontSize=6;
    
    %lateral heatmap
    ax24=nexttile([nr 3]);
    
    hm=heatmap(ind_lat,'Colormap',col1,'ColorLimits',[0 ceil(max_val6)],'ColorbarVisible','off',...
        'GridVisible','off', 'CellLabelColor','none'); 
    %hm.XDisplayLabels = nan(size(hm.XDisplayData));
    hm.YDisplayLabels = nan(size(hm.YDisplayData));
    hm.FontSize=6;
    
    %labels
    ax25=nexttile([nr 1]);
    curr_lab=newlabels(par_seq==6);
    text(repmat(0.06,1,nr),(1:nr)-0.5,flipud(curr_lab), 'HorizontalAlignment', 'left','FontSize',7);
    set(ax25, 'color', 'none');
    axis off;
    
    linkaxes([ax22 ax25], 'y');
end
%colorbar
colormap(gca,col1);
cl=colorbar('AxisLocation', 'in', 'Location', 'west', 'Position', [0.8 0.8 0.045 0.15]);
pos=get(cl, 'Position');
cl.TickLength=pos(3)/1.5;
cl.Ticks = [0 1/3 2/3 1];
cl.TickLabels = {'0', [], [], 'max'};
cl.YAxisLocation='right';
cl.Label.String="% input";
cl.Label.FontSize=8;
cl.Label.Position=[1.3 0.5 0];

par_areas=["CTX"; "TH"; "HYP"; "MB"];%; "HB"];
max_vals=ceil([max_val1; max_val3; max_val4; max_val5]); %max_val6]);
tbl=table(par_areas, max_vals)

set(f,'renderer','Painters');

