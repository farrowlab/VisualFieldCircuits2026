%% --- --- prepare data -

clc;clear data;clear ttd;clear ttm; clear OD;
%load(fullfile(spath,'OD.mat')); load(fullfile(spath,'Coord.mat'));
celltypes = {'grp'};
list_id = [];myexp_labels = {};nn = 0;n_Data=[];OD={}; Coord={};tt=0;
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
        end
    end
end
set(0,'defaultAxesFontSize',10);
%%
n_Data=reshape(n_Data.', 1, []);
countLat={};densityLat={};percentLat={};threshlLat={};errPerLat={};errDenLat={};
countMed={};densityMed={};percentMed={};threshlMed={};errPerMed={};errDenMed={};
labels = BBB(:,[1,3,6]); ll = length(BBB);
% read Lat & Med data
for t = 1:2 %2 if 1 celltype, 4 if 2 celltypes
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
areas=["ZI"];

idSC = find(~ismember(labels(:,1),areas));

%idSC = find(startsWith(labels(:,1),'SCs'));%removes all SC layers from analysis
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
    thresT2 = 0;  nanIDX2 = find((dmLat{c} < thresT2) & (dmMed{c} < thresT2)); %normalized_density > thresT2
    pmLat{c} = mean(percentLat{c},2);pmMed{c} = mean(percentMed{c},2);
    thresT3 = 0; nanIDX3 = find((pmLat{c} < thresT3) & (pmMed{c} < thresT3)); %percentage > thresT3
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
    disp(['Total number of included ', celltypes{c}, ' medial input cells: ' ,num2str(sum(countMed{c},'all'))]);
    disp(['Total number of included ', celltypes{c},' lateral input cells: ' ,num2str(sum(countLat{c},'all'))]);
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



%Med + Lat
%abbr = {};


labelname = '% of total inputs';
property = 'percentage';
d_labelname = '% of sum of densities';
d_property = 'density';
xx={}; yy={}; errY={}; errX={}; individuals_m={}; individuals_l={};
d_xx={}; d_yy={}; d_errY={}; d_errX={}; d_individuals_m={}; d_individuals_l={};
for c=1:length(celltypes)
%percentage 
    xx{c} = mean(percentMed{c},2);
    yy{c} = mean(percentLat{c},2);
    errY{c} = errPerLat{c};errX{c} = errPerMed{c};
    individuals_m{c} = percentMed{c};
    individuals_l{c} = percentLat{c};
    
    %density 
    d_xx{c} = mean(densityMed{c},2);
    d_yy{c} = mean(densityLat{c},2);
    d_errY{c} = errDenLat{c};d_errX{c} = errDenMed{c};
    d_individuals_m{c} = densityMed{c};
    d_individuals_l{c} = densityLat{c};
       
end

zz=xx{1}+yy{1};
plot_seq=[12 13 8 7 2 3 4 1 6 9 14 11 10 5]';
% [ss, plot_seq] = sortrows(zz, 'descend');
newlabels = labels(plot_seq);


for c=1:length(celltypes)
    xx{c}=xx{c}(plot_seq);yy{c}=yy{c}(plot_seq); d_xx{c}=d_xx{c}(plot_seq); d_yy{c}=d_yy{c}(plot_seq);
    errY{c}=errY{c}(plot_seq);errX{c}=errX{c}(plot_seq); d_errY{c}=d_errY{c}(plot_seq);d_errX{c}=d_errX{c}(plot_seq);
    individuals_m{c}=individuals_m{c}(plot_seq,:);individuals_l{c}=individuals_l{c}(plot_seq,:);
    d_individuals_m{c}=d_individuals_m{c}(plot_seq,:);d_individuals_l{c}=d_individuals_l{c}(plot_seq,:);
end
clearvars ss zz;


nGroups = 2;groupwidth = min(0.7, nGroups/(nGroups + 1));
lr = length(yy{1});
%h = bar(1:lr,[xx yy],'FaceColor','flat','BarWidth',0.98);hold on;
%group data so that you can access medial and lateral with for loop
%x positions for medial and lateral
x={};
for k = 1:nGroups
    x{1,k} = (1:lr) - groupwidth/2 + (2*k-1) * groupwidth/(2*nGroups);
end

%percentage
ydata=[xx yy];%y mean positions for medial and lateral
errors=[errX errY];
individuals=cell(1, nGroups); individuals(1:2:end)=individuals_m; individuals(2:2:end)=individuals_l;

%density
d_ydata=[d_xx d_yy];%y mean positions for medial and lateral
d_errors=[d_errX d_errY];
d_individuals=[d_individuals_m d_individuals_l];




% plot med vs lat first
figure;set(gcf,'position',[606.6000 207.4000 229.7333 398.4000]);
clf;clc;
swarm={};%to plot legends in the right order
for i = 1:2 %first plot med then lat
    %typeName = myType{i};
    %plot means
%     h = scatter(x{1,i},ydata(:,i)',40, '_','LineWidth',1.5);
%     h.CData = [0.3 0.3 0.3];
%     uistack(h, 'top');
%     hold on;
    %plot individual experiment data as scatter dots
    coll = col(i,:);
    for jj=1:lr
        ind=individuals{1,i}(jj,:);
        ind=ind(:,not(isoutlier(individuals{1,i}(jj,:),'quartiles')));
        %for plotting the whiskers manually
        h1=plot(repmat(x{1,i}(jj),2,1)', [min(ind') prctile(individuals{1,i}(jj,:)',25)], "Color",[0.4 0.4 0.4],'LineWidth', 2.5);
        hold on;
        h2=plot(repmat(x{1,i}(jj),2,1)', [max(ind') prctile(individuals{1,i}(jj,:)',75)],"Color",[0.4 0.4 0.4],'LineWidth', 2.5);
        hold on;
        %plot boxcharts
        s=boxchart(repmat(x{1,i}(jj),size(individuals{1,i},2),1),individuals{1,i}(jj,:),'BoxFaceColor',coll, 'BoxWidth',0.3,...
            'BoxMedianLineColor','none', 'LineWidth', 1.5,'MarkerColor',[0.7 0.7 0.7], 'MarkerSize', 4,...
            'WhiskerLineStyle','none', 'BoxFaceAlpha',0.35, 'BoxEdgeColor','none');
        uistack(s, 'down');
        hold on;
        %plot swarmchart
        swarm{1,i}=swarmchart(repmat(x{1,i}(jj),size(ind,2),1),ind', 25,...
            'filled','MarkerFaceColor',coll);%plot swarmchart
        swarm{1,i}.XJitterWidth=0.1;
        uistack(swarm{1,i}, 'bottom');
        hold on;
        %plot median lines a bit wider than boxcharts
        h3=scatter(x{1,i}(jj), median(individuals{1,i}(jj,:)), 450,'_','LineWidth',2.0);
        h3.CData = [0.2 0.2 0.2];
        uistack(h3, 'top');
        hold on;
    end
    
end
%plot xgrid lines and shade every other area in between
CurYLim = ylim;
% for jj=lr
%     rec=rectangle('Position',[(jj-0.5) 0 1 CurYLim(2)], 'EdgeColor','none','FaceColor',[0.9 0.9 0.9]);
%     uistack(rec, 'bottom');
%     hold on;
% end
%fill(xgrid(1:2:end), xgrid(2:2:end));
%for 1:length(xgrid{1,1})

%plot lines connecting scatter dots
% for jj=1:lr
%     xvector=[repmat(x{1,1}(jj),size(individuals{1,1},2),1) repmat(x{1,2}(jj),size(individuals{1,1},2),1)];
%     yvector=[individuals{1,1}(jj,:)' individuals{1,2}(jj,:)'];
%     for i=1:size(individuals{1,1},2)
%         l= plot([xvector(i,1) xvector(i,2)],[yvector(i,1) yvector(i,2)],'LineWidth',1.2,'color',[0.6 0.6 0.6]);
%         uistack(l, 'bottom'); %draws the lines behind all the other graph elements
%     end
% end
% 
% %plot sig stars
% for jj=1:lr
%     xticks = jj; yl=max(max(individuals{1,1}(jj,:)),max(individuals{1,2}(jj,:)))+2;
%     [f, pval] = ttest2(individuals{1,1}(jj,:), individuals{1,2}(jj,:));
%     %sigstars on top
%     mysigstar(gca, xticks, yl, pval);
% end


% legend([swarm{1,1} swarm{1,2}],'Medial', 'Lateral', 'AutoUpdate','off');
% get the current Y limits
ylim([0 CurYLim(2)]);

hold off;
set(gca,'XTick',(1:1:lr),'TickLength',[0.0105 0.0105],'TickLabelInterpreter','none',...
    'XTickLabel',newlabels(:,1),'XTickLabelRotation',0,'FontSize',10,'FontName', 'Helvetica',...
    'GridColor',[.5 .5 .5],'GridLineStyle', '-.','XGrid', 'off','color', 'none', 'TickDir', 'out')
lab=ylabel ('Percentage of total inputs','FontName', 'Helvetica');new_lab=get(lab, 'Position');
% lab.Position(1)=new_lab(1)-0.15;
box off;

set(gcf,'renderer','Painters');set(gca,'color','none');

legend([swarm{1,1} swarm{1,2}],'medial', 'lateral')
% print(gcf,(fullfile(savePath,strcat('summary_med+lat','.eps'))),'-depsc','-r300');


%% plot total percentage
figure;set(gcf,'position',[606.6000 207.4000 439.2000 398.4000]);
clf;clc;
swarm={};%to plot legends in the right order
% max_len_data=max(size(individuals_m,2),size(individuals_l,2));
% min_len_data=min(size(individuals_m,2),size(individuals_l,2));
% individuals=[];
individuals={individuals_m individuals_l};
individuals=(individuals{1,1}+individuals{1,2})./2;

% for i=1:max_len_data
%     if i<=min_len_data
%         individuals(:,i)=(individuals_m(:,i)+individuals_l(:,i))./2;
%     elseif size(individuals_m,2)>size(individuals_l,2)
%         individuals(:,i)=individuals_m(:,i);
%     else
%         individuals(:,i)=individuals_l(:,i);
%     end
% end

    %typeName = myType{i};
    %plot means
%     h = scatter(x{1,i},ydata(:,i)',40, '_','LineWidth',1.5);
%     h.CData = [0.3 0.3 0.3];
%     uistack(h, 'top');
%     hold on;
%plot individual experiment data as scatter dots
coll = [0.17 0.43 0.73];
for jj=1:lr
%     if jj==lr %make the last boxchart red
%         coll=[0.69 0.22 0.22];
%     end
    ind=individuals(jj,:);
    ind=ind(:,not(isoutlier(individuals(jj,:),'quartiles')));
    %for plotting the whiskers manually
    h1=plot(repmat(jj,2,1)', [min(ind') prctile(individuals(jj,:)',25)], "Color",[0.4 0.4 0.4],'LineWidth', 2.5);
    hold on;
    h2=plot(repmat(jj,2,1)', [max(ind') prctile(individuals(jj,:)',75)],"Color",[0.4 0.4 0.4],'LineWidth', 2.5);
    hold on;
    %plot boxcharts
    s=boxchart(repmat(jj,size(individuals,2),1),individuals(jj,:),'BoxFaceColor',coll, 'BoxWidth',0.3,...
        'BoxMedianLineColor','none', 'LineWidth', 1.5,'MarkerColor',[0.7 0.7 0.7], 'MarkerSize', 4,...
        'WhiskerLineStyle','none', 'BoxFaceAlpha',0.35, 'BoxEdgeColor','none');
    uistack(s, 'down');
    hold on;
    %plot swarmchart
    swarm=swarmchart(repmat(jj,size(ind,2),1),ind', 25,...
        'filled','MarkerFaceColor',coll);%plot swarmchart
    swarm.XJitterWidth=0.1;
    uistack(swarm, 'bottom');
    hold on;
    %plot median lines a bit wider than boxcharts
    h3=scatter(jj, median(individuals(jj,:)), 450,'_','LineWidth',2.0);
    h3.CData = [0.2 0.2 0.2];
    uistack(h3, 'top');
    hold on;
end


%plot xgrid lines and shade every other area in between
CurYLim = ylim;
% for jj=lr
%     rec=rectangle('Position',[(jj-0.5) 0 1 CurYLim(2)], 'EdgeColor','none','FaceColor',[0.9 0.9 0.9]);
%     uistack(rec, 'bottom');
%     hold on;
% end
%fill(xgrid(1:2:end), xgrid(2:2:end));
%for 1:length(xgrid{1,1})

%plot lines connecting scatter dots
% for jj=1:lr
%     xvector=[repmat(x{1,1}(jj),size(individuals{1,1},2),1) repmat(x{1,2}(jj),size(individuals{1,1},2),1)];
%     yvector=[individuals{1,1}(jj,:)' individuals{1,2}(jj,:)'];
%     for i=1:size(individuals{1,1},2)
%         l= plot([xvector(i,1) xvector(i,2)],[yvector(i,1) yvector(i,2)],'LineWidth',1.2,'color',[0.6 0.6 0.6]);
%         uistack(l, 'bottom'); %draws the lines behind all the other graph elements
%     end
% end

%plot sig stars
% for jj=1:lr
%     xticks = jj; yl=max(max(individuals{1,1}(jj,:)),max(individuals{1,2}(jj,:)))+2;
%     [f, pval] = ttest2(individuals{1,1}(jj,:), individuals{1,2}(jj,:));
%     %sigstars on top
%     mysigstar(gca, xticks, yl, pval);
% end


%legend([swarm{1,1} swarm{1,2}],'Medial', 'Lateral', 'AutoUpdate','off');
% get the current Y limits
ylim([0 CurYLim(2)]);

hold off;
set(gca,'XTick',(1:1:lr),'TickLength',[0.0105 0.0105],'TickLabelInterpreter','none',...
    'XTickLabel',newlabels(:,1),'XTickLabelRotation',0,'FontSize',10,'FontName', 'Helvetica',...
    'GridColor',[.5 .5 .5],'GridLineStyle', '-.','XGrid', 'off','color', 'none', 'TickDir', 'out')
lab=ylabel ('Percentage of total inputs','FontName', 'Helvetica');new_lab=get(lab, 'Position');
% lab.Position(1)=new_lab(1)-0.15;
box off;

set(gcf,'renderer','Painters');set(gca,'color','none');
print(gcf,(fullfile(savePath,strcat('summary_med+lat','.eps'))),'-depsc','-r300');

%% calculate colormap to color vis cortices by percentage
individuals={individuals_m individuals_l};
individuals=(individuals{1,1}+individuals{1,2})./2;
ind_mean=median(individuals,2);
norm_mean=ind_mean;%/max(ind_mean);
colors=cbrewer('seq', 'YlOrRd',300);
lin=linspace(0,ceil(ind_mean(2)),300);
rv_colors={};
for i=1:length(norm_mean)-1
    for j=1:length(lin)-1
        if (norm_mean(i+1) >= lin(j)) && (norm_mean(i+1) <= lin(j+1))
            rv_colors{i}=colors(j,:);
        else
            continue
        end
    end
end
colormap(colors);
for i=1:length(norm_mean)-1
    rec=rectangle('Position',[i-0.5 0 0.5 0.5], 'EdgeColor','none','FaceColor',rv_colors{1,i}); hold on;
end
set(gca, 'XTick',(1:1:lr),'XTickLabel',newlabels(2:end,1));
cl = colorbar;
cl.Ticks = [0 0.2 0.4 0.6 0.8 1];
max_val=ceil(ind_mean(2));
nn=max_val/5;
cl.TickLabels={'0', num2str(nn), num2str(nn*2), num2str(nn*3), num2str(nn*4) ,num2str(max_val)};
