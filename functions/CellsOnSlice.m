%% Plot cell position overlaid with atlas slices

% 
% myType = flip(unique(OD.type),2);
% data.Lat = load(fullfile(spath,'Lat.mat'));
% data.Med = load(fullfile(spath,'Med.mat'));

clc;clear data; clear data2; clear ttd;clear ttm; clear OD;
%load(fullfile(spath,'OD.mat')); load(fullfile(spath,'Coord.mat'));
celltypes = {'grp'};
list_id = [];myexp_labels = {};nn = 0;n_Data=[];OD={}; Coord={};tt=0;
labels={};
for c=1:length(celltypes)
    
    OD{c}=load(fullfile(spath,['OD_',celltypes{c},'.mat']));
    Coord{c}=load(fullfile(spath,['Coord_',celltypes{c},'.mat']));
    
    myType = flip(unique(OD{c}.OD_cell.type),2);
    labels = BBB(:,[1,3,6,8]);
    descendants=BBB(:,4);
        
    for t = 1: length(myType)
        typeName = myType{t};   
        tableid = find(contains(OD{c}.OD_cell.type,typeName));
        expid= OD{c}.OD_cell.exp(tableid);
        coords=Coord{c}.Coord_cell.ori(tableid);
        list_id = [list_id tableid];
        tt=tt+1;

        %data based on brain regions and filtered
        data(tt).type = typeName;
        data(tt).table = load(fullfile(spath,[celltypes{c}, '_', typeName,'.mat']));
        thr = []; thr = data(tt).table.Threshold;

        %raw unfiltered wholebrain data
        data2(tt).type = typeName;
        data2(tt).table.tableid = tableid;
        data2(tt).table.expid = expid;
        data2(tt).table.Coords = cellfun(@(x) table2array(x(:, 3:5)), coords, 'UniformOutput', false);

        n_Data(c,t)=size(data(tt).table.CoordMerge,2);
        for j = 1:length(tableid)
            nn = nn+1;
            myexp_labels{nn,1} = [typeName,OD{c}.OD_cell.group{tableid(j)}];
            myexp_labels{nn,2} = OD{c}.OD_cell.exp{tableid(j)};
        end
    end
end
set(0,'defaultAxesFontSize',10);



bregma = allenCCFbregma; %get bregma coordinates %AP, ML, DV
ML = [3, 2.64, 2.28, 1.92, 1.2, 0.84]; %ML=ML+0.1;
slices_sag = bregma(2)+ML*100; %sliceLabel = {'3.10','2.74','2.38','2.02','1.30','0.94'};
% AP = [-0.9, -1.26, -1.62, -1.98, -2.34, -2.7]; %ZI
%AP = [1.6, 1.16, 0.72, 0.28, -0.16, -0.6]; %ACA
%AP = [-2, -2.5, -3, -3.5, -4, -4.5];
%AP = [-2.8, -2.9, -3, -3.1, -3.2, -3.3];
%AP = [-4.4 -4.55 -4.7 -4.85 -5.0 -5.15];
AP=[0, -1, -2, -3, -4, -5];
DV= [0 -1 -2 -3 -4 -5];
slices_dv=1-DV*100;
slices_cor=bregma(1)-AP*100;
%slice=squeeze(av(940,:,:));%choose slice to plot
%areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl"];
areas=[];
%areas=[];

%colors
colors={};
colors{1,1}=cbrewer('seq','RdPu',8);
colors{1,2}=cbrewer('seq','BuGn',8);
f=figure(1);
figName = 'cells_on_slice';
til=tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');
% til.Title.String='Summary input bias percentage';
% til.Title.FontSize=12;
% til.Title.FontWeight='bold';

for sid = 1:length(slices_cor)
    sa=nexttile;
    s=slices_cor(sid);
    slice=squeeze(av(s,:,:));
    prettyslice(slice,st,gca); hold on;
    %set(gca, 'Visible', 'off');
    for t = 1:length(myType)
        dt=data(t).table;
        typeName = myType{t};
        dt.TT_title;
        labels = dt.TT(:,1:3);
        threCur= dt.Threshold;
        myCoord = dt.CoordMerge;
        % myCoord = dt.Coords;
        if ~isempty(areas)
            id_ar = find(ismember(labels(:,1),areas));
            myCoord=myCoord(id_ar,:);
        end
        myCoord(strcmp(labels(:,1),'SCs')==1,:)={[]}; %remove cells in SCs
        myCoord(ismember(labels(:,1),["MBsen", "MBmot", "MBsta"]),:)={[]}; %remove cells in MB
        myCoord(ismember(labels(:,1),["MH", "LH", "CL"]),:)={[]}; %remove likely inputs from iSC dSC starters

            for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
                if threCur(i,1)==0
                    myCoord(i,:)={[]};
                end
            end
        %-----plot transformed coordinates
        [~, sl] = size(myCoord);coordXYZ = [];
        odd=[0 2 4 6 8 10 12];
        for i = 1:sl
            coordXYZ = cell2mat(myCoord(:,i));
            coordXYZ=coordXYZ+randn(size(coordXYZ))/500;%add some small random noise to visualize scatter points that would be on the same pixel
            coordXYZ=coordXYZ(coordXYZ(:,1)>=-AP(sid)-0.05 & coordXYZ(:,1)<=-AP(sid)+0.05,:); %keep cells that are +-50 nm of the slice
            s=scatter(bregma(3)-coordXYZ(:,3)*100, bregma(2)+coordXYZ(:,2)*100, 5, 'filled',...
                CData=col(t,:), MarkerFaceAlpha=0.7);
            if t==2
                %uistack(s, 'down', odd(sl));
                uistack(s, 'down',2);
            end
            hold on;
        end
    end
    set(sa,'Box','off','Color','none','XColor', 'none','YColor','none');
    xlabel(sa,"AP "+num2str(AP(sid),'%.2f'), 'Color', 'k');
%     xlim([300 550]); ylim([390 580]); %ZI
    %xlim([420 575]); ylim([90 320]); %ZI
end
%set(gca, 'Visible', 'off');
set(gcf,'renderer','Painters', 'WindowState', 'maximized');

%% Plot med and lat independently
clc;clf;clear data;clear ttd;clear ttm;
labels = BBB(:,[1,3,5]); descendants=BBB(:,4);
load(fullfile(spath,'OD.mat'));
load(fullfile(spath,'Coord.mat'));
myType = flip(unique(OD.type),2);
data.Lat = load(fullfile(spath,'Lat.mat'));
data.Med = load(fullfile(spath,'Med.mat'));


bregma = allenCCFbregma; %get bregma coordinates %AP, ML, DV
ML = [3, 2.64, 2.28, 1.92, 1.2, 0.84]; %ML=ML+0.1;
slices_sag = bregma(2)+ML*100; %sliceLabel = {'3.10','2.74','2.38','2.02','1.30','0.94'};
AP = [0, -1, -2, -3, -4, -5];
%AP = [-2, -2.5, -3, -3.5, -4, -4.5];
%AP = [-2.8, -2.9, -3, -3.1, -3.2, -3.3];
slices_cor=bregma(1)-AP*100;
slice=squeeze(av(940,:,:));%choose slice to plot
%areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl"];
areas=[];

%colors
colors=cbrewer('qual', 'Set1',9);
colors=colors([1:5, 7:8], :);

% til.Title.String='Summary input bias percentage';
% til.Title.FontSize=12;
% til.Title.FontWeight='bold';

for t = 1:length(myType)
        f=figure(t);
        figName = 'cells_on_slice';
        til=tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');
        dt=data.(myType{t});
        typeName = myType{t};
        dt.TT_title;
        labels = dt.TT(:,1:3);
        threCur= dt.Threshold;
        myCoord = dt.CoordMerge;
        if ~isempty(areas)
            id_ar = find(ismember(labels(:,1),areas));
            myCoord=myCoord(id_ar,:);
        end
        myCoord(strcmp(labels(:,1),'SCs')==1,:)={[]}; %remove cells in SCs

        %     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
        %         if threCur(i,1)==0
        %             myCoord(i,:)={[]};
        %         end
        %     end
        %-----plot transformed coordinates
        [~, sl] = size(myCoord);coordXYZ = [];
        for sid = 1:length(slices_cor)
            sa=nexttile;
            s=slices_cor(sid);
            slice=squeeze(av(s,:,:));
            prettyslice(slice,st,gca); hold on;
            %set(gca, 'Visible', 'off');

            for i = 1:sl
                coordXYZ = cell2mat(myCoord(:,i));
                coordXYZ=coordXYZ+randn(size(coordXYZ))/500;%add some small random noise to visualize scatter points that would be on the same pixel
                coordXYZ=coordXYZ(coordXYZ(:,1)>=-AP(sid)-0.05 & coordXYZ(:,1)<=-AP(sid)+0.05,:); %keep cells that are +-50 nm of the slice
                s=scatter(bregma(2)-coordXYZ(:,3)*100, -bregma(3)+coordXYZ(:,2)*100, 10, 'filled',...
                    CData=colors(i,:), MarkerFaceAlpha=0.6);
                if t==2
                    %uistack(s, 'down', odd(sl));
                    uistack(s, 'down',2);
                end
                hold on;
            end
            set(sa,'Box','off','Color','none','XColor', 'none','YColor','none');
            xlabel(sa,"AP "+num2str(AP(sid),'%.2f'), 'Color', 'k');
        end
        set(gcf,'renderer','Painters', 'WindowState', 'maximized');
end
%set(gca, 'Visible', 'off');


%% plot horizontal slices

clc;clf;clear data;clear ttd;clear ttm;
labels = BBB(:,[1,3,5]); descendants=BBB(:,4);
load(fullfile(spath,'OD.mat'));
load(fullfile(spath,'Coord.mat'));
myType = flip(unique(OD.type),2);
data.Lat = load(fullfile(spath,'Lat.mat'));
data.Med = load(fullfile(spath,'Med.mat'));


bregma = allenCCFbregma; %get bregma coordinates %AP, ML, DV
ML = [3, 2.64, 2.28, 1.92, 1.2, 0.84]; %ML=ML+0.1;
slices_sag = bregma(2)+ML*100; %sliceLabel = {'3.10','2.74','2.38','2.02','1.30','0.94'};
AP = [0, -1, -2, -3, -4, -5];
%AP = [-2, -2.5, -3, -3.5, -4, -4.5];
%AP = [-2.8, -2.9, -3, -3.1, -3.2, -3.3];
% slices_cor=bregma(1)-AP*100;
% slice=squeeze(av(940,:,:));%choose slice to plot
DV=[0.7 1.1 1.5 1.9 2.3 2.7];
slices_cor=bregma(3)+DV*100;
areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl"];
%areas=[];

%colors
colors={};
colors{1,1}=cbrewer('seq','RdPu',8);
colors{1,2}=cbrewer('seq','BuGn',8);
f=figure(1);
figName = 'cells_on_slice';
til=tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');
% til.Title.String='Summary input bias percentage';
% til.Title.FontSize=12;
% til.Title.FontWeight='bold';

for sid = 1:length(slices_cor)
    sa=nexttile;
    s=round(slices_cor(sid));
    slice=squeeze(av(:,s,:));
    prettyslice(slice,st,gca); hold on;
    %set(gca, 'Visible', 'off');
    for t = 1:length(myType)
        dt=data.(myType{t});
        typeName = myType{t};
        dt.TT_title;
        labels = dt.TT(:,1:3);
        threCur= dt.Threshold;
        myCoord = dt.CoordMerge;
        if ~isempty(areas)
            id_ar = find(ismember(labels(:,1),areas));
            myCoord=myCoord(id_ar,:);
        end
        myCoord(strcmp(labels(:,1),'SCs')==1,:)={[]}; %remove cells in SCs

        %     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
        %         if threCur(i,1)==0
        %             myCoord(i,:)={[]};
        %         end
        %     end
        %-----plot transformed coordinates
        [~, sl] = size(myCoord);coordXYZ = [];
        odd=[0 2 4 6 8 10 12];
        for i = 1:sl
            coordXYZ = cell2mat(myCoord(:,i));
            coordXYZ=coordXYZ+randn(size(coordXYZ))/500;%add some small random noise to visualize scatter points that would be on the same pixel
            %coordXYZ(:,1)=coordXYZ(:,1)+(-0.05 + (0.05+0.05).*rand(length(coordXYZ),1)); %jitter points along AP axis to cover slice spacing effect
            coordXYZ=coordXYZ(coordXYZ(:,2)>=DV(sid)-0.05 & coordXYZ(:,2)<=DV(sid)+0.05 & coordXYZ(:,3)>=0,:); %keep cells that are +-50 nm of the slice
            s=scatter(bregma(2)+coordXYZ(:,3)*100, bregma(1)+coordXYZ(:,1)*100, 5, 'filled',...
                CData=col(t+2,:), MarkerFaceAlpha=0.6);
            if t==2
                %uistack(s, 'down', odd(sl));
                uistack(s, 'down',2);
            end
            hold on;
        end
    end
    set(sa,'Box','off','Color','none','XColor', 'none','YColor','none');
    xlabel(sa,"DV "+num2str(DV(sid),'%.2f'), 'Color', 'k');
end
%set(gca, 'Visible', 'off');
set(gcf,'renderer','Painters', 'WindowState', 'maximized');

%% with rotation (CAN'T ALIGN SCATTER WITH SLICE POSITION)
clc;clf;clear data;clear ttd;clear ttm;
labels = BBB(:,[1,3,5]); descendants=BBB(:,4);
load(fullfile(spath,'OD.mat'));
load(fullfile(spath,'Coord.mat'));
myType = flip(unique(OD.type),2);
data.Lat = load(fullfile(spath,'Lat.mat'));
data.Med = load(fullfile(spath,'Med.mat'));


bregma = allenCCFbregma; %get bregma coordinates %AP, ML, DV
ML = [3, 2.64, 2.28, 1.92, 1.2, 0.84]; %ML=ML+0.1;
slices_sag = bregma(2)+ML*100; %sliceLabel = {'3.10','2.74','2.38','2.02','1.30','0.94'};
AP = [0, -1, -2, -3, -4, -5];
%AP = [-2, -2.5, -3, -3.5, -4, -4.5];
%AP = [-2.8, -2.9, -3, -3.1, -3.2, -3.3];
% slices_cor=bregma(1)-AP*100;
% slice=squeeze(av(940,:,:));%choose slice to plot
DV=[0.5 1.5 2.5 3.5 4.5 5.5];
%DV=[3 3.2 3.4 3.6 3.8 4];

%areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl"];
areas=[];

volume=squeeze(av);
%volshow(volume);
theta = [0, -28.022, -8.29];
transl = [0 0 0];
tform = rigidtform3d(theta,transl);
VolumeRotated = imwarp(volume,tform);


slices_cor=bregma(3)+DV*100;
volshow(VolumeRotated);
volshow(subAtlas.mar.D)
%stats = regionprops3(subAtlas.mar.D,'all')
%colors
colors={};
colors{1,1}=cbrewer('seq','RdPu',8);
colors{1,2}=cbrewer('seq','BuGn',8);
f=figure(1);
figName = 'cells_on_slice';
til=tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');


for sid = 1:length(slices_cor)
    sa=nexttile;
    s=round(slices_cor(sid));
    slice=squeeze(VolumeRotated(:,s,:));
    prettyslice(slice,st,gca); hold on;
    for t = 1:length(myType)
        dt=data.(myType{t});
        typeName = myType{t};
        dt.TT_title;
        labels = dt.TT(:,1:3);
        threCur= dt.Threshold;
        myCoord = dt.CoordMerge;
        if ~isempty(areas)
            id_ar = find(ismember(labels(:,1),areas));
            myCoord=myCoord(id_ar,:);
        end
        myCoord(strcmp(labels(:,1),'SCs')==1,:)={[]}; %remove cells in SCs

        %     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
        %         if threCur(i,1)==0
        %             myCoord(i,:)={[]};
        %         end
        %     end
        %-----plot transformed coordinates
        [~, sl] = size(myCoord);coordXYZ = [];
        odd=[0 2 4 6 8 10 12];
        for i = 1:sl
            coordXYZ = cell2mat(myCoord(:,i));
            coordXYZ=coordXYZ+randn(size(coordXYZ))/500;%add some small random noise to visualize scatter points that would be on the same pixel
            X=coordXYZ(:,1);
            Y=coordXYZ(:,3);
            Z=coordXYZ(:,2);
            Xnew=X*cos(phi)-Y*sin(theta2)*sin(phi)-Z*cos(theta2)*sin(phi);
            Ynew=Y*cos(theta2)-Z*sin(theta2);
            Znew=X*sin(phi)+Y*sin(theta2)*cos(phi)+Z*cos(theta2)*cos(phi);
            coordXYZ=[Xnew Znew Ynew];
            coordXYZ=coordXYZ(coordXYZ(:,2)>=DV(sid)-0.05 & coordXYZ(:,2)<=DV(sid)+0.05,:); %keep cells that are +-50 nm of the slice
            s=scatter(bregma(2)+coordXYZ(:,3)*100, bregma(1)+coordXYZ(:,1)*100, 5, 'filled',...
                CData=col(t+2,:), MarkerFaceAlpha=0.6);
            if t==2
                %uistack(s, 'down', odd(sl));
                uistack(s, 'down',2);
            end
            hold on;
        end
    end
    set(sa,'Box','off','Color','none','XColor', 'none','YColor','none');
    xlabel(sa,"DV "+num2str(DV(sid),'%.2f'), 'Color', 'k');
end
%set(gca, 'Visible', 'off');
set(gcf,'renderer','Painters', 'WindowState', 'maximized');

%% Calculate brain area volumes


stats_all = regionprops3(av(~ismember(av,1)),'Volume');

stats_all.Acronyms=st.acronym(1:size(stats_all,1));


