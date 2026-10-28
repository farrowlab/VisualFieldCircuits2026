%% rotate 3D coordinates by two angles
clc;
theta=deg2rad(28.022);
phi=deg2rad(8.29);
x1=-70; y1=249; z1=-35941;
%compute new coordinates
y2=y1*cos(theta)-z1*sin(theta);
x2=x1*cos(phi)-y1*sin(theta)*sin(phi)-z1*cos(theta)*sin(phi);
z2=x1*sin(phi)+y1*sin(theta)*cos(phi)+z1*cos(theta)*cos(phi);
disp(['x2= ',num2str(x2)]); disp(['y2= ',num2str(y2)]); disp(['z2= ',num2str(z2)]);


%% transform brainrender coordinates to MATLAB-generated coordinates
clc;
X=7830; Y=5694; Z=4296;
x=(X-5400)./(-1000);
y=(Y-5700)./1000;
z=Z./1000;
disp(['x= ',num2str(x)]); disp(['y= ',num2str(y)]); disp(['z= ',num2str(z)]);

%% Plot visual cortex scatter before and after transformation
clc;clf;clear data;clear ttd;clear ttm;
load(fullfile(spath,'OD.mat'));
load(fullfile(spath,'Coord.mat'));
myType = flip(unique(OD.type),2);
Coordinates_tot1=[];
Coordinates_tot2=[];
set(gcf, 'WindowState','maximized');
for t = 1:length(myType)
    coll = col(t,:);
    typeName = myType{t}; 
    dt = load(fullfile(spath,[typeName,'.mat']));
    dt.TT_title;
    labels = dt.TT(:,1:3);
    vis_areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl", "RSPd", "RSPagl"];
    idSC = find(contains(labels(:,1),vis_areas));
    threCur= dt.Threshold;
    myCoord = dt.CoordMerge;
    myCoord=myCoord(idSC,:);
%     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
%         if threCur(i,1)==0
%             myCoord(i,:)={[]};
%         end
%     end
    [sa, sl] = size(myCoord);coordXYZ = [];
    hAx(1)=subplot(1,2,1);   
    for i = 1:sl
        coordXYZ = cell2mat(myCoord(:,i));
        coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:);
        Coordinates_tot1=[Coordinates_tot1; coordXYZ];
        s=scatter(-coordXYZ(:,3), -coordXYZ(:,1), 10, 'filled', CData=coll, MarkerFaceAlpha=0.3);
        hold on;
    end
    
    xlabel ('ML (mm)'); ylabel ('AP (mm)');
    %plot with transformed coordinates
    hAx(2)=subplot(1,2,2);
    theta=deg2rad(-28.022);
    phi=deg2rad(-8.29);
    for i = 1:sl
        coordXYZ = cell2mat(myCoord(:,i));
        coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:);
        X=coordXYZ(:,1); 
        Y=coordXYZ(:,3); 
        Z=coordXYZ(:,2);
        Xnew=X*cos(phi)-Y*sin(theta)*sin(phi)-Z*cos(theta)*sin(phi);
        Ynew=Y*cos(theta)-Z*sin(theta);
        Znew=X*sin(phi)+Y*sin(theta)*cos(phi)+Z*cos(theta)*cos(phi);
        coordXYZ=[Xnew Znew Ynew];
        Coordinates_tot2=[Coordinates_tot2; coordXYZ];
        s=scatter(-coordXYZ(:,3), -coordXYZ(:,1), 10, 'filled', CData=coll, MarkerFaceAlpha=0.3);
        hold on;
    end
    
end

axis(hAx, 'equal');
%linkaxes(hAx)
xlabel ('ML (mm)'); ylabel ('AP (mm)');
xlim(hAx(1),[-(max(Coordinates_tot2(:,3))+0.3) -(min(Coordinates_tot2(:,3))-0.7)]); 
ylim(hAx(1),[-(max(Coordinates_tot2(:,1))+0.3) -(min(Coordinates_tot2(:,1))-0.3)]);
xlim(hAx(2),[-(max(Coordinates_tot2(:,3))+0.3) -(min(Coordinates_tot2(:,3))-0.7)]); 
ylim(hAx(2),[-(max(Coordinates_tot2(:,1))+0.3) -(min(Coordinates_tot2(:,1))-0.3)]);
set(hAx,'color','none','box','off');
% set(hAx(1), 'Position', [0.05, 0.69, 0.4, 0.4]);
% set(hAx(2), 'Position', [0.05, 0.37, 0.92, 0.27]);
%set(gcf, "Position", [488 105 391.4000 657]);

%% Plot visual cortices with contours
clc;clf;clear data;clear ttd;clear ttm;
celltypes = {'ntsr'};
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
Coordinates_tot=[];
f=figure(1);
type_len=[];
ZZs={};
for t = 1:length(myType)
    coll = col(t+2,:);
    typeName = myType{t}; 
    dt=data(t).table;
    labels = dt.TT(:,1:3);
    % vis_areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl"];
    vis_areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl", "RSPd", "RSPagl", "RSPv"];
    idSC = find(contains(labels(:,1),vis_areas));
    threCur= dt.Threshold;
    myCoord = dt.CoordMerge;
    myCoord=myCoord(idSC,:);
%     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
%         if threCur(i,1)==0
%             myCoord(i,:)={[]};
%         end
%     end  
    %plot transformed coordinates
    %figure(f);
    [sa, sl] = size(myCoord);coordXYZ = [];
    theta=deg2rad(-28.022);
    phi=deg2rad(-8.29);
    Coord_type=[];
    for i = 1:sl
        coordXYZ = cell2mat(myCoord(:,i));
        coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:);
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
%         s=scatter(-coordXYZ(:,3), -coordXYZ(:,1), 10, 'filled', CData=coll, MarkerFaceAlpha=0.3);  
%         hold on;
    end
    type_len(1,t)=length(Coord_type);
    
    z = ksdensity([-Coord_type(:,3),-Coord_type(:,1)], [-Coord_type(:,3),-Coord_type(:,1)]);
    x=-Coord_type(:,3);
    y=-Coord_type(:,1);
    newx=x+randn(size(x))/100000;%add some small random noise to not get duplicate x-y values for interpolation
    newy=y+randn(size(y))/100000;
    N = 100;
    [X, Y] = meshgrid(linspace(min(newx)+0, max(newx)-0, N), linspace(min(newy)+0, max(newy)-0, N));
    [Xnan, Ynan] = meshgrid(linspace(min(newx)+0.23, max(newx)-0.23, N), linspace(min(newy)+0.23, max(newy)-0.23, N));
    %Zi=griddata(newx,newy,z,X,Y, 'linear');
    %cutoff_distance = 6;
%     DT = delaunayTriangulation(newx, newy);
%     [vi, d] = nearestNeighbor(DT, newx(:), newy(:));
    %mask = d > cutoff_distance;
    %Zi(mask) = nan;
    Fnan=scatteredInterpolant(newx, newy, z, 'linear', 'none');
    F = scatteredInterpolant(newx, newy, z, 'linear', 'linear');
    ZZnan{1,t}= Fnan(Xnan, Ynan);
    %ZZs{1,t} = Zi;
    ZZs{1,t} = F(X, Y);
    v=0.05:0.06:1; %keep only contour levels with these values
    %v=Z(Z>=0);
    %smooth contour
    n = 7;%higher values smooth more
    kernel = ones(n)/n.^2 ;
    Zs = conv2(ZZs{1,t}, kernel,'same') ;
%     [C,h]=contour(X,Y,Zs,v);
%     h.LineWidth=1.5;
%     h.EdgeColor=coll;
    hold on;
end
%uncomment these strings if lines if you want density normalized to med+lat together
% z = ksdensity([-Coordinates_tot(:,3),-Coordinates_tot(:,1)], [-Coordinates_tot(:,3),-Coordinates_tot(:,1)]);
% z(type_len(1,1):end, :)=-z(type_len(1,1):end, :);
% x=-Coordinates_tot(:,3);
% y=-Coordinates_tot(:,1);
% z=ZZ';
% N = 100;
% [X, Y] = meshgrid(linspace(min(x)-0.2, max(x)+0.2, N), linspace(min(y)-0.2, max(y)+0.2, N));
% F = scatteredInterpolant(x, y, z, 'linear', 'none');
% Z = F(X, Y);
% [C,h] = contour(X, Y, Z);
% h.LineWidth=1.5;


%plot contours
X=-Coordinates_tot(:,3);
Y=-Coordinates_tot(:,1);
newX=X+randn(size(X))/100000;%add some small random noise to not get duplicate x-y values for interpolation
newY=Y+randn(size(Y))/100000;
dt = delaunayTriangulation(newX,newY);
k = convexHull(dt);
ex=newX(k);
ey=newY(k);
p = polyshape(ex,ey); %create a polygon around the scattered data
[cx, cy] = centroid(p); %find the centroid
m = scale(p, 1.0, [cx cy]); %scale it by a factor (1.01 for countour) 
scaled =m.Vertices;
% plot(newX(k), newY(k), 'r'); hold on;
N = 100;
Zs=ZZs{1,1}-ZZs{1,2};
[XX,YY]=meshgrid(linspace(min(newX), max(newX), N), linspace(min(newY), max(newY), N));
X=linspace(min(newX), max(newX), N); Y=linspace(min(newY), max(newY), N);
[in,~] = inpolygon(XX,YY,scaled(:,1),scaled(:,2)); %find indices of XX and YY that are inside the polygon
Zs(~in)=0; %put nan when plotting with contour, 0 when plotting with imagesc
%imagesc(Zs);
n = 6;%higher values smooth more
kernel = ones(n)/n.^2 ;
Zfin = conv2(Zs, kernel,'same');
% Zfin(~in)=nan; %replace Z values outside of polygon with nan
% v=-0.4:0.01:0.4;
% % 
% [C,h]=contour(XX,YY,Zfin,v); %try to fix open contours with rmoutliers
% h.LineWidth=1.5;


% colors = cbrewer('div', 'PiYG', 300); colors=flipud(colors);
%colors=div_cmap;
%colors=CustomColormap;
% colors=div_cmap2;
colors=div_cmap_grp2;
%colors=div_cmap2;
Zmax=max(abs(Zfin),[],'all');
colormap(colors); colorbar; 
% clim([-0.5 +0.5]);
clim([-0.3 +0.3]);
im=imagesc(XX(1,:),YY(:,1),Zfin);
uistack(im, 'bottom');

% ZZ=peaks(XX,YY);
% dim=length(ZZ{1,1}) + length(ZZ{1,2});
% Zmed=repmat(ZZ{1,1}(1,:),dim,1);
% Zlat=repmat(ZZ{1,2}(1,:),dim,1);
% Z=cat(2, Zmed, Zlat);
% 
% contour(XX,YY,Zs,20);
% 
% v=-3:1:5; % contour levels

xlabel ('ML (mm)'); ylabel ('AP (mm)');

ax=gca;
axis(ax, 'equal');
xlim(ax,[-(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3)]); 
ylim(ax,[-(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3)]);
set(ax,'box','on', 'Color', 'none');%put 'color' zerocolor when plotting imagesc
set(f, "Position", [359.4000 170.6000 743.2000 559.2000], 'renderer','Painters');


%% plot single vis area by color
clc;clf;clear data;clear ttd;clear ttm;
celltypes = {'ntsr'};
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
colors=cbrewer('qual','Set1', 13); colors=colors(randperm(size(colors, 1)), :);
% colors{1,1}=cbrewer('seq','RdPu',12);
% colors{1,2}=cbrewer('seq','BuGn',12);
Coordinates_tot=[];
f=figure(1);
type_len=[];
ZZs={};
vis_areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl", "RSPd", "RSPagl", "RSPv"];

for a=1:length(vis_areas)
    curr_area=vis_areas(a);
    coord_area=[];
    for t = 1:length(myType)
        %coll = col(t+2,:);
        typeName = myType{t};
        dt=data(t).table;
        labels = dt.TT(:,1:3);
        idSC = find(ismember(labels(:,1),curr_area));
        %threCur= dt.Threshold;
        myCoord = dt.CoordMerge;
        myCoord=myCoord(idSC,:);
        %     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
        %         if threCur(i,1)==0
        %             myCoord(i,:)={[]};
        %         end
        %     end
        %plot transformed coordinates
        %figure(f);
        [sa, sl] = size(myCoord);coordXYZ = [];
        theta=deg2rad(-28.022);
        phi=deg2rad(-8.29);
        Coord_type=[];
        for i = 1:sl
            coordXYZ = cell2mat(myCoord(:,i));
            coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:);
            X=coordXYZ(:,1);
            Y=coordXYZ(:,3);
            Z=coordXYZ(:,2);
            Xnew=X*cos(phi)-Y*sin(theta)*sin(phi)-Z*cos(theta)*sin(phi);
            Ynew=Y*cos(theta)-Z*sin(theta);
            Znew=X*sin(phi)+Y*sin(theta)*cos(phi)+Z*cos(theta)*cos(phi);
            coordXYZ=[Xnew Znew Ynew];
            Coordinates_tot=[Coordinates_tot; coordXYZ];
            Coord_type=[Coord_type; coordXYZ];
            s=scatter(-coordXYZ(:,3), -coordXYZ(:,1), 10, 'filled', CData=colors(a,:), MarkerFaceAlpha=0.5);
            hold on;
            
        end
        coord_area=[coord_area; Coord_type];
       
        type_len(1,t)=length(Coord_type);

        z = ksdensity([-Coord_type(:,3),-Coord_type(:,1)], [-Coord_type(:,3),-Coord_type(:,1)]);
        x=-Coord_type(:,3);
        y=-Coord_type(:,1);
        newx=x+randn(size(x))/100000;%add some small random noise to not get duplicate x-y values for interpolation
        newy=y+randn(size(y))/100000;
        N = 100;
        [X, Y] = meshgrid(linspace(min(newx)+0, max(newx)-0, N), linspace(min(newy)+0, max(newy)-0, N));
        [Xnan, Ynan] = meshgrid(linspace(min(newx)+0.23, max(newx)-0.23, N), linspace(min(newy)+0.23, max(newy)-0.23, N));
        %Zi=griddata(newx,newy,z,X,Y, 'linear');
        %cutoff_distance = 6;
        %     DT = delaunayTriangulation(newx, newy);
        %     [vi, d] = nearestNeighbor(DT, newx(:), newy(:));
        %mask = d > cutoff_distance;
        %Zi(mask) = nan;
        Fnan=scatteredInterpolant(newx, newy, z, 'linear', 'none');
        F = scatteredInterpolant(newx, newy, z, 'linear', 'linear');
        ZZnan{1,t}= Fnan(Xnan, Ynan);
        %ZZs{1,t} = Zi;
        ZZs{1,t} = F(X, Y);
        %     v=0.05:0.06:1; %keep only contour levels with these values
        %v=Z(Z>=0);
        %smooth contour
        %     n = 6;%higher values smooth more
        %     kernel = ones(n)/n.^2 ;
        %     Zs = conv2(Z, kernel,'same') ;
        %     %[C,h]=contour(X,Y,Zs,v);
        %     h.LineWidth=1.5;
        %     h.EdgeColor=coll;
        %     hold on;
    end
    k=boundary(-coord_area(:,3), -coord_area(:,1), 0.2);
    plot(-coord_area(k,3), -coord_area(k,1),'color', colors(a,:), 'LineWidth', 1.5);
    hold on
end
xlabel ('ML (mm)'); ylabel ('AP (mm)');

ax=gca;
axis(ax, 'equal');
xlim(ax,[-(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3)]); 
ylim(ax,[-(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3)]);
set(ax,'box','on');%put 'color' zerocolor when plotting imagesc
set(f, "Position", [359.4000 170.6000 743.2000 559.2000]);
%uncomment these strings if lines if you want density normalized to med+lat together
% z = ksdensity([-Coordinates_tot(:,3),-Coordinates_tot(:,1)], [-Coordinates_tot(:,3),-Coordinates_tot(:,1)]);
% z(type_len(1,1):end, :)=-z(type_len(1,1):end, :);
% x=-Coordinates_tot(:,3);
% y=-Coordinates_tot(:,1);
% z=ZZ';
% N = 100;
% [X, Y] = meshgrid(linspace(min(x)-0.2, max(x)+0.2, N), linspace(min(y)-0.2, max(y)+0.2, N));
% F = scatteredInterpolant(x, y, z, 'linear', 'none');
% Z = F(X, Y);
% [C,h] = contour(X, Y, Z);
% h.LineWidth=1.5;

%% plot scatterplothistogram %version below is better

clc;clf;clear data;clear ttd;clear ttm;
load(fullfile(spath,'OD.mat'));
load(fullfile(spath,'Coord.mat'));
myType = flip(unique(OD.type),2);
Coordinates_tot=[];
coord_tot={};
f=figure(1);
set(gca,'box', 'off', 'Color', 'none');
type_len=[];
for t = 1:length(myType)
    coll = col(t+2,:);
    typeName = myType{t}; 
    dt = load(fullfile(spath,[typeName,'.mat']));
    dt.TT_title;
    labels = dt.TT(:,1:3);
    vis_areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl"];
    idSC = find(contains(labels(:,1),vis_areas));
    threCur= dt.Threshold;
    myCoord = dt.CoordMerge;
    myCoord=myCoord(idSC,:);
%     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
%         if threCur(i,1)==0
%             myCoord(i,:)={[]};
%         end
%     end  
    %plot transformed coordinates
    %figure(f);
    [sa, sl] = size(myCoord);coordXYZ = [];
    theta=deg2rad(-28.022);
    phi=deg2rad(-8.29);
    Coord_type=[];
    for i = 1:sl
        coordXYZ = cell2mat(myCoord(:,i));
        coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:);
        X=coordXYZ(:,1); 
        Y=coordXYZ(:,3); 
        Z=coordXYZ(:,2);
        Xnew=X*cos(phi)-Y*sin(theta)*sin(phi)-Z*cos(theta)*sin(phi);
        Ynew=Y*cos(theta)-Z*sin(theta);
        Znew=X*sin(phi)+Y*sin(theta)*cos(phi)+Z*cos(theta)*cos(phi);
        coordXYZ=[Xnew Znew Ynew];
        Coordinates_tot=[Coordinates_tot; coordXYZ];
        Coord_type=[Coord_type; coordXYZ];
%         s=scatter(-coordXYZ(:,3), -coordXYZ(:,1), 10, 'filled', CData=coll, MarkerFaceAlpha=0.3);  
%         hold on;
    end
    coord_tot{1,t}=Coord_type;
end
%create table for medial points
arraymed=[-coord_tot{1,1}(:,3), -coord_tot{1,1}(:,1)];
tblmed=array2table(arraymed,"VariableNames", ["xval" "yval"]);
tblmed.Position=repelem(string(myType{1,1}), length(coord_tot{1,1}))';

%create table for lateral points
arraylat=[-coord_tot{1,2}(:,3), -coord_tot{1,2}(:,1)];
tbllat=array2table(arraylat,"VariableNames", ["xval" "yval"]);
tbllat.Position=repelem(string(myType{1,2}), length(coord_tot{1,2}))';

%concatenate tables
tbl=[tblmed; tbllat];
s=scatterhistogram(f,tbl,"xval","yval", 'GroupVariable','Position', 'LegendVisible', 'off', 'HistogramDisplayStyle','smooth',...
    'LineStyle', '-', 'LineWidth', 1.5, 'MarkerSize', 10, 'MarkerAlpha', 0.3, 'ScatterPlotProportion', 0.75); 
s.Color={col(3,:), col(4,:)}; s.XLimits=[-(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3)];
s.YLimits=[-(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3)];

xlabel ('ML (mm)'); ylabel ('AP (mm)');

% ax=gca;
% axis(ax, 'equal');
% xlim(ax,[-(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3)]); 
% ylim(ax,[-(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3)]);
%set(gca,'box','on', 'Color', 'none');%put 'color' zerocolor when plotting imagesc
set(f, "Position", [359.4000 170.6000 743.2000 559.2000]);

%% plot ksdensity with scatter
clc;clf;clear data;clear ttd;clear ttm;
list_id = [];myexp_labels = {};nn = 0;n_Data=[];OD={}; Coord={};tt=0;
clr=2; %0 for ntsr, 2 for grp
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
Coordinates_tot=[];
coord_tot={};


type_len=[];
ZZs={};
for t = 1:length(myType)
    coll = col(t,:);
    typeName = myType{t}; 
    dt=data(t).table;
    labels = dt.TT(:,1:3);
    % vis_areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl"];
    vis_areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl", "RSPd", "RSPagl", "RSPv"];
    idSC = find(contains(labels(:,1),vis_areas));
    threCur= dt.Threshold;
    myCoord = dt.CoordMerge;
    myCoord=myCoord(idSC,:);
%     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
%         if threCur(i,1)==0
%             myCoord(i,:)={[]};
%         end
%     end  
    %plot transformed coordinates
    %figure(f);
    [sa, sl] = size(myCoord);coordXYZ = [];
    theta=deg2rad(-28.022);
    phi=deg2rad(-8.29);
    Coord_type=[];
    for i = 1:sl
        coordXYZ = cell2mat(myCoord(:,i));
        coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:);
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
%         s=scatter(-coordXYZ(:,3), -coordXYZ(:,1), 10, 'filled', CData=coll, MarkerFaceAlpha=0.3);  
%         hold on;
    end
    type_len(1,t)=length(Coord_type);
    
    z = ksdensity([-Coord_type(:,3),-Coord_type(:,1)], [-Coord_type(:,3),-Coord_type(:,1)]);
    x=-Coord_type(:,3);
    y=-Coord_type(:,1);
    newx=x+randn(size(x))/100000;%add some small random noise to not get duplicate x-y values for interpolation
    newy=y+randn(size(y))/100000;
    N = 100;
    [X, Y] = meshgrid(linspace(min(newx)+0, max(newx)-0, N), linspace(min(newy)+0, max(newy)-0, N));
    [Xnan, Ynan] = meshgrid(linspace(min(newx)+0.23, max(newx)-0.23, N), linspace(min(newy)+0.23, max(newy)-0.23, N));
    %Zi=griddata(newx,newy,z,X,Y, 'linear');
    %cutoff_distance = 6;
%     DT = delaunayTriangulation(newx, newy);
%     [vi, d] = nearestNeighbor(DT, newx(:), newy(:));
    %mask = d > cutoff_distance;
    %Zi(mask) = nan;
    Fnan=scatteredInterpolant(newx, newy, z, 'linear', 'none');
    F = scatteredInterpolant(newx, newy, z, 'linear', 'linear');
    ZZnan{1,t}= Fnan(Xnan, Ynan);
    %ZZs{1,t} = Zi;
    ZZs{1,t} = F(X, Y);
    v=0.05:0.06:1; %keep only contour levels with these values
    %v=Z(Z>=0);
    %smooth contour
    n = 7;%higher values smooth more
    kernel = ones(n)/n.^2 ;
    Zs = conv2(ZZs{1,t}, kernel,'same') ;
%     [C,h]=contour(X,Y,Zs,v);
%     h.LineWidth=1.5;
%     h.EdgeColor=coll;
    coord_tot{1,t}=Coord_type;
    hold on;
end


%plot
f=figure(1);
fact=2;
til=tiledlayout(3*fact,3*fact, 'TileSpacing','none', 'Padding','tight');
%plot horizontal density curve
nexttile([1 2*fact+1]);
% h=histfit(-coord_tot{1,1}(:,3),100,'kernel');delete(h(1));h(2).Visible='off'; %medial
% med_x=h(2).XData; med_y=h(2).YData;
% h=histfit(-coord_tot{1,2}(:,3),100,'kernel');delete(h(1));h(2).Visible='off'; %lateral
% lat_x=h(2).XData; lat_y=h(2).YData;
% plot(med_x, med_y/max(max(med_y), max(lat_y)),'LineWidth',2, 'Color',col(3,:)); hold on;
% plot(lat_x, lat_y/max(max(med_y), max(lat_y)),'LineWidth',2, 'Color',col(4,:));
t=1; %1 for grp, 3 for ntsr
[yi,~] = ksdensity(-coord_tot{1,1}(:,3),'Support', [min(-coord_tot{1,1}(:,3))-0.1 ...
    max(-coord_tot{1,1}(:,3))+0.1], 'NumPoints',300);
xi=linspace(min(-coord_tot{1,1}(:,3)), max(-coord_tot{1,1}(:,3)),300);
plot(xi,yi,'LineWidth',2, 'Color',col(t+clr,:)); hold on;
[yi,~] = ksdensity(-coord_tot{1,2}(:,3),'Support', [min(-coord_tot{1,2}(:,3))-0.1 ...
    max(-coord_tot{1,2}(:,3))+0.1], 'NumPoints',300);
xi=linspace(min(-coord_tot{1,2}(:,3)), max(-coord_tot{1,2}(:,3)),300);
plot(xi,yi,'LineWidth',2, 'Color',col(t+1+clr,:));

xlim(gca,[-(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3)]); 
set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');

ax1=nexttile(4*fact-1,[2*fact+1 2*fact+1]);

%plot scatter points
s_med=scatter(-coord_tot{1,1}(:,3), -coord_tot{1,1}(:,1), 10, 'filled', CData=col(t+clr,:), MarkerFaceAlpha=0.35); hold on;
s_lat=scatter(-coord_tot{1,2}(:,3), -coord_tot{1,2}(:,1), 10, 'filled', CData=col(t+1+clr,:), MarkerFaceAlpha=0.35);hold on;


axis(ax1, 'equal');
xlim(ax1,[-(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3)]); 
ylim(ax1,[-(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3)]);
set(ax1,'box','on', 'Color', 'none');%put 'color' zerocolor when plotting imagesc
ax1.FontSize=12;
xlabel ('M-L', 'FontSize', 14); ylabel ('A-P','FontSize', 14);


%plot vertical density curve
nexttile(6*fact, [2*fact+1 1])
% h=histfit(-coord_tot{1,1}(:,1),100,'kernel');delete(h(1));h(2).Visible='off';
% med_y = h(2).XData; med_x = h(2).YData;
% h=histfit(-coord_tot{1,2}(:,1),100,'kernel');delete(h(1));h(2).Visible='off';
% lat_y = h(2).XData; lat_x = h(2).YData;
% plot(med_x/max(med_x), med_y, 'LineWidth',2, 'Color',col(3,:)); hold on;
% plot(lat_x/max(lat_x), lat_y, 'LineWidth',2, 'Color',col(4,:));

% [yi,~] = ksdensity(-coord_tot{1,1}(:,3),'Support', [min(-coord_tot{1,1}(:,3))-0.1 ...
%     max(-coord_tot{1,1}(:,3))+0.1], 'NumPoints',300);
% xi=linspace(min(-coord_tot{1,1}(:,3)), max(-coord_tot{1,1}(:,3)),300);
% plot(xi,yi,'LineWidth',2, 'Color',col(3,:)); hold on;
% [yi,~] = ksdensity(-coord_tot{1,2}(:,3),'Support', [min(-coord_tot{1,2}(:,3))-0.1 ...
%     max(-coord_tot{1,2}(:,3))+0.1], 'NumPoints',300);
% xi=linspace(min(-coord_tot{1,2}(:,3)), max(-coord_tot{1,2}(:,3)),300);
% plot(xi,yi,'LineWidth',2, 'Color',col(4,:));

[yi,~] = ksdensity(-coord_tot{1,1}(:,1),'Support', [min(-coord_tot{1,1}(:,1))-0.1 ...
    max(-coord_tot{1,1}(:,1))+0.1], 'NumPoints',300, 'BandWidth', 0.2);
xi=linspace(min(-coord_tot{1,1}(:,1)), max(-coord_tot{1,1}(:,1)),300);
plot(yi,xi,'LineWidth',2, 'Color',col(t+clr,:)); hold on;
[yi,~] = ksdensity(-coord_tot{1,2}(:,1),'Support', [min(-coord_tot{1,2}(:,1))-0.1 ...
    max(-coord_tot{1,2}(:,1))+0.1], 'NumPoints',300);
xi=linspace(min(-coord_tot{1,2}(:,1)), max(-coord_tot{1,2}(:,1)),300);
plot(yi,xi,'LineWidth',2, 'Color',col(t+1+clr,:));

ylim(gca,[-(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3)]);
set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');

%set(gca,'box','on', 'Color', 'none');%put 'color' zerocolor when plotting imagesc
%make figure square so that scatter axes are equal
x_range=(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3);
y_range=(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3);
xy_ratio=y_range./x_range;
fig_width=582;
set(f, "Position", [572.2000 61 fig_width fig_width*xy_ratio], 'renderer','Painters');

%% plot ksdensity with scatter - ONE COLOR


clc;clf;clear data;clear ttd;clear ttm;
load(fullfile(spath,'OD.mat'));
load(fullfile(spath,'Coord.mat'));
myType = flip(unique(OD.type),2);
Coordinates_tot=[];
coord_tot={};


type_len=[];
for t = 1:length(myType)
    coll = col(t+2,:);
    typeName = myType{t}; 
    dt = load(fullfile(spath,[typeName,'.mat']));
    dt.TT_title;
    labels = dt.TT(:,1:3);
    vis_areas=["VISal", "VISam", "VISl", "VISli", "VISpl", "VISpor", "VISp", "VISpm", "VISa", "VISrl"];
    idSC = find(ismember(labels(:,1),vis_areas));
    threCur= dt.Threshold;
    myCoord = dt.CoordMerge;
    myCoord=myCoord(idSC,:);
%     for i=1:length(myCoord) %this for loop applies the thresholding computed in LoadData to the data before plotting the distributions
%         if threCur(i,1)==0
%             myCoord(i,:)={[]};
%         end
%     end  
    %plot transformed coordinates
    %figure(f);
    [sa, sl] = size(myCoord);coordXYZ = [];
    theta=deg2rad(-28.022);
    phi=deg2rad(-8.29);
    Coord_type=[];
    for i = 1:sl
        coordXYZ = cell2mat(myCoord(:,i));
        coordXYZ=coordXYZ(coordXYZ(:,3)>=0,:);
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
%         s=scatter(-coordXYZ(:,3), -coordXYZ(:,1), 10, 'filled', CData=coll, MarkerFaceAlpha=0.3);  
%         hold on;
    end
    coord_tot{1,t}=Coord_type;
end

%plot
f=figure(1);
fact=2;
til=tiledlayout(3*fact,3*fact, 'TileSpacing','none', 'Padding','tight');
%plot horizontal density curve
nexttile([1 2*fact+1]);


[yi,~] = ksdensity(-coord_tot{1,1}(:,3),'Support', [min(-coord_tot{1,1}(:,3))-0.1 ...
    max(-coord_tot{1,1}(:,3))+0.1], 'NumPoints',300);
xi=linspace(min(-coord_tot{1,1}(:,3)), max(-coord_tot{1,1}(:,3)),300);
plot(xi,yi,'LineWidth',2, 'Color',col(3,:)); hold on;
% [yi,~] = ksdensity(-coord_tot{1,2}(:,3),'Support', [min(-coord_tot{1,2}(:,3))-0.1 ...
%     max(-coord_tot{1,2}(:,3))+0.1], 'NumPoints',300);
% xi=linspace(min(-coord_tot{1,2}(:,3)), max(-coord_tot{1,2}(:,3)),300);
% plot(xi,yi,'LineWidth',2, 'Color',col(4,:));

xlim(gca,[-(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3)]); 
set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');

ax1=nexttile(4*fact-1,[2*fact+1 2*fact+1]);
%plot scatter points
s_med=scatter(-coord_tot{1,1}(:,3), -coord_tot{1,1}(:,1), 10, 'filled', CData=col(3,:), MarkerFaceAlpha=0.35); hold on;
% s_lat=scatter(-coord_tot{1,2}(:,3), -coord_tot{1,2}(:,1), 10, 'filled', CData=col(4,:), MarkerFaceAlpha=0.35);hold on;

set(gca,'box', 'on', 'Color', 'none');
xlim(gca,[-(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3)]); 
ylim(gca,[-(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3)]);
ax=gca; ax.FontSize=10;
xlabel ('M-L', 'FontSize', 12); ylabel ('A-P','FontSize', 12);


%plot vertical density curve
nexttile(6*fact, [2*fact+1 1])
% h=histfit(-coord_tot{1,1}(:,1),100,'kernel');delete(h(1));h(2).Visible='off';
% med_y = h(2).XData; med_x = h(2).YData;
% h=histfit(-coord_tot{1,2}(:,1),100,'kernel');delete(h(1));h(2).Visible='off';
% lat_y = h(2).XData; lat_x = h(2).YData;
% plot(med_x/max(med_x), med_y, 'LineWidth',2, 'Color',col(3,:)); hold on;
% plot(lat_x/max(lat_x), lat_y, 'LineWidth',2, 'Color',col(4,:));

% [yi,~] = ksdensity(-coord_tot{1,1}(:,3),'Support', [min(-coord_tot{1,1}(:,3))-0.1 ...
%     max(-coord_tot{1,1}(:,3))+0.1], 'NumPoints',300);
% xi=linspace(min(-coord_tot{1,1}(:,3)), max(-coord_tot{1,1}(:,3)),300);
% plot(xi,yi,'LineWidth',2, 'Color',col(3,:)); hold on;
% [yi,~] = ksdensity(-coord_tot{1,2}(:,3),'Support', [min(-coord_tot{1,2}(:,3))-0.1 ...
%     max(-coord_tot{1,2}(:,3))+0.1], 'NumPoints',300);
% xi=linspace(min(-coord_tot{1,2}(:,3)), max(-coord_tot{1,2}(:,3)),300);
% plot(xi,yi,'LineWidth',2, 'Color',col(4,:));

[yi,~] = ksdensity(-coord_tot{1,1}(:,1),'Support', [min(-coord_tot{1,1}(:,1))-0.1 ...
    max(-coord_tot{1,1}(:,1))+0.1], 'NumPoints',300, 'BandWidth', 0.2);
xi=linspace(min(-coord_tot{1,1}(:,1)), max(-coord_tot{1,1}(:,1)),300);
plot(yi,xi,'LineWidth',2, 'Color',col(3,:)); hold on;
% [yi,~] = ksdensity(-coord_tot{1,2}(:,1),'Support', [min(-coord_tot{1,2}(:,1))-0.1 ...
%     max(-coord_tot{1,2}(:,1))+0.1], 'NumPoints',300);
% xi=linspace(min(-coord_tot{1,2}(:,1)), max(-coord_tot{1,2}(:,1)),300);
% plot(yi,xi,'LineWidth',2, 'Color',col(4,:));

ylim(gca,[-(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3)]);
set(gca,'box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');

%set(gca,'box','on', 'Color', 'none');%put 'color' zerocolor when plotting imagesc
%make figure square so that scatter axes are equal
x_range=(max(Coordinates_tot(:,3))+0.3) -(min(Coordinates_tot(:,3))-0.3);
y_range=(max(Coordinates_tot(:,1))+0.3) -(min(Coordinates_tot(:,1))-0.3);
xy_ratio=y_range./x_range;
fig_width=582;
set(f, "Position", [572.2000 61 fig_width fig_width*xy_ratio], 'renderer','Painters');
%% calculations
% P21 =[0.18 0.18 0.2 0.25 0.25 0.26 0.35 0.35 0.41 0.45 0.47 0.5];
% P90 =[0.1 0.18 0.36 0.45 0.47 0.5 0.55 0.58 0.65 0.77];
 P21 =[0.8 0.9 1 1.15 1.4 1.55 1.8 1.9 2.1 2.4 3];
 P90 = [0 0.4 0.4 0.5 0.55 0.9 1.3 1.4 1.45];
mean_21=mean(P21, 2);
mean_90=mean(P90 ,2);
%mean_90=21.5;
std_21=std(P21);
std_90=std(P90);
pooled_std=sqrt(((std_21.^2)+(std_90.^2))./2);
d=(mean_90-mean_21)./pooled_std; %effect size
d=(0.45-0.32)/0.166

%%
test=readNPY("brainGridData.npy");

%%
lulz=[];
for i=1:length(BBB(:,4))
    for j=1:length(BBB{i,4})
        curr=BBB{i,4}(j);
        if curr==870
            lulz=i;
        end
    end
end
