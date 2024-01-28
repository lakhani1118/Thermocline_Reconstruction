%% Sections a-d
addpath('Functions');addpath('Data')
load('Data/Z20.mat')
it_list={'{\it G. ruber} (white)','{\it G. tumida}','{\it N. dutertrei}','{\it P. obliquiloculata}'};
%Data is in 28x7 structure with columns Latitude, Longitude, G. ruber, T. sacculifer, G.
%tumida, N. dutertrei, and P. obliq d18O

%Remove rows that do not have species data
goods=~any(isnan(good_data(:,[3,5,6,7])),2);
good_data=good_data(goods,[1,2,3,5,6,7]);
clear goods
%Data now doesn't include the T. sacculifer column because method does not
%rely on that species d18O

%Remove rows that has ________
latitudes=[-8.4,-10.667];
for i=1:2
    row=good_data(:,1)==latitudes(i);
    good_data(row,:)=[];
end
clear latitudes row i

% Because all data is on the same scale, create scale range here
d18Os=good_data(:,[3:end]);
cax=[min(d18Os(:)),max(d18Os(:))];

figure('Position',[0,-114.2,955.2,840])
for i=1:4 %Loop to plot all species data into separate subplots
    subplot(3,2,i)
    hold on
    world_coastp
    scatter(good_data(:,2),good_data(:,1),55,d18Os(:,i),'filled','MarkerEdgeColor','k')
    axis([130,290,-20,20])
    if i==1
        colorbar
    end
    colormap(colorblind_jet)
    caxis(cax)
    title(it_list{i},'FontSize',16)
    xticks([160, 200, 240, 280]);
    xticklabels({'160\circE','160\circW','120\circW','80\circW'});
    yticks([-15,0,15]);
    yticklabels({'15\circS','0\circ','15\circN'});
end
%% Section e data loading
load('d18O_paper.mat') 
%d18O field based on WOA Temp climatology and Legrande and Schmidt 2006 d18Osw climatology
%Used in Lakhani et al. 2022

Lat=squeeze(X3(:,:,1));
Lon=squeeze(Y3(:,:,1));
[cycleLon,cycleLat,d18O]=cycle_data_3D(Lon,Lat,d18O1,30);
%function to change longitude field from (-180 to 180) to be (30 to 390)
%This keeps the Pacific clinmatology together

depthvec=squeeze(Z3(1,1,:));
latvec=cycleLat(1,:);
lonvec=cycleLon(:,1);
[X3,Y3,Z3]=meshgrid(latvec,lonvec,depthvec);
%Reset the X, Y, and Z 3-dimensional matrices
clear Lat Lon d18O1

[~,Depth_600]=min(abs(depthvec-600));
plot_domain_lon=[130,300];
plot_domain_lat=[-20,20];
for i=1:length(plot_domain_lon)
    [~,Iplot_lon(i)]=min(abs(lonvec-plot_domain_lon(i)));
    [~,Iplot_lat(i)]=min(abs(latvec-plot_domain_lat(i)));
end
%Find index in 2-d field for the plotting domain

[~,Ilat]=min(abs(latvec--1.25));
[~,Ilon]=min(abs(lonvec-(-89.68+360)));
%Core location in this grid
deviance_600=d18O(:,:,Depth_600)-d18O(Ilon+2,Ilat+2,Depth_600);
%Adjust core location to a location that has data in this field
%% Section e plotting
subplot(3,1,3)
hold on
[~,I1]=min(abs(lonvec-130));
[~,I2]=min(abs(lonvec-290));
[~,I3]=min(abs(latvec--20));
[~,I4]=min(abs(latvec-20));
[~,I5]=min(abs((lonvec-180)));
res=3;
for i=I3:I4 %Making the field continuous over the International Date Line (Lon=-180 / 180)
    deviance_600(I5-res:I5+res,i)=interp1([lonvec(I5-res),lonvec(I5+res)],[deviance_600(I5-res,i),deviance_600(I5+res,i)],lonvec(I5-res:I5+res));
end
contourf(cycleLon(I1:I2,I3:I4),cycleLat(I1:I2,I3:I4),deviance_600(I1:I2,I3:I4),9,'edgecolor','none')

scatter(good_data(:,2),good_data(:,1),55,'MarkerEdgeColor','k','LineWidth',1.5)
world_coastp
axis([130,290,-20,20])
colorbar

caxis([-0.35,0.35])
scatter(-89.68+360,-1.25,250,'pentagram','b','filled','MarkerEdgeColor','k','LineWidth',1.5)
xticks([160, 200, 240, 280]);
xticklabels({'160\circE','160\circW','120\circW','80\circW'});
yticks([-15,0,15]);
yticklabels({'15\circS','0\circ','15\circN'});