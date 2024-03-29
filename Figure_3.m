%% Z20 equivalent depth and TP=0.8 depth
addpath('Functions/')
addpath('Data/')
load('Z20.mat','good_data','Z20')
good_list=~any(isnan(good_data(:,[3,5,6,7])),2);
Z20=Z20(good_list);
good_data=good_data(good_list,:);
clear good_list

a=[-8.4,-10.667];
%Cores VM28-227 and VM28-229 excluded due to questionable d18O
%stratigraphy, see Karim Lakhani's Thesis, chapter 3.2.1
for i=1:length(a)
    b(i)=find(good_data(:,1)==a(i));
end
good_data(b,:)=[];Z20(b)=[];
%Clean up

%% Z20 map loading part 1
load('d18O_paper.mat')
Lat=squeeze(X3(:,:,1));
Lon=squeeze(Y3(:,:,1));

[newLon,newLat,newd18O]=cycle_data_3D(Lon,Lat,d18O1,30);
latvec=squeeze(newLat(1,:,1));
lonvec=squeeze(newLon(:,1,1));
depthvec=squeeze(Z3(1,1,:));
clear X3 Y3 Z3 Lat Lon d18O1

[~,I5]=min(abs((lonvec-180)));
res=5;
for j=1:30
    for i=234:394
        newd18O(I5-res:I5+res,i,j)=interp1([lonvec(I5-res),lonvec(I5+res)],[squeeze(newd18O(I5-res,i,j)),squeeze(newd18O(I5+res,i,j))],lonvec(I5-res:I5+res));
    end
end
%% Z20 map loading part 2
d18O_Z20=-0.6568;
[X3,Y3,Z3]=meshgrid(latvec,lonvec,depthvec);

for i=400:1040
    for j=234:394
        Ilon=i;
        Ilat=j;

        range=2;
        data=newd18O(mindata(Ilon-range):maxdata(Ilon+range,1440),mindata(Ilat-range):maxdata(Ilat+range,674),1:30);
        data2=data;
        data(data==-900)=900;

        data3=newd18O(mindata(Ilon-range):maxdata(Ilon+range,1440),mindata(Ilat-range):maxdata(Ilat+range,674),:);
        X3sub=X3(mindata(Ilon-range):maxdata(Ilon+range,1440),mindata(Ilat-range):maxdata(Ilat+range,674),:);
        Y3sub=Y3(mindata(Ilon-range):maxdata(Ilon+range,1440),mindata(Ilat-range):maxdata(Ilat+range,674),:);
        Z3sub=Z3(mindata(Ilon-range):maxdata(Ilon+range,1440),mindata(Ilat-range):maxdata(Ilat+range,674),:);        

        if data ==data2
            %expectedd18O=interp3(X3sub,Y3sub,Z3sub,data3,final_lat_lon(i,1),final_lat_lon(i,2),species_depth(j),'linear');
            profile=squeeze(interp3(X3sub,Y3sub,Z3sub,data3,latvec(j),lonvec(i),depthvec,'linear'));
            profile=profile+fliplr(1:length(profile))'*0.00000001;
            expectedDepth=interp1(squeeze(profile(3:30)),squeeze(depthvec(3:30)),d18O_Z20);
            Z20_field(i,j)=expectedDepth;
        else
            Z20_field(i,j)=nan;
        end
    end
end
uncertainty=nanstd(Z20_field(:,294:334),'',2);

LGM=0;
bott=40;
[~,I1]=min(abs(latvec+5));
[~,I2]=min(abs(latvec-5));
XZ_d18O=squeeze(nanmean(newd18O(:,I1:I2,1:bott),2));
[Depth,Long]=meshgrid(depthvec(1:bott),lonvec);
for i=1:length(depthvec(1:bott))
    XZ_d18O(598:603,i)=interp1([Long(598,i),Long(603,i)],[XZ_d18O(598,i),XZ_d18O(603,i)],Long(598:603,i));
end

Pacific=[150,281];
[~,I3]=min(abs(lonvec-Pacific(1)));
[~,I4]=min(abs(lonvec-Pacific(2)));
for i=I3:I4
    profile=XZ_d18O(i,:);
    dep=depthvec(1:bott);
    if ~isnan(profile(end))
        Therm_depth(i)=interp1(profile,dep,-0.66);
    end
end
%% Generate average profile, calculate TP=0.8 from model, and get average d18Osw and Z20 equivalent in d18Oc
depths=[0,nan,210,114,96,610];
benthic_d18O=1.75;
model_TP80=ones([length(good_data),1])*nan;
log_func=@(depth,beta) -1*beta(1).^(-1.*(depth+beta(2)))+beta(3);
for i=1:length(good_data)
    [beta,MLD]=run_thermocline_model3(depths([1,3:end]),[good_data(i,[3,5,6,7]),benthic_d18O]);
    depth_plot=linspace(0,610,1000);
    profile=depth_plot*nan;
    profile(depth_plot<MLD)=good_data(i,3);
    profile(depth_plot>=MLD)=log_func(depth_plot(depth_plot>=MLD),beta);
    d80=profile(end)-0.8*abs(profile(1)-profile(end));
    model_TP80(i)=interp1(profile(depth_plot>=MLD),depth_plot(depth_plot>=MLD),d80);
    model_Z20_equivalent(i)=interp1(profile(depth_plot>=MLD),depth_plot(depth_plot>=MLD),d18O_Z20);
end

%%
model_Z20_equivalent=model_Z20_equivalent(:);
cax2=[0,225];
figure('Position',[2708.2,136.2,842.4,596])
subplot(2,1,1)
hold on
world_coastp
c=contourf(Y3(400:1040,234:394,1),X3(400:1040,234:394,1),Z20_field(400:1040,234:394),[0,25,50,75,100,125,150,175,200,225],'EdgeColor','none');
scatter(good_data(:,2),good_data(:,1),55,model_Z20_equivalent,'filled','MarkerEdgeColor','k','LineWidth',1);
caxis(cax2)
colormap(colorblind_jet)
colorbar('Ticks',[0,25,50,75,100,125,150,175,200,225],'TickLabels',{'0','25','50','75','100','125','150','175','200','225'})

axis([130,290,-20,20])
title('Holocene Z20 equivalent from regression')
xticks([160, 200, 240, 280]);
xticklabels({'160\circE','160\circW','120\circW','80\circW'});
yticks([-15,0,15]);
yticklabels({'15\circS','0\circ','15\circN'});

%%
benthic_d18O=1.75;
depth_plot=linspace(0,610,1000);
MC_num=100;
MC_profiles=ones([length(good_data),MC_num,length(depth_plot)])*nan;
Z20_MC=ones([length(good_data),MC_num])*nan;
for k=1:length(good_data)
    data=good_data(k,:);
    data(:,4)=nan;
    for i=1:MC_num
        disp(strcat(num2str(k),"   ",num2str(i)))
        MC_depths=[0];
        for j=2:5
            MC_depths(j)=depth_realization(xi_tot(j,:),ksdens(j,:),minmax(j,1),minmax(j,2));
        end
        MC_depths(6)=610;
        [MC_beta,MC_MLD]=run_thermocline_model3(MC_depths,[data(3:7),benthic_d18O]);
        MC_profiles(k,i,depth_plot<MC_MLD)=data(3);
        MC_profiles(k,i,depth_plot>=MC_MLD)=log_func(depth_plot(depth_plot>=MC_MLD),MC_beta);
        [~,I]=min(abs(squeeze(MC_profiles(k,i,:))--0.66));
        Z20_MC(k,i)=depth_plot(I);
    end
end
means_Hol=mean(Z20_MC,2);
stds_Hol=std(Z20_MC,'',2);
%%
uncertainty(uncertainty<5)=5;
lat_lon=good_data(:,[1,2]);
Equatorial=find(abs(lat_lon(:,1))<5);
subplot(2,1,2)
hold on
plot(lonvec(I3:I4),Therm_depth(I3:I4),'-k')
Therm_depth=Therm_depth(:);
fill([lonvec(I3:I4);flipud(lonvec(I3:I4))],[Therm_depth(I3:I4)-uncertainty(I3:I4);flipud(Therm_depth(I3:I4)+uncertainty(I3:I4))],[0.5,0.5,0.5],'FaceAlpha',0.5)
errorbar(good_data(Equatorial,2),means_Hol(Equatorial),stds_Hol(Equatorial),'.k','LineWidth',2,'MarkerSize',15)

axis([xlim,0,300])
xticks([160, 200, 240, 280]);
xticklabels({'160\circE','160\circW','120\circW','80\circW'});
set(gca,'YDir','reverse')
title('Holocene EW transect')
ylabel('Depth (m)')
legend('Z20 depth','Uncertainty','Regression results','Location','SouthEast')