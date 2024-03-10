addpath('Functions')
addpath('Data')
load('Z20.mat')
good_data(:,4)=good_data(:,4)*nan; %Remove T. sacculifer from code
ODP849=find(good_data(:,1)==0.1830);
ODP849=good_data(ODP849,:);
%% Get profile for ODP849
load('d18O_paper.mat')
Lat=squeeze(X3(:,:,1));
Lon=squeeze(Y3(:,:,1));
depthvec=squeeze(Z3(1,1,:));
[cycleLon,cycleLat,cycled18O]=cycle_data_3D(Lon,Lat,d18O1,30);
latvec=cycleLat(1,:);
latvec=latvec(:);
lonvec=cycleLon(:,1);
[~,I1]=min(abs(lonvec-ODP849(2)));
[~,I2]=min(abs(latvec-ODP849(1)));
true_climatology=squeeze(cycled18O(I1,I2,:));
clear cycled18O cycleLat cycleLon d18O1 Lat Lon X3 Y3 Z3
%% Holocene
load('Depth_Distributions.mat')
benthic_d18O=1.75;
depths=[0,32,210,114,96,610];
%Get average profile
[beta,MLD]=run_thermocline_model3(depths,[ODP849(3:7),benthic_d18O]);

depth_plot=linspace(0,610,1000);
Hol_profile=depth_plot*nan;
Hol_profile(depth_plot<MLD)=ODP849(3);
log_func=@(depth,beta) -1*beta(1).^(-1.*(depth+beta(2)))+beta(3);
Hol_profile(depth_plot>=MLD)=log_func(depth_plot(depth_plot>=MLD),beta);
MC_num=500;
MC_profiles=ones([MC_num,length(depth_plot)])*nan;
%Get realizations of profile based on ACD error
for i=1:MC_num
    disp(i)
    MC_depths=[0];
    for j=2:5
        MC_depths(j)=depth_realization(xi_tot(j,:),ksdens(j,:),minmax(j,1),minmax(j,2));
    end
    MC_depths(6)=depths(end);
    [MC_beta,MC_MLD]=run_thermocline_model3(MC_depths,[ODP849(3:7),benthic_d18O]);
    MC_profiles(i,depth_plot<MC_MLD)=ODP849(3);
    MC_profiles(i,depth_plot>=MC_MLD)=log_func(depth_plot(depth_plot>=MC_MLD),MC_beta);
end

%% Holocene Figure code
col=[215/255,48/255,39/255];
figure()
hold on
plot([ODP849(3:7),benthic_d18O],depths(1:6),'ok','MarkerFaceColor',col)
plot(Hol_profile,depth_plot,'-','LineWidth',1.5,'Color',col)
plot(true_climatology(1:40),depthvec(1:40),'--k','LineWidth',1.5)

for i=1:MC_num
    plot(MC_profiles(i,:),depth_plot,'-','LineWidth',1.5,'Color',[col,5/(MC_num)])
end

set(gca,'YDir','reverse')
xlabel('\delta^{18}O_c')
ylabel('Depth (m)')
axis([-2,2,0,610])

%% LGM
load('Depth_Distributions.mat')
LGMODP849=[ODP849(1:2),-0.07,0.22,1.44,1.33,0.79];
LGMODP849(4)=nan;
LGMbenthic=2.797;
depths=[0,32,210,114,96,500];
%Get average profile
[beta,MLD]=run_thermocline_model3(depths,[LGMODP849(3:7),LGMbenthic]);

depth_plot=linspace(0,610,1000);
LGM_profile=depth_plot*nan;
LGM_profile(depth_plot<MLD)=LGMODP849(3);
log_func=@(depth,beta) -1*beta(1).^(-1.*(depth+beta(2)))+beta(3);
LGM_profile(depth_plot>=MLD)=log_func(depth_plot(depth_plot>=MLD),beta);
MC_num=500;
LGM_MC_profiles=ones([MC_num,length(depth_plot)])*nan;
%Get realizations of profile based on ACD error
for i=1:MC_num
    disp(i)
    LGM_MC_depths=[0];
    for j=2:5
        LGM_MC_depths(j)=depth_realization(xi_tot(j,:),ksdens(j,:),minmax(j,1),minmax(j,2));
    end
    LGM_MC_depths(6)=depths(end);
    [LGM_MC_beta,LGM_MC_MLD]=run_thermocline_model3(LGM_MC_depths,[LGMODP849(3:7),LGMbenthic]);
    LGM_MC_profiles(i,depth_plot<LGM_MC_MLD)=LGMODP849(3);
    LGM_MC_profiles(i,depth_plot>=LGM_MC_MLD)=log_func(depth_plot(depth_plot>=LGM_MC_MLD),LGM_MC_beta);
end

%%
Tshift=0.21*-2.5; %Applied to Holocene climatology to plot it with LGM data and results
col2=[69/255,117/255,180/255];
figure()
hold on
plot([LGMODP849(3:7),LGMbenthic],depths(1:6),'ok','MarkerFaceColor',col2)
plot(LGM_profile,depth_plot,'-','LineWidth',1.5,'Color',col2)
plot(true_climatology(1:40)+1-Tshift,depthvec(1:40),'--k','LineWidth',1.5)

for i=1:MC_num
    plot(LGM_MC_profiles(i,:),depth_plot,'-','LineWidth',1.5,'Color',[col2,5/(MC_num)])
end
set(gca,'YDir','reverse')
xlabel('\delta^{18}O_c')
ylabel('Depth (m)')
axis([-2+1-Tshift,2+1-Tshift,0,610])

%% Combined figure
figure()
hold on
%Plot Holocene results first
plot([ODP849(3:7)],depths(1:5),'ok','MarkerFaceColor',col)
plot(benthic_d18O,610,'ok','MarkerFaceColor',col)
plot(Hol_profile,depth_plot,'-','LineWidth',1.5,'Color',col)

plot(true_climatology(1:40),depthvec(1:40),'--k','LineWidth',1.5)

for i=1:MC_num
    plot(MC_profiles(i,:),depth_plot,'-','LineWidth',1.5,'Color',[col,5/(MC_num)])
end

set(gca,'YDir','reverse')
xlabel('\delta^{18}O_c')
ylabel('Depth (m)')

%Plot LGM results
Tshift=0.21*-2.5;
plot([LGMODP849(3:7),LGMbenthic]-1+Tshift,depths(1:6),'ok','MarkerFaceColor',col2)
plot(LGM_profile-1+Tshift,depth_plot,'-','LineWidth',1.5,'Color',col2)
plot(true_climatology(1:40),depthvec(1:40),'--k','LineWidth',1.5)

for i=1:MC_num
    plot(LGM_MC_profiles(i,:)-1+Tshift,depth_plot,'-','LineWidth',1.5,'Color',[col2,5/(MC_num)])
end
set(gca,'YDir','reverse')
xlabel('\delta^{18}O_c')
ylabel('Depth (m)')
axis([-2,2,0,610])
