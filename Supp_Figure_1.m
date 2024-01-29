addpath('Data')
addpath('Functions')
load('Z20.mat')
load('Depth_Distributions.mat')
sample=good_data(good_data(:,1)==0.83,:);
%Take ML1208-19GC for sample location
sample(4)=nan;
benthic_d18O=1.75;
%clear good_data goods
ACD=[18,32,210,114,96];
[opt_beta,opt_MLD]=run_thermocline_model3([ACD,617],[sample(3:7),1.6]);
beta_ranges=[1.0015,1.0105;      -50,-550;       20,140];
log_func=@(depth,beta) -1*beta(1).^(-1.*(depth+beta(2)))+beta(3);
%Function to get beta3 given other parameters
beta3=@(beta1,beta2,MLD) sample(3)+beta1^(-MLD-beta2);
num=20;

load('Data/d18O_paper.mat')
Lat=squeeze(X3(:,:,1));
Lon=squeeze(Y3(:,:,1));
depthvec=squeeze(Z3(1,1,:));
[cycleLon,cycleLat,cycled18O]=cycle_data_3D(Lon,Lat,d18O1,30);
latvec=cycleLat(1,:);
latvec=latvec(:);
lonvec=cycleLon(:,1);
s=size(sample);
[~,I1]=min(abs(lonvec-sample(2)));
[~,I2]=min(abs(latvec-sample(1)));
true_climatology=squeeze(cycled18O(I1,I2,:));

figure('Position',[66.6,52.2,982.4,710.8])
%Subplot 1: main thermocline model result using ACD and individual Monte Carlo runs of the
%thermocline model based on ACD uncertainty
subplot(2,2,1)
hold on
plot([sample(3:7),benthic_d18O],[ACD,610],'ok','Color','k','LineWidth',1.5)
plot(true_climatology,depthvec,'--k','LineWidth',1.5)
set(gca,'YDir','reverse')
depth_plot=linspace(0,610,1000);
[~,I1]=min(abs(depth_plot-opt_MLD));
profile=log_func(depth_plot,opt_beta);
profile(1:I1)=sample(3);
plot(profile,depth_plot,'-k','LineWidth',2)
MC_profiles=[];
MC_num=10;
for i=1:MC_num
    MC_depths=[0];
    for j=2:5
        MC_depths(j)=depth_realization(xi_tot(j,:),ksdens(j,:),minmax(j,1),minmax(j,2));
    end
    MC_depths(6)=610;
    [MC_beta,MC_MLD]=run_thermocline_model3(MC_depths,[sample(3:7),benthic_d18O]);
    MC_profiles(i,depth_plot<MC_MLD)=sample(3);
    MC_profiles(i,depth_plot>=MC_MLD)=log_func(depth_plot(depth_plot>=MC_MLD),MC_beta);    
    plot(MC_profiles(i,:),depth_plot,'-','LineWidth',1.5,'Color',[0,0,1,5/(MC_num)])
end
axis([-2.5,3,0,610])
xlabel('\delta^{18}O_c (‰)')
ylabel('Depth (m)')

for var=1:3
    %Subplots 2-4: Main thermocline model result using average ACD and
    %varying one of the parameters in the thermocline model function
    %(beta1,beta2,MLD)
    subplot(2,2,var+1)
    hold on
    plot(true_climatology,depthvec,'--k','LineWidth',1.5)
    plot([sample(3:7),benthic_d18O],[ACD,610],'ok','Color','k')
    set(gca,'YDir','reverse')
    depth_plot=linspace(0,610,1000);
    [~,I1]=min(abs(depth_plot-opt_MLD));
    profile=log_func(depth_plot,opt_beta);
    profile(1:I1)=sample(3);
    
    plot(profile,depth_plot,'-k','LineWidth',2)
    for i=1:num
        range=linspace(beta_ranges(var,1),beta_ranges(var,2),num);
        switch var
            case 1
                MLD=opt_MLD;
                beta=[range(i),opt_beta(2),beta3(range(i),opt_beta(2),MLD)];
            case 2
                MLD=opt_MLD;
                beta=[opt_beta(1),range(i),beta3(opt_beta(1),range(i),MLD)];
            case 3
                MLD=range(i);
                beta=[opt_beta(1:2),beta3(opt_beta(1),opt_beta(2),MLD)];
        end
        [~,I1]=min(abs(depth_plot-MLD));
        profile=log_func(depth_plot,beta);
        profile(1:I1)=sample(3);
        plot(profile,depth_plot,'Color',[0,0,1,0.45],'LineWidth',1)
    end
    axis([-2.5,3,0,610])
    xlabel('\delta^{18}O_c (‰)')
    ylabel('Depth (m)')
end