addpath('Functions')
addpath('Data')
load('Z20.mat')
load('Depth_Distributions.mat')
a=find(good_data(:,1)==-7.133);
b=find(good_data(:,1)==0.83);
c=find(good_data(:,1)==0.1830);
good_data(:,4)=good_data(:,4)*nan; %Remove T. sacculifer from code
sample=good_data([a,b,c],:); %Extract only the 3 example locations
MC_num=10 ;

load('d18O_paper.mat')
Lat=squeeze(X3(:,:,1));
Lon=squeeze(Y3(:,:,1));
depthvec=squeeze(Z3(1,1,:));
[cycleLon,cycleLat,cycled18O]=cycle_data_3D(Lon,Lat,d18O1,30);
latvec=cycleLat(1,:);
latvec=latvec(:);
lonvec=cycleLon(:,1);
s=size(sample);
true_climatology=[];

for i=1:3
    [~,I1]=min(abs(lonvec-sample(i,2)));
    [~,I2]=min(abs(latvec-sample(i,1)));
    true_climatology(i,:)=squeeze(cycled18O(I1,I2,:));
    %Get the predicted climatology that the model is compared to
end
ACD=[18,32,210,114,96];
benthic_d18O=1.75;

%With includes benthic data, without excludes benthic data
profile_with=[];
profile_without=[];
MC_with=[];
MC_without=[];
for i=1:3 %For each basin, get profiles that will be plotted
    %Get profile for average ACD for each species
    [beta,MLD]=run_thermocline_model3([ACD,610],[sample(i,3:7),benthic_d18O]);
    depth_plot=linspace(0,610,1000);
    profile=depth_plot*nan;
    profile(depth_plot<MLD)=sample(i,3);
    log_func=@(depth,beta) -1*beta(1).^(-1.*(depth+beta(2)))+beta(3);
    profile(depth_plot>=MLD)=log_func(depth_plot(depth_plot>=MLD),beta);
    profile_with(i,:)=profile; %Including benthic data
    
    %Get profiles for number of MC runs for random realization of ACD for
    %each species
    for k=1:MC_num 
            MC_depths=[0];
        for j=2:5
            MC_depths(j)=depth_realization(xi_tot(j,:),ksdens(j,:),minmax(j,1),minmax(j,2));
        end
        MC_depths(6)=610;
        [MC_beta,MC_MLD]=run_thermocline_model3(MC_depths,[sample(i,3:7),benthic_d18O]);
        MC_with(i,k,depth_plot<MC_MLD)=sample(i,3);
        MC_with(i,k,depth_plot>=MC_MLD)=log_func(depth_plot(depth_plot>=MC_MLD),MC_beta);
    end 
    
    %Get profile for average ACD for each species
    [beta,MLD]=run_thermocline_model3([ACD],[sample(i,3:7)]);
    depth_plot=linspace(0,610,1000);
    profile=depth_plot*nan;
    profile(depth_plot<MLD)=sample(i,3);
    log_func=@(depth,beta) -1*beta(1).^(-1.*(depth+beta(2)))+beta(3);
    profile(depth_plot>=MLD)=log_func(depth_plot(depth_plot>=MLD),beta);
    profile_without(i,:)=profile; %Excluding benthic data
    
    %Get profiles for number of MC runs for random realization of ACD for
    %each species
    for k=1:MC_num
        MC_depths=[0];
        for j=2:5
            MC_depths(j)=depth_realization(xi_tot(j,:),ksdens(j,:),minmax(j,1),minmax(j,2));
        end
        %MC_depths(6)=610;
        [MC_beta,MC_MLD]=run_thermocline_model3(MC_depths,[sample(i,3:7)]);
        MC_without(i,k,depth_plot<MC_MLD)=sample(i,3);
        MC_without(i,k,depth_plot>=MC_MLD)=log_func(depth_plot(depth_plot>=MC_MLD),MC_beta);
    end
end
%% Plot the results
col_list=[215/255,48/255,39/255;
                27/255, 120/255, 55/255;
                69/255,117/255,180/255];
            
figure('Position',[43.4,78.6,1217.6,684.4])
corenames={'VM28-234','ML1208-19GC','ODP 849'};

for i=1:s(1)
    color=col_list(i,:);
    subplot(3,3,i)
    hold on
    plot(true_climatology(i,:),depthvec,'--k','LineWidth',1.5)
    plot([sample(i,3:7)],[ACD],'ok','Color','k')
    plot(profile_without(i,:),depth_plot,'-','LineWidth',1.5,'Color',color)
    
    for j=1:MC_num
        plot(squeeze(MC_without(i,j,:)),depth_plot,'-','LineWidth',1,'Color',[color,5/(MC_num)])
    end
    
    xlabel('\delta^{18}O_c (‰)')
    ylabel('Depth (m)')
    %title(strcat("Lat: ",num2str(sample(i,1)),"; Lon: ",num2str(sample(i,2))))
    title(corenames{i})
    axis([-3,2,0,610])
    set(gca,'YDir','reverse')
    
    subplot(3,3,i+3)
    hold on
    plot(true_climatology(i,:),depthvec,'--k','LineWidth',1.5)
    plot([sample(i,3:7),benthic_d18O],[ACD,610],'ok','Color','k')
    plot(profile_with(i,:),depth_plot,'-','Color',color,'LineWidth',1.5)
    
    for j=1:MC_num
        plot(squeeze(MC_with(i,j,:)),depth_plot,'-','LineWidth',1,'Color',[color,5/(MC_num)])
    end
    
    xlabel('\delta^{18}O_c (‰)')
    ylabel('Depth (m)')
    axis([-3,2,0,610])
    set(gca,'YDir','reverse')
end


subplot(3,4,10)
hold on
for i=1:3
    color=col_list(i,:);
    plot(profile_with(i,:),depth_plot,'-','Color',color,'LineWidth',2)
    
end
xlabel('\delta^{18}O_c (‰)')
ylabel('Depth (m)')
axis([-3,2,0,610])
set(gca,'YDir','reverse')
legend({'Western Pacific','Central Pacific','Eastern Pacific'},'Location','SouthWest')

subplot(3,4,11)
hold on
for i=1:3
    color=col_list(i,:);
    plot(true_climatology(i,:),depthvec,'--','LineWidth',2,'Color',color)
end
xlabel('\delta^{18}O_c (‰)')
ylabel('Depth (m)')
axis([-3,2,0,610])
set(gca,'YDir','reverse')