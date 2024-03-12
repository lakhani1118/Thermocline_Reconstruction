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

load('d18O_paper.mat')
Lat=squeeze(X3(:,:,1));
Lon=squeeze(Y3(:,:,1));
depthvec=squeeze(Z3(1,1,:));
[cycleLon,cycleLat,cycled18O]=cycle_data_3D(Lon,Lat,d18O1,30);
clear d18O1 Lat Lon
latvec=cycleLat(1,:);
latvec=latvec(:);
lonvec=cycleLon(:,1);
s=size(sample);
true_climatology=[];

for i=1:3
    [~,I1]=min(abs(lonvec-sample(i,2)));
    [~,I2]=min(abs(latvec-sample(i,1)));
    true_climatology(i,:)=squeeze(cycled18O(I1,I2,:));
    true_climatology(i,:)=true_climatology(i,:)+fliplr(1:length(true_climatology(i,:)))*0.00000001;
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
            
figure('Position',[490.6,85.8,981.6,886.4])
corenames={'VM28-234','ML1208-19GC','ODP 849'};

for i=1:s(1)
    color=col_list(i,:);
    subplot(3,3,i)
    hold on
    plot(true_climatology(i,:),depthvec,'--k','LineWidth',1.5)
    plot([sample(i,3:7)],[ACD],'ok','Color','k')
    plot(profile_without(i,:),depth_plot,'-','LineWidth',1.5,'Color',color)
    
    n=find(isnan(true_climatology(i,:)),5,'first')-1;
    Z20_equiv=interp1(true_climatology(i,1:n),depthvec(1:n),-0.66);
    plot([-0.66,-0.66],[Z20_equiv-100,Z20_equiv+100],'--','Color',0.15*[1,1,1],'LineWidth',1)
    
    for j=1:MC_num
        plot(squeeze(MC_without(i,j,:)),depth_plot,'-','LineWidth',1,'Color',[color,5/(MC_num)])
    end
    
    xlabel('\delta^{18}O_c (‰)')
    ylabel('Depth (m)')
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
    
    n=find(isnan(true_climatology(i,:)),5,'first')-1;
    Z20_equiv=interp1(true_climatology(i,1:n),depthvec(1:n),-0.66);
    plot([-0.66,-0.66],[Z20_equiv-100,Z20_equiv+100],'--','Color',0.15*[1,1,1],'LineWidth',1)
    
    xlabel('\delta^{18}O_c (‰)')
    ylabel('Depth (m)')
    axis([-3,2,0,610])
    set(gca,'YDir','reverse')
end

%% Plot contour of d18O for East-West Pacific between 5°S and 5°N
bott=40;
[~,I1]=min(abs(latvec+5));
[~,I2]=min(abs(latvec-5));
XZ_d18O=squeeze(nanmean(cycled18O(:,I1:I2,1:bott),2));
[Depth,Long]=meshgrid(depthvec(1:bott),lonvec);
for i=1:length(depthvec(1:bott))
    XZ_d18O(598:603,i)=interp1([Long(598,i),Long(603,i)],[XZ_d18O(598,i),XZ_d18O(603,i)],Long(598:603,i));
end


Pacific=[150,281];
[~,I3]=min(abs(lonvec-Pacific(1)));
[~,I4]=min(abs(lonvec-Pacific(2)));

subplot(3,1,3)
[~,hContour] = contourf(Long(I3:I4,:),Depth(I3:I4,:),XZ_d18O(I3:I4,:),10,'EdgeColor','none');
drawnow
hFills=hContour.FacePrims;
num=120;
eventFcn = @(srcObj, e) updateTransparency(srcObj,num);
addlistener(hContour, 'MarkedClean', eventFcn);

xticks([160, 200, 240, 280]);
xticklabels({'160\circE','160\circW','120\circW','80\circW'});
set(gca,'YDir','reverse')
colormap(colorblind_jet)
h=colorbar;
set(h,'YDir','reverse')

%% Generate average profile, calculate TP=0.8 from model, and get average d18Osw and Z20 equivalent in d18Oc

depths=[0,nan,210,114,96,610];
benthic_d18O=1.75;
log_func=@(depth,beta) -1*beta(1).^(-1.*(depth+beta(2)))+beta(3);
for i=1:length(good_data)
    [beta,MLD]=run_thermocline_model3(depths,[good_data(i,3:7),benthic_d18O]);
    depth_plot=linspace(0,610,1000);
    profile=depth_plot*nan;
    profile(depth_plot<MLD)=good_data(i,3);
    profile(depth_plot>=MLD)=log_func(depth_plot(depth_plot>=MLD),beta);
    mean_profile(i,:)=profile;
end

%% Generate profiles for Holocene data
%EP_profile
EP=find(and(abs(good_data(:,1))<5,good_data(:,2)>210));
EP_average=mean(mean_profile(EP,:),1);

width=2;
EP_loc=260;
xvec=linspace(EP_loc-width,EP_loc+width,3);
[X,Y]=meshgrid(xvec,depth_model);
%Z=squeeze(repmat(EP_average,[3,1]))';
Z=repmat(EP_average,[3,1])';
hold on
contourf(X,Y,Z,10,'EdgeColor','none')

%CP_profile
CP=find(and(and(abs(good_data(:,1))<5,good_data(:,2)<210),good_data(:,2)>180));
CP_average=mean(mean_profile(CP,:),1);

CP_loc=200;
xvec=linspace(CP_loc-width,CP_loc+width,3);
[X,Y]=meshgrid(xvec,depth_model);
Z=squeeze(repmat(CP_average,[3,1]))';
%Z=repmat(CP_average,[3,1])';
hold on
contourf(X,Y,Z,10,'EdgeColor','none')

%WP_profile
WP=find(and(abs(good_data(:,1))<5,good_data(:,2)<180));
WP_average=mean(mean_profile(WP,:),1);

WP_loc=155;
xvec=linspace(WP_loc-width,WP_loc+width,3);
[X,Y]=meshgrid(xvec,depth_model);
Z=squeeze(repmat(WP_average,[3,1]))';
%Z=repmat(WP_average,[3,1])';
hold on
contourf(X,Y,Z,10,'EdgeColor','none')