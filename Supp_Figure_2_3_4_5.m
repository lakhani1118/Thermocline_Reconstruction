%% 
addpath('Functions')
addpath('Data')
load('Z20.mat')
load('Depth_Distributions.mat')
% a=find(good_data(:,1)==-7.133);
% b=find(good_data(:,1)==0.83);
% c=find(good_data(:,1)==0.1830);
good_data(:,4)=good_data(:,4)*nan; %Remove T. sacculifer from code
MC_num=100;

load('d18O_paper.mat')
Lat=squeeze(X3(:,:,1));
Lon=squeeze(Y3(:,:,1));
depthvec=squeeze(Z3(1,1,:));
[cycleLon,cycleLat,cycled18O]=cycle_data_3D(Lon,Lat,d18O1,30);
latvec=cycleLat(1,:);
latvec=latvec(:);
lonvec=cycleLon(:,1);
true_climatology=[];
ACD=[18,32,210,114,96];
benthic_d18O=1.75;

s2=size(good_data);
goods=find(~any(isnan(good_data(:,[3,5,6,7])),2));
good_data2=good_data(goods,:);
% good_data2(good_data2(:,1)==good_data(a,1),:)=[];
% good_data2(good_data2(:,1)==good_data(b,1),:)=[];
% good_data2(good_data2(:,1)==good_data(c,1),:)=[];
s2=size(good_data2);

groups=[1,5; %WPWP cores
                6,9; %WP cores
                12,16; %CP cores
                17,22]; %EP cores
s2=size(groups);
%%
for h=2:s2(1)
    for i=groups(h,1):groups(h,2)
        depth_plot=linspace(0,610,1000);
        [~,I1]=min(abs(lonvec-good_data2(i,2)));
        [~,I2]=min(abs(latvec-good_data2(i,1)));
        if h==1
            %WPWP cores are too close to the mask, so have to be slightly
            %adjusted to get the closest climatology value
            shift=[3,3;
                       0,2;
                       5,0;
                       5,0;
                       -10,-5];
            I1=I1+shift(i,2);
            I2=I2+shift(i,1);
        else
            
        end
        true_climatology=squeeze(cycled18O(I1,I2,:));
        true_climatology=true_climatology+fliplr(1:length(true_climatology))*0.00000001;
    
    
        figure('Position',[731.4,317.8,529.6,445.2])
        hold on
        plot(true_climatology,depthvec,'--k','LineWidth',1.5)
        plot([good_data2(i,3:7),benthic_d18O],[ACD,610],'ok','Color','k')

        [beta,MLD]=run_thermocline_model3([ACD,610],[good_data2(i,3:7),benthic_d18O]);
        profile=depth_plot*nan;
        profile(depth_plot<MLD)=good_data2(i,3);
        log_func=@(depth,beta) -1*beta(1).^(-1.*(depth+beta(2)))+beta(3);
        profile(depth_plot>=MLD)=log_func(depth_plot(depth_plot>=MLD),beta);
        profile_with(i,:)=profile;

        plot(profile_with(i,:),depth_plot,'-b','LineWidth',1.5)
        
        axis([-3.5,2,0,610])
        set(gca,'YDir','reverse')
        MC_with=[];
        for k=1:MC_num
            MC_depths=[0];
            for j=2:5
                MC_depths(j)=depth_realization(xi_tot(j,:),ksdens(j,:),minmax(j,1),minmax(j,2));
            end
            MC_depths(6)=610;
            [MC_beta,MC_MLD]=run_thermocline_model3(MC_depths,[good_data2(i,3:7),benthic_d18O]);
            MC_with(depth_plot<MC_MLD)=good_data2(i,3);
            MC_with(depth_plot>=MC_MLD)=log_func(depth_plot(depth_plot>=MC_MLD),MC_beta);
            plot(squeeze(MC_with(:)),depth_plot,'-','LineWidth',1,'Color',[0,0,1,5/MC_num]);
        end
        
        n=find(isnan(true_climatology),5,'first')-1;
        Z20_equiv=interp1(true_climatology(1:n),depthvec(1:n),-0.66);
        plot([-0.66,-0.66],[Z20_equiv-100,Z20_equiv+100],'--','Color',0.15*[1,1,1],'LineWidth',1)
        
        %title(strcat("Lat: ",num2str(good_data2(i,1))," ; Lon: ",num2str(good_data2(i,2))));
        xlabel('\delta^{18}O_c (‰)')
        ylabel('Depth (m)')
    end
    pause()
end
