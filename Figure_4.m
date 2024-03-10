%% Z20 equivalent depth and TP=0.8 depth
addpath('Functions/')
addpath('Data/')
load('WOA_04_temp.mat') 
%World Ocean Atlas 2013 0.25 degree mean annual temperature climatology

[Lon,Lat,Temps]=cycle_data_3D(Lon',Lat',Temps,30);
latvec=Lat(1,:);lonvec=Lon(:,1);
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
%% Get the TP=0.8 value for each core location
load('TP_80.mat')
i=1;
latvec=Lat(1,:);
lonvec=Lon(:,1);
for i=1:length(good_data)
    [~,I1]=min(abs(latvec-good_data(i,1)));
    [~,I2]=min(abs(lonvec-good_data(i,2)));
    TP80(i)=TP_data(I2,I1);
    if isnan(TP80(i))
        TP80(i)=TP_data(I2+1,I1+1);
    end
end
%% Calculate d18O TP=0.8 from climatology
load('d18O_paper.mat')
Lat=squeeze(X3(:,:,1));
Lon=squeeze(Y3(:,:,1));
depthvec=squeeze(Z3(1,1,:));
[cycleLon,cycleLat,cycled18O]=cycle_data_3D(Lon,Lat,d18O1,30);
latvec=cycleLat(1,:);
latvec=latvec(:);
lonvec=cycleLon(:,1);
%%
i=5;
Lat_Lon_shift=[3,3;1,1;0,3;0,3;-4,-9];
for i=1:length(good_data)
    [~,I1]=min(abs(latvec-good_data(i,1)));
    [~,I2]=min(abs(lonvec-good_data(i,2)));
    top=cycled18O(I2,I1,1);
    bot1000=cycled18O(I2,I1,47);
    bot650=cycled18O(I2,I1,40);
    profile=squeeze(cycled18O(I2,I1,1:47));
    if and(~isnan(top),~isnan(bot650))
        d18TP1000(i)=bot1000+0.8*(top-bot1000);
        d18TP650(i)=bot650+0.8*(top-bot650);
        z1000(i)=interp1(profile(1:47),depthvec(1:47),d18TP1000(i));
        z650(i)=interp1(profile(1:47),depthvec(1:47),d18TP650(i));
    end
    if i<6
        I1=I1+Lat_Lon_shift(i,2);
        I2=I2+Lat_Lon_shift(i,1);
        top=cycled18O(I2,I1,1);
        bot1000=cycled18O(I2,I1,47);
        bot650=cycled18O(I2,I1,40);
        profile=squeeze(cycled18O(I2,I1,1:47));
        d18TP1000(i)=bot1000+0.8*(top-bot1000);
        d18TP650(i)=bot650+0.8*(top-bot650);
        z1000(i)=interp1(profile(1:47),depthvec(1:47),d18TP1000(i));
        z650(i)=interp1(profile(1:47),depthvec(1:47),d18TP650(i));
    end
end
%z650(5)=nan;
%% Generate average profile, calculate TP=0.8 from model, and get average d18Osw and Z20 equivalent in d18Oc
load('Legrande_Schmidt.mat')
Pacific=or(LegrandeLon>150,LegrandeLon<-100);Tropics=abs(LegrandeLat)<20;subset=Legranded18Osw(Pacific,Tropics,9);subset=subset(:);subset(subset<-1e5)=[];
d18Osw=mean(subset);
Paleotemperature=@(Temp,d18Osw) (d18Osw-.27)-.2*Temp+3.25;
d18O_Z20=Paleotemperature(20,d18Osw);

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
%% Error bars for plot a 
load('Depth_Distributions.mat')
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

%% Plot Figure 3
figure()
hold on
yneg=abs(prctile(Z20_MC,2.5,2)-model_Z20_equivalent');
ypos=abs(prctile(Z20_MC,97.5,2)-model_Z20_equivalent');
errorbar(Z20,model_Z20_equivalent,yneg,ypos,"o")
maxmin=xlim;x=linspace(0,250,100);
[p,S]=polyfit(Z20,model_Z20_equivalent,1);
y=polyval(p,x);
plot(x,y)
xlabel('Climatology 20°C isotherm')
ylabel(strcat('Depth where \delta^{18}O = ',(num2str(round(d18O_Z20,2)))))
[ccorr,P]=corrcoef(Z20,model_Z20_equivalent);
P=P(1,2);
if P<0.0001
    title(strcat('Correlation: ',num2str(round(ccorr(1,2),2)),"  P-value: <0.0001"))
else
    title(strcat('Correlation: ',num2str(round(ccorr(1,2),2)),"  P-value: ",num2str(round(P,4))))
end
axis([0,250,0,250])
plot([0,250],[0,250],'--k')
disp(p)
