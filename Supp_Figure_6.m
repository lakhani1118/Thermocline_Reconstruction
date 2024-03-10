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
%% Calculate d18O TP=0.8 from climatology
load('d18O_paper.mat')
Lat=squeeze(X3(:,:,1));
Lon=squeeze(Y3(:,:,1));
depthvec=squeeze(Z3(1,1,:));
[cycleLon,cycleLat,cycled18O]=cycle_data_3D(Lon,Lat,d18O1,30);
latvec=cycleLat(1,:);
latvec=latvec(:);
lonvec=cycleLon(:,1);
Lat_Lon_shift=[3,3;1,1;0,3;0,3;-4,-9];
for i=1:length(good_data)
    [~,I1]=min(abs(latvec-good_data(i,1)));
    [~,I2]=min(abs(lonvec-good_data(i,2)));
    top=cycled18O(I2,I1,1);
    bot650=cycled18O(I2,I1,40);
    profile=squeeze(cycled18O(I2,I1,1:47));
    if and(~isnan(top),~isnan(bot650))
        z20_equiv(i)=interp1(profile(1:47),depthvec(1:47),-0.66);
    end
    if i<6
        I1=I1+Lat_Lon_shift(i,2);
        I2=I2+Lat_Lon_shift(i,1);
        profile=squeeze(cycled18O(I2,I1,1:47));
        z20_equiv(i)=interp1(profile(1:47),depthvec(1:47),-0.66);
    end
end
%% Plot z20 vs z where d18O=-0.66
%figure('Position',[0,0,980,903.2])
figure()
hold on
scatter(Z20,z20_equiv,"o")
maxmin=xlim;x=linspace(0,250,100);
[p,S]=polyfit(Z20,z20_equiv,1);
y=polyval(p,x);
plot(x,y)
xlabel('Climatology 20°C isotherm')
ylabel(strcat('Depth where \delta^{18}O = ',(num2str(round(-0.66,2)))))
[ccorr,P]=corrcoef(Z20,z20_equiv);
P=P(1,2);
if P<0.0001
    title(strcat('Correlation: ',num2str(round(ccorr(1,2),2)),"  P-value: <0.0001"))
else
    title(strcat('Correlation: ',num2str(round(ccorr(1,2),2)),"  P-value: ",num2str(round(P,4))))
end
axis([0,250,0,250])
plot([0,250],[0,250],'--k')
disp(p)