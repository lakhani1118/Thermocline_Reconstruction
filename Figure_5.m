%% Individual species minus benthic
addpath('Functions/')
addpath('Data/')
load('Data/Z20.mat','good_data','Z20')
benthic_d18O=1.75;

%Remove rows that do not have species data
goods=~any(isnan(good_data(:,[3,5,6,7])),2);
good_data=good_data(goods,[1,2,3,5,6,7]);Z20=Z20(goods);
clear goods

%Cores VM28-227 and VM28-229 excluded due to questionable d18O
%stratigraphy, see Karim Lakhani's Thesis, chapter 3.2.1
latitudes=[-8.4,-10.667];
for i=1:2
    row=good_data(:,1)==latitudes(i);
    good_data(row,:)=[];Z20(row)=[];
end
clear latitudes row i
%%
figure('Position',[0,-166.2,980,903.2])
names={'{\it G. tumida} - benthic  \delta^{18}O_c (‰)','{\it N. dutertrei} - benthic  \delta^{18}O_c (‰)','{\it P. obliquiloculata} - benthic  \delta^{18}O_c (‰)'};
for i=1:3
    subplot(2,2,i)
    hold on
    scatter(Z20,good_data(:,i+3)-benthic_d18O,36,'HandleVisibility','off')
    x=linspace(0,250,100);
    [p,S]=polyfit(Z20,good_data(:,i+3)-benthic_d18O,1);
    y=polyval(p,x);
    plot(x,y)
    [ccorr,P]=corrcoef(Z20,good_data(:,i+3)-benthic_d18O);
    P=P(1,2);
    correl(i)=ccorr(1,2);
    xlabel('Climatology 20°C isotherm')
    ylabel(names{i})
    if P<0.0001
        title(strcat('Correlation: ',num2str(round(ccorr(1,2),2)),"  P-value: <0.0001"))
    else
        title(strcat('Correlation: ',num2str(round(ccorr(1,2),2)),"  P-value: ",num2str(round(P,4))))
    end
end

subplot(2,2,4)
hold on
scatter(Z20,mean(good_data(:,[4,5,6]),2)-benthic_d18O)
x=linspace(0,250,100);
[p,S]=polyfit(Z20,mean(good_data(:,[4,5,6]),2)-benthic_d18O,1);
y=polyval(p,x);
plot(x,y)
xlabel('Climatology 20°C isotherm')
ylabel('Subsurface mean - benthic  \delta^{18}O_c(‰)')
[ccorr,P]=corrcoef(Z20,mean(good_data(:,[4,5,6]),2)-benthic_d18O);
P=P(1,2);
if P<0.0001
    title(strcat('Correlation: ',num2str(round(ccorr(1,2),2)),"  P-value: <0.0001"))
else
    title(strcat('Correlation: ',num2str(round(ccorr(1,2),2)),"  P-value: ",num2str(round(P,4))))
end
axis([0,250,ylim])
