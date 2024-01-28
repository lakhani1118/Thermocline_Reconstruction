function [newlon,newlat,newdata]=cycle_data(Lon,Lat,data,pos)
s=size(Lon);
[A,I]=min(abs(pos-Lon(:,1)));

cyclelon=Lon(1:I,:);
cyclelat=Lat(1:I,:);
cycledata=data(1:I,:,:);

cyclelon=cyclelon+360;
d=size(cyclelon);

newlon(1:s(1)-d(1),:)=Lon(I+1:s(1),:);
newlon(s(1)-d(1)+1:s(1),:)=cyclelon(:,:);
newlat(1:s(1)-d(1),:)=Lat(I+1:s(1),:);
newlat(s(1)-d(1)+1:s(1),:)=cyclelat(:,:);
newdata(1:s(1)-d(1),:,:)=data(I+1:s(1),:,:);
newdata(s(1)-d(1)+1:s(1),:,:)=cycledata(:,:,:);

end
