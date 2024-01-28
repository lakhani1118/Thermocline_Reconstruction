function [value]=mindata(value)
%If value goes below the size of the dataset, change it to 1
if value<=0
    value=1;
end
end
