function [value]=maxdata(value,max)
%If value goes above the size of the dataset, change it to the maximum
%value
if value>max
    value=max;
end
end