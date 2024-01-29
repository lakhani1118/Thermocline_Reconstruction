function [depth]=depth_realization(xi,f,mini,maxi)
k=0;
spacing=xi(2)-xi(1);
cumulate=cumsum(f.*spacing);
while k==0
    a=rand;
    [~,I]=min(abs(cumulate-a));
    depth=xi(I);
    if and(depth >=mini,depth<=maxi)
        k=1;
    end
end


end