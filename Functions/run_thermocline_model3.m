function [beta, MLD]=run_thermocline_model3(depths,profile)
    %%% depths is a vector with placeholder MLD depth, depths for subsurface
    %%% species used in model fitting, and the depth of the benthic point. profile 
    %%% is a vector of the MLD d18O (from surface species), the d18O of the 
    %%% subsurface species, and the d18O of the benthic point. 

    %Functional form of the profile below the mixed layer
    log_func=@(depth,beta) -1*beta(1).^(-1.*(depth+beta(2)))+beta(3);

    %Evaluate every parameter combination and the error resulting from it.
    %Coarse search
    beta1_range=linspace(1.0015,1.0105,20);
    beta2_range=linspace(-50,-550,20);
    MLD_range=linspace(20,140,9);
    RMS=ones(length(beta1_range),length(beta2_range),length(MLD_range))*nan;
    for i=1:length(beta1_range)
        for j=1:length(beta2_range)
            for m=1:length(MLD_range)
                depth_profile=linspace(0,max(depths)*1.05,1000);
                d18O_profile=depth_profile.*nan;
                d18O_profile(depth_profile<=MLD_range(m))=profile(1);
                                
                beta3=profile(1)+beta1_range(i)^(-MLD_range(m)-beta2_range(j));
                %Given beta1, beta2, and MLD and continuity at the base of
                %the mixed layer, beta3 can be determined
                beta=[beta1_range(i), beta2_range(j), beta3];
                d18O_profile(depth_profile>MLD_range(m))=log_func(depth_profile(depth_profile>MLD_range(m)),beta);
                
                pred(1)=profile(1); %Make placeholder not affect RMS error
                for k=2:length(depths)
                    [~,I1]=min(abs(depth_profile-depths(k)));
                    frac=-(depths(k)-depth_profile(mindata(I1-1)))/(depth_profile(I1+1)-depth_profile(mindata(I1-1)));
                    pred(k)=d18O_profile(mindata(I1-1))+frac*(d18O_profile(mindata(I1-1))-d18O_profile(I1+1));
                end
                %Calculate error and store
                RMS(i,j,m)=nanmean((pred-profile).^2);
            end
        end
    end
    %Find minimum over this parameter set
    s2=size(RMS);
    [~,I]=min(RMS(:));[x,y,a]=ind2sub(s2,I);
    
    %Find where the parameter combination with the lowest error must be,
    %define a new grid and redo the exhaustive search for the minimum error
    %parameter combination
    %Fine search
    beta1_range=logspace(log10(beta1_range(mindata(x-1))),log10(beta1_range(maxdata(x+1,length(beta1_range)))),20);
    beta2_range=linspace(beta2_range(mindata(y-1)),beta2_range(maxdata(y+1,length(beta2_range))),20);
    MLD_range=linspace(MLD_range(mindata(a-1)),MLD_range(maxdata(a+1,length(MLD_range))),6);
    RMS=ones(length(beta1_range),length(beta2_range),length(MLD_range))*nan;
    for i=1:length(beta1_range)
        for j=1:length(beta2_range)
            for m=1:length(MLD_range)
                depth_profile=linspace(0,max(depths)*1.05,1000);
                d18O_profile=depth_profile.*nan;
                d18O_profile(depth_profile<=MLD_range(m))=profile(1);
                                
                beta3=profile(1)+beta1_range(i)^(-MLD_range(m)-beta2_range(j));
                %Given beta1, beta2, and MLD and continuity at the base of
                %the mixed layer, beta3 can be determined
                beta=[beta1_range(i), beta2_range(j), beta3];
                d18O_profile(depth_profile>MLD_range(m))=log_func(depth_profile(depth_profile>MLD_range(m)),beta);
                
                pred(1)=profile(1); %Make placeholder not affect RMS error
                for k=2:length(depths) %Get prediction for each depth based on current beta parameters
                    [~,I1]=min(abs(depth_profile-depths(k)));
                    frac=-(depths(k)-depth_profile(mindata(I1-1)))/(depth_profile(I1+1)-depth_profile(mindata(I1-1)));
                    pred(k)=d18O_profile(mindata(I1-1))+frac*(d18O_profile(mindata(I1-1))-d18O_profile(I1+1));
                end
                %Calculate error and store
                RMS(i,j,m)=nanmean((pred-profile).^2);
            end
        end
    end
    %Find minimum over this parameter set
    s2=size(RMS);
    [~,I]=min(RMS(:));[x,y,a]=ind2sub(s2,I);
    beta3=profile(1)+beta1_range(x)^(-MLD_range(a)-beta2_range(y));
    %Output results
    beta=[beta1_range(x), beta2_range(y), beta3];
    MLD=MLD_range(a);

end