%{
For each image in imgs fit its average cross-sectional profile to a pair of
Gaussians.

inputs:
    imgs: cell array of matrices or a single matrix representing an
        TRIF-SIM image
    pix: pixel size in nm
    cx: array of x positions marking the center of the structure
    cy: array of y positions marking the center of the structure

outputs:
    results: data structure containing the following fields for each input
        image.
        profile: the calculated average cross-sectional profile
        pos: the radial position of each point in 'profile' in nm
            ie to plot(pos,profile) to plot the profile
        gfit: the fit object from fitting the profile to the sum of two
            Gaussians.
        gof: the R2 coefficient of the fit
        
Author: Nathan Willy (willy.2@osu.edu)
%}

function results = donut_img_analysis(imgs,pix,cx,cy)

scale1 = 10;
div = 5*pi/180;
results = struct('profile',[],'pos',[],'fit',[],'mask',[],'area',[]);

if(~iscell(imgs))
    temp = {};
    temp{1} = imgs;
    imgs = temp;
end

% loop through the input images
for(i=1:length(imgs))
    % scale the image by factor scale1
    curr1 = double(imresize(imgs{i},scale1,'bicubic'));
    
    % if cx is not provoided calculate the central position of the img
    % in either case convert that position to a position on the scaled up
    % image
    if(isempty(cx))
        cx2 = sum(sum(curr1,2).*[1:size(curr1,2)]')/sum(curr1(:));
        cy2 = sum(sum(curr1,1).*[1:size(curr1,1)])/sum(curr1(:));
    else

        cx2 = double((cx(i)-1)/(length(imgs{1})-1)*(length(imgs{1})*scale1-1)+1);
        cy2 = double((cy(i)-1)/(length(imgs{1})-1)*(length(imgs{1})*scale1-1)+1);
    end

    
    %calculate the average profile by averaging the image cross-sections
    %running through the central point but at various angles.
    profile = zeros(min(size(curr1)),1);
    r = floor(min(size(curr1))/2);
    
    count=0;
    for(th = 0:div:pi-div/10)
        [px,py,p] = improfile(curr1,...
            [cy2 - r*sin(th), cy2 + r*sin(th)],...
            [cx2 - r*cos(th), cx2 + r*cos(th)],...
            length(profile),'bicubic');
        profile = profile + p;
        count=count+1;
        
        if(th==0)
            pos = py-cy2;
        end
    end
    
    profile = profile/count;
    
    
    % fit the profile to the sum of two Gassuians which have a common
    % standard deviation
    F = @(back, amp1, amp2, x1, x2, std, x) ...
        back + ...
        amp1*exp(-((x-x1).^2/(2*std^2))) + ...
        amp2*exp(-((x-x2).^2/(2*std^2)));

    x = pos*pix/scale1;
    y = profile;
    
    x = x(~isnan(y));
    y = y(~isnan(y));

    c0 =  double([min(y),  max(y),     max(y),    quantile(x,1/3), quantile(x,2/3), length(x/10)]);
    low = double([0,       mean(y),    mean(y),   x(1),            quantile(x,1/2), 0]);
    up =  double([mean(y), 1.2*max(y), 1.2*max(y) quantile(x,1/2)  x(end),          inf]);

    
    
    [gfit,gof,out] = fit(x, y, F, 'StartPoint', c0, 'Lower', low, 'Upper', up);
    
    
    results(i).profile = profile;
    results(i).pos = pos*pix/scale1;
    results(i).fit = gfit;
    results(i).gof = gof.rsquare;
    
    
end

end
