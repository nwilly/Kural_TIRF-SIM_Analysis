%{
For each image in imgs,  fit the pixel intensity vs the radial distance
from the structure's center to a gaussian.  Return whether the image is
of a donut structure based on the position of the Gaussian peak position.

inputs:
    imgs: cell array of matricies representing TIRF-SIM images
    cx: vector of the central position of the structures in imgs
    cy: vector of the central position of the structures in imgs
    x_grid: vector of x coordinates of each pixel in imgs
    y_grid: vector of y corrdinates of each pixel in imgs

outputs:
    donut: vector of 1 or 0 indicating whether each structure in imgs is a
            donut.

Author: Nathan Willy
    nmwilly@gmail.com
%}
function [donut] = g1fit(imgs,cx,cy,x_grid,y_grid)

    r_thresh_max = 6; % structure is not a donut if peak > r_thresh_max
    r_thresh_min = 2; % structure is not a donut if peak < r_thresh_min
    
    donut = zeros(length(imgs),1,'single'); % preallocate output

    %loop through all images
    for(i=1:length(imgs))

        % calculate radial distance of each pixel from [cx,cy]
        r = double(sqrt((y_grid-cy(i)).^2 + (x_grid-cx(i)).^2));
        % get pixel intensity values
        int = double(imgs{i}(:));

        % fit [r,int] to a Gaussian distribution
        f = fit(r,int,'gauss1','lower',[0,0,2],'upper',[inf,10,10]);
        r_max = f.b1*sqrt(2);
        
        % exclude structures found to be dim (small gaussian ampliture)
        % or flat (large standard deviation)
        if(f.a1<1000 || f.c1>10)
            continue;
        end
         
        % exclude structures where the peak is not within the interval
        if(r_max<r_thresh_min | r_max>r_thresh_max)
            continue;
        end
        
        % the remaining structures are donuts
        donut(i)=1;
        
    end
end