%{
Calculates the area of each structure and appends the results to
tracking_data.  The area is calculated by calculting the points of maxmimum
and minimum slope in the average cross-sectional profile.  The area is then
given as the area of a circle with diameter equal to the distance between
these points.

inputs:
    tracking_data, data structure containing cropped image and fit data
    pixel, pixel size in nm

outputs:
    tracking_data, input structure with area and max_rmax fields appended
        area, calculated area of the structure in nm^2
        max_rmax, the maximum observed radial extent of the structure

Author: Nathan Willy (willy.2@osu.edu)


%}

function tracking_data = calc_dcanny(tracking_data, pixel)

analysis = donut_img_analysis(tracking_data(1).imgs{1},pixel,6,6);
g2fit = analysis(1).fit;

tracking_data(1).area = [];
tracking_data(1).max_rmax = 0;
parfor(i=1:length(tracking_data))
    if(mod(i,100)==0)
        disp(num2str(i))
    end
    
    F = fittype('back + amp1*exp(-((x-x1).^2/(2*std^2))) + amp2*exp(-((x-x2).^2/(2*std^2)))');

    if(isempty(tracking_data(i).g2fit))
        disp('skipped')
        continue;
        
    end
    
    for(j=1:length(tracking_data(i).imgs))
        
        temp = tracking_data(i).g2fit{j};
        
        g2fit = cfit(F,temp(2),temp(3),temp(1),temp(6),temp(4),temp(5));
                
        X = [-250:0.1:250];
        Y = g2fit(X);

        x1 = min(g2fit.x1,g2fit.x2);
        x2 = max(g2fit.x1,g2fit.x2);

        Left = find(X<x1);
        Right = find(X>x2);
        [dY] = differentiate(g2fit,X);
        [MY,Mi] = max(dY(Left));
        [mY,mi] = min(dY(Right));
        
        tracking_data(i).area(j) = pi/4*(X(Right(mi))-X(Left(Mi)))^2;
    end
    tracking_data(i).max_rmax = sqrt(max(tracking_data(i).area)*4/pi);
end