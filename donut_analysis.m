%{

Perform a 2 Gaussian fit on the average cross-sectional profile of all
structure images in tracking_data. Return tracking_data with the fit
information appended.  The fit is performed in donut_img_analysis.

inputs:
    tracking_data: tracking data structures with sub images
   
outputs:
    tracking_data: tracking data structures with sub images and fit
        information appended.

Author: Nathan Willy (willy.2@osu.edu)
%}

function tracking_data = donut_analysis(tracking_data)

tic
tracking_data(end).g2fit=[];
tracking_data(end).gof=[];
num_imgs = length(vertcat([tracking_data.imgs]))
step = 500;

for(i=1:step:length(tracking_data))
    disp([num2str(i),': ',num2str(length(vertcat([tracking_data.g2fit]))),'/',num2str(num_imgs)])
    i2 = min(length(tracking_data),i+step);
    parfor(pari=i:i2)

        if(length(tracking_data(pari).imgs) == length(tracking_data(pari).g2fit))
            continue;
        end


        cx = double(6 + (tracking_data(pari).xpos - round(tracking_data(pari).xpos)));
        cy = double(6 + (tracking_data(pari).ypos - round(tracking_data(pari).ypos)));

        analysis = donut_img_analysis(tracking_data(pari).imgs,33.1,cx,cy);

        tracking_data(pari).g2fit = {};
        tracking_data(pari).gof = [];
        for(j=1:length(analysis))
            tracking_data(pari).g2fit{j} = single(coeffvalues(analysis(j).fit));
            tracking_data(pari).gof(j) = single(analysis(j).gof);
        end
    end
    disp([num2str(toc/60),' min'])
end
disp('done')