%{

Determines whether each SIM image found in tracking_data is a hollow, donut
structure or not. This is done in the function g1fit.

inputs:
    tracking_data: structure containing 11x11 SIM reconstructed images
    cropped around he center of a structure.

outputs:
    tracking_data: input structure with the added field donut, which is a
    boolean list of whether the structure is a donut in each of its frames.

Author: Nathan Willy (willy.2@osu.edu)

%}

function tracking_data = is_donut(tracking_data)

[x,y] = meshgrid(1:11,1:11);
x = double(x(:));
y = double(y(:));

num_imgs = length(vertcat([tracking_data.imgs]));
step = 500/1;

tic
for(i=1:step:length(tracking_data))
    disp([num2str(i),': ',num2str(length(vertcat(tracking_data.donut))),'/',num2str(num_imgs)])
    toc
    i2 = min(length(tracking_data),i+step);
    parfor(pari=i:i2)
        
        % if donut has already been calculated, skip it
        if(length(tracking_data(pari).donut)==length(tracking_data(pari).imgs))
            continue;
        end
        
        cx = 6 +  tracking_data(pari).xpos -  round(tracking_data(pari).xpos);
        cy = 6 +  tracking_data(pari).ypos -  round(tracking_data(pari).ypos);

        donut2 = g1fit(tracking_data(pari).imgs,cx,cy,x,y);

        for(j=1:length(donut2))
            tracking_data(pari).donut(j,1) = donut2(j);

        end

    end
end