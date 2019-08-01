%{ 
Remove structures from analysis based on the following criteria:

%1: trace is fewer than 3 frames
%2: there are fewer than 2 donut frames
%3: the first frame of the trace is a donut
%4: the first donut is more than 60s from the end of the trace
%5: fewer than 50% of frames after the first donut are donuts
%6: the profile fit is quite bad on at least 1 frame and/or sorta bad on
    %first donut frame
%7: there is high intensity pixel at the edge of an image

the matrix categories (not returned) enumerates the reasons for a structure
being excluded

inputs:
    tracking_data: data structure containing tracking and sub_image data
    frame_rate: time between frames (seconds)

outputs:
    tracking_data: tracking_data data structure with offending structures
            removed.

Author: Nathan Willy (willy.2@osu.edu)
%}
function tracking_data = select_good_structures(tracking_data, frame_rate)

list = [];

% keep track of which of the 7 criteria is responsible for a structure's
% exclusion
categories = zeros(length(tracking_data),7,'logical');

% loop through the structures
for(i=1:length(tracking_data))

    %1: trace is fewer than 3 frames
    if(length(tracking_data(i).imgs)<3)
        categories(i,1)=1;
    end
    
    %2: there are fewer than 2 donut frames
    is = find(tracking_data(i).donut);
    if(length(is)<2)
        categories(i,2)=1;
        continue
    end
    
    %3: the first frame of the trace is a donut
    if(is(1)==1)
        categories(i,3)=1;
    end
    
    %4: the first donut is more than 60s from the end of the trace
    if(length(tracking_data(i).donut)-is(1)*frame_rate(tracking_data(i).movie)>60)
        categories(i,4)=1;
    end
    
    %5: fewer than 50% of frames after the first donut are donuts
    subset = max(1,is(1)-5):min(length(tracking_data(i).donut),is(1)+5);
    is2 = find(tracking_data(i).donut(subset));
    if(isempty(is2) || mean(tracking_data(i).donut(subset(is2(1):end)))<.5)
        categories(i,5)=1;
    end
    
    %6: the profile fit is quite bad on at least 1 frame and/or sorta bad on
    %first donut frame
    if(min(tracking_data(i).gof)<.93 || tracking_data(i).gof(is(1))<.98)
        categories(i,6)=1;
    end
    
    %7: there is high intensity pixel at the edge of an image
    edge_problems=0;
    for(j=1:length(tracking_data(i).imgs))
        [a,b] = find(tracking_data(i).imgs{j}>10000);
        list_edge = a==1 | a==11 | b==1 | b==11;
        D = squareform(pdist([a(list_edge),b(list_edge)]));
        for(k=1:length(D))
            if(length(find(D(k,:)==1))>=2)
                edge_problems = 1;
                break;
            end
        end
        if(edge_problems==1)
            break;
        end
    end
    
    if(edge_problems==1)
        categories(i,7)=1;
        continue;
    end
    
    % skip if any of the exclusion criteria were met
    if(sum(categories(i,:))>0)
        continue;
    end

    
    timgs = tracking_data(i).imgs;
    tg2fit = tracking_data(i).g2fit;
    
    tracking_data(i).imgs = {};
    tracking_data(i).g2fit = {};
    for j = subset
        tracking_data(i).imgs{end+1} = timgs{j};
        tracking_data(i).g2fit{end+1} = tg2fit{j};
    end

    tracking_data(i).gof = tracking_data(i).gof(subset);
    
    
    tracking_data(i).xpos = tracking_data(i).xpos(subset);
    tracking_data(i).ypos = tracking_data(i).ypos(subset);
    tracking_data(i).frame = tracking_data(i).frame(subset);
    tracking_data(i).donut = tracking_data(i).donut(subset);
    
    list(end+1) = i;
end
% 
tracking_data = tracking_data(list);