%{

Split the structures by movie and then by size in order to calculuate
statistics about average structures.

Each movie is analyzed serparately.  Structures are then split into three
approximately equal groups based on size (max_rmax).  Area data is then
combined based on these groupings

inputs:
    tracking_data, data structure containing image and area data for
        structures
    frame_rate, list of frame rates for the movies (seconds)

outputs:
    results: data structures containing info about the structure groupings
        results will be 3 x m long where m is the number of movies
        represented in tracking data.  The 3 comes from the structures of
        each movie being split into small medium and large structures.
        e.g. results(4) would give data on medium structures from the
        second movie.

        Each group of structures has 11 frames of data centered such that
        frame 6 is the first frame where the structure had resolvable
        curvature. eg frame 5 gives data on the average of structures 1
        frame before the first resolvable curvature.  Some structures may
        not have 5 full frames before and after the first curved frame and
        so the counts field enumerates how many structures are included in
        the average.

        imgs: average sim images of all structures in the grouping
        imgs2: average of sim images scaled 2x
        imgs10: average of sim images scaled 10x
        
        time: time relative to the first occurance of resolvable curvature
            (seconds)
        count: the number of structures included in the avregae
        area: the calculated average area nm^2
        dpeak: the average distance between peaks in the 2 Gaussian fit of
            the average cross-sectional profile nm
        area_list: the list of all area values in the group (in case you
            want to do further analysis) nm^2
        dpeak_list: the list of all dpeak values in the group (in case you
            want to do further analysis) nm
        data_index: list of indices to allow back reference to
            tracking_data. It gives the index of the structure in
            tracking_data and the frame number in tracking_data that the
            average frame corresponds to.

    splits: cell array, gives the thresholds for small medium and large
        structures as calculated for each movie

Author: Nathan Willy (willy.2@osu.edu)
%}

function [results,splits] = calculate_structure_averages(tracking_data, frame_rate)

list=[];
for(i=1:length(tracking_data))
    first = find(tracking_data(i).donut);
    if(isempty(first))
        list(end+1) = i;
        continue
    end
    tracking_data(i).time = ([1:length(tracking_data(i).donut)]-first(1))*frame_rate(tracking_data(i).movie);
    tracking_data(i).first_time = tracking_data(i).time(1);
end

tracking_data(list)=[];
    
g_list = {};
splits={};
g_ind=1;
for(m=1:max([tracking_data.movie]))
    splits{m}=[60,0,0,120]+100;
    temp_list = find([tracking_data.movie]==m &  [tracking_data.max_rmax]< splits{m}(end) & [tracking_data.max_rmax]>=splits{m}(1) & [tracking_data.first_time]<0);
    splits{m}(2) = quantile([tracking_data(temp_list).max_rmax],.333);
    splits{m}(3) = quantile([tracking_data(temp_list).max_rmax],.666);
    
    splits{m}(2) = round(splits{m}(2)*1)/1;
    splits{m}(3) = round(splits{m}(3)*1)/1;


    for(s=1:3)
        g_list{g_ind} = temp_list(find([tracking_data(temp_list).max_rmax]< splits{m}(s+1) & [tracking_data(temp_list).max_rmax]>=splits{m}(s)));        
        g_ind = g_ind+1;
    end
end

results = struct('imgs',{},'imgs2',{},'imgs10',{},'time',[],'count',[],'area',[],'dpeak',[],'area_list',[],'dpeak_list',[],'data_index',[]);

for(i=1:length(g_list))
    
    disp(num2str(i));
    
    results(i).imgs = cell(11,1);
    results(i).imgs2 = cell(11,1);
    results(i).imgs10 = cell(11,1);
    
    results(i).count = zeros(11,1);
    results(i).area = zeros(11,1);
    results(i).dpeak = zeros(11,1);
    
    results(i).area_list = cell(11,1);
    results(i).dpeak_list = cell(11,1);

    results(i).data_index = cell(11,1);

    
    if(isempty(g_list{i}))
        continue;
    end
    
    curr_fr = frame_rate(tracking_data(g_list{i}(1)).movie);
    
    results(i).time = [-5:5]*curr_fr;
    
    for(j=1:11)
        results(i).imgs{j} = zeros(11);
        results(i).imgs2{j} = zeros(11*2);
        results(i).imgs10{j} = zeros(11*10);
    end
    
    for(j=1:length(g_list{i}))
        
        ind = g_list{i}(j);

        
        for(t=1:length(tracking_data(ind).time))
            
            index = tracking_data(ind).time(t)/curr_fr + 6;
            if(index<1 || index>11)
                continue;
            end
            
            i2 = imresize(tracking_data(ind).imgs{t},2,'bicubic');
            i10 = imresize(tracking_data(ind).imgs{t},10,'bicubic');
            
            results(i).imgs{index} = results(i).imgs{index} + tracking_data(ind).imgs{t};
            results(i).imgs2{index} = results(i).imgs2{index} + i2;
            results(i).imgs10{index} = results(i).imgs10{index} + i10;
            
            results(i).count(index) = results(i).count(index) + 1;
            results(i).area(index) = results(i).area(index) + tracking_data(ind).area(t);
            
            temp = tracking_data(ind).g2fit{t};
            p2p = abs(temp(5)-temp(4));
            
            results(i).dpeak(index) = results(i).dpeak(index) + p2p;
            
            results(i).area_list{index}(end+1) = tracking_data(ind).area(t);
            results(i).dpeak_list{index}(end+1) = p2p;
            results(i).data_index{index}(end+1,:) = [ind,t];
            

        end
    end
    
    for(j=1:11)
        results(i).imgs{j} = results(i).imgs{j}/results(i).count(j);
        results(i).imgs2{j} = results(i).imgs2{j}/results(i).count(j);
        results(i).imgs10{j} = results(i).imgs10{j}/results(i).count(j);

        results(i).area(j) = results(i).area(j)/results(i).count(j);
        results(i).dpeak(j) = results(i).dpeak(j)/results(i).count(j);

    end
end
