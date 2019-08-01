%{

Goes through each trace in each mat file and gets an 11x11 pixel cropped
image centered on the structure in both the SIM and raw movies.

inputs:
    mat_files: cell array of .mat files containing trace information
    sim_mov: cell array of the reconstructed SIM movie(s)
    raw_movs: cell array of the raw SIM movie(s)
    frame_rates: list of movie frame rates (seconds)

output:
    tracking_data: structure containing the combined trace information of
        all the mat files in mat_files as well as the cropped sub images
        (reconstructed and raw).

    tracking_data has the following fields:
        
        frame: the array of frames in which the structure appears
        
        xpos: the x pixel coordinates in which the structure appears
        
        ypos: the y pixel coordinates in which the structure appears
        
        int: array of pixel intensities of the structure
        
        lt: lifetime, the number of frames over which the structure is tracked
        
        donut: boolean whether the structure in the frame appears as a donut
        
        movie: index of the movie the trace came from
        
        imgs: cell array of cropped images centered on the TIRF-SIM
            structure image
        
        imgs_raw: cell array of cropped images centered on the raw
            structure image

Author: Nathan Willy (willy.2@osu.edu)

%}

function tracking_data = load_sub_imgs(mat_files, sim_movs, raw_movs, frame_rates)

frames=zeros(length(frame_rates),1);

% combine tracking data from mat files into a single structure
disp('loading tracking data')
tracking_data = [];
for(i=1:length(mat_files))
    load(mat_files{i});
    first_frame = imread(sim_movs{i},1);

    frames(i) = length(imfinfo(sim_movs{i}));
        
    nsta = rmfield(nsta, 'sl');
    nsta = rmfield(nsta, 'class');
    if(isfield(nsta,'donut')==0)
%         nsta = rmfield(nsta, 'donut');
        nsta(end).donut= [];
    end
    
    if(isfield(nsta,'r2'))
        nsta = rmfield(nsta, 'r2');
    end
    
    for(j=1:length(nsta))
        nsta(j).movie = uint8(i);
        
        nsta(j).frame = uint16(nsta(j).frame);
        nsta(j).lt = uint16(nsta(j).lt);
        
        nsta(j).xpos = single(nsta(j).xpos);
        nsta(j).ypos = single(nsta(j).ypos);
        nsta(j).int = single(nsta(j).int);
        
        if(min(nsta(j).xpos)<10 || min(nsta(j).ypos)<10)
            nsta(j).lt = 0;
        elseif(max(nsta(j).xpos)>size(first_frame,2)-10 || ...
                max(nsta(j).ypos)>size(first_frame,1)-10)
            nsta(j).lt = 0;
        end
        
    end
    
    if(i==1)
        tracking_data = nsta;
    else
        tracking_data(end+1:end+length(nsta)) = nsta;
    end
    clear nsta
end
tracking_data([tracking_data.lt]==0) = [];

% go through the traces and collect cropped sub-images
disp('loading sub images')
% loop through the different movies
for(file=1:length(sim_movs))
    disp(['file: ',num2str(file)])
    first_frame = imread(sim_movs{file},1);
    
    % load entire SIM movie
    movie = zeros(size(first_frame,1),size(first_frame,2),frames(file));
    for(t=1:frames(file))
        movie(:,:,t) = imread(sim_movs{file},t);
    end
    
    % load entire raw movie
    movie_raw = zeros(size(first_frame,1),size(first_frame,2),frames(file));
    for(t=1:frames(file))
        movie_raw(:,:,t) = imread(raw_movs{file},t);
    end
    
    % loop through th tracks for this movie
    for(i=1:length(tracking_data))
        if(tracking_data(i).movie~=file)
            continue;
        end
        
        % crop sub-images
        tracking_data(i).imgs={};
        for(j=1:length(tracking_data(i).frame))
        
            x = tracking_data(i).xpos(j);
            y = tracking_data(i).ypos(j);

            sub = movie(round(y-5:y+5),round(x-5:x+5),tracking_data(i).frame(j));
            sub_raw = movie_raw(round(y-5:y+5),round(x-5:x+5),tracking_data(i).frame(j));
            tracking_data(i).imgs{j} = single(sub);
            tracking_data(i).imgs_raw{j} = single(sub_raw);

        end
    end
    clear movie;
    clear movie_raw;
end