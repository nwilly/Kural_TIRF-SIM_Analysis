%{

Creates a movie in which every trace is marked with a black pixel.  This is
used to visually confirm that the tracking worked as intended.

inputs:
    movie_path: .tif file to be marked
    fxyc_path: .mat file containing the traces

outputs:
    saved a movie in the same location as movie_path but appended with
    _traced e.g. c:/my_folder/movie.tif -> c:/my_folder/movie_traced.tif

Author: Nathan Willy (willy.2@osu.edu)

%}

function [] = trace_movie(movie_path,fxyc_path)

% get number of frames in the movie
info = imfinfo([movie_path]);
[path,file,ext] = fileparts(info(1).Filename);
frames = length(info);

% load trace data
traces = load([fxyc_path]);
traces=traces.nsta;

% loop through each frame
for(t=1:frames)

    % preallocate space for trace data of traces which exist at time t
    tr_x = zeros(1000000,1);% x position
    tr_y = zeros(1000000,1);% y position
    new_trace = zeros(1000000,1);% mark the beginning of a trace
    index=1;
    
    % loop through the traces and find the traces which exists during
    % current time t
    for(i=1:length(traces))

        try
        if(traces(i).frame(1)<=t && traces(i).frame(end)>=t)
            ind = find(traces(i).frame==t);
            tr_x(index:index+ind-1) = traces(i).xpos(1:ind);
            tr_y(index:index+ind-1) = traces(i).ypos(1:ind);
            
            new_trace(index+ind-1) = 1;
            
            index=index+ind;
        end
        catch E
           disp('error!') 
        end
    end

    
    % read current movie frame
    img = double(imread(movie_path,t));
    
    % ensure that positions are within the dimensions of the frame
    tr_x = min(size(img,2),max(1,round(tr_x)));
    tr_y = min(size(img,1),max(1,round(tr_y)));

    
    for(i=1:length(tr_x))
        % mark each trace position with a 0
        img(tr_y(i),tr_x(i)) = 0;

        % if it is not a trace's beginning, mark locations between previous
        % positions with a 0
        if(new_trace(i)~=1)
            for(j=0:.2:1)
                try
                img(round(tr_y(i)*(1-j) + tr_y(i+1)*j),...
                    round(tr_x(i)*(1-j) + tr_x(i+1)*j)) = 0.;
                catch E
                    continue
                end
            end
        end
    end

    img = uint16(img);
    
    % if first frame, create/overwrite file.  Otherwise append img.
    if(t==1)
        imwrite(img,strjoin({path,strcat(file,'_traced.tif')},'\'));
    else
        imwrite(img,strjoin({path,strcat(file,'_traced.tif')},'\'),'WriteMode','append');
    end
    
    % plot the final frame for a quick check if the tracking looks good
%     if(t==frames)
%         imagesc(img)
%     end
end    
end

