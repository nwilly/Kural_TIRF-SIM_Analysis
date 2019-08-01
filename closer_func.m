%{
Preprocesses a TIRF-SIM movie, filling in hollow structures to make the
movie track-able by TraCKer.

inputs:
    folder: the path where the movie can be found
    movie_name: the name of the tiff file to be closed

outputs:
    saves the resultant movie into folder with the movie_name appended with
    _closed. e.g. c:/my_folder/movie.tif -> c:/my_folder/movie_closed.tif

author: Nathan Willy (willy.2@osu.edu)
%}

function [] = closer_func(folder,movie_name)

se = strel('disk',1);
k = fspecial('disk',2);

% get the number of frames in the movie
info = imfinfo(fullfile(folder,movie_name));
frames = numel(info);

% loop through the frames
for(t=1:frames)
    
    %close by average blur followed by an erosion operation
    img = double(imread(fullfile(folder,movie_name),t));
    img2 = conv2(img,k,'same');
    img3 = imerode(img2,se);
    img3 = uint16(img3);
    
    % create/overwrite file if the first frame, otherwise append frame
    if(t==1)
        imwrite(img3,fullfile(folder,[movie_name(1:end-4),'_closed.tif']));
    else
        imwrite(img3,fullfile(folder,[movie_name(1:end-4),'_closed.tif']),'writemode','append');
    end
end
end