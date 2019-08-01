%{
    This will run through the entire analysis pipeline and generate figures
    similar to those which appear in the paper.  The estimated runtime is 
    about 10 minutes.

	Included with the distribution of the software should be files 
    test_movie_sim.tif and test_movie_raw.tif.  These are an excerpt from 
    a movie taken of a COS-7 cell expressing clathrin-mEmerald using a high
    NA TIRF-SIM microscope.  test_movie_sim is the SIM reconstruction from 
    the data, and test_movie_raw is the average of the 9 fluorescent images 
    used to make the SIM reconstruction.  Frame rate is 1frame/7seconds.

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
        
        g2fit: coeffiecients of the fit of the sum of two Gaussians to the
            average cross-sectional profile of the structure
            F = fittype('back + amp1*exp(-((x-x1).^2/(2*std^2))) + amp2*exp(-((x-x2).^2/(2*std^2)))');
            temp = tracking_data(i).g2fit{j};
            g2fit = cfit(F,temp(2),temp(3),temp(1),temp(6),temp(4),temp(5));
        
        gof: R2 coeficient of the 2 gaussian fit to the average
            cross-sectional profile
        
        area: calculated area in nm^2
        
        max_rmax: the structures largest peak distance of the single gaussian fit

Author: Nathan Willy (willy.2@osu.edu)

%}

% analysis parameters
folder = '.'; % folder where the data is kept
movie = {'test_movie_sim.tif'}; % cell array of SIM movie(s) to be analyzed
raw_movie = {'test_movie_raw.tif'}; % cell array of raw movie(s) to be analyzed
frame_rate = [7]; % array of movie frame rates (seconds)
pixel = 33.1; % pixel size in nm


% preprocess the SIM movie by filling any hollow structures.  This allows
% our tracking software TraCKer to work much better.  This will save the
% resulting movie as *_closed.tif.
disp('Closing Movies')
for i = 1:length(movie)
    closer_func(folder,movie{i});
end

% Run particle tracking program on the closed movies.  294 is a hard coded
% intensity threshold which we found worked well for the demo movie.  That
% value should be modified/removed before running on other data.
disp('TraCKing')
mat_files = {};
for i = 1:length(movie)
    mat_files{i} = TraCKer(fullfile(folder,[movie{i}(1:end-4),'_closed.tif']),9,294);
    
    % Convert the output of TraCKer to a more convenient data structure.
    trace_to_struct(mat_files{i},frame_rate(i),100);
    mat_files{i} = [mat_files{i}(1:end-4),'_struct.mat'];
    
    % Create a movie in which each particle's trace is marked with a black
    % pixel.  This step is ooptional and only serves to confirm that TraCKer
    % work as intended.
    load(mat_files{i});
    trace_movie(fullfile(folder,[movie{i}(1:end-4),'_closed.tif']),mat_files{i});
end





% Get cropped images of each structure found in the traces.
disp('creating sub images')
tracking_data = load_sub_imgs(mat_files, movie, raw_movie, frame_rate);

% The code aggressively removes 
% any found to be in close proximity to another in order to guarantee that 
% results are not misleading. In order that a small test movie can have 
% significant number of traces be analyzed, for demonstration purposed I have
% skipped this step.  Simply uncomment the line to re-include it.
% disp('remove proximal traces')
% tracking_data = stitch_tracks(tracking_data);

% Determine whether each cropped image is of a donut (a curved clathrin
% coat) or not.
disp('determining donuts')
tracking_data = is_donut(tracking_data);

% Perform further analysis of the sub images calculating the
% cross-sectional fits
disp('running donut analysis')
tracking_data = donut_analysis(tracking_data);

% remove structures which are unanalyzable
disp('selecting good structures')
tracking_data = select_good_structures(tracking_data,frame_rate);

% calculate areas of the remaining structures
disp('calculating dcanny')
tracking_data = calc_dcanny(tracking_data, pixel);

% Calculate structure averages
disp('calculating average structures')
[results, splits] = calculate_structure_averages(tracking_data,frame_rate);

% Plot the average images and profiles
% This will create a warning about setting coefficient values.  This is
% normal.
graph_average_structures(results, splits, frame_rate, 1, pixel,'COS7')



