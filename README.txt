This software package documents the analysis we used to chacterize clathrin coated structures imaged using TIRF-SIM.  The main goal of the analysis is to determine when during the maturation of the clathrin coat curvature can be resolved.  By showing that the maximum projected area of the structure occurs at or after the first image with resolvable curvature, we show that the model in which clathrin coats first develop flat and then curve is atypical.

Author: Nathan Willy
contact: willy.2@osu.edu

System Requirements: 
	None

Software Requirements: 
	MATLAB.  I believe it will run on anything as early as MATLAB 2008a without modification. Software was written and tested in MATLAB 2017b.
	https://www.mathworks.com/help/install/ug/install-mathworks-software.html

Software Demonstrations:

1)
	Included with the distribution of the software should be files test_movie_sim.tif and test_movie_raw.tif.  These are an excerpt from a movie taken of a COS-7 cell 		expressing clathrin-mEmerald using a high NA TIRF-SIM microscope.  test_movie_sim is the SIM reconstruction from the data, and test_movie_raw is the average of the 9 		fluorescent images used to make the SIM reconstruction.  Frame rate is 1frame/7seconds.

	To see the code in action, simply run the script donut_analysis_demo_script.m.  This will run through the entire analysis pipeline and generate figures similar to 		those which appear in the paper.  The estimated runtime is about 15 minutes.  This script can be easily modified to run on other data.  See the script's comments for 		details.

2)
	Included with distribution of the software should be the file ccs_coords.mat, which contains sample data for use in our fluorophore simulation. The sample data is 		the hemisphere of a 60 vertex polyhedron with icosahedral symmetry.  The radius is made to be 80nm.  This models a half completed clathrin coat.

	Running the script fluor_sim_demo_script.m will load the sample coordinates and create and display an image using the fluorescence simulation.
		

