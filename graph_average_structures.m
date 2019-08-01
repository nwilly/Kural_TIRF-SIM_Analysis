%{

Plot the average structures and their cross sectional profiles of a movie.
The sub graphs are saved to the working directory as .png images

Note: This script produces a warning:
Warning: Setting coefficient values clears confidence bounds information.
This is expected and normal.

inputs:
    results: data structure containing data on the average structures
    splites: cell array of size thresholds used for separating structures
        into small medium and large.
    fr: list of movie frame rates (seconds)
    movie_index: which movie should be plotted
    pix: pixels size (nm)
    cell_type: str, label used for graphs

outputs:
    saves the sub images to the working directory.  The image naming is as
    follows:
        i-(1/2/3): small/medium/large
        j-(1->11): time in frames where frame 6 is the first donut frame
        num-x: number of structures represented in the average
        low-x: the lowest max_rmax threshold (nm)
        high-x: the highest max_rmax threshold (nm)
        _profile: the image is of the average cross sectional profile

        _areas: plots the progression of the average area for small (red)
            medium (green) and large (blue) structures (nm^2)

        _dpeaks: plots the progression of the average dpeak value for 
            small (red) medium (green) and large (blue) structures (nm)

Author: Nathan Willy (willy.2@osu.edu)

%}


function [] = graph_average_structures(results, splits, fr, movie_index, pix, cell_type)

close all
figure
areas = zeros(3,11);
d_peak=zeros(3,11);


for(j=(1:3)+(movie_index-1)*3)
    
    first_sub=1;
    
    for(i=1:11)
            
        if(sum(isnan(results(j).imgs{i}(:)))>0)
            areas(mod(j-1,3)+1,i) = nan;
            d_peak(mod(j-1,3)+1,i) = nan;
            continue;
        end
        
        img2 = results(j).imgs10{i};
        mask = edge(img2,'canny',[.5 .7],10)>0;
        mask = imdilate(mask,ones(8));
        mask = bwmorph(mask,'thin',12);
        mask = imfill(mask,'holes')>0;
        area = sum(mask(:))*pix*pix/100;
        
        results(j).mask{i} = mask;
        areas(mod(j-1,3)+1,i) = area;
        
        img_analysis = donut_img_analysis(results(j).imgs{i},pix,[],[]);
        dpeak = img_analysis(1).fit.x2-img_analysis(1).fit.x1;
        
        
        d_peak(mod(j-1,3)+1,i) = dpeak;
                
        subplot(3,11,mod((j-1),3)*11 + i)

        imagesc(results(j).imgs10{i})
        draw_mask_outline(mask,1,1);

        axis equal
        title(strcat('t=',num2str((i-6)*fr(movie_index)),'s, #=',num2str(results(j).count(i))))%,' area=',num2str(floor(area/100)*100),'nm^2'))
        
        if(first_sub==1)
            if(mod(j-1,3)==0)
                ylabel([num2str(splits{movie_index}(1)),'<rmax_{max}<',num2str(splits{movie_index}(2)),'nm'])
            elseif(mod(j-1,3)==1)
                ylabel([num2str(splits{movie_index}(2)),'<rmax_{max}<',num2str(splits{movie_index}(3)),'nm'])
            else
                ylabel([num2str(splits{movie_index}(3)),'<rmax_{max}<',num2str(splits{movie_index}(4)),'nm'])
            end
        end
        
        first_sub=0;

            
        figure('innerposition',[512,512,512,512],'visible','off');
        imagesc(results(j).imgs10{i})
        draw_mask_outline(mask,1,1);
        axis equal
        
        saveas(gca,['./i-',num2str(i),'_j-',num2str(j),'_num-',num2str(results(j).count(i)),...
            '_low-',num2str(splits{movie_index}(1)),'_high-',num2str(splits{movie_index}(2)),'.png']);

        close

        temp = imread(['./i-',num2str(i),'_j-',num2str(j),'_num-',num2str(results(j).count(i)),...
            '_low-',num2str(splits{movie_index}(1)),'_high-',num2str(splits{movie_index}(2)),'.png']);
        temp2 = temp(78:695,106:723,:);

        imwrite(temp2,['./i-',num2str(i),'_j-',num2str(j),'_num-',num2str(results(j).count(i)),...
            '_low-',num2str(splits{movie_index}(1)),'_high-',num2str(splits{movie_index}(2)),'.png']);

        figure('innerposition',[512,512,512,512],'visible','off');

        plot(img_analysis.pos,img_analysis.profile/10^4,'b','linewidth',5)
        hold on
        f1 = img_analysis.fit;
        f1.back = f1.back/10^4;
        f2 = f1;
        f1.amp2=0;
        f1.amp1=f1.amp1/10^4;
        f2.amp1=0;
        f2.amp2=f2.amp2/10^4;

        x_vals = img_analysis.pos(1):.1:img_analysis.pos(end);
        
        plot(x_vals,f1(x_vals),'color','r','linewidth',5)
        plot(x_vals,f2(x_vals),'color','r','linewidth',5)
        plot(img_analysis.pos,img_analysis.profile/10^4,'b','linewidth',5)

        legend off
        
        xlabel('Radial Position (nm)')
        ylabel('Intensity (AU)')
        
        set(gca,'fontsize',30)
        
        saveas(gca,['./i-',num2str(i),'_j-',num2str(j),'_num-',num2str(results(j).count(i)),...
            '_low-',num2str(splits{movie_index}(1)),'_high-',num2str(splits{movie_index}(2)),'_profile.png']);
        
        close
    end
end

a = annotation('textbox',[.5,.9,.1,.1],'String',cell_type,'FitBoxToText','on');
a.FontSize=20;


figure('units','pixel','innerposition',[0 0 1542*4/3/2 1542/2],'visible','off')

plot((-5:5)*fr(movie_index),areas(1,:),'color',[.7 .05 .05],'linewidth',5)
hold on
plot((-5:5)*fr(movie_index),areas(2,:),'color',[.05 .7 .05],'linewidth',5)
plot((-5:5)*fr(movie_index),areas(3,:),'color',[.05 .05 .7],'linewidth',5)
title([cell_type,' Area'])

t = ylabel('Projected Area (nm^2)');

xlabel('Time (s)')
set(gca,'fontsize',50/2)

saveas(gca,['./j-',num2str(j),...
    '_low-',num2str(splits{movie_index}(1)),'_high-',num2str(splits{movie_index}(2)),'_areas.png']);

figure('units','pixel','innerposition',[0 0 1542*4/3/2 1542/2])
plot((-5:5)*fr(movie_index),d_peak(1,:),'color',[.7 .05 .05],'linewidth',5)
hold on
plot((-5:5)*fr(movie_index),d_peak(2,:),'color',[.05 .7 .05],'linewidth',5)
plot((-5:5)*fr(movie_index),d_peak(3,:),'color',[.05 .05 .7],'linewidth',5)
title([cell_type,' delta peak'])

t = ylabel('Distance peak-to-peak (nm)');


xlabel('Time (s)')                       
set(gca,'fontsize',50/2)

saveas(gca,['./j-',num2str(j),...
    '_low-',num2str(splits{movie_index}(1)),'_high-',num2str(splits{movie_index}(2)),'_dpeaks.png']);

