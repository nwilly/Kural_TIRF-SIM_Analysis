%{

Draw lines around the perimeter of a structure mask.

inputs:
    mask: boolean matrix showing the pixels which include the structure
    d_scale: the current image scaling
    m_scale: mask scaling

draws onto the current axis

Author: Nathan Willy (willy.2@osu.edu)

%}

function [] = draw_mask_outline(mask,d_scale,m_scale)

    mask = imresize(mask,d_scale,'bilinear');
    hold on
    
    for(i=1:size(mask,1)-1)
        for(j=1:size(mask,2)-1)
            
            %right
            if(mask(i,j)~=mask(i,j+1))
                line(([j,j]+.5)*m_scale+.5*0,([i,i+1]-.5)*m_scale+.5*0,'Color','black','linewidth',3,'marker','s','markersize',.5)
            end
            
            %down
            if(mask(i,j)~=mask(i+1,j))
                line(([j,j+1]-.5)*m_scale+.5*0,([i,i]+.5)*m_scale+.5*0,'Color','black','linewidth',3,'marker','s','markersize',.5)
            end
        end
    end

end