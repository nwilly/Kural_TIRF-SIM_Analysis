%{
Converts the data format that comes out of TraCKer to our struct data
structure.  

inputs:
    file: .mat file to be converted (output from TraCKer)
    frame_rate: time between frames (seconds)
    bkgrd: average background intensity

Author: Josh Ferguson
Contact: Nathan Willy (willy.2@osu.edu)
%}

function trace_to_struct(file,frame_rate,bkgrd)

    load(file)
    nsta = struct('frame',[],'xpos',[],'ypos',[],'int',[],'class',3);
    for(i=1:size(TraceX,1))
        list = find(TraceX(i,:));
        nsta(i).frame = list;
        nsta(i).xpos = TraceX(i,list);
        nsta(i).ypos = TraceY(i,list);
        nsta(i).int = TraceINT(i,list)';
        nsta(i).class=3;
        nsta(i).lt = length(list);
    end
    
    nsta = slope_finding(nsta,frame_rate,bkgrd);
    save(strcat(file(1:end-4),'_struct.mat'),'nsta');

end

