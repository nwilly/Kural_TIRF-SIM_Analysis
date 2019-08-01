%{
    Removes any traces which are ever within 5 pixels of another trace.

    inputs:
        tracking_data: structure with trace information.

    outputs:
        tracking_data: trace data structure with proximal traces removed

Author: Nathan Willy (willy.2@osu.edu)
%}

function tracking_data = remove_proximal_tracks(tracking_data)

for m=1:max([tracking_data.movie])
    
    disp(num2str(m))
    
    list = find([tracking_data.movie]==m);
    X = zeros(length(vertcat([tracking_data(list).xpos])),1);
    Y = zeros(size(X));
    T = zeros(size(X));
    index = zeros(size(X));
    
    ind=1.;
    for(i=1:length(list))
        X(ind:ind+length(tracking_data(list(i)).xpos)-1) = tracking_data(list(i)).xpos;
        Y(ind:ind+length(tracking_data(list(i)).xpos)-1) = tracking_data(list(i)).ypos;
        T(ind:ind+length(tracking_data(list(i)).xpos)-1) = tracking_data(list(i)).frame;
        index(ind:ind+length(tracking_data(list(i)).xpos)-1) = ones(length(tracking_data(list(i)).xpos),1)*i;
        
        ind = ind + length(tracking_data(list(i)).xpos);
    end
    
    conflicted = zeros(length(list),1,'logical');
    for(t=1:max(T))
        list2 = find(abs(T-t)<=1);
        D = squareform(pdist([X(list2),Y(list2)]));
        [a,b] = find(D<5);
        
        conflicted(a)=1;
        conflicted(b)=1;
    end
            
    tracking_data(list(conflicted)) = [];
end
end