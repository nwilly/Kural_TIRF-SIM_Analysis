function fxyc_struct = slope_finding(fxyc_struct,frame_rate,bkgrd,varargin)
%SLOPE_FINDING Create slope field in trace structure.
%   Inputs:
%       fxyc_struct: Structure output by fxyc_to_struct function.
%       frame_rate: Frame rate of movie in seconds. Important for
%           generating as correct as possible 12s windows.
%       bkgrd: Background. This is usually the same as the input from
%           comb_run. It is subtracted from all intensity traces, since it
%           is judged to be the lowest possible recordable signal.
%   A fourth argument can be added if you want to monitor the process, but
%   the code is very fast and rarely warrants monitoring.
%
% Also note that non analyzed slopes are set to zero. Due to the extreme
% unlikelihood of getting a real slope of exactly zero, all zeroes are
% excluded during analysis. If you don't like this, feel free to alter the
% code and use NaNs. For some reason I can't remember, this plan wasn't
% effective for us.
%
% Josh Ferguson, Kural Group, Ohio State University, ferguson.621@osu.edu
%
if nargin == 4
    monitor = varargin{1};
else
    monitor = false;
end
ints = {fxyc_struct.int}; %reduce text in code.
prange = 12/frame_rate; %develop a frame range of 12 seconds
%we want at least points to fit the line, so if the frame rate is greater
%than 4 seconds, the frame range will increase beyond 12 seconds.
if prange<3, prange=3; end 
%originally designed to analyze data using no future knowledge, you can
%adjust the percent of the frame range that looks forward in time.
%adjusting this should change the shape of histograms, so maintain
%consistency while comparing analyses.
forwardp = .25;
front = ceil(forwardp*(prange-1)); %frames forward in time
rear = floor((1-forwardp)*(prange-1)); %frame backward in time
if monitor, fprintf('Percent Complete: %3i%%',0); end
for i = 1:length(ints)
    if isempty(ints{i})
        fxyc_struct(i).sl = [];
        fxyc_struct(i).r2 = [];
        continue;
    end
    int = ints{i}; %reduce text in code.
    lint = length(int);
    intdif = zeros(lint,1); 
    r2 = zeros(lint,1);
    if lint<=prange %if the trace is too short, skip it.
        fxyc_struct(i).sl = single(intdif);
        fxyc_struct(i).r2 = single(intdif);
        continue;
    end
    %this is where we perform manual least-squares fitting
    for j = (rear+1):(lint-front)
        sub = (j-rear):(j+front);
        curmax = max(int)-bkgrd;
        tmp = (int(sub)-bkgrd)/curmax;
        tmpx = sub*frame_rate;
        tmpy = tmp';
        numer = length(tmpx)*sum(tmpx.*tmpy)-sum(tmpx)*sum(tmpy);
        denom = length(tmpx)*sum(tmpx.^2)-sum(tmpx)^2;
        numer2 = sum(tmpy)*sum(tmpx.^2)-sum(tmpx)*sum(tmpx.*tmpy);
        intdif(j) = numer/denom;
        intercept = numer2/denom;

        R = corrcoef((intdif(j)*sub + intercept),tmpy);
        r2(j) = R(1,2);
    end
    fxyc_struct(i).sl = single(intdif); %add data to existing structure.
    fxyc_struct(i).r2 = single(r2);
    if monitor, fprintf('\b\b\b\b%3i%%',ceil(100*i/length(ints))); end
end
if monitor, fprintf('\b\b\b\b%3i%%\n',100); end
end