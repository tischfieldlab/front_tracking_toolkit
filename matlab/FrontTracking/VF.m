function [vx,vy,x,y,t,fi,th]=VF(vfronts,framerange,ROI,noisy)
% Usage: [vx,vy,x,y,t,fi,th]=VF(vfronts,[framerange],[ROI],[noisy])
% Working from the tracked fronts in "vfronts", VF plots the fronts and
% their velocities in frames "framerange" in the area "ROI". Velocity
% values are returned in vectors "vx" and "vy"; corresponding positions,
% times, and track indices are returned in "x", "y", and "t", respectively.
% Front indices are returned in "fi", and thicknesses are returned in "th".
% If noisy==0, no plot is produced. The input "vfronts" must be a structure
% of the form produced by TheFrontTracker.m. Specify "framerange" as a
% two-element vector: [starttime endtime]. Specify the region of interest,
% "ROI" as a four element vector, [xmin ymin width height]. 

% Modified March 21 2017 by Thomas Nevins from Velocities.m.

framerangedefault = [-inf inf]; % all frames
ROIdefault=[0 0 inf inf]; % all space
noisydefault=1;

if nargin<1
    error(['Usage: [vx,vy,x,y,t,fi,th] = ' mfilename ...
        '(vfronts,[framerange],[ROI],[noisy])'])
end
if ~exist('framerange','var') || isempty(framerange)
    framerange=framerangedefault;
elseif numel(framerange)==1
    framerange=framerange*[1 1];
end
if ~exist('ROI', 'var') || isempty(ROI)
    ROI=ROIdefault;
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisydefault;
end
fn=fieldnames(vfronts);
if ~any(strcmp(fn,'X')) || ~any(strcmp(fn,'Y')) || ~any(strcmp(fn,'T')) ...
        || ~any(strcmp(fn,'Vx')) || ~any(strcmp(fn,'Vy')) || isempty(vfronts)
    error('Sorry, the input does not appear to contain tracks.')
end

if isempty(vertcat(vfronts.Thickness))
    hasthick=0;
else
    hasthick=1;
end

Nall=sum([vfronts.len]);
vx=NaN(Nall,1); % initialize arrays
vy=NaN(Nall,1);
x=NaN(Nall,1);
y=NaN(Nall,1);
t=NaN(Nall,1);
fi=NaN(Nall,1);
if hasthick
    th=NaN(Nall,1);
else
    th=[];
end
pos=1; % marker to current position in arrays
for ii=1:numel(vfronts) % assemble arrays from vfronts structure
    ind = (vfronts(ii).T>=min(framerange)) & (vfronts(ii).T<=max(framerange));
    if ind
        Np=vfronts(ii).len;
        vx(pos:pos+Np-1)=vfronts(ii).Vx;
        vy(pos:pos+Np-1)=vfronts(ii).Vy;
        x(pos:pos+Np-1)=vfronts(ii).X;
        y(pos:pos+Np-1)=vfronts(ii).Y;
        if hasthick
            th(pos:pos+Np-1)=vfronts(ii).Thickness;
        end
        t(pos:pos+Np-1)=vfronts(ii).T;
        fi(pos:pos+Np-1)=ii;
        pos=pos+Np;
    end
end
ind=~isnan(vx); % remove unwanted frames
vx=vx(ind);
vy=vy(ind);
x=x(ind);
y=y(ind);
t=t(ind);
fi=fi(ind);
if hasthick
    th=th(ind);
end
[t,ind]=sort(t(:)); % sort by time
vx=vx(ind);
vy=vy(ind);
x=x(ind);
y=y(ind);
fi=fi(ind);
if hasthick
    th=th(ind);
end
ind= ( (x>=ROI(1)) & (x<=ROI(1)+ROI(3)) & ...
        (y>=ROI(2)) & (y<=ROI(2)+ROI(4)) ); % exclude points near edges
vx=vx(ind);
vy=vy(ind);
x=x(ind);
y=y(ind);
t=t(ind);
fi=fi(ind);
if hasthick
    th=th(ind);
end
% If requested, plot the particles
if noisy
    figure;
    quiver(x,y,vx,vy,'r');
    set(gca,'dataaspectratio',[1 1 1],'nextplot','add','ydir','reverse');
    axis tight
    [~,ends,~]=unique(fi,'legacy');
    begins=circshift(ends,1)+1;
    begins(1)=1;
    Nfr=numel(ends);
    for ii=1:Nfr
        ind=begins(ii):ends(ii);
        plot(x(ind),y(ind),'k');
    end
    if framerange(1)==framerange(2)
        title([num2str(numel(x)) ' fronts in frame ' ...
            num2str(framerange(1))]);
    else
        if isinf(framerange(1))
            framerange(1)=t(1);
        end
        if isinf(framerange(2))
            framerange(2)=t(end);
        end
        title([num2str(numel(x)) ' fronts in frames ' ...
            num2str(min(framerange)) ' to ' num2str(max(framerange)) ]);
    end
    if exist('region','var') && ~isempty(region)
        xlim(ROI(1)+[0 ROI(3)]);
        ylim(ROI(2)+[0 ROI(4)]);
    end
end
end


% Copyright 2017 Thomas Nevins and Douglas H. Kelley
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%      http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
