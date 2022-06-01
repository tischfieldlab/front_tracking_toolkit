function [vx,vy,x,y,t,tr]=Velocities(vtracks,framerange,noisy)
% Usage: [u,v,x,y,t,tr]=Velocities(vtracks,[framerange],[noisy])
% Working from the velocity tracks in "vtracks", Velocities plots the
% velocity of each particle in frames "framerange". Velocity values are 
% returned in vectors "u" and "v"; corresponding positions, times, and 
% track indices are returned in "x", "y", and "t". If noisy==0, no plot is 
% produced. The input must be a structure of the form produced by 
% PredictiveTracker.m. Specify "framerange" as a two-element vector: 
% [starttime endtime]. This file can be downloaded from 
% http://leviathan.eng.yale.edu/software.

% Written 13 April 2011 by Douglas H. Kelley.
% Fixed frame range in title when using default 1 August 2011. 
% Included tr output 26 March 2012. 
% Pre-allocated for speed (making a real difference) 27 March 2012. 

framerangedefault = [-inf inf]; % all frames
noisydefault=1;

if nargin<1
    error(['Usage: [u,v,x,y,t,tr] = ' mfilename ...
        '(vtracks,[framerange],[noisy])'])
end
if ~exist('framerange','var') || isempty(framerange)
    framerange=framerangedefault;
elseif numel(framerange)==1
    framerange=framerange*[1 1];
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisydefault;
end
fn=fieldnames(vtracks);
if ~any(strcmp(fn,'X')) || ~any(strcmp(fn,'Y')) || ~any(strcmp(fn,'T')) ...
        || ~any(strcmp(fn,'U')) || ~any(strcmp(fn,'V')) || isempty(vtracks)
    error('Sorry, the input does not appear to contain tracks.')
end
Nall=sum([vtracks.len]);
vx=NaN(Nall,1); % initialize arrays
vy=NaN(Nall,1);
x=NaN(Nall,1);
y=NaN(Nall,1);
t=NaN(Nall,1);
tr=NaN(Nall,1);
pos=1; % marker to current position in arrays
for ii=1:numel(vtracks) % assemble arrays from vtracks structure
    ind = (vtracks(ii).T>=min(framerange)) & (vtracks(ii).T<=max(framerange));
    if any(ind)
        Np=sum(ind);
        vx(pos:pos+Np-1)=vtracks(ii).U(ind(:));
        vy(pos:pos+Np-1)=vtracks(ii).V(ind(:));
        x(pos:pos+Np-1)=vtracks(ii).X(ind(:));
        y(pos:pos+Np-1)=vtracks(ii).Y(ind(:));
        t(pos:pos+Np-1)=vtracks(ii).T(ind(:));
        tr(pos:pos+Np-1)=ii;
        pos=pos+Np;
    end
end
ind=~isnan(vx); % remove unwanted frames
vx=vx(ind);
vy=vy(ind);
x=x(ind);
y=y(ind);
t=t(ind);
tr=tr(ind);
[t,ind]=sort(t(:)); % sort by time
vx=vx(ind);
vy=vy(ind);
x=x(ind);
y=y(ind);
tr=tr(ind);

% If requested, plot the particles
if noisy
    figure;
    quiver(x,y,vx,vy);
    set(gca,'dataaspectratio',[1 1 1]);
    axis tight
    if framerange(1)==framerange(2)
        title([num2str(numel(x)) ' particles in frame ' ...
            num2str(framerange(1))]);
    else
        if isinf(framerange(1))
            framerange(1)=t(1);
        end
        if isinf(framerange(2))
            framerange(2)=t(end);
        end
        title([num2str(numel(x)) ' particles in frames ' ...
            num2str(min(framerange)) ' to ' num2str(max(framerange)) ]);
    end
end

%Copyright 2011 Douglas H. Kelley, and Nicholas Ouellette
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

