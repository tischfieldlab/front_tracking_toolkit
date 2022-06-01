function [wx,wy,vx,vy,x,y,t,fi,th]=VF(filename,timerange,region,noisy)
% Usage: [wx,wy,vx,vy,x,y,t,fi,th]=VF(filename,[timerange],[region],[noisy])
% Working from the reaction fronts recorded in "filename", vf plots the
% fronts and their velocities at times "timerange" in the area "region".
% Velocities are returned in "wx" and "wy"; positions, "x" and "y"; times,
% "t". The velocity of the front without fluid flow (v=w-u) is also
% returned in "vx" and "vy". If noisy==0, no plot is produced. Specify
% "region" in pixels, as a vector of the form [left top width height]. The
% input must be a .gdf file, sorted framewise with x, y, t, wx, wy, vx, and
% vy in columns 2-8, respectively. For improved speed , first run
% index_gdf.m. Requires getmaxlen.m for plotting. See also vp.m and
% FrontTracker.m. 

% Written by Doug Kelley 22 March 2015, based largely on vp.m. 
% Updated 6 May 2015 for compatibility with new version of FrontTracker.m.
% Updated 12 May 2015 to plot both chemical and total velocity.
% Bug fix 14 May 2015: now dropping front indices for points outside
% region.
% Upated 25 July 2015 to scale chemical and total velocity quivers by same
% factor.

magicnum = 82991;
indexsuffix = '_ind'; % for naming files output by index_gdf.m 
timerangedefault = [0 inf]; % all frames
noisydefault = true;

% Format of input data:
ficol = 1;
xcol = 2;
ycol = 3;
tcol = 4;
wxcol=5; % front velocity w/ flow; w=v+u
wycol=6;
vxcol=7; % front velocity w/o flow
vycol=8;
thcol=9;

% Parse input arguments
if nargin < 1
    error(['Usage: [wx,wy,vx,vy,x,y,t,fi] = ' mfilename ...
        '(filename,[timerange],[region],[noisy])']);
end
if ~exist('timerange','var') || isempty(timerange)
    timerange=timerangedefault;
elseif numel(timerange)==1
    timerange=timerange*[1 1];
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisydefault;
end

if ~exist(filename,'file')
    error(['Sorry, cannot find ' filename '.']);
end

[path,shortname]=fileparts(filename);
indfile=fullfile(path,[shortname indexsuffix '.mat']); % specify full path to stay in same directory
if exist(indfile,'file') % with an index file, can do this the fast way.
    frameind=load(indfile);
    tmin=max([timerange(1) frameind.n(1)]);
    tmax=min([timerange(2) frameind.n(end)]);
    if tmin>tmax
        error(['Sorry, ' filename ...
            ' does not contain times between ' ...
            num2str(timerange(1)) ' and ' ...
            num2str(timerange(2)) '.']);
    end
    fidin=fopen(filename);
    header=fread(fidin,6,'int32');
    if header(1)~=magicnum
        error(['Sorry, ' filename ...
            ' does not appear to be a .gdf file.']);
    end
    ncol=header(3);
    indmin=find(frameind.n>=tmin,1);
    indmax=find(frameind.n>tmax,1)-1;
    if isempty(indmax) % we want the whole file
        indmax=numel(frameind.n);
    end
    startloc=frameind.nstarts(indmin);
    readsize=frameind.nstarts(indmax) + ...
        frameind.nsizes(indmax) - ...
        frameind.nstarts(indmin); % in bytes
    fseek(fidin,startloc,'bof');
    bt=fread(fidin,readsize/4,'float32'); % 4 bytes per float32
    bt=reshape(bt,ncol,numel(bt)/ncol)';
    fclose(fidin);
    fi=bt(:,ficol);
    x=bt(:,xcol);
    y=bt(:,ycol);
    t=bt(:,tcol);
    wx=bt(:,wxcol);
    wy=bt(:,wycol);
    vx=bt(:,vxcol);
    vy=bt(:,vycol);
    th=bt(:,thcol);

else % Open the .gdf input file, read the header, get initial time
    fidin=fopen(filename);
    header=fread(fidin,6,'int32');
    if header(1)~=magicnum
        error(['Sorry, ' filename ...
            ' does not appear to be a .gdf file.']);
    end
    ncol=header(3);
    buffer=fread(fidin,ncol,'float32')';
    tmin=max(buffer(tcol),timerange(1));

    fi=[];
    x=[];
    y=[];
    t=[];
    wx=[];
    wy=[];
    vx=[];
    vy=[];
    th=[];
    % Loop over frames, keeping all front points in each
    while ~feof(fidin)

        % Read in the next frame
        tt=buffer(tcol);
        fseek(fidin,-4*ncol,'cof'); % back up one sample
        bt=[];
        while buffer(tcol)==tt
            buffer=fread(fidin,ncol,'float32')';
            bt=[bt;buffer];
            if feof(fidin)
                break
            end
        end
        if ~feof(fidin)
            bt(end,:)=[]; % remove last sample, which has wrong time
        end    
        if tt < timerange(1) % before firstframe
            continue
        end
        if tt > timerange(2) % after lastframe
            break
        end
        fi=[fi;bt(:,ficol)];
        x=[x;bt(:,xcol)];
        y=[y;bt(:,ycol)];
        t=[t;tt*ones(size(bt,1),1)];
        wx=[wx;bt(:,wxcol)];
        wy=[wy;bt(:,wycol)];
        vx=[vx;bt(:,vxcol)];
        vy=[vy;bt(:,vycol)];
        th=[th;bt(:,thcol)];
    end % loop over frames
    tmax=tt-1;
    fclose(fidin);
    if numel(x)==0
        error(['Sorry, ' filename ...
            ' does not contain times between ' ...
            num2str(timerange(1)) ' and ' ...
            num2str(timerange(2)) '.']);
    end
end % if exist(indfile,'file')

if exist('region','var') && ~isempty(region)
    keepers=( (x>=region(1)) & (x<=region(1)+region(3)) & ...
        (y>=region(2)) & (y<=region(2)+region(4)) ); % exclude points near edges
    x=x(keepers);
    y=y(keepers);
    wx=wx(keepers);
    wy=wy(keepers);
    vx=vx(keepers);
    vy=vy(keepers);
    t=t(keepers);
    fi=fi(keepers);
    th=th(keepers);
end

% If requested, plot the front points
if noisy
    figure;
    ind=~isnan(wx);
    if sum(ind)>0
        maxlen=max( [getmaxlen(x(ind),y(ind),vx(ind),vy(ind)) ...
            getmaxlen(x(ind),y(ind),wx(ind),wy(ind)) ]); % single scale for both chemical and total
    else
        maxlen=1; % anything will do
    end
    quiver(x(ind),y(ind),vx(ind)/maxlen,vy(ind)/maxlen,0,'color','r');
    set(gca,'dataaspectratio',[1 1 1], ...
        'nextplot','add','ydir','reverse');
    quiver(x(ind),y(ind),wx(ind)/maxlen,wy(ind)/maxlen,0,'b');
    legend('chemical','total')
    [~,ends,~]=unique(fi,'legacy');
    begins=circshift(ends,1)+1;
    begins(1)=1;
    Nfr=numel(ends);
    for ii=1:Nfr
        ind=begins(ii):ends(ii);
        plot(x(ind),y(ind),'k');
    end
    if timerange(1)==timerange(2)
        title([num2str(numel(x)) ' front points in ' ...
            strrep(filename,'_','\_') ': time ' num2str(timerange(1))]);
    else
        title([num2str(numel(x)) ' front points in ' ...
            strrep(filename,'_','\_') ': times ' ...
            num2str(tmin) ' to ' num2str(tmax) ]);
    end
    if exist('region','var') && ~isempty(region)
        xlim(region(1)+[0 region(3)]);
        ylim(region(2)+[0 region(4)]);
    end
end

%Copyright 2018 Thomas Nevins and Douglas H. Kelley
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
