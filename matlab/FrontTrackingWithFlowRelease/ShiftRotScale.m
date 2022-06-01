function tracksSRS=ShiftRotScale(tracks,A,phi,res,frr,outfile)
% Usage: tracksSRS=ShiftRotScale(tracks,A,phi,res,frr,[outfile])
% Given particle tracks recorded by one camera, the shift and rotation
% required to align it with a second camera, and the resolution (in
% mm/pixel) and frame rate (in frames/s), ShiftRotScale shifts and rotates
% the tracks into the reference frame of the second camera, then converts
% the units to mm and s. The tracks must be contained in "tracks", which
% can be either a struct array of the sort produced by PredictiveTracker.m,
% or the name of a .gdf file with eight columns (track number, x, y, t, u,
% v, ax ay) of the sort produced by difftracks. Provide the x and y shift
% (in mm) in the two-element vector "A", the rotation (in radians) in
% "phi", the resolution (in mm/pixel) in "res", and the frame rate (in
% frames/s) in "frr". Alternately, provide in "frr" the name of the .seq
% movie from which the tracks were made, and ShiftRotScale will use the
% original timestamp from each frame (more precise). The shifted, rotated,
% and scaled tracks are output as the struct array "tracksSRS". If a
% filename is provided in "outfile", they are also saved to that file (in
% .gdf format and sorted by time, not track number). Requires read_gdf.m,
% Velocities.m, and write_gdf.m. Requires read_seq.m to pull timestamps
% from .seq file. See also camAlign.m.

% Written 13 March 2015 by Doug Kelley.
% Fixed sign errors 22 March 2015.
% Updated 23 March 2015 for compatibility with new version of camAlign.m
% (first scale, then rotate, then shift). 
% Updated 25 March 2015 for compatibility with PredictiveTracker output
% (struct arrays) in addition to difftracks output (.gdf files). 
% Updated 13 July 2015 to copy exact timestamps if movie is .seq. 

% Next: Make compatible with 3D data.

% -=- gdf file format -=-
Ncols=8;
trcol=1;
xcol=2;
ycol=3;
tcol=4;
ucol=5;
vcol=6;
axcol=7;
aycol=8;

% -=- Check inputs -=-
if nargin<1
    error(['Usage: tracksSRS = ' mfilename ...
        '(tracks,A,phi,res,frr,[outfile])']);
end
if numel(res)>1
    error('Sorry, res must be a scalar.')
end
if ischar(frr) % expect filename of .seq
    seqname=frr;
    if exist(seqname,'file')~=2
        error(['Sorry, cannot find file ' seqname '.'])
    end
    [~,~,info]=read_seq(seqname,0,0);
    frr=info.FrameRate;
    hasseq=true;
else
    hasseq=false;
end % if ischar(frr)

% -=- Read data -=-
if isstruct(tracks) % tracks from PredictiveTracker
    [u0,v0,x0,y0,t0,tr]=Velocities(tracks,[0 inf],0);
else % tracks from difftracks: .gdf file
    dd=read_gdf(tracks);
    if size(dd,2)~=Ncols
        error(['Sorry, ' tracks ' must have exactly ' num2str(Ncols) ...
            ' columns.'])
    end
    tr=dd(:,trcol);
    x0=dd(:,xcol);
    y0=dd(:,ycol);
    t0=dd(:,tcol);
    u0=dd(:,ucol);
    v0=dd(:,vcol);
    ax0=dd(:,axcol);
    ay0=dd(:,aycol);
end % if isstruct(tracks) 

% -=- Scale, rotate, and shift -=-
x0=x0*res; % change px --> mm
y0=y0*res;
x=x0*cos(phi)+y0*sin(phi)-A(1); % shift & rotate positions
y=-x0*sin(phi)+y0*cos(phi)-A(2);
if hasseq % use individual timestamps from .seq file
    Np=numel(t0); % particle count
    t=NaN(Np,1); % particle time
    ind=NaN(Np,1);
    frlist=(min(t0)-1):(max(t0)+1); % assumes no frames are skipped!
    disp(['Now reading ' seqname '.'])
    [~,tlist,~]=read_seq(seqname,[frlist(1) frlist(end)],0);
    for ii=1:Np
        ind(ii) = find(t0(ii)==frlist);
        t(ii)=tlist(ind(ii));
    end
    for ii=1:Np
        u0(ii)=u0(ii)*res/( tlist(ind(ii)+1) - tlist(ind(ii)-1) )/2; % two-frame central difference
        v0(ii)=v0(ii)*res/( tlist(ind(ii)+1) - tlist(ind(ii)-1) )/2;
    end
else % no .seq; use mean frame rate for everything
    t=t0/frr; % change fr --> s
    u0=u0*res*frr; % change px/fr --> mm/s
    v0=v0*res*frr;
end % if hasseq
u=u0*cos(phi)+v0*sin(phi); % rotate velocities
v=-u0*sin(phi)+v0*cos(phi);
if exist('ax0','var')
    if hasseq
        ax0(ii)=ax0(ii)*res/( tlist(ind(ii)+1) - tlist(ind(ii)-1) )^2; % two-frame central difference
        ay0(ii)=ay0(ii)*res/( tlist(ind(ii)+1) - tlist(ind(ii)-1) )^2;
    else
        ax0=ax0*res*frr^2; % change px/fr^2 --> mm/s^2
        ay0=ay0*res*frr^2;
    end % if hasseq
    ax=ax0*cos(phi)+ay0*sin(phi); % rotate accelerations
    ay=-ax0*sin(phi)+ay0*cos(phi);
end

% -=- Save to disk if requested -=-
if exist('outfile','var') && ~isempty(outfile)
    Np=numel(x);
    if isstruct(tracks) % need to initialize the data array
        if exist('ax0','var')
            dd=NaN(Np,8); % 6 columns: tr, x, y, t, u, v, ax, ay
        else
            dd=NaN(Np,6); % 6 columns: tr, x, y, t, u, v
        end
    end
    [dd(:,tcol),ind]=sort(t); % Sort by time.
    dd(:,trcol)=tr(ind);
    dd(:,xcol)=x(ind);
    dd(:,ycol)=y(ind);
    dd(:,ucol)=u(ind);
    dd(:,vcol)=v(ind);
    if exist('ax0','var')
        dd(:,axcol)=ax(ind);
        dd(:,aycol)=ay(ind);
    end
    write_gdf(dd,outfile);
end

% -=- Assemble struct array for output if requested -=-
if nargout>0
    [tr,ind]=sort(tr); % Sort by track number.
    x=x(ind);
    y=y(ind);
    t=t(ind);
    u=u(ind);
    v=v(ind);
    if exist('ax0','var')
        ax=ax(ind);
        ay=ay(ind);
    end
    [ttr,ends]=unique(tr,'legacy'); % split track-by-track
    begins=circshift(ends,[1 0])+1;
    begins(1)=1;
    Ntr=numel(ttr);
    if exist('ax0','var')
        tracksSRS=repmat( ...
            struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[], ...
            'AX',[],'AY',[]) ,[Ntr,1]);
    else
        tracksSRS=repmat( ...
            struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[]) ,[Ntr,1]);
    end
    for ii=1:Ntr
        ind=begins(ii):ends(ii);
        tracksSRS(ii).len=ends(ii)-begins(ii)+1;
        tracksSRS(ii).X=x(ind);
        tracksSRS(ii).Y=y(ind);
        tracksSRS(ii).T=t(ind);
        tracksSRS(ii).U=u(ind);
        tracksSRS(ii).V=v(ind);
        if exist('ax0','var')
            tracksSRS(ii).AX=ax(ind);
            tracksSRS(ii).AY=ay(ind);
        end
    end
end % if nargout>0

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
