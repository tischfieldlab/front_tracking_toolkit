function vfronts = TheFrontTracker(fname,thr,step,minarea,uname,span,ROI,framerange,res,invert,outfile,noisy,thick,bac,maxSpeed,maxPx)
% Usage: vfronts =
% TheFrontTracker(fname,thr,[step],[minarea],[uname],[span],[ROI],[framerange],[res],[invert],[outfile],[noisy],[thick],[bac],[maxSpeed],[maxPx])
% Given a movie of moving bright regions, FrontTracker locates the front at
% the edge of each region and finds its local velocity. If corresponding
% flow velocities are also given, both the chemical velocity of the front
% (due to molecular diffusion and reaction kinetics, presumed perpendicular
% to the front) and the total velocity of the front (which also includes
% advection by the underlying flow) are returned. If no flow velocities are
% given, the fronts are assumed to move in a stagnant fluid, and the total
% velocity equals the chemical velocity. 
%
% The movie must be saved in one of the following formats: - a series of
% image files with sequential numerical filenames, - an image stack in .tif
% or .gif format, - a Norpix .seq movie, or - a video format readable by
% Matlab's VideoReader.m.  Specify the movie in "fname" (e.g., '0*.png' or
% 'stack.tif', or 'movie.avi'). To be identified as a bright region, part
% of the movie must have brightness greater than "thr+bac" and area larger
% than "minarea".  The "bac" is a background image of the dish before
% reaction. If invert==1, the movie is first inverted (white to black). If
% thick==1, the front thickness will be calculated at all points also. The
% tracker will skip frames between each boundary based on the "step" value.
% Whatever integer is entered is the number of frames to step. Parts of the
% movie outside the region of interest "ROI" (specified as [x y width
% height], in pixels) are ignored.  Frames outside "framerange" (specified
% as [firstframe lastframe]) are also ignored.  Region boundaries are
% smoothed using a local line fit with span width "span". Corresponding
% velocity measurements, specified with filename "uname", must be in a .gdf
% file of the sort output by ShiftRotScale.m. The parameter "maxSpeed"
% controls how far the tracker looks for nearby fronts when determining
% front speed. The input is a distance in mm representing how far you want
% the tracker to look for another front if it were looking in the very next
% frame (regardless of "step"). The parameter "maxPx" is the mm distance
% you want to use for thickness curve fitting.
% 
% Results are returned in the structure "vfronts", whose fields "len", "T",
% "X", "Y", "Wx", "Wy", "Vx", "Vy", and "Thick" contain the length, time,
% horizontal coordinates, vertical coordinates, horizontal total
% velocities, vertical total velocities, horizontal chemical velocities,
% vertical chemical velocities, and front thickness of each point on each
% front, respectively. In the outputs, times are measured in seconds. If
% the camera resolution (in mm/pixel) is provided in "res", lengths are
% measured in mm; otherwise lengths are measured in pixels. If a filename
% is provided in "outfile", results are written to disk in .gdf format
% readable by vf.m. Parameters are always saved in a .mat to allow easier
% repeatability. If noisy~=0, the movie is repeated onscreen with overlaid
% region boundaries. If noisy==2, each movie frame is also saved to disk as
% an image. Requires Norpix2MATLAB.m, write_gdf.m, vp.m,
% logisticcdf2left.m, and logisticcdf2right.m. See also
% PredictiveTracker.m, ShiftRotScale.m, and vf.m.

% Written 11 March 2015 by Doug Kelley. 
% Updated 13 March 2015 to save results in mm and s, and to run faster with
% flat loops. 
% Accelerated 18 March 2015 by checking distance to old boundary before
% running polyxpoly and by including velocity_downsample. Also fixed
% plotting.
% Updated 22 March 2015 to save as .gdf and keep original origin (top left
% of image). 
% Fixed bugs 2 May 2015: default ROI was transposed, and ROI later too big.
% Updated 6 May 2015 to output chemical velocity and total velocity (not
% apparent velocity), and to allow movie formats other than .seq. Also
% fixed a bug in the plotting.
% Overhaul 10 May 2015: More intuitive algorithm, proven to work with dummy
% data. 
% Updated 11 May 2015 to handle frames that have no fronts.
% Updated 29 May 2015 to do Gaussian image smoothing before thresholding. 
% Updated 15 June 2015 to skip frames to avoid pixel jitter. Also fixed bug
% that was producing incorrect front speeds where boundaries were not
% closed.
% Updated 14 July 2015 to choose length of perpendiculars based on a
% reasonable speed limit. 
% Fixed advection bug 24 July 2015.
% Updated 16 November 2015 by Thomas Nevins to save parameters out in a Matlab file. 
% Updated 02 June 2016 by Thomas Nevins to use a background image.
% Updated 30 June 2016 by Thomas Nevins to calculate front thickness.
% Updated 01 August 2016 by Thomas Nevins. Fixed front tracking edges of
% ROI, lowered maxSpeed, and bugfix on front thicker complete.
% January 23 2017 by Thomas Nevins, version used on 2017 Chaos paper is finalized.
% March 16 2017, Thomas changed imaging in noisy to imagesc from image. Added
% user control for maxPx and maxSpeed.
% August 28 2017, Thomas found bug in calculation of WX and WY. Now using a
% changing interpolant as we loop through shifting, and using displacement
% for velocity instead of the interpolant.
% October 19 2017, Thomas added error handling for step variable. 
% Minor updates 27 November 2017: input processing, comments. 
% January 9 2018 Thomas prepared initial open access release corresponding to another paper submission.

if nargin<2
    error(['Usage: vfronts = ' mfilename ...
        '(fname,thr,step, [minarea],[uname],[span],[ROI],[framerange],[res],[invert],[outfile],[noisy],[thick],[bac])']);
end

% -=- Set defaults -=-----------------------------------------------------
step_default=1; %Don't skip if no skip given.
minarea_default=0; % px^2
span_default=50; % points along boundary
res_default=1; % mm/px
invert_default=0;%Default to not inverting
noisy_default=0; %Default to no plot
thick_default=0;%Default to no thickness
framerange_default=[-inf inf]; %Default to all the data
bac_default=0; %Default to subtracting nothing
maxSpeed_default=1.5; %Default maximum speed.
frr_default=1; %Default to 1 fr/s.

% -=- Parameters -=-------------------------------------------------------
hsize=6; % Gaussian filter size 6 px
sigma=2; % Gaussian filter width (standard deviation) 2 px
pausetime=0.2; % seconds to pause between frames when plotting
savedirname='frontsmovie';
figsize=[800 600]; % figure width and height, in pixels
velocity_downsample=3; % calculate velocity only at every 5th point
ficol=1; % front index in column 1 of .gdf file
xcol=2; % x in column 2 of .gdf file
ycol=3; % y in column 3 of .gdf file
tcol=4; % t in column 4 of .gdf file
wxcol=5; % wx in column 5 of .gdf file (front velocity w/ flow; w=v+u)
wycol=6; % wy in column 6 of .gdf file
vxcol=7; % vx in column 7 of .gdf file (front velocity w/o flow)
vycol=8; % vy in column 8 of .gdf file

% -=- Parse inputs -=-----------------------------------------------------
if ~exist('step','var') || isempty(step)
    step=step_default;
end
if ~exist('framerange','var') || isempty(framerange)
    framerange=framerange_default;
elseif numel(framerange)==1
    framerange=framerange*[1 1];
end
if ~exist('minarea','var') || isempty(minarea)
    minarea=minarea_default;
end
if ~exist('uname','var') || isempty(uname)
    warning('No flow file specified; assuming zero flow.');
    has_flow=false;
    uname=[];
elseif ~exist(uname,'file')
    error(['Sorry, cannot find ' uname '.'])
else
    has_flow=true;
end
if ~exist('span','var') || isempty(span)
    span=span_default;
end
if ~exist('res','var') || isempty(res)
    if has_flow
        error('Sorry, if uname is given, res must also be given.')
    else
        res=res_default;
    end
end
maxPx_default=ceil(1.5/res); %Default maxPx in pixels
if ~exist('invert','var') || isempty(invert)
    invert=invert_default;
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisy_default;
end
if ~exist('thick','var') || isempty(thick)
    thick=thick_default;
end
if thick ==1
    thickcol=9; % Thickness in column 9 of .gdf file
    Ncols=9; %Make this many columns
else
    Ncols=8; % number of columns in .gdf file
end
if ~exist('bac','var') || isempty(noisy)
    bac=bac_default;
end
if ~exist('maxPx','var') || isempty(maxPx)
    maxPx=maxPx_default;
else
    maxPx=ceil(maxPx/res);
end
if ~exist('maxSpeed','var') || isempty(maxSpeed)
    maxSpeed=maxSpeed_default;
end

% -=- Determine movie format, get info, get set up -=---------------------
gcp; %Starts parallel pool
[filepath,~,ext]=fileparts(fname);
names=dir(fname);
thickguess=maxPx*0.05;
hv=fspecial('gaussian',hsize,sigma);
ftl=fittype('logisticcdf2left(x,a,b,mu,I)'); %Custom fittype
ftr=fittype('logisticcdf2right(x,a,b,mu,I)'); 

%Code to handle various types of supported movies
if numel(names)==0
    error(['Sorry, cannot find ' fname '.']);
elseif ( strcmpi(ext,'.avi') || ...
    strcmpi(ext,'.mp4') || strcmpi(ext,'.mov') || strcmpi(ext,'.mj2') || ...
    strcmpi(ext,'.wmv') || strcmpi(ext,'.mpg') || strcmpi(ext,'.asf') || ...
    strcmpi(ext,'.asx') || strcmpi(ext,'.m4v') || strcmpi(ext,'.ogg') ) % single movie file
    movtype='avi';
    if isempty(which('VideoReader.m')) % check for req'd helper function
        error(['Sorry, reading .avi files requires ' ...
            'VideoReader.m.'])
    end
    vobj=VideoReader(fullfile(filepath,names.name));
    frrate=vobj.FrameRate;
    if strcmp(vobj.VideoFormat(1:3),'RGB')
        color_depth=2^(vobj.BitsPerPixel/3); % RGB
    elseif strcmp(vobj.VideoFormat,'Indexed')
        error('Sorry, indexed color is not supported.')
    else
        color_depth=2^vobj.BitsPerPixel; % grayscale
    end
    ht=vobj.Height;
    wd=vobj.Width;
    tmin=max([framerange(1) 1]);
    tmax=min([framerange(2) vobj.NumberOfFrames]);
    t0=tmin/frrate;
elseif strcmpi(ext,'.seq')
    movtype='seq';
    [~,t0,movinfo]=read_seq(fname,0,0); %See the read_seq.m function
    frrate=movinfo.FrameRate;
    color_depth=2^(movinfo.ImageBitDepthReal);
    ht=movinfo.ImageHeight;
    wd=movinfo.ImageWidth;
    tmin=max([framerange(1) 0]);
    tmax=min([framerange(2) movinfo.AllocatedFrames-1]);
elseif numel(names)==1 && ( strcmpi(ext,'.tif') || ...
        strcmpi(ext,'.tiff') || strcmpi(ext,'.gif') ) % single file, looks like an image stack
    movtype='stack';
    movinfo=imfinfo(fullfile(filepath,names.name));
    warning('Cannot determine frame rate; assuming 1 fr/s.')
    frrate=1;
    color_depth=2^(movinfo(1).BitDepth);
    ht=movinfo(1).Height;
    wd=movinfo(1).Width;
    tmin=max([framerange(1) 1]);
    tmax=min([framerange(2) numel(movinfo)]);
    t0=tmin/frrate;
else % multiple files probably a series of images
    movtype='images';
    movinfo=imfinfo(fullfile(filepath,names(1).name));
    warning('Cannot determine frame rate; assuming 1 fr/s.')
    frrate=1;
    color_depth=2^(movinfo.BitDepth);
    ht=movinfo.Height;
    wd=movinfo.Width;
    tt=NaN(size(names));
    for ii=1:numel(names)
        [~,myname]=fileparts(names(ii).name);
        tt(ii)=str2double(myname);
    end
    if any(isnan(tt)) % at least one filename was not a number, so just make ordinals
        tt=1:numel(tt);
    end
    tmin=max([framerange(1) min(tt)]);
    tmax=min([framerange(2) max(tt)]);
    if tmax<tmin
        error(['Sorry, ' fname ...
            ' does not contain frames in range [' num2str(framerange) ...
            '].'])
    end
    names = names( (tt>=tmin) & (tt<=tmax) );
    t0=tmin/frrate;
end % if numel(names)==0
if ~exist('ROI','var') || isempty(ROI) % if no ROI specified...
    ROI=[1 1 wd ht]; % ... use whole image
end
tt=tmin:step:tmax; % frames
Nf=numel(tt); % frame count
vfronts=cell(Nf,1);
limc=[inf -inf]; % for caxis

% -=- Loop through frames -=----------------------------------------------
% The plan: in each frame, locate fronts by finding regions brighter than
% thr, getting their boundaries, and smoothing. To get chemical velocity,
% consider frames in pairs. Use known flow velocity to remove flow
% displacement from the second frame, so the remaining displacement is due
% to chemical velocity alone. Then find chemical velocity by drawing lines
% locally perpendicular to the earlier front, through the (shifted) later
% front. 

for ii=0:Nf

    if ii>0
        % -=- Already have boundaries for current frame -=-
        Nbo=Nbn; % count of boundaries
        bdszold=bdsz-1; % count of points on each boundary; we'll drop repeated final point
        bdryold=bdrynew; % front boundaries
        txold=tx; % unit tangent vector at each boundary point, used in thick
        tyold=ty;
        
        if thick==1
            bdrypxold=bdrypx;
        end
        for jj=1:Nbo
            bdryold{jj}=bdryold{jj}(1:bdszold(jj),:); % drop repeated final point
            if thick == 1
                bdrypxold{jj}=bdrypxold{jj}(1:bdszold(jj),:);
            end
            txold{jj}=txold{jj}(1:bdszold(jj));
            tyold{jj}=tyold{jj}(1:bdszold(jj));
        end
        oldt=t1; % time
        Npo=Npn-Nbo; % total count of points on boundaries (dropping one from each)
        bold=b2;
        tx_shift=txold; % unit tangent vector at each boundary point, used in velocity
        ty_shift=tyold;
    else
        Nbo=0;
    end

    if ii<Nf
        % -=- Read fronts movie, find boundaries for next frame -=-
        switch movtype
            case('seq')
                [b,t1]=read_seq(fname,tt(ii+1),0); % read frame
                b=im2double(b)*255;
            case('avi')
                b = double(read(vobj,tmin+ii));
                t1=(tmin+ii)/frrate;
            case('stack')
                b = double(imread(fullfile(filepath,names.name),tmin+ii));
                t1=(tmin+ii)/frrate;
            case('images')
                b = double(imread(fullfile(filepath,names(ii+1).name)));
                t1=(tmin+ii)/frrate;
        end % switch movtype
        b=b-bac;
        b2=b;
        b=imfilter(b,hv,'replicate'); % Gaussian low-pass filter
        b=b(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
        b2=b2(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
        if invert
            b=color_depth-1-b;
            b2=color_depth-1-b2;
        end
        limc1=[prctile(b(:),1) prctile(b(:),99)];
        if limc1(1)<limc(1)
            limc(1)=limc1(1);
        end
        if limc1(2)>limc(2)
            limc(2)=limc1(2);
        end
        [bdrynew,lab]=bwboundaries(b>thr); % region boundaries and labels
        if thick==1
            bdrypx=bdrynew; %Thick needs to work in pixels until the very end.
        end
        Nbn=numel(bdrynew); % boundary count
        for jj=1:Nbn
            bdrynew{jj}(:,2)=(bdrynew{jj}(:,2)+ROI(1))*res; % mm from original origin
            bdrynew{jj}(:,1)=(bdrynew{jj}(:,1)+ROI(2))*res; % mm from original origin
        end
        
        %Cut out the small fronts.
        bdsz=NaN(Nbn,1);
        intlab=NaN(Nbn,1);
        for jj=1:Nbn % check for boundaries that are too short
            bdsz(jj)=size(bdrynew{jj},1);
            intlab(jj)=jj;
        end
        ind =(bdsz>=span+1);
        bdrynew=bdrynew(ind); % drop 'em
        if thick == 1
            bdrypx=bdrypx(ind);
        end
        bdsz=bdsz(ind);
        intlab=intlab(ind);
        Nbn=numel(bdrynew); % boundary count, again
        % *By doing the span check first we get major speed increase.
        % Area takes much longer than length, so the spancut first limits
        % how many area calculations we need.
        a=NaN(Nbn,1); 
        for jj=1:Nbn %Check for boundaries that enclose small area.
            a(jj)=sum(lab(:)==intlab(jj));
        end
        ind=(a>=minarea);
        bdrynew=bdrynew(ind);
        bdsz=bdsz(ind);
        Nbn=numel(bdrynew);
        Npn=sum(bdsz);
        
        if ii>0
            disp(['Tracking ' num2str(Nbn) ' fronts in frame ' ...
                num2str(tt(ii+1)) ', ' num2str(t1-t0) ' seconds (' num2str(ii+1) ...
                ' of ' num2str(Nf) ').'])
        end

% -=- Smooth fronts and find perpendicular angle at each point -=-
        tx=cell(Nbn,1); % tangential vector, x component
        ty=cell(Nbn,1); % tangential vector, y component
        
        if thick==1 %If thick, need to process bdrypx separately. 
            for jj=1:Nbn 
                bdrynew{jj}=[bdrynew{jj}(end-span:end-1,:) ; bdrynew{jj} ; ...
                    bdrynew{jj}(2:span,:)]; % append for circular smoothing
                bdrynew{jj}(:,1)=smooth(bdrynew{jj}(:,1),span,'lowess'); % smooth
                bdrynew{jj}(:,2)=smooth(bdrynew{jj}(:,2),span,'lowess');
                bdrynew{jj}=bdrynew{jj}(span+1:end-span+1,:); % drop appended points
                bdrypx{jj}=[bdrypx{jj}(end-span:end-1,:) ; bdrypx{jj} ; bdrypx{jj}(2:span,:)];
                bdrypx{jj}(:,1)=smooth(bdrypx{jj}(:,1),span,'lowess');
                bdrypx{jj}(:,2)=smooth(bdrypx{jj}(:,2),span,'lowess');
                bdrypx{jj}=bdrypx{jj}(span+1:end-span+1,:);
                tx{jj} = circshift(bdrynew{jj}(:,2),[1 0]) - ...
                    circshift(bdrynew{jj}(:,2),[-1 0]);
                ty{jj} = circshift(bdrynew{jj}(:,1),[1 0]) - ...
                    circshift(bdrynew{jj}(:,1),[-1 0]);
                ds=sqrt(tx{jj}.^2+ty{jj}.^2);
                tx{jj}=tx{jj}./ds;
                ty{jj}=ty{jj}./ds;
            end % for jj=1:Nbn
        else % Don't do the thick stuff
            for jj=1:Nbn 
                bdrynew{jj}=[bdrynew{jj}(end-span:end-1,:) ; bdrynew{jj} ; ...
                    bdrynew{jj}(2:span,:)]; % append for circular smoothing
                bdrynew{jj}(:,1)=smooth(bdrynew{jj}(:,1),span,'lowess'); % smooth
                bdrynew{jj}(:,2)=smooth(bdrynew{jj}(:,2),span,'lowess');
                bdrynew{jj}=bdrynew{jj}(span+1:end-span+1,:); % drop appended points
                tx{jj} = circshift(bdrynew{jj}(:,2),[1 0]) - ...
                    circshift(bdrynew{jj}(:,2),[-1 0]);
                ty{jj} = circshift(bdrynew{jj}(:,1),[1 0]) - ...
                    circshift(bdrynew{jj}(:,1),[-1 0]);
                ds=sqrt(tx{jj}.^2+ty{jj}.^2);
                tx{jj}=tx{jj}./ds;
                ty{jj}=ty{jj}./ds;
            end % for jj=1:Nbn
        end
        
        
% The prior frame has fronts, so measure their velocity.        
        if Nbo>0 
            bdry_shifted=bdryold;
 % -=- Subtract displacement due to flow -=-
 % Plan: to advect fronts, use all flow data from skipped
 % frames. Work backward in time from final "particle"
 % locations. tu0 is earlier time and tu1 is later time.
            if has_flow
                [ux0,uy0,x0,y0,tu0]=vp(uname, ... % setting up for the loop
                    (tt(ii)+[-1/2 1/2])/frrate,[],0);
                [~,ind]=unique(x0+sqrt(-1)*y0); % unique in both spatial coordinates
                x0=x0(ind); % drop duplicate points
                y0=y0(ind);
                ux0=ux0(ind);
                uy0=uy0(ind);
                tu0=tu0(ind);
                tu0=tu0(1); % assuming all times same and framerate matches!
                for jj=tt(ii)+1:1:tt(ii+1) % consider skipped frames individually; jj is frame number
                    ux1=ux0; % use leftovers from last iteration: what was earlier becomes later
                    uy1=uy0;
                    x1=x0;
                    y1=y0;
                    tu1=tu0;
                    UX=scatteredInterpolant(x1,y1,ux1);
                    UY=scatteredInterpolant(x1,y1,uy1);
                    [ux0,uy0,x0,y0,tu0]=vp(uname, ...
                        (jj+[-1/2 1/2])/frrate,[],0);
                    [~,ind]=unique(x0+sqrt(-1)*y0); % unique in both spatial coordinates
                    x0=x0(ind); % drop duplicate points
                    y0=y0(ind);
                    ux0=ux0(ind);
                    uy0=uy0(ind);
                    tu0=tu0(ind);
                    tu0=tu0(1); % assuming all times same and framerate matches!
                    for ll=1:Nbo
                        bdry_shifted{ll}(:,2) = bdry_shifted{ll}(:,2) + ...
                            UX(bdry_shifted{ll}(:,2),bdry_shifted{ll}(:,1))*(tu0-tu1); % shifted new front
                        bdry_shifted{ll}(:,1) = bdry_shifted{ll}(:,1) + ...
                            UY(bdry_shifted{ll}(:,2),bdry_shifted{ll}(:,1))*(tu0-tu1); % shifted new front
                    end
                end % for jj=tt(ii)+1:1:tt(ii+1)
                %Calculate tangents to shifted front, will superseed
                %tangents in first section of computation loop.
                for ll=1:Nbo
                    tx_shift{ll} = circshift(bdry_shifted{ll}(:,2),[1 0]) - ...
                        circshift(bdry_shifted{ll}(:,2),[-1 0]);
                    ty_shift{ll} = circshift(bdry_shifted{ll}(:,1),[1 0]) - ...
                        circshift(bdry_shifted{ll}(:,1),[-1 0]);
                    ds=sqrt(tx_shift{ll}.^2+ty_shift{ll}.^2);
                    tx_shift{ll}=tx_shift{ll}./ds;
                    ty_shift{ll}=ty_shift{ll}./ds;
                end %Loop through all old boundaries
                %Store flow displacement velocities.
                bdryprev=vertcat(bdryold{:});
                bdry_shifted=vertcat(bdry_shifted{:});
                ux_large=(bdry_shifted(:,2)-bdryprev(:,2))./(t1-oldt);
                uy_large=(bdry_shifted(:,1)-bdryprev(:,1))./(t1-oldt);
            else
                bdry_shifted=vertcat(bdry_shifted{:});
            end % if has_flow

% -=- Calculate chemical velocity, point-by-point -=- 
            tx_shift=vertcat(tx_shift{:});
            ty_shift=vertcat(ty_shift{:});
            ends=cumsum(bdszold); % indices where bdrys end
            begins=circshift(ends,[1 0])+1; % indices where bdrys start
            begins(1)=1;
            vdi=1:velocity_downsample:Npo;
            Npds=numel(vdi); % number of points after downsampling
            vx_small=NaN(Npds,1); % chemical velocity components, from perpendicular distance
            vy_small=NaN(Npds,1);
            bdryold1_small=bdry_shifted(vdi,:);
            tx_small=tx_shift(vdi);
            ty_small=ty_shift(vdi);
            maxDist=(maxSpeed/frrate)*step; % max distance a front can advance            
            parfor jj=1:Npds %Loop over points we want data at
                perpx=bdryold1_small(jj,2)*[1 1]+maxDist.*ty_small(jj)*[1 -1]; % ends of perpendicular segment
                perpy=bdryold1_small(jj,1)*[1 1]-maxDist.*tx_small(jj)*[1 -1];
                if (min(perpx)<=ROI(1)*res || max(perpx)>=(ROI(1)+ROI(3))*res || ...
                        min(perpy)<=ROI(2)*res || max(perpy)>=(ROI(2)+ROI(4))*res)
                    vx_small(jj)=NaN; %Drop boundary points
                    vy_small(jj)=NaN;
                else
                    ddmin=inf;
                    for ll=1:Nbn % loop over (shifted) new boundaries
                        minsqrdist=min( ...
                            (bdrynew{ll}(:,2)-bdryold1_small(jj,2)).^2 + ...
                            (bdrynew{ll}(:,1)-bdryold1_small(jj,1)).^2 );
                        if minsqrdist<=maxDist^2 % close enough to intersect
                            [ix,iy]=polyxpoly(perpx,perpy, ...
                                bdrynew{ll}(:,2),bdrynew{ll}(:,1)); % find intersections
                            dx=ix-bdryold1_small(jj,2);
                            dy=iy-bdryold1_small(jj,1);
                            [dd,ind]=min(sqrt(dx.^2+dy.^2)); % shortest distance from bdry pt to an old bdry
                            if dd<ddmin % closer than the existing closest intersection?
                                ddmin=dd; % remember this intersection
                                vx_small(jj)=dx(ind)/(t1-oldt); % mm/s
                                vy_small(jj)=dy(ind)/(t1-oldt);
                            end
                        end % if minsqrdist<=maxDist.^2
                    end % for ll=1:Nbn (loop over old boundaries)
                end % Edge Check
            end % parfor ii=1:Npds (loop over points)
            
% -=- Calculate Front Thickness, if asked -=-
            if thick ==1
                bdrypxold_small=vertcat(bdrypxold{:});
                bdrypxold_small=bdrypxold_small(vdi,:);
                fthick=NaN(Npds,1); %Will eventually store the thickness value
                txold=vertcat(txold{:});
                tyold=vertcat(tyold{:});
                tx_small=txold(vdi);
                ty_small=tyold(vdi);
                parfor jj=1:Npds 
                    perpxt=bdrypxold_small(jj,2)*[1 1]+maxPx.*ty_small(jj)*[1 -1];
                    perpyt=bdrypxold_small(jj,1)*[1 1]-maxPx.*tx_small(jj)*[1 -1];
                    if (min(perpxt)<= 0 || max(perpxt)>=ROI(3) || ...
                        min(perpyt)<= 0 || max(perpyt)>=ROI(4)) %Get rid of stuff too close to boundary
                        fthick(jj)=NaN;
                    elseif (isnan(perpxt(1)) || isnan(perpxt(2)) || isnan(perpyt(1)) ...
                            || isnan(perpyt(2)) || isinf(perpxt(1)) || isinf(perpxt(2)) ...
                            || isinf(perpyt(1)) || isinf(perpyt(2)))
                        fthick(jj)=NaN;
                    else   
                        c=improfile(bold,perpxt,perpyt,'bicubic'); %Grabs image profile in c.
                        if ~isnan(sum(c))
                            NProf=numel(c);
                            pxdist=(1/(NProf-1))*res*sqrt((perpxt(2)-perpxt(1))^2+...
                            (perpyt(2)-perpyt(1))^2); %mm distance between two indices of c. Not trivial due to improfile.
                            pxs=(1:numel(c))'-1; %x coords for eventual fit
                            mguess=(numel(c)-1)/2;% Making the guesses for the fitter.
                            high_c=max(c);
                            low_c=min(c);
                            topleft=mean(c(floor(mguess)+1:end))<=mean(c(1:floor(mguess))); %Which of the two forms to use
                            startpt=[thickguess, low_c, high_c-low_c,mguess];
                            if topleft
                                f0=fit(pxs,c,ftl,'StartPoint',startpt); %do a curve fit
                            else
                                f0=fit(pxs,c,ftr,'StartPoint',startpt); %do a curve fit
                            end %if topleft
                            
%                             if f0.b>0 %Diagnostic
%                                 plot(mean(perpxt),mean(perpyt),'bo');
%                             else
%                                 plot(mean(perpxt),mean(perpyt),'ro');
%                             end
                            
%                             if mod(jj,10)==0 %Diagnostic
%                                 figure;
%                                 hold on
%                                 plot(pxs,c,'bo');
%                                 plot(pxs,f0(pxs),'r');
%                             end
                            fthick(jj)=f0.I*pxdist; %Finally turn thickness into metric.
                            if abs(fthick(jj))>=maxPx*2*res
                                fthick(jj)=NaN;
                            end
                        else
                            fthick(jj)=NaN;
                        end %Checking for NaN in C. Inelegant solution. 
                    end %if perpendicular is inbounds
                end %parfor jj=1:Npds
                fthick=abs(fthick); 
            end %if thick ==1
            

% -=- Calculate total velocity too, and consolidate -=-
            vx=NaN(Npo,1);
            vy=NaN(Npo,1);
            wx=NaN(Npo,1);
            wy=NaN(Npo,1);
            if thick==1
                thickness=NaN(Npo,1);
                thickness(vdi)=fthick;
            end
            vx(vdi)=vx_small;
            vy(vdi)=vy_small;
            if has_flow
                wx(vdi) = vx_small+ux_large(vdi); %UX(bdryold1_small(:,2),bdryold1_small(:,1));
                wy(vdi) = vy_small+uy_large(vdi); %UY(bdryold1_small(:,2),bdryold1_small(:,1));
            else
                wx(vdi)=vx_small;
                wy(vdi)=vy_small;
            end

        end % if Nbo>0

        if ii>0
% -=- Arrange as structure for output -=-
            bdrystruct=repmat( struct('len',[],'T',oldt-t0, ...
                'X',[],'Y',[],'Wx',[],'Wy',[],'Vx',[],'Vy',[],'Thickness',[]) ,[Nbo 1]);
            bdryold=vertcat(bdryold{:});
            for jj=1:Nbo
                bdrystruct(jj).len=bdszold(jj); % unitless (point count)
                ind=begins(jj):ends(jj);
                bdrystruct(jj).X=(bdryold(ind,2)); % mm from original origin
                bdrystruct(jj).Y=(bdryold(ind,1)); % mm from original origin
                bdrystruct(jj).Wx=wx(ind); % mm/s
                bdrystruct(jj).Wy=wy(ind); % mm/s
                bdrystruct(jj).Vx=vx(ind); % mm/s
                bdrystruct(jj).Vy=vy(ind); % mm/s
                if thick==1
                    bdrystruct(jj).Thickness=thickness(ind); % mm
                end
            end
            vfronts{ii}=bdrystruct;
        end % if ii>0     
        
    elseif ii==Nf % if ii<Nf

% -=- Handle the final frame differently -=-
        ends=cumsum(bdszold); % indices where bdrys end
        begins=circshift(ends,[1 0])+1; % indices where bdrys start
        begins(1)=1;
        wx=NaN(Npn,1);
        wy=NaN(Npn,1);
        vx=NaN(Npn,1);
        vy=NaN(Npn,1);
        if thick==1
            thickness=NaN(Npn,1);
        end
        bdrystruct=repmat( struct('len',[],'T',oldt-t0, ...
            'X',[],'Y',[],'Wx',[],'Wy',[],'Vx',[],'Vy',[],'Thickness',[]) ,[Nbo 1]);
        bdryold=vertcat(bdryold{:});
        for jj=1:Nbo
            bdrystruct(jj).len=bdszold(jj); % unitless (point count)
            ind=begins(jj):ends(jj);
            bdrystruct(jj).X=(bdryold(ind,2)); % mm from original origin
            bdrystruct(jj).Y=(bdryold(ind,1)); % mm from original origin
            bdrystruct(jj).Wx=wx(ind); % mm/s
            bdrystruct(jj).Wy=wy(ind); % mm/s
            bdrystruct(jj).Vx=vx(ind); % mm/s
            bdrystruct(jj).Vy=vy(ind); % mm/s
            if thick==1
                bdrystruct(jj).Thickness=thickness(ind); % mm
            end
        end
        vfronts{ii}=bdrystruct;

    end % if ii<Nf
    
end % for ii=1:Nf (loop over time)

% -=- Consolidate everything and produce output -=-
vfronts=vertcat(vfronts{:}); % turn cell array of structs into single struct
if exist('outfile','var') && ~isempty(outfile)
    Npall=sum(vertcat(vfronts.len)); % total number of points, all times
    dd=NaN(Npall,Ncols); % big array for all vfronts information
    dd(:,xcol)=vertcat(vfronts.X);
    dd(:,ycol)=vertcat(vfronts.Y);
    dd(:,wxcol)=vertcat(vfronts.Wx);
    dd(:,wycol)=vertcat(vfronts.Wy);
    dd(:,vxcol)=vertcat(vfronts.Vx);
    dd(:,vycol)=vertcat(vfronts.Vy);
    if thick==1
        dd(:,thickcol)=vertcat(vfronts.Thickness);
    end
    startind=0;
    for ii=1:numel(vfronts); % expand T and front index
        ind=startind+(1:vfronts(ii).len);
        dd(ind,ficol)=ii*ones(vfronts(ii).len,1);
        dd(ind,tcol)=vfronts(ii).T*ones(vfronts(ii).len,1);
        startind=startind+vfronts(ii).len;
    end
    disp(['Writing fronts to ' outfile '.'])
    write_gdf(dd,outfile); % write a .gdf
end % if exist('outfile','var') && ~isempty(outfile)

% -=- Plot if requested -=------------------------------------------------
if noisy
    if isnumeric(noisy) && noisy>1
        disp(['Plotting and saving frames. ' ...
            'Please do not cover the figure window!'])
        if exist(fullfile('.',savedirname),'file')~=7
            mkdir(savedirname)
        end
    else
        disp('Plotting...')
    end
    defaultpos=get(0,'DefaultFigurePosition');
    figure('units','pixels','position',[defaultpos(1:2) figsize]);
    axes('nextplot','add','dataaspectratio',[1 1 1],'ydir','reverse','clim',limc);
    xlabel(['thr = ' num2str(thr) ...
        ', span = ' num2str(span)  ', minarea = ' num2str(minarea) ...
        ', invert = ' num2str(invert)], 'interpreter','none');
    switch movtype
        case('seq')
            b=read_seq(fname,0,0);
        case('avi')
            b = double(read(vobj,tmin));
        case('stack')
            b = double(imread(fullfile(filepath,names.name),tmin));
        case('images')
            b = double(imread(fullfile(filepath,names(1).name)));
    end % switch movtype
    b=b(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
    hi=imagesc((ROI(1)+(0:ROI(3))),(ROI(2)+(0:ROI(4))),b);
    xlim(0.5+ROI(1)+[0 ROI(3)]);
    ylim(0.5+ROI(2)+[0 ROI(4)]);
    colormap(gray); % has no effect if image is color b/c RGB values override
    allT=vertcat(vfronts(:).T);
    for ii=1:Nf
        if ii>1
            delete(hl)
        end
        switch movtype
            case('seq')
                [b,t1]=read_seq(fname,tt(ii),0);
                b=im2double(b)*255;
            case('avi')
                b = double(read(vobj,tmin+ii-1));
                t1=(tmin+ii-1)/frrate;
            case('stack')
                b = double(imread(fullfile(filepath,names.name),tmin+ii-1));
                t1=(tmin+ii-1)/frrate;
            case('images')
                b = double(imread(fullfile(filepath,names(ii).name)));
                t1=(tmin+ii-1)/frrate;
        end % switch movtype
        b=b-bac;
        b=b(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
        set(hi,'cdata',b);
        title([fname ' frame ' num2str(tt(ii)) ', ' num2str(t1-t0) ...
            ' seconds (' num2str(ii) ' of ' num2str(Nf) ')'], ...
            'interpreter','none');
        list=find(allT==(t1-t0));
        hl=NaN(size(list));
        for jj=1:numel(list)
            hl(jj)=plot(vfronts(list(jj)).X/res, ... % plot in px, not mm
                vfronts(list(jj)).Y/res,'r', ...
                'linewidth',2,'displayname',['front ' num2str(list(jj))]);
        end
        drawnow
        pause(pausetime); 
        if isnumeric(noisy) && noisy>1
            snap=getframe(gca);
            imwrite(snap.cdata, ...
                fullfile('.',savedirname,[num2str(ii,'%05.0f') '.png'])); 
        end 
    end % for ii=1:Nf
    disp('Done.')
end % if noisy

% Code for saving the input parameters.
s.fname=fname;
s.thr=thr;
s.step=step;
s.minarea=minarea;
s.uname=uname;
s.span=span;
s.ROI=ROI;
s.framerange=framerange;
s.res=res;
s.invert=invert;
s.outfile=outfile;
s.noisy=noisy;
s.bac=bac;
s.thick=thick;
s.maxPx=maxPx*res;
s.maxSpeed=maxSpeed;
save([outfile(1:end-4) '.mat'],'-struct','s');

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
