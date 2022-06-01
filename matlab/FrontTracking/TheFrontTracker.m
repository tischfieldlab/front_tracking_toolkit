function vfronts = TheFrontTracker(mov,thr,step,minarea,span,ROI,framerange,res,frrate,invert,noisy,bac,maxSpeed,maxPx)
% Usage: vfronts =
% TheFrontTracker(mov,thr,step,[minarea],[span],[ROI],[framerange],[res],[frrate],[invert],[noisy],[bac],[maxSpeed],[maxPx])
% Given a movie "mov" of moving bright regions, TheFrontTracker locates the
% fronts at the edges of each region and tracks them over time, measuring
% local velocity and thickness at many points along each front. Specify
% "mov" as a three-dimensional matrix whose elements are pixel
% brightnesses, whose first two dimensions vary with spatial coordinates,
% and whose third dimension varies with time (x,y,t).
%
% To be identified as a bright region, part of the movie must have
% brightness greater than the background image by an amount "thr" and have
% area larger than "minarea" (in square pixels). Optionally provide a
% background image in "bac".  If invert==1, the movie is first inverted
% (white to black). The tracker will skip "step" movie frames between
% velocity and thickness measurements, for smoothing. 

% Parts of the movie outside the region of interest "ROI" (specified as [x
% y width height], in pixels) are ignored.  Frames outside "framerange"
% (specified as [firstframe lastframe]) are also ignored.  Region
% boundaries are smoothed using a local line fit with span width "span".
% The distance "maxSpeed" (in pixels per frame) controls how far the
% tracker looks for nearby fronts when determining front speed.  The
% parameter "maxPx" is the pixel distance you want to use for thickness
% curve fitting.
% 
% Results are returned in the structure "vfronts", whose fields "len", "T",
% "X", "Y", "Vx", "Vy", and "Thickness" contain the length (point count),
% time, horizontal coordinates, vertical coordinates, horizontal chemical
% velocities, vertical chemical velocities, and front thickness of each
% point on each front, respectively. If "frrate" (in frames/second) and
% "res" (in mm/pixel) are provided, outputs are reported in units of
% seconds and mm; otherwise, frames and pixels are used.  If noisy~=0, the
% movie is repeated onscreen with overlaid region boundaries. If noisy==2,
% each movie frame is also saved to disk as an image.
%
% See also VF.m.

% Written March 2017 by Thomas D. Nevins and Douglas H. Kelley, based
% largely on FrontTracker.m. 

% -=- Set defaults -=-----------------------------------------------------
minarea_default=0; % px^2
span_default=50; % points along boundary
res_default=1; % mm/px
invert_default=0; %Default to not inverting
noisy_default=0; %Default to no plot
framerange_default=[-inf inf]; %Default to all the data
bac_default=0; %Default to subtracting nothing
maxPx_default=20; %Default to 20 pixel length for thickness
if exist('res','var')
    maxSpeed_default=20*res; % Default to 20 pixels as max speed jump.
else
    maxSpeed_default=20*res_default;
end
frr_default=1; %Default to 1 fr/s.

% -=- Parameters -=-------------------------------------------------------
thickguess=1; % guess of thickness in pixels.
pausetime=0.2; % seconds to pause between frames when plotting
savedirname='frontsmovie';
figsize=[800 600]; % figure width and height, in pixels
velocity_downsample=5; % calculate velocity only at every 2nd point

% -=- Parse inputs -=-----------------------------------------------------
if nargin<3
    error(['Usage: vfronts = ' mfilename ...
        '(mov,thr,step,[minarea],[span],[ROI],[framerange],[res],[frrate],[invert],[noisy],[bac],[maxSpeed],[maxPx])']);
end
if ~exist('framerange','var') || isempty(framerange)
    framerange=framerange_default;
elseif numel(framerange)==1
    framerange=framerange*[1 1];
end
if ~exist('minarea','var') || isempty(minarea)
    minarea=minarea_default;
end
if ~exist('span','var') || isempty(span)
    span=span_default;
end
if ~exist('res','var') || isempty(res)
    res=res_default;
end
if ~exist('frrate','var') || isempty(frrate)
    frrate=frr_default;
end
if ~exist('invert','var') || isempty(invert)
    invert=invert_default;
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisy_default;
end
if ~exist('bac','var') || isempty(bac)
    bac=bac_default;
end
if ~exist('maxPx','var') || isempty(maxPx)
    maxPx=maxPx_default;
end
if ~exist('maxSpeed','var') || isempty(maxSpeed)
    maxSpeed=maxSpeed_default;
else
    maxSpeed=maxSpeed*res;
end

% -=- Get info, get set up -=---------------------
gcp; % Opens a parallel pool
if ~exist('mov','var') || isempty(mov)
    error('No input data provided.');
end

ftl=fittype('logisticcdf2left(x,a,b,mu,I)'); %Custom fittypes
ftr=fittype('logisticcdf2right(x,a,b,mu,I)'); 
movsize=size(mov); %Getting array sizes.
ht=movsize(1);
wd=movsize(2);
tmin=max([framerange(1) 1]); %Frames
tmax=min([framerange(2) movsize(3)]); %Frames
t0=tmin/frrate; %Seconds
color_depth=256;

if ~exist('ROI','var') || isempty(ROI) % if no ROI specified...
    ROI=[1 1 wd ht]; % ... use whole image
end
tt=tmin:step:tmax; % frames
tts=tt/frrate; %seconds
Nf=numel(tt); % frame count
vfronts=cell(Nf,1);
limc=[inf -inf]; % for caxis

% -=- Loop through frames -=----------------------------------------------
% The plan: in each frame, locate fronts by finding regions brighter than
% thr, getting their boundaries, and smoothing. To get chemical velocity,
% consider frames in pairs. Then find chemical velocity by drawing lines
% locally perpendicular to the earlier front, through the later
% front. Use those lines to grab image profiles, and fit curves for
% thickness values.

for ii=0:(Nf-1)

    if ii>0 %ii actually records data from the previous loop iteration.
        % -=- Already have boundaries for current frame -=-
        Nbo=Nbn; % count of boundaries
        bdszold=bdsz-1; % count of points on each boundary; we'll drop repeated final point
        bdryold=bdrynew; % front boundaries
        txold=tx; % unit tangent vector at each boundary point
        tyold=ty;
        bdrypxold=bdrypx;
        for jj=1:Nbo
            bdryold{jj}=bdryold{jj}(1:bdszold(jj),:); % drop repeated final point
            bdrypxold{jj}=bdrypxold{jj}(1:bdszold(jj),:);
            txold{jj}=txold{jj}(1:bdszold(jj));
            tyold{jj}=tyold{jj}(1:bdszold(jj));
        end
        oldt=t1; % seconds
        Npo=Npn-Nbo; % total count of points on boundaries (dropping one from each)
        bold=b2;
    else
        Nbo=0;
    end

    if ii<Nf 
        % -=- Read fronts movie, find boundaries for next frame -=-
        b=mov(:,:,tt(ii+1));
        t1=tts(ii+1); %Grab the seconds into the film we are.
        b=b-bac;
        b2=b;
        %b=imfilter(b,hv,'replicate'); % Gaussian low-pass filter
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
        bdrypx=bdrynew;
        Nbn=numel(bdrynew); % boundary count
        for jj=1:Nbn
            bdrynew{jj}(:,2)=(bdrynew{jj}(:,2)+ROI(1))*res; % mm from original origin
            bdrynew{jj}(:,1)=(bdrynew{jj}(:,1)+ROI(2))*res; % mm from original origin
        end
        a=NaN(Nbn,1);
        bdsz=NaN(Nbn,1);
        for jj=1:Nbn % check for boundaries that are too short
            a(jj)=sum(lab(:)==jj);
            bdsz(jj)=size(bdrynew{jj},1);
        end
        ind = (a>=minarea) & (bdsz>=span+1);
        bdrynew=bdrynew(ind); % drop 'em
        bdrypx=bdrypx(ind);
        bdsz=bdsz(ind);
        Nbn=numel(bdrynew); % boundary count, again
        Npn=sum(bdsz); % points count
        if ii>0
            disp(['Tracking ' num2str(Nbn) ' fronts in frame ' ...
                num2str(tt(ii)) ', ' num2str(t1-t0) ' seconds (' num2str(ii) ...
                ' of ' num2str(Nf) ').'])
        end

        % -=- Smooth fronts and find perpendicular angle at each point -=-
        tx=cell(Nbn,1); % tangential vector, x component
        ty=cell(Nbn,1); % tangential vector, y component
        parfor jj=1:Nbn
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
        end % parfor jj=1:Nbn
        
        if Nbo>0 % The prior frame has fronts, so measure their velocity.
            bdry_shifted=bdrynew;

            % -=- Calculate chemical velocity, point-by-point -=-
            txold=vertcat(txold{:});
            tyold=vertcat(tyold{:});
            ends=cumsum(bdszold); % indices where bdrys end
            begins=circshift(ends,[1 0])+1; % indices where bdrys start
            begins(1)=1;
            vdi=1:velocity_downsample:Npo;
            Npds=numel(vdi); % number of points after downsampling
            vx_small=NaN(Npds,1); % chemical velocity components, from perpendicular distance
            vy_small=NaN(Npds,1);
            bdryold1_small=vertcat(bdryold{:});
            bdryold1_small=bdryold1_small(vdi,:);
            tx_small=txold(vdi);
            ty_small=tyold(vdi);
            maxDist=(maxSpeed/frrate)*step; % max distance a front can advance            
            parfor jj=1:Npds
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
                            (bdry_shifted{ll}(:,2)-bdryold1_small(jj,2)).^2 + ...
                            (bdry_shifted{ll}(:,1)-bdryold1_small(jj,1)).^2 );
                        if minsqrdist<=maxDist^2 % close enough to intersect
                            [ix,iy]=polyxpoly(perpx,perpy, ...
                                bdry_shifted{ll}(:,2),bdry_shifted{ll}(:,1)); % find intersections
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
                 end % Edge check
            end % parfor ii=1:Npds (loop over points)
            
            % -=- Calculate Front Thickness, if asked -=-
            bdrypxold_small=vertcat(bdrypxold{:});
            bdrypxold_small=bdrypxold_small(vdi,:);
            fthick=NaN(Npds,1); %Will eventually store the thickness values
            parfor jj=1:Npds
                perpxt=bdrypxold_small(jj,2)*[1 1]+maxPx.*ty_small(jj)*[1 -1];
                perpyt=bdrypxold_small(jj,1)*[1 1]-maxPx.*tx_small(jj)*[1 -1];
                if (min(perpxt)<=0 || max(perpxt)>=ROI(3) || ...
                    min(perpyt)<= 0 || max(perpyt)>=ROI(4)) %Removes points too close to edge.
                    fthick(jj)=NaN;
                else
                    c=improfile(bold,perpxt,perpyt,'bicubic'); %Grab pixel values.
                    if ~isnan(sum(c))
                        NProf=numel(c);
                        pxdist=(1/(NProf-1))*res*sqrt((perpxt(2)-perpxt(1))^2+...
                        (perpyt(2)-perpyt(1))^2); %mm distance between two indices of c. Not trivial due to improfile.
                        pxs=(1:numel(c))'-1; %x coords for eventual fit. 
                        mguess=(numel(c)-1)/2; %Making guesses for the fitter.
                        high_c=max(c); %Used in fit guess
                        low_c=min(c); %Used in fit guess
                        topleft=mean(c(floor(mguess)+1:end))<=mean(c(1:floor(mguess))); %Which of two general forms do we have? 
                        startpt=[thickguess, low_c, high_c-low_c,mguess];
                        if topleft %Pick which type of fit to use
                            f0=fit(pxs,c,ftl,'StartPoint',startpt); %do a curve fit
                        else
                            f0=fit(pxs,c,ftr,'StartPoint',startpt); %do a curve fit
                        end %if topleft
                        fthick(jj)=f0.I*pxdist; %Scale the fit to metric units
                        if abs(fthick(jj))>=maxPx*2*res
                            fthick(jj)=NaN;
                        end
                    else
                        fthick(jj)=NaN;
                    end %Check for a NaN in c. Inelegant solution
                end %if perpendicular is inbounds
            end %for jj=1:Npds
            fthick=abs(fthick);
            
            % -=- Calculate and consolidate -=-
            vx=NaN(Npo,1);
            vy=NaN(Npo,1);
            thickness=NaN(Npo,1);
            thickness(vdi)=fthick;
            vx(vdi)=vx_small;
            vy(vdi)=vy_small;

        end % if Nbo>0

        if ii>0
            % -=- Arrange as structure for output -=-
            bdrystruct=repmat( struct('len',[],'T',oldt-t0, ...
                'X',[],'Y',[],'Vx',[],'Vy',[],'Thickness',[]) ,[Nbo 1]);
            bdryold=vertcat(bdryold{:});
            for jj=1:Nbo
                bdrystruct(jj).len=bdszold(jj); % unitless (point count)
                ind=begins(jj):ends(jj);
                bdrystruct(jj).X=(bdryold(ind,2)); % mm from original origin
                bdrystruct(jj).Y=(bdryold(ind,1)); % mm from original origin
                bdrystruct(jj).Vx=vx(ind); % mm/s
                bdrystruct(jj).Vy=vy(ind); % mm/s
                bdrystruct(jj).Thickness=thickness(ind); % mm
            end
            vfronts{ii}=bdrystruct;
        end % if ii>0     
        
    elseif ii==Nf % if ii<Nf

        % -=- Handle the final frame differently -=-
        ends=cumsum(bdszold); % indices where bdrys end
        begins=circshift(ends,[1 0])+1; % indices where bdrys start
        begins(1)=1;
        vx=NaN(Npn,1);
        vy=NaN(Npn,1);
        thickness=NaN(Npn,1);
        bdrystruct=repmat( struct('len',[],'T',oldt-t0, ...
            'X',[],'Y',[],'Vx',[],'Vy',[],'Fi',[],'Thickness',[]) ,[Nbo 1]);
        bdryold=vertcat(bdryold{:});
        for jj=1:Nbo
            bdrystruct(jj).len=bdszold(jj); % unitless (point count)
            ind=begins(jj):ends(jj);
            bdrystruct(jj).X=(bdryold(ind,2)); % mm from original origin
            bdrystruct(jj).Y=(bdryold(ind,1)); % mm from original origin
            bdrystruct(jj).Vx=vx(ind); % mm/s
            bdrystruct(jj).Vy=vy(ind); % mm/s
            bdrystruct(jj).Thickness=thickness(ind); % mm
        end
        vfronts{ii}=bdrystruct;

    end % if ii<Nf
    
end % for ii=1:Nf (loop over time)

% -=- Consolidate everything and produce output -=-
vfronts=vertcat(vfronts{:}); % turn cell array of structs into single struct

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
    axes('nextplot','add','dataaspectratio',[1 1 1],'ydir','reverse', ...
        'clim',limc);
    xlabel(['thr = ' num2str(thr) ...
        ', span = ' num2str(span)  ', minarea = ' num2str(minarea) ...
        ', invert = ' num2str(invert)], 'interpreter','none');
    b=mov(:,:,tt(1));
    b=b(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
    hi=imagesc((ROI(1)+(0:ROI(3))),(ROI(2)+(0:ROI(4))),b);
    xlim(0.5+ROI(1)+[0 ROI(3)]);
    ylim(0.5+ROI(2)+[0 ROI(4)]);
    colormap(gray); % has no effect if image is color b/c RGB values override
    allT=vertcat(vfronts(:).T);
    tyme=unique(allT);
    for ii=1:Nf-1
        if ii>1
            delete(hl)
        end
        b=mov(:,:,tt(ii));
        b=b(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
        set(hi,'cdata',b);
        title(['frame ' num2str(tt(ii)) ', ' num2str(tyme(ii)) ...
            ' seconds (' num2str(ii) ' of ' num2str(Nf) ')'], ...
            'interpreter','none');
        list=find(allT==(tyme(ii)));
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
