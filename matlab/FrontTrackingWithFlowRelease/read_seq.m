function [img,t,info]=read_seq(fname,framerange,noisy,BytesAlignment)
% Usage: [img,t,info]=read_seq(fname,[framerange],[noisy],[BytesAlignment])
% Given the Norpix *.seq movie "fname", read_seq reads the file, returning
% movie contents in "img", an array of size H x W x N, where each of N
% frames is H x W pixels. The times at which the frames were recorded (in
% seconds) are returned in the N x 1 array "timestamp".  If a two-element
% vector "framerange" is provided, frames before framerange(1) and after
% framerange(2) are ignored.  Frame numbering begins at zero. If noisy~=0,
% plots are produced with a pausetime equal to noisy. 
%
% Information about the movie is returned in the struct "info", with these
% fields:
%   info.Version - version number of the .seq file
%   info.Description - user-provided description
%   info.ImageWidth - width of each image (pixels)
%   info.ImageHeight - height of each image (pixels)
%   info.ImageBitDepthReal - bit depth of each image
%   info.ImageFormat - format of each image
%   info.AllocatedFrames - number of frames in the movie
%   info.Origin - 0 if not pre/post recorded
%   info.FrameRate - movie frame rate (frames/second)
%   info.ReferenceFrame
%   info.TimeOffsetMicroseconds - offset of image timestamps (microseconds)
%   info.CustomReferenceTime - custom reference time (seconds)
%   info.OldestFrameIndex - Zero unless recording in a loop.
%   info.BytesAlignment - See below.
%   info.DatenumStart - 1st frame time, accurate to seconds, Matlab datenum
%   info.MicrosecondStart - 1st frame time, microseconds
% 
% The Matlab datenum format has insufficient precision to track fractions
% of a second, so info.DatenumStart and info.MicrosecondStart are stored
% separately. The BytesAlignment of a .seq file specifies the positions
% within the file where each frame starts. Unfortunately it cannot always
% be deduced from the file header, so read_seq attempts to guess the
% correct value, returning the result in "info.BytesAlignment". If images
% are returned incorrectly, try setting the "BytesAlignment" input to 1024
% or 8192. 

% Written 4 June 2015 by Douglas H. Kelley.
% Fixed info.Description 5 June 2015.
% Added info.DatenumStart and info.MicrosecondsStart 15 July 2015. 
% Updated 22 July 2015 to deal gracefully with non-integer framerange.
% Updated 1 Oct 2015 to take an input pausetime.

% Next: Read 10-bit and packed 10-bit formats. 

% -=- More about versions of the .seq format -=-
% StreamPix 5.0 to 5.19: .seq version 4, BytesAlignment value in header is
% 0, actual alignment is 512.
% StreamPix 5.19: .seq version 4 or 5, BytesAlignment value in header is
% 0, actual alignment is 512 or 8192. 
% StreamPix 6.0 to 6.3.X: .seq version 5, BytesAlignment value in header is
% 0, actual alignment is 8192.
% StreamPix 6.4: .seq version 5, BytesAlignment value in header is
% 8192, actual alignment is 8192.

% -=- Defaults -=-
framerange_default = [0 inf]; % all frames
noisy_default = false; % do not plot unless requested
pausetime_default=0.5; % seconds between frames, when plotting

% -=- File format details -=-
magicnum = 65261; % or, in hex, 0xFEED
HeaderSize = 1024; % expect 1024-byte header
DescriptionFormatOffset = 592; % bytes from beginning of file
ImageFormatCode = [0,100,101,200,300,400,500,600,610,700,800,900,905, ...
    906,112,113,114,115,123,124,125,131,132,133,134,135,136,1000,1001,1002];
ImageFormatString = {'Unknown','Monochrome','Raw Bayer','BGR','Planar', ...
    'RGB','BGRx','YUV422','YUV422','UVY422','UVY411','UVY444','PhynxRGB', ...
    'PhynxRGB','Mono MSB','Mono bayer MSB','Mono MSB Swap', ...
    'Mono Bayer MSB Swap','BGR10 packed','BGR10 packed Phoenix', ...
    'RGB10 packed Phoenix','Mono packed','Mono bayer packed', ...
    'Mono packed 8448','BGR10 v1 packed','BGR10 v2 packed','Basler', ...
    'Euresys JPEG','ISG JPEG'};
CompressionFormatCode = [0 1 2 3 4 5 6 7 8 9 16 17 18];
CompressionFormatString = {'none','JPEG','RLE','Huffman','LZ','RLE fast', ...
    'Huffman fast','LZ fast','H264','wavelet','H264 RGB','lossless JPEG', ...
    'Cayer'};
epochStart='1 Jan 1970';
updateDate='16 May 2015'; % Dirty hack for guessing BytesAlignment

% -=- Parse inputs -=-
if nargin<1
    error(['Usage: [img,t,info] = ' mfilename ...
        '(fname,[framerange],[noisy],[BytesAlignment])']);
end

if ~exist('framerange','var') || isempty(framerange)
    framerange = framerange_default;
elseif numel(framerange)==1
    framerange=framerange*[1 1];
end

if ~exist('noisy','var') || isempty(noisy)
    plotit = noisy_default;
elseif noisy>0
    plotit=true;
    pausetime=noisy;
elseif noisy==0
    plotit=false;
else
    plotit=true;
    pausetime=pausetime_default;
end

% -=- Open file for reading -=-
if ~exist(fname,'file')
    error(['Sorry, cannot find ' fname '.']);
end
fid = fopen(fname,'r');
if fid == -1 
    error(['Sorry, cannot find file ' fname '.'])
end

% -=- Read header -=-
if fread(fid,1,'long')~=magicnum
    fclose(fid);
    error(['Sorry, ' fname ' is not a Norpix .seq file.']);
end
fread(fid,12,'uint16'); % expect 'Norpix seq\n'
info=struct;
info.Version=fread(fid,1,'long');
if fread(fid,1,'long')~=HeaderSize
    fclose(fid);
    error(['Sorry, ' fname ' is not a Norpix .seq file.']);
end
mypos=ftell(fid);
fseek(fid,DescriptionFormatOffset,'bof');
DescriptionFormat=fread(fid,1,'long');
fseek(fid,mypos,'bof');
info.Description=fread(fid,512,'uint8')'; % don't read char or Matlab will go 513 bytes instead of 512
switch DescriptionFormat
    case 0 % unicode
        info.Description=info.Description(1:2:end);
        info.Description(find(info.Description==0,1):end)=[]; % drop first ascii 0 and the rest
        info.Description=char(info.Description);
    case 1 % ascii
        info.Description(find(info.Description==0,1):end)=[]; % drop first ascii 0 and the rest
        info.Description=char(info.Description);
    case 2 % data
end
info.ImageWidth=fread(fid,1,'ulong');
info.ImageHeight=fread(fid,1,'ulong');
bitstr=['uint' num2str(fread(fid,1,'ulong'),'%.0f')]; % 8, 16, 24, or 32
info.ImageBitDepthReal=fread(fid,1,'ulong');
fread(fid,1,'ulong'); % ImageSizeBytes
ind=fread(fid,1,'ulong')==ImageFormatCode;
if sum(ind)==0
    fclose(fid);
    error(['Sorry, ' fname ' contains an unsupported image format.']);
else
    info.ImageFormat=ImageFormatString{ind};
end
info.AllocatedFrames=fread(fid,1,'ulong');
info.Origin=fread(fid,1,'ulong'); % 0 if not pre/post recorded
TrueImageSize=fread(fid,1,'ulong'); % number of bytes between first pixel of subsequent frames
info.FrameRate=fread(fid,1,'double');
fread(fid,1,'long'); % DescriptionFormat, again
info.ReferenceFrame=fread(fid,1,'ulong');
if fread(fid,1,'ulong')~=0 % FixedSize. Should be 0 for uncompressed data.
    fclose(fid);
    error('Sorry, compressed .seq files are not supported.')
end
fread(fid,1,'ulong'); % Flags (reserved for Norpix use)
fread(fid,1,'long'); % BayerPattern
info.TimeOffsetMicroseconds=fread(fid,1,'long');
fread(fid,1,'long'); % ExtendedHeaderSize
ind=fread(fid,1,'long')==CompressionFormatCode;
if CompressionFormatCode(ind)~=0 % Could later support compressed formats.
    fclose(fid);
    error('Sorry, compressed .seq files are not supported.')
end
info.CustomReferenceTime=fread(fid,1,'long'); % seconds since 1 Jan 1970 (C time_t)
info.CustomReferenceTime=info.CustomReferenceTime+fread(fid,1,'ushort')/10^3; % include millisecond count
info.CustomReferenceTime=info.CustomReferenceTime+fread(fid,1,'ushort')/10^6; % include microsecond count
fread(fid,1,'ulong'); % H264GOP (Group of Picture value)
fread(fid,1,'ulong'); % H264Bitrate
fread(fid,1,'ulong'); % JPEGQualityInfo
fread(fid,1,'long'); % H264DecodeFormat
fread(fid,1,'long'); % IndexOffset: Offset of compression index data
info.OldestFrameIndex=fread(fid,1,'ulong');
tmp=fread(fid,1,'ulong');
if ~exist('BytesAlignment','var') || isempty(BytesAlignment) % not specified by user
    info.BytesAlignment=tmp;
else
    info.BytesAlignment=BytesAlignment;
end
if info.BytesAlignment==0
    if info.Version==4
        info.BytesAlignment=1024;
    elseif info.Version==5
        modDate=dir(fname);
        modDate=datenum(modDate.date); % modification date of .seq file (Matlab datenum format)
        updateDate=datenum(updateDate); % date of StreamPix update (Matlab datenum format)
        if modDate > updateDate % newer file; guess StreamPix 6.3
            info.BytesAlignment=8192;
        else
            info.BytesAlignment=1024; % older file; guess StreamPix 5.19
        end
    else
        fclose(fid);
        error('Sorry, only v4 and v5 .seq files are supported.')
    end
end

% -=- Set up to read images -=-
framerange=[ceil(max(0,framerange(1))) ...
    floor(min(info.AllocatedFrames-1,framerange(2)))];
Nf=diff(framerange)+1;
img=zeros(info.ImageHeight,info.ImageWidth,Nf,bitstr);
t=NaN(Nf,1);
fseek(fid,info.BytesAlignment,'bof'); % start of first image
fread(fid,info.ImageWidth*info.ImageHeight,bitstr); % read first image
t0=fread(fid,1,'long')-info.CustomReferenceTime; % first frame time, seconds
info.DatenumStart=t0/3600/24+datenum(epochStart);
info.MicrosecondsStart=fread(fid,1,'ushort')*10^3; % millisecond count, converted to microseconds
info.MicrosecondsStart=info.MicrosecondsStart+fread(fid,1,'ushort'); % microsecond count
t0=t0+info.MicrosecondsStart/10^6; % add microseconds
frlist=framerange(1):framerange(2);

% -=- Loop through framerange, reading images and times, perhaps plotting -=-
for ii=1:Nf
    fseek(fid,info.BytesAlignment+frlist(ii)*TrueImageSize,'bof');
    img(:,:,ii)=permute(reshape( ...
        fread(fid,info.ImageWidth*info.ImageHeight,bitstr), ...
        info.ImageWidth,info.ImageHeight,[]),[2 1 3]);
    t(ii)=fread(fid,1,'long'); % seconds
    t(ii)=t(ii)+fread(fid,1,'ushort')/10^3; % add milliseconds
    t(ii)=t(ii)+fread(fid,1,'ushort')/10^6-t0; % add microseconds, remove offset

    % -=- Plot if requested -=-
    if plotit
        if ii==1
            figure;
            hi=imagesc(squeeze(img(:,:,ii)));
            caxis([0 2^info.ImageBitDepthReal-1]); % keep color scale of raw image
            set(gca,'DataAspectRatio',[1 1 1],'xtick',[],'ytick',[]);
            set(gcf,'colormap',gray)
        else
            set(hi,'cdata',squeeze(img(:,:,ii)));
        end
        title([fname ': ' num2str(t(ii)) ' s (' num2str(ii) ' of ' ...
            num2str(Nf) ')'],'interpreter','none');
        drawnow
        pause(pausetime);
    end % if noisy


end % for ii=1:Nf
fclose(fid);

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

