function [n,nstarts,nsizes]=index_gdf(filename,colnum,savename)
% Usage: [n,nstarts,nsizes]=index_gdf(filename,[colnum],[savename])
% Given the name of a tracks file in .gdf format, index_gdf indexes that 
% file based on column "colnum", returning the unique values occuring in that
% column in "n", the locations on disk at which each value first occurs in 
% "nstarts", and the size on disk of the data associated with each value in
% "nsizes". Unless savename==0, the results are also saved to disk as a .mat
% file. The input file must be sorted by column "colnum". index_gdf was 
% written for indexing framewise-sorted particle tracks by frame number 
% (time); its output file makes analysis codes like vp.m and vortp.m run 
% much faster. 

% Written 4 January 2011 by Doug Kelley.
% Added extension check 22 June 2011. 
% Re-wrote to make compatible with non-integer times (seconds, not frames!) 
% 16 November 2011. 
% Fixed bug for last few frames of large files, 5 December 2011. 

colnum11 = 5; % if .gdf has 11 columns (3D data), guess time in column 5
colnumdefault = 4; % otherwise, guess time in column 4
magicnum = 82991; % identifier for .gdf files
suffix = '_ind'; % suffix for filename of output file; used only if savename empty

if nargin<1
    error(['Usage: [n,nstarts,nsizes] = ' mfilename ...
        '(filename,[colnum],[savename])'])
end

if ~exist('savename','var') || isempty(savename)
    savename=strrep(filename,'.gdf',suffix);
end
[pathstr,name,ext]=fileparts(filename);
if isempty(ext)
    filename=fullfile(pathstr,[name '.gdf']); % append extension if necessary
end
[pathstr,name,ext]=fileparts(savename);
if isempty(ext) || ~strcmp(ext,'.mat')
    savename=fullfile(pathstr,[name '.mat']); % append extension if necessary
end

fidin=fopen(filename);
header=fread(fidin,6,'int32');
if header(1)~=magicnum
    error(['Sorry, ' filename ' does not appear to be a .gdf file.']);
end
ncol=header(3);
if ~exist('colnum','var') || isempty(colnum)
    if ncol==11
        colnum=colnum11;
    else
        colnum=colnumdefault;
    end
end
mypos=ftell(fidin); % remember where we were...
buffer=fread(fidin,ncol,'float32')';
nmin=buffer(colnum); % first frame
while buffer(colnum)==nmin % find second frame
    if feof(fidin)
        error('Sorry, the .gdf must have more than one frame.')
    end
    buffer=fread(fidin,ncol,'float32')';
end
dn=buffer(colnum)-nmin; % frame duration
fseek(fidin,-ncol*4,'eof'); % go to last row of data
buffer=fread(fidin,ncol,'float32')';
nmax=buffer(colnum);
n=NaN(1,round((nmax-nmin)/dn)); % preallocate time vector (guessing the size)
nstarts=NaN(1,numel(n)); % addresses in the input file
fseek(fidin,mypos,'bof'); % go back to first row of data 
ii=1;
n(ii)=nmin; % first frame has time nmin
nstarts(ii)=mypos; % first frame starts here
buffer=fread(fidin,ncol,'float32')';
while ~feof(fidin)
    if buffer(colnum)>n(ii) % if the data is for a new frame...
        ii=ii+1;
        n(ii)=buffer(colnum); % ...write down the new time...
        nstarts(ii)=mypos; % ...and note our position.
    end
    mypos=ftell(fidin);
    buffer=fread(fidin,ncol,'float32')';
end
fclose(fidin);
nstarts(ii+1)=mypos; % end of file is here
nsizes=diff(nstarts); % size on disk
n(ii+1:end)=[]; % now trim to actual data range
nstarts(ii+1:end)=[];
nsizes(ii+1:end)=[];
if savename
    disp(['Saving ' savename '.']);
    save(savename,'n','nstarts','nsizes');
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

