function a = read_gdf(filename)
% Usage:  a = read_gdf(filename) 

% Updated 15 October 2009 by Doug Kelley for better error handling. 
% And again 29 September 2010. 

magic = 082991;

if nargin<1
	error(['Usage: a = ' mfilename '(filename)'])
end

  % guess first that it is not natively written on a mac
  fid = fopen(filename,'r','l');
  if fid == -1 
    error(['Sorry, cannot find file ' filename '.'])
  end
  % first check the endian (I think this is correct)
  a = fread(fid,1,'ushort');
  stat = fseek(fid,0,-1);
  if a(1) == 256  % written on the mac possibly native 
    fid = fopen(filename,'r'); 
    head = fread(fid,4,'ushort'); % indeed mac native reread it 
    if head(1)~=magic
        error(['Sorry, ' filename ' does not appear to be a .gdf file.']);
    end
    dim = head(4);                % figure out the dimensionality
    fseek(fid,0,'eof');
    f_size = (ftell(fid)-24)/4.0;
    stat = fseek(fid,0,-1);       % takes you back to the beginning of
    
    if dim == 2                   % the file
       head = fread(fid,12,'ushort');  % read the entire header
       n = head(6);
       m = f_size/n; 
       stat = fseek(fid,0,-1);            % back again to the front
       stat = fseek(fid,24,-1);
       a = fread(fid,[n,m],'single');    % finally read in the data
       a = a';
    elseif dim == 3
       head = char(fread(fid,12,'ushort'))';    % entire header
       n = head(6);                      % main dimensions
       m = head(8);
       l = head(10);                     % loop dimension
       % complete this after figuring out how to make 3d arrays
    end
    
  else       % not made on the mac

    head = fread(fid,20,'long');
    if head(1)~=magic
        error(['Sorry, ' filename ' does not appear to be a .gdf file.']);
    end
    dim = head(2);
    if dim == 2
      n = head(3);
      m = head(4);
      if n==0 || m==0
          warning(['According to its header, ' filename ' contains no data.'])
      end
      stat = fseek(fid,0,-1);
      stat = fseek(fid,24,-1);
      a = fread(fid,[n,m],'single');
      a = a';
    elseif dim == 3
      n = head(6);                      % main dimensions
      m = head(8);
      l = head(10);     
      % complete this after figuring out how to make 3d arrays
    end
  end 

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

