function write_gdf(d,filename)
% Usage: write_gdf(d,filename)
% For writing to .gdf format; roughly equivalent to IDL program
% write_gdf.pro, but without ASCII support. See also read_gdf.m.

% Written 27 Apr 2009 by Doug Kelley.
% Updated 15 October 2009 to handle headers better (but only in two
% dimensions). 
% Updated 10 November 2009 to correct usage of "code" and header(6).

magicnum = 082991;
Ndim = 2; % two-dimensional arrays only
code = 4; % IDL flag for float32

if nargin<2
    error(['Usage: ' mfilename '(d,filename)']);
end

if ndims(d)>2
    error(['Sorry, ' mfilename ' supports only two-dimensional arrays.'])
end

[Nrow,Ncol]=size(d);
header=[magicnum Ndim Ncol Nrow code Nrow*Ncol];

fid=fopen(filename,'w');
if fid==-1
    error(['Sorry, cannot create ' filename ' for output.']);
end
fwrite(fid,header,'int32');
fwrite(fid,d','float32');
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

