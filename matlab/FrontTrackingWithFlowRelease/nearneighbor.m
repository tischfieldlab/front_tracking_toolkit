function [dnn,ind]=nearneighbor(x,S)
% Usage: [dnn,ind]=nearneighbor(x,[S])
% Given a set of points "x", nearneighbor calculates the distance from each
% point to its nearest neighbor, accounting for edge effects, returning
% those distances in the vector "dnn". The indices of neighbor points are
% also returned in "ind". Provide "x" as an Nx2 or Nx3 array, where N is
% the number of points. For any particle nearer to the edge than to its
% nearest neighbor, its nearest-neighbor distance and the corresponding
% index are returned as NaN. The edge is determined from "S", which can be
% either a structure produced by alphavol or an alpha-radius for use with
% alphavol. By default, S=inf is used, producing a convex edge. Requires
% alphavol by Jonas Lundgren. See also dpart.m. 

% Written 20 July 2012 by Doug Kelley. 

S_default=inf; % by default, use convex boundary

if nargin<1
    error(['Usage: [dnn,ind] = ' mfilename '(x,[S])']);
end
if ~exist('S','var') || isempty(S)
    S=S_default;
end

[Np,Ndim]=size(x);
if ~any(Ndim==[2 3])
    error(['Sorry, nearneighbor supports only two- and ' ...
        'three-dimensional data.'])
end
if ~isstruct(S)
    [~,S]=alphavol(x,S); % input was alpha-radius, so compute struct
end

dist=zeros(Np);
for ii=1:Ndim
    xx=repmat(x(:,ii),1,Np);
    dist=dist+(xx-xx').^2;
end
dist=sqrt(dist); % all interparticle distances
dist(dist==0)=NaN; % ignore self-distances
[dnn,ind]=min(dist); % distances to, and indices of, nearest neighbors
if Ndim==3
    % For a plane given by ax+by+cz+d=0 and a point (x0,y0,z0), the
    % distance from the point to the plane is 
    % | ax0+by0+cz0-d | / (a^2+b^2+c^2)^(1/2). 
    warnstate=warning('off','MATLAB:TriRep:PtsNotInTriWarnId'); % keep it quiet
    TR=TriRep(S.bnd,x); % triangulation of edge facets
    warning(warnstate); % put things back the way we found them
    fn=faceNormals(TR); % normals to edge facets: a, b, and c for a plane given by ax+by+cz+d=0
    d=-sum(fn.*x(S.bnd(:,1),:),2); % d, for a plane given by ax+by+cz+d=0
    edgedist=abs(fn*x'+repmat(d,1,Np))./repmat(sqrt(sum(fn.^2,2)),1,Np); % distances from all particles to all edges
    edgedist=min(edgedist,[],1); % distance from each particle to closest edge
else
    % In two dimensions, the distance between a point (x0,y0) and the line
    % through points (x1,y1) and (x2,y2) is
    % | (x2-x1)(y2-y0) - (x1-x0)(y2-y1) | / ( (x2-x1)^2 + (y2-y1)^2 )^(1/2)
    Ne=size(S.bnd,1);
    x2x1 = repmat(x(S.bnd(:,2),1)-x(S.bnd(:,1),1),1,Np);
    y2y0 = ( repmat(x(S.bnd(:,2),2),1,Np)' - repmat(x(:,2),1,Ne) )';
    x1x0 = repmat(x(S.bnd(:,2),2)-x(S.bnd(:,1),2),1,Np);
    y2y1 = ( repmat(x(S.bnd(:,1),1),1,Np)' - repmat(x(:,1),1,Ne) )';
    denom=repmat(sqrt( sum((x(S.bnd(:,2),:)-x(S.bnd(:,1),:)).^2,2) ),1,Np);
    edgedist = abs( x2x1 .* y2y0 - x1x0 .* y2y1 ) ./ denom; % distances from all particles to all edges
    edgedist=min(edgedist,[],1); % distance from each particle to closest edge
end
pastedge=dnn>edgedist;
dnn(pastedge)=NaN; % neighbor distances further than the edge are invalid
ind(pastedge)=NaN; % nearest neighbors further than the edge are invalid


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
