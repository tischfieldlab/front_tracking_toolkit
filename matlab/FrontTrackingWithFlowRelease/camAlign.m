function [phi,A,E,res]=camAlign(fnames,L,noisy)
% Usage: [phi,A,E,res]=camAlign(fnames,[L],[noisy]);
% Given two images from two different cameras of the same checkerboard
% pattern, camAlign calculates the best rotation and shift to express data
% from the second camera in the reference frame of the first camera.
% Provide image filenames in the two-element cell array "fnames"; the best
% rotation angle is returned in "phi", the best shift is returned in "A",
% and the total squared error is returned in "E". If the checkerboard
% length "L" (in mm) is provided, each image resolution (in mm/pixel) is
% calculated and returned in "res". Unless noisy==0, a plot is also
% produced. Requires nearneighbor.m. See also ShiftRotScale.m.

% Written 13 March 2015 by Doug Kelley. 
% Updated 23 March 2015 to scale by the resolution before calculating the
% best-fit rotation and shift. 
% Updated 8 April 2015 to output resolution if only one image is provided.
% Updated 22 April 2015 to deal more gracefully with single-image case. 

% Next: What if the number of points detected is not the same in both
% cameras? 

outlier_deviations=3; % ignore nearest neighbor distances more than 3 deviations from mean
noisy_default=1;

if nargin<1
    error(['Usage: [phi,A,E,res] = ' mfilename '(fnames,[L],[noisy])'])
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisy_default;
end
if ~iscell(fnames) % for single file given as string
    fnames={fnames};
end

% -=- First grab mask points from calibration images -=-
Ncam=numel(fnames);
if Ncam==1
    warning('Only one image provided. Not calculating displacement.')
end
pos=cell(Ncam,1);
res=NaN(Ncam,2); 
for ii=1:Ncam
    d=imread(fnames{ii});
    disp(['Finding checkerboard points in ' fnames{ii} '.']);
    pos{ii}=detectCheckerboardPoints(d);
    if exist('L','var') && ~isempty(L)
        dnn=nearneighbor(pos{ii});
        dnn(isnan(dnn))=[];
        dm=mean(dnn);
        dd=std(dnn);
        ind = (dnn>dm-outlier_deviations*dd) & ...
            (dnn<dm+outlier_deviations*dd);
        allres=L./dnn(ind);
        res(ii,1)=mean(allres);
        res(ii,2)=std(allres);
        pos{ii}=pos{ii}*res(ii,1); % scale to mm
    else
        res=[];
    end
end
N=size(pos{1},1);

% -=- Solve the least-squares fit -=-
% We'll minimize the squared error E = sum_i( (x_i - R*X_i + A)^2 ), where
% x_i are vectors giving point coordinates from the reference mask, X_i are
% vectors giving point coordinates from the other mask, A = [a b] is an
% unknown vector offset, and R = [ cos(phi) sin(phi) ; -sin(phi) cos(phi) ]
% is an unknown rotation matrix. The three equations solved below follow
% from dE/da = 0, dE/db = 0, and dE/dphi = 0, respectively.
if Ncam>1
    disp('Calculating shift and rotation.')
    syms phi a b
    S=solve( ...
        0 == sum(pos{1}(:,1)) - sum(pos{2}(:,1))*cos(phi) ...
            - sum(pos{2}(:,2))*sin(phi) + a*N , ...
        0 == sum(pos{1}(:,2)) + sum(pos{2}(:,1))*sin(phi) ...
            - sum(pos{2}(:,2))*cos(phi) + b*N , ...
        0 == sin(phi)*sum( pos{1}(:,1).*pos{2}(:,1) ...
            + pos{1}(:,2).*pos{2}(:,2) ...
            + a*pos{2}(:,1) + b*pos{2}(:,2) ) ...
            - cos(phi)*sum( pos{1}(:,1).*pos{2}(:,2) ...
            -pos{1}(:,2).*pos{2}(:,1) ...
            + a*pos{2}(:,2) - b*pos{2}(:,1) ) );
    A=[double(S.a) double(S.b)];
    phi=double(S.phi);

    % -=- Choose solution with smallest error -=-
    E=NaN(size(phi));
    for ii=1:numel(phi)
        E(ii)=sum( pos{1}(:,1).^2+pos{1}(:,2).^2 ...
            +A(ii,1).^2+A(ii,2).^2 ...
            +pos{2}(:,1).^2+pos{2}(:,2).^2 ...
            +2*pos{1}(:,1)*A(ii,1)+2*pos{1}(:,2)*A(ii,2)+ ...
            -2*pos{1}(:,1).*pos{2}(:,1)*cos(phi(ii)) ...
            -2*pos{1}(:,2).*pos{2}(:,2)*cos(phi(ii)) ...
            -2*pos{1}(:,1).*pos{2}(:,2)*sin(phi(ii)) ...
            +2*pos{1}(:,2).*pos{2}(:,1)*sin(phi(ii)) ...
            -2*pos{2}(:,1)*A(ii,1)*cos(phi(ii)) ...
            +2*pos{2}(:,1)*A(ii,2)*sin(phi(ii)) ...
            -2*pos{2}(:,2)*A(ii,1)*sin(phi(ii)) ...
            -2*pos{2}(:,2)*A(ii,2)*cos(phi(ii)) );
    end
    [E,ind]=min(E);
    A=A(ind,:);
    phi=phi(ind);
else
    E=NaN;
    A=[0 0];
    phi=0;
end % if Ncam>1

% -=- Plot results if requested-=-
if noisy
    figure;
    hold on
    box on
    daspect([1 1 1])
    res_string=cell(2,1);
    if exist('L','var') && ~isempty(L)
        shiftunits=' mm';
        errorunits=' mm^2';
        for ii=1:Ncam
            res_string{ii}=[': ' num2str(res(ii,1)) '\pm' ...
                num2str(res(ii,2)) ' mm/px'];
        end
    else
        shiftunits=' px';
        errorunits=' px^2';
        res_string{1}=[];
        res_string{2}=[];
    end
    plot(pos{1}(:,1),pos{1}(:,2),'.','displayname', [fnames{1} ...
        res_string{1}])
    if Ncam>1
        plot(pos{2}(:,1),pos{2}(:,2),'r.','displayname', [fnames{2} ...
            res_string{2}])
        plot(pos{2}(:,1)*cos(phi)+pos{2}(:,2)*sin(phi)-A(1), ...
            -pos{2}(:,1)*sin(phi)+pos{2}(:,2)*cos(phi)-A(2),'ro', ...
            'displayname',['\phi = ' num2str(phi*180/pi) '^\circ, A = [' ...
            num2str(A(1)) ',' num2str(A(2)) ']' shiftunits ...
            ', E = ' num2str(E(1)) errorunits])
    end % if Ncam>1
    legend toggle
    set(legend,'location','northoutside')
    axis ij
    if exist('L','var') && ~isempty(L)
        xlabel('x (mm)')
        ylabel('y (mm)')
    else
        xlabel('x (px)')
        ylabel('y (px)')
    end
end % if noisy

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
