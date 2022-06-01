Thanks for taking an interest in FrontTracker. The software is freely available for redistribution, modification, and use. However with any use or modification and redistribution, include the original information (this document, and the original copyright), along with a statement about what you have modified. See the copyright.txt document for the copyright statement.

If your use of this code leads to a publication, please cite "T. D. Nevins and D. H. Kelley. Front Tracking Velocimetry in Advection-Reaction-Diffusion Systems. Under review."

Below you should find all the basic help information for the codes contained. 
TheFrontTracker.m - takes reaction data organized as time series formats, such as .tif stack, .avi, and .seq (x,y,t) and obtains velocity and thickness measurements. 
TheFrontTrackerArrays.m - takes reaction data organized as time series as a matlab array or cell, and does exactly what TheFrontTracker does. These formats require some special handling not compatible with the methods of TheFrontTracker.m.
VF.m - can be used to plot the resulting output. 
PredictiveTracker.m - is a software for tracking tracer particles in experiments and obtaining flow data.
vp.m and Velocities.m - Codes for quickly reading the output of PredictiveTracker.m with slightly different operation styles.
write_gdf.m,read_gdf.m, and index_gdf.m - helper codes for working with .gdf format. The read and write commands are called by FrontTracker and VF, but index_gdf can be used on any gdf file to create an index file for faster reading.
read_seq.m - read data from .seq files, a format used by out camera system.  
camAlign.m - Code for aligning two cameras using an image of a checkerboard taken in both cameras. Included for completeness to show our process.
ShiftRotScale.m - Uses the data from camAlign.m to adjust flow tracks to match reaction camera coordinates.
nearneighbor.m - Helper function for camAlign.m. 

All of this information can also be obtained by typing help [corresponding code] to the command line in MATLAB, but here it is for convenience. If you still have further questions or bug reports, email tnevins@ur.rochester.edu or d.h.kelley@me.rochester.edu. 

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: vfronts =
TheFrontTracker(fname,thr,[step],[minarea],[uname],[span],[ROI],[framerange],[res],[invert],[outfile],[noisy],[thick],[bac],[maxSpeed],[maxPx])
Given a movie of moving bright regions, FrontTracker locates the front at
the edge of each region and finds its local velocity. If corresponding
flow velocities are also given, both the chemical velocity of the front
(due to molecular diffusion and reaction kinetics, presumed perpendicular
to the front) and the total velocity of the front (which also includes
advection by the underlying flow) are returned. If no flow velocities are
given, the fronts are assumed to move in a stagnant fluid, and the total
velocity equals the chemical velocity. 

The movie must be saved in one of the following formats: - a series of
image files with sequential numerical filenames, - an image stack in .tif
or .gif format, - a Norpix .seq movie, or - a video format readable by
Matlab's VideoReader.m.  Specify the movie in "fname" (e.g., '0*.png' or
'stack.tif', or 'movie.avi'). To be identified as a bright region, part
of the movie must have brightness greater than "thr+bac" and area larger
than "minarea".  The "bac" is a background image of the dish before
reaction. If invert==1, the movie is first inverted (white to black). If
thick==1, the front thickness will be calculated at all points also. The
tracker will skip frames between each boundary based on the "step" value.
Whatever integer is entered is the number of frames to step. Parts of the
movie outside the region of interest "ROI" (specified as [x y width
height], in pixels) are ignored.  Frames outside "framerange" (specified
as [firstframe lastframe]) are also ignored.  Region boundaries are
smoothed using a local line fit with span width "span". Corresponding
velocity measurements, specified with filename "uname", must be in a .gdf
file of the sort output by ShiftRotScale.m. The parameter "maxSpeed"
controls how far the tracker looks for nearby fronts when determining
front speed. The input is a distance in mm representing how far you want
the tracker to look for another front if it were looking in the very next
frame (regardless of "step"). The parameter "maxPx" is the mm distance
you want to use for thickness curve fitting.
 
Results are returned in the structure "vfronts", whose fields "len", "T",
"X", "Y", "Wx", "Wy", "Vx", "Vy", and "Thick" contain the length, time,
horizontal coordinates, vertical coordinates, horizontal total
velocities, vertical total velocities, horizontal chemical velocities,
vertical chemical velocities, and front thickness of each point on each
front, respectively. In the outputs, times are measured in seconds. If
the camera resolution (in mm/pixel) is provided in "res", lengths are
measured in mm; otherwise lengths are measured in pixels. If a filename
is provided in "outfile", results are written to disk in .gdf format
readable by vf.m. Parameters are always saved in a .mat to allow easier
repeatability. If noisy~=0, the movie is repeated onscreen with overlaid
region boundaries. If noisy==2, each movie frame is also saved to disk as
an image. Requires Norpix2MATLAB.m, write_gdf.m, vp.m,
logisticcdf2left.m, and logisticcdf2right.m. See also
PredictiveTracker.m, Fronts.m, ShiftRotScale.m, and VF.m.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: vfronts =
TheFrontTrackerArrays(fname,thr,step,[minarea],[uname],[span],[ROI],[framerange],[res],[frrate],[invert],[outfile],[noisy],[thick],[bac],[maxSpeed],[maxPx])
Given a 3D matrix of moving bright regions (x,y,t), ToyFrontTracker
locates the front at the edge of each region and finds its local
velocity. If corresponding flow velocities are also given, both the
chemical velocity of the front (due to molecular diffusion and reaction
kinetics, presumed perpendicular to the front) and the total velocity of
the front (which also includes advection by the underlying flow) are
returned. If no flow velocities are given, the fronts are assumed to move
in a stagnant fluid, and the total velocity equals the chemical velocity.
Flowing reactions have not been extensively tested in this code.

Specify the "movie variable" in "toy". To be identified as a bright
region, part of the movie must have brightness greater than "thr+bac" and
area larger than "minarea".  The "bac" is a background image before
reaction. If invert==1, the movie is first inverted (white to black). If
thick==1, the front thickness will be calculated at all points also. The
tracker will skip frames between each boundary based on the "step" value,
which improves velocity measurement.  Whatever integer is entered is the
number of frames to step. Parts of the movie outside the region of
interest "ROI" (specified as [x y width height], in pixels) are ignored.
Frames outside "framerange" (specified as [firstframe lastframe]) are
also ignored.  Region boundaries are smoothed using a local line fit with
span width "span". Corresponding velocity measurements, specified with
filename "uname", must be in a .gdf file of the sort output by
ShiftRotScale.m. The parameter "maxSpeed" controls how far the tracker
looks for nearby fronts when determining front speed. The input is a
distance in pixels representing how far you want the tracker to look for
another front if it were looking in the very next frame (regardless of
"step"). The parameter "maxPx" is the pixel distance you want to use for
thickness curve fitting.
 
Results are returned in the structure "vfronts", whose fields "len", "T",
"X", "Y", "Wx", "Wy", "Vx", "Vy", and "Thick" contain the length, time,
horizontal coordinates, vertical coordinates, horizontal total
velocities, vertical total velocities, horizontal chemical velocities,
vertical chemical velocities, and front thickness of each point on each
front, respectively. In the outputs, times are measured in seconds. If
the camera resolution (in mm/pixel) is provided in "res", frame rate in
"frrate", lengths are measured in mm; otherwise lengths are measured in
pixels. If a filename is provided in "outfile", results are written to
disk in .gdf format readable by vft.m. Parameters are always saved in a
.mat to allow easier repeatability. If noisy~=0, the movie is repeated
onscreen with overlaid region boundaries. If noisy==2, each movie frame
is also saved to disk as an image. Requires vp.m, logisticcdf2right, and
logisticcdf2left.  See also PredictiveTracker.m, Fronts.m,
ShiftRotScale.m, and vf.m.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: [wx,wy,vx,vy,x,y,t,fi,th]=VF(filename,[timerange],[region],[noisy])
Working from the reaction fronts recorded in "filename", vf plots the
fronts and their velocities at times "timerange" in the area "region".
Velocities are returned in "wx" and "wy"; positions, "x" and "y"; times,
"t". The velocity of the front without fluid flow (v=w-u) is also
returned in "vx" and "vy". If noisy==0, no plot is produced. Specify
"region" in pixels, as a vector of the form [left top width height]. The
input must be a .gdf file, sorted framewise with x, y, t, wx, wy, vx, and
vy in columns 2-8, respectively. For improved speed , first run
index_gdf.m. Requires getmaxlen.m for plotting. See also vp.m and
TheFrontTracker.m. 

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: [img,t,info]=read_seq(fname,[framerange],[noisy],[BytesAlignment])
Given the Norpix *.seq movie "fname", read_seq reads the file, returning
movie contents in "img", an array of size H x W x N, where each of N
frames is H x W pixels. The times at which the frames were recorded (in
seconds) are returned in the N x 1 array "timestamp".  If a two-element
vector "framerange" is provided, frames before framerange(1) and after
framerange(2) are ignored.  Frame numbering begins at zero. If noisy~=0,
plots are produced with a pausetime equal to noisy. 
Information about the movie is returned in the struct "info", with these
fields:
  info.Version - version number of the .seq file
  info.Description - user-provided description
  info.ImageWidth - width of each image (pixels)
  info.ImageHeight - height of each image (pixels)
  info.ImageBitDepthReal - bit depth of each image
  info.ImageFormat - format of each image
  info.AllocatedFrames - number of frames in the movie
  info.Origin - 0 if not pre/post recorded
  info.FrameRate - movie frame rate (frames/second)
  info.ReferenceFrame
  info.TimeOffsetMicroseconds - offset of image timestamps (microseconds)
  info.CustomReferenceTime - custom reference time (seconds)
  info.OldestFrameIndex - Zero unless recording in a loop.
  info.BytesAlignment - See below.
  info.DatenumStart - 1st frame time, accurate to seconds, Matlab datenum
  info.MicrosecondStart - 1st frame time, microseconds

The Matlab datenum format has insufficient precision to track fractions
of a second, so info.DatenumStart and info.MicrosecondStart are stored
separately. The BytesAlignment of a .seq file specifies the positions
within the file where each frame starts. Unfortunately it cannot always
be deduced from the file header, so read_seq attempts to guess the
correct value, returning the result in "info.BytesAlignment". If images
are returned incorrectly, try setting the "BytesAlignment" input to 1024
or 8192. 

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: write_gdf(d,filename)
For writing to .gdf format; roughly equivalent to IDL program
write_gdf.pro, but without ASCII support. See also read_gdf.m.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: [n,nstarts,nsizes]=index_gdf(filename,[colnum],[savename])
Given the name of a tracks file in .gdf format, index_gdf indexes that 
file based on column "colnum", returning the unique values occuring in that
column in "n", the locations on disk at which each value first occurs in 
"nstarts", and the size on disk of the data associated with each value in
"nsizes". Unless savename==0, the results are also saved to disk as a .mat
file. The input file must be sorted by column "colnum". index_gdf was 
written for indexing framewise-sorted particle tracks by frame number 
(time); its output file makes analysis codes like VF.m run 
much faster.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage:  a = read_gdf(filename) 

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: [phi,A,E,res]=camAlign(fnames,[L],[noisy]);
Given two images from two different cameras of the same checkerboard
pattern, camAlign calculates the best rotation and shift to express data
from the second camera in the reference frame of the first camera.
Provide image filenames in the two-element cell array "fnames"; the best
rotation angle is returned in "phi", the best shift is returned in "A",
and the total squared error is returned in "E". If the checkerboard
length "L" (in mm) is provided, each image resolution (in mm/pixel) is
calculated and returned in "res". Unless noisy==0, a plot is also
produced. Requires nearneighbor.m. See also ShiftRotScale.m.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: tracksSRS=ShiftRotScale(tracks,A,phi,res,frr,[outfile])
Given particle tracks recorded by one camera, the shift and rotation
required to align it with a second camera, and the resolution (in
mm/pixel) and frame rate (in frames/s), ShiftRotScale shifts and rotates
the tracks into the reference frame of the second camera, then converts
the units to mm and s. The tracks must be contained in "tracks", which
can be either a struct array of the sort produced by PredictiveTracker.m,
or the name of a .gdf file with eight columns (track number, x, y, t, u,
v, ax ay) of the sort produced by difftracks. Provide the x and y shift
(in mm) in the two-element vector "A", the rotation (in radians) in
"phi", the resolution (in mm/pixel) in "res", and the frame rate (in
frames/s) in "frr". Alternately, provide in "frr" the name of the .seq
movie from which the tracks were made, and ShiftRotScale will use the
original timestamp from each frame (more precise). The shifted, rotated,
and scaled tracks are output as the struct array "tracksSRS". If a
filename is provided in "outfile", they are also saved to that file (in
.gdf format and sorted by time, not track number). Requires read_gdf.m,
Velocities.m, and write_gdf.m. Requires read_seq.m to pull timestamps
from .seq file. See also camAlign.m.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: [dnn,ind]=nearneighbor(x,[S])
Given a set of points "x", nearneighbor calculates the distance from each
point to its nearest neighbor, accounting for edge effects, returning
those distances in the vector "dnn". The indices of neighbor points are
also returned in "ind". Provide "x" as an Nx2 or Nx3 array, where N is
the number of points. For any particle nearer to the edge than to its
nearest neighbor, its nearest-neighbor distance and the corresponding
index are returned as NaN. The edge is determined from "S", which can be
either a structure produced by alphavol or an alpha-radius for use with
alphavol. By default, S=inf is used, producing a convex edge. Requires
alphavol by Jonas Lundgren. See also dpart.m. 

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Copyright 2018 Thomas Nevins and Douglas H. Kelley

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: [vtracks,ntracks,meanlength,rmslength] = PredictiveTracker(inputnames,threshold,max_disp,[bground_name],[minarea],[invert],[noisy])
Given a movie of particle motions, PredictiveTracker produces Lagrangian
particle tracks using a predictive three-frame best-estimate algorithm.
The movie must be saved in one of the following formats: 
- a series of image files with sequential numerical filenames,
- an image stack in .tif or .gif format,
- a Norpix .seq movie, or
- a video format readable by Matlab's VideoReader.m.
Specify the movie in "inputnames" (e.g., '0*.png' or 'stack.tif', or
'movie.avi'). To be identified as a particle, a part of the image must
have brightness that differs from the background by at least "threshold".
If invert==0, PredictiveTracker seeks particles brighter than the
background; if invert==1, PredictiveTracker seeks particles darker than
the background; and if invert==-1, PredictiveTracker seeks any sort of
contrast. The background is read from the file "bground_name"; see
BackgroundImage. If minarea==1, PredictiveTracker seeks single-pixel
particles by comparing brightness to adjacent pixels (fast and good for
small particles); otherwise PredictiveTracker seeks particles having
areas larger than "minarea" (in square pixels; this method is better for
tracking large particles). Once identified, each particle is tracked
using a kinematic prediction, and a track is broken when no particle lies
within "max_disp" pixels of the predicted location. The results are
returned in the structure "vtracks", whose fields "len", "X", "Y", "T",
"U", and "V" contain the length, horizontal coordinates, vertical
coordinates, times, horizontal velocities, and vertical velocities of
each track, respectively. If minarea~=1, "vtracks" is returned with an
additonal field, "Theta", giving the orientation of the major axis of the
particle with respect to the x-axis, in radians. The total number of
tracks is returned as "ntracks"; the mean and root-mean-square track
lengths are returned in "meanlength" and "rmslength", respectively. If
noisy~=0, the movie is repeated with overlaid velocity quivers and the
tracks are plotted. If noisy==2, each movie frame is also saved to disk
as an image.  Requires ParticleFinder.m, and also requires read_seq.m for
use with .seq movies. If publications result from your use of this code,
please cite "Using particle tracking to measure flow instabilities in an
undergraduate laboratory experiment." D. H. Kelley and N. T.  Ouellette.
Am. J. Phys. 79 (3) 267-273 (2011). 

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: [pvx,pvy,px,py,pt,ptr]=vp(filename,[timerange],[region],[noisy])
Working from the velocity tracks recorded in "filename", vp plots the
velocity at each particle in times "timerange" in the area "region".
Velocities are returned in "pvx" and "pvy"; positions, "px" and "py";
times, "pt"; track numbers, "ptr". If noisy==0, no plot is produced.
Specify "region" in pixels, as a vector of the form [left bottom width
height]. The input can either be a .gdf file (sorted framewise with x, y,
t, vx, and vy in columns 2-6, respectively) or a directory containing
.mat files created by streamprojbN. Requires getmaxlen.m for plotting.
For improved speed with .gdf files, first run index_gdf.m. See also
vframe.m, aframe.m, vrms.m, divframe.m, vortframe.m.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: [u,v,x,y,t,tr]=Velocities(vtracks,[framerange],[noisy])
Working from the velocity tracks in "vtracks", Velocities plots the
velocity of each particle in frames "framerange". Velocity values are 
returned in vectors "u" and "v"; corresponding positions, times, and 
track indices are returned in "x", "y", and "t". If noisy==0, no plot is 
produced. The input must be a structure of the form produced by 
PredictiveTracker.m. Specify "framerange" as a two-element vector: 
[starttime endtime]. This file can be downloaded from 
http://leviathan.eng.yale.edu/software.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Copyright 2011 Douglas H. Kelley and Nicholas Ouellette

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

