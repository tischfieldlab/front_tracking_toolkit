Thanks for taking an interest in FrontTracker. The software is freely available for redistribution, modification, and use. However with any use or modification and redistribution, include the original information (this document, and the original copyright), along with a statement about what you have modified. See the copyright.txt document for the copyright statement.

If your use of this code leads to a publication, please cite "T. D. Nevins and D. H. Kelley. Front tracking for quantifying advection-reaction-diffusion. Under review."

Below you should find all the basic help information for the codes contained. "TheFrontTracker.m" takes reaction data organized as a normal 3D Matlab array (x,y,t) and obtains velocity and thickness measurements. "VF.m" can be used to plot the resulting output.

All of this information can also be obtained by typing help [corresponding code] to the command line in MATLAB, but here it is for convenience. If you still have further questions or bug reports, email tnevins@ur.rochester.edu or d.h.kelley@me.rochester.edu. 

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: vfronts =
TheFrontTracker(mov,thr,step,[minarea],[span],[ROI],[framerange],[res],[frrate],[invert],[noisy],[bac],[maxSpeed],[maxPx])
Given a movie "mov" of moving bright regions, TheFrontTracker locates the
fronts at the edges of each region and tracks them over time, measuring
local velocity and thickness at many points along each front. Specify
"mov" as a three-dimensional matrix whose elements are pixel
brightnesses, whose first two dimensions vary with spatial coordinates,
and whose third dimension varies with time (x,y,t).
%
To be identified as a bright region, part of the movie must have
brightness greater than the background image by an amount "thr" and have
area larger than "minarea" (in square pixels). Optionally provide a
background image in "bac".  If invert==1, the movie is first inverted
(white to black). The tracker will skip "step" movie frames between
velocity and thickness measurements, for smoothing. 

Parts of the movie outside the region of interest "ROI" (specified as [x
y width height], in pixels) are ignored.  Frames outside "framerange"
(specified as [firstframe lastframe]) are also ignored.  Region
boundaries are smoothed using a local line fit with span width "span".
The distance "maxSpeed" (in pixels per frame) controls how far the
tracker looks for nearby fronts when determining front speed.  The
parameter "maxPx" is the pixel distance you want to use for thickness
curve fitting.

Results are returned in the structure "vfronts", whose fields "len", "T",
"X", "Y", "Vx", "Vy", and "Thickness" contain the length (point count),
time, horizontal coordinates, vertical coordinates, horizontal chemical
velocities, vertical chemical velocities, and front thickness of each
point on each front, respectively. If "frrate" (in frames/second) and
"res" (in mm/pixel) are provided, outputs are reported in units of
seconds and mm; otherwise, frames and pixels are used.  If noisy~=0, the
movie is repeated onscreen with overlaid region boundaries. If noisy==2,
each movie frame is also saved to disk as an image.

See also VF.m.

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Usage: [vx,vy,x,y,t,fi,th]=VF(vfronts,[framerange],[ROI],[noisy])
Working from the tracked fronts in "vfronts", VF plots the fronts and
their velocities in frames "framerange" in the area "ROI". Velocity
values are returned in vectors "vx" and "vy"; corresponding positions,
times, and track indices are returned in "x", "y", and "t", respectively.
Front indices are returned in "fi", and thicknesses are returned in "th".
If noisy==0, no plot is produced. The input "vfronts" must be a structure
of the form produced by TheFrontTracker.m. Specify "framerange" as a
two-element vector: [starttime endtime]. Specify the region of interest,
"ROI" as a four element vector, [xmin ymin width height]. 


Copyright 2017 Thomas Nevins and Douglas H. Kelley

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
