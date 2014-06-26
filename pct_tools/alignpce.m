%ALIGNPCE  Converts a pointcloud to normalized egg length and standard DV
%rotation.
%
%   pc = alignpce(pc, [dv_rot])
%
%   PC is the pointcloud to align
%   dv_rot is the current DV rotation of the pointcloud (defaults to 0, 
%          which is correct for atlases)
%
% Tara Martin, 2014-03-04


function pc = alignpce(pc, dv_rot)

if nargin==1, 
  % Set default DV rotation to 0 (correct for atlases)
  dv_rot = 0;
end
  
pc = egglengthnormalize(pc);    %normalize length to 1
pc = rotation(pc,dv_rot+pi/2);  %rotate so ventral side is at 0 degrees, rather than 90 degrees.
pc = align(pc);                 %apply rotation to coordinate system
