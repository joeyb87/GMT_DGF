function		scaleaxis(h,sc,shft)
%
% Scales the size of an axis.
%
% function		scaleaxis(h,sc)
%
% input:
%	h: handle of the axis to scale
%	sc: scaling factor
%	shft: [x y] shift
%
% Daniel Sodickson 11/1/99

if nargin < 3,
   shft = [0 0];
end
pos = get(h,'position');
set(h,'position',[pos(1)-(sc-1)*pos(3)/2+shft(1) pos(2)-(sc-1)*pos(4)/2+shft(2) sc*pos(3:4)]);
