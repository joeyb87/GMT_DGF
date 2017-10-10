function	plotcarr_sphere(vec,nc,nseg,lw,ls,lc)

% Plots coil arrays with geometry specified in array vec on current axis.
% Used by coilsim.m.
%
% function	plotcarr(vec,nc,nseg,lw,ls,lc)
% 
% input:
%	vec: nc*nseg x 3 array of component coil starting and ending positions
%	nc: number of coils in array
%	nseg: 1 x nc vector containing number of line segments in each component coil
%		(if nseg is input as a scalar, it is converted to a 1 x nc vector)
%	lw: linewidth for line segments
%	ls: linestyle for line segments
%	lc: linecolor for line segments
%
% output:
%	hlset: handles to lines
%
% Michael Ohliger 5/2000
% nseg is now treated so that any sized coil can be plotted
%
% Daniel Sodickson	7/22/96
% vector nseg, lw, ls, and lc added 12/18/97
% vector lw,ls,lc input added 12/23/97

if nargin < 4,
	lw = 0.5;
	ls = '-';
	lc = 'black';
elseif nargin < 5,
	ls = '-';
	lc = 'black';
elseif nargin < 6,
	lc = 'black';
end
if length(nseg) == 1,
	nseg = nseg*ones(1,nc);
end
if length(lw) == 1,
	lw = repmat(lw,[nc 1]);
end
if length(ls(:,1)) == 1,
	ls = repmat(ls,[nc 1]);
end
if length(lc(:,1)) == 1,
	lc = repmat(lc,[nc 1]);
end
hlset = [];

nseg;

for ic = 1:nc,
	closevec = [vec(sum(nseg(1:ic-1))+1:sum(nseg(1:ic)),[1 2 3]); vec(sum(nseg(1:ic-1))+1,[1 2 3])];	% DKS/MAO 5/00
	%for iseg = 1:nseg(ic),
% 		hl = line(closevec(:,1),closevec(:,2),closevec(:,3),...
% 			'linewidth',lw(ic,:),'linestyle',ls(ic,:),'color',lc(ic,:),'erasemode','normal');
		hl = line(closevec(:,1),closevec(:,2),closevec(:,3),...
			'linewidth',lw(ic,:),'linestyle',ls(ic,:),'color',lc(ic,:));
		hlset = [hlset; hl];
        %end
end
