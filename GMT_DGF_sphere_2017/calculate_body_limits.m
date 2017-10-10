function [coords, idx] = calculate_body_limits(r, idxS)
	% geom.r, geom.idxS
	[L,M,P,~] = size(r);
	[IDX{1}, IDX{2}, IDX{3}] = ind2sub([L,M,P], idxS);
	idx.x = minmax(IDX{1}.');
	idx.y = minmax(IDX{2}.');
	idx.z = minmax(IDX{3}.');
	coords.x = squeeze(r(idx.x,1,1,1));
	coords.y = squeeze(r(1,idx.y,1,2));
	coords.z = squeeze(r(1,1,idx.z,3));
end