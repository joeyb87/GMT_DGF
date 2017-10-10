function fig = slice_view(scat,mapcell,titles,lims)
	default_lims = [0,120];
	if nargin == 3
		lims = default_lims;
	end
	if isnumeric(lims)
		lims = repelem({lims},size(mapcell,1),size(mapcell,2));
	end
	lims = cellfunNU(@(x)x(:).',lims);
	
	% set gui spacing properties
	ax = struct('h',480,'w',480,'hgap',20,'pad',30);
	cb = struct('w',20,'hgap',5);
	sl = struct('w',20,'hgap',10,'rgap',5);
	sl.c = struct('side',sl.w,'tgap',5);
	vgap = struct('top',20,'bot',20,'mid',20);
	
	% instantiate and store properties in map struct
	map.dims = size(mapcell);
	if numel(map.dims) < 2
		map.dims(2) = 1;
	end
	if numel(map.dims) < 3
		map.dims(3) = 1;
	end
	mapcell = reshape(mapcell, map.dims);
	if numel(map.dims) > 3
		error('Only 3 dimensional cell arrays are handled');
	end
	
	% extract relevant info
	map.mask = ~cellfun(@isempty,mapcell);
	map.T = nnz(map.mask);
	map.idxT = find(map.mask);
	map.idx = find(map.mask(:,:,1));
	map.sub = cell(1,ndims(mapcell));
	[map.sub{:}] = ind2sub(map.dims,map.idx);
	map.lims = cellfun(@max,map.sub);
	map.N = numel(map.idx);
	map.page = 1;
	map.idxP = bsxfun(@plus,map.idx(:),prod(map.dims(1:2))*(0:map.dims(3)-1));
	
	% size figure so as to fit all components
	[M,N] = deal(map.lims(1),map.lims(2));
	dims(1) = 2*N*ax.pad+N*ax.w+2*N*ax.hgap+N*cb.w+N*cb.hgap+sl.hgap+sl.w+sl.rgap;
	dims(2) = 2*M*ax.pad+vgap.bot+M*ax.h+vgap.top+(M-1)*vgap.mid;
	
	% calculate limits
	[~,axlim] = calculate_body_limits(scat.dom.r,scat.ind);

	% create figure window
	pos = get(0,'screensize');
	pos = pos(3:4)/2-dims/2;
	fig = figure('Visible','off','Position',[pos,dims],'Name','slice_view','MenuBar','none');
	
	% ensure that all pages share the same contiguous shape
	contiguous = numel(map.idx) == max(map.idx);
	shape_match = bsxfun(@xor,map.mask(:,:,:),map.mask(:,:,1));
	if ~contiguous
		error('Pages must have contiguous elements (left to right, top to bottom)');
	end
	if any(shape_match(:))
		error('Pages must have the same two dimensions');
	end
	
	% instantiate axes and colorbars
	map.Axes = cell(map.dims);
	map.CB = cell(map.dims);
	options = {'Units','pixels','PlotBoxAspectRatioMode','manual','PlotBoxAspectRatio',[1,1,1],'NextPlot','replacechildren'};
	map.Axes(map.idxT) = arrayfunNU(@(~)axes('Parent',fig,options{:}),map.idxT);
	map.CB(map.idxT)   = cellfunNU(@(a)colorbar(a,'Units','pixels'),map.Axes(map.idxT));
	cellfun(@(a,t) title(a,t), map.Axes(map.idxT), titles(map.idxT));
	map.titles = titles;
	
	% instantiate menu items
	slider  = uicontrol('Style','slider','Tag','slider','Parent',fig);
	cbutton = uicontrol('Style','pushbutton','String','C','Tag','center','Parent',fig);
	map.data = num2cell(zeros([scat.dom.dims,map.N]),1:numel(scat.dom.dims));
	for k=1:map.N
		map.data{k}(scat.ind) = mapcell{k};
	end
	slider.Callback  = {@update_plots, map, slider, lims};
	cbutton.Callback = {@center_slider, slider, map, lims};
	
	sel.dims = scat.dom.dims;
	sel.lim = axlim;
	
	% ui menu items
	sel.F = uimenu(fig,'Label','File');
	sel.S = uimenu(sel.F,'Label','Save all','Accelerator','s');
	sel.S.Callback = {@save_all,sel,slider,map,lims};
	
	sel.P = uimenu(fig,'Label','Axis');
	sel_ob = cellfunNU(@(a)uimenu(sel.P,'Label',a,'Accelerator',a),{'X','Y','Z'});
	[sel.X,sel.Y,sel.Z] = sel_ob{:};
	sel_cb = {@set_normal_axis,sel,slider,map,lims};
	set([sel_ob{:}],'Callback',sel_cb);
		
	% set positions and callback function for resizing, and enable figure
	initialize_positions(fig,[],ax,vgap,cb,map);
	set_normal_axis(sel.Z,[],sel,slider,map,lims);
	fig.SizeChangedFcn = {@update_positions,sl,cb,map};
	op = map.idxP(:,2:end);
	set([map.Axes{op(:)}], 'Visible', 'off');
	set([map.CB{op(:)}], 'Visible', 'off');
	
	% normalize axis units for proper rescaling
	set([map.Axes{map.idxT}],'Units','normalized');
	
	% resize figure window if necessary
	pos = get(0,'screensize');
	dpos = pos(3:4)-100;
	if fig.Position(3) > dpos(1)
		fig.Position(4) = round(dpos(1)*fig.Position(4)/fig.Position(3));
		fig.Position(3) = dpos(1);
	end
	if fig.Position(4) > dpos(2)
		fig.Position(3) = round(dpos(2)*fig.Position(3)/fig.Position(4));
		fig.Position(4) = dpos(2);
	end
	pos = pos(3:4)/2-fig.Position(3:4)/2;
	fig.Position(1:2) = pos(1:2);	
	fig.Visible = 'on';
end

% one-time function
function initialize_positions(fig,~,ax,vgap,cb,map)
	D = map.dims;
	axpos = zeros([4,D]);
	cbpos = zeros([4,D]);
	ins = cell(D);
	ins(map.idxT) = cellfun(@(a)a.TightInset(:),map.Axes(map.idxT),'UniformOutput',false);
	
	hins = cell(D);
	vins = cell(D);
	hins(map.idxT) = cellfun(@(in)in(1)+in(3),ins(map.idxT),'UniformOutput',false);
	vins(map.idxT) = cellfun(@(in)in(2)+in(4),ins(map.idxT),'UniformOutput',false);
	hins = max([hins{map.idxT}]);
	vins = max([vins{map.idxT}]);
	
	% horizontal placement (l2r)
	axpos([1,3],1,1) = [ax.hgap+ax.pad,ax.w].';
	cbpos([1,3],1,1) = [sum(axpos([1,3],1,1))+hins+cb.hgap,cb.w].';
	for n=2:map.dims(2)
		axpos([1,3],1,n) = [cbpos(1,1,n-1)+cb.w+ax.hgap+2*ax.pad,ax.w].';
		cbpos([1,3],1,n) = [sum(axpos([1,3],1,n))+hins+cb.hgap,cb.w].';
	end
	% copy to other rows
	axpos([1,3],2:map.lims(1),:) = repelem(axpos([1,3],1,:),1,map.lims(1)-1,1);
	cbpos([1,3],2:map.lims(1),:) = repelem(cbpos([1,3],1,:),1,map.lims(1)-1,1);
	
	% vertical placement (b2u)
	axpos([2,4],map.lims(1),1) = [vgap.bot+ax.pad,ax.h].';
	cbpos([2,4],map.lims(1),1) = [vgap.bot+ax.pad,ax.h].';
	for m=map.lims(1)-1:-1:1
		axpos([2,4],m,1) = [sum(axpos([2,4],m+1,1))+vins+2*ax.pad+vgap.mid,ax.h].';
		cbpos([2,4],m,1) = axpos([2,4],m,1);
	end
	axpos([2,4],:,2:map.lims(2)) = repelem(axpos([2,4],:,1),1,1,map.lims(2)-1);
	cbpos([2,4],:,2:map.lims(2)) = repelem(cbpos([2,4],:,1),1,1,map.lims(2)-1);
	
	[I,J] = ndgrid(1:map.dims(1),1:map.dims(2));
	arrayfun(@(k,l) set([map.Axes{k,l,:}],'Position',axpos(:,k,l),'ActivePositionProperty','outerposition'), I, J);
	arrayfun(@(k,l) set([map.CB{k,l,:}],'Position',cbpos(:,k,l)), I, J);
end

% function that updates visible ui elements
function update_positions(fig,~,sl,cb,map)
	slider  = findobj(fig,'Tag','slider');
	cbutton = findobj(fig,'Tag','center');
	
	fpos  = fig.Position;
	axesv = [map.Axes{map.idxP(:,map.page)}];
	set(axesv,'Units','pixels');
	pos = get(axesv,'Position');
	ins = get(axesv,'TightInset');
	set(axesv,'Units','normalized');
	if ~iscell(pos)
		pos = {pos};
	end
	if ~iscell(ins)
		ins = {ins};
	end
	
	laxpos = pos{map.lims(1)};
	faxpos = pos{1};
	voff = max(0,laxpos(4)-laxpos(3))/2;
	hoff = min(0,laxpos(4)-laxpos(3))/2;

	CBpos = cellfun(@(pos,ins)[pos(1)+pos(3)+ins(3)+cb.hgap+hoff,pos(2)+voff,cb.w,min(pos(3:4))].',pos,ins,'UniformOutput',false);
	arrayfun(@(k) set(map.CB{k},'Position',CBpos{k}), 1:map.N, 'UniformOutput', false);
	slider.Position = [fpos(3)-sl.rgap-sl.w,laxpos(2)+voff,sl.w,faxpos(2)+faxpos(4)-laxpos(2)];
	cbutton.Position = [slider.Position(1),slider.Position(2)-sl.c.tgap-sl.c.side,sl.c.side,sl.c.side];
end

function save_all(~,~,sel,slider,map,lims)
	folder = uigetdir('');
	base = inputdlg('Enter base filename for image files');
	pos = slider.Value;
	fig = slider.Parent;
	switch fig.UserData.axis
		case 'x'
			sel = @(d) squeeze(d(pos,:,:)).';
		case 'y'
			sel = @(d) squeeze(d(:,pos,:)).';
		case 'z'
			sel = @(d) squeeze(d(:,:,pos)).';
	end
	comp = @real;
	p = cellfun(sel, map.data(map.idxT), 'UniformOutput', false);
	p = cellfun(comp, p, 'UniformOutput', false);
	f = cellfun(@(~)figure('Visible','off'),p,'UniformOutput',false);
	ax = cellfun(@(fi)axes('Parent',fi,'Visible','off','PlotBoxAspectRatioMode','manual','PlotBoxAspectRatio',[1,1,1]),f,'UniformOutput',false);
	im = arrayfun(@(k) imagesc(p{k},'Parent',ax{k},lims{map.idxT(k)}), 1:map.T, 'UniformOutput', false);
	cb = cellfun(@(ax0) colorbar(ax0), ax, 'UniformOutput', false);
	cellfun(@(ax0,t) title(ax0,t), ax(:), map.titles(map.idxT));
	%cellfun(@(f0,t) print('-r300', f0, strjoin([folder,filesep,base,'_',t],''), '-depsc'), f(:), map.titles(map.idxT));
	cellfun(@(f0,t) print('-r300', f0, strjoin([folder,filesep,base,'_',t],''), '-dpng'), f(:), map.titles(map.idxT));
end

function update_plots(~,~,map,slider,lims)
	fig = slider.Parent;
	slider.Value = round(slider.Value);
	pos = slider.Value;
	switch fig.UserData.axis
		case 'x'
			sel = @(d) squeeze(d(pos,:,:)).';
		case 'y'
			sel = @(d) squeeze(d(:,pos,:)).';
		case 'z'
			sel = @(d) squeeze(d(:,:,pos)).';
	end
	comp = @real;
	p = cell(map.dims);
	I = map.idxP(:,map.page);
	p(I) = cellfunNU(sel, map.data(I));
	p(I) = cellfunNU(comp, p(I));
	cellfunNU(@(d,a,l) imagesc(d,'Parent',a,l), p(I), map.Axes(I), lims(I));
end

function set_normal_axis(menu,~,sel,slider,map,lims)
	fig = slider.Parent;
	axes = {'X','Y','Z'};
	tangent_x = {'Y','X','Z'};
	tangent_y = {'Z','Z','Y'};
	normal = strcmp(axes, menu.Label);
	
	sel.(axes{normal}).Checked = 'on';
	cellfun(@(a) set(sel.(a), 'Checked', 'off'), axes(~normal));
	sel.P.Label = ['Axis: ', lower(axes{normal})];
	nlims = sel.lim.(lower(axes{normal}));
	slider.Min = nlims(1);
	slider.Max = nlims(2);
	XLim = [1,sel.dims(strcmp(tangent_x(normal),axes))];
	YLim = [1,sel.dims(strcmp(tangent_y(normal),axes))];
	fig.UserData.axis = lower(axes{normal});

	set([map.Axes{:}],'XLim', XLim);
	set([map.Axes{:}],'YLim', YLim);
	slider.SliderStep = [1/(slider.Max-slider.Min), 5/(slider.Max-slider.Min)];
	center_slider([],[],slider,map,lims);
end

function center_slider(~,~,slider,map,lims)
	slider.Value = round((slider.Min+slider.Max)/2);
	update_plots([],[],map,slider,lims);
end

function set_scale(menu,~,scale)
	switch menu.Label
		case 'Log'
			scale.log.Checked = 'on';
			scale.lin.Checked = 'off';
		case 'Linear'
			scale.log.Checked = 'off';
			scale.lin.Checked = 'on';
	end
end

function [varargout] = cellfunNU(f, C1, varargin)
	[varargout{1:nargout}] = cellfun(f,C1,varargin{:},'UniformOutput',false);
end

function [varargout] = arrayfunNU(f, A1, varargin)
	[varargout{1:nargout}] = arrayfun(f,A1,varargin{:},'UniformOutput',false);
end