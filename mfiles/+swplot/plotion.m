function varargout = plotion(varargin)
% plots magnetic ion properties
%
% SWPLOT.PLOTION('option1', value1, ...)
%
% hFigure = SWPLOT.PLOTION(...)
%
% The function plots selected properties of magnetic ions stored in a SpinW
% object onto an swplot figure.
%
% Input:
%
% Options:
%
% obj       SpinW object.
% range     Plotting range of the lattice parameters in lattice units,
%           dimensions are [3 2]. For example to plot the first unit cell,
%           use: [0 1;0 1;0 1]. Also the number unit cells can be given
%           along the a, b and c directions: [2 1 2], that is equivalent to
%           [0 2;0 1;0 2]. Default is the single unit cell.
% rangeunit Unit in which the range is defined. It can be the following
%           string:
%               'lu'        Lattice units (default).
%               'xyz'       Cartesian coordinate system in Angstrom units.
% mode      String, defines how the bond is plotted
%               'aniso'     Ellipsoid is plotted for single ion anisotropy.
%               'g'     	Ellipsoid is drawn for g-tensor.
% scale     Scaling factor for the size of the ellipsoid relative to the 
%           shortest bond length. Default value is 1/3.
% alpha     Transparency (alpha value) of the ellipsoid, default value is 
%           0.3.
% radius1   Minimum radius of the ellipsoid, default value is 0.08.
% lineWidth Line width in pt of the main circles surrounding the ellipsoid, 
%           if zero no circles are drawn. Default is 0.5.
% figure    Handle of the swplot figure. Default is the selected figure.
% legend    Whether to add the plot to the legend, default is true.
% color     Color of the ellipsoid:
%               'auto'      All ellipsoid gets the color of the ion.
%               'colorname' All ellipsoid will have the same given color.
%               [R G B]     RGB code of the color that fix the color of all
%                           ellipsoid.
% color2    Color of the main circles, default is [0 0 0] for black. Can be
%           either a row vector of RGB code or string of a color name.
% nMesh     Resolution of the ellipse surface mesh. Integer number that is
%           used to generate an icosahedron mesh with #mesh number of
%           additional triangulation, default value is stored in
%           swpref.getpref('nmesh')
% nPatch    Number of points on the curve for the arrows, default
%           value is stored in swpref.getpref('npatch').
% tooltip   If true, the tooltips will be shown when clicking on atoms.
%           Default is true.
% shift     Column vector with 3 elements, all atomic positions will be
%           shifted by the given value. Default value is [0;0;0].
% replace   Replace previous atom plot if true. Default is true.
%
% Output:
%
% hFigure           Handle of the swplot figure.
%
% The name of the objects that are created called 'bond'. To find the
% handles and the stored data on these objects, use e.g.
%
%   sObject = swplot.findobj(hFigure,'name','bond')
%


% default values
%fontSize0 = swpref.getpref('fontsize',[]);
nMesh0    = swpref.getpref('nmesh',[]);
nPatch0   = swpref.getpref('npatch',[]);
range0    = [0 1;0 1;0 1];

inpForm.fname  = {'range' 'legend' 'label' 'scale' 'linewidth' 'alpha'};
inpForm.defval = {range0  true     true    1/3     'fix'       0.3    };
inpForm.size   = {[-1 -2] [1 1]    [1 1]   [1 1]   [1 3]       [1 1]  };
inpForm.soft   = {false   false    false   false   false       false  };

inpForm.fname  = [inpForm.fname  {'mode' 'color' 'nmesh' 'npatch'}];
inpForm.defval = [inpForm.defval {[]     'auto'  nMesh0  nPatch0 }];
inpForm.size   = [inpForm.size   {[1 -4] [1 -5]  [1 1]   [1 1]   }];
inpForm.soft   = [inpForm.soft   {true  false   false   false    }];

inpForm.fname  = [inpForm.fname  {'figure' 'obj' 'rangeunit' 'tooltip'}];
inpForm.defval = [inpForm.defval {[]       []    'lu'        true     }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1] [1 -6]      [1 1]    }];
inpForm.soft   = [inpForm.soft   {true     true  false       false    }];

inpForm.fname  = [inpForm.fname  {'shift' 'replace' 'radius1' }];
inpForm.defval = [inpForm.defval {[0;0;0] true      0.08      }];
inpForm.size   = [inpForm.size   {[3 1]   [1 1]     [1 1]     }];
inpForm.soft   = [inpForm.soft   {false   false     false     }];

param = sw_readparam(inpForm, varargin{:});

if isempty(param.figure)
    hFigure  = swplot.activefigure('plot');
else
    hFigure = param.figure;
end

% takes care of spinw object saved/loaded in/from figure
if isempty(param.obj)
    obj = getappdata(hFigure,'obj');
else
    setappdata(hFigure,'obj',copy(param.obj));
    obj = param.obj;
    setappdata(hFigure,'base',obj.basisvector);
end

% the basis vectors in columns.
BV = obj.basisvector;

% set figure title
set(hFigure,'Name', 'SpinW: Single ion properties');

% change range, if the number of unit cells are given
if numel(param.range) == 3
    param.range = [ zeros(3,1) param.range(:)];
elseif numel(param.range) ~=6
    error('plotion:WrongInput','The given plotting range is invalid!');
end

range = param.range;

switch param.rangeunit
    case 'lu'
        rangelu = [floor(range(:,1)) ceil(range(:,2))];
    case 'xyz'
        % corners of the box
        corners = BV\[range(1,[1 2 2 1 1 2 2 1]);range(2,[1 1 2 2 1 1 2 2]);range(3,[1 1 1 1 2 2 2 2])];
        rangelu = [min(corners,[],2) max(corners,[],2)];
        rangelu = [floor(rangelu(:,1)) ceil(rangelu(:,2))];
    otherwise
        error('plotion:WrongInput','The given rangeunit string is invalid!');
end

% atom data
mAtom  = obj.matom;

% generate bonds, but don't sort the bonds on DM interactions
[SS, ~] = intmatrix(obj,'plotmode',true,'extend',false,'sortDM',false,'zeroC',param.zero,'nExt',[1 1 1]);

if isempty(SS.all)
    warning('plotion:EmptyPlot','No bonds to plot!')
end

coupling        = struct;
coupling.dl     = double(SS.all(1:3,:));
coupling.atom1  = SS.all(4,:);
coupling.atom2  = SS.all(5,:);
coupling.matidx = SS.all(15,:);
coupling.idx    = SS.all(16,:);

% matrix values
mat   = reshape(SS.all(6:14,:),3,3,[]);
% DM interaction
matDM = (mat-permute(mat,[2 1 3]))/2;
% keep symmetric part of the interactions
matSym = (mat+permute(mat,[2 1 3]))/2;
matDM  = permute(cat(2,matDM(2,3,:),matDM(3,1,:),matDM(1,2,:)),[2 3 1]);

% are there non-zero DM vectors
isDM  = any(matDM(:));
isSym = any(matSym(:));

% select bond type to plot
if isempty(param.mode)
    if isDM
        param.mode = 'arrow';
    else
        param.mode = 'cylinder';
    end
end

if isempty(param.mode2)
    if isDM
        param.mode2 = 'antisym';
    elseif isSym
        param.mode2 = 'sym';
    else
        param.mode2 = 'none';
    end
end
  

% generate the  positions of the bonds
% generate positions from rangelu the inclusive range
pos = cell(1,3);
[pos{:}] = ndgrid(rangelu(1,1):rangelu(1,2),rangelu(2,1):rangelu(2,2),rangelu(3,1):rangelu(3,2));
pos  = reshape(cat(4,pos{:}),[],3)';
pos1 = bsxfun(@plus,pos,permute(mAtom.r(:,coupling.atom1),[1 3 2]));
pos2 = bsxfun(@plus,pos,permute(mAtom.r(:,coupling.atom2)+coupling.dl,[1 3 2]));

% number of unit cells
nCell = size(pos,2);

% generate matrix indices
matidx = repmat(coupling.matidx,[nCell 1]);

% generate matrix values
matSym = reshape(repmat(permute(matSym,[1 2 4 3]),[1 1 nCell 1]),3,3,[]);
mat    = reshape(repmat(permute(mat,[1 2 4 3]),[1 1 nCell 1]),3,3,[]);
matDM  = reshape(repmat(permute(matDM,[1 3 2]),[1 nCell 1]),3,[]);

pos1  = reshape(pos1,3,[]);
pos2  = reshape(pos2,3,[]);

% cut out the bonds that are out of range
switch param.rangeunit
    case 'lu'
        % lower range<=L<= upper range
        pIdx1 = all(bsxfun(@ge,pos1,range(:,1)) & bsxfun(@le,pos1,range(:,2)),1);
        pIdx2 = all(bsxfun(@ge,pos2,range(:,1)) & bsxfun(@le,pos2,range(:,2)),1);
        pIdx  = all([pIdx1;pIdx2],1);
    case 'xyz'
        % convert to xyz
        posxyz1 = BV*pos1;
        posxyz2 = BV*pos2;
        pIdx1   = all(bsxfun(@ge,posxyz1,range(:,1)) & bsxfun(@le,posxyz1,range(:,2)),1);
        pIdx2   = all(bsxfun(@ge,posxyz2,range(:,1)) & bsxfun(@le,posxyz2,range(:,2)),1);
        pIdx    = all([pIdx1;pIdx2],1);
end

if ~any(pIdx)
    warning('plotion:EmptyPlot','There are no bonds in the plotting range!')
    return
end

pos1   = pos1(:,pIdx);
pos2   = pos2(:,pIdx);
matidx = matidx(pIdx);
matSym = matSym(:,:,pIdx);
mat    = mat(:,:,pIdx);
matDM  = matDM(:,pIdx);

% bond vector in rlu
dpos = pos2-pos1;
% length of bond in Angstrom
lxyz = sqrt(sum((BV*dpos).^2,1));

% shift positions
pos1 = bsxfun(@plus,pos1,param.shift);
pos2 = bsxfun(@plus,pos2,param.shift);

% color
if strcmp(param.color,'auto')
    color = double(obj.matrix.color(:,matidx));
else
    color = swplot.color(param.color);
end

% save original matrix values into data
%mat0   = obj.matrix.mat(:,:,matidx);
matDat = mat2cell(mat,3,3,ones(1,numel(matidx)));

% legend label
lLabel = obj.matrix.label(matidx);

% plot ellipse/arrow on top of bond
switch param.mode2
    case 'none'
        % do nothing
    case 'antisym'
        % plot arrows
        % remove zero DM vectors
        matDM2 = sum(matDM.^2,1);
        zeroDM = matDM2==0;
        matDM  = matDM(:,~zeroDM);
        posDM  = (pos1(:,~zeroDM)+pos2(:,~zeroDM))/2;
        % maximum value of DM
        maxDM = sqrt(max(matDM2));
        % scale and convert vectors to RLU
        vecDM = BV\((matDM/maxDM)*param.scale*min(lxyz));
        colDM = color(:,~zeroDM);
        
        swplot.plot('type','arrow','name','bond_DM','position',cat(3,posDM,posDM+vecDM),'text','',...
            'figure',hFigure,'legend',lLabel,'color',colDM,'R',param.radius1,...
            'tooltip',false,'replace',param.replace,'npatch',param.npatch,...
            'data',{},'label',{},'nmesh',param.nmesh,'ang',param.ang,...
            'lhead',param.lhead);
        
    case 'sym'
        % plot ellipsoid
        % remove zero ellipsoids
        rmSym = permute(sumn(abs(matSym),[1 2])==0,[1 3 2]);
        % remove Heisenberg interactions
        diagSym = zeros(size(matSym));
        diagSym(1,1,:) = matSym(1,1,:);
        diagSym(2,2,:) = matSym(1,1,:);
        diagSym(3,3,:) = matSym(1,1,:);
        rmSym = rmSym | ~any(reshape(bsxfun(@minus,matSym,diagSym),9,[]),1);
        
        matSym = matSym(:,:,~rmSym);
        posSym = (pos1(:,~rmSym)+pos2(:,~rmSym))/2;
        colSym = color(:,~rmSym);
        % calculating the main radiuses of the ellipsoid.
        [V, Rell] = eigorth(matSym,1e-5);
        % creating positive definite matrix by adding constant to all
        % eigenvalues.
        maxR  = sqrt(max(sum(Rell.^2,1)));
        Rell = ((Rell+param.radius1)/(maxR+param.radius1))*param.scale*min(lxyz);
        % V*diag(R) vectorized
        V = bsxfun(@times,V,permute(Rell,[3 1 2]));
        
        swplot.plot('type','ellipsoid','name','bond_sym','position',posSym,'text','',...
            'figure',hFigure,'legend',lLabel,'color',colSym,'T',V,...
            'tooltip',false,'replace',param.replace,'npatch',param.npatch,...
            'data',{},'label',{},'nmesh',param.nmesh);

    otherwise
        error('plotion:WrongInput','The given mode2 string is invalid!');
end


lineStyle0 = '-';
switch param.mode
    case 'cylinder'
        type0 = 'cylinder';
    case 'arrow'
        type0 = 'arrow';
        % shorten the bonds
        % radius
        switch param.radius
            case 'auto'
                radius = sw_atomdata(obj.unit_cell.label,'radius');
                aidx1 = repmat(mAtom.idx(1,coupling.atom1),[nCell 1]);
                aidx1 = aidx1(pIdx);
                aidx2 = repmat(mAtom.idx(1,coupling.atom2),[nCell 1]);
                aidx2 = aidx2(pIdx);
                rad1 = radius(aidx1)*param.radius2;
                rad2 = radius(aidx2)*param.radius2;
            case 'fix'
                rad1 = param.radius2;
                rad2 = param.radius2;
            otherwise
                error('plotmag:WrongInput','The given radius string is invalid!');
        end
        pos1 = pos1 + bsxfun(@times,rad1./lxyz,dpos);
        pos2 = pos2 - bsxfun(@times,rad2./lxyz,dpos);
    case 'line'
        type0 = 'line';
        % change linestyle based on the matrix label
        % 1 for continuous line
        % 2 for dashed line
        switch param.linestyle
            case 'auto'
                lineStyle0 = (cellfun(@(C)C(end),lLabel)=='-')+1;
            case {'--' '-'}
                lineStyle0 = param.linestyle;
            otherwise
                error('plotion:WrongInput','The given lineStyle string is illegal!')
        end
    case 'empty'
        type0 = [];
    otherwise
        error('plotion:WrongInput','The given mode string is illegal!')
end
    
switch param.linewidth
    case 'fix'
        lineWidth = param.linewidth0;
    case 'lin'
        absmat = permute(sumn(abs(mat),[1 2]),[1 3 2]);
        lineWidth = absmat/max(absmat)*param.linewidth0;
    case 'pow'
        absmat = permute(sumn(abs(mat),[1 2]),[1 3 2]);
        lineWidth = (absmat/max(absmat)).^param.widthpow*param.linewidth0;
    otherwise
        error('plotion:WrongInput','The given linewidth string is illegal!')
end

% plot bond vectors
if ~isempty(type0)
    swplot.plot('type',type0,'name','bond','position',cat(3,pos1,pos2),'text','',...
        'figure',hFigure,'legend',[],'color',color,'R',param.radius0,...
        'tooltip',false,'replace',param.replace,'lineStyle',lineStyle0,...
        'data',matDat,'label',lLabel,'nmesh',param.nmesh,'ang',param.ang,...
        'lineWidth',lineWidth,'lhead',param.lhead);
end

% legend data
lDat = getappdata(hFigure,'legend');

if param.replace
    % remove old legend entries
    lIdx = ~ismember(lDat.name,'bond');
    lDat.color = lDat.color(:,lIdx);
    lDat.type  = lDat.type(:,lIdx);
    lDat.name  = lDat.name(:,lIdx);
    lDat.text  = lDat.text(:,lIdx);
    % redraw legend
    swplot.legend('refresh',hFigure);
end

if param.legend
    % append color
    lDat.color = [lDat.color double(obj.matrix.color)/255];
    % append type
    lDat.type = [lDat.type (cellfun(@(C)C(end),obj.matrix.label) == '-')+1];
    % append name
    lDat.name = [lDat.name repmat({'bond'},1,numel(obj.matrix.label))];
    % append text
    lDat.text = [lDat.text obj.matrix.label];
    
    setappdata(hFigure,'legend',lDat);
    swplot.legend('on',hFigure);
end

if nargout > 0
    varargout{1} = hFigure;
end

if param.tooltip
    swplot.tooltip('on',hFigure);
end

end