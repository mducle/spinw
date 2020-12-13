function brille_init(obj, varargin)
% Initialises and constructs a brille grid for interpolation
% 
% ### Syntax
% 
% `brille_init(obj, varargin)`
% 
% ### Description
% 
% This function sets up a Brille BZTrellisQ grid for interpolating the magnon
% energies (eigenvalues) and either the spin-spin correlation function Sab or
% the eigenvectors of the spin Hamiltonian. It then uses `spinwave` to 
% calculate these quantities at specified points in the grid (`fill` the grid).
%
% This function is meant to be called from `spinwave` when the option 
% `use_brille` is set. 
%
% There is a check that if the grid has been filled and the state of the spinw
% object has not changed then it will not be refilled.
% 
% ### Input Arguments
%
% `k`
% : A propagation vector - if specified will override the propagation vector
%   stored in the .mag_str field of the spinW object.
% `nExt`
% : The super cell size - if specified will override the size in .mag_str.
%
% In order to generate the correct first (irreducible) _magnetic_ Brillouin zone,
% we need to use either the magnetic propagation vector or supercell size.
% Normally, the values in the .mag_str field is used but this can be overriden
% by the user if it turns out that the interpolation is incorrect.
%
% `use_primitive`
% : If set to true, brille will find the (irreducible) first Brillouin zone
%   of the primitive rather than conventional unit cell. default is true
%
% `wedge_search`
% : If set to true, brille will find the irreducible first Brillouin zone,
%   otherwise just the first BZ will be used. default is true
%
% `search_length`
% : An integer to control how-far the vertex-finding algorithm should
%   search in τ-index. The default (==1) indicates that (1̄1̄1̄), (1̄1̄0), (1̄1̄1),
%   (1̄0̄1), ..., (111) are included. 
% 
% `node_volume_fraction`
% : The fractional volume of a single tetrahedron in the grid. Smaller numbers
%   will result in more accurate interpolation at the cost of greate computing
%   time. default is 1e-5
%
% `use_vectors`
% : If true will interpolate the eigenvectors of the Hamiltonian rather than
%   the spin-spin correlation functions Sab. default is false
%
% In addition, all parameters accepted by spinwave is accepted here and will
% be passed to `spinwave` to `fill` the grid.
%
% ### Output Arguments
%
% None. The function sets the hidden .brille field of the spinW object which
% is then used by `spinwave` to interpolate.
% 
% ### See Also
% 
% [spinw.spinwave]
%
% (C) 2020 Greg Tucker and Duc Le

% TODO: 

inpForm.fname  = {'k'          'nExt'      'use_primitive' 'search_length'};
inpForm.defval = {NaN*[0;0;0]  NaN*[0;0;0] true            1              };
inpForm.size   = {[1 -1]       [1 -1]      [1 1]           [1 1]          };

inpForm.fname  = [inpForm.fname  {'wedge_search' 'node_volume_fraction'}];
inpForm.defval = [inpForm.defval {true           1e-5                  }];
inpForm.size   = [inpForm.size   {[1 1]          [1 1]                 }];

inpForm.fname  = [inpForm.fname  {'use_vectors' 'fid' 'hermit'}];
inpForm.defval = [inpForm.defval {false         -1    true}   ];
inpForm.size   = [inpForm.size   {[1 1]         [1 1] [1 1]   }];

[kwds, passthrough] = sw_readparam(inpForm, varargin{:});
passthrough = [passthrough {'hermit' kwds.hermit}];

pref = swpref;
if kwds.fid == -1
    kwds.fid = pref.fid;
end
fid = kwds.fid;

% Check if we have done this before
obj_to_hash = struct(obj);
obj_to_hash.kwds = kwds;
obj_to_hash.passthrough = passthrough;
new_hash = get_hash(obj_to_hash);
if isfield(obj.brille, 'hash') && all(obj.brille.hash == new_hash)
    return;
end

% In SpinW notation, nExt and k define separate functionalities and need not be related.
% SpinW uses nExt to define the supercell and k to define a single-k magnetic structure
% which can be represented by a single rotating coordinate frame (so not all single-k
% structures can be represented). It is possible in SpinW to have both a supercell
% and a single-k structure by specifying non-unity values for nExt and non-zero
% values for k. In this case the spins of adjacent supercells will have their
% angles rotated by the amount determined by k.

if all(isnan(kwds.k)) && all(isnan(kwds.nExt))
    % For incommensurate structure, should use the full (structural) BZ
    % rather than the reduced BZ corresponding to the magnetic k-vector
    % So here we ignore the SpinW k-vector and just use the supercell size.
    nExt = double(obj.mag_str.nExt);
    k = 1 ./ nExt;
elseif ~all(isnan(kwds.k)) && ~all(isnan(kwds.nExt))
    % User specified nExt and k
    k = kwds.k;
    nExt = kwds.nExt;
elseif ~all(isnan(kwds.k))
    % User specified k only
    warning('Using user specified k, ignoring SpinW nExt and k values');
    k = kwds.k;
    [~, nExt] = rat(k(:), 1e-5);
else
    % User specified nExt only
    warning('Using user specified nExt, ignoring SpinW nExt and k values');
    nExt = kwds.nExt;
    k = 1 ./ nExt;
end

assert(numel(nExt)==3, 'the number of unit cell extensions, nExt, must be (1,3) or (3,1)')
nExt = nExt(:);

if numel(k)==3
    % nExt and k should be compatible, with nExt a direct lattice "vector" and
    % k a reciprocal lattice vector
    [kn, kd] = rat(k(:), 1e-5);
    % the rationalized denominator of k should (normally) be nExt
    if sum(abs(kd - nExt)) > sum(abs(kd + nExt))*eps()
        warning('k=(%d/%d,%d/%d,%d/%d) and nExt=(%d,%d,%d) are not compatible',...
            kn(1), kd(1), kn(2), kd(2), kn(3), kd(3), ...
            nExt(1), nExt(2), nExt(3))
        nExt = max( kd, nExt);
        warning('Returning supercell lattice for nExt=(%d,%d,%d)',nExt(1),nExt(2),nExt(3))
    end
end

lens = obj.lattice.lat_const(:) .* nExt;
angs = obj.lattice.angle(:); % SpinW stores angles in radian
spg = obj.lattice.label;

if kwds.use_vectors
    obj.brille.Qtrans = diag(nExt);
else
    obj.brille.Qtrans = eye(3);
    kwds.node_volume_fraction = kwds.node_volume_fraction * prod(nExt);
end

% Construct the brille grids and fill it.
fprintf0(fid, 'Creating Brille grid\n');
obj.brille.bz = brille.create_bz(lens, angs, spg, ...
                                 'use_primitive', kwds.use_primitive, ...
                                 'search_length', kwds.search_length, ...
                                 'time_reversal_symmetry', false, ...
                                 'wedge_search', kwds.wedge_search);
obj.brille.grid = brille.create_grid(obj.brille.bz, ...
                                     'node_volume_fraction', kwds.node_volume_fraction, ...
                                     'complex_vectors', kwds.use_vectors, ...
                                     'complex_values', ~kwds.hermit);
hkl = permute(brille.p2m(obj.brille.grid.rlu), [2 1]);
if sum(sum(obj.brille.Qtrans - eye(3))) > 1e-6
    hkl = obj.brille.Qtrans \ hkl;
end
fprintf0(fid, 'Filling Brille grid\n');
if kwds.use_vectors
    spec = obj.spinwave(hkl, passthrough{:}, 'saveV', true, 'sortMode', false, 'optmem', 0);
    [omega, V] = parse_twin(spec);
    if (size(omega, 1) / size(V, 1)) == 3 && (size(V, 3) / size(omega, 2)) == 3
        % Incommensurate
        kmIdx = repmat(sort(repmat([1 2 3],1,size(omega, 2))),1,1);
        eigvec = permute(cat(1, V(:,:,kmIdx==1), V(:,:,kmIdx==2), V(:,:,kmIdx==3)), [3 1 2]);
    else
        eigvec = permute(V, [3 1 2]);
    end
else
    spec = obj.spinwave(hkl, passthrough{:}, 'sortMode', false, 'formfact', false, 'optmem', 0);
    [omega, Sab] = parse_twin(spec, 'Sab');
    eigvec = permute(real(Sab), [4 3 1 2]);
end
eigval = permute(real(omega), [2 1]);
eigvec_sz = size(eigvec);
eigvec_span = prod(eigvec_sz(3:end));
obj.brille.grid.fill(eigval, 1, eigvec, eigvec_span);

% Update hash so we don't refill it unecessarily
obj.brille.hash = new_hash;

end

function out = get_hash(obj)
    % Calculates a hash for an object or struct using undocumented built-ins
    % Based on DataHash (https://uk.mathworks.com/matlabcentral/fileexchange/31272-datahash)  
    Engine = java.security.MessageDigest.getInstance('MD5');
    Engine.update(getByteStreamFromArray(obj));
    out = typecast(Engine.digest, 'uint8');
end

function [omega, V] = parse_twin(spec, use_Sab)
    if iscell(spec.omega)
        % Has twins
        omega = spec.omega{1};
        if nargin > 1 && use_Sab
            V = spec.Sab{1};
        else
            V = spec.V{1};
        end
    else
        omega = spec.omega;
        if nargin > 1 && all(logical(use_Sab))
            V = spec.Sab;
        else
            V = spec.V;
        end
    end
end

