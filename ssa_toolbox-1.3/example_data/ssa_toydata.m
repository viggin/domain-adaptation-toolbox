function [ X, A, cov_sources, mean_sources ] = ssa_toydata2(n, ds, dn, varargin)
%SSA_TOYDATA Generate toy data for the SSA algorithm.
%
%usage 
%  [ X, A, cov_sources, mean_sources ] = ssa_toydata2(n, ds, dn, <options>)
%
%input
%  n          Number of epochs. 
%  ds         Number of stationary sources. 
%  dn         Number of non-stationary sources. 
%  <options>  List of key-value pairs to set the following options.
%    n_samples    Number of samples in each epoch; can be either a vector 
%                 of length n or a scalar (default: 500)
%    min_n_logvar Minimum absolute log-variance of the n-sources (default: 0) 
%    n_logvar     Mean absolute log-variance above min_n_logvar of the n-sources (default: 1)
%    corr_min     Minimum canonical correlation between s/n-sources (default: 0)
%    corr_max     Maximum canonical correlation between s/n-sources (default: 0.5) 
%    mean_nonstat Non-stationarity in the mean relative to the non-stationarity 
%                 in the covariance matrix (default: 0) 
%    rand_ndir    Randomize basis of the non-stationary sources (default: true) 
%    p_nv_larger  Probability that the variance of the n-sources is 
%                 larger than the variance of the s-sources (default: 0.5)
%    s_var        Scaling of the source variance (default: 1)
%
%output 
%  X                  Cell array of epoch datasets.
%  A                  Mixing matrix. 
%  cov_epo,           Covariance matrix and mean of each epoch. These are 
%    mean_epo         the true moments, from which X has been sampled. 
%
%example 
%  [ X, A ] = ssa_toydata(10, 2, 2); % generate data 
%  [ Ps, Pn, As, An ] = ssa(X, 2);   % apply SSA 
%  subspace_error(An, A(:,[3 4]))    % measure the error
% 
%author 
%  paul.buenau@tu-berlin.de
%
%license
%  This software is distributed under the BSD license. See COPYING for
%  details.

% Copyright (c) 2010, Jan Saputra Müller, Paul von Bünau, Frank C. Meinecke,
% Franz J. Kiraly and Klaus-Robert Müller.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
% list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice, this
% list of conditions and the following disclaimer in the documentation and/or other
%  materials provided with the distribution.
% 
% * Neither the name of the Berlin Institute of Technology (Technische Universität
% Berlin) nor the names of its contributors may be used to endorse or promote
% products derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
%  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
% SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, ...
						'min_n_logvar', 0, ...
						'n_logvar', 1, ...
						'corr_min', 0, ...
						'corr_max', 0.5, ...
						'p_nv_larger', 0.5, ...
						'rand_ndir', true, ...
						'orth_mixing', false, ...
						'mean_nonstat', 0, ...
						'n_samples', 500, ...
					  's_var', 1 ...
						 );

d = ds + dn;

% Check parameters.
assert(opt.p_nv_larger >= 0 && opt.p_nv_larger <= 1);
assert(opt.corr_min >= 0 && opt.corr_max > opt.corr_min && opt.corr_max <= 1);
assert(n > 0);
assert(ds >= 0 && dn >= 0 && d>0);
assert(length(opt.n_samples) == 1 || length(opt.n_samples) == n);

if length(opt.n_samples) == 1 
	opt.n_samples = repmat(opt.n_samples, [1 n]);
end

E = rand(dn,n);
E = opt.n_logvar * (E ./ repmat(mean(E, 1), [dn 1]));
E = E + opt.min_n_logvar

% Sample the sign of the log-variances. 
ii = find(rand(dn, n) > opt.p_nv_larger);
E(ii) = -E(ii);
E = exp(E);

% Sample the canonical correlation (max. CCA eigenvalue) for each epoch. 
corrs = opt.corr_min + (opt.corr_max-opt.corr_min)*rand(1,n);

% Determine the scaling of the cross-covariance; this controls the 
% canonical correlation between s- and n-sources (corrs).
cc_scaling = sqrt(E ./ repmat((corrs.^-2) - 1, [dn 1]));

% Set direction of correlations (randomized later, if requested). 
Bn = eye(dn);
Bs = eye(ds);

% Choose random mixing matrix. 
if opt.orth_mixing
	% Random orthonormal. 
	A = randrot(d);
else
	% With random entries, normalized. 
	A = rand(d,d) - 0.5;
	A = A*diag(sum(A.^2).^-0.5);
end

% Compute the source covariance matrices for each epoch. 
cov_sources = repmat(eye(d), [1 1 n]);
for i=1:n
	if opt.rand_ndir
		Bn = randrot(dn);
		Bs = randrot(ds);
	end
	
	% Determine the matrix C which controls the correlation between 
  % s- and n-sources. 
	C = diag(cc_scaling(1:min(ds,dn),i));
	if ds>dn, C = [ C; zeros(ds-dn, dn) ];
	else
		if dn>ds, C = [ C zeros(ds, dn-ds) ]; end
	end
	C = Bs*C*Bn';

	cov_sources(1:ds,(ds+1):end, i) = C;
	cov_sources((ds+1):end,1:ds, i) = C';
	cov_sources((ds+1):end,(ds+1):end, i) = C'*C + Bn*diag(E(:,i))*Bn';

	% Overall scaling of the variance. 
	cov_sources(:,:,i) = opt.s_var*cov_sources(:,:,i);
end

% Mix the covariance matrices. 
cov_epo = zeros(d, d, n);
for i=1:n
	cov_epo(:,:,i) = A*cov_sources(:,:,i)*A';
end

% Now compute the level of non-stationarity in the covariance matrices. 
W = sqrtm(inv(mean(cov_epo, 3)));

cov_nonstat = 0;
for i=1:n 
	[foo, eigvals] = eig(W*cov_epo(:,:,i)*W');
	cov_nonstat = cov_nonstat + sum(abs(log(diag(eigvals))));
end

% Compute the epoch-means of the n-sources. Scale the the sum of the squared 
% norms relative to sum of the absolute log-eigenvalues. 

% Scale the sum of the squared norm of the mean vectors \mu_1, ..., \mu_n 
% relative to the total non-stationarity in the covariance matrix. 
norm_n = rand(1, n);
norm_n = norm_n/sum(norm_n) * opt.mean_nonstat*cov_nonstat;

% Compute mean vectors with random entries in [-0.5, 0.5].
mean_n = rand(dn, n) - 0.5;

Q = (W*A)'*(W*A);
Q = Q(ds+1:end,ds+1:end);

for i=1:n
	mean_n(:,i) = sqrt(norm_n(i))/sqrt(mean_n(:,i)'*Q*mean_n(:,i))*mean_n(:,i);
end

mean_sources = [ zeros(ds,n); mean_n ];
mean_epo = zeros(d, n);

for i=1:n
	mean_epo(:,i) = A*mean_sources(:,i);
end

% Generate samples. 
X = cell(1, n);
for i=1:n
	X{i} = mvnrnd(mean_epo(:,i)', cov_epo(:,:,i), opt.n_samples(i))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = plot_cov(C, col) 

[V,D] = eig(C);

N = 500;
alpha = linspace(0,2*pi,N);

c = 1;

E = [ cos(alpha)*sqrt(D(1,1)); sin(alpha)*sqrt(D(2,2)) ];
XY = V*E;

h = plot(XY(1,:), XY(2,:), col);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ R ] = randrot(d)
%RANDROT        Generate random orthogonal matrix. 
%
%usage
%  [R,M] = randrot(d)

M = 100*(rand(d,d)-0.5);
M = 0.5*(M-M');
R = expm(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = propertylist2struct(varargin)
% PROPERTYLIST2STRUCT - Make options structure from parameter/value list
%
%   OPT = PROPERTYLIST2STRUCT('param1', VALUE1, 'param2', VALUE2, ...)
%   Generate a structure OPT with fields 'param1' set to value VALUE1, field
%   'param2' set to value VALUE2, and so forth.
%
%   OPT has an additional field 'isPropertyStruct' that is meant to identify
%   OPT is a structure containing options data. Only in the case of missing
%   input arguments, no such identification field is written, that is,
%   PROPERTYLIST2STRUCT() returns [].
%
%   OPT2 = PROPERTYLIST2STRUCT(OPT1, 'param', VALUE, ...) takes the options
%   structure OPT1 and adds new fields 'param' with according VALUE.
%
%   See also SET_DEFAULTS
%

if nargin==0,
  % Return an empty struct without identification tag
  opt= [];
  return;
end

if isstruct(varargin{1}) | isempty(varargin{1}),
  % First input argument is already a structure: Start with that, write
  % the additional fields
  opt= varargin{1};
  iListOffset= 1;
else
  % First argument is not a structure: Assume this is the start of the
  % parameter/value list
  opt = [];
  iListOffset = 0;
end
% Write the identification field. ID field contains a 'version number' of
% how parameters are passed.
opt.isPropertyStruct = 1;

nFields= (nargin-iListOffset)/2;
if nFields~=round(nFields),
  error('Invalid parameter/value list');
end

for ff= 1:nFields,
  fld = varargin{iListOffset+2*ff-1};
  if ~ischar(fld),
    error(sprintf('String required on position %i of the parameter/value list', ...
                  iListOffset+2*ff-1));
  end
  prp= varargin{iListOffset+2*ff};
  opt= setfield(opt, fld, prp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt, isdefault]= set_defaults(opt, varargin)
%[opt, isdefault]= set_defaults(opt, defopt)
%[opt, isdefault]= set_defaults(opt, field/value list)
%
% This functions fills in the given struct opt some new fields with
% default values, but only when these fields DO NOT exist before in opt.
% Existing fields are kept with their original values.
% There are two forms in which you can can specify the default values,
% (1) as struct, 
%   opt= set_defaults(opt, struct('color','g', 'linewidth',3));
%
% (2) as property/value list, e.g.,
%   opt= set_defaults(opt, 'color','g', 'linewidth',3);
%
% The second output argument isdefault is a struct with the same fields
% as the returned opt, where each field has a boolean value indicating
% whether or not the default value was inserted in opt for that field.
%
% The default values should be given for ALL VALID property names, i.e. the
% set of fields in 'opt' should be a subset of 'defopt' or the field/value
% list. A warning will be issued for all fields in 'opt' that are not present
% in 'defopt', thus possibly avoiding a silent setting of options that are
% not understood by the receiving functions. 
%
% $Id$
% 
% Authors: Frank Meinecke (meinecke@first.fhg.de)
%          Benjamin Blankertz (blanker@first.fhg.de)
%          Pavel Laskov (laskov@first.fhg.de)

if length(opt)>1,
  error('first argument must be a 1x1 struct');
end

% Set 'isdefault' to ones for the field already present in 'opt'
isdefault= [];
if ~isempty(opt),
  for Fld=fieldnames(opt)',
    isdefault= setfield(isdefault, Fld{1}, 0);
  end
end

% Check if we have a  field/value list
if length(varargin) > 1
  
  % If the target is a propertylist structure use propertylist2struct to
  % convert the property list to a defopt structure.
  if (ispropertystruct(opt))
    defopt = propertylist2struct(varargin{:});
      
  else  % otherwise construct defopt from scratch
    
    
    % Create a dummy defopt structure: a terrible Matlab hack to overcome
    % impossibility of incremental update of an empty structure.
    defopt = struct('matlabsucks','foo');
  
    % Check consistency of a field/value list: even number of arguments
    nArgs= length(varargin)/2;
    if nArgs~=round(nArgs) & length(varargin~=1),
      error('inconsistent field/value list');
    end
    
    % Write a temporary defopt structure
    for ii= 1:nArgs,
      defopt= setfield(defopt, varargin{ii*2-1}, varargin{ii*2});
    end
    
    % Remove the dummy field from defopt
    defopt = rmfield(defopt,'matlabsucks');
  end
  
else  
  
  % If varargin has only one element, it must be a defopt structure.
  defopt = varargin{1};
  
end
  
% Replace the missing fields in 'opt' from their 'defopt' counterparts. 
for Fld=fieldnames(defopt)',
  fld= Fld{1};
  if ~isfield(opt, fld),
    opt= setfield(opt, fld, getfield(defopt, fld));
    isdefault= setfield(isdefault, fld, 1);
  end
end

% Check if some fields in 'opt' are missing in 'defopt': possibly wrong
% options.
for Fld=fieldnames(opt)',
  fld= Fld{1};
  if ~isfield(defopt,fld)
    error('set_defaults:DEFAULT_FLD',['field ''' fld ''' does not have a valid default option']);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = ispropertystruct(opts)
% ISPROPERTYSTRUCT - Check whether a structure contains optional parameters
%
%   T = ISPROPERTYSTRUCT(OPTS)
%   returns 1 if OPTS is a structure generated by PROPERTYLIST2STRUCT.
%   
%   
%   See also PROPERTYLIST2STRUCT
%

error(nargchk(1, 1, nargin));
% Currently, we do not check the version number. Existence of the field
% is enough to identify the opts structure as a property list
t = isfield(opts, 'isPropertyStruct');
