function [Ps, Pn, As, An, ssa_results] = ssa(X, d, varargin)
%SSA Stationary Subspace Analysis
%usage 
%  [Ps, Pn, As, An, ssa_results] = ssa(X, d, <options>)
%
%input
%  <no input>     Show version information
%  X              Data in one of two possible formats:
%                  * D x n matrix with data in the columns
%                  * cell array where each X{i} is a D x n_i matrix
%                    and n_i the number of samples in epoch i
%  d              Dimensionality of stationary subspace
%  <options>      List of key-value pairs to set the following options.
%    reps                Number of restarts (the one with the lowest
%                         objective function value is returned). Default: 5
%    equal_epochs        Number of equally sized epochs. equal_epochs=0 means, that
%                         the number of epochs is chosen by a heuristic. Default: 0 (chose by heuristic)
%    use_mean            Set this to false to ignore changes in the mean
%                         (for example if your dataset ensures you that no changes
%                         in the mean occur). Default: true
%    use_covariance      Set this to false to ignore changes in the
%                         covariance matrices. Default: true
%    matrix_library      matrix library to use. Has to be 'colt' or 'jblas'.
%                         Default: 'colt'
%    random_seed         Random seed, to make runs repeatable.
%                         Default: 0 (which means: no fixed random seed)
%    quiet               Set this to true in order to suppress all output.
%                         Default: false
%    ignore_determinacy  Set this to true, if the determinacy bounds should
%                         be ignored. Default: false
%                     
%
%output
%  Ps             Projection matrix to stationary subspace (d x D)
%  Pn             Projection matrix to non-stationary subspace ((D-d) x D)
%  As             Basis of stationary subspace (D x d)
%  An             Basis of non-stationary subspace (D x (D-d))
%  ssa_results    SSA results structure as described in the manual
%
%author
%  Jan Saputra Mueller, saputra@cs.tu-berlin.de
%
%license
%  This software is distributed under the BSD license. See COPYING for
%  details.
% http://mloss.org/revision/download/851/
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

% set dynamic java path to needed libraries
basedir = fileparts(mfilename('fullpath'));
javaclasspath({[basedir filesep 'ssa.jar']});

% set default parameters
opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, ...
						'reps', 5, ...
					    'equal_epochs', 0, ...
						'use_mean', true, ...
						'use_covariance', true, ...
						'matrix_library', 'colt', ...
						'random_seed', 0,  ...
					    'quiet', false, ...
                        'ignore_determinacy', false ...
						 );

% instantiate classes
if ~opt.quiet
	cl = ssatoolbox.ConsoleLogger;
else
	cl = ssatoolbox.QuietLogger;
end
ssamain = ssatoolbox.Main(false, cl);
ssamain.toolboxMode = ssatoolbox.Main.TOOLBOX_MODE_MATLAB;

% no parameter?
if ~exist('X', 'var')
    version = ssamain.data.getClass.getPackage.getImplementationVersion;
    fprintf(['SSA Toolbox version ' char(version) '\n']);
    fprintf('This software is distributed under the BSD license. See COPYING for details.\n');
    return;
end

if ~exist('d', 'var')
    error('You have to specify the number of stationary sources\n');
end

if strcmp(opt.matrix_library, 'colt')
    if ~opt.quiet, fprintf('Using Colt library...\n'); end
    ssatoolbox.SSAMatrix.setGlobalLib(ssatoolbox.SSAMatrix.COLT);
elseif strcmp(opt.matrix_library, 'jblas')
    if ~opt.quiet
			fprintf('Using jBlas library...\n'); 
    	fprintf('If you get problems with jBlas, try the Colt library.\n');
		end
    ssatoolbox.SSAMatrix.setGlobalLib(ssatoolbox.SSAMatrix.JBLAS);
else
    error('Error: Unknown matrix library %s.\n', opt.matrix_library);
end

% detect whether to use equal or custom epochization
try
 if iscell(X)
    % create custom epochization
    if ~opt.quiet, fprintf('Custom epochization found...\n'); end
    Xwoeps = [X{:}];
    Xdm = ssatoolbox.SSAMatrix(Xwoeps);
    ssamain.data.setTimeSeries(Xdm, []);
    epdef = zeros(1, size(X, 2));
    epochs = length(X);
    min_ep_size = Inf;
    last_pos = 1;
    for i=1:epochs
        ep_size = length(X{i});
        epdef(last_pos:(last_pos+ep_size-1)) = i*ones(1, ep_size);
        last_pos = last_pos + ep_size;
        if min_ep_size > ep_size
            min_ep_size = ep_size;
        end
    end
    fakefile = java.io.File('');
    ssamain.data.setCustomEpochDefinition(epdef, epochs, min_ep_size, fakefile);
    ssamain.data.setEpochType(ssamain.data.EPOCHS_CUSTOM);
 else
    % epochize equally
    if ~opt.quiet, fprintf('No custom epochization found. Using equally sized epochs.\n'); end
    Xdm = ssatoolbox.SSAMatrix(X);
    ssamain.data.setTimeSeries(Xdm, []);
    if opt.equal_epochs == 0
        % use heuristic
        ssamain.data.setEpochType(ssamain.data.EPOCHS_EQUALLY_HEURISTIC);
    else
        ssamain.data.setNumberOfEqualSizeEpochs(opt.equal_epochs);
        ssamain.data.setEpochType(ssamain.data.EPOCHS_EQUALLY);
    end
 end
catch
 e = lasterror;
 if ~isempty(strfind(e.message, 'OutOfMemoryError'))
     printJavaHeapSpaceError;
 else
     disp(e.message);
 end
 return;
end

% set random seed
if opt.random_seed ~= 0
    ssatoolbox.SSAMatrix.setRandomSeed(opt.random_seed);
end

% set SSA parameters
ssamain.parameters.setNumberOfStationarySources(d);
ssamain.parameters.setNumberOfRestarts(opt.reps);
ssamain.parameters.setUseMean(opt.use_mean);
ssamain.parameters.setUseCovariance(opt.use_covariance);
ssamain.parameters.setIgnoreDeterminacy(opt.ignore_determinacy);

% run optimization
try
%ssaresult = ssaopt.optimize(ssaparam, ssadata);
ret = ssamain.runSSA(false);
catch
 e = lasterror;
 if ~isempty(strfind(e.message, 'OutOfMemoryError'))
     printJavaHeapSpaceError;
 else
     disp(e.message);
 end
 return;
end

if ret == false
    return;
end

% return results
Ps = ssamain.results.Ps.getArray;
Pn = ssamain.results.Pn.getArray;
As = ssamain.results.Bs.getArray;
An = ssamain.results.Bn.getArray;

% ssa_results structure as described in the manual
ssa_results = struct;
ssa_results.Ps = Ps;
ssa_results.Pn = Pn;
ssa_results.As = As;
ssa_results.An = An;
if iscell(X)
    ssa_results.s_src = Ps * Xwoeps;
    ssa_results.n_src = Pn * Xwoeps;
else
    ssa_results.s_src = Ps * X;
    ssa_results.n_src = Pn * X;
end
parameters = struct;
parameters.input_file = '';
parameters.epoch_file = '';
parameters.no_s_src = d;
parameters.no_restarts = opt.reps;
parameters.use_mean = opt.use_mean;
parameters.use_covariance = opt.use_covariance;
if iscell(X)
    parameters.eq_epochs = 0;
else
    parameters.eq_epochs = opt.equal_epochs;
end
ssa_results.parameters = parameters;
ssa_results.loss_s = ssamain.results.loss_s;
ssa_results.loss_n = ssamain.results.loss_n;
ssa_results.iterations_s = ssamain.results.iterations_s;
ssa_results.iterations_n = ssamain.results.iterations_n;

ssa_results.description = ['SSA results (' datestr(now) ')'];

function printJavaHeapSpaceError
    fprintf('ERROR: Not enough Java heap space.\n');
    fprintf('To increase the Java heap space in Matlab, have a look at this website:\n\n');
    fprintf('http://www.mathworks.com/support/solutions/en/data/1-18I2C/\n\n');
    fprintf('In case you are using Matlab 2010a or later,\nthis can be easily done using Matlab''s preferences dialog.\n');

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
    warning('set_defaults:DEFAULT_FLD',['field ''' fld ''' does not have a valid default option']);
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
