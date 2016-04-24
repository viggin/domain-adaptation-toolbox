function options = svmlopt(varargin)
% SVMLOPT - Generate/alter options structure for SVM light
% 
%   OPTIONS = SVMLOPT('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates an
%   optimization options structure OPTIONS in which the named parameters
%   have the specified values.  Any unspecified parameters are set to []
%   (parameters with value [] indicate to use the default value for that
%   parameter when OPTIONS is passed to SVM light. It is sufficient to type
%   only the leading characters that uniquely identify the parameter.  Case
%   is ignored for parameter names.
%   OPTIONS = SVMLOPT(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of OLDOPTS 
%   with the named parameters altered with the specified values.
%   OPTIONS = SVMLOPT(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS.  Any parameters in NEWOPTS
%   with non-empty values overwrite the corresponding old parameters in 
%   OLDOPTS. 
%   OPTIONS = SVMLOPT (with no input arguments) creates an options structure
%   OPTIONS where all the fields are set to [].
%
%   The correspondence of the SVM light options to the field in the
%   OPTIONS structure is as follows:
%   Field          SVM light  Range, description
%   'Verbosity'      -v       {0 .. 3}, default value 1
%                             Verbosity level
%   'Regression'     -z       {0, 1}, default value 0
%                             Switch between regression [1] and
%                             classification [0]
%   'C'              -c       (0, Inf), default value (avg. x*x)^-1
%                             Trade-off between error and margin
%   'TubeWidth'      -w       (0, Inf), default value 0.1
%                             Epsilon width of tube for regression
%   'CostFactor'     -j       (0, Inf), default value 1
%                             Cost-Factor by which training errors on
%                             positive examples outweight errors on
%                             negative examples
%   'Biased'         -b       {0, 1}, default value 1
%                             Use biased hyperplane x*w+b0 [1] instead of
%                             unbiased x*w0 [0]
%   'RemoveIncons'   -i       {0, 1}, default value 0
%                             Remove inconsistent training examples and
%                             retrain
%   'ComputeLOO'     -x       {0, 1}, default value 0
%                             Compute leave-one-out estimates [1]
%   'XialphaRho'     -o       )0, 2), default value 1.0
%                             Value of rho for XiAlpha-estimator and for
%                             pruning leave-one-out computation
%   'XialphaDepth'   -k       {0..100}, default value 0
%                             Search depth for extended XiAlpha-estimator 
%   'TransPosFrac'   -p       (0..1), default value ratio of
%                             positive and negative examples in the
%                             training data. Fraction of unlabeled
%                             examples to be classified into the positive
%                             class
%   'Kernel'         -t       {0..4}, default value 1
%                             Type of kernel function:
%                             0: linear
%                             1: polynomial (s a*b+c)^d
%                             2: radial basis function exp(-gamma ||a-b||^2)
%                             3: sigmoid tanh(s a*b + c)
%                             4: user defined kernel from kernel.h
%   'KernelParam'    -d, -g, -s, -r, -u
%                             Depending on the kernel, this vector
%                             contains [d] for polynomial kernel, [gamma]
%                             for RBF, [s, c] for tanh kernel, string for
%                             user-defined kernel
%   'MaximumQP'      -q       {2..}, default value 10
%                             Maximum size of QP-subproblems
%   'NewVariables'   -n       {2..}, default value is the value chosen
%                             for 'MaximumQP'. Number of new variables
%                             entering the working set in each
%                             iteration. Use smaller values to prevent
%                             zig-zagging
%   'CacheSize'      -m       (5..Inf), default value 40.
%                             Size of cache for kernel evaluations in MB
%   'EpsTermin'      -e       (0..Inf), default value 0.001
%                             Allow that error for termination criterion
%                             [y [w*x+b] - 1] < eps
%   'ShrinkIter'     -h       {5..Inf}, default value 100.
%                             Number of iterations a variable needs to be
%                             optimal before considered for shrinking
%   'ShrinkCheck'    -f       {0, 1}, default value 1
%                             Do final optimality check for variables
%                             removed by shrinking. Although this test is
%                             usually positive, there is no guarantee
%                             that the optimum was found if the test is
%                             omitted.
%   'TransLabelFile' -l       String. File to write predicted labels of
%                             unlabeled examples into after transductive
%                             learning.
%   'AlphaFile'      -a       String. Write all alphas to this file after
%                             learning (in the same order as in the
%                             training set).
%
%   See also SVML,SVM_LEARN,SVM_CLASSIFY
%

% 
% Copyright (c) by Anton Schwaighofer (2001)
% $Revision: 1.7 $ $Date: 2002/05/30 10:33:28 $
% mailto:anton.schwaighofer@gmx.net
% 
% This program is released unter the GNU General Public License.
% 

options = struct('ExecPath', [], ...
                 'Verbosity', [], ...
                 'Regression', [], ...
                 'C', [], ...
                 'TubeWidth', [], ...
                 'CostFactor', [], ...
                 'Biased', [], ...
                 'RemoveIncons', [], ...
                 'ComputeLOO', [], ...
                 'XialphaRho', [], ...
                 'XialphaDepth', [], ...
                 'TransPosFrac', [], ...
                 'Kernel', [], ...
                 'KernelParam', [], ...
                 'MaximumQP', [], ...
                 'NewVariables', [], ...
                 'CacheSize', [], ...
                 'EpsTermin', [], ...
                 'ShrinkIter', [], ...
                 'ShrinkCheck', [], ...
                 'TransLabelFile', [], ...
                 'AlphaFile', []);

numberargs = nargin;

Names = fieldnames(options);
[m,n] = size(Names);
names = lower(Names);

i = 1;
while i <= numberargs
  arg = varargin{i};
  if isstr(arg)
    break;
  end
  if ~isempty(arg)
    if ~isa(arg,'struct')
      error(sprintf('Expected argument %d to be a string parameter name or an options structure.', i));
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),Names{j,:}))
        val = getfield(arg, Names{j,:});
      else
        val = [];
      end
      if ~isempty(val)
        [valid, errmsg] = checkfield(Names{j,:},val);
        if valid
          options = setfield(options, Names{j,:},val);
        else
          error(errmsg);
        end
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
  error('Arguments must occur in name-value pairs.');
end
expectval = 0;
while i <= numberargs
  arg = varargin{i};
  if ~expectval
    if ~isstr(arg)
      error(sprintf('Expected argument %d to be a string parameter name.', i));
    end
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)
      error(sprintf('Unrecognized parameter name ''%s''.', arg));
    elseif length(j) > 1 
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
        msg = [msg '(' Names{j(1),:}];
        for k = j(2:length(j))'
          msg = [msg ', ' Names{k,:}];
        end
        msg = sprintf('%s).', msg);
        error(msg);
      end
    end
    expectval = 1;
  else           
    [valid, errmsg] = checkfield(Names{j,:}, arg);
    if valid
      options = setfield(options, Names{j,:}, arg);
    else
      error(errmsg);
    end
    expectval = 0;
  end
  i = i + 1;
end
  

function [valid, errmsg] = checkfield(field,value)
% CHECKFIELD Check validity of structure field contents.
%   [VALID, MSG] = CHECKFIELD('field',V) checks the contents of the specified
%   value V to be valid for the field 'field'.
%

valid = 1;
errmsg = '';
if isempty(value)
  return
end
isFloat = length(value==1) & isa(value, 'double');
isPositive = isFloat & (value>=0);
isString = isa(value, 'char');
range = [];
requireInt = 0;
switch field
  case {'Biased', 'ComputeLOO', 'ShrinkCheck'}
    if isFloat,
      if all(value~=[0 1]);
        valid = 0;
        errmsg = sprintf('Invalid value for %s parameter: Must be either 0 or 1', field);
      end
    end
  case 'XialphaRho'
    if ~isFloat | (value<=0) | (value>2),
      valid = 0;
      errmsg = sprintf('Invalid value for %s parameter: Must be scalar in the range (0, 2]', field);
    end
  case {'ExecPath', 'TransLabelFile', 'AlphaFile'}
    if ~isString,
      valid = 0;
      errmsg = sprintf('Invalid value for %s parameter: Must be a string', field);
    end
  case 'RemoveIncons'
    requireInt = 1;
    range = [0 3];
    % does not comply to the SVMlight help text, but is mentioned in the code
  case 'Verbosity'
    requireInt = 1;
    range = [0 3];
  case 'Regression'
    requireInt = 1;
    range = [0 1];
  case {'C', 'CostFactor', 'EpsTermin', 'TubeWidth'}
    range = [0 Inf];
  case 'Kernel'
    requireInt = 1;
    range = [0 4];
  case 'KernelParam'
    valid = 1;
  case 'XialphaDepth'
    requireInt = 1;
    range = [0 100];
  case 'TransPosFrac'
    range = [0 1];
  case 'MaximumQP'
    requireInt = 1;
    range = [2, Inf];
  case 'NewVariables'
    requireInt = 1;
    range = [2, Inf];
  case 'CacheSize'
    requireInt = 1;
    range = [5, Inf];
  case 'ShrinkIter'
    requireInt = 1;
    range = [5, Inf];
  otherwise
    valid = 0;
    error('Unknown field name for Options structure.')
end

if ~isempty(range),
  if (value<range(1)) | (value>range(2)),
    valid = 0;
    errmsg = sprintf('Invalid value for %s parameter: Must be scalar in the range [%g..%g]', ...
                     field, range(1), range(2));
  end
end

if requireInt & ((value-round(value))~=0),
  valid = 0;
  errmsg = sprintf('Invalid value for %s parameter: Must be integer', ...
                   field);
end
