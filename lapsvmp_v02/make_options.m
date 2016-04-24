function options = make_options(varargin)
% {make_options} generate a complete option data structure.
%     
%      options = make_options()
%      options = make_options('Param1','StringValue1','Param2',NumValue2)
%
%      options: a structure with a wide set of fields required by the
%               functions provided in this library. On each function you
%               can find a description of the parameters that it requires
%               and of their default values.
%               
% Author: Stefano Melacci (2009)
%         mela@dii.unisi.it

% options default values
options = struct('Kernel', 'rbf', ...
                 'KernelParam', 1, ...
                 'NN', 6, ...
                 'GraphDistanceFunction', 'euclidean', ... 
                 'GraphWeights', 'heat', ...
                 'GraphWeightParam', 0, ...
                 'LaplacianNormalize', 1, ...
                 'LaplacianDegree', 1, ...
                 'gamma_A', 1e-6, ... 
                 'gamma_I', 1.0, ...
                 'Verbose', true, ...
                 'Hinge', true, ...
                 'Cg', false, ...
                 'UseBias', false, ...
                 'NewtonLineSearch', false, ...
                 'NewtonCholesky', false, ...
                 'MaxIter', 200, ...
                 'InitAlpha', false, ...
                 'InitBias', 0, ...
                 'CgStopType', 0, ...
                 'CgStopParam', [], ...
                 'CgStopIter', [], ...
                 'UseOneClass', 0);

numberargs = nargin;
if numberargs == 0
    return
end
if rem(numberargs,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
for i = 1:2:numberargs
    if ~ischar(varargin{i})
        error('Arguments name must be strings.');
    end
    [valid, errmsg] = checkfield(varargin{i},varargin{i+1});
    if valid
        options.(varargin{i}) = varargin{i+1};
    else
        error(errmsg);
    end
end
return


function [valid,errmsg] = checkfield(field,value)
% {checkfield} checks validity of structure field contents.

valid = 1;
errmsg = '';

if isempty(value)
    return
end

isString = isa(value, 'char');
range = [];
requireInt = 0;
requireScalar = 0;
requireString = 0;

switch field
    case 'NN'
        requireInt = 1;
        requireScalar = 1;
        range = [1 Inf];        
    case 'LaplacianNormalize'
        requireInt = 1;
        requireScalar = 1;
        range = [0 1];      
    case 'LaplacianDegree'
        requireInt = 1;
        range = [1 Inf];
        requireScalar = 1;         
    case {'Kernel', 'GraphWeights', 'GraphDistanceFunction'}
        requireString = 1;
    case {'gamma_A', 'gamma_I', 'GraphWeightParam'}
        range = [0 Inf];
        requireScalar = 1;       
    case 'KernelParam'
        requireScalar = 1; 
        
    case {'Verbose', 'Hinge', 'Cg', 'UseBias', 'NewtonLineSearch', ...
          'NewtonCholesky'}
        requireInt = 1;
        requireScalar = 1;
        range = [0 1];
    case 'MaxIter'
        requireInt = 1;
        requireScalar = 1;
        range = [1 +Inf]; 
    case 'InitAlpha'        
    case 'InitBias'
        requireScalar = 1;
        range = [-Inf +Inf];   
    case 'CgStopType'
        requireInt = 1;
        requireScalar = 1;
        range = [-1 7];     
    case 'CgStopParam'
        requireScalar = 1;
        range = [0 1];
    case 'CgStopIter'
        requireInt = 1;
        requireScalar = 1;
        range = [1 +Inf];   
    case 'UseOneClass'    
        requireScalar = 1;
        range = [0 0.5];
    otherwise
        valid = 0;
        errmsg = ['Unknown field ' field ' for options structure.'];
end

if valid==1 && requireString && ~isString
    valid = 0;
    errmsg = (['Invalid value for' field ...
               ' parameter: Must be a string.']);
end

if valid==1 && requireScalar && isString
    valid = 0;
    errmsg = (['Invalid value for' field ...
               ' parameter: Must be a scalar.']);    
end

if valid==1 && ~isempty(range),
    if (value<range(1)) || (value>range(2))
        valid = 0;
        errmsg = sprintf('Invalid value for %s parameter: ', field);
        errmsg = strcat(errmsg, ...
                        sprintf('Must be in the range [%g..%g]', ...
                        range(1), range(2)));
    end
end

if valid==1 && requireInt && ((value-round(value))~=0)
    valid = 0;
    errmsg = (['Invalid value for' field ...
               ' parameter: Must be integer.']);
end

