function [ varargout ] = colorbarlbl( varargin )
  %colorbarlbl Wrapper for colorbar that accepts a label as first argument
  %
  %
  % This script and its functions follow the coding style that can be
  % sumarized in:
  % * Variables have lower camel case
  % * Functions upper camel case
  % * Constants all upper case
  % * Spaces around operators
  %
  % Authors:  Néstor Uribe-Patarroyo
  %
  % NUP: 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA;
  % <uribepatarroyo.nestor@mgh.harvard.edu>

  % MGH Flow Measurement project (v1.0)
  %
  % Changelog:
  %
  % V1.1 (2015-01-13): Some bugfixes, allows supplying property-value pairs
  % V1.0 (2014-07-03): Initial version released
  %
  % Copyright Néstor Uribe-Patarroyo (2014)
  
  
  % Look
  % http://stackoverflow.com/questions/4895556/how-to-wrap-a-function-using-varargin-and-varargout
  % for info on why this works (weird Matlab behavior!)
  
  cmds = {'off', 'hide', 'delete'}; % commands to ignore
  % list of arguments that should not be passed to colorbar, but to xlabel
  % instead
  textArgs = {'Interpreter', 'Fontsize'};
  
  % Identify handles
  handles = false(nargin, 1);
  for k = 1:nargin
    handles(k) = isscalar(ishghandle(varargin{k})) &&...
      all(ishghandle(varargin{k})) && ~strcmpi(varargin{k - 1}, 'fontsize') &&...
      ~strcmpi(varargin{k - 1}, 'linewidth');
  end
  anyHandle = any(handles);
  varArgIn = varargin;
  if anyHandle % Remove handle
    inputHandle = varArgIn{handles};
    varArgIn(handles) = [];
  end
  nArgin = numel(varArgIn);
  
  % Check if a single command is given as only argument
  if nArgin >= 1
    iscmd = any(strcmpi(varArgIn{1}, cmds));
  end
  
  % Indentify label, if supplied
  if (mod(nArgin, 2) == 0) || iscmd % Only pairs, no label, or one of the commands
    if ~anyHandle
      [varargout{1:nargout}] = colorbar(varArgIn{:});
    else
      [varargout{1:nargout}] = colorbar(inputHandle, varArgIn{:});
    end
  else
    % Get label and remove from args
    label = varArgIn{1};
    varArgIn(1) = [];
    % First identify argument pairs that belong to xlabel
    labelArgs = cellfun(@any,...
      cellfun(@strcmpi, repmat({textArgs}, [1, nArgin - 1]), varArgIn,...
      'UniformOutput', false));
    % Because they come in pairs, extend the "true" to rightmost neighbor
    if numel(labelArgs) >= 2
      labelArgs(2:2:end) = labelArgs(1:2:end - 1);
    end
    
    % Now create colorbar with proper args
    if ~anyHandle
      [colorbarHandle] = colorbar(varArgIn{~labelArgs});
    else
      [colorbarHandle] = colorbar(inputHandle, varArgIn{~labelArgs});
    end
    % Return colorbar handle if required
    if nargout > 0
      varargout{1} = colorbarHandle;
    end
    
    % And now label with proper args
    if sum(labelArgs) == 0
      ylabel(colorbarHandle, label);
    else
      ylabel(colorbarHandle, label, varArgIn{labelArgs});
    end
  end
end

