function [ varargout ] = imagescnan(varargin )
  %imagescnan Wrapper for imagesc which sets nan values as transparent
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
  
  % MGH Therapy monitoring project (v1.0)
  %
  % Changelog:
  %
  % V1.0 (2017-11-01): Initial version released
  %
  % Copyright Néstor Uribe-Patarroyo (2017)
  
  
  % Look
  % http://stackoverflow.com/questions/4895556/how-to-wrap-a-function-using-varargin-and-varargout
  % for info on why this works (weird Matlab behavior!)
  
  % Firs check if we gave as last - 1 and last parameters a 'nancolor' and
  % value pair
  setNanColor = false;
  nanAlpha = 0;
  nArgs = nargin;
  if nargin >= 3
    if ischar(varargin{end - 1})
      if strcmpi(varargin{end - 1}, 'nancolor')
        setNanColor = true;
        nanColor = varargin{end};
        nArgs = nArgs - 2;
      elseif strcmpi(varargin{end - 1}, 'nanalpha')
        nanAlpha = varargin{end};
        nArgs = nArgs - 2;
      end
    end
    if ischar(varargin{end - 3})
      if strcmpi(varargin{end - 3}, 'nancolor')
        setNanColor = true;
        nanColor = varargin{end - 2};
        nArgs = nArgs - 2;
      elseif strcmpi(varargin{end - 3}, 'nanalpha')
        nanAlpha = varargin{end - 2};
        nArgs = nArgs - 2;
      end
    end
  end
  
  % Indentify image
  if (nArgs == 1) || (nArgs == 2) || (nArgs == 4) % Only image given or image and CLim
    if nArgs == 2
      im = varargin{1}; % Only image and crange
      alphaMat = single(~isnan(im));
      alphaMat(~alphaMat) = nanAlpha;
      [varargout{1:nargout}] = imagesc(varargin{1},...
        'AlphaData', alphaMat, 'AlphaDataMapping', 'none', varargin{2});
    elseif nArgs == 4 % two axes, image and range
      im = varargin{3};
      alphaMat = ~isnan(im);
      alphaMat(~alphaMat) = nanAlpha;
      [varargout{1:nargout}] = imagesc(varargin{1:3},...
        'AlphaData', alphaMat, 'AlphaDataMapping', 'none', varargin{4});
    else
      im = varargin{1}; % Only image given
      alphaMat = ~isnan(im);
      alphaMat(~alphaMat) = nanAlpha;
      [varargout{1:nargout}] = imagesc(varargin{1},...
        'AlphaData', alphaMat, 'AlphaDataMapping', 'none');
    end
  else
    im = varargin{3}; % Two axes must be given, image is third
    alphaMat = ~isnan(im);
    alphaMat(~alphaMat) = nanAlpha;
    [varargout{1:nargout}] = imagesc(varargin{1:nArgs},...
      'AlphaData', alphaMat, 'AlphaDataMapping', 'none');
  end

  if setNanColor
    set(gca, 'color', nanColor)
  end
end

