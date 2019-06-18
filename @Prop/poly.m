function out = poly(varargin)

out = @(x)sum([varargin{2:2:end}].*x.^[varargin{1:2:end}]);
    
end

