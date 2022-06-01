function maxlen = getmaxlen(varargin)
    maxlen = 0;
    for ii=1:size(varargin, 2)
        maxlen = max(maxlen, size(varargin{ii}, 1));
    end
end

