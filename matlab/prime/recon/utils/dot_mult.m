function [ out ] = dot_mult( x, y )
%DOT_MULT Summary of this function goes here
%   Detailed explanation goes here

x_size = size(x);
y_size = size(y);


if ndims(x) == ndims(y)

    div_size = x_size ./ y_size;
    
    if sum(div_size > 1)
        Y = repmat(y, div_size);

        out = x .* Y; 
    else
        
        X = repmat(x, y_size ./ x_size);

        out = X .* y; 
    end
    
else
    
    if ndims(x) > ndims(y)


        y_size = padarray(y_size, [0, ndims(x) - ndims(y)], 1, 'post');

        Y = repmat(y, x_size ./ y_size);

        out = x .* Y; 


    else

        x_size = padarray(x_size, [0, ndims(y) - ndims(x)], 1, 'post');

        X = repmat(x, y_size ./ x_size);

        out = X .* y; 

    end
end


end

