% ------------------------------------------------------------------------
% function i = MidpointCircle(i, radius, xc, yc, value) draws a circle in
% a matrix using the integer midpoint circle algorithm.
% Does not miss or repeat pixels
%
%   Output: - i: image with circle (matrix mXn)
%
%   Input:  - i: original image (matrix mXn)
%           - radius: circle radius (scalar)
%           - xc: x coordinate of the center (scalar)
%           - yc: y coordinate of the center (scalar)
%           - values: intensity of the circle (scalar)
%
% Created by : Peter Bone
% Created : 19th March 2007

% 17th October 2014
% Modified by Riccardo Lattanzi to maintain symmetry in the case that matrix dimensions are even numbers

function i = MidpointCircle(i, radius, xc, yc, value)

if rem(size(i,1),2) == 0
    xc = int16(xc);
    yc = int16(yc);
    
    x = int16(0);
    y = int16(radius);
    d = int16(1 - radius);
    
    i(xc, yc+y) = value;
    i(xc, yc-y+1) = value;
    i(xc+y, yc) = value;
    i(xc-y+1, yc) = value;
    
    while ( x < y - 1 )
        x = x + 1;
        if ( d < 0 )
            d = d + x + x + 1;
        else
            y = y - 1;
            a = x - y + 1;
            d = d + a + a;
        end
        i( x+xc,  y+yc) = value;
        i( y+xc,  x+yc) = value;
        i( y+xc, -x+yc+1) = value;
        i( x+xc, -y+yc+1) = value;
        i(-x+xc+1, -y+yc+1) = value;
        i(-y+xc+1, -x+yc+1) = value;
        i(-y+xc+1,  x+yc) = value;
        i(-x+xc+1,  y+yc) = value;
    end
    % end
    
else
    
    xc = int16(xc);
    yc = int16(yc);
    
    x = int16(0);
    y = int16(radius);
    d = int16(1 - radius);
    
    i(xc, yc+y) = value;
    i(xc, yc-y) = value;
    i(xc+y, yc) = value;
    i(xc-y, yc) = value;
    
    while ( x < y - 1 )
        x = x + 1;
        if ( d < 0 )
            d = d + x + x + 1;
        else
            y = y - 1;
            a = x - y + 1;
            d = d + a + a;
        end
        i( x+xc,  y+yc) = value;
        i( y+xc,  x+yc) = value;
        i( y+xc, -x+yc) = value;
        i( x+xc, -y+yc) = value;
        i(-x+xc, -y+yc) = value;
        i(-y+xc, -x+yc) = value;
        i(-y+xc,  x+yc) = value;
        i(-x+xc,  y+yc) = value;
    end
end
