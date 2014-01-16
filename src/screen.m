% Display Function
%   Input is any variable that is compatible with the disp function
%   Will only display if the flag variable is non-zero
function screen(message,flag)
    if(flag)
        disp(message);
    end
