function feat = Feature(x1, x2)
[row1, col1] = FromElement(x1);
[row2, col2] = FromElement(x2);

feat = [row1+row2 row1^2+row2^2 row1*row2 row1^2*row2^2 ...
    col1+col2 col1^2+col2^2 col1*col2 col1^2*col2^2, ...
    row1*x1+row2*x2, ...
    sqrt(x1)+sqrt(x2) sqrt(x1*x2) x1+x2 x1^2+x2^2 x1*x2, ...
    1];
end

function [row1, col1] = FromElement(element)
if(element<=2)
    row1 = 1;
    col1 = element;
elseif(element<=10)
    row1 = 2;
    col1 = element-2;
elseif(element<=18)
    row1 = 3;
    col1 = element-10;
else
    row1 = 4;
    col1 = element-18;
end
end