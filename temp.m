function cartesian = ZMatrix2Cartesian(zMatrix)
x = zeros(size(zMatrix,1) + 2, 1);
y = zeros(size(zMatrix,1) + 2, 1);
z = zeros(size(zMatrix,1) + 2, 1);
x(1) = 1; y(1) = 1; z(1) = 0;
x(2) = 0; y(2) = 1; z(2) = 0;
x(3) = 0; y(3) = 0; z(3) = 0;

xyz = zeros(size(zMatrix,1) + 2, 3);
for i = 2:size(zMatrix,1)
    k1 = zMatrix(i, 2);
    k2 = zMatrix(i, 4);
    k3 = zMatrix(i, 6);
    b = zMatrix(i, 3);
    if(i == 2)
        k2 = i - 2;
        a = 90;
    else
        a = zMatrix(i, 5);
    end
    if(i <= 3)
        k3 = i - 3;
        t = 0;
    else
        t = zMatrix(i, 7);
    end
    x(i + 2) = x(k1 + 2); y(i + 2) = y(k1 + 2); z(i + 2) = z(k1 + 2);
    
    x1 = x(k2 + 2) - x(k1 + 2); y1 = y(k2 + 2) - y(k1 + 2); z1 =z(k2 + 2) - z(k1 + 2);
    norm = sqrt(x1^2 + y1^2 + z1^2);
    x1 = x1 / norm; y1 = y1/norm; z1 = z1 / norm;
    
    x2 = x(k3 + 2) - x(k2 + 2); y2 = y(k3 + 2) - y(k2 + 2); z2 = z(k3 + 2) - z(k2 + 2);
    norm = x1 * x2 + y1 * y2 + z1 * z2;
    x2 = x2 - norm * x1; y2 = y2 - norm * y1; z2 = z2 - norm * z1;
    norm = sqrt(x2^2 + y2^2 + z2^2);
    if(norm > 0)
        x2 = x2 / norm; y2 = y2 / norm; z2 = z2 / norm;
    end
    
    x3 = y1 * z2 - y2 * z1;
    y3 = z1 * x2 - z2 * x1;
    z3 = x1 * y2 - x2 * y1;
    
    x(i + 2) = x(i + 2) + b * (cosd(a) * x1 + sind(a) * ( cosd(t) * x2 - sind(t) * x3));
    y(i + 2) = y(i + 2) + b * (cosd(a) * y1 + sind(a) * ( cosd(t) * y2 - sind(t) * y3));
    z(i + 2) = z(i + 2) + b * (cosd(a) * z1 + sind(a) * ( cosd(t) * z2 - sind(t) * z3));
    
    
end
% Gaussian convention is apparently different, this fixes it
cartesian = [y(3:end) z(3:end) x(3:end)];
% coord = [x(3:end) y(3:end) z(3:end)];
end