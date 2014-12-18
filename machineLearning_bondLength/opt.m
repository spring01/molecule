
% options = optimoptions('lsqnonlin', 'MaxFunEvals', 1e4, 'TolX', 1e-30, 'TolFun', 1e-30);

for i = 1:100
    [coeffs{i}, err(i)] = lsqnonlin(@error1, 10*rand(8,1));
    disp(i);
%     coeffs{i} = 5*rand(5,1);
%     err(i) = sum(error1(coeffs{i}).^2);
end

