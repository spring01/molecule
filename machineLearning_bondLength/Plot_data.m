function Plot(obj)

allShapes = 4000 .* [ ...
    1             1 ...
    3 3 2 2 2 2 2 2 ...
    4 4 3 3 3 3 3 3 ...
    5 5];
allColors = 0.12 .* [ ...
    [8 8 8];                                                       [8 7 7]; ...
    [0 8 8]; [8 0 8]; [8 8 0]; [5 5 5]; [0 0 8]; [8 0 0]; [0 8 0]; [6 5 5];
    [0 6 6]; [6 0 6]; [6 6 0]; [3 3 3]; [0 0 6]; [6 0 0]; [0 6 0]; [4 3 3];
    [0 4 4]; [4 0 4];];

geom = obj.cartesian(:, 2:end);
shapes = allShapes(obj.cartesian(:,1));
% colors = zeros(size(obj.cartesian,1), 1);
colors = allColors(obj.cartesian(:,1), :);
bondLines = {};

BondLengthTensor = InitializeBondLengthTensor();
for iAtom = 1:size(obj.cartesian,1)
        
    % determine bonds 
    for jAtom = iAtom+1:size(obj.cartesian,1)
        numBonds = NumberOfBonds(obj.cartesian, iAtom, jAtom, BondLengthTensor);
        vector = geom(iAtom, :) - geom(jAtom, :);
        normal_vector = ones(1, 3);
        normal_vector(3) = -(vector(1) + vector(2))/vector(3);
        normal_vector = normal_vector ./ norm(normal_vector);
        switch(numBonds)
            case(1)
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)];
            case(2)
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)] ...
                    - 0.03.*[normal_vector; normal_vector];
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)] ...
                    + 0.03.*[normal_vector; normal_vector];
            case(3)
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)];
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)] ...
                    - 0.05.*[normal_vector; normal_vector];
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)] ...
                    + 0.05.*[normal_vector; normal_vector];
        end
    end
end
figure();
axis vis3d;
hold;
for oneLine = bondLines
    line(oneLine{1}(:, 1), oneLine{1}(:, 2), oneLine{1}(:, 3), 'LineWidth', 4);
end
scatter3(geom(:,1), geom(:,2), geom(:,3), shapes, colors, 'fill');

end


function numBonds = NumberOfBonds(cartesian, atom1, atom2, BondLengthTensor)

distance = norm(cartesian(atom1, 2:end) - cartesian(atom2, 2:end));

indices = sort([cartesian(atom1, 1), cartesian(atom2, 1)]);

criterions = reshape(BondLengthTensor(indices(1), indices(2), :), [], 1);
criterions(2:3) = (criterions(1:2) + 4.*criterions(2:3)) ./ 5;
criterions(1) = 1.2 .* criterions(1);

numBonds = 3 - sum(distance>criterions);

end

function bondLengthTensor = InitializeBondLengthTensor()
singleBondMat = ...
    [...H    He    Li    Be     B     C     N     O     F    Ne    Na    Mg    Al    Si     P     S    Cl    Ar     K    Ca
    0.741 0.000 1.595 1.343 1.320 1.140 1.090 1.000 1.014 0.000 1.887 1.730 1.648 1.520 1.435 1.374 1.315 1.292 2.243 2.003 % H
    0.000 1.081 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 % He
    0.000 0.000 2.673 0.000 0.000 0.000 0.000 1.688 1.564 0.000 2.889 0.000 0.000 0.000 0.000 0.000 2.180 0.000 3.270 0.000 % Li
    0.000 0.000 0.000 2.460 0.000 0.000 0.000 1.331 1.400 0.000 0.000 0.000 0.000 0.000 0.000 1.742 1.797 0.000 0.000 0.000 % Be
    0.000 0.000 0.000 0.000 1.763 1.491 1.645 1.265 1.317 0.000 0.000 0.000 0.000 0.000 0.000 1.609 1.750 0.000 0.000 0.000 % B
    0.000 0.000 0.000 0.000 0.000 1.580 1.492 1.448 1.398 0.000 0.000 0.000 0.000 1.875 1.858 1.849 1.813 0.000 0.000 0.000 % C
    0.000 0.000 0.000 0.000 0.000 0.000 1.864 1.507 1.512 0.000 0.000 0.000 1.786 1.575 1.491 1.497 1.975 0.000 0.000 0.000 % N
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 1.516 1.421 0.000 1.950 1.767 1.618 1.510 1.512 1.574 1.696 0.000 0.000 1.976 % O
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 1.412 0.000 1.926 1.770 1.654 1.604 1.593 1.646 1.697 0.000 2.171 1.967 % F
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 3.100 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 % Ne
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 3.079 0.000 0.000 0.000 0.000 0.000 2.584 0.000 3.589 0.000 % Na
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 3.891 0.000 0.000 0.000 2.143 2.199 0.000 0.000 0.000 % Mg
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 2.701 0.000 0.000 2.029 2.240 0.000 0.000 0.000 % Al
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 2.320 0.000 1.929 2.076 0.000 0.000 0.000 % Si
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 2.219 1.900 2.214 0.000 0.000 0.000 % P
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 2.059 2.076 0.000 0.000 0.000 % S
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 1.988 0.000 2.667 2.437 % Cl
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 % Ar
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 3.905 0.000 % K
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 % Ca
    ];

doubleBondMat = zeros(20);
doubleBondMat(5, 7) = 1.325; % B=N
doubleBondMat(5, 8) = 1.265; % B=O
doubleBondMat(6, 6) = 1.382; % C=C
doubleBondMat(6, 7) = 1.338; % C=N
doubleBondMat(6, 8) = 1.245; % C=O
doubleBondMat(7, 7) = 1.252; % N=N
doubleBondMat(7, 8) = 1.258; % N=O

tripleBondMat = zeros(20);
tripleBondMat(6, 6) = 1.246; % C#C
tripleBondMat(6, 7) = 1.177; % C#N
tripleBondMat(7, 7) = 1.133; % N#N

bondLengthTensor = zeros(20,20,3);
bondLengthTensor(:,:,1) = singleBondMat;
bondLengthTensor(:,:,2) = doubleBondMat;
bondLengthTensor(:,:,3) = tripleBondMat;
end
