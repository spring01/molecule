rng(100);


bigMat1 = rand(100);
bigMat2 = rand(100);
bigMat1 = bigMat1*bigMat1';
bigMat2 = bigMat2*bigMat2';

tic
for i = 1:100000
    a = trace(bigMat1*bigMat2);
end
toc

tic
for i = 1:100000
    b = sum(sum(bigMat1.*bigMat2));
end
toc

tic
for i = 1:100000
    c = reshape(bigMat1, 1, []) * reshape(bigMat2, [], 1);
end
toc

tic
for i = 1:100000
    c = reshape(bigMat1, [], 1)' * reshape(bigMat2, [], 1);
end
toc

tic
for i = 1:100000
    c = ones(1,numel(bigMat1)) * reshape(bigMat1.*bigMat2, [], 1);
end
toc


a
b
c
