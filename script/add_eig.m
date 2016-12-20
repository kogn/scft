function res = test(output)
%load ./data/output
[m n] = size(output);
res = zeros(m,3);
A = zeros(3,3);
for i=1:m
    A(1,1) = output(i,2);
    A(2,2) = output(i,3);
    A(3,3) = output(i,4);
    A(1,2) = output(i,5);
    A(1,3) = output(i,6);
    A(2,3) = output(i,7);
    A(2,1) = output(i,5);
    A(3,1) = output(i,6);
    A(3,2) = output(i,7);
    res(i,:) = eig(A);
end
