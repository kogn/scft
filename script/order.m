function test(output)
%load ./data/output
A = zeros(3,3);
A(1,1) = output(1);
A(2,2) = output(2);
A(3,3) = output(3);
A(1,2) = output(4);
A(1,3) = output(5);
A(2,3) = output(6);
A(2,1) = output(4);
A(3,1) = output(5);
A(3,2) = output(6);
eig(A)
B = zeros(3,3);
B(1,1) = output(7);
B(2,2) = output(8);
B(3,3) = output(9);
B(1,2) = output(10);
B(1,3) = output(11);
B(2,3) = output(12);
B(2,1) = output(10);
B(3,1) = output(11);
B(3,2) = output(12);
eig(B)
