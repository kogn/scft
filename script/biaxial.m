function res = biaxial(output)
[m n] = size(output);
res = zeros(m,1);
tmp = add_eig(output);
for i = 1:m
    res(i) = 1-6*(sum(tmp(i,:).^3))^2/(sum(tmp(i,:).^2))^3;
    %res(i) = (sum(tmp(i,:).^2))^3-6*(sum(tmp(i,:).^3))^2;
end

