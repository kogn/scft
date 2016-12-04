function A = sym(l)
A = zeros(2*l+1,2*l+1);
for i=1:(2*l)
  A(i,i+1) = 0.5*c(i-l-1,l);
  A(i+1,i) = -0.5*c(i-l-1,l);
end
end

function res = c(j,l);
res = sqrt((l-j)*(l+j+1));
end
