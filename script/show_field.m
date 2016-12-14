function show_field(data)
x = 32;
y = 128;
st = reshape(data,y,x*2);
surf(st(:,1:x)-st(:,x+1:end));
hold;
surf(st(:,1:x)+st(:,x+1:end));
end
