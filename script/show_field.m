function show_field(data)
st = reshape(data,64,128);
surf(st(:,1:64)-st(:,65:end));
hold;
surf(st(:,1:64)+st(:,65:end));
end
