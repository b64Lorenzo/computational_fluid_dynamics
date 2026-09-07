function mesh = create_uniform_mesh(L,N)

x = linspace(0,L,N+1)';

dx = diff(x);

mesh.x = x;

mesh.dx_w = [dx(1); dx];
mesh.dx_e = [dx; dx(end)];

end