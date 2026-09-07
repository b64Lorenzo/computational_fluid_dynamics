function mesh = create_nonuniform_mesh(L,N,beta)

eta = linspace(0,1,N+1)';

x = L*(exp(beta*eta)-1)/(exp(beta)-1);

dx = diff(x);

mesh.x = x;

mesh.dx_w = [dx(1); dx];
mesh.dx_e = [dx; dx(end)];

end