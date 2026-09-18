
function mesh = nonuniform_121(L,k)

x = [ ...
    0
    0.2*L - 1*k
    0.4*L - 2*k
    0.6*L - 3*k
    0.8*L - 4*k
    1.0*L - 5*k
    1.0*L - 4*k
    1.0*L - 3*k
    1.0*L - 2*k
    1.0*L - 1*k
    1.0*L ];

dx = diff(x);

mesh.x    = x;
mesh.dx_w = [dx(1); dx];
mesh.dx_e = [dx; dx(end)];

end