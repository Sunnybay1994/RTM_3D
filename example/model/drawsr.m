X = 5;
Y = 5;
Z = 3;
T = 60 * 1e-9;%s

dx = 0.02; %m
dy = dx;
dz = dx;

nx = round(X/dx);
ny = round(Y/dy);
nz_air = 10;
nz = round(Z/dz) + nz_air;

dx_src = 0.4;
dx_rec = 0.2;
dx_nmo = dx_rec/2;
dnx_src = dx_src / dx;
dnx_rec = dx_rec / dx;
nx_src = round((nx-2*round(0.4/dx))/dnx_src);
if iseven(nx_src)
    nx_src = nx_src-1;
end
nx_rec = round((nx-2*round(0.2/dx))/dnx_rec);
if iseven(nx_rec)
    nx_rec = nx_rec-1;
end
srcx = ((-floor(nx_src/2):floor(nx_src/2)) * dx_src) + nx/2*dx;
recx = ((-floor(nx_rec/2):floor(nx_rec/2)) * dx_rec) + nx/2*dx;
nmox = recx(1)+dx_nmo:dx_nmo:recx(end)-dx_nmo;
srcz = 0-dz;
[Xs,Ys,Zs] = meshgrid(srcx,srcx,srcz);
[Xr,Yr,Zr] = meshgrid(recx,recx,srcz);
[Xn,Yn,Zn] = meshgrid(nmox,nmox,srcz);
figure(1)
daspect([1,1,1])
hold on
plot3(Xs,Ys,Zs,'r^','MarkerSize',12)
plot3(Xr,Yr,Zr,'bo')
plot3(Xn,Yn,Zn,'k.')
hold off