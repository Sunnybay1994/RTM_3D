%%
nx = 250;
ny = 250;

loc = load('merge_gather_loc_0060.dat');
gather = load('merge_gather_0060.dat');
loc_nmo = load('merge_gather_loc_NMO_0060.dat');
gather_nmo = load('merge_gather_NMO_0060.dat');

%%
img = zeros(nx,ny,size(gather,2));
img_nmo = zeros(nx,ny,size(gather_nmo,2));

for i=1:size(loc,1)
    img(loc(i,1),loc(i,2),:) = gather(i,:);
end

for i=1:size(loc_nmo,1)
    img_nmo(loc_nmo(i,1),loc_nmo(i,2),:) = gather_nmo(i,:);
end

%%
X = 5;
Y = 5;
Z = 3;
T = 60 * 1e-9;%s

dx = 0.02; %m
dy = dx;
dz = dx;
dt = 0.03 *1e-9;

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

%%
nt = size(img,3);
x=(1:nx)*dx;
tt = (1:nt)*dt;
for i=125
    subplot(1,2,1)
    imagesc(x,tt,squeeze(img(i,:,:))');colorbar
    hold on
    plot(srcx(6),dt,'r^','MarkerSize',12)
    plot(recx,dt*ones(size(recx)),'bo')
    hold off
    title(['image at y\_grid=' num2str(i)])
    caxis([-6e-3,6e-3])
    subplot(1,2,2)
    imagesc(x,tt,squeeze(img_nmo(i,:,:))');colorbar
    hold on
    plot(srcx(6),dt,'r^','MarkerSize',12)
    plot(recx,dt*ones(size(recx)),'bo')
    hold off
    title(['image after nmo at y\_grid=' num2str(i)])
    caxis([-6e-3,6e-3])
%     pause()
end