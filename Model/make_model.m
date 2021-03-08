%% modelname
modelname = '2vs3_fault';
fn = 'model';
fn_save = [fn '.mat'];
fig_save = [fn '.png'];

%% Grid parameter
epr_max = 15;
epr_min = 1;
miur_max = 1;
miur_min = 1;
freq_src = 800*1e6;
dxmax = finddx(epr_max, miur_max, freq_src);
disp(['dxmax=' num2str(dxmax) 'm'])

X = 2;
Y = 1.6;
Z = 1;
T = 30 * 1e-9;%s

dx = 0.01; %m
dy = dx;
dz = dx;

nx = round(X/dx);
ny = round(Y/dy);
nz_air = 10;
nz = round(Z/dz) + nz_air;

npmlx = 8;
npmly = npmlx;
npmlz = 8;

outstep_t_wavefield = 5;
outstep_x_wavefield = 2;
outstep_slice = 5;

% dtmax = finddt(epr_min, miur_min, dx, dy, dz);
% disp(['dt_max=' num2str(dtmax/1e-9) 'ns']);
dt = 0.01 *1e-9;
nt = T/dt;

%% background model
ep_bg = ones(nx,ny,nz) * 9;
ep_bg(:,:,1:nz_air) = 1; % air layer

ep = ep_bg;

%% slice
slicez_z = 0.5;
slicex = [nx/2];
slicey = [ny/2];
slicez = [round(slicez_z/dz) + nz_air];

%%% the parameter names above should be changed togather with those in 'model_em.py' %%%
%% layers
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;
surf_ep = [12];
surf_pos = [0.3]; %m
for i = 1:length(surf_pos)
    layer_z_begin = surf_pos(i);
    if i ~= length(surf_pos)
        layer_z_end = surf_pos(i+1);
    else
        layer_z_end = z(end);
    end
    ep(:,:,z >= layer_z_begin & z<=layer_z_end) = surf_ep(i);
end
%% faults model
dot1 = [1.8 0 1];
dot2 = [0.2 1.6 0];
dot3 = [0.8 0 0];
% 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0

dh = 0.5; %m
idh = round(dh / dz);

ep1 = ep;
for iz = (1+nz_air):nz
    if iz-idh <= nz_air
        ep1(:,:,iz) = ep(:,:,nz_air+1);
    else
        ep1(:,:,iz) = ep(:,:,iz-idh);
    end
end

for ix = 1:nx
    xi = (ix-1)*dx;
    for iy = 1:ny
        yi = (iy-1)*dy;
        for iz = nz:-1:(nz_air+1)
            zi = (iz-nz_air)*dz;
            doti = [xi,yi,zi];
            if dot(cross((dot3-dot1),(dot2-dot1)),(doti-dot1)) > 0
                ep(ix,iy,iz) = ep1(ix,iy,iz);
            end
        end
    end
end

%% src and rec para
dx_src = 0.2;
dx_rec = 0.04;
dy_src = dx_src;
dy_rec = dx_rec;
src_margin_nx = npmlx;
src_margin_ny = npmlx;
rec_margin_nx = npmlx;
rec_margin_ny = npmlx;
src_span = 0;

% place src and rec
dnx_src = dx_src / dx;
dnx_rec = dx_rec / dx;
nx_src = round((nx-2*src_margin_nx)/dnx_src);
if mod(nx_src,2) == 0
    nx_src = nx_src-1;
end
nx_rec = round((nx-2*rec_margin_nx)/dnx_rec);
if mod(nx_rec,2) == 0
    nx_rec = nx_rec-1;
end
srcx = ((-floor(nx_src/2):floor(nx_src/2)) * dx_src) + nx/2*dx;
recx = ((-floor(nx_rec/2):floor(nx_rec/2)) * dx_rec) + nx/2*dx;

dny_src = dy_src / dy;
dny_rec = dy_rec / dy;
ny_src = round((ny-2*src_margin_ny)/dny_src);
if mod(ny_src,2) == 0
    ny_src = ny_src-1;
end
ny_rec = round((ny-2*rec_margin_ny)/dny_rec);
if mod(ny_rec,2) == 0
    ny_rec = ny_rec-1;
end
srcy = ((-floor(ny_src/2):floor(ny_src/2)) * dy_src) + ny/2*dy;
recy = ((-floor(ny_rec/2):floor(ny_rec/2)) * dy_rec) + ny/2*dy;

srcz = 0-dz;
[Xs,Ys,Zs] = meshgrid(srcx,srcy,srcz);
[Xr,Yr,Zr] = meshgrid(recx,recy,srcz);
srcx = reshape(Xs,1,[]);
srcy = reshape(Ys,1,[]);
recx = reshape(Xr,1,[]);
recy = reshape(Yr,1,[]);
recz = srcz;
%%
save(fn_save,'modelname','dx','dy','dz','nx','ny','nz','nz_air','T','dt','freq_src',...
    'dx_src','dx_rec','dy_src','dy_rec','srcx','srcy','srcz','recx','recy','recz',...
    'src_margin_nx','src_margin_ny','rec_margin_nx','rec_margin_ny','src_span',...
    'npmlx','npmly','npmlz','outstep_t_wavefield','outstep_x_wavefield',...
    'outstep_slice','slicex','slicey','slicez','ep_bg','ep');

% %% y slice view
% for i = 1:length(y)
%     imagesc(x,z,squeeze(ep(:,i,:))');colorbar
%     pause(0.1)
% end
% 
% %% z slice view
% for i = 1:length(z)
%     imagesc(x,y,squeeze(ep(:,:,i))');colorbar
%     pause(0.1)
% end

%% 3D-view
figure(10)
clf;
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])
set(gca,'fontsize',30,'fontname','Times')
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;
[X,Y,Z] = meshgrid(y,x,z);
p0 = patch(isosurface(Y,X,Z,ep,8.9));
isonormals(X,Y,Z,ep,p0)
p0.FaceColor = 'black';
p0.EdgeColor = 'none';
% p3 = patch(isosurface(X,Y,Z,ep,6));
% isonormals(X,Y,Z,ep,p3)
% p3.FaceColor = 'magenta';
% p3.EdgeColor = 'none';
p5 = patch(isosurface(Y,X,Z,ep,11.9));
isonormals(X,Y,Z,ep,p5)
p5.FaceColor = 'b';
p5.EdgeColor = 'none';
daspect([1,1,1])
view(3); alpha(.3);axis tight
set(gca,'ZDir','reverse')
set(gca,'YDir','reverse')
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
zlim([z(1),z(end)]);
% title('layers with fault')
set(gca,'fontsize',20);
xlabel('x(m)','Fontsize',24);ylabel('y(m)','Fontsize',24);zlabel('depth(m)','Fontsize',24);
camlight('headlight') 
lighting gouraud
camlight('right')
lighting gouraud

hold on
plot3(Xs,Ys,Zs,'r^')
plot3(Xr,Yr,Zr,'b.')
% xlim([0 10]);ylim([0 10]);zlim([-0.5 5])
hold off
set(gca,'fontsize',20,'fontname','Times')
%%
export_fig(fig_save)

%%
function dtmax = finddt(epr_min, miur_min, dx, dy, dz)
    mu0 = 1.2566370614e-6;
    ep0 = 8.8541878176e-12;
    epmin = epr_min * ep0;
    mumin = miur_min * mu0;
    dtmax = 6.0 / 7.0 * sqrt(epmin * mumin / (dx^-2 + dy^-2 + dz^-2));
end

function dxmax = finddx(epr_max, miur_max, fmax)
    mu0 = 1.2566370614e-6;
    ep0 = 8.8541878176e-12;
    epmax = epr_max * ep0;
    mumax = miur_max * mu0;
    wlmin = 1 / (fmax * sqrt(epmax * mumax));
    dxmax = wlmin;
end