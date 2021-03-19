%% modelname
filename = mfilename;
modelname = '3layers_with_0.2m_fault_0o';
fn = 'model';
fn_save = [fn '.mat'];
fig_save = [fn '.png'];

%% Grid parameter
epr_max = 15;
epr_min = 1;
miur_max = 1;
miur_min = 1;
freq_src = 300*1e6;
dxmax = finddx(epr_max, miur_max, freq_src);
disp(['dxmax=' num2str(dxmax) 'm'])

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

npmlx = 8;
npmly = npmlx;
npmlz = 8;

outstep_t_wavefield = 5;
outstep_x_wavefield = 2;
outstep_slice = 5;

% dtmax = finddt(epr_min, miur_min, dx, dy, dz);
% disp(['dt_max=' num2str(dtmax/1e-9) 'ns']);
dt = 0.03 *1e-9;
nt = T/dt;

%% background model
ep_bg = ones(nx,ny,nz) * 9;
ep_bg(:,:,1:nz_air) = 1; % air layer

ep = ep_bg;

%% slice
slicex = [nx/2];
slicey = [ny/2];
slicez = [round(1/dz) + nz_air];

%%% the parameter names above should be changed togather with those in 'model_em.py' %%%

%% layers
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;
surf_ep = [6,9,12];
surf_pos = [0.9,1.5,2.4]; %m
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
dot1 = [2 0 0.8];
dot2 = [2.5 0 2.3];
dot3 = [3 5 2.3];
% 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0

dh = 0.2; %m
idh = round(dh / dz);

ep1 = ep;
for iz = (1+nz_air):nz
    ep1(:,:,iz) = ep(:,:,iz-idh);
end

for ix = 1:nx
    xi = (ix-1)*dx;
    for iy = 1:ny
        yi = (iy-1)*dy;
        for iz = nz:-1:(nz_air+1+idh)
            zi = (iz-nz_air)*dz;
            doti = [xi,yi,zi];
            if dot(cross((dot3-dot1),(dot2-dot1)),(doti-dot1)) > 0
                ep(ix,iy,iz) = ep1(ix,iy,iz);
            end
        end
    end
end

%% src&rec
dy_src = 0.5;
dx_src = 0.1;
src_margin_nx = 0.5/dx;
src_margin_ny = 0.5/dy;

src_span = floor(dx_src/2);

dx_rec = dx_src;
dy_rec = dy_src;
rec_margin_nx = src_margin_nx;
rec_margin_ny = src_margin_ny;

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
save(fn_save,'filename','modelname','dx','dy','dz','nx','ny','nz','nz_air','T','dt','freq_src',...
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
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;
[X,Y,Z] = meshgrid(y,x,z);
p0 = patch(isosurface(X,Y,Z,ep,8.9));
isonormals(X,Y,Z,ep,p0)
p0.FaceColor = 'black';
p0.EdgeColor = 'none';
% p3 = patch(isosurface(X,Y,Z,ep,6));
% isonormals(X,Y,Z,ep,p3)
% p3.FaceColor = 'magenta';
% p3.EdgeColor = 'none';
p5 = patch(isosurface(X,Y,Z,ep,11.9));
isonormals(X,Y,Z,ep,p5)
p5.FaceColor = 'black';
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
xlabel('y(m)','Fontsize',24);ylabel('x(m)','Fontsize',24);zlabel('depth(m)','Fontsize',24);
% ����װ�ƣ��Ҹо��������ƣ��������飩�ȽϺ��ʣ�һ��̫���ˡ�
% �����Է�8���ƣ�headlight��ʾͷ�ƣ�����left��right��
camlight('headlight') 
lighting gouraud
camlight('right')
lighting gouraud

% place src and rec
hold on
plot3(Ys,Xs,Zs,'r.')
% plot3(Yr,Xr,Zr,'b.')
% xlim([0 10]);ylim([0 10]);zlim([-0.5 5])
hold off

% saveas(gcf,fig_save)