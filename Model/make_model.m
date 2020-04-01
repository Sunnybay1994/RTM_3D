%% modelname
modelname = 'layers_with_fault';
fn = 'model';
fn_save = [fn '.mat'];
fig_save = [fn '.png'];

%% Grid parameter
epr_max = 15;
epr_min = 1;
miur_max = 1;
miur_min = 1;
fmax = 300*1e6;
dxmax = finddx(epr_max, miur_max, fmax);
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

dtmax = finddt(epr_min, miur_min, dx, dy, dz);
disp(['dt_max=' num2str(dtmax/1e-9) 'ns']);
dt = 0.03 *1e-9;
nt = T/dt;

%% background model
ep_bg = ones(nx,ny,nz) * 9;
% ep_bg(:,:,1:nz_air) = 1; % air layer

ep = ep_bg;

%% slice
slicex = [nx/2];
slicey = [ny/2];
slicez = [round(1.2/dz) + nz_air];

%%% the parameter names above should be changed togather with those in 'model_em.py' %%%

%% layers
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;
surf_ep = [5,9,13];
surf_pos = [0.8,1.3,2.3]; %m
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
% -6*(x-2)+3*(y-2)+2*z=0
dh = 0.3; %m
idh = round(dh / dz);
for ix = 1:nx
    xi = (ix-1)*dx;
    for iy = 1:ny
        yi = (iy-1)*dy;
        for iz = nz:-1:(nz_air+1+idh)
            zi = (iz-nz_air)*dz;
            doti = [xi,yi,zi];
            if dot(cross((dot3-dot1),(dot2-dot1)),(doti-dot1)) > 0
                ep(ix,iy,iz) = ep(ix,iy,iz-idh);
            end
        end
    end
end


%%
save(fn_save,'modelname','dx','dy','dz','nx','ny','nz','nz_air','T','dt',...
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
p0 = patch(isosurface(X,Y,Z,ep,1));
isonormals(X,Y,Z,ep,p0)
p0.FaceColor = 'black';
p0.EdgeColor = 'none';
% p1 = patch(isosurface(X,Y,Z,ep,14.9));
% isonormals(X,Y,Z,ep,p1)
% p1.FaceColor = 'yellow';
% p1.EdgeColor = 'none';
% p2 = patch(isosurface(X,Y,Z,ep,3));
% isonormals(X,Y,Z,ep,p2)
% p2.FaceColor = 'blue';
% p2.EdgeColor = 'none';
p3 = patch(isosurface(X,Y,Z,ep,6));
isonormals(X,Y,Z,ep,p3)
p3.FaceColor = 'magenta';
p3.EdgeColor = 'none';
% p4 = patch(isosurface(X,Y,Z,ep,8.9));
% isonormals(X,Y,Z,ep,p4)
% p4.FaceColor = 'black';
% p4.EdgeColor = 'none';
p5 = patch(isosurface(X,Y,Z,ep,11.9));
isonormals(X,Y,Z,ep,p5)
p5.FaceColor = 'blue';
p5.EdgeColor = 'none';
daspect([1,1,1])
view(3); alpha(.3);axis tight
set(gca,'ZDir','reverse')
set(gca,'YDir','reverse')
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
zlim([z(1),z(end)]);
title('layers with fault')
xlabel('y(m)');ylabel('x(m)');zlabel('depth(m)');
% 下面装灯，我感觉放两个灯（拷贝两遍）比较合适，一个太暗了。
% 最多可以放8个灯，headlight表示头灯，还有left和right。
camlight('headlight') 
lighting gouraud
camlight('right')
lighting gouraud

% place src and rec
hold on
dx_src = 0.4;
dx_rec = 0.2;
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
srcz = 0-dx;
[Xs,Ys,Zs] = meshgrid(srcx,srcx,srcz);
[Xr,Yr,Zr] = meshgrid(recx,recx,srcz);
plot3(Xs,Ys,Zs,'r^')
plot3(Xr,Yr,Zr,'b.')
% xlim([0 10]);ylim([0 10]);zlim([-0.5 5])
hold off

saveas(gcf,fig_save)