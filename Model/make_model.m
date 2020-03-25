%% modelname
modelname = 'test3';
fn = 'test3_bigger_T';
fn_save = [fn '.mat'];
fig_save = [fn '.png'];

%% Grid parameter
epr_max = 15;
epr_min = 1;
miur_max = 1;
miur_min = 1;
fmax = 400*1e6;
dxmax = finddx(epr_max, miur_max, fmax);
disp(['dxmax=' num2str(dxmax) 'm'])

X = 4;
Y = 4;
Z = 3;
T = 120 *1e-9;%s

dx = 0.05; %m
dy = dx;
dz = dx;

nx = round(X/dx);
ny = round(Y/dy);
nz_air = 6;
nz = round(Z/dz) + nz_air;

npmlx = 8;
npmly = npmlx;
npmlz = 5;

outstep_t_wavefield = 5;
outstep_x_wavefield = 1;
outstep_slice = 5;

dtmax = finddt(epr_min, miur_min, dx, dy, dz);
disp(['dt_max=' num2str(dtmax/1e-9) 'ns']);
dt = 0.06 *1e-9;

%% background model
ep_bg = ones(nx,ny,nz) * 9;
% ep_bg(:,:,1:nz_air) = 1; % air layer

ep = ep_bg;

%% slice
slicex = [nx/2];
slicey = [ny/2];
slicez = [round(1.2/dz) + nz_air];

%%% the parameter names above should be changed togather with those in 'model_em.py' %%%

% sphere
r = 0.5;
xm = X/2;
ym = Y/2;
zm = 1.5;
for ix = floor((xm-r)/dx):ceil((xm+r)/dx)
    for iy = floor((ym-r)/dy):ceil((ym+r)/dy)
        for iz = floor((zm-r)/dz):ceil((zm+r)/dz)
            if norm([ix*dx,iy*dx,iz*dx] - [xm,ym,zm]) < r
                ep(ix,iy,iz + nz_air) = 6;
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
dx_src = 0.6;
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