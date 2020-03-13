%% modelname
modelname = 'obstacles';

%% Grid parameter
dx = 0.05; %m
dy = dx;
dz = dx;

nx = 200;
ny = 200;
nz_air = 10;
nz = 120 + nz_air;

%% background model
ep_bg = ones(nx,ny,nz) * 9;
ep_bg(:,:,1:nz_air) = 1; % air layer
ep = ep_bg;

%% slice
slicex = [50,150];
slicey = [50,150];
slicez = [50 + nz_air];

%%% the parameter names above should be changed togather with those in 'model_em.py' %%%

%% obstacles
x = (0:(nx-1)) * dx;
y = (0:(ny-1)) * dy;
z = ((0:(nz-1)) - nz_air) * dz;

% cuboid
xm = 2.5;
ym = 2.5;
zm = 2.5;
lx = 2;
ly = 2;
lz = 2;
ep(x>(xm-lx/2)&x<(xm+lx/2),y>(ym-ly/2)&y<(ym+ly/2),z>(zm-lz/2)&z<(zm+lz/2)) = 12;

% sphere
r = 1;
xm = 7.5;
ym = 2.5;
zm = 2.5;
for ix = floor((xm-r)/dx):ceil((xm+r)/dx)
    for iy = floor((ym-r)/dy):ceil((ym+r)/dy)
        for iz = floor((zm-r)/dz):ceil((zm+r)/dz)
            if norm([ix*dx,iy*dx,iz*dx] - [xm,ym,zm]) < r
                ep(ix,iy,iz + nz_air) = 6;
            end
        end
    end
end

% Spherical cavity
r = 1.2;
r2 = 0.8;
xm = 2.5;
ym = 7.5;
zm = 2.5;
for ix = floor((xm-r)/dx):ceil((xm+r)/dx)
    for iy = floor((ym-r)/dy):ceil((ym+r)/dy)
        for iz = floor((zm-r)/dz):ceil((zm+r)/dz)
            if norm([ix*dx,iy*dx,iz*dx] - [xm,ym,zm]) < r2
                ep(ix,iy,iz + nz_air) = 1;
            elseif norm([ix*dx,iy*dx,iz*dx] - [xm,ym,zm]) < r
                ep(ix,iy,iz + nz_air) = 12;
            end
        end
    end
end

% % regular tetrahedron
% z_top = 1.5;
% h = sqrt(6)*2/3;
% xm = 11;
% ym = 3;
% k_yz = sqrt(3)/(sqrt(6)*2/3);
% k_xy = 2/sqrt(3);
% for iz = 0:round(h/dz)
%     y_max = round((iz*dz*k_yz)/dy);
%     zi = iz + nz_air + round(z_top/dz);
%     for iy = 0:round(y_max)
%         yi = iy + round(ym/dy) + round(-2/3*y_max);
%         x_max = round((k_xy*iy*dy)/dx);
%         ep((round(-0.5*x_max):round(x_max/2)) + round(xm/dx),yi,zi) = 15;
%     end
% end
% reverse regular tetrahedron
z_top = 1.5;
h = sqrt(6)*2/3;
xm = 7.5;
ym = 7.5;
k_yz = sqrt(3)/(sqrt(6)*2/3);
k_xy = 2/sqrt(3);
for iz = 0:round(h/dz)
    y_max = round(((round(h/dz)-iz)*dz*k_yz)/dy);
    zi = iz + nz_air + round(z_top/dz);
    for iy = 0:round(y_max)
        yi = iy + round(ym/dy) + round(-2/3*y_max);
        x_max = round((k_xy*iy*dy)/dx);
        ep((round(-0.5*x_max):round(x_max/2)) + round(xm/dx),yi,zi) = 12;
    end
end

%%
save('model.mat','modelname','dx','dy','dz','nx','ny','nz','nz_air','slicex','slicey','slicez','ep_bg','ep');

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
[X,Y,Z] = meshgrid(y,x,z);
p0 = patch(isosurface(X,Y,Z,ep,1));
isonormals(X,Y,Z,ep,p0)
p0.FaceColor = 'red';
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
p4 = patch(isosurface(X,Y,Z,ep,9));
isonormals(X,Y,Z,ep,p4)
p4.FaceColor = 'black';
p4.EdgeColor = 'none';
p5 = patch(isosurface(X,Y,Z,ep,11.9));
isonormals(X,Y,Z,ep,p5)
p5.FaceColor = 'cyan';
p5.EdgeColor = 'none';
daspect([1,1,1])
view(3); alpha(.3);axis tight
set(gca,'ZDir','reverse')
set(gca,'YDir','reverse')
zlim([z(1),z(end)]);
title('A Complex Model')
xlabel('y(m)');ylabel('x(m)');zlabel('depth(m)');
% 下面装灯，我感觉放两个灯（拷贝两遍）比较合适，一个太暗了。
% 最多可以放8个灯，headlight表示头灯，还有left和right。
% camlight('headlight') 
% lighting gouraud
% camlight('right')
% lighting gouraud