%% basic para
mu0 = 1.2566370614e-6;
ep0 = 8.8541878176e-12;
c = 299792458;

%% modelname
dx = 0.025;
freq_src = 400*1e6;
dot_zPosition = 2;

filename = mfilename;
% modelname = sprintf('dm%g_g%gm',dot_zPosition,dx);
modelname = 'test_dot';
fprintf('\nmodelname=%s\n',modelname);
% modelname = 'test_0o';
% modelname = 'test_2d';
fn = modelname;
fn_save = [fn '.mat'];
fig_save = [fn '.png'];
src_save = [fn '_src.png'];

% wavelet
[wl, twl, f_thr, wl_tw] = ricker_src(0.01 *1e-9,freq_src,0.85);
% export_fig(src_save)

% Grid parameter
epr_max = 10;
epr_min = 1;
miur_max = 1;
miur_min = 1;
epr0 = 1;
wvl = 3e8/freq_src/sqrt(epr0);
wl_w = wl_tw * 3e8/sqrt(epr0);
fprintf('wavelet width=%gm\n',wl_w);

disp('----------Grid Parameter----------')
disp('<Space Grid>')
X = 5;
Y = 5;
Z = 5;

dx = dx; %m
dy = dx;
dz = dx;

dxmax = finddx(epr_max, miur_max, freq_src);
fprintf('dxmax=%.2gm,dx=%gm(%.2g%%),dy=%gm(%.2g%%),dz=%gm(%.2g%%)\n',dxmax,dx,dx/dxmax*100,dy,dy/dxmax*100,dz,dz/dxmax*100);
fprintf('wavelength=%gm(%gdx)\n',wvl,wvl/dx);

npmlx = 8;
npmly = npmlx;
npmlz = npmlx;

nx = round(X/dx) + 2 * npmlx;
ny = round(Y/dy) + 2 * npmly;
nz_air = 2 + npmlz;
nz = round(Z/dz) + npmlz + nz_air;

X = nx*dx;
Y = ny*dy;
Z = (nz-nz_air)*dz;


% time prarameter
disp('<Time Grid>')
T_ref = sqrt(X^2+Y^2+(2*Z)^2)/(3e8/sqrt(epr0));%s
fprintf('Max 2-way traveltime=%gns\n',T_ref/1e-9)
T = 1.1*T_ref;%s

dtmax = finddt(epr_min, miur_min, dx, dy, dz);
dt = 0.7 * dtmax;
nt = T/dt;
fprintf('t_total=%gns.dt=%gns,nt=%g\n',T/1e-9,dt/1e-9,nt);
fprintf('T=%gns(%gdt)\n',1e9/freq_src,1/(dt*freq_src));

check_dispersion(dx,dt,f_thr,c/sqrt(epr_max*miur_max));
check_stability(dx,dt);

outstep_x_wavefield = ceil((nx*ny*nz/150^3)^(1/3));
fprintf('outstep_x_wavefield=%g\n',outstep_x_wavefield)
outstep_t_wavefield = 5;
outstep_slice = outstep_t_wavefield;

disp('----------Grid Parameter End----------')

%% background model
ep_bg = ones(nx,ny,nz) * epr0;
ep_bg(:,:,1:nz_air) = 1; % air layer

ep = ep_bg;

%% slice
slicez_z = dot_zPosition;
slicex = [nx/2];
slicey = [ny/2];
slicez = [round(slicez_z/dz) + nz_air];

%%% the parameter names above should be changed togather with those in 'model_em.py' %%%
%% dot model
r = 0.5;
dots{1} = [X/2,Y/2,slicez_z];
% dots{2} = [X/2+r,Y/2,slicez_z];

for ix = 1:nx
    xi = ix * dx;
        for iy = 1:ny
            yi = iy * dy;
            for iz = 1:nz
                zi = (iz-nz_air) * dz;
                    dot = dots{1};
                    if (xi - dot(1))^2 + (yi - dot(2))^2 + (zi - dot(3))^2 <= r^2
                        ep(ix,iy,iz) = 6;
                    end
%                     dot = dots{2};
%                     if (xi - dot(1))^2 + (yi - dot(2))^2 + (zi - dot(3))^2 <= r^2
%                         ep(ix,iy,iz) =ep0-3;
%                     end
            end
        end
end

%% src and rec para
dx_src = 0.5;
dx_rec = 0.1;
% dx_rec = 0.2;
dy_src = dx_src;
dy_rec = dx_rec;


src_margin_nx = npmlx;
src_margin_ny = npmly;
rec_margin_nx = npmlx;
rec_margin_ny = npmly;
src_span = 0;

% place src and rec
dnx_src = dx_src / dx;
dnx_rec = dx_rec / dx;
nx_src = floor((nx-2*src_margin_nx)/dnx_src);
if mod(nx_src,2) == 0
    nx_src = nx_src-1;
end
nx_rec = floor((nx-2*rec_margin_nx)/dnx_rec);
if mod(nx_rec,2) == 0
    nx_rec = nx_rec-1;
end
srcx = ((-floor(nx_src/2):floor(nx_src/2)) * dx_src) + nx/2*dx;
recx = ((-floor(nx_rec/2):floor(nx_rec/2)) * dx_rec) + nx/2*dx;

dny_src = dy_src / dy;
dny_rec = dy_rec / dy;
ny_src = floor((ny-2*src_margin_ny)/dny_src);
if mod(ny_src,2) == 0
    ny_src = ny_src-1;
end
ny_rec = floor((ny-2*rec_margin_ny)/dny_rec);
if mod(ny_rec,2) == 0
    ny_rec = ny_rec-1;
end
srcy = ((-floor(ny_src/2):floor(ny_src/2)) * dy_src) + ny/2*dy;
recy = ((-floor(ny_rec/2):floor(ny_rec/2)) * dy_rec) + ny/2*dy;

srcz = -0.1;
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
figure(10)
clf;
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])
set(gca,'fontsize',30,'fontname','Times')
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;
[X,Y,Z] = meshgrid(y,x,z);
p0 = patch(isosurface(Y,X,Z,ep,1));
isonormals(X,Y,Z,ep,p0)
p0.FaceColor = 'black';
p0.EdgeColor = 'none';
% p3 = patch(isosurface(X,Y,Z,ep,ep0-3));
% isonormals(X,Y,Z,ep,p3)
% p3.FaceColor = 'magenta';
% p3.EdgeColor = 'none';
% p5 = patch(isosurface(Y,X,Z,ep,ep0+3));
% isonormals(X,Y,Z,ep,p5)
% p5.FaceColor = 'blue';
% p5.EdgeColor = 'none';
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
function [wavelet,t_wavelet,f_thr,width] = ricker_src(dt,f,thr,tlength)
    disp('<Ricker Wavelet>')
    if(nargin<4)
       tlength=4*1/f;
    end
    fs = 1/dt;
    [~,t_wavelet] = ricker(dt,f,tlength);
    wavelet = (1 - 2 * (pi*f.*t_wavelet).^2) .* exp(-(pi*f.*t_wavelet).^2);
    subplot(2,1,1)
    plot(t_wavelet,wavelet)
    width = 2*0.88521/(2*pi*f);
    fprintf('width:%gm\n',width);
    hold on
    plot([-width/2,width/2],[0.5,0.5],'r-')
    plot([t_wavelet(1),-width/2],[0.5,0.5],'r--')
    hold off
    % f ana
    wly = fft(wavelet);
    n = length(wavelet);
    wly0 = fftshift(wly);
    wlf0 = (-n/2:n/2-1)*(fs/n);
    power = abs(wly0).^2/n;
    subplot(2,1,2)
    wlf0p = wlf0(wlf0>0&wlf0<(4*f));
    powerp = power(wlf0>0&wlf0<(4*f));
    plot(wlf0p/1e6,powerp)
    xlabel('Frequency(MHz)')
    ylabel('Power')
    f_thr = max(wlf0(power>((1-thr)*max(power))));
    fprintf('freq of %d%% max freq of riker wavelet:%gMHz\n',thr*100,f_thr/1e6);
    hold on
    plot([f_thr/1e6,f_thr/1e6],[0,max(power)],'--')
    hold off
end

function [cx,ct] = check_dispersion(dx,dt,fmax,cmin)
    disp('-----Dispersion Check-----')
    w = 2*pi*fmax;
    k0 = w/cmin;
    ax = k0*dx/2;
    at = w*dt/2;
    bx = sin(ax)/ax;
    bt = sin(at)/at;
    cx = 1/bx;
    ct = bt;
    fprintf("kdx/2=%g, wdt/2=%g (should be both << 1)\n",ax,at);
    fprintf("k/k0≈%g(for dx), k/k0=%g(for dt)\n",cx,ct);
    disp('-----Dispersion Check End-----')
end

function alpha = check_stability(dx, dt, cmax, dim)
    disp('-----Stability Check-----')
    if nargin < 3
        cmax = 299792458;%m/s
    end
    if length(dx) > 1
        dx = 1/sqrt(sum(dx.^-2));
    else
        if nargin < 4
            dim = 3;
        end
        dx = dx/sqrt(dim);
    end
    
    alpha = cmax*dt/dx;
    fprintf("cdt/dx=%g, should ≤ 2/pi(PSTD) or 1(FDTD_Yee) or 6/7(FDTD_Zhu)\n",alpha);
    if alpha > 1
        disp('!!!Failed.!!!')
    elseif alpha > 6/7
        disp('!!Failed for FDTD_Zhu.!!')
    elseif alpha > 2/pi
        disp('!!Failed for PSTD.!!')
    else
        disp('Success.')
    end
    disp('-----Stability Check End-----')
end

function dxmax = finddx(epr_max, miur_max, fmax)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    mu0 = 1.2566370614e-6;
    ep0 = 8.8541878176e-12;
    epmax = epr_max * ep0;
    mumax = miur_max * mu0;
    wlmin = 1 / (fmax * sqrt(epmax * mumax));
    dxmax = wlmin;
end

function dtmax = finddt(epr_min, miur_min, dx, dy, dz)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    mu0 = 1.2566370614e-6;
    ep0 = 8.8541878176e-12;
    epmin = epr_min * ep0;
    mumin = miur_min * mu0;
    dtmax = 6.0 / 7.0 * sqrt(epmin * mumin / (dx^-2 + dy^-2 + dz^-2));
end
