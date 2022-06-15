clear all
%% basic para
mu0 = 1.2566370614e-6;
ep0 = 8.8541878176e-12;
c = 299792458;

% modelname
zmode = true;
freq = 1
00;
freq_src = freq * 1e6;
filename = mfilename;

suffix = '';
suffix = [suffix, sprintf('f%g',freq/100)];
    
modelname = [ sprintf('glyc'),suffix ];
fprintf('\nmodelname=%s\n',modelname);
fn = modelname;
fn_save = [fn '.mat'];
fig_save = [fn '.png'];
src_save = [fn '_src.png'];


% Grid parameter
epr_bg = 3 * 4;% half v in media

epr_max = epr_bg;
epr_min = 1;
miur_max = 1;
miur_min = 1;


% source wavelet
[wl, twl, f_thr, wl_tw] = ricker_src(0.01 *1e-9,freq_src,0.85);
wvl = 3e8/freq_src/sqrt(epr_bg);
wl_w = wl_tw * 3e8/sqrt(epr_bg);
fprintf('wavelet width in medium=%gm\n',wl_w);

disp('↓↓↓↓↓↓↓↓↓↓Grid Parameter↓↓↓↓↓↓↓↓↓↓')
disp('<Space Grid>')
dx = 0.1; %m
dy = dx;

X = 60;
% nx = 1210;
if freq == 100
    Y = 5; % 0.5:0.5:4.5
%     ny = 115;
    Z = 12;
    dz = 0.04;
    dt = 0.1466*10^-9/4;
    nt = 3693;
elseif freq == 400
    Y = 6; % 0.5:0.5:5.5
%     ny = 135;
    Z = 6;
    dz = 0.02;
    nt = 1903;
    dt = 0.0587*10^-9/2;
end

dxmax = finddx(epr_max, miur_max, freq_src);
fprintf('dxmax=%.2gm,dx=%gm(%.2g%%),dy=%gm(%.2g%%),dz=%gm(%.2g%%)\n',dxmax,dx,dx/dxmax*100,dy,dy/dxmax*100,dz,dz/dxmax*100);
fprintf('wavelength=%gm(%gdx)\n',wvl,wvl/dx);

npmlx = 8;
npmly = npmlx;
npmlz = npmlx;

nx = round(X/dx) + 2 * npmlx;
ny = round(Y/dy) + 2 * npmly;
nz_air = 2 + npmlz;
nz = round(Z/dz) + nz_air;
% dx/y/z,dt,nt,npmlx/y/z,nz_air, shold be same with prep_3d


% time prarameter
disp('<Time Grid>')
T_ref = sqrt(X^2+Y^2+(2*Z)^2)/(3e8/sqrt(epr_bg));%s
fprintf('Max 2-way traveltime=%gns\n',T_ref/1e-9)

dtmax = finddt(epr_min, miur_min, dx, dy, dz);
% dt = 0.9*dtmax/2;
% nt = ceil(T_ref/dt);
T = dt*nt;%s
fprintf('t_total=%gns.dt=%gns,nt=%g\n',T/1e-9,dt/1e-9,nt);
fprintf('Period=%gns(%gdt)\n',1e9/freq_src,1/(dt*freq_src));

check_dispersion(dx,dt,f_thr,c/sqrt(epr_max*miur_max));
check_stability(min([dx,dy,dz]),dt,c,1);

outstep_x_wavefield = 1;
fprintf('outstep_x_wavefield=%g\n',outstep_x_wavefield)
outstep_t_wavefield = 4;
outstep_slice = outstep_t_wavefield;

disp('↑↑↑↑↑↑↑↑↑↑Grid Parameter End↑↑↑↑↑↑↑↑↑↑')

%% model
%% background model
ep_bg = ones(nx,ny,nz) * epr_bg;
ep_bg(:,:,1:nz_air) = 1; % air layer

ep = ep_bg;

%% slice, src&rec
slicex_x = X/2;
slicex = round(slicex_x/dx+npmlx);
slicey_y = Y/2;
slicey = round((slicey_y)/dy+npmly);
% slicez_z = [0.1,0.5,1];
slicez_z = [0.9];
slicez = round(slicez_z/dz+nz_air);


%%
srclenth = 4*1/freq_src;
[src_wl, src_wl_t, ~, ~] = ricker_src(dt,freq_src,0.85,srclenth);
src_wl_rtm = src_wl;
src_wl_rtm_t = src_wl_t;

figure(13)
clf;
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])

src_t_min = src_wl_rtm_t(1);
if src_wl_t(1) < src_t_min
    src_t_min = src_wl_t(1);
end
src_t_min = src_t_min/1e-9;

subplot(2,1,2)
plot(src_wl_rtm_t/1e-9,src_wl_rtm)
xlim([src_t_min,-src_t_min]);
title('RTM source wavelet')
xlabel('t(ns)')
subplot(2,1,1)
plot(src_wl_t/1e-9,src_wl)
xlim([src_t_min,-src_t_min]);
title('real source wavelet')

export_fig(sprintf("%s_src.png",fn),'-transparent')

%% src and rec para
if zmode
    dy_src = Y/2;
    dx_src = X/2;
    src_margin_nx = 0.5/dx;
    src_margin_ny = 0.5/dy;
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
    ny_src = round((ny-4*src_margin_ny)/dny_src);
    if mod(ny_src,2) == 0
        ny_src = ny_src-1;
    end
    ny_rec = round((ny-4*rec_margin_ny)/dny_rec);
    if mod(ny_rec,2) == 0
        ny_rec = ny_rec-1;
    end
    srcy = ((-floor(ny_src/2):floor(ny_src/2)) * dy_src) + ny/2*dy;
    recy = ((-floor(ny_rec/2):floor(ny_rec/2)) * dy_rec) + ny/2*dy;
    srcz = 0-dz;
    recz = srcz;
else
    dx_src = 0.4;
    dx_rec = 0.2;
    dnx_src = dx_src / dx;
    dnx_rec = dx_rec / dx;
    nx_src = round((nx-2*(round(dx_src/dx)+npmlx))/dnx_src);
    if iseven(nx_src)
        nx_src = nx_src-1;
    end
    nx_rec = round((nx-2*round(dx_rec/dx+npmlx))/dnx_rec);
    if iseven(nx_rec)
        nx_rec = nx_rec-1;
    end
    srcx = ((-floor(nx_src/2):floor(nx_src/2)) * dx_src) + nx/2*dx;
    srcy = srcx;
    recx = ((-floor(nx_rec/2):floor(nx_rec/2)) * dx_rec) + nx/2*dx;
    recy = recx;
    srcz = 0-dz;
    recz = srcz;
end

%%%useless
dx_src=5;
dy_src=5;
dx_rec=0.1;
dy_rec=0.5;
src_margin_nx=npmlx;
src_margin_ny=npmly;
rec_margin_nx=npmlx;
rec_margin_ny=npmly;
src_span = 0;
%%%

[Xs,Ys,Zs] = meshgrid(srcx,srcy,srcz);
srcx = reshape(Xs,1,[]);
srcy = reshape(Ys,1,[]);
[Xr,Yr,Zr] = meshgrid(recx,recy,recz);
% if length(srcy) == 1
%     for j = 1:length(recy)
%         if recy(j) == srcy(1)
%             for i = 1:length(recx)
%                 if abs(srcx(1)-recx(i)) <= 0.2
%                     Xr(j,i) = nan;
%                     Yr(j,i) = nan;
%                     Zr(j,i) = nan;
%                 end
%             end
%         end
%     end
% end
recx = reshape(Xr,1,[]);
recx = recx(~isnan(recx));
recy = reshape(Yr,1,[]);
recy = recy(~isnan(recy));


srcinfo = round([srcx/dx+npmlx;srcy/dy+npmly;srcz*ones(size(srcx))/dz+nz_air]);
recinfo = round([recx/dx+npmlx;recy/dy+npmly;recz*ones(size(recx))/dz+nz_air]);
%%

save(fn_save,'filename','modelname','dx','dy','dz','nx','ny','nz','nz_air',...
    'npmlx','npmly','npmlz','T','dt','freq_src',...
    'outstep_t_wavefield','outstep_x_wavefield','outstep_slice',...
    'slicex','slicey','slicez','ep_bg','ep','srcinfo','recinfo',...
    'src_wl', 'src_wl_t','src_wl_rtm', 'src_wl_rtm_t');
%     'dx_src','dx_rec','dy_src','dy_rec','srcx','srcy','srcz','recx','recy','recz',...
%     'src_margin_nx','src_margin_ny','rec_margin_nx','rec_margin_ny','src_span',...

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

% %% 3D-view
% figure(10)
% clf;
% set(gcf,'Unit','centimeters')
% set(gcf,'Position',[0,0,29.7,21])
% set(gca,'fontsize',30,'fontname','Times')
% x = ((1:nx)-npmlx)*dx;
% y = ((1:ny)-npmlx)*dy;
% z = ((1:nz)-nz_air)*dz;
% [X,Y,Z] = meshgrid(y,x,z);
% % p0 = patch(isosurface(Y,X,Z,ep,epr1));
% % isonormals(X,Y,Z,ep,p0)
% % p0.FaceColor = 'black';
% % p0.EdgeColor = 'none';
% p1 = patch(isosurface(Y,X,Z,ep,epr_bg-0.1));
% isonormals(X,Y,Z,ep,p1)
% p1.FaceColor = 'black';
% p1.EdgeColor = 'none';
% % p3 = patch(isosurface(X,Y,Z,ep,ep0-3));
% % isonormals(X,Y,Z,ep,p3)
% % p3.FaceColor = 'magenta';
% % p3.EdgeColor = 'none';
% p5 = patch(isosurface(Y,X,Z,ep,11.9));
% isonormals(X,Y,Z,ep,p5)
% p5.FaceColor = 'black';
% p5.EdgeColor = 'none';
% daspect([1,1,1])
% view(3); alpha(.3);axis tight
% set(gca,'ZDir','reverse')
% set(gca,'YDir','reverse')
% xlim([x(npmlx),x(end-npmlx)]);
% ylim([y(npmly),y(end-npmly)]);
% zlim([z(npmlz),z(end-npmlz)]);
% % title('layers with fault')
% set(gca,'fontsize',20);
% xlabel('x(m)','Fontsize',24);ylabel('y(m)','Fontsize',24);zlabel('depth(m)','Fontsize',24);
% camlight('headlight') 
% lighting gouraud
% camlight('right')
% lighting gouraud
% 
% hold on
% if zmode
%     plot3(Xs,Ys,Zs,'r.')
% else
% plot3(Xs,Ys,Zs,'r^')
% plot3(Xr,Yr,Zr,'b.')
% end
% % plot3(Xr2,Yr2,Zr2,'b.')
% % xlim([0 10]);ylim([0 10]);zlim([-0.5 5])
% hold off
% set(gca,'fontsize',20,'fontname','Times')
% %%
% export_fig(fig_save,'-transparent')


% %%
% fsrc = fopen('src.in_0000','w');
% fprintf(fsrc,'%d,%d,%d\n',nrec,span,length(tt_m));
% for i =1:nrec
%     for ixspan = (1:xspan)-1
%         for iyspan = -yspan_half:yspan_half
%             fprintf(fsrc,'%d,%d,%d,Ey\n',recinfo(1,i)+ixspan,recinfo(2,i)+iyspan,recinfo(3,i));
%         end
%     end
% end
% for i =1:nrec
%     gather = interp1(tt_zc,gatherdata(:,i),tt_m,'linear',0);
%     gatherdatam(:,i) = gather;
%     fprintf(fsrc,'%g ',gather(end:-1:1));
%     fprintf(fsrc,'\n');
% end
% fclose(fsrc);
% %%
% start = 0;
% for i = 1:length(pdata.y)
%     ind{i} = (1:length(pdata.x{i})) + start;
%     start = start + length(pdata.x{i});
%     figure(i);imagesc(pdata.x{i},tt_m,gatherdatam(:,ind{i}))
%     pause(0.1)
% end


%%
function [wavelet,t_wavelet,f_thr,width] = ricker_src(dt,f,thr,tlength)
    disp('<Ricker Wavelet>')
    if(nargin<3)
       thr=0.85;
    end
    if(nargin<4)
       tlength=4*1/f;
    end
    fs = 1/dt;
    [~,t_wavelet] = ricker(dt,f,tlength);
    wavelet = (1 - 2 * (pi*f.*t_wavelet).^2) .* exp(-(pi*f.*t_wavelet).^2);
    figure(1);subplot(2,1,1)
    plot(t_wavelet,wavelet)
    width = 2*0.88521/(2*pi*f);
    fprintf('width:%gns\n',width/1e-9);
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
    disp('↓↓↓↓↓↓Dispersion Check↓↓↓↓↓↓')
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
    disp('↑↑↑↑↑↑Dispersion Check End↑↑↑↑↑↑')
end

function alpha = check_stability(dx, dt, cmax, dim)
    disp('↓↓↓↓↓↓Stability Check↓↓↓↓↓↓')
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
    disp('↑↑↑↑↑↑Stability Check End↑↑↑↑↑↑')
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