workdir = 'checkbox1.0_0.4_0.2';
ab_type = 'eps';

nx = 100;
ny = 100;
nz = 70;
nz_air = 10;
dx = 0.1;
dy = 0.1;
dz = 0.1;
order = 2;
x = dx: dx : (dx * nx);
y = (1:ny) * dy;
z = ((1:nz) - nz_air) * dz;

% %% plot src & rec
% nz0=20;
% dx = 0.02;
% 
% npmlx = 12;
% dx_src = 0.2;
% dy_src = 0.5;
% 
% dnx_src = dx_src / dx;
% dny_src = dy_src / dx;
% % nx_src = floor((300-2*npmlx)/dnx_src);
% % ny_src = floor((300-2*npmlx)/dny_src);
% nx_src = 25
% ny_src = 7
% if iseven(nx_src)
%     nx_src = nx_src - 1;
% end
% if iseven(ny_src)
%     ny_src = ny_src - 1;
% end
% 
% srcx = ((-(dnx_src*(nx_src-1)/2):dnx_src:(dnx_src*(nx_src-1)/2)) + nx/2)*dx;
% srcy = ((-(dny_src*(ny_src-1)/2):dny_src:(dny_src*(ny_src-1)/2)) + ny/2)*dx;
% srcx = repmat(srcx,[ny_src,1]);
% srcx = reshape(srcx,[nx_src*ny_src,1]);
% srcy = repmat(srcy,[1,nx_src]);
% plot(srcx,srcy,'xb','markersize',16)
% hold on
% 
% zslice = [0];
% x = (1:nx) * 0.02;
% y = (1:ny) * 0.02;
% z = ((1:nz) - 20)*0.02;
% slice(x,y,z,permute(eps,[2,1,3]),[],[],zslice,'nearest')

%%
fn_infos = dir(fullfile(workdir,['Input/' ab_type '*']));
folder = fn_infos(1).folder;
eps = [];
for i = 1:length(fn_infos)
    fn = fullfile(folder,[ab_type '.in_' num2str(i-1,'%03d')]);
    disp(fn)
    dat = load(fn);
    dat = permute(reshape(dat,ny,[],nz),[2,1,3]);
%     if i == 1
%         dat = dat(1:end-2,:,:);
%     elseif i == length(fn_infos)
%         dat = dat(3:end,:,:);
%     else
%         dat = dat(3:end-2,:,:);
%     end
    eps = cat(1,eps,dat(3:end-2,:,:));
end

%% 2d slice view
% for i = 1:size(eps,3)
%     imagesc(squeeze(eps(:,:,i)));colorbar;
%     title(['y_n=' num2str(i)])
%     pause(0.1)
% end
for i = 15:30:45
    figure()
    imagesc(x,z,squeeze(eps(:,i,:))');colorbar;
    title(['yslice at y=' num2str(i/10) 'm'])
    xlabel('x(m)');ylabel('depth(m)')
    saveas(gca,['yslice at y=' num2str(i/10) 'm.png'])
    imagesc(x,y,squeeze(eps(:,:,(i+nz_air)))');colorbar;
    title(['zslice at z=' num2str(i/10) 'm'])
    xlabel('x(m)');ylabel('y(m)');
    saveas(gca,['zslice at z=' num2str(i/10) 'm.png'])
end

% %% 3d view
% ind1 = find(eps > 10 );
% % ind2 = find(eps == 4);
% for k=1:length(ind1)
%     [I1(k),J1(k),K1(k)] = ind2sub(size(eps),ind1(k));
% end
% I1 = I1*0.02;J1=J1*0.02;K1=((K1-20)*0.02);
% % for k=1:length(ind2)
% %     [I2(k),J2(k),K2(k)] = ind2sub(size(eps),ind2(k));
% % end
% % I2 = I2*0.02;J2=J2*0.02;K2=((K2-20)*0.02);
% 
% plot3(I1,J1,K1,'.r')%,I2,J2,K2,'.r')
% xlabel('x(m)');ylabel('y(m)');zlabel('depth(m)');
% xlim([0,6]);ylim([0,4.8]);zlim([-0.4,2])
% set(gca,'zdir','reverse','FontSize',24)

%% 3D slice view
x = dx: dx : (dx * nx);
y = (1:ny) * dy;
z = ((1:nz) - nz_air) * dz;
result = eps;
figure(21)
slice(x,y,z,result,[],[10],[1.5,4.5]);
shading interp;alpha(0.6);
colorbar;daspect([1 1 1]);
% caxis([6,12])
zlim([-1,6]);xlim([0,10]);ylim([0,10]);
set(gca,'fontsize',20);
xlabel('x(m)');ylabel('y(m)');zlabel('depth(m)');
set(gca,'zdir','reverse','FontSize',24)
