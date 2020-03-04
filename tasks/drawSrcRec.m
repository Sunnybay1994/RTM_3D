modelname='EM300';

nx = 300;
dx = 0.02;
npmlx = 12;
dx_src = 0.6;
dx_rec = 0.2;


pic=imread([modelname '.png']);
m=squeeze(pic(:,:,1));
m(m~=255)=1;
m(m==255)=0.001;
figure(10)
x=(1:300)*dx;
imagesc(x,x,m);colorbar;title('model');xlabel('x(m)');ylabel('y(m)');
hold on

%% plot src
dnx_src = dx_src / dx;
dnx_rec = dx_rec / dx;
nx_src = floor((300-2*npmlx)/dnx_src) -1;
nx_rec = floor((300-2*npmlx)/dnx_rec) -1;

figure(10)

srcx = (npmlx+dnx_src:dnx_src:(nx-npmlx-dnx_src))*0.02;
srcy = srcx;
srcx = repmat(srcx,[nx_src,1]);
srcx = reshape(srcx,[nx_src^2,1]);
srcy = repmat(srcy,[1,nx_src]);
plot(srcx,srcy,'r^','markersize',6)
saveas(gcf,fullfile('Zero_offset',[modelname '_s.png']));
hold on

recx = (npmlx+dnx_rec:dnx_rec:(nx-npmlx-dnx_rec))*0.02;
recy = recx;
recx = repmat(recx,[nx_rec,1]);
recx = reshape(recx,[nx_rec^2,1]);
recy = repmat(recy,[1,nx_rec]);
plot(recx,recy,'cx')
hold off

saveas(gcf,[modelname '_sr.png']);