%%
mu0 = 4*pi*1e-7;
epsl0 = 8.85*1e-12;
dx = 0.02;
dy = 0.02;
dz = 0.02;

nx = 128;
ny=128;
nz=64;
nrec=25;
nt=2006;
xstep = 2;
tstep = 5;
%%
slx1={};
sly1={};
slz1={};
slx2={};
sly2={};
slz2={};
workdir1 = 'pstdtest_800MHz_2.0m_0.5m_fdtd_16';
workdir2 = 'pstdtest_800MHz_2.0m_0.5m_pstd_16';
xlist1 = dir(fullfile(workdir1,'Output', 'slx_Ey_0000*.bin'));
ylist1 = dir(fullfile(workdir1,'Output', 'sly_Ey_0000*.bin'));
zlist1 = dir(fullfile(workdir1,'Output', 'slz_Ey_0000*.bin'));
xlist2 = dir(fullfile(workdir2,'Output', 'slx_Ey_0000*.bin'));
ylist2 = dir(fullfile(workdir2,'Output', 'sly_Ey_0000*.bin'));
zlist2 = dir(fullfile(workdir2,'Output', 'slz_Ey_0000*.bin'));

%%
for i=1:length(xlist1)
    xname1 = fullfile(xlist1(i).folder,xlist1(i).name);
    yname1 = fullfile(ylist1(i).folder,ylist1(i).name);
    zname1 = fullfile(zlist1(i).folder,zlist1(i).name);
    
    fo = fopen(xname1);
    slx1{i} = fread(fo,[ny,nz],'single');
    fclose(fo);
    fo = fopen(yname1);
    sly1{i} = fread(fo,[nx,nz],'single');
    fclose(fo);
    fo = fopen(zname1);
    slz1{i} = fread(fo,[nx,ny],'single');
    fclose(fo);
end

for i=1:length(xlist2)
    xname2 = fullfile(xlist2(i).folder,xlist2(i).name);
    yname2 = fullfile(ylist2(i).folder,ylist2(i).name);
    zname2 = fullfile(zlist2(i).folder,zlist2(i).name);
    
    fo = fopen(xname2);
    slx2{i} = fread(fo,[ny,nz],'single');
    fclose(fo);
    fo = fopen(yname2);
    sly2{i} = fread(fo,[nx,nz],'single');
    fclose(fo);
    fo = fopen(zname2);
    slz2{i} = fread(fo,[nx,ny],'single');
    fclose(fo);
end

fo = fopen(fullfile(workdir2,'output','gather_0000_00000.bin'));
gather = fread(fo,[nrec,nt],'single');
fclose(fo);
%%
figure(10)
for i = 1:length(slx1)
    subplot(2,3,1)
    imagesc(slx1{i}(:,:)');colorbar
    title(['slicex fdtd at it=' num2str(i*tstep)]);xlabel('y');ylabel('z')
    subplot(2,3,2)
    imagesc(sly1{i}(:,:)');colorbar
    title(['slicex fdtd at it=' num2str(i*tstep)]);xlabel('x');ylabel('z')
    subplot(2,3,3)
    imagesc(slz1{i}');colorbar
    title(['slicex fdtd at it=' num2str(i*tstep)]);xlabel('x');ylabel('y')
    pause(0.1)
    subplot(2,3,4)
    imagesc(slx2{i}(:,:)');colorbar
    title(['slicex pstd at it=' num2str(i*tstep)]);xlabel('y');ylabel('z')
    subplot(2,3,5)
    imagesc(sly2{i}(:,:)');colorbar
    title(['slicex pstd at it=' num2str(i*tstep)]);xlabel('x');ylabel('z')
    subplot(2,3,6)
    imagesc(slz2{i}');colorbar
    title(['slicex pstd at it=' num2str(i*tstep)]);xlabel('x');ylabel('y')
    pause(0.1)
end

%%
figure(11)
for i = 1:nrec
    plot(gather(i,:))
    title(num2str(i))
    pause(0.1)
end

%% time:
np_max = 24;
np = 1:np_max;
n_txt = 6;

t_fdtd = 0;
t_cal_fdtd = 0;
t_io_fdtd = 0;

for i=1:n_txt
    fid = fopen(sprintf('time_fdtd%d.txt',i));
    c1 = '%*s %*s %*s %f\n';
    c2 = '%*s %*s %*s %f\n';
    c3 = '%*s %*s %f\n';
    formatSpec = [c1,c2,c3]; 
    A = textscan(fid, formatSpec,np_max);
    fclose(fid);
    t_fdtd = t_fdtd + A{3};
    t_cal_fdtd = t_cal_fdtd + A{2};
    t_io_fdtd = t_io_fdtd + A{1};
end
t_fdtd = t_fdtd / n_txt;
t_cal_fdtd = t_cal_fdtd / n_txt;
t_io_fdtd = t_io_fdtd / n_txt;

v_fdtd = t_fdtd(1)./t_fdtd;
v_cal_fdtd = t_cal_fdtd(1)./t_cal_fdtd;
v_io_fdtd = t_io_fdtd(1)./t_io_fdtd;

figure(21)
clf
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])
set(gca,'fontsize',30,'fontname','Times')
    
p = plot(np,t_fdtd,'-');
% p = plot(np,t_fdtd,'-',np,t_cal_fdtd,'--',np,t_io_fdtd,'-.');
% cr = p(2).Color;
% p(2).Color = p(1).Color; 
% p(3).Color = p(1).Color; 
% title('Time usage VS cores/threads','Fontsize',36);
xlabel('num of cores','Fontsize',20);ylabel('time(s)','Fontsize',20)
hold on
yyaxis right
pr = plot(np,v_fdtd);
% pr = plot(np,v_fdtd,'-',np,v_cal_fdtd,'--',np,v_io_fdtd,'-.');
% pr(1).Color = cr; 
% pr(2).Color = cr; 
% pr(3).Color = cr; 
ylabel('calculation speed','Fontsize',20)
% xlim([0,20])
% ylim(xlim/2)
hold off
% legend('Total time','Calculation time','I/O time','Total speed','Calculation speed','I/O speed','Location','north')
set(gca,'fontsize',20,'fontname','Times')%,'ycolor',cr)
yyaxis left
set(gca,'ycolor',p(1).Color)

%%
export_fig("time vs cores of FDTD.png")

%%
fid = fopen('time_pstd.txt');
c1 = 'Total time: %fs, Total calculate time: %fs, Total I/O time: %fs\n';
formatSpec = [c1]; 
B = textscan(fid, formatSpec);
fclose(fid);
t_pstd = B{1};
t_cal_pstd = B{2};
t_io_pstd = B{3};
%%
figure(21)
semilogy(np,t_fdtd,'r',np,t_cal_fdtd,'r--',np,t_io_fdtd,'r-.',np,t_pstd,'b-',np,t_cal_pstd,'b--',np,t_io_pstd,'b-.')
set(gca,'fontsize',24);
title('Time usage VS cores/threads','Fontsize',36);xlabel('num of cores/threads','Fontsize',32);ylabel('time(s)','Fontsize',32)
legend('Total time(FDTD)','Calculation time(FDTD)','I/O time(FDTD)','Total time(PSTD)','Calculation time(PSTD)','I/O time(PSTD)')
export_fig("time vs cores(threads)_all_in_one.png")
%%
figure(20)
subplot(2,2,1)
plot(np,t_cal_fdtd)
title('Total calculate time of FDTD');xlabel('num of cores');ylabel('time(s)')
set(gca,'fontsize',18);
subplot(2,2,2)
plot(np,t_cal_pstd)
title('Total calculate time of PSTD');xlabel('num of threads');ylabel('time(s)')
set(gca,'fontsize',18);
subplot(2,2,3)
plot(np,t_io_fdtd)
title('Total I/O time of FDTD');xlabel('num of cores');ylabel('time(s)')
set(gca,'fontsize',18);
subplot(2,2,4)
plot(np,t_io_pstd)
title('Total I/O time of PSTD');xlabel('num of threads');ylabel('time(s)')
set(gca,'fontsize',18);
export_fig("time vs cores(threads).png")