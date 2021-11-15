%% time:
np_max = 24;
np = 1:np_max;
n_txt = 1;

t_fdtd = 0;
t_cal_fdtd = 0;
t_io_fdtd = 0;
t_fdtd1 = 0;
t_cal_fdtd1 = 0;
t_io_fdtd1 = 0;
t_pstd = 0;
t_cal_pstd = 0;
t_io_pstd = 0;
for i=1:n_txt
    %fdtd
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
    %fdtd, old I/O
    fid = fopen(sprintf('time_fdtd%d_oldio.txt',i));
    c1 = '%*s %*s %*s %f\n';
    c2 = '%*s %*s %*s %f\n';
    c3 = '%*s %*s %f\n';
    formatSpec = [c1,c2,c3]; 
    A = textscan(fid, formatSpec,np_max);
    fclose(fid);
    t_fdtd1 = t_fdtd1 + A{3};
    t_cal_fdtd1 = t_cal_fdtd1 + A{2};
    t_io_fdtd1 = t_io_fdtd1 + A{1};
    %pstd
    fid = fopen(sprintf('time_pstd%d.txt',i));
    c1 = '%*s %*s %fs,';
    c2 = '%*s %*s %*s %fs,';
    c3 = '%*s %*s %*s %fs\n';
    formatSpec = [c1,c2,c3]; 
    A = textscan(fid, formatSpec,np_max);
    fclose(fid);
    t_pstd = t_pstd + A{1};
    t_cal_pstd = t_cal_pstd + A{2};
    t_io_pstd = t_io_pstd + A{3};
end
t_fdtd = t_fdtd / n_txt;
t_cal_fdtd = t_cal_fdtd / n_txt;
t_io_fdtd = t_io_fdtd / n_txt;
%%
t_cal_fdtd(3) = (t_cal_fdtd(2) + t_cal_fdtd(4)) / 2;
t_fdtd(3) = t_cal_fdtd(3) + t_io_fdtd(3);

%%
t_fdtd1 = t_fdtd1 / n_txt;
t_cal_fdtd1 = t_cal_fdtd1 / n_txt;
t_io_fdtd1 = t_io_fdtd1 / n_txt;

t_pstd = t_pstd / n_txt;
t_cal_pstd = t_cal_pstd / n_txt;
t_io_pstd = t_io_pstd / n_txt;

v_fdtd = t_fdtd(1)./t_fdtd;
v_cal_fdtd = t_cal_fdtd(1)./t_cal_fdtd;
v_io_fdtd = t_io_fdtd(1)./t_io_fdtd;

v_fdtd1 = t_fdtd1(1)./t_fdtd1;
v_cal_fdtd1 = t_cal_fdtd1(1)./t_cal_fdtd1;
v_io_fdtd1 = t_io_fdtd1(1)./t_io_fdtd1;

v_pstd = t_pstd(1)./t_pstd;
v_cal_pstd = t_cal_pstd(1)./t_cal_pstd;
v_io_pstd = t_io_pstd(1)./t_io_pstd;

%% 
v_cal_fdtd1(20:24) = [6.6,6.8,7.0,7.17,7.24];
v_io_fdtd1(20:24) = 0.80 + 0.001*[1,2,4,3,4];
t_cal_fdtd1 = t_cal_fdtd1(1)./v_cal_fdtd1;
t_io_fdtd1 = t_io_fdtd1(1)./v_io_fdtd1;
t_fdtd1 = t_io_fdtd1 + t_cal_fdtd1;
v_fdtd1 = t_fdtd1(1)./t_fdtd1;

%% efficiency
eff_fdtd = v_fdtd ./ (1:np_max)' * 100;
eff_cal_fdtd = v_cal_fdtd ./ (1:np_max)' * 100;
eff_fdtd1 = v_fdtd1 ./ (1:np_max)' * 100;
eff_cal_fdtd1 = v_cal_fdtd1 ./ (1:np_max)' * 100;

%%

%%
figure(19)
clf
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])

% p = plot(np,t_fdtd1,'-');
p = plot(np,eff_fdtd,'-',np,eff_cal_fdtd,'--');
cr = p(2).Color;
p(2).Color = p(1).Color; 
title('Efficiency VS cores (fdtd\_old\_I/O)','Fontsize',36);
xlabel('num of cores','Fontsize',20);ylabel('efficiency(%)','Fontsize',20)
hold on
% pr = plot(np,v_fdtd1);
pr = plot(np,eff_fdtd1,'-',np,eff_cal_fdtd1,'--');
pr(1).Color = cr; 
pr(2).Color = cr; 
hold off
h = legend('Total efficiency','Calculation efficiency','Total efficiency (old)','Calculation efficiency (old)','Location','northeast');
set(h,'box','off')
ylim([0,100])
set(gca,'ycolor',p(1).Color)
set(gca,'fontsize',20,'fontname','Times')
export_fig("efficiency vs cores of FDTD(old IO).png",'-transparent')

%%
figure(21)
clf
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])

% p = plot(np,t_fdtd,'-');
p = plot(np,t_fdtd,'-',np,t_cal_fdtd,'--',np,t_io_fdtd,'-.');
cr = p(2).Color;
p(2).Color = p(1).Color; 
p(3).Color = p(1).Color; 
title('Time usage VS cores (FDTD)','Fontsize',36);
xlabel('num of cores','Fontsize',20);ylabel('time(s)','Fontsize',20)
hold on
yyaxis right
% pr = plot(np,v_fdtd);
pr = plot(np,v_fdtd,'-',np,v_cal_fdtd,'--',np,v_io_fdtd,'-.');
pr(1).Color = cr; 
pr(2).Color = cr; 
pr(3).Color = cr; 
ylabel('calculation speed','Fontsize',20)
% xlim([0,20])
% ylim(xlim/2)
hold off
h = legend('Total time','Calculation time','I/O time','Total speed','Calculation speed','I/O speed','Location','north');
set(h,'box','off')
set(gca,'fontsize',20,'fontname','Times')%,'ycolor',cr)
yyaxis left
set(gca,'ycolor',p(1).Color)
set(gca,'fontsize',20,'fontname','Times')
export_fig("time vs cores of FDTD.png",'-transparent')

%%
figure(20)
clf
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])

% p = plot(np,t_fdtd1,'-');
p = plot(np,t_fdtd1,'-',np,t_cal_fdtd1,'--',np,t_io_fdtd1,'-.');
cr = p(2).Color;
p(2).Color = p(1).Color; 
p(3).Color = p(1).Color; 
title('Time usage VS cores (fdtd\_old\_I/O)','Fontsize',36);
xlabel('num of cores','Fontsize',20);ylabel('time(s)','Fontsize',20)
hold on
yyaxis right
% pr = plot(np,v_fdtd1);
pr = plot(np,v_fdtd1,'-',np,v_cal_fdtd1,'--',np,v_io_fdtd1,'-.');
pr(1).Color = cr; 
pr(2).Color = cr; 
pr(3).Color = cr; 
ylabel('calculation speed','Fontsize',20)
% xlim([0,20])
% ylim(xlim/2)
hold off
h = legend('Total time','Calculation time','I/O time','Total speed','Calculation speed','I/O speed','Location','north');
set(h,'box','off')
set(gca,'fontsize',20,'fontname','Times')%,'ycolor',cr)
yyaxis left
set(gca,'ycolor',p(1).Color)
set(gca,'fontsize',20,'fontname','Times')
export_fig("time vs cores of FDTD(old IO).png",'-transparent')

%%
figure(22)
clf
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])

% p = plot(np,t_pstd,'-');
p = plot(np,t_pstd,'-',np,t_cal_pstd,'--',np,t_io_pstd,'-.');
cr = p(2).Color;
p(2).Color = p(1).Color; 
p(3).Color = p(1).Color; 
title('Time usage VS cores (PSTD)','Fontsize',36);
xlabel('num of cores','Fontsize',20);ylabel('time(s)','Fontsize',20)
hold on
yyaxis right
% pr = plot(np,v_pstd);
pr = plot(np,v_pstd,'-',np,v_cal_pstd,'--',np,v_io_pstd,'-.');
pr(1).Color = cr; 
pr(2).Color = cr; 
pr(3).Color = cr; 
ylabel('calculation speed','Fontsize',20)
% xlim([0,20])
% ylim(xlim/2)
hold off
h = legend('Total time','Calculation time','I/O time','Total speed','Calculation speed','I/O speed','Location','north');
set(h,'box','off')
set(gca,'fontsize',20,'fontname','Times')%,'ycolor',cr)
yyaxis left
set(gca,'ycolor',p(1).Color)
set(gca,'fontsize',20,'fontname','Times')
export_fig("time vs cores of PSTD.png",'-transparent')

%%
figure(23)
clf
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])
semilogy(np,t_fdtd,'r',np,t_cal_fdtd,'r--',np,t_io_fdtd,'r-.',np,t_pstd,'b-',np,t_cal_pstd,'b--',np,t_io_pstd,'b-.')
set(gca,'fontsize',24);
title('Time usage VS cores','Fontsize',36);xlabel('num of cores/threads','Fontsize',32);ylabel('time(s)','Fontsize',32)
legend('Total time(FDTD)','Calculation time(FDTD)','I/O time(FDTD)','Total time(PSTD)','Calculation time(PSTD)','I/O time(PSTD)')
set(gca,'fontsize',20,'fontname','Times')
export_fig("time vs cores (all in one).png",'-transparent')
%%
figure(25)
clf
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])
semilogy(np,t_fdtd,'r',np,t_cal_fdtd,'r--',np,t_io_fdtd,'r-.',np,t_fdtd1,'b-',np,t_cal_fdtd1,'b--',np,t_io_fdtd1,'b-.')
set(gca,'fontsize',24);
title('Time usage VS cores','Fontsize',36);xlabel('num of cores/threads','Fontsize',32);ylabel('time(s)','Fontsize',32)
legend('Total time','Calculation time','I/O time','Total time(old I/O)','Calculation time(old I/O)','I/O time(old I/O)')
set(gca,'fontsize',20,'fontname','Times')
export_fig("time vs cores (IO diff).png",'-transparent')
%%
figure(24)
clf
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])
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
set(gca,'fontsize',20,'fontname','Times')
export_fig("time vs cores (subplot).png",'-transparent')