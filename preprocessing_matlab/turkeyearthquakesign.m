clear all
clc
xi=0.05;%阻尼比
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%对单条数据进行处理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path=['C:\DataFiles\Data_files\turkeyearthquake\'];
namelist1=dir([path,'20230206011734_2104_ap_RawAcc_U.asc']);
% namelist2=dir([path,'20230206102447_0122_mp_RawAcc_E.asc']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fullname1=[path,namelist1.name];
% fullname2=[path,namelist2.name];
%%%%%%%%%%%%%%%%%%%%%%%%%%5
A=importdata(fullname1);
% C=importdata(fullname2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=length(A.textdata);
for j=1:L
    tline2=A.textdata(j);
    if strncmp(tline2,'PGA_CM/S^2',8)
        pga0 = regexp(tline2,'\d*\.?\d*','match');
        gbit0=str2num(pga0{1,1}{1,2});
    end
    if strncmp(tline2,'STATION_LATITUDE_DEGREE',23)
        latitude = regexp(tline2,'\d*\.?\d*','match');
        gbit1=str2num(latitude{1,1}{1,1});
    end
    if strncmp(tline2,'STATION_LONGITUDE_DEGREE',24)
        longitude = regexp(tline2,'\d*\.?\d*','match');
        gbit2=str2num(longitude{1,1}{1,1});
    end
    if strncmp(tline2,'SAMPLING_INTERVAL_S',19)
        sampt = regexp(tline2,'\d*\.?\d*','match');
        gbit3=str2num(sampt{1,1}{1,1});
    end
end
BB=A.data';
BBB=BB(:);
z=find(isnan(BBB));
BBB(z)=0;
B={BBB,namelist1.name,gbit1,gbit2,gbit3,gbit0};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M=length(C.textdata);
% for j=1:M
%     tline2=C.textdata(j);
%     if strncmp(tline2,'PGA_CM/S^2',8)
%         pga0 = regexp(tline2,'\d*\.?\d*','match');
%         gbit0=str2num(pga0{1,1}{1,2});
%     end
%     if strncmp(tline2,'STATION_LATITUDE_DEGREE',23)
%         latitude = regexp(tline2,'\d*\.?\d*','match');
%         gbit1=str2num(latitude{1,1}{1,1});
%     end
%     if strncmp(tline2,'STATION_LONGITUDE_DEGREE',24)
%         longitude = regexp(tline2,'\d*\.?\d*','match');
%         gbit2=str2num(longitude{1,1}{1,1});
%     end
%     if strncmp(tline2,'SAMPLING_INTERVAL_S',19)
%         sampt = regexp(tline2,'\d*\.?\d*','match');
%         gbit3=str2num(sampt{1,1}{1,1});
%     end
% end
% DD=C.data';
% DDD=DD(:);
% u=find(isnan(DDD));
% DDD(u)=0;
% D={DDD,namelist2.name,gbit1,gbit2,gbit3,gbit0};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
data=B{1,1};%数据
result_name=B{1,2};%用于新建文件夹命名
pre_event_time=20;%事件发生前时间
sp=B{1,5};%采样周期
pa=B{1,6};
a=length(data);
c=pre_event_time/sp+1;
base_line_start=mean(data(1:pre_event_time/sp+1));%基线计算
data_start=data(pre_event_time/sp+2:a)-base_line_start;%数据前处理
wp=[0.025/sp 0.5/sp-1];%自动调节
order=2;
type='bandpass';
[c,d]=butter(order,2*wp*sp,type);%带通滤波
data_start=filter(c,d,data_start);
X=zeros(size(data_start));
b=length(data_start);
t=zeros(b,1);%时间
for i=1:b
    t(i)=(i-1)*sp;
end
YY=sort(abs(data_start));
a0=max(YY(ceil(length(YY)*0.8)),25);
dnum=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstderaccel = diff(data_start)/sp; % First derivative of record.
signchange = firstderaccel(1:end-1).*firstderaccel(2:end);
secndderaccel = diff(data_start,2)/sp^2; % Second derivative of record.
iextreme1 = find(firstderaccel==0)+1; % Finds flat portions in record.
iextreme2 = find(signchange < 0)+1; % Finds changes of sign in first derivative.
iextreme = [iextreme1' iextreme2'];
iextreme = sort(iextreme);
maxomin = sign(secndderaccel); % Gets sign of second derivative.
% Creates vectors ‘maximo’ and ‘minimo’ with local maximum and minimum values.
iextremez=[];
iextremef=[];
countz=1;
countf=1;
for j=1:length(iextreme)
    if(iextreme(j)>length(secndderaccel))
        break;
    else
        if ((maxomin(iextreme(j)-1) == -1) && (data_start(iextreme(j)) > 0))
            iextremez(countz,1)=iextreme(j);
            countz=countz+1;
        else if ((maxomin(iextreme(j)-1) == 1) && (data_start(iextreme(j)) < 0))
                iextremef(countf,1)=iextreme(j);
                countf=countf+1;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=dnum/2+1:length(iextremez)-dnum/2
    xz0=data_start(iextremez(k));
    if abs(xz0)>=a0
        xz=zeros(dnum,1);
        bz=zeros(dnum,1);
        for l=1:dnum/2
            xz(l)=data_start(iextremez(k-l));
            xz(l+dnum/2)=data_start(iextremez(k+l));
            bz(l)=abs(log(abs(xz(l)/xz0)));
            bz(l+dnum/2)=abs(log(abs(xz(l+dnum/2)/xz0)));
        end
        theta=log(50)/log(2);
        bmin=min(bz);
        %bzmean=mean(bz);
        fz0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(xz0/a0)));
        fz=(1-fz0)*exp(-bz)/sum(exp(-bz));
        xz0_new=fz0*xz0+fz'*xz;
    else
        xz0_new=xz0;
    end
    X(iextremez(k))=xz0_new;
end
for k=dnum/2+1:length(iextremef)-dnum/2
    xf0=data_start(iextremef(k));
    if abs(xf0)>=a0
        xf=zeros(dnum,1);
        bf=zeros(dnum,1);
        for l=1:dnum/2
            xf(l,1)=data_start(iextremef(k-l));
            xf(l+dnum/2,1)=data_start(iextremef(k+l));
            bf(l,1)=abs(log(abs(xf(l)/xf0)));
            bf(l+dnum/2,1)=abs(log(abs(xf(l+dnum/2)/xf0)));
        end
        theta=log(50)/log(2);
        bmin=min(bf);
        %bfmean=mean(bf);
        ff0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(xf0/a0)));
        ff=(1-ff0)*(exp(-bf))/sum(exp(-bf));
        xf0_new=ff0*xf0+ff'*xf;
    else
        xf0_new=xf0;
    end
    X(iextremef(k))=xf0_new;
    s=[iextremez;iextremef];
    s=sort(s);
    s=s';
    datagoal=data_start;
    dzs = X(s);
    dl=length(dzs);
    tgoal=t;
    c1=dl;
    C1{i,1}=datagoal;
    C1{i,2}=dzs;
    C1{i,3}=tgoal;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data2=D{1,1};%数据
% sp2=D{1,5};%采样周期
% a2=length(data);
% c2=pre_event_time/sp2+1;
% base_line_start2=mean(data2(1:pre_event_time/sp2+1));%基线计算
% data_start2=data2(pre_event_time/sp2+2:a2)-base_line_start2;%数据前处理
% data_start2=data_start2(1:50/sp2);%数据前处理
% b2=length(data_start2);
% t2=zeros(b2,1);%时间
% for i=1:b2
%     t2(i)=(i-1)*sp2;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [peakdata,f]=(max(abs(data_start)));
% pal=peakdata-mod(peakdata,10)+10;
% [peakdata2,f2]=(max(abs(data_start2)));
% pal2=peakdata2-mod(peakdata2,10)+10;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('visible','on');
% set(gcf,'Units','centimeter','Position',[0 1 54 60]);
% ax1=subplot('Position',[3/54 33/60 50/54 25/60],'XAxisLocation','bottom','YAxisLocation','left');
% plot(t,data_start,'-k','LineWidth',2);
% legend('20230206102447-0131-mp-N','Location','southeast','fontname','Times New Roman','fontsize',18,'Box','off');
% axis([t(1) t(end) -pal pal])
% set(gca,'YTick',[-pal:200:pal],'fontname','Times New Roman','fontsize',14,'linewidth',1) %改变y轴坐标间隔显示
% box off
% ax12=axes('Position',[3/54 33/60 50/54 25/60],'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k','linewidth',1);
% set(ax12,'YTick', [],'XTick', []);
% ax2=subplot('Position',[3/54 4/60 50/54 25/60]);
% plot(t,data_start2,'-b','LineWidth',2);
% legend('20230206102447-0127-mp-E','Location','southeast','fontname','Times New Roman','fontsize',18,'Box','off');
% axis([t(1) t(end) -pal2 pal2]);
% pos=axis(ax2);
% ylabel('Acceleration(cm·s^{-2})','position',[pos(1)-pos(2)/30,pos(4)+pos(4)/10],'fontsize',24,'fontname','Times New Roman');
% xlabel('Time(s)','fontsize',20,'fontname','Times New Roman');
% set(gca,'YTick',[-pal2:30:pal2],'fontname','Times New Roman','fontsize',14,'linewidth',1) %改变y轴坐标间隔显示
% box off
% ax22=axes('Position',[3/54 4/60 50/54 25/60],'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k','linewidth',1);
% set(ax22,'YTick', [],'XTick', []);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('visible','on');
% set(gcf,'Units','centimeter','Position',[0 1 54 31]);
% ax3=subplot('Position',[3/54 4/31 50/54 25/31]);
% plot(t,data_start,'-k','linewidth',2);
% hold on
% plot(t,data_start2,'-b','linewidth',2);
% legend('20230206102447-0131-mp-N','20230206102447-0127-mp-E','Location','southeast','fontsize',18,'fontname','Times New Roman','Box','off');
% ylabel('Acceleration(cm·s^{-2})','fontsize',24,'fontname','Times New Roman');
% xlabel('Time(s)','fontsize',20,'fontname','Times New Roman');
% axis([t(1) t(end) -pal pal]);
% set(gca,'YTick',[-pal:200:pal],'fontname','Times New Roman','fontsize',14,'linewidth',1) %改变y轴坐标间隔显示
% box off
% ax32=axes('Position',[3/54 4/31 50/54 25/31],'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k','linewidth',1);
% set(ax32,'YTick', [],'XTick', []);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wp=[0.025/sp 0.5/sp];%自动调节
% order=2;
% type='bandpass';
% [c,d]=butter(order,wp*sp,type);%带通滤波
% X=filter(c,d,data_start);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[peakdata,f]=(max(abs(data_start)));
pal=peakdata-mod(peakdata,10)+10;
peak_t_data_satrt=t(f);
peakdata=sprintf('%8.2f',peakdata);
peakdata=num2str(peakdata);
peak_t_data_satrt=num2str(peak_t_data_satrt);
[peakX,g]=(max(abs(X)));
pdl=peakX-mod(peakX,10)+10;
peak_t_X=t(g);
peakX=sprintf('%8.1f',peakX);
peakX=num2str(peakX);
peak_t_X=num2str(peak_t_X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('visible','on');
set(gcf,'Units','centimeter','Position',[0 1 54 60]);
ax1=subplot('Position',[3/54 33/60 50/54 25/60],'XAxisLocation','bottom','YAxisLocation','left');
plot(t,data_start,'-k','LineWidth',2);
legend(['Peak.raw: ',peakdata,' cm·s^{-2} at ',peak_t_data_satrt,' s'],'Location','southeast','fontname','Times New Roman','fontsize',18,'Box','off');
axis([t(1) t(end) -pal pal])
set(gca,'YTick',[-pal:20:pal],'fontname','Times New Roman','fontsize',14,'linewidth',1) %改变y轴坐标间隔显示
box off
ax12=axes('Position',[3/54 33/60 50/54 25/60],'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k','linewidth',1);
set(ax12,'YTick', [],'XTick', []);
ax2=subplot('Position',[3/54 4/60 50/54 25/60]);
plot(t,X,'-b','LineWidth',2);
legend(['Peak.fil: ',peakX,' cm·s^{-2} at ',peak_t_X,' s'],'Location','southeast','fontname','Times New Roman','fontsize',18,'Box','off');
axis([t(1) t(end) -pal pal]);
pos=axis(ax2);
ylabel('Acceleration(cm·s^{-2})','position',[pos(1)-pos(2)/30,pos(4)+pos(4)/10],'fontsize',24,'fontname','Times New Roman');
xlabel('Time(s)','fontsize',20,'fontname','Times New Roman');
set(gca,'YTick',[-pal:20:pal],'fontname','Times New Roman','fontsize',14,'linewidth',1) %改变y轴坐标间隔显示
box off
ax22=axes('Position',[3/54 4/60 50/54 25/60],'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k','linewidth',1);
set(ax22,'YTick', [],'XTick', []);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
second_folder = sprintf('%s/%s', path, 'Results2');%制作app后results可编辑
results_filename=strcat(path,'Results2','\');
mkdir(second_folder);
for i=1:1
    wordname=result_name(1:end-4);
    word_filename=strcat(wordname,'.xls');
    filespec_user = [results_filename wordname];
    try
        Excel = actxGetRunningServer('Excel.Application');
    catch
        Excel = actxserver('Excel.Application');
    end
    
    Excel.Visible = 1;    % set(Excel, 'Visible', 1);
    
    if exist(filespec_user,'file')
        Workbook = Excel.Workbooks.Open(filespec_user);
    else
        Workbook = Excel.Workbooks.Add;
        Workbook.SaveAs(filespec_user);
    end
    
    Sheets = Excel.ActiveWorkbook.Sheets;    % Sheets = Workbook.Sheets;
    Sheet1 = Sheets.Item(1);
    Sheet1.Activate;
    lengthd=length(t);
    for n=1:lengthd
        sl=num2str(n);
        start_str1=['A' sl];
        start_str2=['B' sl];
        start_str3=['C' sl];
        Sheet1.Range(start_str1).Value=t(n);
        Sheet1.Range(start_str2).Value=data_start(n);
        Sheet1.Range(start_str3).Value=X(n);
    end
    Workbook.Save; % 保存表格
    Workbook.Close; % 关闭表格
    Excel.Quit; % 退出Excel服务器
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pga=max(abs(data_start));
% data_start=data_start*0.6*980/pga;
% dd=[sp;data_start];
% cname=strcat(final_name,result_name,'.txt');
% fid = fopen(cname,'wt');
% fprintf(fid,'%g\n',dd);
% vx=(cumtrapz(data_start)*sp);
% dx=(cumtrapz(data_start)*sp);
% tt=zeros(c,1);
% for j=1:c
%     tt(j)=(j-1)*sp;
% end
% data_prevent=data(1:c)-base_line_start;
% s=nextpow2(c);
% M=2^s;
% ddata_prevent=[data_prevent' zeros(1,M-c)];
% fpy=fft(ddata_prevent)*sp;%快速傅里叶变换
% Fpy=abs(fpy(1:M/2))/(sp*M);%只取低于采样频率一半的频率作图
% ffpreq=(0:M/2-1)*(1/sp)/M;
% k=find(max(Fpy));
% frl=ffpreq(k(1));
% if frl>=0.1
%     lt=0.1;
% elseif frl<=0.05
%     lt=0.05;
% else
%     lt=frl;
% end
% m=1;
% m=num2str(m);
% XX=zeros(size(t));
% v=mod(10,pa);
% pal=pa-v+10;
% wp=[lt 35];%自动调节
% order=2;
% type='bandpass';
% [c,d]=butter(order,wp*sp,type);%带通滤波
% X=filter(c,d,data_start);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perlow=0.05;%反应谱绘制最低周期
% perupp=20;%反应谱绘制最高周期
% t_freq_ver=perupp-perlow;
% ver=t_freq_ver*100;%周期分割
% ys=zeros(size(t));%存放所计算相对位移
% yv=zeros(size(t));%存放所计算相对速度
% Xya=zeros(size(t));%存放所计算绝对加速度
% ysmax=zeros(ver,1);%位移反应谱
% yvmax=zeros(ver,1);%速度反应谱
% Xyamax=zeros(ver,1);%加速度反应谱
% delta_t=sp;%采样时间间隔
% for m=0:ver
%     omega=2*pi/(t_freq_ver*m/ver+0.05);
% a11=exp(-xi*omega*delta_t)*(cos(sqrt(-xi^2 + 1)*omega*delta_t) + xi*sin(sqrt(-xi^2 + 1)*omega*delta_t)/sqrt(-xi^2 + 1));
% a12=exp(-xi*omega*delta_t)*sin(sqrt(-xi^2 + 1)*omega*delta_t)/(omega*sqrt(-xi^2 + 1));
% a21=omega*exp(xi*omega*(-delta_t))*sin((-delta_t)*omega*sqrt(-xi^2 + 1))/sqrt(-xi^2 + 1);
% a22=-xi*exp(-xi*omega*(delta_t))*sin(sqrt(-xi^2 + 1)*omega*(delta_t))/sqrt(-xi^2 + 1) + exp(-xi*omega*(delta_t))*cos(sqrt(-xi^2 + 1)*omega*(delta_t));
% b11=exp(-xi*omega*(delta_t))*((1/omega^2 + 2*xi/((delta_t)*omega^3))*cos(sqrt(-xi^2 + 1)*omega*(delta_t)) + (xi/omega + (2*xi^2 - 1)/((delta_t)*omega^2))*sin(sqrt(-xi^2 + 1)*omega*(delta_t))/(sqrt(-xi^2 + 1)*omega)) - 2*xi/((delta_t)*omega^3);
% b12=exp(-xi*omega*(delta_t))*(-2*xi*cos(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)*omega^3) - (2*xi^2 - 1)*sin(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)*omega^3*sqrt(-xi^2 + 1))) + 2*xi/((delta_t)*omega^3) - 1/omega^2;
% b21=-xi*omega*exp(-xi*omega*(delta_t))*((1/omega^2 + 2*xi/((delta_t)*omega^3))*cos(sqrt(-xi^2 + 1)*omega*(delta_t)) + (xi/omega + (2*xi^2 - 1)/((delta_t)*omega^2))*sin(sqrt(-xi^2 + 1)*omega*(delta_t))/(sqrt(-xi^2 + 1)*omega)) + exp(-xi*omega*(delta_t))*(-2*xi*cos(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)^2*omega^3) - (1/omega^2 + 2*xi/((delta_t)*omega^3))*sqrt(-xi^2 + 1)*omega*sin(sqrt(-xi^2 + 1)*omega*(delta_t)) - (2*xi^2 - 1)*sin(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)^2*omega^3*sqrt(-xi^2 + 1)) + (xi/omega + (2*xi^2 - 1)/((delta_t)*omega^2))*cos(sqrt(-xi^2 + 1)*omega*(delta_t))) + 2*xi/((delta_t)^2*omega^3);
% b22=-xi*omega*exp(-xi*omega*(delta_t))*(-2*xi*cos(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)*omega^3) - (2*xi^2 - 1)*sin(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)*omega^3*sqrt(-xi^2 + 1))) + exp(-xi*omega*(delta_t))*(2*xi*cos(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)^2*omega^3) + 2*xi*sqrt(-xi^2 + 1)*sin(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)*omega^2) + (2*xi^2 - 1)*sin(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)^2*omega^3*sqrt(-xi^2 + 1)) - (2*xi^2 - 1)*cos(sqrt(-xi^2 + 1)*omega*(delta_t))/((delta_t)*omega^2)) - 2*xi/((delta_t)^2*omega^3);
% for j=1:b-1
%     ys(j+1)=a11*ys(j)+a12*yv(j)+b11*data_start(j)+b12*data_start(j+1);
%     yv(j+1)=a21*ys(j)+a22*yv(j)+b21*data_start(j)+b22*data_start(j+1);
% end
% for s=1:b
%     Xya(s)=-(2*xi*omega*yv(s)+omega^2*ys(s));
% end
% ysmax(m+1)=max(abs(ys));
% yvmax(m+1)=max(abs(yv));
% Xyamax(m+1)=max(abs(Xya));
% end
% s_ysmax=ysmax/max(dx);
% s_yvmax=yvmax/max(vx);
% s_Xyamax=Xyamax/max(data_start);
% Data_acca=readtable([final_name,'1.csv']);
% Data_accv=readtable([final_name,'2.csv']);
% Data_accd=readtable([final_name,'3.csv']);
% [row,col]=size(Data_acca);
% period=table2array(Data_acca(3:row,1));
% vwa=table2array(Data_acca(3:row,2))/max(abs(data_start));
% vwv=table2array(Data_accv(3:row,2))/max(abs(vx));
% vwd=table2array(Data_accd(3:row,2))/max(abs(dx));
% t_omega_v=perlow:t_freq_ver/ver:perupp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,1,1);
% loglog(t_omega_v,s_Xyamax,period,vwa);
% xlabel('周期T/s');
% ylabel('绝对加速度反应谱值/cm·s^{-2}');
% legend('精确解','viewwave解');
% title('绝对加速度反应谱');
% subplot(3,1,2);
% loglog(t_omega_v,s_yvmax,period,vwv);
% xlabel('周期T/s');
% ylabel('相对速度反应谱值/cm·s^{-1}');
% legend('精确解','viewwave解');
% title('相对速度反应谱');
% subplot(3,1,3);
% loglog(t_omega_v,s_ysmax,period,vwd);
% xlabel('周期T/s');
% ylabel('相对位移反应谱值/cm');
% legend('精确解','viewwave解');
% title('相对位移反应谱');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [peakdata,f]=(max(abs(data_start)));
% peak_t_data_satrt=t(f);
% peakdata=sprintf('%8.1f',peakdata);
% peakdata=num2str(peakdata);
% peak_t_data_satrt=num2str(peak_t_data_satrt);
% [peakX,g]=(max(abs(X)));
% peak_t_X=t(g);
% peakX=sprintf('%8.1f',peakX);
% peakX=num2str(peakX);
% peak_t_X=num2str(peak_t_X);
% figure('visible','on');
% subplot(2,1,1);
% set(gcf,'Units','centimeter','Position',[0 0 50 20])
% plot(t,data_start,'-k');
% hold on;
% plot(t,XX,'--r');
% ylabel('Acceleration(cm·s^{-2})','position',[-6.5,-700],'fontsize',16,'fontname','Times New Roman');
% legend(['Peak.raw: ',peakdata,' cm·s^{-2} at ',peak_t_data_satrt,' s'],'zeroline','Location','southwest','fontname','Times New Roman','Box','off');
% axis([t(1) t(end) -pal pal]);
% subplot(2,1,2);
% set(gcf,'Units','centimeter','Position',[0 0 50 20]);
% plot(t,X,'-b');
% hold on;
% plot(t,XX,'--r');
% xlabel('Time(s)','position',[46.245,-700],'fontsize',16,'fontname','Times New Roman');
% legend(['Peak.cor: ',peakX,' cm·s^{-2} at ',peak_t_X,' s'],'zeroline','Location','southwest','fontname','Times New Roman','Box','off');
% axis([0 t(end) -pal pal]);
% saveas(gcf,[final_name,[m,'.png']]);
%'position',[-15,3.5],
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r=nextpow2(b);
% N=2^r;
% ddata_start=[data_start' zeros(1,N-b)];
% fy=fft(ddata_start)*sp;%快速傅里叶变换
% Fy=abs(fy(1:N/2))/(sp*N);%只取低于采样频率一半的频率作图
% ffreq=(0:N/2-1)*(1/sp)/N;
% figure('visible','off');
% set(gcf,'Units','centimeter','Position',[5 5 50 20]);
% loglog(ffreq,Fy);
% axis([0 100 0 Inf]);%可编辑
% xlabel('频率/Hz');
% ylabel('F(f)');
% title([result_name,'傅氏谱(滤波前)']);
% saveas(gcf,[final_name,'傅氏谱(未滤波).png']);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('visible','off');
% set(gcf,'Units','centimeter','Position',[5 5 50 20]);
% V=cumtrapz(t,data_start);
% VV=abs(V);
% [pv,w]=(max(VV));
% v=mod(10,pv);
% pvl=pv-v+10;
% plot(t,V);
% hold on;
% plot(t,XX);
% xlabel('时间/s');
% ylabel('地面速度值/cm·s^{-2}');
% axis([0 t(end) -pvl pvl]);
% title([result_name,'时间历程图(滤波前)']);
% saveas(gcf,[final_name,'速度时间历程图(未滤波).png']);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('visible','off');
% set(gcf,'Units','centimeter','Position',[5 5 50 20]);
% W=cumtrapz(t,X);
% plot(t,W);
% hold on;
% plot(t,XX);
% xlabel('时间/s');
% ylabel('地面速度值/cm·s^{-2}');
% axis([0 t(end) -pvl pvl]);
% title([result_name,'时间历程图(滤波后)']);
% saveas(gcf,[final_name,'速度时间历程图(已滤波).png']);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xx=[X' zeros(1,N-b)];
% bfy=fft(Xx)*sp;%快速傅里叶变换
% Fby=abs(bfy(1:N/2))/(sp*N);%只取低于采样频率一半的频率作图
% fbfreq=(0:N/2-1)*(1/sp)/N;
% figure('visible','off');
% set(gcf,'Units','centimeter','Position',[5 5 50 20]);
% loglog(fbfreq,Fby);
% axis([0 100 0 Inf]);%可编辑
% xlabel('频率/Hz');
% ylabel('F(f)');
% title([result_name,'傅氏谱(滤波后)']);
% saveas(gcf,[final_name,'傅氏谱(已滤波).png']);

% Document.Save; % 保存文档
% Document.Close; % 关闭文档
% Word.Quit; % 退出word服务器