clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输出尖刺校正及未矫正图像
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(gcp('nocreate'));
path=['C:\DataFiles\Data_files\turkeyearthquake\'];
second_folder = sprintf('%s/%s', path, 'Results4');%制作app后results可编辑
results_filename=strcat(path,'Results4','\');
mkdir(second_folder);
namelist=dir([path,'*.asc']);
wordname='all6';
word_filename=strcat(wordname,'.doc');
filespec_user = [results_filename wordname];
try
    % 若Word服务器已经打开，返回其句柄Word
    Word = actxGetRunningServer('Word.Application');
catch
    % 否则，创建一个Microsoft Word服务器，返回句柄Word
    Word = actxserver('Word.Application');
end
Word.Visible = 1; % 或set(Word, 'Visible', 1);
% 若测试文件存在，打开该测试文件，否则，新建一个文件，并保存，文件名为测试.doc
if exist(filespec_user,'file')
    Document = Word.Documents.Open(filespec_user);
    % Document = invoke(Word.Documents,'Open',filespec_user);
else
    Document = Word.Documents.Add;
    % Document = invoke(Word.Documents, 'Add');
    Document.SaveAs2(filespec_user);
end
Content = Document.Content;
Selection = Word.Selection;
Paragraphformat = Selection.ParagraphFormat;
Document.PageSetup.TopMargin = 60;
Document.PageSetup.BottomMargin = 45;
Document.PageSetup.LeftMargin = 45;
Document.PageSetup.RightMargin = 45;
headline = wordname;
Content.Start = 0; % 起始点为0，即表示每次写入覆盖之前资料
Content.Text = headline;
Content.Font.Size = 10; % 字体大小
Content.Font.Bold = 1; % 字体加粗
Content.Paragraphs.Alignment = 'wdAlignParagraphCenter'; % 居中,wdAlignParagraphLeft/Center/Rig
parpool(8);
l=length(namelist);
A=cell(l,1);
parfor i=1:l
    fullname2=[path,namelist(i).name];
    A{i,1}=importdata(fullname2);
end
B={};
countb=1;
for i=1:l
    L=length(A{i,1}.textdata);
    for j=1:L
        tline2=A{i,1}.textdata(j);
        if strncmp(tline2,'PGA_CM/S^2',8)
            pga = regexp(tline2,'\d*\.?\d*','match');
            gbit0=str2num(pga{1,1}{1,2});
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
    if gbit0>=50
    BB=A{i,1}.data';
    BBB=BB(:);
    z=find(isnan(BBB));
    BBB(z)=0;
    B{countb,1}={BBB,namelist(i).name,gbit1,gbit2,gbit3,gbit0};
    countb=countb+1;
    end
end
% mm=1;
% ll=length(C);
% B={};
% for i=1:ll
%     pp=C{1,i}{1,6};
%     if  pp>=50%查找str中是否有pattern，返回出现位置，没有出现返回空数组
%         B{mm,1}=C{i,1};
%         mm=mm+1;
%     end
% end
dataname={};
count=1;
for w=1:length(B)
    data=B{w,1}{1,1};%数据
    result_name=B{w,1}{1,2};%用于新建文件夹命名
    pre_event_time=20;%事件发生前时间
    sp=B{w,1}{1,5};%采样周期
    pva=B{w,1}{1,6};
    a=length(data);
    c=pre_event_time/sp+1;
    base_line_start=mean(data(1:pre_event_time/sp+1));%基线计算
    data_start=data(pre_event_time/sp+2:a)-base_line_start;
    %     data_start=400.*data_start./max(abs(data_start));%数据前处理
    b=length(data_start);
    t=zeros(b,1);%时间
    wp=[0.01/sp 0.5/sp-1];%自动调节
    order=2;
    type='bandpass';
    [c,d]=butter(order,2*wp*sp,type);%带通滤波
    data_start=filter(c,d,data_start);
    X=data_start;
    for i=1:b
        t(i)=(i-1)*sp;
    end
    %     v=cumtrapz(data_start)*sp;
    %     d=cumtrapz(v)*sp;
    %     dd=d/max(abs(d));
    YY=sort(abs(data_start));
    if max(YY)>50
        dataname{count,1}=result_name;
        count=count+1;
        Paragraphformat.Alignment = 'wdAlignParagraphCenter'; % 居中
        Word.Selection.TypeParagraph;
        Selection.Start = Content.end; % 开始的地方在上一个的结尾
        Selection.Text = result_name; % 当前时间作为输出
        Selection.Font.Size = 12; % 字号
        Selection.MoveDown; %将所选内容向下移动，并返回移动距位数离的单
        a0=max(YY(ceil(length(YY)*0.75)),25);
        %     a0=40;
        dnum=2;
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     firstderaccel = diff(data_start)/sp; % First derivative of record.
        %     signchange = firstderaccel(1:end-1).*firstderaccel(2:end);
        %     secndderaccel = diff(data_start,2)/sp^2; % Second derivative of record.
        %     iextreme1 = find(firstderaccel==0)+1; % Finds flat portions in record.
        %     iextreme2 = find(signchange < 0)+1; % Finds changes of sign in first derivative.
        %     iextreme = [iextreme1' iextreme2'];
        %     iextreme = sort(iextreme);
        %     iextremez=[];
        %     iextremef=[];
        %     data_startz=[];
        %     data_startf=[];
        %     countz=1;
        %     countf=1;
        %     for j=1:length(iextreme)
        %         x=data_start(iextreme(j));
        %         if x>=0
        %             iextremez(countz,1)=iextreme(j);
        %             countz=countz+1;
        %         else
        %             iextremef(countf,1)=iextreme(j);
        %             countf=countf+1;
        %         end
        %     end
        %     for k=dnum/2+1:length(iextremez)-dnum/2
        %         xz0=data_start(iextremez(k));
        %         if abs(xz0)>=a0
        %             xz=zeros(dnum,1);
        %             bz=zeros(dnum,1);
        %             for l=1:dnum/2
        %                 xz(l)=data_start(iextremez(k-l));
        %                 xz(l+dnum/2)=data_start(iextremez(k+l));
        %                 bz(l)=abs(log(abs(xz(l)/xz0)));
        %                 bz(l+dnum/2)=abs(log(abs(xz(l+dnum/2)/xz0)));
        %             end
        %             theta=log(50)/log(4);
        %             bmin=min(bz);
        %             %bzmean=mean(bz);
        %             fz0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(xz0/a0)));
        %             fz=(1-fz0)*exp(-bz)/sum(exp(-bz));
        %             xz0_new=fz0*xz0+fz'*xz;
        %         else
        %             xz0_new=xz0;
        %         end
        %         data_startz(k-dnum/2,1)=xz0_new;
        %     end
        %     for k=dnum/2+1:length(iextremef)-dnum/2
        %         xf0=data_start(iextremef(k));
        %         if abs(xf0)>=a0
        %             xf=zeros(dnum,1);
        %             bf=zeros(dnum,1);
        %             for l=1:dnum/2
        %                 xf(l,1)=data_start(iextremef(k-l));
        %                 xf(l+dnum/2,1)=data_start(iextremef(k+l));
        %                 bf(l,1)=abs(log(abs(xf(l)/xf0)));
        %                 bf(l+dnum/2,1)=abs(log(abs(xf(l+dnum/2)/xf0)));
        %             end
        %             theta=log(50)/log(4);
        %             bmin=min(bf);
        %             %bfmean=mean(bf);
        %             ff0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(xf0/a0)));
        %             ff=(1-ff0)*(exp(-bf))/sum(exp(-bf));
        %             xf0_new=ff0*xf0+ff'*xf;
        %         else
        %             xf0_new=xf0;
        %         end
        %         data_startf(k-dnum/2,1)=xf0_new;
        %     end
        %     for k=dnum/2+1:length(iextremez)-dnum/2
        %         data_start(iextremez(k))=data_startz(k-dnum/2,1);
        %     end
        %     for k=dnum/2+1:length(iextremef)-dnum/2
        %         data_start(iextremef(k))=data_startf(k-dnum/2,1);
        %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     YY=sort(abs(data_start));
        %     a0=max(YY(ceil(length(YY)*0.85)),40);
        %     dnum=4;
        %     pv=mod(pva,10);
        %     pal=pva-pv+10;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        firstderaccel = diff(data_start)/sp; % First derivative of record.
        signchange = firstderaccel(1:end-1).*firstderaccel(2:end);
        secndderaccel = diff(data_start,2)/sp^2; % Second derivative of record.
        iextreme1 = find(firstderaccel==0)+1; % Finds flat portions in record.
        iextreme2 = find(signchange < 0)+1; % Finds changes of sign in first derivative.
        iextreme = [iextreme1' iextreme2'];
        iextreme = sort(iextreme);
        ddata=data_start(iextreme);
        tt=t(iextreme);
        XX=zeros(size(tt));
        maxomin = sign(secndderaccel); % Gets sign of second derivative.
        iextremez=[];
        iextremef=[];
        countz=1;
        countf=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j=1:length(iextreme)
            if(iextreme(j)>length(secndderaccel))
                break;
            else
                if ((maxomin(iextreme(j)-1) == -1) && data_start(iextreme(j)) > 0)
                    iextremez(countz,1)=iextreme(j);
                    countz=countz+1;
                else if ((maxomin(iextreme(j)-1) == 1) && data_start(iextreme(j)) < 0)
                        iextremef(countf,1)=iextreme(j);
                        countf=countf+1;
                        
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        end
        dzs = X(iextreme);
        %     for k=dnum/2+1:length(iextremez)-dnum/2
        %         data_start(iextremez(k))=data_startz(k-dnum/2,1);
        %     end
        %     for k=dnum/2+1:length(iextremef)-dnum/2
        %         data_start(iextremef(k))=data_startf(k-dnum/2,1);
        %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [peakdata,f]=(max(abs(dzs)));
        pvl=peakdata-mod(peakdata,10)+10;
        peak_t_data_satrt=tt(f);
        peakdata=sprintf('%8.2f',peakdata);
        peakdata=num2str(peakdata);
        peak_t_data_satrt=num2str(peak_t_data_satrt);
        [peakX,g]=(max(abs(ddata)));
        pdl=peakX-mod(peakX,10)+10;
        peak_t_X=tt(g);
        peakX=sprintf('%8.1f',peakX);
        peakX=num2str(peakX);
        peak_t_X=num2str(peak_t_X);
        figure('visible','off');
        set(gcf,'Units','centimeter','Position',[0 1 54 60]);
        ax1=subplot('Position',[3/54 33/60 50/54 25/60]);
        plot(tt,ddata,'-b','LineWidth',2);
        hold on;
        plot(tt,XX,'--r');
        legend(['Peak.raw: ',peakX,' cm·s^{-2} at ',peak_t_X,' s'],'zeroline','Location','southwest','fontname','Time New Roman','Box','off');
        axis([t(1) t(end) -pdl pdl]);
        ylabel('Acceleration(cm·s^{-2})','fontsize',16,'fontname','Time New Roman');
        xlabel('Time(s)','fontsize',16,'fontname','Time New Roman');
        ax2=subplot('Position',[3/54 4/60 50/54 25/60]);
        plot(tt,dzs,'-k','LineWidth',2);
        hold on;
        plot(tt,XX,'--r');
        legend(['Peak.cor: ',peakdata,' cm·s^{-2} at ',peak_t_data_satrt,' s'],'zeroline','Location','southwest','fontname','Time New Roman','Box','off');
        axis([t(1) t(end) -pdl pdl]);
        ylabel('Acceleration(cm·s^{-2})','fontsize',16,'fontname','Time New Roman');
        xlabel('Time(s)','fontsize',16,'fontname','Time New Roman');
        %     xlabel('Time(s)','fontsize',16,'fontname','Time New Roman');
        %     ax2=subplot('Position',[3/54 26/48 50/54 20/48]);
        %     plot(t,dd,'-b',LineWidth=3);
        %     hold on;
        %     plot(t,XX,'--r');
        %     pos2=axis(ax2);
        %     ylabel('Acceleration(cm·s^{-2})','position',[pos2(1)-pos2(2)/40,pos2(4)],'fontsize',16,'fontname','Time New Roman');
        %     xlabel('Time(s)','fontsize',16,'fontname','Time New Roman');
        %     legend(['Peak.raw: ',peakX,' cm·s^{-2} at ',peak_t_X,' s'],'zeroline','Location','southwest','fontname','Time New Roman','Box','off');
        %     axis([0 t(end) -1 1]);
        %     Selection.Start = Content.end;
        %     Selection=Word.Selection;
        %     print -dmeta
        %     invoke(Selection, 'Paste');
        %     Selection.MoveDown;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5傅氏谱
        % r=nextpow2(b);
        %     N=2^r;
        %     ddata_start=[data_start' zeros(1,N-b)];
        %     fy=fft(ddata_start)*sp;%快速傅里叶变换
        %     Fy=abs(fy(1:N/2));%只取低于采样频率一半的频率作图
        %     ffreq=(0:N/2-1)*(1/sp)/N;
        %     [peakX,g]=(max(abs(Fy)));
        %     pdl=peakX-mod(peakX,10)+10;
        %     peak_t_X=ffreq(g);
        %     peakX=sprintf('%8.1f',peakX);
        %     peakX=num2str(peakX);
        %     peak_t_X=num2str(peak_t_X);
        %     ax2=subplot('Position',[3/54 4/60 50/54 25/60]);
        %     loglog(ffreq,Fy);
        %     xlabel('频率/Hz');
        %     ylabel('F(f)');
        %     legend(['Peak.raw: ',peakX,' at ',peak_t_X,' Hz'],'zeroline','Location','southwest','fontname','Time New Roman','Box','off');
        %     axis([0 100 0 Inf]);%可编辑
        Selection.Start = Content.end;
        Selection=Word.Selection;
        print -dmeta
        invoke(Selection,'Paste');
        Selection.MoveDown;
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%小波系数图
        %     figure('visible','off')
        %     [wt,f,coi] = cwt(data_start,'bump',1/sp);
        %     pcolor(t,f,abs(wt));shading interp
        %     colormap(bone);
        %     set(gca,'ytick',[],'xtick',[],'xcolor','w','ycolor','w');
        %     set(gca, 'LooseInset', [0,0,0,0]);
        %     xlabel('时间 t/s');
        %     ylabel('频率 f/Hz');
        %     title('小波时频图');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    end
end
Document.Save; % 保存文档
Document.Close; % 关闭文档
Word.Quit; % 退出word服务器
%%%%%%%%%%%%%%%%%%%%%%%%
% excelname='all';
% excel_filename=strcat(excelname,'.xlsx');
% filespec_user = [results_filename excelname];
% try
%     Excel = actxGetRunningServer('Excel.Application');
% catch
%     Excel = actxserver('Excel.Application');
% end
% 
% Excel.Visible = 1;    % set(Excel, 'Visible', 1);
% 
% if exist(filespec_user,'file')
%     Workbook = Excel.Workbooks.Open(filespec_user);
% else
%     Workbook = Excel.Workbooks.Add;
%     Workbook.SaveAs(filespec_user);
% end
% 
% Sheets = Excel.ActiveWorkbook.Sheets;    % Sheets = Workbook.Sheets;
% Sheet1 = Sheets.Item(1);
% Sheet1.Activate;
% for n=1:length(dataname)
%     sl=num2str(n);
%     start_str=['A' sl];
%     Sheet1.Range(start_str).Value=dataname{n,1};
% end
% Workbook.Save; % 保存表格
% Workbook.Close; % 关闭表格
% Excel.Quit; % 退出Excel服务器