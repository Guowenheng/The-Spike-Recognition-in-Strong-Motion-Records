clear all
clc
delete(gcp('nocreate'));
path=['C:\DataFiles\Data_files\turkeyearthquake\'];
second_folder = sprintf('%s/%s', path, 'Results_2');%制作app后results可编辑
results_filename=strcat(path,'Results_2','\');
mkdir(second_folder);
[num,txt,raw]=xlsread([path,'all.xlsx']);
namelist1={};
namelist0={};
% namelist2={};
r=length(num);
aa=1;
bb=1;
cll=0.05:0.05:0.8;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% cc=1;
for i=1:r
    if raw{i,5}==0
        namelist0{aa}=raw{i,1};
        aa=aa+1;
    elseif raw{i,5}==1
        namelist1{bb}=raw{i,1};
        bb=bb+1;
%     elseif raw{i,5}==3
%         namelist2{cc}=raw{i,1};
%         cc=cc+1;
    end
end
parpool(8);
l0=length(namelist0);
A0=cell(l0,1);
parfor i=1:l0
    fullname0=[path,namelist0{i}];
    A0{i,1}=importdata(fullname0);
end
l1=length(namelist1);
A1=cell(l1,1);
parfor i=1:l1
    fullname1=[path,namelist1{i}];
    A1{i,1}=importdata(fullname1);
end
% l2=length(namelist2);
% A2=cell(l2,1);
% for i=1:l2
%     fullname2=[path,namelist2{i}];
%     A2{i,1}=importdata(fullname2);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:l0
    L=length(A0{i,1}.textdata);
    for j=1:L
        tline0=A0{i,1}.textdata(j);
        if strncmp(tline0,'PGA_CM/S^2',8)
            pga = regexp(tline0,'\d*\.?\d*','match');
            gbit0=str2num(pga{1,1}{1,2});
        end
        if strncmp(tline0,'STATION_LATITUDE_DEGREE',23)
            latitude = regexp(tline0,'\d*\.?\d*','match');
            gbit1=str2num(latitude{1,1}{1,1});
        end
        if strncmp(tline0,'STATION_LONGITUDE_DEGREE',24)
            longitude = regexp(tline0,'\d*\.?\d*','match');
            gbit2=str2num(longitude{1,1}{1,1});
        end
        if strncmp(tline0,'SAMPLING_INTERVAL_S',19)
            sampt = regexp(tline0,'\d*\.?\d*','match');
            gbit3=str2num(sampt{1,1}{1,1});
        end
    end
    gbit4=erase(namelist0(i),'.asc');
    BB=A0{i,1}.data';
    BBB=BB(:);
    z=find(isnan(BBB));
    BBB(z)=0;
    B0{i,1}={BBB,gbit4,gbit1,gbit2,gbit3,gbit0};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:l1
    L=length(A1{i,1}.textdata);
    for j=1:L
        tline1=A1{i,1}.textdata(j);
        if strncmp(tline1,'PGA_CM/S^2',8)
            pga = regexp(tline1,'\d*\.?\d*','match');
            gbit0=str2num(pga{1,1}{1,2});
        end
        if strncmp(tline1,'STATION_LATITUDE_DEGREE',23)
            latitude = regexp(tline1,'\d*\.?\d*','match');
            gbit1=str2num(latitude{1,1}{1,1});
        end
        if strncmp(tline1,'STATION_LONGITUDE_DEGREE',24)
            longitude = regexp(tline1,'\d*\.?\d*','match');
            gbit2=str2num(longitude{1,1}{1,1});
        end
        if strncmp(tline1,'SAMPLING_INTERVAL_S',19)
            sampt = regexp(tline1,'\d*\.?\d*','match');
            gbit3=str2num(sampt{1,1}{1,1});
        end
    end
    gbit4=erase(namelist1(i),'.asc');
    BB=A1{i,1}.data';
    BBB=BB(:);
    z=find(isnan(BBB));
    BBB(z)=0;
    B1{i,1}={BBB,gbit4,gbit1,gbit2,gbit3,gbit0};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% for i=1:l2
%     L=length(A2{i,1}.textdata);
%     for j=1:L
%         tline2=A2{i,1}.textdata(j);
%         if strncmp(tline2,'PGA_CM/S^2',8)
%             pga = regexp(tline2,'\d*\.?\d*','match');
%             gbit0=str2num(pga{1,1}{1,2});
%         end
%         if strncmp(tline2,'STATION_LATITUDE_DEGREE',23)
%             latitude = regexp(tline2,'\d*\.?\d*','match');
%             gbit1=str2num(latitude{1,1}{1,1});
%         end
%         if strncmp(tline2,'STATION_LONGITUDE_DEGREE',24)
%             longitude = regexp(tline2,'\d*\.?\d*','match');
%             gbit2=str2num(longitude{1,1}{1,1});
%         end
%         if strncmp(tline2,'SAMPLING_INTERVAL_S',19)
%             sampt = regexp(tline2,'\d*\.?\d*','match');
%             gbit3=str2num(sampt{1,1}{1,1});
%         end
%     end
%     gbit4=erase(namelist2(i),'.asc');
%     BB=A2{i,1}.data';
%     BBB=BB(:);
%     z=find(isnan(BBB));
%     BBB(z)=0;
%     B2{i,1}={BBB,gbit4,gbit1,gbit2,gbit3,gbit0};
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for v=1:length(cll)
cc=num2str(cll(v));
results_filename2=strcat(results_filename,cc,'\');
normal=strcat('normal',cc);
abnormal=strcat('abnormal',cc);
C0=cell(l0,3);
c0=zeros(l0,1);
for i=1:l0
    data=B0{i,1}{1,1};
    pre_event_time=20;%事件发生前时间
    sp=B0{i,1}{1,5};%采样周期
    a=length(data);
    base_line_start=mean(data(1:pre_event_time/sp+1));%基线计算
    data_start=(data(pre_event_time/sp+2:a)-base_line_start);%数据前处理
    wp=[0.01/sp 0.5/sp-1];%自动调节%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    order=2;
    type='bandpass';
    [c,d]=butter(order,2*wp*sp,type);%带通滤波
    data_start=filter(c,d,data_start);
    X=zeros(size(data_start));
    b=length(data_start);
    t=zeros(b,1);%时间
    for j=1:b
        t(j)=(j-1)*sp;
    end
    YY=sort(abs(data_start));
    a0=max(YY(ceil(length(YY)*0.8)),50);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ccl=a0*cll(v);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    iextremez1=[];
    iextremef1=[];
    countz=1;
    countf=1;
    countz1=1;
    countf1=1;
    for j=1:length(iextreme)
        if(iextreme(j)>length(secndderaccel))
            break;
        else
            if ((maxomin(iextreme(j)-1) == -1) && (data_start(iextreme(j)) > 0))% && (abs(data_start(iextreme(j)))>=ccl) )%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                iextremez(countz,1)=iextreme(j);
                countz=countz+1;
                if abs(data_start(iextreme(j)))>=ccl
                    iextremez1(countz1,1)=iextreme(j);
                    countz1=countz1+1;
                end
            elseif ((maxomin(iextreme(j)-1) == 1) && (data_start(iextreme(j)) < 0))% && (abs(data_start(iextreme(j)))>=ccl) )%%%%%%%%%%%%%%%%%%%%%%%%%5
                iextremef(countf,1)=iextreme(j);
                countf=countf+1;
                if abs(data_start(iextreme(j)))>=ccl
                    iextremef1(countf1,1)=iextreme(j);
                    countf1=countf1+1;
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=dnum/2+1:length(iextremez)-dnum/2
        xz0=data_start(iextremez(k));
        if abs(xz0)>=ccl%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            fz0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(2*xz0/a0)));
            fz=(1-fz0)*exp(-bz)/sum(exp(-bz));
            xz0_new=fz0*xz0+fz'*xz;
        else
            xz0_new=xz0;
        end
        X(iextremez(k))=xz0_new;
    end
    for k=dnum/2+1:length(iextremef)-dnum/2
        xf0=data_start(iextremef(k));
        if abs(xf0)>=ccl%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            ff0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(2*xf0/a0)));
            ff=(1-ff0)*(exp(-bf))/sum(exp(-bf));
            xf0_new=ff0*xf0+ff'*xf;
        else
            xf0_new=xf0;
        end
        X(iextremef(k))=xf0_new;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s=[iextremez1;iextremef1];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    s=sort(s);
    s=s';
    datagoal=data_start;
    dzs = X(s);
    dl=length(dzs);
    tgoal=t;
    c0(i,1)=dl;
    C0{i,1}=datagoal;
    C0{i,2}=dzs;
    C0{i,3}=tgoal;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1=cell(l1,3);
c1=zeros(l1,1);
for i=1:l1
    data=B1{i,1}{1,1};
    pre_event_time=20;%事件发生前时间
    sp=B1{i,1}{1,5};%采样周期
    a=length(data);
    base_line_start=mean(data(1:pre_event_time/sp+1));%基线计算
    data_start=(data(pre_event_time/sp+2:a)-base_line_start);%数据前处理
    wp=[0.01/sp 0.5/sp-1];%自动调节%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    order=2;
    type='bandpass';
    [c,d]=butter(order,2*wp*sp,type);%带通滤波
    data_start=filter(c,d,data_start);
    X=zeros(size(data_start));
    b=length(data_start);
    t=zeros(b,1);%时间
    for j=1:b
        t(j)=(j-1)*sp;
    end
    YY=sort(abs(data_start));
    a0=max(YY(ceil(length(YY)*0.8)),50);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    ccl=cll(v)*a0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    iextremez1=[];
    iextremef1=[];
    countz=1;
    countf=1;
    countz1=1;
    countf1=1;
    for j=1:length(iextreme)
        if(iextreme(j)>length(secndderaccel))
            break;
        else
            if ((maxomin(iextreme(j)-1) == -1) && (data_start(iextreme(j)) > 0))% && abs(data_start(iextreme(j)))>=ccl )%%%%%%%%%%%%%%%%5
                iextremez(countz,1)=iextreme(j);
                countz=countz+1;
                if abs(data_start(iextreme(j)))>=ccl
                    iextremez1(countz1,1)=iextreme(j);
                    countz1=countz1+1;
                end
            elseif ((maxomin(iextreme(j)-1) == 1) && (data_start(iextreme(j)) < 0))% && abs(data_start(iextreme(j)))>=ccl )%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                iextremef(countf,1)=iextreme(j);
                countf=countf+1;
                if abs(data_start(iextreme(j)))>=ccl
                    iextremef1(countf1,1)=iextreme(j);
                    countf1=countf1+1;
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=dnum/2+1:length(iextremez)-dnum/2
        xz0=data_start(iextremez(k));
        if abs(xz0)>=ccl%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
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
            fz0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(2*xz0/a0)));
            fz=(1-fz0)*exp(-bz)/sum(exp(-bz));
            xz0_new=fz0*xz0+fz'*xz;
        else
            xz0_new=xz0;
        end
        X(iextremez(k))=xz0_new;
    end
    for k=dnum/2+1:length(iextremef)-dnum/2
        xf0=data_start(iextremef(k));
        if abs(xf0)>=ccl%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            ff0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(2*xf0/a0)));
            ff=(1-ff0)*(exp(-bf))/sum(exp(-bf));
            xf0_new=ff0*xf0+ff'*xf;
        else
            xf0_new=xf0;
        end
        X(iextremef(k))=xf0_new;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    s=[iextremez1;iextremef1];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s=sort(s);
    s=s';
    datagoal=data_start;
    dzs = X(s);
    dl=length(dzs);
    tgoal=t;
    c1(i,1)=dl;
    C1{i,1}=datagoal;
    C1{i,2}=dzs;
    C1{i,3}=tgoal;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C2=cell(l2,3);
% c2=zeros(l2,1);
% for i=1:l2
%     data=B2{i,1}{1,1};
%     pre_event_time=20;%事件发生前时间
%     sp=B2{i,1}{1,5};%采样周期
%     a=length(data);
%     base_line_start=mean(data(1:pre_event_time/sp+1));%基线计算
%     data_start=(data(pre_event_time/sp+2:a)-base_line_start);%数据前处理
%     wp=[0.025/sp 0.5/sp-1];%自动调节
%     order=2;
%     type='bandpass';
%     [c,d]=butter(order,2*wp*sp,type);%带通滤波
%     data_start=filter(c,d,data_start);
%     X=zeros(size(data_start));
%     b=length(data_start);
%     t=zeros(b,1);%时间
%     for j=1:b
%         t(j)=(j-1)*sp;
%     end
%     YY=sort(abs(data_start));
%     a0=max(YY(ceil(length(YY)*0.75)),25);
%     dnum=2;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     firstderaccel = diff(data_start)/sp; % First derivative of record.
%     signchange = firstderaccel(1:end-1).*firstderaccel(2:end);
%     secndderaccel = diff(data_start,2)/sp^2; % Second derivative of record.
%     iextreme1 = find(firstderaccel==0)+1; % Finds flat portions in record.
%     iextreme2 = find(signchange < 0)+1; % Finds changes of sign in first derivative.
%     iextreme = [iextreme1' iextreme2'];
%     iextreme = sort(iextreme);
%     maxomin = sign(secndderaccel); % Gets sign of second derivative.
%     % Creates vectors ‘maximo’ and ‘minimo’ with local maximum and minimum values.
%     iextremez=[];
%     iextremef=[];
%     countz=1;
%     countf=1;
%     for j=1:length(iextreme)
%         if(iextreme(j)>length(secndderaccel))
%             break;
%         else
%             if ((maxomin(iextreme(j)-1) == -1) && (data_start(iextreme(j)) > 0))
%                 iextremez(countz,1)=iextreme(j);
%                 countz=countz+1;
%             else if ((maxomin(iextreme(j)-1) == 1) && (data_start(iextreme(j)) < 0))
%                     iextremef(countf,1)=iextreme(j);
%                     countf=countf+1;
%                 end
%             end
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%             theta=log(50)/log(2);
%             bmin=min(bz);
%             %bzmean=mean(bz);
%             fz0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(xz0/a0)));
%             fz=(1-fz0)*exp(-bz)/sum(exp(-bz));
%             xz0_new=fz0*xz0+fz'*xz;
%         else
%             xz0_new=xz0;
%         end
%         X(iextremez(k))=xz0_new;
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
%             theta=log(50)/log(2);
%             bmin=min(bf);
%             %bfmean=mean(bf);
%             ff0=(1-exp(-(theta*bmin)))/(1+exp(1-abs(xf0/a0)));
%             ff=(1-ff0)*(exp(-bf))/sum(exp(-bf));
%             xf0_new=ff0*xf0+ff'*xf;
%         else
%             xf0_new=xf0;
%         end
%         X(iextremef(k))=xf0_new;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     s=[iextremez;iextremef];
%     s=sort(s);
%     s=s';
%     datagoal=data_start;
%     dzs = X(s);
%     d2=length(dzs);
%     tgoal=t;
%     c2(i,1)=d2;
%     C2{i,1}=datagoal;
%     C2{i,2}=dzs;
%     C2{i,3}=tgoal;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
All=max([c0;c1]);
all=ceil(All/100)*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_name0=strcat(results_filename2,normal,'\');
mkdir(final_name0);
for i=1:l0
    excelname=num2str(i);
    excel_filename=strcat(excelname,'.xlsx');
    filespec_user = [final_name0 excelname];
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
    %     data1=C0{i,1};
    data2=C0{i,2};
    %     datafinal1=zeros(all,1);
    datafinal2=zeros(all,1);
    lengthdata=c0(i,1);
    %     st=B0{i,1}{1,5};
    %     t=C0{i,3};
    %     tfinal=[];%时间
    %     lt=t(end)/st;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for j=1:lengthdata
    %         datafinal1(j,1)=data1(j,1);
    %     end
    for j=1:lengthdata
        datafinal2(j,1)=data2(j,1);
    end
    %     for j=1:all-length(t)
    %         tfinal=[tfinal;(lt+1)*st];
    %         lt=lt+1;
    %     end
    %     tt=[t;tfinal];
    %     for n=1:all
    %         sl=num2str(n);
    %         start_str=['A' sl];
    %         Sheet1.Range(start_str).Value=datafinal1(n);
    %     end
    for n=1:all
        sl=num2str(n);
        start_str=['A' sl];
        Sheet1.Range(start_str).Value=datafinal2(n);
    end
    %     for n=1:all
    %         sl=num2str(n);
    %         start_str=['C' sl];
    %         Sheet1.Range(start_str).Value=tt(n);
    %     end
    Workbook.Save; % 保存表格
    Workbook.Close; % 关闭表格
    Excel.Quit; % 退出Excel服务器
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_name1=strcat(results_filename2,abnormal,'\');
mkdir(final_name1);
for i=1:l1
    excelname=num2str(i);
    excel_filename=strcat(excelname,'.xlsx');
    filespec_user = [final_name1 excelname];
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
    %     data1=C1{i,1};
    data2=C1{i,2};
    %     datafinal1=zeros(all,1);
    datafinal2=zeros(all,1);
    lengthdata=c1(i,1);
    %     st=B1{i,1}{1,5};
    %     t=C1{i,3};
    %     tfinal=[];%时间
    %     lt=t(end)/st;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     for j=1:lengthdata
    %         datafinal1(j,1)=data1(j,1);
    %     end
    for j=1:lengthdata
        datafinal2(j,1)=data2(j,1);
    end
    %     for j=1:all-length(t)
    %         tfinal=[tfinal;(lt+1)*st];
    %         lt=lt+1;
    %     end
    %     tt=[t;tfinal];
    %     for n=1:all
    %         sl=num2str(n);
    %         start_str=['A' sl];
    %         Sheet1.Range(start_str).Value=datafinal1(n);
    %     end
    for n=1:all
        sl=num2str(n);
        start_str=['A' sl];
        Sheet1.Range(start_str).Value=datafinal2(n);
    end
    %     for n=1:all
    %         sl=num2str(n);
    %         start_str=['B' sl];
    %         Sheet1.Range(start_str).Value=tt(n);
    %     end
    Workbook.Save; % 保存表格
    Workbook.Close; % 关闭表格
    Excel.Quit; % 退出Excel服务器0
end
% final_name2=strcat(results_filename,'other','\');
% mkdir(final_name2);
% for i=1:l2
%     excelname=num2str(i);
%     excel_filename=strcat(excelname,'.xlsx');
%     filespec_user = [final_name2 excelname];
%     try
%         Excel = actxGetRunningServer('Excel.Application');
%     catch
%         Excel = actxserver('Excel.Application');
%     end
%
%     Excel.Visible = 1;    % set(Excel, 'Visible', 1);
%
%     if exist(filespec_user,'file')
%         Workbook = Excel.Workbooks.Open(filespec_user);
%     else
%         Workbook = Excel.Workbooks.Add;
%         Workbook.SaveAs(filespec_user);
%     end
%
%     Sheets = Excel.ActiveWorkbook.Sheets;    % Sheets = Workbook.Sheets;
%     Sheet1 = Sheets.Item(1);
%     Sheet1.Activate;
%     %     data1=C0{i,1};
%     data2=C2{i,2};
%     %     datafinal1=zeros(all,1);
%     datafinal2=zeros(all,1);
%     lengthdata=c2(i,1);
%     %     st=B0{i,1}{1,5};
%     %     t=C0{i,3};
%     %     tfinal=[];%时间
%     %     lt=t(end)/st;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %     for j=1:lengthdata
%     %         datafinal1(j,1)=data1(j,1);
%     %     end
%     for j=1:lengthdata
%         datafinal2(j,1)=data2(j,1);
%     end
%     %     for j=1:all-length(t)
%     %         tfinal=[tfinal;(lt+1)*st];
%     %         lt=lt+1;
%     %     end
%     %     tt=[t;tfinal];
%     %     for n=1:all
%     %         sl=num2str(n);
%     %         start_str=['A' sl];
%     %         Sheet1.Range(start_str).Value=datafinal1(n);
%     %     end
%     for n=1:all
%         sl=num2str(n);
%         start_str=['A' sl];
%         Sheet1.Range(start_str).Value=datafinal2(n);
%     end
%     %     for n=1:all
%     %         sl=num2str(n);
%     %         start_str=['C' sl];
%     %         Sheet1.Range(start_str).Value=tt(n);
%     %     end
%     Workbook.Save; % 保存表格
%     Workbook.Close; % 关闭表格
%     Excel.Quit; % 退出Excel服务器
% end
end