clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%修改Atfd格式使其能够被importdata读取
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path=['D:\Data and Files\3132\'];
namelist=dir([path,'*.asc']);
l=length(namelist);
for i=1:l
    fullname=[path,namelist(i).name];
    fid=fopen(fullname);
    j = 0;
    while ~feof(fid)
        tline = fgetl(fid);                            %逐行读取原始文件
        j = j+1;
        newline{j} = tline;                               %创建新对象接受原始文件每行数据                           %想增加的数据                         %想增加的数据
        if j == 51                                     %判断是否到达待修改的行
            newline{j} =strrep(tline,'.,',' ');	                  %新增数据	  							%将待修改的行全部替换为新内容							%如需要一行替换为2行
        end
    end
    fclose(fid);
    filename_new=fullname;
    fileID = fopen(filename_new,'w+');                    %以可读写的方式打开输出文件，文件若存在则清空文件内容从文件头部开始写，若不存在则根据文件名创建新文件并只写打开
    for k=1:j
        fprintf(fileID,'%s\t\n',newline{k});              %将newline内的内容逐行写出
    end
    fclose(fileID);
end
