% function main_blood(var)
% N=size(var,2);
% for i=1:N
%     main_blood2(var(1,i),var(2,i),var(3,i),var(4,i));
% end
% end

% function out=main_blood(var)
% out=main_blood2(var(1),var(2),var(3),var(4));
% end
% 
% function out=main_blood2(Q_b,Q_f,V_t2,V_t4)
% 一阶Euler方法

clear;
clc;

Path='D:\data\2012-9-16-21：15\';
time=fix(clock);
SubPath=[num2str(time(1)),'-',num2str(time(2)),'-',num2str(time(3)),'-',num2str(time(4)),'-',num2str(time(5)),'-',num2str(time(6)),'=',num2str(randi(10000,1,1)),'\'];
mkdir([Path SubPath(1:end-1)]);

% temp=[];
% Q_b_str=num2str(Q_b);
% Q_f_str=num2str(Q_f);
% V_t2_str=num2str(V_t2);
% V_t4_str=num2str(V_t4);
% save([Path,SubPath,Q_b_str,'-',Q_f_str,'-',V_t2_str,'-',V_t4_str,'.mat'],'temp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   setting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_element=2000;
N_combine=10;% 每次选取10个体积元组成一个大体积元进入管道 保持外部溶液浓度相同
N_big_element=N_element/N_combine;
% A_c L_p P_s V_iso V_bc sigma alpha 可随机 但必须以列向量出现
Constant=struct('A_c',135*10^-12*ones(N_element,1),'L_p',1.74*10^-12*ones(N_element,1),'P_s',6.61*10^-8*ones(N_element,1),'V_iso',98.3*10^-18*ones(N_element,1),...
    'V_bc',0.283*98.3*10^-18*ones(N_element,1),'V_s',71*10^-6,'R',8.31,'alpha',10^7*ones(N_element,1),'T',298,'sigma',0.841*ones(N_element,1),...
    'cv',1*10^-8,'V_t1',5,'V_t3',85,'v1',[],'m_s1',1000,'m_s2',1000,'m_n1',290,'m_n2',290);
Constant.v1=Constant.V_iso+Constant.V_s.*Constant.m_s2.*(Constant.V_iso-Constant.V_bc);

Save_Data=1;% 存储 Data
Save_Bag=1;% 存储 Bag
Save_Random_Flag=1;% 存储 Random_Flag

%%%%% 管子底面积如下：s1,2,4=1.25*10^-5;s3=5*10^-4;s5=5*10^-3;所以长度为：l1,2,4=40cm;l3=17cm;体积为V1,2,4=5ml;V3=85ml;%%%%%%%%%%%%
V_t1=Constant.V_t1;
V_t2=5;
V_t3=Constant.V_t3;
V_t4=5;


Constant.V_t1=V_t1;
Constant.V_t2=V_t2;
Constant.V_t3=V_t3;
Constant.V_t4=V_t4;


Q_b=40;
Q_f=20;


N_Compact=100;
N_bag=N_element;% 小体积元个数
N_outcycle=0;% 小体积元个数
flag_bag=1:N_element;
flag_outcycle=zeros(1,N_element);
cycle_number=40*N_big_element;
store_size=40*N_Compact;
store_size_2=N_big_element+1;% 缓冲池大小
V_iso=Constant.V_iso;
store_flag=0;
Time=ones(1,N_element);
Time_Bag=1;
Gate=1;% 血袋进入管道的阀门开

%%%%%%%%%%%%%%%%%%%%%%%%%   initializaion   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 存储个体积元的数据
if Save_Data==1
    Data=struct('T',zeros(N_element,store_size),'v1',zeros(N_element,store_size),'m_s1',zeros(N_element,store_size),...
        'm_s2',zeros(N_element,store_size),'m_n1',zeros(N_element,store_size),'m_n2',zeros(N_element,store_size),'flag',zeros(N_element,store_size));
end
% 缓冲池
Data_2=struct('T',zeros(N_element,store_size_2),'v1',zeros(N_element,store_size_2),'m_s1',zeros(N_element,store_size_2),...
                'm_s2',zeros(N_element,store_size_2),'m_n1',zeros(N_element,store_size_2),'m_n2',zeros(N_element,store_size_2),'flag',zeros(N_element,store_size_2));
Data_2.T(:,1)=0;
Data_2.v1(:,1)=Constant.v1;
Data_2.m_s1(:,1)=Constant.m_s1;
Data_2.m_s2(:,1)=Constant.m_s2;
Data_2.m_n1(:,1)=Constant.m_n1;
Data_2.m_n2(:,1)=Constant.m_n2;
Data_2.flag(:,1)=5;

% 存储血袋内溶液参数的变化
if Save_Bag==1
    Bag=struct('T',zeros(1,store_size),'m_s2',zeros(1,store_size),'m_n2',zeros(1,store_size),'m_s1_mean',zeros(1,store_size),...
        'm_s1_std',zeros(1,store_size),'m_n1_mean',zeros(1,store_size),'m_n1_std',zeros(1,store_size));
end
% 缓冲池
Bag_2=struct('T',zeros(1,store_size_2),'m_s2',zeros(1,store_size_2),'m_n2',zeros(1,store_size_2),'m_s1_mean',zeros(1,store_size_2),...
            'm_s1_std',zeros(1,store_size_2),'m_n1_mean',zeros(1,store_size_2),'m_n1_std',zeros(1,store_size_2));
Bag_2.m_s2(1)=Constant.m_s2;
Bag_2.m_n2(1)=Constant.m_n2;
Bag_2.m_s1_mean(1)=Constant.m_s1;
Bag_2.m_s1_std(1)=0;
Bag_2.m_n1_mean(1)=Constant.m_n1;
Bag_2.m_n1_std(1)=0;

% 存储每次被随机选取进入管道的体积元的编号 不压缩
if Save_Random_Flag==1
    Random_Flag=zeros(N_combine,41*N_big_element);
    N_Random_Flag=0;
end

% 存储各体积元在管道中的具体位置 不保存
Position=zeros(3,N_element);% 管内分为4段，第一行表示管内段号1~4，第二行表示段内位置0~V_ti 第三行表示体积元的体积(用于第2、3段的数据处理)

Result.V_max=[];
Result.V_min=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%    calculation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parallel computing
matlabpool local 2

Current_Time=0;
for i0=1:100000%cycle_number
    if mod(i0,100)==0
        fprintf(['i0 = %d, time =  ',num2str(fix(clock)),'\n'],i0);
    end

    if N_bag>0 && i0>1 && Data_2.m_s1(flag_bag(1),Time(flag_bag(1)))<=50&&Data_2.m_s2(flag_bag(1),Time(flag_bag(1)))<=50
        Gate=0;
    end
    %%%%%%%%%%%%%% 初始化动态参数 %%%%%%%%%%%%%
%     Q_b=200;
%     Q_f=50;
    delta_t=Constant.cv*N_combine/(Q_b/10^6/60);
    Current_Time=Current_Time+delta_t;
    %%%%%%%%%%%%%%%%% store %%%%%%%%%%%%%%%%%
    % 将缓冲池Data_2中的数据压缩后放入Data 对Bag_2做相同处理
    if min(Time)==store_size_2 || (Gate==0&&N_outcycle==0)
        Result.V_max=max([Result.V_max,Data_2.v1],[],2);
        Result.V_min=min([Result.V_min,Data_2.v1],[],2);
        
        store_flag=store_flag+1;
        Time(:)=1;
        if Save_Data==1
            temp=Compact(Data_2,N_Compact,store_size_2-1);
            temp2=(store_flag-1)*N_Compact+1:store_flag*N_Compact;
            Data.T(:,temp2)=temp.T;
            Data.v1(:,temp2)=temp.v1;
            Data.m_s1(:,temp2)=temp.m_s1;
            Data.m_s2(:,temp2)=temp.m_s2;
            Data.m_n1(:,temp2)=temp.m_n1;
            Data.m_n2(:,temp2)=temp.m_n2;
            Data.flag(:,temp2)=temp.flag;
        end
        Data_2.T(:,1)=Data_2.T(:,end);Data_2.T(:,2:end)=0;
        Data_2.v1(:,1)=Data_2.v1(:,end);Data_2.v1(:,2:end)=0;
        Data_2.m_s1(:,1)=Data_2.m_s1(:,end);Data_2.m_s1(:,2:end)=0;
        Data_2.m_s2(:,1)=Data_2.m_s2(:,end);Data_2.m_s2(:,2:end)=0;
        Data_2.m_n1(:,1)=Data_2.m_n1(:,end);Data_2.m_n1(:,2:end)=0;
        Data_2.m_n2(:,1)=Data_2.m_n2(:,end);Data_2.m_n2(:,2:end)=0;
        Data_2.flag(:,1)=Data_2.flag(:,end);Data_2.flag(:,2:end)=0;
            
        Time_Bag=1;
        if Save_Bag==1
            temp2=(store_flag-1)*N_Compact+1:store_flag*N_Compact;
            temp=Compact_Bag(Bag_2,N_Compact,store_size_2-1);
            Bag.T(1,temp2)=temp.T;
            Bag.m_s2(1,temp2)=temp.m_s2;
            Bag.m_n2(1,temp2)=temp.m_n2;
            Bag.m_s1_mean(1,temp2)=temp.m_s1_mean;
            Bag.m_s1_std(1,temp2)=temp.m_s1_std;
            Bag.m_n1_mean(1,temp2)=temp.m_n1_mean;
            Bag.m_n1_std(1,temp2)=temp.m_n1_std;
        end
        Bag_2.T(1)=Bag_2.T(end);Bag_2.T(2:end)=0;
        Bag_2.m_s2(1)=Bag_2.m_s2(end);Bag_2.m_s2(2:end)=0;
        Bag_2.m_n2(1)=Bag_2.m_n2(end);Bag_2.m_n2(2:end)=0;
        Bag_2.m_s1_mean(1)=Bag_2.m_s1_mean(end);Bag_2.m_s1_mean(2:end)=0;
        Bag_2.m_s1_std(1)=Bag_2.m_s1_std(end);Bag_2.m_s1_std(2:end)=0;
        Bag_2.m_n1_mean(1)=Bag_2.m_n1_mean(end);Bag_2.m_n1_mean(2:end)=0;
        Bag_2.m_n1_std(1)=Bag_2.m_n1_std(end);Bag_2.m_n1_std(2:end)=0;
    end
    %%%%%%%%%%%  Break  %%%%%%%%%%%
    if Gate==0 && N_outcycle==0
        break
    end
    %%%%%%%%%%%%    bag   %%%%%%%%%%%%
    if N_bag>0
        temp=flag_bag(1:N_bag);
        temp2=Time(1);% Time各元相同
        x1=Data_2.v1(temp,temp2);
        x2=Data_2.m_s1(temp,temp2);
        x3=Data_2.m_s2(temp,temp2);
        x4=Data_2.m_n1(temp,temp2);
        x5=Data_2.m_n2(temp,temp2);
        data=bagcycle(x1,x2,x3,x4,x5,delta_t,Class_Choose(Constant,temp));
        Time(temp)=Time(temp)+1;
        temp2=temp2+1;
        Data_2.T(temp,temp2)=delta_t+Data_2.T(temp,temp2-1);
        Data_2.v1(temp,temp2)=data.v1;
        Data_2.m_s1(temp,temp2)=data.m_s1;
        Data_2.m_s2(temp,temp2)=data.m_s2;
        Data_2.m_n1(temp,temp2)=data.m_n1;
        Data_2.m_n2(temp,temp2)=data.m_n2;
        Data_2.flag(temp,temp2)=5;
        
        Time_Bag=Time_Bag+1;
        Bag_2.T(Time_Bag)=Bag_2.T(Time_Bag-1)+delta_t;
        Bag_2.m_s2(Time_Bag)=data.m_s2(1,1);
        Bag_2.m_n2(Time_Bag)=data.m_n2(1,1);
        Bag_2.m_s1_mean(Time_Bag)=mean(data.m_s1);
        Bag_2.m_s1_std(Time_Bag)=std(data.m_s1);
        Bag_2.m_n1_mean(Time_Bag)=mean(data.m_n1);
        Bag_2.m_n1_std(Time_Bag)=std(data.m_n1);
    end
    %%%%%%%%%%%%    outcycle   %%%%%%%%%%%%
    if N_outcycle>0
        temp=flag_outcycle(1:N_outcycle);
        temp2=min(Time);% Time各元不相同
        x1=Data_2.v1(temp,temp2);
        x2=Data_2.m_s1(temp,temp2);
        x3=Data_2.m_s2(temp,temp2);
        x4=Data_2.m_n1(temp,temp2);
        x5=Data_2.m_n2(temp,temp2);
        [data,position]=outcycle(x1,x2,x3,x4,x5,delta_t,Q_b,Q_f,V_t2,V_t4,Class_Choose(Constant,temp),Position(:,temp),N_combine);
        Position(:,temp)=position;
        Time(temp)=Time(temp)+1;
        temp2=temp2+1;
        Data_2.T(temp,temp2)=delta_t+Data_2.T(temp,temp2-1);
        Data_2.v1(temp,temp2)=data.v1;
        Data_2.m_s1(temp,temp2)=data.m_s1;
        Data_2.m_s2(temp,temp2)=data.m_s2;
        Data_2.m_n1(temp,temp2)=data.m_n1;
        Data_2.m_n2(temp,temp2)=data.m_n2;
        Data_2.flag(temp,temp2)=data.flag;
    end
    %%%%%%%%%%%%    in-out   %%%%%%%%%%%%
    flag_input=-ones(1,N_combine);% 进入管道
    flag_output=-ones(1,N_combine);% 进入血袋 可能长为2*N_combine
    
    if N_outcycle>0 && (Position(1,flag_outcycle(N_outcycle))==4&&Position(2,flag_outcycle(N_outcycle))+Constant.cv/3>=V_t4/10^6)
        flag_output=flag_outcycle(N_outcycle-N_combine+1:N_outcycle);
        Position(:,flag_output)=[zeros(1,N_combine);zeros(1,N_combine);Constant.cv*ones(1,N_combine)];
        flag_outcycle(N_outcycle-N_combine+1:N_outcycle)=0;
        N_outcycle=N_outcycle-N_combine;
        % 可能出现2个大体积元进入bag
        if N_outcycle>0 && (Position(1,flag_outcycle(N_outcycle))==4&&Position(2,flag_outcycle(N_outcycle))+Constant.cv/3>=V_t4/10^6)
            flag_output=flag_outcycle(N_outcycle-N_combine+1:N_outcycle);
            Position(:,flag_output)=[zeros(1,N_combine);zeros(1,N_combine);Constant.cv*ones(1,N_combine)];
            flag_outcycle(N_outcycle-N_combine+1:N_outcycle)=0;
            N_outcycle=N_outcycle-N_combine;
        end
    end
    if Gate==1 && N_bag>0
        % 顺序进入管道
        choose=1:N_combine;
        flag_input=flag_bag(1:N_combine);
        Position(:,flag_input)=[ones(1,N_combine);zeros(1,N_combine);Constant.cv*ones(1,N_combine)];
        flag_bag=flag_bag(N_combine+1:end);
        N_bag=N_bag-N_combine;
        
%         % 随机进入管道
%         choose=Random_Choose(N_bag,N_combine);
%         flag_input=flag_bag(choose);
%         Position(:,flag_input)=[ones(1,N_combine);zeros(1,N_combine);Constant.cv*ones(1,N_combine)];
%         flag_bag(choose)=0;
%         N_bag=N_bag-N_combine;
%         temp=sort(flag_bag);
%         flag_bag(:)=0;
%         flag_bag(1:N_bag)=temp(end-N_bag+1:end);
    end
    
    if flag_output(1)>-1
        flag_bag(N_bag+1:N_bag+length(flag_output))=flag_output;
        N_bag=N_bag+length(flag_output);
    end
    
    if flag_input(1)>-1
        if Save_Random_Flag==1
            N_Random_Flag=N_Random_Flag+1;
            Random_Flag(:,N_Random_Flag)=flag_input.';
        end
        N_outcycle=N_outcycle+N_combine;
        if N_outcycle>1
            flag_outcycle(N_combine+(1:N_outcycle))=flag_outcycle(1:N_outcycle);
        end
        flag_outcycle(1:N_combine)=flag_input;
    end
end
return
if i0==cycle_number
    fprintf('out of time\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Store Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Save_Random_Flag==1
    Random_Flag=Random_Flag(1:N_Random_Flag);
    save([Path,SubPath,'Random_Flag.mat'],'Random_Flag','-v7.3');
end
if Save_Data==1
    Time=sum(Data.T(1,2:end)>0)+1;
    Data.T=Data.T(:,1:Time);
    Data.v1=Data.v1(:,1:Time);
    Data.m_s1=Data.m_s1(:,1:Time);
    Data.m_s2=Data.m_s2(:,1:Time);
    Data.m_n1=Data.m_n1(:,1:Time);
    Data.m_n2=Data.m_n2(:,1:Time);
    Data.flag=Data.flag(:,1:Time);
    save([Path,SubPath,'Data.mat'],'Data','-v7.3');
end
if Save_Bag==1
    Bag.m_s2=Bag.m_s2(1:Time);
    Bag.m_n2=Bag.m_n2(1:Time);
    Bag.m_s1_mean=Bag.m_s1_mean(1:Time);
    Bag.m_s1_std=Bag.m_s1_std(1:Time);
    Bag.m_n1_mean=Bag.m_n1_mean(1:Time);
    Bag.m_n1_std=Bag.m_n1_std(1:Time);
    save([Path,SubPath,'Bag.mat'],'Bag','-v7.3');
end

Result.v1_end=Data_2.v1(:,1);
Result.m_s1_end=Data_2.m_s1(:,1);
Result.m_s2_end=Data_2.m_s2(:,1);
Result.m_n1_end=Data_2.m_n1(:,1);
Result.m_n2_end=Data_2.m_n2(:,1);
Result.T_end=Current_Time;

save([Path,SubPath,'Result.mat'],'Result','-v7.3');
save([Path,SubPath,'Constant.mat'],'Constant','-v7.3');

matlabpool close


%***********
% V_max_tolerance=1.15*Constant.V_iso;
% V_min_tolerance=0.95*Constant.V_iso;
% T_end_tolerance=20*60;
% out=exp(50*(mean(Result.V_max)-V_max_tolerance)/Constant.V_iso)+...
%     exp(50*(V_min_tolerance-Constant.V_iso)/Constant.V_iso)+...
%     exp(50*(Result.T_end(1)-T_end_tolerance)/T_end_tolerance);
% save([Path,SubPath,'out.mat'],'out','-v7.3');
% end