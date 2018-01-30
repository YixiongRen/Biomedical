% Compact Bag
% 采用类递归手法实现最小误差下的降采样
% 通过取子段线性误差最大点实现，并保留flag变化的点和不同参数的极值点
% 输入数据长度不超过1也不小于N 以防数据变化太剧烈
function out=Compact_Bag(Data,N,Length)
out=struct('T',zeros(1,N),'m_s2',zeros(1,N),'m_n2',zeros(1,N),'m_s1_mean',zeros(1,N),...
            'm_s1_std',zeros(1,N),'m_n1_mean',zeros(1,N),'m_n1_std',zeros(1,N));
Time=sum(Data.T(1,2:end)>0)+1;
Time=min([Time,Length]);
if Time>=N
    t=Data.T(1,1:Time);
    m_s2=Data.m_s2(1,1:Time);
    m_n2=Data.m_n2(1,1:Time);
    m_s1_mean=Data.m_s1_mean(1,1:Time);
    m_s1_std=Data.m_s1_std(1,1:Time);
    m_n1_mean=Data.m_n1_mean(1,1:Time);
    m_n1_std=Data.m_n1_std(1,1:Time);
    ps=Compact2(t,m_s2,m_n2,m_s1_mean,m_s1_std,m_n1_mean,m_n1_std,N);
    %display(ps);
    out.T(1,1:N)=t(ps);
    out.m_s2(1,1:N)=m_s2(ps);
    out.m_n2(1,1:N)=m_n2(ps);
    out.m_s1_mean(1,1:N)=m_s1_mean(ps);
    out.m_s1_std(1,1:N)=m_s1_std(ps);
    out.m_n1_mean(1,1:N)=m_n1_mean(ps);
    out.m_n1_std(1,1:N)=m_n1_std(ps);
else
    out.T(1,1:N)=Data.T(1,1:N);
    out.m_s2(1,1:N)=Data.m_s2(1,1:N);
    out.m_n2(1,1:N)=Data.m_n2(1,1:N);
    out.m_s1_mean(1,1:N)=Data.m_s1_mean(1,1:N);
    out.m_s1_std(1,1:N)=Data.m_s1_std(1,1:N);
    out.m_n1_mean(1,1:N)=Data.m_n1_mean(1,1:N);
    out.m_n1_std(1,1:N)=Data.m_n1_std(1,1:N);
end
end

function out=Compact2(t,m_s2,m_n2,m_s1_mean,m_s1_std,m_n1_mean,m_n1_std,N)
Length=size(t,2);
Flag=zeros(1,Length);
Flag(1)=1;
Flag(end)=1;

temp2=Flag+Peak(m_s2)+Peak(m_n2)+Peak(m_s1_mean)+Peak(m_s1_std)+Peak(m_n1_mean)+Peak(m_n1_std);
if sum(temp2>0)<=N
    Flag=temp2;
end

Flag=find(Flag>0);
Count=length(Flag);
while Count>N
    temp=randi([1,Count],1,1);
    Flag=Flag([1:temp-1,temp+1:Count]);
    Count=Count-1;
end

m_s2=Linear(m_s2/1000,Flag);
m_n2=Linear(m_n2/300,Flag);
m_s1_mean=Linear(m_s1_mean/1000,Flag);
m_s1_std=Linear(m_s1_std/100,Flag);
m_n1_mean=Linear(m_n1_mean/300,Flag);
m_n1_std=Linear(m_n1_std/30,Flag);
list=m_s2+m_n2+m_s1_mean+m_s1_std+m_n1_mean+m_n1_std;

Count=length(Flag);
while(Count<N)
%     list=Linear(list,Flag);
%     figure;
%     plot(list);hold on
    [~,ps]=max(list);
%     plot(ps,0,'r*');
    if ps==1
        break;
    end
    Flag=sort([Flag,ps]);
    Count=Count+1;
    ps=find(Flag==ps);
    temp1=abs(list(Flag(ps-1):Flag(ps))-((list(Flag(ps))-list(Flag(ps-1)))/(Flag(ps)-Flag(ps-1))*(0:Flag(ps)-Flag(ps-1))+list(Flag(ps-1))));
    temp2=abs(list(Flag(ps):Flag(ps+1))-((list(Flag(ps+1))-list(Flag(ps)))/(Flag(ps+1)-Flag(ps))*(0:Flag(ps+1)-Flag(ps))+list(Flag(ps))));
    list(Flag(ps-1):Flag(ps))=temp1;
    list(Flag(ps):Flag(ps+1))=temp2;
end
while Count<N
    temp=randi([1,N],1,1);
    Flag=union(temp,Flag);
    Count=length(Flag);
end
out=Flag;
end

function out=Peak(in)
temp=in(2:end)-in(1:end-1);
out=zeros(1,size(in,2));
for i=2:length(temp)-1
    if temp(i-1)*temp(i)<0
        out(i+1)=1;
    end
end
end


function out=Linear(in,flag)
N=length(flag);
out=zeros(1,length(in));
for i=1:N-1
    cut=flag(i+1)-flag(i);
    if cut>1
        out(flag(i):flag(i+1))=abs(in(flag(i):flag(i+1))-((in(flag(i+1))-in(flag(i)))/cut*(0:cut)+in(flag(i))));
    end
end
end