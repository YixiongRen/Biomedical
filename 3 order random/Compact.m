% Compact Data
% 采用类递归手法实现最小误差下的降采样
% 通过取子段线性误差最大点实现，并保留flag变化的点和不同参数的极值点
% 输入数据长度不超过N_element也不小于N 以防数据变化太剧烈
function out=Compact(Data,N,Length)
N_element=size(Data.T,1);
out=struct('T',zeros(N_element,N),'v1',zeros(N_element,N),'m_s1',zeros(N_element,N),'m_s2',zeros(N_element,N),'m_n1',zeros(N_element,N),...
            'm_n2',zeros(N_element,N),'flag',zeros(N_element,N));
parfor i=1:N_element
    Time=sum(Data.T(i,2:end)>0)+1;
    Time=min([Time,Length]);
    if Time>=N
        t=Data.T(i,1:Time);
        v1=Data.v1(i,1:Time);
        m_s1=Data.m_s1(i,1:Time);
        m_s2=Data.m_s2(i,1:Time);
        m_n1=Data.m_n1(i,1:Time);
        m_n2=Data.m_n2(i,1:Time);
        flag=Data.flag(i,1:Time);
        ps=Compact2(t,v1,m_s1,m_s2,m_n1,m_n2,flag,N);
        %display(ps);
        out.T(i,1:N)=t(ps);
        out.v1(i,1:N)=v1(ps);
        out.m_s1(i,1:N)=m_s1(ps);
        out.m_s2(i,1:N)=m_s2(ps);
        out.m_n1(i,1:N)=m_n1(ps);
        out.m_n2(i,1:N)=m_n2(ps);
        out.flag(i,1:N)=flag(ps);
    else
        out.T(i,1:N)=Data.T(i,1:N);
        out.v1(i,1:N)=Data.v1(i,1:N);
        out.m_s1(i,1:N)=Data.m_s1(i,1:N);
        out.m_s2(i,1:N)=Data.m_s2(i,1:N);
        out.m_n1(i,1:N)=Data.m_n1(i,1:N);
        out.m_n2(i,1:N)=Data.m_n2(i,1:N);
        out.flag(i,1:N)=Data.flag(i,1:N);
    end
end
end

function out=Compact2(t,v1,m_s1,m_s2,m_n1,m_n2,flag,N)
Length=size(t,2);
Flag=zeros(1,Length);
Flag(1)=1;
Flag(end)=1;

temp=abs(flag(2:end)-flag(1:end-1));
Flag=Flag+[temp,0]+[0,temp];
temp2=Flag+Peak(v1)+Peak(m_s1)+Peak(m_s2)+Peak(m_n1)+Peak(m_n2);
if sum(Flag>0)<N&&sum(temp2>0)<=N
    Flag=temp2;
end

Flag=find(Flag>0);
Count=length(Flag);
while Count>N
    temp=randi([1,Count],1,1);
    Flag=Flag([1:temp-1,temp+1:Count]);
    Count=Count-1;
end


v1=Linear(v1/10^-15,Flag);
m_s1=Linear(m_s1/1000,Flag);
m_s2=Linear(m_s2/1000,Flag);
m_n1=Linear(m_n1/300,Flag);
m_n2=Linear(m_n2/300,Flag);
list=v1+m_s1+m_s2+m_n1+m_n2;

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