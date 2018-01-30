% 一阶Euler方法
function [data,position]=outcycle(x1,x2,x3,x4,x5,delta_t,Q_b,Q_f,V_t2,V_t4,constant,position,N_combine)
N_element=size(x1,1);
V_t1=constant.V_t1;
V_t3=constant.V_t3;

data=struct('v1',zeros(N_element,1),'m_s1',zeros(N_element,1),'m_s2',zeros(N_element,1),...
    'm_n1',zeros(N_element,1),'m_n2',zeros(N_element,1),'flag',zeros(N_element,1));
%%%%%%%%%%%%%%% part 1 %%%%%%%%%%%%%%%
ps=find(position(1,:)==1);
if ~isempty(ps)
    % 选出进入 part 2
    ps2=ps(end-N_combine+1:end);
    if position(2,ps2(1))+constant.cv/3>=V_t1/10^6
        n_s20=(constant.cv-constant.alpha(ps2).*x1(ps2))./(1+constant.V_s*x3(ps2)).*x3(ps2);
        a=((1+Q_f/Q_b)*constant.cv-constant.alpha(ps2).*x1(ps2)-constant.V_s*n_s20);
        x3(ps2)=n_s20./a;
        position(:,ps2)=[ones(1,N_combine)*2;position(2,ps2)-V_t1/10^6;ones(1,N_combine)*constant.cv*(Q_b+Q_f)/Q_b];
        ps=ps(1:end-N_combine);
    end
end
if ~isempty(ps)
    data2=outcycle_2(x1(ps),x2(ps),x3(ps),x4(ps),x5(ps),delta_t,constant,50,N_combine,ps,position(3,ps));
    position(2,ps)=position(2,ps)+constant.cv*N_combine;
    data.v1(ps)=data2.v1;
    data.m_s1(ps)=data2.m_s1;
    data.m_s2(ps)=data2.m_s2;
    data.m_n1(ps)=data2.m_n1;
    data.m_n2(ps)=data2.m_n2;
    data.flag(ps)=1;
end
%%%%%%%%%%%%%%% part 2 %%%%%%%%%%%%%%%
ps=find(position(1,:)==2);
while ~isempty(ps) && position(2,ps(end))+constant.cv/3>=V_t2/10^6
    % 选出进入 part 3 可能有2个进入 part 3
    ps2=ps(end-N_combine+1:end);
    position(1:2,ps2)=[3*ones(1,N_combine);position(2,ps2)-V_t2/10^6];
    ps=ps(1:end-N_combine);
end
if ~isempty(ps)
    data2=outcycle_2(x1(ps),x2(ps),x3(ps),x4(ps),x5(ps),delta_t,constant,100,N_combine,ps,position(3,ps));
    position(2,ps)=position(2,ps)+position(3,ps)*N_combine;
    data.v1(ps)=data2.v1;
    data.m_s1(ps)=data2.m_s1;
    data.m_s2(ps)=data2.m_s2;
    data.m_n1(ps)=data2.m_n1;
    data.m_n2(ps)=data2.m_n2;
    data.flag(ps)=2;
end
%%%%%%%%%%%%%%% part 3 %%%%%%%%%%%%%%%
ps=find(position(1,:)==3);
while ~isempty(ps) && position(2,ps(end))+constant.cv/3>=V_t3/10^6
    % 选出进入 part 4 可能有2个进入 part 4
    ps2=ps(end-N_combine+1:end);
    position(:,ps2)=[4*ones(1,N_combine);position(2,ps2)-V_t3/10^6;constant.cv*ones(1,N_combine)];
    ps=ps(1:end-N_combine);
end
if ~isempty(ps)
% 保证每个体积元在出超滤管时体积尽可能接近constant.cv
    N=length(ps);
    cut=50;
    v20=zeros(cut,N);
    delta_t_cut=delta_t/cut;
    Q=constant.cv*(Q_b+Q_f)/Q_b/cut*N_combine;% 进口流量
    for i0=1:cut
        v2=position(3,ps);% 体积元大小
        v_residue=v2-constant.cv;% 残余的多余溶液体积 向量
        v_h2=V_t3/10^6-position(2,ps);% 残余的超滤管体积 向量
        v_move=v_residue(1)*(1-(1-Q/v_h2(1)));% 超滤量
        ps2=ps(1:N_combine);
        position(2,ps2)=position(2,ps2)+Q;
        position(3,ps2)=position(3,ps2)-v_move;
        for i=N_combine+1:N_combine:N-N_combine+1
            ps2=ps(i:i+N_combine-1);
            ps3=ps2(1);
            v_temp=position(2,ps3-1)+position(3,ps3-1)*N_combine-position(2,ps3);% 位移量
            position(2,ps2)=position(2,ps2)+v_temp;
            v_move=v_residue(i)*(1-(1-v_temp/v_h2(i))^2);
            position(3,ps2)=position(3,ps2)-v_move;
        end
        v20(i0,:)=position(3,ps);
    end
    data2=outcycle_3(x1(ps),x2(ps),x3(ps),x4(ps),x5(ps),delta_t_cut,constant,cut,N_combine,ps,v20);
    data.v1(ps)=data2.v1;
    data.m_s1(ps)=data2.m_s1;
    data.m_s2(ps)=data2.m_s2;
    data.m_n1(ps)=data2.m_n1;
    data.m_n2(ps)=data2.m_n2;
    data.flag(ps)=3;
end
%%%%%%%%%%%%%%% part 4 %%%%%%%%%%%%%%%
ps=find(position(1,:)==4);
if ~isempty(ps)
    data2=outcycle_2(x1(ps),x2(ps),x3(ps),x4(ps),x5(ps),delta_t,constant,50,N_combine,ps,position(3,ps));
    position(2,ps)=position(2,ps)+constant.cv*N_combine;
    data.v1(ps)=data2.v1;
    data.m_s1(ps)=data2.m_s1;
    data.m_s2(ps)=data2.m_s2;
    data.m_n1(ps)=data2.m_n1;
    data.m_n2(ps)=data2.m_n2;
    data.flag(ps)=4;
end

end

function data=outcycle_2(x1,x2,x3,x4,x5,delta_t,constant,cut,N_combine,ps,v2)
N_element=size(x1,1);
N_big_element=N_element/N_combine;
T_t1=delta_t/cut;
data=struct('v1',zeros(N_element,1),'m_s1',zeros(N_element,1),'m_s2',zeros(N_element,1),...
    'm_n1',zeros(N_element,1),'m_n2',zeros(N_element,1));

A_c=constant.A_c(ps);
L_p=constant.L_p(ps);
P_s=constant.P_s(ps);

V_bc=constant.V_bc(ps);
V_s=constant.V_s;
R=constant.R;
alpha=constant.alpha(ps);
T=constant.T;
sigma=constant.sigma(ps);

v2=v2.';

v1=x1;
m_s1=x2;
m_s2=x3;
m_n1=x4;
m_n2=x5;
v1_0=x1;
m_s1_0=x2;
m_s2_0=x3;
m_n1_0=x4;
m_n2_0=x5;

for i=1:cut
    for j=1:N_big_element
        ps2=N_combine*(j-1)+(1:N_combine);
        temp=(v2(ps2)-alpha(ps2).*v1(ps2))./(1+V_s.*m_s2(ps2));
        x3=sum(temp.*m_s2(ps2))/sum(temp);
        x5=sum(temp.*m_n2(ps2))/sum(temp);
        m_s2(ps2)=x3;
        m_n2(ps2)=x5;
    end
    
    dy1=A_c.*L_p.*R.*T.*(m_n1.*(v1-V_bc)./(1+V_s.*m_s1).*(1+V_s.*m_s1)./(v1-V_bc)-m_n2.*(v2-alpha.*v1)./(1+V_s.*m_s2).*(1+V_s.*m_s2)./(v2-alpha.*v1)+sigma.*(m_s1-m_s2));
    dy2=(1+V_s.*m_s1).^2./(v1-V_bc).*(((1-sigma).*((m_s2-m_s1)./(log(m_s2)-log(m_s1)+10^-20))-m_s1./(1+V_s.*m_s1)).*dy1+P_s.*A_c.*(m_s2-m_s1));
    dy3=(1+V_s.*m_s2).^2./(v2-alpha.*v1).*(alpha.*(m_s2./(1+V_s.*m_s2)-m_s1./(1+V_s.*m_s1)).*dy1-alpha.*(v1-V_bc)./(1+V_s.*m_s1).^2.*dy2);
    v1=v1+dy1*T_t1;
    m_s1=m_s1+dy2*T_t1;
    m_s2=m_s2+dy3*T_t1;

    data.v1=v1;
    data.m_s1=m_s1;
    data.m_s2=m_s2;
    data.m_n1=m_n1_0.*(v1_0-V_bc)./(1+V_s.*m_s1_0).*(1+V_s.*m_s1)./(v1-V_bc);
    data.m_n2=m_n2_0.*(v2-alpha.*v1_0)./(1+V_s.*m_s2_0).*(1+V_s.*m_s2)./(v2-alpha.*v1);
    
    m_n1=data.m_n1;
    m_n2=data.m_n2;
end
end

function data=outcycle_3(x1,x2,x3,x4,x5,delta_t,constant,cut,N_combine,ps,v20)% 超滤管内
N_element=size(x1,1);
N_big_element=N_element/N_combine;
T_t1=delta_t/cut;
data=struct('v1',zeros(N_element,1),'m_s1',zeros(N_element,1),'m_s2',zeros(N_element,1),...
    'm_n1',zeros(N_element,1),'m_n2',zeros(N_element,1));

A_c=constant.A_c(ps);
L_p=constant.L_p(ps);
P_s=constant.P_s(ps);

V_bc=constant.V_bc(ps);
V_s=constant.V_s;
R=constant.R;
alpha=constant.alpha(ps);
T=constant.T;
sigma=constant.sigma(ps);

v1=x1;
m_s1=x2;
m_s2=x3;
m_n1=x4;
m_n2=x5;
v1_0=x1;
m_s1_0=x2;
m_s2_0=x3;
m_n1_0=x4;
m_n2_0=x5;

for i=1:cut
    v2=v20(i,:).';% 数据均采用列向量计算
    
    for j=1:N_big_element
        ps2=N_combine*(j-1)+(1:N_combine);
        temp=(v2(ps2)-alpha(ps2).*v1(ps2))./(1+V_s.*m_s2(ps2));
        x3=sum(temp.*m_s2(ps2))/sum(temp);
        x5=sum(temp.*m_n2(ps2))/sum(temp);
        m_s2(ps2)=x3;
        m_n2(ps2)=x5;
    end
    
    dy1=A_c.*L_p.*R.*T.*(m_n1.*(v1-V_bc)./(1+V_s.*m_s1).*(1+V_s.*m_s1)./(v1-V_bc)-m_n2.*(v2-alpha.*v1)./(1+V_s.*m_s2).*(1+V_s.*m_s2)./(v2-alpha.*v1)+sigma.*(m_s1-m_s2));
    dy2=(1+V_s.*m_s1).^2./(v1-V_bc).*(((1-sigma).*((m_s2-m_s1)./(log(m_s2)-log(m_s1)+10^-20))-m_s1./(1+V_s.*m_s1)).*dy1+P_s.*A_c.*(m_s2-m_s1));
    dy3=(1+V_s.*m_s2).^2./(v2-alpha.*v1).*(alpha.*(m_s2./(1+V_s.*m_s2)-m_s1./(1+V_s.*m_s1)).*dy1-alpha.*(v1-V_bc)./(1+V_s.*m_s1).^2.*dy2);
    v1=v1+dy1*T_t1;
    m_s1=m_s1+dy2*T_t1;
    m_s2=m_s2+dy3*T_t1;
    
    data.v1=v1;
    data.m_s1=m_s1;
    data.m_s2=m_s2;
    data.m_n1=m_n1_0.*(v1_0-V_bc)./(1+V_s.*m_s1_0).*(1+V_s.*m_s1)./(v1-V_bc);
    data.m_n2=m_n2_0.*(v2-alpha.*v1_0)./(1+V_s.*m_s2_0).*(1+V_s.*m_s2)./(v2-alpha.*v1);
    
    m_n1=data.m_n1;
    m_n2=data.m_n2;
end
end
