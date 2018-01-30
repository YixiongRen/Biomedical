% 一阶Euler方法
function data=bagcycle(x1,x2,x3,x4,x5,delta_t,constant)%(v1,m_s1,m_s2,m_n1,m_n2)
N_element=size(x1,1);

V_s=constant.V_s;
alpha=constant.alpha;

N_T=50;%每一个单元分为6*N_big_element次混合；
T_t1=delta_t/N_T;
v2=constant.cv;
%V_iso=constant.V_iso;

data=struct('v1',x1,'m_s1',x2,'m_s2',x3,'m_n1',x4,'m_n2',x5);

for i=1:N_T
    x1=data.v1;
    x2=data.m_s1;
    x3=data.m_s2;
    x4=data.m_n1;
    x5=data.m_n2;
    
    temp=(v2-alpha.*x1)./(1+V_s.*x3);
    m_s2=sum(temp.*x3)/sum(temp);
    m_n2=sum(temp.*x5)/sum(temp);
    x3(:)=m_s2;
    x5(:)=m_n2;
    
    data=bagcycle_part2(x1,x2,x3,x4,x5,T_t1,constant);
end
end

function data=bagcycle_part2(x1,x2,x3,x4,x5,T_t1,constant)
N_element=size(x1,1);
data=struct('T',T_t1*ones(N_element,1),'v1',x1,'m_s1',x2,'m_s2',x3,'m_n1',x4,'m_n2',x5,'flag',5*ones(N_element,1));

A_c=constant.A_c;
L_p=constant.L_p;
P_s=constant.P_s;
%V_iso=constant.V_iso;
V_bc=constant.V_bc;
V_s=constant.V_s;
R=constant.R;
alpha=constant.alpha;
T=constant.T;
sigma=constant.sigma;

v2=constant.cv;

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

dy1=A_c.*L_p.*R.*T.*(m_n1.*(v1-V_bc)./(1+V_s.*m_s1).*(1+V_s.*m_s1)./(v1-V_bc)-m_n2.*(v2-alpha.*v1)./(1+V_s.*m_s2).*(1+V_s.*m_s2)./(v2-alpha.*v1)+sigma.*(m_s1-m_s2));
dy2=(1+V_s.*m_s1).^2./(v1-V_bc).*(((1-sigma).*((m_s2-m_s1)./(log(m_s2)-log(m_s1)+10^-20))-m_s1./(1+V_s.*m_s1)).*dy1+P_s.*A_c.*(m_s2-m_s1));
dy3=(1+V_s.*m_s2).^2./(v2-alpha.*v1).*(alpha.*(m_s2./(1+V_s.*m_s2)-m_s1./(1+V_s.*m_s1)).*dy1-alpha.*(v1-V_bc)./(1+V_s.*m_s1).^2.*dy2);
v1=v1+dy1*T_t1;
m_s1=m_s1+dy2*T_t1;
m_s2=m_s2+dy3*T_t1;

data.T=T_t1;
data.v1=v1;
data.m_s1=m_s1;
data.m_s2=m_s2;
data.m_n1=m_n1_0.*(v1_0-V_bc)./(1+V_s.*m_s1_0).*(1+V_s.*m_s1)./(v1-V_bc);
data.m_n2=m_n2_0.*(v2-alpha.*v1_0)./(1+V_s.*m_s2_0).*(1+V_s.*m_s2)./(v2-alpha.*v1);

end