% use Data
N_point=sum(Data.T(1,2:end)>0)+1;
N_element=size(Data.T,1);
Cycle=10;
D1=2;
N1=N_element/D1;
Resize=[500,Cycle*250];
Field='m_n2';

Limit=3/4;% ֻѡȡHSV ��Hֵ0~0.75 ����ֵΪ1
Image=zeros(N1,Cycle*N_element);
delta_t=0.03;
V_max=max(max(getfield(Data,Field)));
V_min=min(min(getfield(Data,Field)));
N_Compact=100; %ÿ������ѹ����100���� ���1�����ڲ���



for i=D1:D1:N_element
    ps=1;
    data1=Data.T(i,:);
    data2=getfield(Data,Field,{i,1:size(Data.T,2)});
    j=i/D1;
    
    while  ps<Cycle*N_Compact
        ps2=ps+1;
        t1=data1(ps);
        t2=data1(ps2);
        v1=data2(ps);
        v2=data2(ps2);
        Image(j,1+round(t1/delta_t):1+round(t2/delta_t))=v1+(v2-v1)/(t2-t1)*delta_t*(-round(t1/delta_t)+1+round(t2/delta_t));
        ps=ps+1;
        
    end

end
Image=Limit-(Image-V_min)/(V_max-V_min)*Limit;
Image=imresize(Image, Resize);
% imshow(Image,'InitialMagnification','fit');
figure;
Image2=ones([Resize,3]);
Image2(:,:,1)=Image;
Image2=hsv2rgb(Image2);
imshow(Image2,'InitialMagnification','fit');
















