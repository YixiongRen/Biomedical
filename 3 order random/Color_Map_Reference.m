% �ο�ɫ��
% Color_Map_Reference0

Width=100;
Height=500;
Limit=3/4;% ֻѡȡHSV ��Hֵ0~0.75 ����ֵΪ1
Color_Map_Reference0=ones([round(Height*Limit),Width,3]);
for i=1:round(Height*Limit)
    Color_Map_Reference0(i,:,1)=Limit-i/Height;
end
figure;
Color_Map_Reference0=hsv2rgb(Color_Map_Reference0);
imshow(Color_Map_Reference0,'InitialMagnification','fit');




