% 参考色带
% Color_Map_Reference0

Width=100;
Height=500;
Limit=3/4;% 只选取HSV 中H值0~0.75 其余值为1
Color_Map_Reference0=ones([round(Height*Limit),Width,3]);
for i=1:round(Height*Limit)
    Color_Map_Reference0(i,:,1)=Limit-i/Height;
end
figure;
Color_Map_Reference0=hsv2rgb(Color_Map_Reference0);
imshow(Color_Map_Reference0,'InitialMagnification','fit');




