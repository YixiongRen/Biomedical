function out=Random(mean,std,n)
% ȡ������Ϊ6std ������̬�ֲ�99.73%
% out ��һ�б�ʾn���Ӷεľ�ֵ �ڶ��б�ʾ��Ӧ�ĸ���ֵ
out=zeros(2,n);
d=50;
temp=zeros(2,d+1);
for i=1:n
    temp(1,:)=-3*std+6*std/n*(i-1)+6*std/n/d*(0:d);
    temp(2,:)=1/sqrt(2*pi)/std*exp(-temp(1,:).^2/2/std^2);
    temp(2,1)=temp(2,1)/2;
    temp(2,end)=temp(2,end)/2;
%     display(temp);
    out(1,i)=temp(1,:)*temp(2,:).'/sum(temp(2,:))+mean;
    out(2,i)=sum(temp(2,:));
end
out(2,:)=out(2,:)/sum(out(2,:));
end
