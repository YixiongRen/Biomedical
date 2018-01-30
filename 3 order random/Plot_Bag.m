function Plot_Bag( Bag )
% Use Bag
%Plot_Bag(Bag)
Time=sum(Bag.T(1,2:end)>0)+1;
subplot(2,2,1);
time=Bag.T(1,1:Time);
plot(time,Bag.m_s2(1,1:Time),'b-');
xlabel('t');ylabel('m_s2');
subplot(2,2,2);
plot(time,Bag.m_n2(1,1:Time),'b-');
xlabel('t');ylabel('m_n2');
subplot(2,2,3);
plot(time,Bag.m_s1_mean(1,1:Time),'b-');hold on
plot(time,Bag.m_s1_mean(1,1:Time)+Bag.m_s1_std(1,1:Time),'r-');hold on
plot(time,Bag.m_s1_mean(1,1:Time)-Bag.m_s1_std(1,1:Time),'r-');
xlabel('t');ylabel('m_s1');
subplot(2,2,4);
plot(time,Bag.m_n1_mean(1,1:Time),'b-');hold on
plot(time,Bag.m_n1_mean(1,1:Time)+Bag.m_n1_std(1,1:Time),'r-');hold on
plot(time,Bag.m_n1_mean(1,1:Time)-Bag.m_n1_std(1,1:Time),'r-');
xlabel('t');ylabel('m_n1');
hold on
end
