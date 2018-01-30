function Plot_Element( Data,element )
% Use Data
%Plot_Element(Data,2)
Time=sum(Data.T(element,2:end)>0)+1;
subplot(3,2,1);
plot(Data.T(element,1:Time),Data.v1(element,1:Time));hold on
xlabel('t');ylabel('v1');
subplot(3,2,2);
plot(Data.T(element,1:Time),Data.m_s1(element,1:Time));hold on
xlabel('t');ylabel('m_s1');
subplot(3,2,3);
plot(Data.T(element,1:Time),Data.m_s2(element,1:Time));hold on
xlabel('t');ylabel('m_s2');
subplot(3,2,4);
plot(Data.T(element,1:Time),Data.m_n1(element,1:Time));hold on
xlabel('t');ylabel('m_n1');
subplot(3,2,5);
plot(Data.T(element,1:Time),Data.m_n2(element,1:Time));hold on
xlabel('t');ylabel('m_n2');
subplot(3,2,6);
plot(Data.T(element,1:Time),Data.flag(element,1:Time));hold on
xlabel('t');ylabel('flag');
hold on
end

