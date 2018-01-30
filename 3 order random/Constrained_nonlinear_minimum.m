clear;
clc;
% matlabpool('open','local',4);
nvars=4;
A=[];
b=[];
Aeq=[];
beq=[];
LB=[200,20,2,2];
UB=[300,50,10,10];
x0=[200,25,5,5];
nonlcon=[];
options=optimset('UseParallel','always','Display','iter','MaxIter',400,'Algorithm','active-set');

[x,fval,exitflag] = fmincon(@main_blood,x0,A,b,Aeq,beq,LB,UB,nonlcon,options);

if exitflag==-2
    fprintf('No feasible point found.\n');
end
fprintf('best coef is Qb=%f Qf=%f V_t2=%f V_t4=%f \n',x);
fprintf('best value is %f\n',fval);

