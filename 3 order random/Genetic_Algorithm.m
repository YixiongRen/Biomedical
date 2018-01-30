
clear;
clc;
matlabpool('open','local',4);
nvars=4;
A=[];
b=[];
Aeq=[];
beq=[];
LB=[200,20,2,2];
UB=[300,60,10,10];
nonlcon=[];
options=gaoptimset('Generations',30,'PopulationSize',28,'UseParallel','always','Display','iter');

[x,fval,exitflag] = ga(@main_blood,nvars,A,b,Aeq,beq,LB,UB,nonlcon,options);

if exitflag==-2
    fprintf('No feasible point found.\n');
end
fprintf('best coef is Qb=%f Qf=%f V_t2=%f V_t4=%f \n',x);
fprintf('best value is %f\n',fval);

