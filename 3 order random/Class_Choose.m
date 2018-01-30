function Constant=Class_Choose(Constant,flag)
% 将 Constant中的 A_c L_p P_s V_iso V_bc sigma alpha 按flag进行选取
Constant.A_c=Constant.A_c(flag);
Constant.L_p=Constant.L_p(flag);
Constant.P_s=Constant.P_s(flag);
Constant.V_iso=Constant.V_iso(flag);
Constant.V_bc=Constant.V_bc(flag);
Constant.sigma=Constant.sigma(flag);
Constant.alpha=Constant.alpha(flag);
end