clear all; clc
    
mpc=pglib_opf_case200_activ; type=mpc.bus(:,2);
P=-mpc.bus(:,3); P(mpc.gen(:,1))=P(mpc.gen(:,1))+mpc.gen(:,2);
Q=mpc.bus(:,4)*0; Q(mpc.bus(:,2)==1,1)=-mpc.bus(mpc.bus(:,2)==1,4);
for i=1:1:size(mpc.gen,1)
if mpc.bus(mpc.gen(i,1),2)==1
    Q(mpc.gen(i,1))=Q(mpc.gen(i,1))+mpc.gen(i,3);
    end
end
Ym=makeYbus(mpc); Ym=full(Ym); V=mpc.bus(:,8); Vdegs=deg2rad(mpc.bus(:,9));

tol = 1e-5; max_iter = 10;
[V1,Vdegs1, P1, Q1] = NTRaph(V, Vdegs, type, P/mpc.baseMVA, Q/mpc.baseMVA, full(Ym), tol, max_iter);

[V2,Vdegs2, slackPQ] = vectorized_newton(type, V, full(Ym), P/mpc.baseMVA, Q/mpc.baseMVA, max_iter, tol);