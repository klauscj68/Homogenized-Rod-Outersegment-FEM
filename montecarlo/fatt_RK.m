% Calcolo del fattore RK/RK_tot nell'equazione A12 of NLP, J.Gen.Physiol. 116 (2000) 795--824
function [RK_rel]=fatt_RK(v, RK_tot, Rec_tot, M_bind, K_bind)

% v è la concentrazione del Ca2+

C_1=(v/K_bind(1))^2 * (1/K_bind(3)+M_bind/(K_bind(4)*K_bind(2))) * Rec_tot;
C_2=1 + (v/K_bind(1))^2 * (1+M_bind/K_bind(2));

A=C_1 * C_2;
B=C_1*(RK_tot/Rec_tot-1) + C_2;

Rec_rel=(sqrt(B^2+4*A) -B)/(2*A);
RK_rel=1/(1+C_1*Rec_rel);

return
