function [rhs_psi, rhs_Pi_psi, rhs_Phi_psi] = compute_RHS
Globals1D;
GlobalsGR;

%constraints
%C0 = H0+gamma0
%C1 = H1+gamma1

%source terms
src_psi = -lapse.*Pi_psi - paragamma1.*shift.*Phi_psi
src_Pi_psi = lapse.*gamma11.*gamma1.*Phi_psi + lapse.*Pi_psi.*(invg00.*gamma0.*(-lapse)+invg01.*gamma1.*(-lapse)) ...
        - paragamma1.*paragamma2.*shift.*Phi_psi - lapse.*(gamma11.*Phi_psi.*(normal0.*Pi01+normal1.*Pi11) ...
        + 0.5.*Pi_psi.*(normal0.*normal0.*Pi00 + 2.*normal0.*normal1.*Pi01 + normal1.*normal1.*Pi11))
src_Phi_psi = lapse.*(0.5.*Pi_psi.*(normal0.*normal0.*Phi00 + 2.*normal0.*normal1.*Phi01 + normal1.*normal1.*Phi11) ...
       + gamma11.*Phi_psi.*(normal0.*Phi01 + normal1.*Phi11) - paragamma2.*Phi_psi)

%Hhat terms
[deri_plus_psi, deri_minus_psi] = deri(psi, psi, psi);
[deri_plus_Pi_psi, deri_minus_Pi_psi] = deri(Pi_psi, Pi_psi, Pi_psi);
[deri_plus_Phi_psi, deri_minus_Phi_psi] = deri(Phi_psi, Phi_psi, Phi_psi);

max_alpha = 1.0;
avg_deri_psi = 0.5.*(deri_plus_psi+deri_minus_psi);
avg_deri_Pi_psi = 0.5.*(deri_plus_Pi_psi+deri_minus_Pi_psi);
avg_deri_Phi_psi = 0.5.*(deri_plus_Phi_psi+deri_minus_Phi_psi);

Hhat_psi = -(1+paragamma1).*shift.*avg_deri_psi - max_alpha.*0.5.*(deri_plus_psi-deri_minus_psi)
Hhat_Pi_psi = -paragamma1.*paragamma2.*shift.*avg_deri_psi - shift.*avg_deri_Pi_psi ...
        + lapse.*gamma11.*avg_deri_Phi_psi - max_alpha.*0.5.*(deri_plus_Pi_psi-deri_minus_Pi_psi)
Hhat_Phi_psi = -paragamma2.*lapse.*avg_deri_psi + lapse.*avg_deri_Pi_psi - shift.*avg_deri_Phi_psi ...
        - max_alpha.*0.5.*(deri_plus_Phi_psi-deri_minus_Phi_psi)

%RHS
rhs_psi = src_psi - Hhat_psi;
rhs_Pi_psi = src_Pi_psi - Hhat_Pi_psi;
rhs_Phi_psi = src_Phi_psi - Hhat_Phi_psi;

%constraints
Cr_psi = avg_deri_psi - Phi_psi;

return
