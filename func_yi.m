function [theory_y0y0,theory_y1y1] = func_yi(p,n0,n1,z,alpha0,beta0,alpha1,beta1,c,c0,c1,g0_,g1_,g0,g1,Sigma0,Sigma1,Sigma0_tilde_1,Sigma1_tilde_1,Sigma0_tilde,Sigma1_tilde,H0_bar,H1_bar,mu)

%% i = 0, compute asmyptotic equivalence of H1*Sigma0*H1
C = [c0,c1];
g_=[g0_,g1_];
Sigma_tilde_ = cell(2,1);Sigma_tilde_{1,1} = Sigma0_tilde_1;Sigma_tilde_{2,1} = Sigma1_tilde_1;
Omega1 = zeros(2,2); 
for a =1:2
    for b =1:2
        Omega1(a, b) = c*C(b)*z^2*g_(a)^2*1/p*trace(Sigma_tilde_{a,1}*H1_bar*Sigma_tilde_{b,1}*H1_bar); 
    end
end

R_ = inv(eye(2)-Omega1)*Omega1;
R_00 = c0/c0*R_(1,1);
R_10 = c1/c0*R_(2,1);
H1S0H1 = H1_bar*Sigma0_tilde_1*H1_bar+R_00*H1_bar*Sigma0_tilde_1*H1_bar+R_10*H1_bar*Sigma1_tilde_1*H1_bar;
if alpha1 == 0
    Q0_bar = H1_bar*Sigma0*H1_bar;
else
    Q0_bar = H1S0H1/(p*alpha1/n0);
end

R_01 = c0/c1*R_(1,2);
R_11 = c1/c1*R_(2,2);
H1S1H1 = H1_bar*Sigma1_tilde_1*H1_bar+R_01*H1_bar*Sigma0_tilde_1*H1_bar+R_11*H1_bar*Sigma1_tilde_1*H1_bar;
H1S1H1 = H1S1H1/(p*beta1/n1);
theory_y0 = 1/p*mu'*Q0_bar*mu;


%% compute asmyptotic equivalence of H0*Sigma1*H0
C = [c0,c1];
g=[g0,g1];
Sigma_tilde = cell(2,1);Sigma_tilde{1,1} = Sigma0_tilde;Sigma_tilde{2,1} = Sigma1_tilde;
Omega = zeros(2,2); 
for a =1:2
    for b =1:2
        Omega(a, b) = c*C(b)*z^2*g(a)^2*1/p*trace(Sigma_tilde{a,1}*H0_bar*Sigma_tilde{b,1}*H0_bar); 
    end
end

R = inv(eye(2)-Omega)*Omega;
R01 = c0/c1*R(1,2);
R11 = c1/c1*R(2,2);
H0S1H0 = H0_bar*Sigma1_tilde*H0_bar+R01*H0_bar*Sigma0_tilde*H0_bar+R11*H0_bar*Sigma1_tilde*H0_bar;
if beta0 == 0
    Q1_bar = H0_bar*Sigma1*H0_bar;
else
    Q1_bar = H0S1H0/(p*beta0/n1);
end

R00 = c0/c0*R(1,1);
R10 = c1/c0*R(2,1);
HS0H = H0_bar*Sigma0_tilde*H0_bar+R00*H0_bar*Sigma0_tilde*H0_bar+R10*H0_bar*Sigma1_tilde*H0_bar;
H0S0H0 = HS0H/(p*alpha0/n0);
theory_y1 = 1/p*mu'*Q1_bar*mu;
theory_y0y0 = theory_y0+1/p*1/n1*trace(Sigma1*Q0_bar)+1/p*1/n0*trace(Sigma0*H0S0H0);
theory_y1y1 = theory_y1+1/p*1/n1*trace(Sigma1*H1S1H1)+1/p*1/n0*trace(Sigma0*Q1_bar);

end
