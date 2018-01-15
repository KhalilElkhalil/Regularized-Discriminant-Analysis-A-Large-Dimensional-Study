function [theory_trace_B0B0,theory_trace_B1B1] = func_traceBB(p,n0,n1,n,alpha0,beta0,alpha1,beta1,Sigma0,Sigma1,H0,H1,H0_bar,H1_bar)

delta0 = 1/n0*trace(Sigma0*H0_bar);
delta1 = 1/n0*trace(Sigma0*H1_bar);
delta0_bar = 1/n1*trace(Sigma1*H0_bar);
delta1_bar = 1/n1*trace(Sigma1*H1_bar);
delta0_tilde = alpha1/(1+alpha1*delta1);
delta1_tilde = beta1/(1+beta1*delta1_bar);
delta00_tilde = alpha0/(1+alpha0*delta0);
delta11_tilde = beta0/(1+beta0*delta0_bar);

M0 = delta00_tilde^2/(p*n0);
T0 = delta11_tilde^2/(p*n1);

% test_Sig0_H0_bar_H0_bar = (trace(Sigma0*H0_bar*Sigma0*H0_bar)*(1-p*T0*trace(Sigma1*H0_bar*Sigma1*H0_bar))+...
%     p*T0*trace(Sigma1*H0_bar*Sigma0*H0_bar)*trace(Sigma1*H0_bar*Sigma0*H0_bar))/(p-p^2*T0*trace(Sigma1*H0_bar*Sigma1*H0_bar)*...
%     (1+p*M0*trace(Sigma0*H0_bar*Sigma0*H0_bar))-p^3*T0*M0*trace(Sigma1*H0_bar*Sigma0*H0_bar*trace(Sigma0*H0_bar*Sigma1*H0_bar)));
% 
% test_Sig1_H0_bar_H0_bar=(trace(Sigma1*H0_bar*Sigma1*H0_bar)*(1-p*T0*trace(Sigma1*H0_bar*Sigma1*H0_bar))+...
%     p*T0*trace(Sigma1*H0_bar*Sigma0*H0_bar)*trace(Sigma0*H0_bar*Sigma1*H0_bar))/(p-p^2*T0*trace(Sigma1*H0_bar*Sigma1*H0_bar)*...
%     (1+p*M0*trace(Sigma0*H0_bar*Sigma0*H0_bar))-p^3*T0*M0*trace(Sigma1*H0_bar*Sigma0*H0_bar*trace(Sigma0*H0_bar*Sigma1*H0_bar)));

test_Sig0_H0_bar_H0_bar = 1/p*(1/p*trace(Sigma0*H0_bar*Sigma0*H0_bar)+T0*trace(Sigma1*H0_bar*Sigma0*H0_bar)*...
1/p*trace(Sigma1*H0_bar*Sigma0*H0_bar)/(1/p-T0*trace(Sigma1*H0_bar*Sigma1*H0_bar)))/(1/p-M0*trace(Sigma0*H0_bar*Sigma0*H0_bar)...
-T0*trace(Sigma1*H0_bar*Sigma0*H0_bar)*M0*trace(Sigma0*H0_bar*Sigma1*H0_bar)/(1/p-T0*trace(Sigma1*H0_bar*Sigma1*H0_bar)));

test_Sig1_H0_bar_H0_bar = 1/p*(1/p*trace(Sigma1*H0_bar*Sigma1*H0_bar)+M0*trace(Sigma0*H0_bar*Sigma1*H0_bar)*...
1/p*trace(Sigma0*H0_bar*Sigma1*H0_bar)/(1/p-M0*trace(Sigma0*H0_bar*Sigma0*H0_bar)))/(1/p-T0*trace(Sigma1*H0_bar*Sigma1*H0_bar)...
-M0*trace(Sigma0*H0_bar*Sigma1*H0_bar)*T0*trace(Sigma1*H0_bar*Sigma0*H0_bar)/(1/p-M0*trace(Sigma0*H0_bar*Sigma0*H0_bar)));

%% compute asmyptotic equivalence of H1*Sigma0*H1
M1 = delta0_tilde^2/(p*n0);
T1 = delta1_tilde^2/(p*n1);

% test_Sig0_H1_bar_H1_bar = (trace(Sigma0*H1_bar*Sigma0*H1_bar)*(1-p*T1*trace(Sigma1*H1_bar*Sigma1*H1_bar))+...
%     p*T1*trace(Sigma1*H1_bar*Sigma0*H1_bar)*trace(Sigma1*H1_bar*Sigma0*H1_bar))/(p-p^2*T1*trace(Sigma1*H1_bar*Sigma1*H1_bar)*...
%     (1+p*M1*trace(Sigma0*H1_bar*Sigma0*H1_bar))-p^3*T1*M1*trace(Sigma1*H1_bar*Sigma0*H1_bar*trace(Sigma0*H1_bar*Sigma1*H1_bar)));

test_Sig0_H1_bar_H1_bar = 1/p*(1/p*trace(Sigma0*H1_bar*Sigma0*H1_bar)+T1*trace(Sigma1*H1_bar*Sigma0*H1_bar)*...
1/p*trace(Sigma1*H1_bar*Sigma0*H1_bar)/(1/p-T1*trace(Sigma1*H1_bar*Sigma1*H1_bar)))/(1/p-M1*trace(Sigma0*H1_bar*Sigma0*H1_bar)...
-T1*trace(Sigma1*H1_bar*Sigma0*H1_bar)*M1*trace(Sigma0*H1_bar*Sigma1*H1_bar)/(1/p-T1*trace(Sigma1*H1_bar*Sigma1*H1_bar)));

% test_Sig1_H1_bar_H1_bar = (trace(Sigma1*H1_bar*Sigma1*H1_bar)*(1-p*T1*trace(Sigma1*H1_bar*Sigma1*H1_bar))+...
%     p*T1*trace(Sigma0*H1_bar*Sigma1*H1_bar)*trace(Sigma1*H1_bar*Sigma0*H1_bar))/(p-p^2*T1*trace(Sigma1*H1_bar*Sigma1*H1_bar)*...
%     (1+p*M1*trace(Sigma0*H1_bar*Sigma0*H1_bar))-p^3*T1*M1*trace(Sigma1*H1_bar*Sigma0*H1_bar*trace(Sigma0*H1_bar*Sigma1*H1_bar)));

test_Sig1_H1_bar_H1_bar = 1/p*(1/p*trace(Sigma1*H1_bar*Sigma1*H1_bar)+M1*trace(Sigma0*H1_bar*Sigma1*H1_bar)*...
1/p*trace(Sigma0*H1_bar*Sigma1*H1_bar)/(1/p-M1*trace(Sigma0*H1_bar*Sigma0*H1_bar)))/(1/p-T1*trace(Sigma1*H1_bar*Sigma1*H1_bar)...
-M1*trace(Sigma0*H1_bar*Sigma1*H1_bar)*T1*trace(Sigma1*H1_bar*Sigma0*H1_bar)/(1/p-M1*trace(Sigma0*H1_bar*Sigma0*H1_bar)));

%% compute asmyptotic equivalence of H0*Sigma*H1

M = alpha0*delta0_tilde/(p*n0*(1+alpha0*delta0));
T = beta0*delta1_tilde/(p*n1*(1+beta0*delta0_bar));

% test_Sig0_H0_bar_H1_bar = (trace(Sigma0*H0_bar*Sigma0*H1_bar)*(1-p*T*trace(Sigma1*H1_bar*Sigma1*H0_bar))+...
%     p*T*trace(Sigma1*H1_bar*Sigma0*H0_bar)*trace(Sigma1*H0_bar*Sigma0*H1_bar))/(p-p^2*T*trace(Sigma1*H1_bar*Sigma1*H0_bar)*...
%     (1+p*M*trace(Sigma0*H1_bar*Sigma0*H0_bar))-p^3*T*M*trace(Sigma1*H1_bar*Sigma0*H0_bar*trace(Sigma0*H1_bar*Sigma1*H0_bar)));

test_Sig0_H0_bar_H1_bar = 1/p*(1/p*trace(Sigma0*H0_bar*Sigma0*H1_bar)+T*trace(Sigma1*H1_bar*Sigma0*H0_bar)*...
1/p*trace(Sigma1*H0_bar*Sigma0*H1_bar)/(1/p-T*trace(Sigma1*H1_bar*Sigma1*H0_bar)))/(1/p-M*trace(Sigma0*H1_bar*Sigma0*H0_bar)...
-T*trace(Sigma1*H1_bar*Sigma0*H0_bar)*M*trace(Sigma0*H1_bar*Sigma1*H0_bar)/(1/p-T*trace(Sigma1*H1_bar*Sigma1*H0_bar)));

test_H0H1 = 1/p*trace(Sigma0*H0*Sigma0*H1);
% disp(['empirical sigma0*H0*sigma0*H1:',num2str(test_H0H1)]);
% disp(['test sigma0*H0_bar*sigma0*H1_bar:',num2str(test_Sig0_H0_bar_H1_bar)]);

% test_Sig1_H0_bar_H1_bar = (trace(Sigma1*H0_bar*Sigma1*H1_bar)*(1-p*T*trace(Sigma1*H1_bar*Sigma1*H0_bar))+...
%     p*T*trace(Sigma1*H1_bar*Sigma0*H0_bar)*trace(Sigma0*H0_bar*Sigma1*H1_bar))/(p-p^2*T*trace(Sigma1*H1_bar*Sigma1*H0_bar)*...
%     (1+p*M*trace(Sigma0*H1_bar*Sigma0*H0_bar))-p^3*T*M*trace(Sigma1*H1_bar*Sigma0*H0_bar*trace(Sigma0*H1_bar*Sigma1*H0_bar)));

test_Sig1_H0_bar_H1_bar = 1/p*(1/p*trace(Sigma1*H0_bar*Sigma1*H1_bar)+M*trace(Sigma0*H1_bar*Sigma1*H0_bar)*...
1/p*trace(Sigma0*H0_bar*Sigma1*H1_bar)/(1/p-M*trace(Sigma0*H1_bar*Sigma0*H0_bar)))/(1/p-T*trace(Sigma1*H1_bar*Sigma1*H0_bar)...
-M*trace(Sigma0*H1_bar*Sigma1*H0_bar)*T*trace(Sigma1*H1_bar*Sigma0*H0_bar)/(1/p-M*trace(Sigma0*H1_bar*Sigma0*H0_bar)));

test_Sig1_H0H1 = 1/p*trace(Sigma1*H0*Sigma1*H1);
% disp(['empirical sigma1*H0*sigma1*H1:',num2str(test_Sig1_H0H1)]);
% disp(['test sigma1*H0_bar*sigma0*_H1_bar:',num2str(test_Sig1_H0_bar_H1_bar)]);


theory_trace_B0B0 = test_Sig0_H0_bar_H0_bar+test_Sig0_H1_bar_H1_bar-2*test_Sig0_H0_bar_H1_bar;
emprical_trace_B0B0 = 1/p*trace(Sigma0*(H1-H0)*Sigma0*(H1-H0));

%  disp(['empirical traceB0B0:',num2str(emprical_trace_B0B0)]);
%  disp(['theory traceB0B0:',num2str(theory_trace_B0B0)]);

 theory_trace_B1B1 = test_Sig1_H0_bar_H0_bar+test_Sig1_H1_bar_H1_bar-2*test_Sig1_H0_bar_H1_bar;
%emprical_trace_B1B1 = 1/p*trace(Sigma1*(H1-H0)*Sigma1*(H1-H0));
emprical_trace_B1B1 = 1/p*trace(Sigma1*H1*Sigma1*H1)+1/p*trace(Sigma1*H0*Sigma1*H0)-2/p*trace(Sigma1*H0*Sigma1*H1);


%  disp(['empirical traceB1B1:',num2str(emprical_trace_B1B1)]);
%  disp(['theory traceB1B1:',num2str(theory_trace_B1B1)]);

end