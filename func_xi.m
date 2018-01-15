function [xi_0,xi_1]= func_xi(p,n0,n1,n,alpha0,beta0,alpha1,beta1,Sigma0,Sigma1,H0,H1,H0_bar,H1_bar,mu,gamma,lambda)

delta0 = 1/n0*trace(Sigma0*H0_bar);
delta1 = 1/n0*trace(Sigma0*H1_bar);
delta0_bar = 1/n1*trace(Sigma1*H0_bar);
delta1_bar = 1/n1*trace(Sigma1*H1_bar);

sim_logH0_H1 = n0/sqrt(p)*log((1+alpha0*delta0)/(1+alpha1*delta1))+n1/sqrt(p)*log((1+beta0*delta0_bar)/(1+beta1*delta1_bar))+...
    1/sqrt(p)*log(det(H1_bar)/...
    det(H0_bar))+1/sqrt(p)*(alpha1*delta1*n0/(1+alpha1*delta1)+beta1*delta1_bar*n1/(1+beta1*delta1_bar)-alpha0*delta0*n0/(1+alpha0*delta0)-beta0*delta0_bar*n1/(1+beta0*delta0_bar));


%% i = 0
theory_mu0_term = 1/sqrt(p)*(delta0-delta1_bar-mu'*H1_bar*mu);
% emprical_mu0_term = 1/sqrt(p)*((mu0-x0_)'*H0*(mu0-x0_)-(mu0-x1_)'*H1*(mu0-x1_));
% disp(['empirical mu0 term:',num2str(emprical_mu0_term)]);
% disp(['theory mu0 term:',num2str(theory_mu0_term)]);
%% i = 1
theory_mu1_term = 1/sqrt(p)*(delta0-delta1_bar+mu'*H0_bar*mu);
% emprical_mu1_term = 1/sqrt(p)*((mu1-x0_)'*H0*(mu1-x0_)-(mu1-x1_)'*H1*(mu1-x1_));
% disp(['empirical mu1 term:',num2str(emprical_mu1_term)]);
% disp(['theory mu1 term:',num2str(theory_mu1_term)]);

% theory_traceB0 = 1/sqrt(p)*trace(Sigma0*(H1_bar-H0_bar));
% emprical_traceB0 = 1/sqrt(p)*trace(Sigma0*(H1-H0));
% disp(['empirical traceB0:',num2str(emprical_traceB0)]);
% disp(['theory traceB0:',num2str(theory_traceB0)]);

xi_0 = sim_logH0_H1+theory_mu0_term;
xi_1 = sim_logH0_H1+theory_mu1_term;
end