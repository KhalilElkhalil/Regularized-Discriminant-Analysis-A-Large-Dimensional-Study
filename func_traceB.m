function [theory_traceB0,theory_traceB1] = func_traceB(p,Sigma0,Sigma1,H0,H1,H0_bar,H1_bar)
theory_traceB0 = 1/sqrt(p)*trace(Sigma0*(H1_bar-H0_bar));
emprical_traceB0 = 1/sqrt(p)*trace(Sigma0*(H1-H0));
% disp(['empirical traceB0:',num2str(emprical_traceB0)]);
% disp(['theory traceB0:',num2str(theory_traceB0)]);

theory_traceB1 = 1/sqrt(p)*trace(Sigma1*(H1_bar-H0_bar));
emprical_traceB1 = 1/sqrt(p)*trace(Sigma1*(H1-H0));
% disp(['empirical traceB1:',num2str(emprical_traceB1)]);
% disp(['theory traceB1:',num2str(theory_traceB1)]);
end
