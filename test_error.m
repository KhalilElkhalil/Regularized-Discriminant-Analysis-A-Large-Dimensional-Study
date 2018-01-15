clc;clear all;
tic; 
p_vect = 100:50:500;
Error = zeros(length(p_vect),1);
Num = 100;
true_error_rda = zeros(length(p_vect),1) ; % RDA true error
theory_error_rda = zeros(length(p_vect),1) ; 
Theory_error_rda = zeros(length(p_vect),1);
for i = 1:length(p_vect)
    p = p_vect(i);
%     n0 = 0.5* p;    % c = 1/2
%     n1 = 0.5* p;   

%     n0 = p;    % c = 1
%     n1 = p;

    n0 = 1.5*p;  %c = 3/2
    n1 = 1.5*p;
    n = n0+n1;
    pi0 = n0/n ; pi1 = n1/n ; % priors
    n0_test = 1000 ; n1_test = round(n0_test * n1/n0) ; n_test = n0_test + n1_test ;

    gamma = 0.5; lambda = 0.5;
    
     mu0 = ones(p,1) ; mu1 = mu0 + 2/sqrt(p)*ones(p,1) ; mu = mu0-mu1;  
%    mu0 = ones(p,1) ; mu1 = mu0 + 2/p^4*ones(p,1) ; mu = mu0-mu1;    

    Sigma0 = zeros(p,p);
    for k =1:p
        for iter =1:p
            Sigma0(k, iter) = .6^(abs(k-iter));
        end
    end
    k = round(sqrt(p)); 
    A = [eye(k) zeros(k, p-k); zeros(p-k, k) zeros(p-k, p-k)];
    Sigma1 = Sigma0 +2*A; 
    Sigma_r0 = sqrtm(Sigma0); Sigma_r1 = sqrtm(Sigma1);        
    err = 0;
    Err_rda_theory = 0;
    err_rda_Theory = 0; 
    error_rda = 0;
    parfor j = 1:Num
        
        [i, j]
    Z0_train = randn(p,n0) ; %
    Z1_train = randn(p,n1) ;
    Y0 = Sigma_r0 * Z0_train ; u0 = ones(n0,1)/n0 ;
    Y1 = Sigma_r1 * Z1_train ; u1 = ones(n1,1)/n1 ;
    Y0_train = mu0 * ones(n0,1)' + Sigma_r0 * Z0_train ; % All training data as columns of C0
    Y1_train = mu1 * ones(n1,1)' + Sigma_r1 * Z1_train ; % All training data as columns of C0
    x0_ = Y0_train * u0 ; x1_ = Y1_train * u1 ;   % sample means
    Sigma0_hat = 1/n0*(Y0_train - x0_ * ones(1,n0))*(Y0_train - x0_ * ones(1,n0)).' ; % SCM for C0
    Sigma1_hat = 1/n1*(Y1_train - x1_ * ones(1,n1))*(Y1_train - x1_ * ones(1,n1)).' ; % SCM for C1
    
    alpha0 = gamma*n0/(n0+lambda*n1);
    beta0 = gamma*n1*lambda/(n0+lambda*n1);
    alpha1 = gamma*lambda*n0/(n1+lambda*n0);
    beta1  = gamma*n1/(n1+lambda*n0);
    
    H0 = inv(alpha0*Sigma0_hat+beta0*Sigma1_hat+(1-gamma)*eye(p));
    H1 = inv(beta1*Sigma1_hat+alpha1*Sigma0_hat+(1-gamma)*eye(p));
    %% compute H0_bar, H1_bar
    Sigma_r0_tilde = sqrt(p*alpha0/n0)*Sigma_r0;
    Sigma_r1_tilde = sqrt(p*beta0/n1)*Sigma_r1;
    Sigma0_tilde = Sigma_r0_tilde*Sigma_r0_tilde;
    Sigma1_tilde = Sigma_r1_tilde*Sigma_r1_tilde;
    
    Sigma_r0_tilde_1 = sqrt(p*alpha1/n0)*Sigma_r0;
    Sigma_r1_tilde_1 = sqrt(p*beta1/n1)*Sigma_r1;
    Sigma0_tilde_1 = Sigma_r0_tilde_1*Sigma_r0_tilde_1;
    Sigma1_tilde_1 = Sigma_r1_tilde_1*Sigma_r1_tilde_1;
    
    Z0 = randn(p,n0) ;
    Z1 = randn(p,n1) ;
    Y0_tilde = Sigma_r0_tilde * Z0 ;
    Y1_tilde = Sigma_r1_tilde * Z1 ;
    w_tilde  = 1/sqrt(p)*[Y0_tilde,Y1_tilde];
    
    Y0_tilde_1 = Sigma_r0_tilde_1 * Z0 ;
    Y1_tilde_1 = Sigma_r1_tilde_1 * Z1 ;
    w_tilde_1  = 1/sqrt(p)*[Y0_tilde_1,Y1_tilde_1];
    z = gamma - 1;
    
    c = p/n;
    c0 = n0/n;
    c1 = n1/n;
    
    g0 = 1; g1 = 1;
    g0_tilde = 1; g1_tilde = 1;
    err0 = 1; err1 = 1;
    while  err0>1e-6 || err1>1e-6
        a0 = g0;
        a1 = g1;
        g0_tilde = -1/z*1/p*trace(Sigma0_tilde*inv((eye(p)+c0*g0*Sigma0_tilde+c1*g1*Sigma1_tilde)));
        g1_tilde = -1/z*1/p*trace(Sigma1_tilde*inv((eye(p)+c0*g0*Sigma0_tilde+c1*g1*Sigma1_tilde)));
        g0 = -1/z*1/c*1/(1+g0_tilde);
        g1 = -1/z*1/c*1/(1+g1_tilde);
        err0 = abs(a0 - g0)^2;
        err1 = abs(a1 - g1)^2;
    end
    
    H0_bar = -1/z*inv((eye(p)+c0*g0*Sigma0_tilde+c1*g1*Sigma1_tilde));
    
    g0_ = 1; g1_ = 1;
    g0_tilde_ = 1; g1_tilde_ = 1;
    err0 = 1; err1 = 1;
    while  err0>1e-6 || err1>1e-6
        a0 = g0_;
        a1 = g1_;
        g0_tilde_ = -1/z*1/p*trace(Sigma0_tilde_1*inv((eye(p)+c0*g0_*Sigma0_tilde_1+c1*g1_*Sigma1_tilde_1)));
        g1_tilde_ = -1/z*1/p*trace(Sigma1_tilde_1*inv((eye(p)+c0*g0_*Sigma0_tilde_1+c1*g1_*Sigma1_tilde_1)));
        g0_ = -1/z*1/c*1/(1+g0_tilde_);
        g1_ = -1/z*1/c*1/(1+g1_tilde_);
        err0 = abs(a0 - g0_)^2;
        err1 = abs(a1 - g1_)^2;
    end
    
    H1_bar = -1/z*inv((eye(p)+c0*g0_*Sigma0_tilde_1+c1*g1_*Sigma1_tilde_1));

    Z0_test = randn(p,n0_test) ; % testing data for C0
    Z1_test = randn(p,n1_test) ; % testing data for C1
    Y0_test = mu0 * ones(n0_test,1)' + Sigma_r0 * Z0_test ; % All testing data as columns of C0
    Y1_test = mu1 * ones(n1_test,1)' + Sigma_r1 * Z1_test ; % All testing data as columns of C1

    %% RDA true error
    s0 = eig(H0); 
    s1 = eig(H1);
        for ii = 1:n0_test
            x = Y0_test(:,ii) ;
            %% RDA
            RDA0 = 0.5 * sum(log(s0)) - 0.5 * (x-x0_).' * H0 * (x-x0_) + log(pi0) ; % quadratic disc. for C0
            RDA1 = 0.5 * sum(log(s1)) - 0.5 * (x-x1_).' * H1 * (x-x1_) + log(pi1) ; % quadratic disc. for C1

            if RDA0 < RDA1
                error_rda = error_rda + 1 ;
            end            
             
            
        end
        for ii = 1:n1_test
            x = Y1_test(:,ii) ;
            %% RDA
            RDA0 = 0.5 * sum(log(s0)) - 0.5 * (x-x0_).' * H0 * (x-x0_) + log(pi0) ; % quadratic disc. for C0
            RDA1 = 0.5 * sum(log(s1)) - 0.5 * (x-x1_).' * H1 * (x-x1_) + log(pi1) ; % quadratic disc. for C1

            if RDA1 < RDA0
                error_rda = error_rda + 1 ;
            end
            
        end


 %%  RDA theory error    
%         B0 = Sigma_r0*(H1 - H0)*Sigma_r0; 
%         theory_y0 = Sigma_r0*(H1*(mu0 - x1_) - H0*(mu0 - x0_));
%         xi0 = -log(det(H0)/det(H1)) + (mu0 - x0_)'* H0 * (mu0 - x0_) - (mu0 - x1_)'*H1*(mu0 - x1_)+2*log(pi1/pi0);         
%         B1 = Sigma_r1*(H1 - H0)*Sigma_r1; 
%         theory_y1 = Sigma_r1*(H1*(mu1 - x1_) - H0*(mu1 - x0_));
%         xi1 = -log(det(H0)/det(H1)) + (mu1 - x0_)'* H0 * (mu1 - x0_) - (mu1 - x1_)'*H1*(mu1 - x1_)+2*log(pi1/pi0);      
%         Arg0 = (xi0-trace(B0))/sqrt(2*trace(B0*B0)+4*theory_y0'*theory_y0);
%         Arg1 = (xi1-trace(B1))/sqrt(2*trace(B1*B1)+4*theory_y1'*theory_y1)
%         Eps_0 = 1-qfunc(Arg0);
%         Eps_1 = 1-qfunc(-Arg1);
%         Eps = pi0*Eps_0+pi1*Eps_1;
%         Err_rda_theory = Err_rda_theory+Eps;
       %% RDA asymptotic theory error    
       [trace_B0B0,trace_B1B1] = func_traceBB(p,n0,n1,n,alpha0,beta0,alpha1,beta1,Sigma0,Sigma1,H0,H1,H0_bar,H1_bar);
       [xi_0,xi_1]=func_xi(p,n0,n1,n,alpha0,beta0,alpha1,beta1,Sigma0,Sigma1,H0,H1,H0_bar,H1_bar,mu,gamma,lambda);  
       [traceB0,traceB1] = func_traceB(p,Sigma0,Sigma1,H0,H1,H0_bar,H1_bar);
       [y0y0,y1y1] = func_yi(p,n0,n1,z,alpha0,beta0,alpha1,beta1,c,c0,c1,g0_,g1_,g0,g1,Sigma0,Sigma1,Sigma0_tilde_1,Sigma1_tilde_1,Sigma0_tilde,Sigma1_tilde,H0_bar,H1_bar,mu);
       
        arg0 = (xi_0 - traceB0)/sqrt(2*trace_B0B0+4*y0y0);
        arg1 = (xi_1 - traceB1)/sqrt(2*trace_B1B1+4*y1y1);      

        eps_0 = 1-qfunc(arg0);
        eps_1 = 1-qfunc(-arg1);
        eps = pi0*eps_0 + pi1*eps_1; 
        err_rda_Theory = err_rda_Theory + eps;
    end
theory_error_rda(i) = err_rda_Theory/Num;
true_error_rda(i) = error_rda/n_test/Num;
%Theory_error_rda(i) = Err_rda_theory/Num;
end
toc
figure;
plot(p_vect,theory_error_rda,'b','LineWidth',2) ; hold on
plot(p_vect,true_error_rda,'r','LineWidth',2) ; hold on
%plot(p_vect,Theory_error_rda,'g') 
set(gca,'FontSize',12);

xlabel({'p'},'Fontsize',12)
ylabel({'Classification Error'},'Fontsize',12) ;

legend({'asymptotic error','empirical error'},'Fontsize',12);
title({'c'},'Fontsize',12);

grid on;
