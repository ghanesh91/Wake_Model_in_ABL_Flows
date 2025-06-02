%% ABL velocity model from Narasimhan et al. (BLM, 2024)

function [Ug,Vg,ustar,h,mu,epsilon]=GDL_ABL_model(G,fc,z0,Cr,Theta0,lapse_rate,eps_tol)
%Definition of model constants
%Matching height = 20% of ABL height
cm=0.2;

%von-Karman constant
kappa=0.41;

%Fitting constants from LES (see Narasimhan et al., BLM, 2024)
C_TN  = 0.5;
C_CN  = 1.6;
C_NS  = 0.78;
Cg    = 1.43;
Gamma = 0.83;


g     = 9.81;%m/s^2
N_infty = sqrt(g/Theta0*lapse_rate);
muN   = N_infty/fc;
%Iterative estimation of h, ustar, Ug, Vg

%Step a-Guess values of h, ustar
h0      = 1000;%km
ustar0  = 0.1;%m/s
epsilon = 1;

while epsilon > eps_tol 
    %Step b-Compute next iteration of ustar
    xi_hat_0 = z0*fc/ustar0;
    Ro0       = 1/xi_hat_0;
    h_hat_0  = h0*fc/ustar0;
    mu0      = -(g*Cr*h_hat_0)/(ustar0*fc^2*Theta0);
    
    g_func_cm       = Cg*(1-exp(-cm/Gamma));
    g_prime_func_cm = (Cg/Gamma/h_hat_0)*exp(-cm/Gamma);

    B0  = 3*kappa/(2*h_hat_0);
    A0  = -log(cm*h_hat_0) ...
          - kappa*((5*mu0+0.3*muN)*(cm*h_hat_0-xi_hat_0) ...
                   +g_prime_func_cm*(1-cm)^(3/2)-g_func_cm*(3/2/h_hat_0)*sqrt(1-cm));
    ustar1 = kappa*G/sqrt((log(Ro0)-A0)^2 + B0^2);

    %Step-c: Calculate next iteration of ABL height
    h1 = (ustar1/fc)*(1/C_TN^2 ...
                    +muN/C_CN^2 ...
                    +(1/C_NS^2)*(g/Theta0*(-Cr)*h0)/(ustar1^2*fc))^(-1/2);

    epsilon = max(abs(h1-h0)/h0,abs(ustar1-ustar0)/ustar0);
    if epsilon<eps_tol
        ustar = ustar1;
        h     = h1;
        xi_hat_0 = z0*fc/ustar;
        Ro       = 1/xi_hat_0;
        h_hat    = h*fc/ustar;
        mu       = -(g*Cr*h_hat)/(ustar*fc^2*Theta0);
        g_prime_func_cm = (Cg/Gamma/h_hat)*exp(-cm/Gamma);
        B  = 3*kappa/(2*h_hat);
        A  =  -log(cm*h_hat) ...
              - kappa*((5*mu+0.3*muN)*(cm*h_hat-xi_hat_0) ...
              + g_prime_func_cm*(1-cm)^(3/2)-g_func_cm*(3/2/h_hat)*sqrt(1-cm));
       
        Ug    = (ustar/kappa)*(log(Ro)-A);
        Vg    = -ustar/kappa*B;
    else
        h0=h1;
        ustar0=ustar1;
    end
end
end