function [Cx,A_star,xo,kw]= Cx_func(beta,Ct,R,u_star,Uh,zh,h,x,A1,B1,model_type)
    
    D         = 2*R;
    A_star    = (1+sqrt(1-Ct*cos(beta)^2))./ ...
                (2*sqrt(1-Ct*cos(beta)^2));
    a_ind     = 0.5*(1 - sqrt(1-Ct*(cos(beta))^2));
    xi0_tilde = R*sqrt(A_star);

    Iu_zh   = sqrt(A1*log(zh/h)+B1)*(u_star/Uh);
    if (model_type ==1)
         kw      = 0.05;
    end
    if (model_type==2)
        kw = 0.075;
    end
    if (model_type==3)
        kw=( 0.021^6 + (0.33*Iu_zh)^6 )^(1/6);
    end

    sigma2  = (kw*x+0.4*xi0_tilde).* ...
              (kw*x+0.4*xi0_tilde*cos(beta));

    xo      = -0.4*xi0_tilde/kw+(0.2/kw)*(sqrt(25*D^2*Ct*cos(beta)^3/(8*(1-(1-2*a_ind)^2)) +xi0_tilde^2*(1-cos(beta))^2 )+ ...
              xi0_tilde*(1-cos(beta)));  

    Cx      = 1-sqrt(1-(Ct*(cos(beta))^3)./(8*sigma2./D^2));

    Cx(x<xo)=2*a_ind; 
end