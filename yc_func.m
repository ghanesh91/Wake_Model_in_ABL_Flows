function [yc_majid,yc_CVP_majid,yc_g_majid]=yc_func(t_hat,R,A_star,zh,x,z)
    Nx = numel(x);
    Nz = numel(z);


        yc_hat_majid = zeros(Nx,Nz);
      yc_hat_g_majid = zeros(Nx,Nz);
    yc_hat_CVP_majid = zeros(Nx,Nz);
            yc_majid = zeros(Nx,Nz);
          yc_g_majid = zeros(Nx,Nz);
        yc_CVP_majid = zeros(Nx,Nz);
           xi0_tilde = R*sqrt(A_star);
    for i=1:Nx
        for k=1:Nz
            num = (pi-1)*abs(t_hat(i,k))^3   + 2*sqrt(3)*pi^2*t_hat(i,k)^2    + ...
                  48*(pi-1)^2*abs(t_hat(i,k));
            den = 2*pi*(pi-1)*(t_hat(i,k))^2 + 4*sqrt(3)*pi^2*abs(t_hat(i,k)) + 96*(pi-1)^2;
            yc_hat_CVP_majid(i,k) =  (num./den).*sign(t_hat(i,k));
            yc_hat_g_majid(i,k)   = -(2/pi)*(t_hat(i,k)/(((z(k)+zh)/xi0_tilde)^2-1));   
            yc_hat_majid(i,k)     = yc_hat_CVP_majid(i,k)+yc_hat_g_majid(i,k);   
        end
    end
   
    yc_majid(:,:)     =     real(yc_hat_majid*xi0_tilde);
    yc_g_majid(:,:)   =   real(yc_hat_g_majid*xi0_tilde);
    yc_CVP_majid(:,:) = real(yc_hat_CVP_majid*xi0_tilde);
    
end