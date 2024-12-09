function [t_hat_xz]= t_hat_xz_func(beta,Ct,R,A_star,u_star,Uh,zh,h,x,z,U,V)
    Nx = numel(x);
    Nz = numel(z);
    t_hat_xz = zeros(Nx,Nz);
    U_res = sqrt(U.^2+0*V.^2);
    for i = 1:Nx
        for k = 1:Nz
            t_hat_xz(i,k)=-1.44*(Uh/u_star)*(1/sqrt(A_star))*Ct*(cos(beta)^2)*sin(beta) ...
                          .*(1-exp( -0.35*(u_star/U_res(k))*(x(i)/R)  ) );
        end
    end
end