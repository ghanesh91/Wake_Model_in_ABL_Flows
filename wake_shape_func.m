function [xi_dim,xi0_theta]= wake_shape_func(beta,A_star,R,x,y,z,theta_angle,t_hat)
    
    Nx=numel(x);
    Ny=numel(y);
    Nz=numel(z);
    xi0_tilde=R*sqrt(A_star);
    xi0_theta=zeros(Nx,Ny,Nz);
    t_hat_3D =zeros(Nx,Ny,Nz);
    
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                xi0_theta(i,j,k) = xi0_tilde * abs(cos(beta)) ./ sqrt(1 - (sin(beta) .* sin(theta_angle(i, j, k))).^2);
            end
        end
    end
    
    t_hat_3D(:,:,:) = permute(repmat(t_hat, [1, 1, Ny]), [1, 3, 2]);
    
    xi_hat=zeros(Nx,Ny,Nz);
    xi_dim=zeros(Nx,Ny,Nz);
    alpha=1.263;
    xi_hat(:,:,:) = 1-alpha.*((1/2 )*tanh(t_hat_3D.^2/4 /alpha).*cos(2*theta_angle) ...
                             -(1/4 )*tanh(t_hat_3D.^3/8 /alpha).*cos(3*theta_angle) ...
                             -(5/48)*tanh(t_hat_3D.^4/16/alpha).*cos(2*theta_angle) ...
                             +(7/48)*tanh(t_hat_3D.^4/16/alpha).*cos(4*theta_angle));
    xi_dim(:,:,:) = real(xi_hat.*xi0_theta);
   
   
end