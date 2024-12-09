function [theta_angle]= theta_func(x,y,z,zh,yc_xz)
   Nx=numel(x);
   Ny=numel(y);
   Nz=numel(z);
   theta_angle=zeros(Nx,Ny,Nz);
   for i=1:Nx
      for j=1:Ny
        for k=1:Nz  
            theta_angle(i,j,k)=atan2(z(k)-zh,y(j)-yc_xz(i,k));
            % if theta_angle(i,j,k)<0
            %    theta_angle(i,j,k)=theta_angle(i,j,k)+2*pi;
            % end
        end
       end
   end
end