classdef HGBeam < handle
   properties
      w0x;
      w0y;
      lambda;
      m;
      n;
      zRx;
      zRy;
      v;
   end
   methods
       function obj = HGBeam(m, n, lambda, w0x, w0y)
           obj.m = m;
           obj.n = n;
           obj.lambda = lambda;
           obj.w0x = w0x;
           obj.w0y = w0y;
           obj.zRx = pi*w0x^2/lambda;
           obj.zRy = pi*w0y^2/lambda;
       end
 
       function result = wx(obj, z)
           result = obj.w0x*sqrt(1+(z./obj.zRx).^2);
       end
       
       function result = wy(obj, z)
           result = obj.w0y*sqrt(1+(z./obj.zRy).^2);
       end
       
       function result = u(obj, x, y, z)
           wx = obj.wx(z);
           wy = obj.wy(z);
           result = (2/pi)^(1/2)/sqrt(2^(obj.m+obj.n)*factorial(obj.m)*factorial(obj.n)*wx*wy)...
               .*exp(-x.^2/wx^2-y.^2/wy^2).*hermiteH(obj.m,sqrt(2)*x/wx).*hermiteH(obj.n,sqrt(2)*y/wy);
       end
       
       function result = phi(obj, x, y, z)
           result = exp(-1i*(obj.k*z + obj.k*z/2.*(x.^2/obj.zRx^2 + y.^2/obj.zRy^2) -...
               (obj.m+1/2)*atan(z/obj.zRx) - (obj.n+1/2)*atan(z/obj.zRy)));
       end
       
       function result = FEu(obj, x, y, z, k)
           wx = obj.wx(z);
           wy = obj.wy(z);
           if z == 0
            invert_R = 0;
           else
               invert_R = (z^2+obj.zRx^2)/z;
           end
           xpart = exp(-x.^2/wx^2).*hermiteH(obj.m,sqrt(2)*x/wx);
           ypart = exp(-y.^2/wy^2).*hermiteH(obj.n,sqrt(2)*y/wy);
           Lx = length(x);
           Ly = length(y);
           x = repmat(x, Ly, 1);
           y = repmat(y',1, Lx);
           result = (2/pi)^(1/2)/sqrt(2^(obj.m+obj.n)*factorial(obj.m)*factorial(obj.n)*wx*wy)...
               *kron(xpart, ypart').*exp(-1i*k*(x.^2+y.^2)*invert_R)*exp((obj.m+obj.n+1)*1i*atan(z/obj.zRy));
       end
       
   end
   
   methods(Static)
       function theta = guoyAngle(w0, lambda, z)
           theta = atan(z./(pi*w0^2/lambda));
       end
       
       function result = zR(w0, lambda)
           result = pi*w0^2/lambda;
       end
       
       function result = wz(w0, lambda, z)
           zR = pi*w0^2/lambda;
           result = w0*sqrt(1+(z./zR).^2);
       end
       
       function result = w0(w, R, lambda)
           result = sqrt(w^2/(1+(pi*w^2/(lambda*R))^2));
       end
       
   end
end