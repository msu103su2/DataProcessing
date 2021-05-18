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
       
   end
   
   methods(Static)
       function theta = guoyAngle(w0, lambda, z)
           theta = atan(z/(pi*w0^2/lambda));
       end
       
       function result = zR(w0, lambda)
           result = pi*w0^2/lambda;
       end
   end
end