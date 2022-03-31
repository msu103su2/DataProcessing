classdef OLExperiment < handle
    properties
        sp;
        device;
        F0;
        F1;
        laser;
        PD;
        w0;
        lenf;
        consts = consts;
        T = 300;
        f;
    end
    
    methods
        function obj = OLExperiment(sp, device, laser, w0, lenf, PD)
            obj.sp = sp;
            obj.device = device;
            obj.laser = laser;
            obj.w0 = w0;
            obj.lenf = lenf;
            obj.PD = PD;
        end
        
        function CalF(obj, spdx, d)
            HGB0 = HGBeam(0, 0, obj.laser.lambda, obj.w0, obj.w0);
            HGB1 = HGBeam(0, 1, obj.laser.lambda, obj.w0, obj.w0);
            wz = HGB0.wx(obj.lenf);
            oldoffset = obj.sp.offset;
            oldd = obj.sp.d;
            
            obj.sp.offset = spdx/wz*obj.w0;
            obj.sp.d = d/wz*obj.w0;
            obj.F0 = obj.sp.F_onlyx(HGB0);
            obj.F1 = obj.sp.F_onlyx(HGB1);
            obj.sp.offset = oldoffset;
            obj.sp.d = oldd;
        end
        
        function re = CalFmn(obj, spdx, d, m, n)
            HGB = HGBeam(m, n, obj.laser.lambda, obj.w0, obj.w0);
            wz = HGB.wx(obj.lenf);
            oldoffset = obj.sp.offset;
            oldd = obj.sp.d;
            
            obj.sp.offset = spdx/wz*obj.w0;
            obj.sp.d = d/wz*obj.w0;
            re = obj.sp.F_onlyx(HGB);
            obj.sp.offset = oldoffset;
            obj.sp.d = oldd;
        end
        
        function re = CMcor(obj, dx, theta)           
            re = (obj.F0 + 2*obj.F1*dx/obj.w0*cos(theta)) + ...
                2*obj.laser.Pdc*obj.F1*(obj.device.kmx())^2*obj.laser.k*obj.w0*dx/obj.consts.c*obj.device.response(obj.f)*sin(theta);
            re = abs(re).^2;
        end
        
        function re = CMcor_chi(obj, dx, theta)           
            re = 2*obj.laser.Pdc*obj.F1*(obj.device.kmx())^2*obj.laser.k*obj.w0*dx/obj.consts.c*obj.device.response(obj.f)*sin(theta);
            re = abs(re).^2;
        end
        
        function re = Qnoi(obj, theta)
            re = obj.laser.Pdc*obj.F1^2*obj.consts.hbar*obj.consts.c*obj.laser.k*(abs(2*obj.laser.Pdc/(obj.consts.c*obj.laser.k)*(obj.laser.k*obj.device.kmx()*obj.w0)^2*obj.device.response(obj.f).*sin(theta) - cos(theta)).^2/obj.PD.WtoA + sin(theta).^2/obj.PD.WtoA);
        end
        
        function SetInspectFreqs(obj, f)
            obj.f = f;
        end
        
        function re = thermal(obj, theta)
            re = 1/(2*pi) * 8*obj.consts.kb*obj.T*obj.device.gamma*obj.device.meff()...
                *(sin(theta)*obj.F1*obj.laser.k*obj.device.kmx()*obj.w0)^2*abs(obj.device.response(obj.f)).^2;
        end
        
        function [total, signal, noise, signal_chi, thermal, Qnoi] = CalSig_chirp(obj, dx, theta)
            %assuming laser psd is only composed of a dc and some spectrum
            %around device.fm, whose amplitude is descripbed by Pac in time
            %domain
            signal = obj.CMcor(dx, theta)*obj.laser.Pac^2/(4*obj.laser.deltaf);
            signal_chi = obj.CMcor_chi(dx, theta)*obj.laser.Pac^2/(4*obj.laser.deltaf);  
            
            thermal = obj.thermal(theta)*obj.laser.Pdc^2;
            Qnoi = obj.Qnoi(theta);
            noise = thermal + Qnoi;
            total = signal + noise;    
        end
        
        function re = WtoA(obj, psd)
            re = psd*(obj.PD.WtoA)^2;
        end
        
        function re = AtoV(obj, psd)
            re = psd*(obj.PD.Gain)^2;
        end
        
        function re = MiscConversion(obj, psd)
            re = 10*log10(20) + 10*log10(2*psd);%to single sided dBm/Hz
        end
        
        function re = ConvertPSD(obj, psd)
            re = obj.MiscConversion(obj.AtoV(obj.WtoA(psd)));
        end
        
        function [total, signal, noise, signal_chi, thermal, sideband] = GetPSDs(obj, dx, theta)
            [total, signal, noise, signal_chi, thermal, Qnoi] = obj.CalSig_chirp(dx, theta);
            total = obj.ConvertPSD(total);
            signal = obj.ConvertPSD(signal);
            noise = obj.ConvertPSD(noise);
            sideband = obj.ConvertPSD((thermal + signal_chi)*0.5*(obj.laser.Pac_sideband/obj.laser.Pdc)^2 + Qnoi);
            signal_chi = obj.ConvertPSD(signal_chi);
            thermal = obj.ConvertPSD(thermal);
        end
        
        function A = CoefA(obj, dx, theta)
            A = obj.F0 + 2*obj.F1*dx/obj.w0*cos(theta);
        end
        
        function B = CoefB(obj, dx, theta)
            B = 2*obj.F1*(obj.laser.k*(obj.device.kmx())^2*obj.w0*dx)/obj.consts.c*sin(theta);
        end
        
        function DC = DC_CM(obj, dx, theta)
            DC = obj.CoefA(dx, theta)*obj.laser.Pdc +...
                obj.CoefB(dx, theta)*(obj.laser.Pac)^2*(obj.device.response(-obj.f)+obj.device.response(obj.f))/4;
            DC = DC*obj.PD.WtoA*(obj.PD.Gain*1e3);
        end
        
        function re = twoOmega(obj, dx, theta)
            re = obj.laser.Pac^2/4*obj.CoefB(dx, theta)*obj.device.response(-obj.f);
            re = abs(re).^2;
        end
        
        function [total, signal, noise, thermal, Qnoi, ref] = AOM_omega(obj, theta, DeltaX0)
            F0002 = obj.CalFmn(obj.sp.offset, obj.sp.d, 0, 2);
            F0102 = obj.CalFmn(obj.sp.offset, obj.sp.d, 1, 2);
            re0 = (sqrt(2*pi)*DeltaX0/obj.w0*cos(theta) + sqrt(2*pi)*obj.laser.k*(obj.device.kmx())^2*DeltaX0*obj.w0*obj.laser.Pdc/obj.consts.c*obj.device.response(obj.f)*sin(theta))*obj.F1;
            re = (sqrt(2*pi)*DeltaX0/obj.w0*cos(theta) + sqrt(2*pi)*obj.laser.k*(obj.device.kmx())^2*DeltaX0*obj.w0*obj.laser.Pdc/obj.consts.c*obj.device.response(obj.f)*sin(theta))*obj.F1+...
                (2*sqrt(2)*cos(2*theta)*F0002 + sqrt(2)*sin(theta)*F0102 - sin(theta)*obj.F1)*sqrt(2*pi)*obj.laser.k*(obj.device.kmx())^2*DeltaX0^3*obj.laser.Pdc/(2*obj.w0*obj.consts.c)*(obj.device.response(-obj.f)/2+obj.device.response(obj.f));
            ref = sqrt(2*pi)*DeltaX0/obj.w0*cos(theta)*obj.F1 + obj.f - obj.f;
            re = abs(re).^2;
            re0 = abs(re0).^2;
            ref = abs(ref).^2;
            signal = re*obj.laser.Pdc^2/(4*obj.laser.deltaf);
            ref = ref*obj.laser.Pdc^2/(4*obj.laser.deltaf);
            thermal = obj.thermal(theta)*obj.laser.Pdc^2;
            Qnoi = obj.Qnoi(theta);
            noise = thermal + Qnoi;
            total = signal + noise;
        end
        
        function [total, signal, noise, thermal, Qnoi, ref] = GetPSDs_AOM(obj, theta, DeltaX0)
            [total, signal, noise, thermal, Qnoi, ref] = obj.AOM_omega(theta, DeltaX0);
            total = obj.ConvertPSD(total);
            signal = obj.ConvertPSD(signal);
            ref = obj.ConvertPSD(ref);
            noise = obj.ConvertPSD(noise);
            thermal = obj.ConvertPSD(thermal);
            Qnoi = obj.ConvertPSD(Qnoi);
        end
    end
end