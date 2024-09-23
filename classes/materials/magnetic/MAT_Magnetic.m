classdef MAT_Magnetic < MAT
    properties
        IsLinearMAG = true;
        Mu0 = 4*pi*1e-7
        Sigma
        Mur
        B
        H
        HBspline;
        HBsplineDer;
        HBsplineDer2;
        HBsplineInt;
    end

    methods (Access = public)
        function mu = getMuLinear(obj)
            mu = obj.Mur*obj.Mu0;
        end

        function nu = getNuLinear(obj)
            nu = 1./obj.getMuLinear();
        end

        function setMurLinear(obj, mur)
            obj.Mur = mur;
        end

        function setSigma(obj, sigma)
            obj.Sigma = sigma;
        end

        function sig = getSigma(obj)
            sig = obj.Sigma;
        end
        
        function setLinear(obj, trueFalse)
            assert(trueFalse==0 || trueFalse == 1);
            obj.IsLinearMAG = trueFalse;
        end
        
        function fitHBspline(obj)
            obj.HBspline = pchip(obj.B, obj.H);
            obj.HBsplineDer = fnder(obj.HBspline);
            obj.HBsplineDer2 = fnder(obj.HBsplineDer);
            obj.HBsplineInt = fnint(obj.HBspline);
        end

        function mu = getMuNonlinear(obj, B)
            B(B<=1e-8) = 1e-8; % avoid zero division
            mu = B./ppval(obj.HBspline, B);
            % correct B values higher than measurement (assume behavior like vacuum)
            mu(B>obj.B(end)) = B(B>obj.B(end))./((B(B>obj.B(end))-obj.B(end))./obj.Mu0 + obj.H(end));
        end
        function nu = getNuNonlinear(obj, B)
            nu = 1./obj.getMuNonlinear(B);
        end

        function nuPrime = getNuPrimeNonlinear(obj, B)
            B(B<=1e-8) = 1e-8; % avoid zero division
            nuPrime = (ppval(obj.HBsplineDer, B).*B - ppval(obj.HBspline, B))./(B.^2);
            % correct B values higher than measurement (assume behavior like vacuum)
            B2correct = B(B>obj.B(end));
            % nuPrime(B>obj.B(end)) = (B2correct./obj.Mu0 - obj.getNuNonlinear(B2correct).*B2correct)./B2correct.^2; 
            nuPrime(B>obj.B(end)) = (obj.B(end)-obj.Mu0*obj.H(end))/obj.Mu0./B2correct.^2; 
        end

        function d2nudB2 = getd2nudB2(obj, B)
            B(B<=1e-8) = 1e-8; % avoid zero division
            d2nudB2 = (ppval(obj.HBsplineDer2, B).*B.^2 - 2* ppval(obj.HBsplineDer, B).*B  + 2*ppval(obj.HBspline, B))./(B.^3);
            B2correct = B(B>obj.B(end));
            d2nudB2(B>obj.B(end)) = -2*(obj.B(end)-obj.Mu0*obj.H(end))/obj.Mu0./B2correct.^3; 
        
        end

        function C = getCoenergy(obj, B)
            if obj.IsLinearMAG
                C = 1/(2*obj.getMuLinear())*B.^2;
            else
                C = B.^2./obj.getMuNonlinear(B) - obj.getEnergy(B);
            end
        end
        
        function E = getEnergy(obj, B)
            if obj.IsLinearMAG
                E = 1/(2*obj.getMuLinear())*B.^2;
            else
                E = ppval(obj.HBsplineInt, B);
            end
        end 

        function plotCharacteristics(obj)
            f = figure("Name", "Plots");
            % B-H curve
            subplot(2,2,1)
            obj.plotBH();
            % B-mu-curve
            subplot(2,2,2)
            obj.plotBMu();
            % B-nu-curve
            subplot(2,2,3)
            obj.plotBNu();
            % B-nuprime-curve
            subplot(2,2,4)
            obj.plotBNuPrime();
        end

        function plotBH(obj)
            b = 0:0.001:max(obj.B)*1.1;
            h = b./obj.getMuNonlinear(b);
            % B-H curve
            plot(h, b, "Color", "black", LineWidth=1.5)
            hold on
            scatter(obj.H, obj.B, "filled", "MarkerEdgeColor","black", "Marker","x","LineWidth",2)
            grid on
            xlabel("H (A/m)")
            ylabel("B (T)")
            title("B-H-curve")
        end

        function plotBMu(obj)
            b = 0:0.001:max(obj.B)*1.1;
            mur = obj.getMuNonlinear(b)/obj.Mu0;
            plot(b, mur, "Color", "black", LineWidth=1.5)
            hold on
            scatter(obj.B, obj.B./obj.H/obj.Mu0, "filled", "MarkerEdgeColor","black", "Marker","x","LineWidth",2)
            grid on
            xlabel("B (T)")
            ylabel("Mu_r")
            title("B-mu-curve")
        end

        function plotBNu(obj)
            b = 0:0.001:max(obj.B)*1.1;
            nu = obj.getNuNonlinear(b);
            plot(b, nu, "Color", "black", LineWidth=1.5)
            hold on
            scatter(obj.B, obj.H./obj.B, "filled", "MarkerEdgeColor","black", "Marker","x","LineWidth",2)
            grid on
            xlabel("B (T)")
            ylabel("Nu ")
            title("B-nu-curve")
        end

        function plotBNuPrime(obj)
            b = 0:0.001:max(obj.B)*1.1;
            nuP = obj.getNuPrimeNonlinear(b);
            plot(b, nuP, "Color", "black", LineWidth=1.5)
            grid on
            xlabel("B (T)")
            ylabel("dNu/dB ")
            title("B-nu'-curve")
        end

        function plotBNuPrimeB2(obj)
            % TBD NOT YET FINISHED
            b = 0:0.001:max(obj.B)*1.1;
            nuP = obj.getNuPrimeB2Nonlinear(b);
            plot(b, nuP, "Color", "black", LineWidth=1.5)
            grid on
            xlabel("B (T)")
            ylabel("dNu/dB^2 ")
            title("B^2-nu'-curve")
        end
    end
end