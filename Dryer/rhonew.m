      function [Frho] =  rhonew(Fcomp,Fctype,FtempK)

%  RHO       Density of multi-component streams
%
%            Synopsis: MATLAB function to determine the density of a
%            multi-component stream (kg/m^3).
%
%            Call: rho(Fcomp,Fctype,FtempK)
%
%            Passed Variables:
%             Fcomp - component mass fractions (w/w)
%             Fctype - component types: 1- water
%                                       2- protein
%                                       3- carbohydrate
%                                       4- fiber
%                                       5- fat
%                                       6- ash
%                                       7- air
%             FtempK - temperature (K)
%
%            Returned Variables:
%             Frho - density (kg/m^3)
%
%            Example:
%             The density of a 12% (w/w) carbohydrate and 5% (w/w) protein
%             solution at 22oC may be found by:
%             Fcomp = [0.12 0.05 0.83]; 
%             Fctype = [2 3 1]; 
%             T = 22+tref; 
%             p = rho(Fcomp, Fctype, T)
%             
%             p =  1.0458e+03
%             
%            Reference:
%             Choi, Y & Okos, M.R. 1986.  Effects of Temperature and
%               Composition on the thermal properties of foods. In, "Food
%               Engineering Applications. Vol. 1," Elsevier Applied Science
%               Publishers, N.Y.

%  ************************************************************************

%     Reference temperture (K)  
      Tref = 273.15;  

%     Convert T (K) to T (C)
      TC = FtempK - Tref;

      for I = 1:length(Fcomp)

        if (Fctype(I) == 1)
 
%         Density of water component (kg/m^3)        
          rhoh20 = 997.18 + 0.0031439*TC - (0.0037574)*(TC^2); 

%         Density of ice component (kg/m^3)
        %  rhoice = 9.1689e2 - 1.3071e-1*TC;

%         Unfrozen and frozen water fractions (w/w)
         % ice = unfroz(Fcomp,Fctype,FtempK);
          %Fh2ouf = ice(1);  % unfrozen
          %Fh2of = ice(2);   % frozen

%         Weighted specific volume of water component (m^3/kg)        
       %   Fcsv(I) = (Fh2ouf/rhoh20)+(Fh2of/rhoice);  
          Fcsv(I) = Fcomp(I)/rhoh20;
        elseif (Fctype(I) == 2)

%         Density of protein component (kg/m^3)        
          rhopro = 1329.9 -0.5184*TC;  
%         Weighted specific volume of protein component (m^3/kg)        
          Fcsv(I) = Fcomp(I)/rhopro;  

        elseif (Fctype(I) == 3)

%         Density of carbohydrate component (kg/m^3)        
          rhocarb = 1599.1 - 0.31046*TC;    
%         Weighted specific volume of carbohydrate component (m^3/kg)
          Fcsv(I) = Fcomp(I)/rhocarb;   

        elseif (Fctype(I) == 4)
 
%         Density of fiber component (kg/m^3)
          rhofiber = 1311.5 - 0.36589*TC;
%         Weighted specific volume of fiber component (m^3/kg) 
          Fcsv(I) = Fcomp(I)/rhofiber;   

        elseif (Fctype(I) == 5)

%         Density of fat component (kg/m^3)        
          rhofat = 925.59 - 0.41757*TC;    
%         Weighted specific volume of fat component (m^3/kg) 
          Fcsv(I) = Fcomp(I)/rhofat;   

        elseif (Fctype(I) == 6)

%         Density of ash component (kg/m^3)        
          rhoash = 2423.8 - 0.28063*TC;  
%         Weighted specific volume of ash component (m^3/kg)  
          Fcsv(I) = Fcomp(I)/rhoash;   

        end  
      end 
 
%     Density of stream weighted by compositional breakdown (kg/m^3)
      Frho = 1.0/sum(Fcsv);
      end
      
     
  
