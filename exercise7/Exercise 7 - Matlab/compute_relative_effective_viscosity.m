function [mueffrel]=compute_relative_effective_viscosity(d,htd,invivo)

%Computes the effective viscosity 
%based on the equations from:
%Microvascular blood viscosity in vivo and the endothelial surface layer
%by Pries and Secomb, published in 2005 in Am J Physiol Heart Circ Physiol

if invivo
    %Compute physical ESL (endothelial surface layer) thickness
    %ONLY for diameters from 2.4 to 10.5
    doff =2.4;
    dcrit = 10.5;
    d50 = 100.0;
    eamp = 1.1;
    ewidth = 0.03;
    epeak = 0.6;
    wmax = 2.6;
    ehd=1.18;
    
    was=(d-doff)./(d+d50-2*doff)*wmax;
    wpeak=eamp*(d-doff)./(dcrit-doff);
    wph=was+wpeak*epeak;
    %effective esl thickness
    weff = was+wpeak.*(1+ehd*htd);
    %Compute physical vessle diameter
    dph= d-2*wph;
    deff= d-2*weff;
    
else
    dph=d;
end

nu45=220.0*exp(-1.3*dph)+3.2-2.44*exp(-0.06*dph.^0.645);
c=(0.8+exp(-0.075*dph)).*(-1.0+1./(1+1e-11*dph.^12))+1./(1+1e-11*dph.^12);
nu=1.0+(nu45-1.0).*((1.0-htd).^c-1.0)./((1-0.45).^c-1.0);

if invivo
    mueffrel = nu.*(d./deff).^4.0;
else
    mueffrel = nu;
end
end
    
    
