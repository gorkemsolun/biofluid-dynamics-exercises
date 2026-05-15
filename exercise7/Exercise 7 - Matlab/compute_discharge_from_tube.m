function [htd]=compute_discharge_from_tube(ht,d,invivo)
%Computes the discharge hematocrit from the tube hematocrit.
%based on the equations from:
%Microvascular blood viscosity in vivo and the endothelial surface layer
%by Pries and Secomb, published in 2005 in Am J Physiol Heart Circ Physiol

%The ratio of discharge to tube hematocrit, is used to calculate the
%Fahraeus effect.
%Furthermore it is needed as input to calculate the effective viscosity.

%It can be chosen between the invivo (invivo = 1) and 
%the invitro (invivo = 0) formulation

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
    
    was=(d-doff)./(d+d50-2*doff)*wmax;
    wpeak=eamp*(d-doff)./(dcrit-doff);
    wph=was+wpeak*epeak;
    
    %Compute physical vessel diameter
    dph= d-2*wph;
    
    %Compute discharge hematocrit
    x=1.0+1.7*exp(-0.415*dph)-0.6*exp(-0.011*dph);
    ht_vitro=ht.*(d./dph).^2;
    htd = 0.5*(x-sqrt(-4*ht_vitro.*x+x.^2+4*ht_vitro))./(x-1);
    
else
    %Compute discharge hematocrit
    x=1.0+1.7*exp(-0.415*d)-0.6*exp(-0.011*d);
    htd = 0.5*(x-sqrt(-4*ht.*x+x.^2+4*ht))./(x-1);
end

%Recheck that no value is larger than 1.0:
for i=1:size(htd,2)
    if htd(i) >= 1.0
        htd(i) = 1.0;
    end
end
end