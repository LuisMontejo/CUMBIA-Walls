function [ec,fc] = manderconfWALL(Ec,AstBZ,Dh,clb,sh,fpc,fyh,eco,esm,...
    espall,HBZ,B,ncx,ncy,wi,dels)
        

% confined concrete:

sp  = sh - Dh;
Ash = 0.25*pi*(Dh^2);


bc   = B - 2*clb + Dh;
dc   = HBZ - clb + Dh;
Asx  = ncx*Ash;
Asy  = ncy*Ash;
Ac   = bc*dc;
rocc = AstBZ/Ac;
rox  = Asx/(sh*dc);
roy  = Asy/(sh*bc);
ros  = rox + roy;
ke   = ((1 - sum(wi.^2)/(6*bc*dc)) * (1 - sp/(2*bc)) ...
            * (1-sp/(2*dc))) / (1 - rocc);
ro   = 0.5*ros;
fpl  = ke*ro*fyh;
   

fpcc = (-1.254 + 2.254*sqrt(1 + 7.94*fpl/fpc) - 2*fpl/fpc)*fpc;
ecc  = eco*(1 + 5*(fpcc/fpc-1));
Esec = fpcc/ecc;
r    = Ec/(Ec-Esec);
ecu  = max((0.004 + 0.6*ros*fyh*esm/fpcc), espall);

ec = [0:dels:ecu];
x  = (1/ecc)*ec;
fc = fpcc*x*r./(r-1+x.^r);
