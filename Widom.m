% Widom Particle Insertion Method-for hybrid NVT/muVT
% L means length(x-direction), H means height(y-direction)
% L_1 GC region on the left, L_nvt: NVT region in the middle
function rhs=Widom(lattice,L_1,Lx,el,bw,rT, left, right)
size1=size(lattice); sizeH=size1(1);
pickx=randi([L_1,Lx]); picky=randi(sizeH);
if (lattice(picky,pickx)==1)
      abc = 0;                %rhs=exp(-u/rT)
elseif (lattice(picky,pickx)==0)
    u=single_countE(pickx,picky,lattice,el,bw, left, right);
    abc=exp(-u/rT);
end

rhs=abc;