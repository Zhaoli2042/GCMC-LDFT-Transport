%Initial Configuration builder--version 2
%random build
%designed for putting particles randomly in 2 separated regions, can also
%be used for 1 region
%random section didn't work if n1=1;
%I think it is working...
function [c,d]=build(n1,lattice,matrix,x0,xn,y0,yn)
%here x0 and xn means the size of creating area in x-direction
%y0 and yn mean size in y-direction
size1=size(lattice); sizeL=size1(2);sizeH=size1(1);
%% Put a particle at corner(y0,x0)
    lattice(y0,x0)=1;     
    n11=sum(sum(lattice));
    matrix(n11,1)=y0;
    matrix(n11,2)=x0;
    chose_x=x0;
    chose_y=y0;
%% Randomly put particles
for l=n11+1:1:n11+n1-1
    while (lattice(chose_y,chose_x)==1)
        %randomly choose site until you find an empty one
        chose_x=randi([x0,xn]);
        chose_y=randi([y0,yn]);
    end
    lattice(chose_y,chose_x)=1;
    matrix(l,1)=chose_y;
    matrix(l,2)=chose_x;
end
c=lattice;
d=matrix;