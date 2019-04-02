%Initial Configuration builder
%can choose regular build or random build
%I will first write a regular one
%random section didn't work if n1=1;
function [c,d]=inic(n1,lattice,matrix,random_or_not)
Na_put=0;
size1=size(lattice); sizeL=size1(2);sizeH=size1(1);
if (random_or_not==0)
    for x=1:sizeL
        for y=1:sizeH
            if Na_put<n1
                lattice(y,x)=1;
                Na_put=Na_put+1;
                matrix(Na_put,1)=y;
                matrix(Na_put,2)=x;
            end
        end
    end
else
    lattice(1,1)=1;
    matrix(1,1)=1;
    matrix(1,2)=1;
    chose_x=1;
    chose_y=1;
    for l=1:1:n1-1
        while (lattice(chose_y,chose_x)==1)
            chose_x=randi(sizeL);
            chose_y=randi(sizeH);
        end
        lattice(chose_y,chose_x)=1;
        matrix(l+1,1)=chose_y;
        matrix(l+1,2)=chose_x;
    end
end
c=lattice;
d=matrix;