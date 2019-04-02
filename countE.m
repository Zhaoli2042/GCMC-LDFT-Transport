%To avoid double counting, we multiply 0.5 before lattice occupancy
%So this subroutine gives positive energy(no double counting)

function abc=countE(lattice,matrix,el,bw)
ij=0;
size2=size(matrix);
size1=size(lattice); sizeL=size1(2); sizeH=size1(1);
for m=1:1:size2(1)  %check from first to last site if it is occupied, also checked periodic boundary
    x=matrix(m,2);
    y=matrix(m,1);
    w=0;
    if y==1
        w=1;
    else 
        down=y-1;
    end
    if y==sizeH
        w=1;
    else
        up=y+1;
    end
    if x==1
        left=sizeL;
    else 
        left=x-1;
    end
    if x==sizeL
        right=1;
    else
        right=x+1;
    end
    if w==0
        ij=ij+(lattice(y,left)+lattice(y,right)+lattice(up,x)+lattice(down,x))*0.5*el;
    else
        if y==1
            ij=ij+0.5*el*(lattice(y,left)+lattice(y,right)+lattice(up,x))+w*bw;
        elseif y==sizeH
            ij=ij+0.5*el*(lattice(y,left)+lattice(y,right)+lattice(down,x))+w*bw;
        end
    end
end
abc=-ij;