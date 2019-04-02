%count single point energy,all walls
function e=single_countE(x,y,lattice,el,bw,L_1, L_2)
w=0;
wlr=0;      %non-interacting wall on the left and right
size1=size(lattice); sizeL=size1(2);sizeH=size1(1);
ij=0;
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
    wlr=1;
else 
    left=x-1;
end
if x==sizeL
    wlr=1;
else
    right=x+1;
end
if w==0
    if wlr == 0
        ij=ij+(lattice(y,left)+lattice(y,right)+lattice(up,x)+lattice(down,x))*el;
    else
        if x == 1
            ij=ij+(lattice(y,right)+lattice(up,x)+lattice(down,x))*el;
        elseif x == sizeL
            ij=ij+(lattice(y,left)+lattice(up,x)+lattice(down,x))*el;
        end
    end
else
    if wlr == 0
        if (x > L_1 && x < (sizeL-L_2+1))
            % only wall effect in NVT REGION
            if (y==1)
                ij=ij+el*(lattice(y,left)+lattice(y,right)+lattice(up,x))+w*bw;
            elseif (y==sizeH)
                ij=ij+el*(lattice(y,left)+lattice(y,right)+lattice(down,x))+w*bw;
            end
        else
            % no wall effect in GC regions, periodic
            if (y==1)
                ij=ij+el*(lattice(y,left)+lattice(y,right)+lattice(up,x) + lattice(sizeH,x));
            elseif (y==sizeH)
                ij=ij+el*(lattice(y,left)+lattice(y,right)+lattice(down,x) + lattice(1,x));
            end
        end
    else
        % ignore y direction wall effect in GC REGION--periodic
        if x == 1
            if y==1
                ij=ij+el*(lattice(y,right)+lattice(up,x) + lattice(sizeH,x));
            elseif y==sizeH
                ij=ij+el*(lattice(y,right)+lattice(down,x) + lattice(1,x));
            end
        elseif x == sizeH
            if y==1
                ij=ij+el*(lattice(y,left)+lattice(up,x) + lattice(sizeH,x));
            elseif y==sizeH
                ij=ij+el*(lattice(y,left)+lattice(down,x) + lattice(1,x));
            end
        end
    end
            
end
e=-0.5*ij;