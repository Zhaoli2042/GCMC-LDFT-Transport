%% Introduction
%This script uses neighbor list method. At 1 subcycle, a particle is
%selected to move to its adjacent sites. A neighbor list will be generated
%at the end of the script, size of the neighbor list will be used as the
%random number to make sure that every adjacent sites has equal probability
%of being selected. 
%%
function abc=neigh_list(lattice,x,y)
up=y+1; down=y-1; left=x-1; right=x+1;
size1=size(lattice); sizeL=size1(2); sizeH=size1(1);

neighbors=[up,x;down,x;y,left;y,right];
%% if particle at boundary, eliminate corresponding neighbor from neighbor list
if (up == (sizeH+1))
    neighbors(1,:)=[-1 0];
end
if (down == 0)
    neighbors(2,:)=[-2 0];
end
if (left == 0)
    neighbors(3,:)=[-3 0];
end
if (right == (sizeL+1))
    neighbors(4,:)=[-4 0];
end
minim=min(neighbors(:,1));
while (minim <=0)
    ind=find(neighbors(:,1) == minim);
    neighbors(ind,:)=[]; %#ok<FNDSB>
    minim=min(neighbors(:,1));
end
abc=neighbors;