%x-direction periodic
%y-direction wall effect
%rectangular lattice
%x-column
%y-row
%lattice(y,x)=at coordinate(x,y) on the lattice
clear all
clc        
L=10;              %length in x direction
H=200;               %length in y direction
lattice=zeros(H,L); %lattice matrix, size L*H
Z=4;                %coord #
el=1;
wall=0;
bw=1;
rT=100;    
abcd=size(lattice);
N=L*H;              %number of total sites
n1=10;               %initial number of total particles
matrix=zeros(n1,2); %used to record coordinates of the occupied sites
Na_put=0;           %number of particles already put in the initial configuration
reduced_chemical_potential=-2;
activity=exp(reduced_chemical_potential/rT);
for x = 1:L
    for y = 1:H
        if Na_put<n1
            lattice(y,x)=1;
            Na_put=Na_put+1;
            matrix(Na_put,1)=y;
            matrix(Na_put,2)=x; %record the coordinate of the occupied sites
        end
    end
end  
%count initial energy
ij=0;
size2=size(matrix);
for m=1:1:size2(1)  %check from first to last site if it is occupied, also checked periodic boundary
    x=matrix(m,2);
    y=matrix(m,1);
    w=0;
    if y==1
        w=1;
    else 
        down=y-1;
    end
    if y==H
        w=1;
    else
        up=y+1;
    end
    if x==1
        left=L;
    else 
        left=x-1;
    end
    if x==L
        right=1;
    else
        right=x+1;
    end
    if w==0
        ij=ij+lattice(y,left)+lattice(y,right)+lattice(up,x)+lattice(down,x);
    else
        if y==1
            ij=ij+lattice(y,left)+lattice(y,right)+lattice(up,x)+w;
        elseif y==H
            ij=ij+lattice(y,left)+lattice(y,right)+lattice(down,x)+w;
        end
    end
end
E=-el/2*ij-wall*bw; %initial config energy
fprintf('Initial Energy is %d',E)

%remove a particle
Neq=10000;
Nsam=100000;                  %same as equilibrium steps
Ncyc=Neq+Nsam;                %# of cycles to run this MC
samstep=100;
Energy=zeros(1,Nsam/samstep);
insert=100;
NGC=1;                        %run GC moves every 100 cycle
for a=1:Ncyc
%MC move

pick=randi(n1);
x=matrix(pick,2);
y=matrix(pick,1);
    w=0;
    if y==1
        w=1;
    else
        down=y-1;
    end
    if y==H
        w=1;
    else
        up=y+1;
    end
    if x==1
        left=L;
    else 
        left=x-1;
    end
    if x==L
        right=1;
    else
        right=x+1;
    end
    if w==0
        ij=lattice(y,left)+lattice(y,right)+lattice(up,x)+lattice(down,x);
    else
        if y==1
            ij=lattice(y,left)+lattice(y,right)+lattice(up,x)+w;
        elseif y==H
            ij=lattice(y,left)+lattice(y,right)+lattice(down,x)+w;
        end
    end
    eold=ij*el;    
%randomly choose an adjacent particle
    ra=randi(1000);
    x_neigh=(-1)^ra+x;
    if x_neigh==(L+1)
        x_neigh=1;         %apply boundary conditions to x-direction
    end
    if x_neigh==0
        x_neigh=L;
    end
    y_neigh=(-1)^ra+y;
    if y_neigh==(H+1)
        y_neigh=H-1;       %if the particle is at y-boundaries, force it to choose inner adjacent sites
    end
    if y_neigh==0
        y_neigh=2;
    end
    if (lattice(y_neigh,x_neigh)==0)
            w=0;
            if y_neigh==1
                w=1;
            else
                down=y_neigh-1;
            end
            if y_neigh==H
                w=1;
            else
                up=y_neigh+1;
            end
            if x_neigh==1
                left=L;
            else 
                left=x_neigh-1;
            end
            if x_neigh==L
                right=1;
            else
                right=x_neigh+1;
            end
            if w==0
                ij=lattice(y,left)+lattice(y,right)+lattice(up,x)+lattice(down,x);
            else
                if y==1
                    ij=lattice(y,left)+lattice(y,right)+lattice(up,x)+w;
                elseif y==H
                    ij=lattice(y,left)+lattice(y,right)+lattice(down,x)+w;
                end
            end
            enew=ij*el;


            de=enew-eold;
            if de<0 
                f=1;
            else
                g=rand();
            if exp(-de/rT)>g
                f=1;
            else 
                f=0;
            end
            end
    elseif (lattice(y_neigh,x_neigh)==1)
            f=0;    
    end
    if f==1
        lattice(y,x)=0;
        lattice(y_neigh,x_neigh)=1;
        E=E+de;
        matrix(pick,:)=[y_neigh,x_neigh];
    end


%GC MOVE

    if (a>Neq) && (mod(a,NGC)==0)
        if (rand<0.5)
    %adding particle
            x=randi(L);
            y=randi(H);
            if lattice(y,x)==0
                w=0;
                if y==1
                    w=1;
                else 
                    down=y-1;
                end
                if y==H
                    w=1;
                else
                    up=y+1;
                end
                if x==1
                    left=L;
                else 
                    left=x-1;
                end
                if x==L
                    right=1;
                else
                    right=x+1;
                end
                if w==0
                    ij=lattice(y,left)+lattice(y,right)+lattice(up,x)+lattice(down,x);
                else
                    if y==1
                        ij=lattice(y,left)+lattice(y,right)+lattice(up,x)+w;
                    elseif y==H
                        ij=lattice(y,left)+lattice(y,right)+lattice(down,x)+w;
                    end
                end
                delta_e=-ij*el;
                if (rand<exp(-delta_e/rT)*activity*N/(n1+1))
                    E=E+delta_e;
                    n1=n1+1;
                    new_row=[y,x];
                    lattice(y,x)=1;
                    matrix=[matrix;new_row];
                end
            end
        else
            if (n1>0)      
            size1=size(matrix);
            particle_id=randi(size1(1));
            chosey=matrix(particle_id,1);
            chosex=matrix(particle_id,2);
                w=0;
                if chosey==1
                    down=H;
                else 
                    down=chosey-1;
                end
                if chosey==H
                    up=1;
                else
                    up=chosey+1;
                end
                if chosex==1
                    left=L;
                else 
                    left=chosex-1;
                end
                if chosex==L
                    right=1;
                else
                    right=chosex+1;
                end
                if w==0
                    ij=lattice(chosey,left)+lattice(chosey,right)+lattice(up,chosex)+lattice(down,chosex);
                else
                    if chosey==1
                        ij=lattice(chosey,left)+lattice(chosey,right)+lattice(up,chosex)+w;
                    elseif chosey==H
                        ij=lattice(chosey,left)+lattice(chosey,right)+lattice(down,chosex)+w;
                    end
                end
                delta_e=ij*el; % well, need to ask why they fliped the sign
                if (rand<exp(delta_e/rT)*n1/(N*activity))
                    E=E+delta_e;
                    lattice(chosey,chosex)=0;
                    matrix(particle_id,:)=[];
                    n1=n1-1;
                end
            end
        end
    end 
end
density=n1/N;
disp(density)