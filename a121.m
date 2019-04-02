%x-direction periodic
%y-direction wall effect
%rectangular lattice
%x-column
%y-row
%lattice(y,x)=at coordinate(x,y) on the lattice
clear all
clc        
L=100;              %length in x direction
H=20;               %length in y direction
lattice=zeros(H,L); %lattice matrix, size L*H
Z=4;                %coord #
el=1;
wall=0;
bw=1;
rT=100;    
abcd=size(lattice);
N=L*H;              %number of total sites
n1=100;              %initial number of total particles
matrix=zeros(n1,2); %used to record coordinates of the occupied sites
Na_put=0;           %number of particles already put in the initial configuration
L_1=10; %length of the GC area on the left
reduced_chemical_potential_1=-2;
activity1=exp(reduced_chemical_potential_1/rT);

random_or_not=1;
[lattice,matrix]=build(n1,lattice,matrix,1,20,1,20);
E=countE(lattice,matrix,el,bw);