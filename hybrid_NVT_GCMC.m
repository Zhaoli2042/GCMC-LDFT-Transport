%x-direction wall at 2 GCMC region, bw=0
%y-direction wall effect
%rectangular lattice any size
%x-column
%y-row
%lattice(y,x)=at coordinate(x,y) on the lattice
%% Test cases
% -2&-2:    ~0.5
% -5&-5:    ~0.38
%-10&-10:   ~
%% Initial Settings
clear all
clc        
L=160;              %length in x direction
H=20;               %length in y direction
lattice=zeros(H,L); %lattice matrix, size L*H
Z=4;                %coord #
el=1;
bw=1;
rT=2;    
N=L*H;              %number of total sites
n1=50;              %initial number of total particles
matrix=zeros(n1,2); %used to record coordinates of the occupied sites
Na_put=0;           %number of particles already put in the initial configuration
L_1=40; %length of the GC area on the left
reduced_chemical_potential_1=-1;
activity1=exp(reduced_chemical_potential_1/rT);

L_2=40; %length of GC area on the right
L_nvt=L-L_1-L_2;
reduced_chemical_potential_2=-10;
activity2=exp(reduced_chemical_potential_2/rT);


%% WIDOM INSERTION INPUTS

        LHS_MUid = 0; RHS_MUid = 0;
        LHS_MUex = 0; RHS_MUex = 0;
        LHS_MU = 0; RHS_MU = 0;
%% Building Initial Configuration and Count Initial Energy

[lattice,matrix]=build(25,lattice,matrix,L-L_2+1,L,1,H); 
[lattice,matrix]=build(25,lattice,matrix,1,L_1,1,H); 
%count initial energy
E=countE(lattice,matrix,el,bw);
fprintf('Initial Energy is %d',E)
%% Settings for MC cycles
Neq=300000;
Nsam=100000;                  %same as equilibrium steps
Ncyc=Neq+Nsam;                %# of cycles to run this MC
samstep=100;
% Energy=zeros(1,Nsam/samstep);
NGC=1;                        %run GC moves every NGC cycle(s)
center=zeros(Nsam,L);         %record the center densities
mid=H/2;                      %calculate the middle point
Ninsert=10;                   %Do Widom test every 10 MC cycles
%% Cycle begins
for a=1:Ncyc
%% MC move
    if (n1>0 && (a >= Neq))
        for substeps=1:1:5
            pick=randi(n1);
            x=matrix(pick,2);
            y=matrix(pick,1);
            eold=single_countE(x,y,lattice,el,bw, L_1, L_2);    
            %randomly choose an adjacent site
            neighbors=neigh_list(lattice,x,y);
            aa=size(neighbors);
            ab=aa(1);
            ac=randi(ab);
            y_neigh=neighbors(ac,1); x_neigh=neighbors(ac,2);
            if (lattice(y_neigh,x_neigh)==0)
                enew=single_countE(x_neigh,y_neigh,lattice,el,bw, L_1, L_2);
                    de=enew-eold;
                    if (rand < exp(-de/rT))
                        lattice(y,x)=0;
                        lattice(y_neigh,x_neigh)=1;
                        E=E+de;
                        matrix(pick,:)=[y_neigh,x_neigh];
                    end
            end
        end
    end
    %% GC MOVE

    %if (mod(a,NGC)==0)
        for abcdefg=1:1:40
            LorR=randi(2);            %left or right
            if (LorR==1)
                activity=activity1;
                pic_ran=[1,L_1];      %pick a random site from left
                N=L_1*H;
                n=sum(sum(lattice(1:H,1:L_1)));%count the number of particles 
                %in the left region
            else
                activity=activity2;
                pic_ran=[L-L_2+1,L];    %pick a random site from right
                N=L_2*H;
                n=sum(sum(lattice(1:H,(L-L_2+1):L)));
                %count in the right region
            end
            if (rand<=0.5)
        %% adding particle
                x=randi(pic_ran);
                y=randi(H);
                if lattice(y,x)==0
                    delta_e=single_countE(x,y,lattice,el,bw, L_1, L_2)-0;
                    if (rand<exp(-delta_e/rT)*activity*N/(n+1))
                        E=E+delta_e;
                        n1=n1+1;
                        new_row=[y,x];
                        lattice(y,x)=1;
                        matrix=[matrix;new_row]; %#ok<AGROW>
                    end
                end
            else
            %% Removing particle
            % particle in GC regions
             ngc = sum(sum(lattice(1:H,1:L_1))) + sum(sum(lattice(1:H,L-L_2+1:L)));
                if (n1>0 && ngc > 0)      
                size1=size(matrix);
                chosex=L_1+1;
                    while ((chosex>L_1) && (chosex<(L-L_2+1))) %random until you get 
                        %a value not in NVT region
                        particle_id=randi(size1(1));
                        chosex=matrix(particle_id,2);
                    end
                    chosey=matrix(particle_id,1);
                    w=0;
                    delta_e=0-single_countE(chosex,chosey,lattice,el,bw, L_1, L_2); 
                    % well, need to ask why they fliped the sign
                    if (rand<exp(-delta_e/rT)*n/(N*activity))
                        E=E+delta_e;
                        lattice(chosey,chosex)=0;
                        matrix(particle_id,:)=[];
                        n1=n1-1;
                    end
                end
            end
        end
    %end

    %% record averages of occurance in matrix called "relative_density "
    if a>Neq
        numb=a-Neq;
        if a==(Neq+1)
            re_rou=lattice;
        else
            re_rou=(re_rou*(numb-1)+lattice)/numb;
        end
    end
    
    %% WIDOM PARTS
     if ((a>Neq) && (mod(a,Ninsert)==0))
         rh = 0; lh = 0;

        for times=1:1:200       %run Widom 200 times every test
            rh = Widom(lattice,L-L_2+1,L,el,bw,rT, L_1, L_2);
            lh = Widom(lattice,1,L_1,el,bw,rT, L_1, L_2);
            if times==1
                rhs= rh;
                lhs = lh;
            else
                rhs=[rhs rh]; %#ok<AGROW>
                lhs = [lhs lh]; %#ok<AGROW>
            end
        end
        RHSS=mean(rhs);
        LHSS=mean(lhs);
           RHS = 0;
           LHS =0;
        if RHS == 0
            RHS=RHSS;
            LHS=LHSS;
        else
            RHS=[RHS RHSS]; %#ok<AGROW>
            LHS=[LHS LHSS]; %#ok<AGROW>

        end
        if LHS_MUid ==0
            LHS_MUid = rT * log(sum(sum(lattice(1:H,1:L_1)))/(L_1*H));
            RHS_MUid = rT * log(sum(sum(lattice(1:H,L-L_2+1:L)))/(L_2*H));
            LHS_MUex = -log(LHS)*rT;
            RHS_MUex = -log(RHS)*rT;
            LHS_MU = LHS_MUid + LHS_MUex;
            RHS_MU = RHS_MUid + RHS_MUex;
        else
            LHS_MUid = LHS_MUid + rT * log(sum(sum(lattice(1:H,1:L_1)))/(L_1*H));
            RHS_MUid = RHS_MUid + rT * log(sum(sum(lattice(1:H,L-L_2+1:L)))/(L_2*H));
            LHS_MUex = LHS_MUex -log(LHS)*rT;
            RHS_MUex = RHS_MUex -log(RHS)*rT;
            LHS_MU = LHS_MU + LHS_MUid + LHS_MUex;
            RHS_MU = RHS_MU + RHS_MUid + RHS_MUex;  
            
%             if (LHS_MU == -inf || RHS_MU == -inf)
%                 break;
%             end
        end
     end
     
% WRITE TO TEXT FILES FOR VISUALIZATION%
    if (a > Neq && (mod(a, 100) == 0))
    write_to_file(matrix, L_1, L_2, L, H, a);
    end
end
density=n1/(L*H);
disp(density)
E=0-countE(lattice,matrix,el,bw);
fprintf('Final Energy is %d',E)
fprintf('\nDensity on the Left is %d',sum(sum(lattice(1:H,1:L_1)))/(L_1*H))
fprintf('\nDensity on the Right is %d',sum(sum(lattice(1:H,L-L_2+1:L)))/(L_2*H))
xx=linspace(1,L,L);
xx_nvt = linspace(L_1+1, L-L_2, (L-L_1-L_2));
yy=linspace(1,H,H);
% surf(xx,yy,re_rou)
plot(yy,re_rou(:,1))
%plot density profile(rou_center, density profile in the center)
%start from a closer initial point

%rou_wall(at two walls)
% surf(xx_nvt,yy,re_rou(:,L_1+1: L-L_2))
% xlabel('X-axis distance')
% ylabel('Y-axis distance')
% zlabel('reduced density')
surf(xx,yy,re_rou)
xlabel('X-axis distance')
ylabel('Y-axis distance')
zlabel('reduced density')

%% CALCULATE AVERAGE WIDOMs
LHS_MUid = LHS_MUid / (Nsam/Ninsert);
RHS_MUid = RHS_MUid / (Nsam/Ninsert);
LHS_MUex = LHS_MUex / (Nsam/Ninsert);
RHS_MUex = RHS_MUex / (Nsam/Ninsert);
LHS_MU = LHS_MUid + LHS_MUex;
RHS_MU = RHS_MUid + RHS_MUex;
fprintf("lefthand chem potent is %.2f\n", LHS_MU);
fprintf("righthand chem potent is %.2f\n", RHS_MU);

% LHS_MU = LHS_MU / (Nsam/Ninsert);
% RHS_MU = RHS_MU / (Nsam/Ninsert);


%     %% Widom Particle Insertion Method for NVT region
% %     if ((a>Neq) && (mod(a,Ninsert)==0))
%         for times=1:1:200       %run Widom 5 times every test
%             rh = Widom(lattice,L-L_2+1,L,el,bw,rT, L_1, L_2);
%             lh = Widom(lattice,1,L_1,el,bw,rT, L_1, L_2);
%             if (rh ~= -inf && lh ~= -inf)
%                 if times==1 
%                     rhs= rh;
%                     lhs = lh;
%                 else
%                     rhs=[rhs rh]; %#ok<AGROW>
%                     lhs = [lhs lh]; %#ok<AGROW>
%                 end
%             end
%         end
%         RHSS=mean(rhs);
%         LHSS=mean(lhs);
%            RHS = 0;
%            LHS =0;
%         if RHS == 0
%             RHS=RHSS;
%             LHS=LHSS;
%         else
%             RHS=[RHS RHSS]; %#ok<AGROW>
%             LHS=[LHS LHSS]; %#ok<AGROW>
% 
%         end
%         LHS_MUid = rT * log(sum(sum(lattice(1:H,1:L_1)))/(L_1*H));
%         RHS_MUid = rT * log(sum(sum(lattice(1:H,L-L_2+1:L)))/(L_2*H));
%         LHS_MUex = -log(LHS)*rT;
%         RHS_MUex = -log(RHS)*rT;
%         LHS_MU = LHS_MUid + LHS_MUex;
%         RHS_MU = RHS_MUid + RHS_MUex;
%         fprintf("lefthand chem potential is %.2f\n", LHS_MUid + LHS_MUex);
%         fprintf("fixed LHS chemical potential is %.2f\n", reduced_chemical_potential_1);
%         fprintf("righthand chem potential is %.2f\n", RHS_MUid + RHS_MUex);
%         fprintf("fixed RHS chemical potential is %.2f\n", reduced_chemical_potential_2);

% %     end