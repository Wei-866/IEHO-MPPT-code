% The code has been taken from the study:

%_________________________________________________________________________
%  Elk herd optimization (EHO) source codes
%
%  Author and programmer:: Mohammed Azmi Al‑Betar,Mohammed A.Awadallah,Malik Shehadeh Braik,Sharif Makhadmeh,Iyad Abu Doush

%  Paper: Al-Betar, M. A., Awadallah, M. A., Braik, M. S., Makhadmeh, S., & Doush, I. A. (2024). Elk herd optimizer: a novel nature-inspired metaheuristic algorithm. Artificial Intelligence Review, 57(3), 48.
%_________________________________________________________________________

%  Improved Elk Herd Optimization(IEHO)                                                                   
%                                                                                                     
%  The Intel(R) Core(TM) i9 processor with the primary frequency of 2.80GHz, 32GB memory, and the operating system of 64-bit windows 11 using matlab2024a. 
%
%  Author and programmer: Gang Zheng,Wenchang Wei,Heming Jia,Yiqi Liu,Jiankai Lin                                                                      
%         e-Mail: ericzg@nefu.edu.cn; weiwenchang@nefu.edu.cn; jiaheming@fjsmu.edu.cn                                                                                                                                                                                                                                                

function [duty,iterations,Pos] = IEHO_MPPT(vpv,ipv)
persistent p Pold u dc dbest counter iteration iter_max num num_of dc_old p_old t zero_indices num_elements Mmax Mw Memory Memoryf;

if isempty(num)
    num=6; % population quantity
    MalesRate=0.2;
    num_of = round(num*MalesRate);
end
if isempty(p)
    p=zeros(1,num);
    p_old=zeros(1,num);
    Pold=1;
    counter=1;
    u=1;
    iteration=0;
    iter_max=15;
    t=1;
    zero_indices=[];
    zero_indices=[1 2 3 4 5 6];
    num_elements=num;
    Mw=0;
    Mmax=num*2;
    Memory=zeros(1,num*2);
    Memoryf=zeros(1,num*2);
end
if isempty(dc)
    %% Initialize the population
    dc=zeros(1,num);
    dc = linspace(0, 1, 8);  
    dc = dc(2:end-1);
    dbest=0;
    dc_old=dc;
end


Pos=dc';
iterations=iteration;
if iterations<=iter_max
    if(counter>=1 && counter<=1500)
        duty=dc(u);
        counter=counter+1;
        return;
    end
%% Memory archive construction
    if(u>=1 && u<=num)
        p(u) = vpv*ipv;
        [Mw,Memory,Memoryf]=saveMemory(Mmax,Mw,Memory,Memoryf,dc(u),p(u));
    end

    if num_elements>1 && t~=num_elements
        t=t+1;
        u=zero_indices(t);
        duty=dc(u);
        counter=1;
        return;
    end
    u=1;
    counter=1;
    if iteration>0
        NewFitness=cat(2,p,p_old);
        NewPopulation=cat(2,dc,dc_old);
        [sorted_NewFitness,sorted_indexes]=sort(NewFitness,'descend');
        for i=1:num
            dc(i)=NewPopulation(sorted_indexes(i));
            p(i)= sorted_NewFitness(i);
        end
    end
    iteration=iteration+1;
    p_old=p;
    dc_old=dc;
    [m,j]=sort(p,'descend');
    dbest=dc(j(1));
    [dc,K,p]=IEHO_UpdateDuty(dc',num,p,num_of,iteration,iter_max,Memory,Memoryf);
    zero_indices = find(K==0);
    if isempty(zero_indices)
        u = 1;
        duty=dc(u);
        Pold=vpv*ipv;
        return;
    end
    num_elements = length(zero_indices);
    t=1;
    u=zero_indices(t);
    duty=dc(u);
    Pold=vpv*ipv;
    return;
else
    duty=dbest;
    P=vpv*ipv;
    if abs(P-Pold)>0.05*Pold  
        iteration=0;
        dc = linspace(0, 1, 8);   
        dc = dc(2:end-1);
        Mw=0;
        Memory=zeros(1,num*2);
        Memoryf=zeros(1,num*2);
    end
    Pold=P;
    return;
end
end

%% IEHO-MPPT
function [D,K,p]=IEHO_UpdateDuty(dc,num,p,num_of,iteration,iter_max,Memory,Memoryf)


K=zeros(1,num);
D=zeros(1,num);
TransposeFitness=zeros(1,num_of);

%Sort the ELK positions
    [sorted_ELKS_fitness,sorted_indexes]=sort(p,'descend');

    % Make a copy of population
    NewElkHerd = dc;
    NewElkHerdFitness = p;
    
    dbest=dc(sorted_indexes(1),:); % ELK with best position
    dbest_fit=sorted_ELKS_fitness(1); % the fitness of best position

    % Number of females for each male
    for i=1:num_of
        TransposeFitness(i)=sorted_ELKS_fitness(i);
    end

    %% Predation risk-avoidance migration mechanism
    A=6*exp(-2.6*(randi([1 2]))*((iteration)/iter_max))./(1+0.12);

    if A>3
        for i=1:num
            randomNumber1=randi([1, num]);
            thta=2*pi*rand;
            D2=(dbest-dc(i))*rand+dc(i);
            B=1.2-1.2*cos(pi*(1-iteration/iter_max)/2);
            D(i)=D2+(dc(randomNumber1)-dc(i))*cos(thta)*(1/A)*B*rand+dc(randomNumber1)*sin(thta)*(1/A)*B*rand;
        end
    else
        %%
        Familes=zeros(1,num);
        for i=(num_of+1):num
            FemaleIndex=sorted_indexes(i); % index of female
            randNumber = rand;
            MaleIndex = 0;
            sum_fitness=0;

            for j=1:num_of
                sum_fitness = sum_fitness + (TransposeFitness(j)/sum(TransposeFitness));
                if (sum_fitness > randNumber)
                    MaleIndex =  j;
                    break;
                end
            end
            Familes(FemaleIndex)=sorted_indexes(MaleIndex);
        end

        for i=1:num
            %Male
            if(Familes(i)==0)
                h=fix(rand*num)+1;
                D(i)=dc(i)+rand*(dc(h)-dc(i));
                %% Triangle walk strategy
                rg=2*exp(-(2*iteration/iter_max)^2);
                L = dbest-D(i);
                LP = L*rand;
                alph = L*L+LP*LP-2*LP*L*cos(2*pi*rand);
                D(i)=dbest+rg*alph;
            else
                h=fix(rand*num)+1;
                MaleIndex=Familes(i);
                hh = randperm(size(find(Familes==MaleIndex),2));
                h=round(1+(size(hh,2)-1)*rand);
                rd = -2+4*rand;
                D(i)=dc(i)+(dc(Familes(i))-dc(i))+rd*(dc(h)-dc(i));
            end
        end

    end
    
    for i=1:num
        dup_uplim=(D(i)-1)<0;
        dup_lowlim=D(i)>0;
        D(i)=D(i).*dup_lowlim.*dup_uplim+(~dup_uplim);
    end

   %% Memory-guided redirection strategy
        for i=1:num
            [p,K]=Memorydetection(Memory,Memoryf,D(i),i,p,K,iteration);
        end
        
end

%% Memory-guided redirection strategy

function [p,K,D]=Memorydetection(Memory,Memoryf,D,i,p,K,iteration)
% global githowdes
which1=findMemory(Memory,D(1),iteration);
if isempty(which1)
    p(i)=p(i);
    K(i)=0;
else
    [m,j]=sort(Memoryf,'descend');
    D=Memory(j(1))+rand*(Memory(j(1))-D);
    which2=findMemory(Memory,D(1),iteration);
    if isempty(which2)
        p(i)=p(i);
        K(i)=0;
    else
        p(i)=Memoryf(which2(1));
        K(i)=i;
    end
end
end

function which1=findMemory(Memory,pop,iteration)
if  iteration<=3
    which1=find(abs(Memory-pop(1))<=0.03);
else
     which1=find(abs(Memory-pop(1))<=0.0001);
end
end

function [Mw,Memory,Memoryf]=saveMemory(Mmax,Mw,Memory,Memoryf,pop,popf)
Mw=Mw+1;
if Mw>Mmax
    Mw=1;
    [Memoryf,sy]=sort(Memoryf,"ascend");
    Memory=Memory(:,sy);
end
Memory(Mw)=pop(1);
Memoryf(Mw)=popf(1);
end