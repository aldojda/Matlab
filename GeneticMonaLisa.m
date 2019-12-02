clc;
clear;
tic
imdata = imread('PS4.jpeg');
imdata=rgb2gray(imdata);
imdata2=imdata;
imdata = double(imdata);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
N_shapes=100;                %Should be even number%% 
Popu=30;                     %Should be even number%% 
Chrom=36;                             %Don't change%% 
MD_Success = 0.9;                %Reproduction rate%% 
Gen=1;                   %Start value of Generation%% 
P=0.01;                     %Pobability of Mutation%% 
End_Gen=120000;              %Number of Generations%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
a=1;
Pix=size(imdata,1);
Gioconda=ones([Pix,Pix,Popu]);
Fitness=zeros(Popu,1);
BetterFitness=zeros(1,End_Gen/10);
BetterGio=zeros([Pix,Pix,End_Gen/10]);
%Generate matrix with poblation, shapes and colors;
%[1 0 0]=N_shapes, number of squares, [0 1 0]=Number of chromosomes
%[0 0 1]= Population
%7->x1;7->y1;7->x2;7->y2;Color->8:total 36 chromosomes
%Initial random Matrix
Matrix=randi([0 1], [N_shapes,Chrom,Popu]);
NewGen=Matrix;
Vars=zeros(N_shapes,5,Popu);
BetterFitness(1,1)=1e6;

while (Gen<End_Gen)
    Mutation_P=P;
    %Conversion of binary-> Number in order to get the fitness
    for k=1:Popu
        Vars(:,5,k)=bi2de(Matrix(:,29:Chrom,k));%Color
        Vars(:,1,k)=bi2de(Matrix(:,1:7,k))+1;%X1
        Vars(:,2,k)=bi2de(Matrix(:,8:14,k))+1;%X2
        Vars(:,3,k)=bi2de(Matrix(:,15:21,k))+1;%Y1
        Vars(:,4,k)=bi2de(Matrix(:,22:28,k))+1;%Y2
        %Order Vars in scale of colors
        [~ , index_sort] = sort(Matrix(:,:,k),1,'descend');
        Vars(:,:,k)=Vars(index_sort(:,end),:,k);
        %Reorder x,y min and max
        X(:,1)=min(Vars(:,1:2,k).').';
        X(:,2)=max(Vars(:,1:2,k).').';
        X(:,3)=min(Vars(:,3:4,k).').';
        X(:,4)=max(Vars(:,3:4,k).').';
        X(:,5)=Vars(:,5,k);
        %Make imaging and fitness
        for l=1:N_shapes
            Square=X(l,5)*ones(X(l,2)-X(l,1)+1,X(l,4)-X(l,3)+1);
            Gioconda(X(l,1):X(l,2),X(l,3):X(l,4),k)=(Gioconda(X(l,1):X(l,2),X(l,3):X(l,4),k)+Square)./2;
        end
        %Calculate Fitness and stored
        Fitness(k,1)=sum(sum(abs(imdata-Gioconda(:,:,k))));
       % Fitness(k,1)=abs(sum(sum(imdata(10:end-10,10:end-10)))-sum(sum(Gioconda(10:end-10,10:end-10,k))));
    end
    Gioconda=floor(Gioconda);
    F1=1./Fitness;
    Prob=F1/sum(F1);
    BestFit=find(Fitness==min(Fitness));
    Test1=repmat(cumsum(Prob)',Popu,1);
    Test0=[zeros(Popu,1) Test1(:,1:end-1)];
    %Natural selection
    %Matrix test conditional generate the new population w/
    Test=repmat(rand(Popu,1),1,Popu);
    Final=sum((Test0<Test).*(Test<Test1)).';
    Acom=1;
    %Save old values in order to get control
    %Remember save the one with the best Fitness...
    OldGen=Matrix;
    for k=1:Popu
        for kp=1:Final(k)
            NewGen(:,:,Acom)=Matrix(:,:,k);
            Acom=Acom+1;
        end
    end

    %Process of reproduction
    %Selection of the Father and mother:
    Dad=randi([1,Popu],Popu/2,1);
    Mom=randi([1,Popu],Popu/2,1);
    Dad(1)=1;
    Dad(2)=2;
    %Success of reproduction
    Success_Rep=rand(Popu/2,1)< MD_Success;
    
    for k=1:Popu/2
        NewGen(:,:,2*k-1)=Matrix(:,:,Dad(k));
        NewGen(:,:,2*k)=Matrix(:,:,Mom(k));
        if Dad(k)==Mom(k)
            %Be sure dad and mom aren´t equal
            Success_Rep(k)=0;
        end
        if Success_Rep(k)==1
            %Reproduction and permutation
            S1=randi([1,Chrom],2,1);
            S1(2)=randi([1,N_shapes],1,1);
            NewGen(S1(2):1:end,:,[2*k,2*k-1])=NewGen(S1(2):1:end,:,[2*k-1,2*k]);
            NewGen(S1(2),1:S1(1),[2*k,2*k-1])=NewGen(S1(2),1:S1(1),[2*k-1,2*k]);
            
        end
    end
        %Add those with better fitness
    SortFit=sort(Fitness);
    Find=[find(Fitness==SortFit(1)),find(Fitness==SortFit(2))];
    NewGen(:,:,1)=Matrix(:,:,Find(1));
    NewGen(:,:,2)=Matrix(:,:,Find(2));
    %Mutation

    for kp=1:Popu
        if rand < Mutation_P
            %Occurs the mutation
            Mut_Bit=500;%randi([1,10]);
            MutV=[randi([1,Chrom],Mut_Bit,1),randi([1,N_shapes],Mut_Bit,1)];
            for kpp=1:Mut_Bit
                if NewGen(MutV(kpp,2),MutV(kpp,1),kp)==1
                    NewGen(MutV(kpp,2),MutV(kpp,1),kp)=0;
                else
                    NewGen(MutV(kpp,2),MutV(kpp,1),kp)=1;
                end
            end
        end
    end
 
    %Save results;
    Matrix=NewGen;
    if mod(Gen,10)==0
    BetterFitness(1,a+1)=min(Fitness);
    BetterGio(:,:,a)=Gioconda(:,:,BestFit);
    a=a+1;
    end
    Gen=Gen+1;
end
BetterGio=uint8(BetterGio);
toc
subplot(1,2,1)
imshow(BetterGio(:,:,end-1)); 
subplot(1,2,2)
imshow(imdata2)
figure(2)
plot(BetterFitness(2:end))