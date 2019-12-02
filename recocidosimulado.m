% close all
 clc
N=8;
b=20;
a=0;
kb=1;
t=1;
rest=15;
XY=zeros(2,N);
XY(1,:)=[1,-2,3,-5,9,3,-12,0];
XY(2,:)=[2,4,-6,8,0,9,7,0];
fin=10000;
matrizres=zeros(rest,N);
distres=zeros(1,rest);
for arc=1:rest
    resulti=randperm(N);
for lol=2:1:fin 
    t=t-1e-4; t=t*.85;
    resultaux=randperm(N);
    distact=0; distprue=0;
    for a=1:(N-1)
   distact=norm(XY(:,resulti(a+1))-XY(:,resulti(a)))+distact;
   distprue= norm(XY(:,resultaux(a+1))-XY(:,resultaux(a)))+distprue;
    end
    distact= norm(XY(:,resulti(1))-XY(:,resulti(N)))+distact;
   distprue= norm(XY(:,resultaux(1))-XY(:,resultaux(N)))+distprue;
   D_E= distprue-distact;
   if distprue<distact
       resulti=resultaux;
       resultaux=resulti;
   elseif distprue>=distact
           if rand(1,1)<=exp(D_E/kb/t)
              resulti=resultaux;
              resultaux=resulti; 
           else
               resulti=resulti;
               resultaux=resulti;
          end 
   end
end
distres(arc)= distact;
matrizres(arc,:)= resulti;
end
distres
matrizres