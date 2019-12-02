
% %% Problema 1
% %Molecular Dynamics
tic
% clc;
% clear;
% close all
%%%%%%%%%%%%%%%
N=100;%50
Space=50;%10
h=0.06;
Temp=100;  
Step1=12000;
 Step2=200;
%%%%%%%%%%%%%%
% Initial condition
%Xi position
Xi= Space*rand([1,N,3]);
%Xi= repmat(Space*rand([1,N]),[1,1,3]);
[U1]=Force_U(Xi,N,2);
Uvec=zeros(1,Step1);
Prob=zeros(1,Step1);
Probability=1;
%Create new random cordinate with one diferent position
for jj=1:Step1
    Xi_new=Xi(1,:,:);
    N1=randi([1,N]);
    %Xi_new(1,N1,randi([1,3]))=Space*rand;
     Xi_new(1,N1,:)=Xi_new(1,N1,:)+((0.4).*rand - 0.2);
    for ii=1:3
        if Xi_new(1,N1,ii)>Space
            Xi_new(1,N1,ii)=Xi_new(1,N1,ii)-Space;
        end
    end
    [U2]=Force_U(Xi_new,N,2);
    %Proceso de minimización de energía
    Del_U=U2-U1;
    Temp=Temp*0.9999^(jj);
    if U1>=U2
        Xi=Xi_new;
        U1=U2;

    elseif U1<U2
        Probability=exp(-(Del_U)/Temp);
        if  Probability>= rand
            Xi=Xi_new;
            U1=U2;
        end
    end
    Uvec(1,jj)=U1;
    Prob(jj)=Probability;
end
figure(1)
plot(Uvec)
% % figure(2)
% % plot(Prob)
% % disp(U1)
figure(4)
scatter3(Xi(1,:,1),Xi(1,:,2),Xi(1,:,3),...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
    xlim([0 Space]);
    ylim([0 Space]);
    zlim([0 Space]);
  X_inicio=Xi;
toc
%[Xi]=Config(N,Space);
%% Evolución en el tiempo
%Determinar Xi vec(3 variables)
to=1;
%v_i=Xi*0;
v_i = to*normrnd(0,1,[1,N,3]);
for jj=1:Step2
    T(jj)=1/(3*N)*sum(sum(v_i.^2));
    [F_i]=Force_U(Xi,N,1);
    X_new=Xi+h*v_i+h^2/2*F_i;
    while isempty(X_new( X_new<0))==0
    X_new( X_new<0)=X_new( X_new<0)+Space;
    end
    while isempty(X_new( X_new>Space))==0
    X_new( X_new>Space)=X_new( X_new>Space)-Space;
    end
    [F_new]=Force_U(X_new,N,1);
    v_new=v_i+h*(F_new+F_i)/2;
    figure(2)
    
    scatter3(Xi(1,:,1),Xi(1,:,2),Xi(1,:,3),...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
    xlim([0 Space]);
    ylim([0 Space]);
    zlim([0 Space]);
   figure(3)
    plot(T);
    Xi=X_new;
    v_i=v_new;
end
%Definir velocidades
figure(4)
plot(T)

toc