clc
clear
%Vector de eigenvelores y eigenfunciones
Evl=0.02:3:3600;
Row=length(Evl);
Efc=zeros(Row,1001);
Nod=zeros(1,Row);
Nodnumb=Nod;
Limites=zeros(1,50);
Results=zeros(50,1001);
EResults=zeros(1,25);
%Generar Eigenvectores de los Eigenvalores principales
for i=1:Row
    [x,y,z]=Rungekutta(Evl(i));
    Efc(i,:)=y;
    % Determinar los nodos de cada eigenvector
    Nodos=Efc(i,2:end-1).*Efc(i,3:end);
    Nodos(Nodos(1,:)<0)=0;
    Nodos=find(~Nodos(1,:));
    Nodnumb(i)=length(Nodos)+1;
end

%Generar los limites para encontrar los posibles 25 Eigenvalores en sus respectivos limites
for i= 1:25
    Limites(2*i-1)=find(Nodnumb==i,1,'last'); 
    Limites(2*i)=find(Nodnumb==i+1,1,'first');
end
Evl=Evl(Limites);
Efc=Efc(Limites,:);
error=1;
%Proceso de biseccion para encontrar valor optimo
for i=1:25
    error=1;
    disp(i);
    while error > 1e-7
        
        %Eigen valor medio
        Evlm=(Evl(2*i-1)+Evl(2*i))/2;
        %Valor de la funcion en el punto medio
        [x,y,z]=Rungekutta(Evlm);
        Efcm=y;
        %Metodo de bisección
        if Efcm(end)> 0
            if Efc(2*i-1,end)>Efcm(end)& Efc(2*i,end)<Efcm(end)
                Efc(2*i-1,:)=Efcm;
                Evl(2*i-1)=Evlm;
                
            elseif Efc(2*i-1,end)<Efcm(end)& Efc(2*i,end)>Efcm(end)
                Efc(2*i,:)=Efcm;
                Evl(2*i)=Evlm;
            end       
        elseif Efcm(end)<0
            if Efc(2*i-1,end)<Efcm(end)& Efc(2*i,end)>Efcm(end)
                Efc(2*i-1,:)=Efcm;
                Evl(2*i-1)=Evlm;
            elseif Efc(2*i-1,end)>Efcm(end)& Efc(2*i,end)<Efcm(end)
                Efc(2*i,:)=Efcm;
                Evl(2*i)=Evlm;
            end
        end
        error=abs(Efcm(end));
    end
    %Guardar resultados en Vector
    Results(i,:)=Efcm;
    EResults(i)=Evlm;
end
%Graficar eigenfunciones 1 a 25
hold on
for i=1:25
    
    plot(x,Results(i,:))
    pause(0.9)
end
hold off
figure(2)
plot(EResults*2/pi^2,'o')
