function[Distance]=TSalesman(City,Xv,Yv)
close all;
clc;
MDistance=zeros(City);
Xpath=zeros(City,1);
Ypath=zeros(City,1);
% matrix of distance.
for i=1:City
    for j=1:City
        MDistance(i,j)=sqrt((Xv(i)-Xv(j))^2+(Yv(i)-Yv(j))^2);
    end
end
%Definición del caso inicial, orden aleatorio
CityNum= randperm(City);
CityNum2=[CityNum(2:end),CityNum(1)];
%Calcular la distancia
Distance=0;
for j=1:City
    Distance=Distance+MDistance(CityNum(j),CityNum2(j));
end
disp(Distance);
%%%%%SIGUE DETERMINAR Cual valor cambiar aleatoreamente y la funcion de
%%%%%probabilidad
for k=1:2000
    %Seleecionar particula i Al azar y cambiarla
    ChangeCity=randperm(City,2);
    CityNumNew=CityNum;
    CityNumNew([ChangeCity(1),ChangeCity(2)])=CityNum([ChangeCity(2),ChangeCity(1)]);
    CityNum2New=[CityNumNew(2:end),CityNumNew(1)];
    %   Evaluar distancia Caso Nuevo
    DistanceNew=0;
    for j=1:City
        DistanceNew=DistanceNew+MDistance(CityNumNew(j),CityNum2New(j));
    end
    disp(DistanceNew)
    if Distance>=DistanceNew
        CityNum=CityNumNew;
        CityNum2=CityNum2New;
        Distance=DistanceNew;
    elseif Distance<DistanceNew
        Probability=exp(-(DistanceNew-Distance)/(0.5*City*0.3^(k/2+1)));
        
        if  Probability>= rand
            CityNum=CityNumNew;
            CityNum2=CityNum2New;
            Distance=DistanceNew;
            
        end
    end
end
%Proceso para plotear
for i=1:City 
    Xpath(i)=Xv(CityNum(i));
    Ypath(i)=Yv(CityNum(i));
    %1000 is the length of x_vector and y_vector
    plot(Xv(i), Yv(i), 'o')
     hold on
end
plot([Xpath;Xpath(1)],[Ypath;Ypath(1)])
disp(CityNum)
disp(DistanceNew)
end