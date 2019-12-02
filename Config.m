
function[Xi]=Config(N,Space)
%% Problema 1
%Molecular Dynamics
tic
%%%%%%%%%%%%%%%
N=N;
Space=Space;
pobla=20;
partici=2;
generacionF=5000;
generacion=0;
xop=0.65;% probabilidad de cruce
mutap=0.2;% probabilidad de mutacion
%%%%%%%%%%%%%%%
% Initial condition
%Xi position
Xi= Space*rand([1,N,3,pobla]);
while generacion<generacionF
    Xir= repmat(Xi,N,1);
    for j=1:pobla
        for i=1:3
            Del_Xi(:,:,i,j)=Xir(:,:,i,j)-Xir(:,:,i,j)';
        end
        R2=Del_Xi(:,:,1,j).^2+Del_Xi(:,:,2,j).^2+Del_Xi(:,:,3,j).^2;
        U_mat=1.*R2.^(-6)-4.*R2.^(-3);
        U_mat(isnan(U_mat))=0;
        OutFU(j)=sum(sum(U_mat));
    end
    % for j=1:pobla
    %  OutFU(j)=Force_U(Xi,N,2);
    % end
    
    cand_elite=floor(0.1*pobla);% aproximadamente un 10% de la poblacion pasa a
    %la siguiente generacion por elitismo
    OutFU1=OutFU;
    for el=1:cand_elite
        [min1,ind]=min(OutFU1);
        elites(el)=ind;
        OutFU1(ind)=Inf;
    end
    
    pama=0;% cantidad de padres y madres
    alma=[];% aqui se almacenan aquellos que no son candidatos para cruces (elites}
    %y aquellos que ya fueron escogidos como padre/madre
    alma=[alma elites];
    while pama<pobla/2% la mitad de la poblacion sera padre/madre pero no todos se cruzaran
        parti=randperm(pobla,partici);% se escoge a los participantes del torneo
        while isempty(intersect(parti,alma))==0% se checa si estas no fueron escogidos ya
            parti=randperm(pobla,partici);
        end
        Xir= repmat(Xi,N,1);
        for j=1:length(parti)
            for i=1:3
                Del_Xi(:,:,i,j)=Xir(:,:,i,j)-Xir(:,:,i,j)';
            end
            R2=Del_Xi(:,:,1,j).^2+Del_Xi(:,:,2,j).^2+Del_Xi(:,:,3,j).^2;
            U_mat=1.*R2.^(-6)-4.*R2.^(-3);
            U_mat(isnan(U_mat))=0;
            OutFU(j)=sum(sum(U_mat));
        end
        OutFUpadres=OutFU(parti);
        padreind=find(OutFUpadres==min(OutFUpadres(:)));% se escoge como padre al que tenga
        %menos distancia (mayor fitness)
        pama=pama+1;
        alma=[alma parti(padreind(1))];% el ganador del torneo se incorpora a
        %este vector para que no vuelva a ser escogido
    end
    alma=alma(length(elites)+1:end);
    alma2=[];
    alma2=[alma2 elites];% alma2 almacenara los que si se cruzan y los elites,
    %es decir almacenara todo aquel que no tendra mutaciones
    for m=1:2:length(alma)-1
        xoprob=rand;
        %     if generacion<250% las ´primeras 100 generaciones tendran mayor probabilidad de mutacion
        %         xoprob=xoprob-.5;
        %     end
        if xoprob<xop% se checa si se cruza
            xopointx=randi([1,N],1,1);% se selecciona un punto de corte para la distancia
            xopointy=randi([1,N],1,1);% se selecciona un punto de corte para el angulo de disparo
            xopointz=randi([1,N],1,1);
            % se separa el adn el cual en este caso es la matriz mov
            string1=Xi(1,1:xopointx,1,m);
            string1_1=Xi(1,xopointx+1:N,1,m);
            string2=Xi(1,1:xopointx,1,m+1);
            string2_1=Xi(1,xopointx+1:N,1,m+1);
            
            string3=Xi(1,1:xopointy,2,m);
            string3_1=Xi(1,xopointy+1:N,2,m);
            string4=Xi(1,1:xopointy,2,m+1);
            string4_1=Xi(1,xopointy+1:N,2,m+1);
            
            string5=Xi(1,1:xopointz,3,m);
            string5_1=Xi(1,xopointz+1:N,3,m);
            string6=Xi(1,1:xopointz,3,m+1);
            string6_1=Xi(1,xopointz+1:N,3,m+1);
            
            % se reforma el adn el cual en este caso es la matriz mov
            Xi(:,:,1,alma(m))=cat(2,string1,string2_1);
            Xi(:,:,1,alma(m+1))=cat(2,string2,string1_1);
            Xi(:,:,2,alma(m))=cat(2,string3,string4_1);
            Xi(:,:,2,alma(m+1))=cat(2,string4,string3_1);
            Xi(:,:,3,alma(m))=cat(2,string5,string6_1);
            Xi(:,:,3,alma(m+1))=cat(2,string6,string5_1);
            alma2=[alma2 alma(m) alma(m+1)];% se almacenan los dos padres
        end
    end
    % Mutaciones
    a=1:pobla;
    nocruza=setxor(a,alma2);% se selecciona todo aquel que no se cruzo y que no es elite
    for m=1:length(nocruza)
        mutaprob=rand;
        if generacion<500% las ´primeras 100 generaciones tendran mayor probabilidad de mutacion
            mutaprob=mutaprob-.5;
        end
        if mutaprob<mutap
            i=randi([1,N],1,1);
            Xi(1,i,:,nocruza(m))=-Space+2*Space*rand([1,1,3]);
        end
    end
    Xir= repmat(Xi,N,1);
    for j=1:pobla
        for i=1:3
            Del_Xi(:,:,i,j)=Xir(:,:,i,j)-Xir(:,:,i,j)';
        end
        R2=Del_Xi(:,:,1,j).^2+Del_Xi(:,:,2,j).^2+Del_Xi(:,:,3,j).^2;
        U_mat=1.*R2.^(-6)-4.*R2.^(-3);
        U_mat(isnan(U_mat))=0;
        OutFU(j)=sum(sum(U_mat));
    end
    
    generacion=generacion+1;
    U_mejor(generacion)=unique(min(OutFU));
end
x=1:generacion;
figure(1)
plot(x,U_mejor)
index=find(OutFU==min(OutFU(:)));
mejor_indi=index(1);
minenergy=OutFU(mejor_indi)
figure(2)
scatter3(Xi(:,:,1,mejor_indi),Xi(:,:,2,mejor_indi),Xi(:,:,3,mejor_indi),'filled')
toc
Xi= Xi(:,:,:,mejor_indi);
end