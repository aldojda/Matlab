%Sooting method.....
%Rungekutta
Rungekutta(Ev)
hbar=1;
m=1;
a=1;                     
h=0.001;      
hbar=1;
ka=1;
a=1;      %Tamaño del pozo
%Ev=pi^2/2;
% step size
x = -a/2:h:a/2;              %Vector del espacio que contiene el pozo
y = zeros(1,length(x));      %Funcion
z = zeros(1,length(x));
%Condiciones de frontera
y(1) = 0;         % condición inicial
z(1) = 1;         % condición final
     f = @(x,y,z) z;
    g = @(x,y,z) -2*ka/(hbar^2)*Ev*y;
    %    Runge-Kutta method.
    for i=1:(length(x)-1)
        k1 = f(x(i),y(i),z(i));
        l1 = g(x(i),y(i),z(i));
        k2 = f(x(i)+0.5*h,y(i)+0.5*h*k1,z(i)+0.5*h*l1);
        l2 = g(x(i)+0.5*h,y(i)+0.5*h*k1,z(i)+0.5*h*l1);
        k3 = f((x(i)+0.5*h),(y(i)+0.5*h*k2),(z(i)+0.5*h*l2));
        l3 = g((x(i)+0.5*h),(y(i)+0.5*h*k2),(z(i)+0.5*h*l2));
        k4 = f((x(i)+h),(y(i)+k3*h),(z(i)+l3*h));
        l4 = g((x(i)+h),(y(i)+k3*h),(z(i)+l3*h));
        y(i+1) = y(i) + (1/6)*(k1+2*k2+2*k3+k4)*h;
        z(i+1) = z(i) + (1/6)*(l1+2*l2+2*l3+l4)*h;
   
    plot(x,y)
    %% 






%DE: D´´Z=-(2*m/h)(E-v(x)); v(x)=0 in all space except |a|>x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Declaración de variables y constantes
hbar=1;
ka=1;
a=1;                         %Tamaño del pozo
h=0.01;        % step size
EVnumb=0;
x = -a/2:h:a/2;   %Vector del espacio que contiene el pozo
y = zeros(1,length(x));             %Funcion
z = zeros(1,length(x));
EV=zeros(1,length(x)-2);%Derivada de la función
Evnumb=zeros(1,900);
EVec=zeros(1,26);                   % Vector d Eigenvalores de Salida

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eigf=0.2:3.6:3000%Vector de Eigenvalores
Eigf2=zeros(1,25)
%El siguiente proceso evalua ED para diferentes valores de eigenvalores
%con el fin de encontrar los valores extremos para comenzar con la busqueda
%del eigenvalor con el error mínimo.

for n=1:length(Eigf)
    %Funciónes de ecuaciones diferenciales lineales
   
    end
    %Evaluar si se encuentran dos Eigenvectores consecutivos
    EV(1,:)=y(2:end-1).*y(3:end);
    EV(EV(1,:)<0) = 0;
    Cero=find(~EV(1,:));
    Evnumb(n)=length(Cero)+1 ;%Define the i-eigenvector of the function
    
end
%Generar los limites para encontrar los posibles 25 Eigenvalores en sus respectivos limites
for k=1:25
    EVec(k)=find(Evnumb==k,1,'last'); % Negative value
end
%Vector de Eigenvalores extremos para la busqueda de El eigenvalor con el
%mínimo error.
Eigf=Eigf(EVec(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%El proceso anterior fué para encontrar las posibles ubicaciones de los eigenvalores

%Proceso de bisección para los 25 eigenvalores...

for k=1:2
    Err=1;
    while Err>0.00001 
        disp(k)
            %Funciónes de ecuaciones diferenciales lineales
            f = @(x,y,z) z;
            g = @(x,y,z) -2*m/(hbar^2)*Eigf(k)*y;
            %    Runge-Kutta method.
            for i=1:(length(x)-1)
                k1 = f(x(i),y(i),z(i));
                l1 = g(x(i),y(i),z(i));
                k2 = f(x(i)+0.5*h,y(i)+0.5*h*k1,z(i)+0.5*h*l1);
                l2 = g(x(i)+0.5*h,y(i)+0.5*h*k1,z(i)+0.5*h*l1);
                k3 = f((x(i)+0.5*h),(y(i)+0.5*h*k2),(z(i)+0.5*h*l2));
                l3 = g((x(i)+0.5*h),(y(i)+0.5*h*k2),(z(i)+0.5*h*l2));
                k4 = f((x(i)+h),(y(i)+k3*h),(z(i)+l3*h));
                l4 = g((x(i)+h),(y(i)+k3*h),(z(i)+l3*h));
                y(i+1) = y(i) + (1/6)*(k1+2*k2+2*k3+k4)*h;
                z(i+1) = z(i) + (1/6)*(l1+2*l2+2*l3+l4)*h;
        end
        Err=abs(y(k));
    Eigf(k)=Eigf(k)-y(end)/z(end);
    end
end
%Process to find number of zeros in the function.
y1=y(end);
plot(x,y)

%Count the number of zeros, bisection method
