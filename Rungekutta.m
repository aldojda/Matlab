%Rungekutta Eigenvalores
%   Este codigo despliega genera la solucion de una ED por el método de Rk
function[x,y,z]=Rungekutta(Ev)
hbar=1;
m=1;
a=1;                     
h=0.001;      
hbar=1;
ka=1;
a=1;      %Tamaño del pozo

x = -a/2:h:a/2;              %Vector del espacio que contiene el pozo
y = zeros(1,length(x));      %Funcion
z = zeros(1,length(x));
%Condiciones de frontera
y(1) = 0;         % condición inicial
z(1) = 1;         % condición final
     f = @(x,y,z) z;
    g = @(x,y,z) -2*ka/(hbar^2)*Ev*y;
    %    Runge-Kutta orden 4
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
end
