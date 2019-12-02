function[OutFU]=Force_U(Xi,N,FU)
%Calcular matriz R cuadrado
Xi= repmat(Xi,N,1);
for i=1:3
    Del_Xi(:,:,i)=Xi(:,:,i)'-Xi(:,:,i);
end
R2=Del_Xi(:,:,1).^2+Del_Xi(:,:,2).^2+Del_Xi(:,:,3).^2;
if FU==1
    %Calcular matriz fuerza
    F=-12*R2.^(-7)+24*R2.^(-4);
   % F(F==Inf)=0;
    F(isnan(F))=0;
    Fvec=-F.*Del_Xi;
    %F_tot is -Fi,j
    OutFU=sum(Fvec);
elseif FU==2
    %Calcular Matrix potencial
    U_mat=1.*R2.^(-6)-4.*R2.^(-3);
    U_mat(isnan(U_mat))=0;
    OutFU=sum(sum(U_mat));
end
end
