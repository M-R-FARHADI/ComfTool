function q = pmveq2(x,mat,j,jj)
q(1)=35.7-0.028*(mat(1)-mat(2))-mat(5)...
    *( 3.96e-8*mat(7)*( (x(1)+273)^4 - (jj+273)^4 )...
    + mat(7)*x(2)*(x(1) -j) )  - x(1);
if (x(1)-j)^0.25 >= 12.1*(mat(6))^0.5
    q(2)=2.28*(x(1)-j)^0.25 - x(2);
else
    q(2)=12.1*mat(6)^0.5 -x(2);
end
end