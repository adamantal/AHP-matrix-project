clear all;

matrixmeret = 4;
NemParetokSzama = 0 ;
ma12 = 7;    ma13 = 6;    ma14 = 5;
             ma23 = 1/2;  ma24 = 1; 
                          ma34 = 1/2;
A = [          1          ma12        ma13      ma14        ; ...
              1/ma12       1          ma23      ma24        ; ...
              1/ma13      1/ma23       1        ma34        ; ...
              1/ma14      1/ma24      1/ma34     1           ]; 

w = [0.64737128 , 0.09361541 , 0.11588291 , 0.14313040 ]
I1elemszam = 0;
I0elemszam = 0;
for indexi = 1:matrixmeret
  v(indexi) = log(w(indexi));
for indexj = 1:matrixmeret
  b(indexi,indexj) = log(A(indexi,indexj));
  if real(w(indexi)/w(indexj)) - A(indexi,indexj) > 10^(-6) 
      I1elemszam = I1elemszam + 1 ;
      I1poziciok(I1elemszam,1) = indexi ;
      I1poziciok(I1elemszam,2) = indexj ;
  end;
end;  % for indexj = 1:matrixmeret
end;  % for indexi = 1:matrixmeret
for indexi = 1:matrixmeret-1
for indexj = indexi+1:matrixmeret
  if abs(w(indexi)/w(indexj) - A(indexi,indexj)) < 10^(-8)   
      I0elemszam = I0elemszam + 1  ;
      I0poziciok(I0elemszam,1) = indexi ;
      I0poziciok(I0elemszam,2) = indexj ;
  end;
end;  % for indexj = 1:matrixmeret
end;  % for indexi = 1:matrixmeret

IneqSor = 0;     
for I1index = 1:I1elemszam
    IneqSor = IneqSor + 1 ;
    indexi = I1poziciok(I1index,1);
    indexj = I1poziciok(I1index,2);    
    Ineq(IneqSor,indexi) = -1 ;
    Ineq(IneqSor,indexj) =  1 ;  
    RHSineq(IneqSor) = -b(indexi,indexj) ;
    IneqSor = IneqSor + 1 ;
    Ineq(IneqSor,indexi) =   1 ;
    Ineq(IneqSor,indexj) =  -1 ;  
    Ineq(IneqSor,matrixmeret+I1index) =  1 ;      
    RHSineq(IneqSor) = v(indexi)-v(indexj) ;
end;   % for I1index = 1:I1elemszam
for indexi = 1:matrixmeret
    c(indexi) = 0 ;
end;    % for indexi = 1:matrixmeret
for indexi = matrixmeret+1:matrixmeret+I1elemszam
    c(indexi) = -1 ;
end;    % for indexi = matrixmeret+1:matrixmeret+I1elemszam
for indexi = 1:matrixmeret
    lb(indexi) = -10 ;
end;    % for indexi = 1:matrixmeret
for indexi = matrixmeret+1:matrixmeret+I1elemszam
    lb(indexi) = 0 ;
end;    % for indexi = matrixmeret+1:matrixmeret+I1elemszam
for indexi = 1:matrixmeret+I1elemszam
    ub(indexi) = 100 ;
end;    % for indexi = matrixmeret+1:matrixmeret+I1elemszam
for indexi = 1:matrixmeret
   start(indexi) = v(indexi);
end;   % for indexi = 1:matrixmeret-1
for indexi = matrixmeret+1:matrixmeret+I1elemszam
   start(indexi) = 0;
end;   % for indexi = matrixmeret+1:matrixmeret+I1elemszam
for Eqindex = 1:matrixmeret+I1elemszam
    Eq(1,Eqindex) = 0;
end;   % for Eindex = matrixmeret+1:matrixmeret+I1elemszam
Eq(1,1) = 1;
RHSeq(1) = 0 ;
for I0index = 1:I0elemszam
    indexi = I0poziciok(I0index,1);
    indexj = I0poziciok(I0index,2);    
    Eq(I0index+1,indexi) =  1 ;
    Eq(I0index+1,indexj) =  -1 ;  
    RHSeq(I0index+1) = -b(indexj,indexi) ;
end;   % for I0index = 1:I0elemszam
if I1elemszam > 0
[xopt,fopt]=linprog(c,Ineq,RHSineq,Eq,RHSeq,lb,ub,start);
else
[xopt,fopt]=linprog(c,[],[],Eq,RHSeq,lb,ub,start);    
end;    
if fopt < -10^(-6) 
   NemParetokSzama = NemParetokSzama + 1 ;
   disp('Nem Pareto optimalis!');
   for indexi = 1:matrixmeret
    dw(indexi) = exp(xopt(indexi));
   end;   % for indexi = 1:matrixmeret
   disp('Az alabbi Pareto optimalis sulyvektor dominalja:');
   disp(dw);
   NemPareto(NemParetokSzama,1) = NemParetokSzama;
   NemPareto(NemParetokSzama,2) = ma12;
   NemPareto(NemParetokSzama,3) = ma13;
   NemPareto(NemParetokSzama,4) = ma14;
   NemPareto(NemParetokSzama,5) = ma23;
   NemPareto(NemParetokSzama,6) = ma24;
   NemPareto(NemParetokSzama,7) = ma34;
   NemPareto(NemParetokSzama,8) = -1;
   %NemPareto(NemParetokSzama,9) = CR;
   NemPareto(NemParetokSzama,10) = w(1);
   NemPareto(NemParetokSzama,11) = w(2);
   NemPareto(NemParetokSzama,12) = w(3);
   NemPareto(NemParetokSzama,13) = w(4);
   NemPareto(NemParetokSzama,14) = -8;
   NemPareto(NemParetokSzama,15) = dw(1);
   NemPareto(NemParetokSzama,16) = dw(2);
   NemPareto(NemParetokSzama,17) = dw(3);
   NemPareto(NemParetokSzama,18) = dw(4);
end;
