clear;
A  = [ ...  // Kou Lin Example 1
1   4    3   1   3    4    ;  ...
0   1    7   3   1/5  1    ;  ...
0   0    1   1/5 1/5 1/6   ;  ...
0   0    0   1   1   1/3   ;  ...
0   0    0   0   1   3     ;  ...
0   0    0   0   0   1     ];


matrixsize = max(size(A));
for indexi = 1:matrixsize-1
    for indexj = indexi+1:matrixsize
        if A(indexi,indexj) == 0
        else
           A(indexj,indexi) = 1/A(indexi,indexj) ;
        end;    
    end;       
end;        // for indexi = 1:matrixsize-1

B = zeros(matrixsize,matrixsize);  
for indexi = 1:matrixsize
    s = 0 ;
    for indexk = 1:matrixsize
        s = s + (A(indexk,indexi))^2
    end;  //      for indexj = 1:matrixsize
    Aoszlopnegyzetosszeg(indexi) = s^(1/2) ;
end;     // for indexi = 1:matrixsize-1
for indexi = 1:matrixsize
    for indexj = 1:matrixsize
        B(indexi,indexj) = A(indexi,indexj)/Aoszlopnegyzetosszeg(indexj) ;
    end;  //      for indexj = 1:matrixsize
end;     // for indexi = 1:matrixsize-1


for indexi = 1:matrixsize
    s = 0 ;
    for indexj = 1:matrixsize
        s = s + B(indexi,indexj)
    end;  //      for indexj = 1:matrixsize
    Bsorosszeg(indexi) = s ;
end;     // for indexi = 1:matrixsize-1
Bsorosszegnegyzetosszeg = 0 ;
for indexi = 1:matrixsize
    Bsorosszegnegyzetosszeg = Bsorosszegnegyzetosszeg + (Bsorosszeg(indexi))^2;
end;     // for indexi = 1:matrixsize-1
for indexi = 1:matrixsize
    s=0 ;
    for indexj = 1:matrixsize
        s = s + B(indexi,indexj) ; 
    end; 
    w(indexi) = s/((Bsorosszegnegyzetosszeg)^(1/2)) ; 
end;     // for indexi = 1:matrixsize-1

disp(w) ;

custom = 1;  // the Pareto optimality of a custom weight vector, given in line 36, is checked

format(16);

L = zeros(matrixsize,matrixsize);             

//w = [0.6514   0.3460 0.0879  0.3262  0.4773  0.3377 ];

I1elemszam = 0;
I0elemszam = 0;
for indexi = 1:matrixsize
  v(indexi) = log(w(indexi));
for indexj = 1:matrixsize
  b(indexi,indexj) = log(A(indexi,indexj));
  if real(w(indexi)/w(indexj)) - A(indexi,indexj) > 10^(-6) 
      I1elemszam = I1elemszam + 1 ;
      I1poziciok(I1elemszam,1) = indexi ;
      I1poziciok(I1elemszam,2) = indexj ;
  end;
end;  // for indexj = 1:matrixsize
end;  // for indexi = 1:matrixsize
for indexi = 1:matrixsize-1
for indexj = indexi+1:matrixsize
  if abs(w(indexi)/w(indexj) - A(indexi,indexj)) < 10^(-8)   
      I0elemszam = I0elemszam + 1  ;
      I0poziciok(I0elemszam,1) = indexi ;
      I0poziciok(I0elemszam,2) = indexj ;
      L(indexi,indexj)=-1 ;                   
      L(indexj,indexi)=-1 ;                   
  end;
end;  // for indexj = 1:matrixsize
end;  // for indexi = 1:matrixsize 


if I0elemszam > matrixsize-2                  
   if I0elemszam == matrixsize*(matrixsize-1)/2
       if custom == 0 
           disp('Pairwise comparison matrix A is consistent, its principal right eigenvector is efficient.');               
           text=strcat(['Pairwise comparison matrix A is consistent, its principal right eigenvector is efficient.']);     
           wout = w;
           custom = 2 ;
       else // if custom == 0    
           disp('Pairwise comparison matrix A is consistent, the weight vector you entered coincides with the principal right eigenvector and it is efficient.'); 
           text=strcat(['Pairwise comparison matrix A is consistent, the weight vector you entered coincides with the principal right eigenvector and it is efficient.']);  
           wout = w;
           custom = 2 ;
       end; // if custom == 0             
   else  
      DIAG = -sum(L,1);                          
      for indexi = 1:matrixsize                 
          L(indexi,indexi) = DIAG(indexi);      
      end;                                       
      L(:,matrixsize) =[];                      
      L(matrixsize,:) =[];                       
      if det(L) == 0 
      else
          disp('The weight vector you entered is efficient, because the set of edges (i,j) fulfilling a_ij = w_i/w_j forms or includes a spanning tree.');                
          text=strcat(['The weight vector you entered is efficient, because the set of edges (i,j) fulfilling a_ij = w_i/w_j forms or includes a spanning tree.']);
          wout = w;    
          custom = 2 ;
      end;                                       
end;   // if I0elemszam == matrixsize*(matrixsize-1)/2
end;   // if I0elemszam > matrixsize-2   
if custom ~= 2                                
    IneqSor = 0;     
    sIneqSor = 0;                                                              
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
        Ineq(IneqSor,matrixsize+I1index) =  1 ;      
        RHSineq(IneqSor) = v(indexi)-v(indexj) ;
    end;   // for I1index = 1:I1elemszam
    for indexi = 1:matrixsize
        c(indexi) = 0 ;
    end;    // for indexi = 1:matrixsize
    for indexi = matrixsize+1:matrixsize+I1elemszam
        c(indexi) = -1 ;
    end;    // for indexi = matrixsize+1:matrixsize+I1elemszam
    for indexi = 1:matrixsize
        lb(indexi) = -10 ;
    end;    // for indexi = 1:matrixsize
    for indexi = matrixsize+1:matrixsize+I1elemszam
        lb(indexi) = 0 ;
    end;    // for indexi = matrixsize+1:matrixsize+I1elemszam
    for indexi = 1:matrixsize+I1elemszam
        ub(indexi) = 30 ;
    end;    // for indexi = matrixsize+1:matrixsize+I1elemszam
    for indexi = 1:matrixsize
       start(indexi) = v(indexi);
    end;   // for indexi = 1:matrixsize-1
    for indexi = matrixsize+1:matrixsize+I1elemszam
       start(indexi) = 0;
    end;   // for indexi = matrixsize+1:matrixsize+I1elemszam
    for Eqindex = 1:matrixsize+I1elemszam
        Eq(1,Eqindex) = 0;
    end;   // for Eindex = matrixsize+1:matrixsize+I1elemszam

    Eq(1,1) = 1;
    RHSeq(1) = 0 ;

    for I0index = 1:I0elemszam
        indexi = I0poziciok(I0index,1);
        indexj = I0poziciok(I0index,2);    
        Eq(I0index+1,indexi) =  1 ;
        Eq(I0index+1,indexj) =  -1 ;  
        RHSeq(I0index+1) = -b(indexj,indexi) ;
    end;   // for I0index = 1:I0elemszam

    if I1elemszam > 0 
       EH = [Eq; Ineq] ;
       rhs = [RHSeq; RHSineq] ;
    else 
       EH = Eq ;
       rhs = RHSeq ;  
    end;
    
    if custom == 1 & I1elemszam == matrixsize*(matrixsize-1)/2                 
        sIneqSor = 0;                                                          
        for I1index = 1:I1elemszam                                             
            sIneqSor = sIneqSor + 1 ;                                          
            indexi = I1poziciok(I1index,1);                                    
            indexj = I1poziciok(I1index,2);                                    
            sIneq(sIneqSor,indexi) = -1 ;                                      
            sIneq(sIneqSor,indexj) =  1 ;                                      
            sRHSineq(sIneqSor) = -b(indexi,indexj) ;                           
            sIneqSor = sIneqSor + 1 ;                                          
            sIneq(sIneqSor,indexi) =   1 ;                                     
            sIneq(sIneqSor,indexj) =  -1 ;                                     
            sIneq(sIneqSor,matrixsize+1) =  1 ;                                
            sRHSineq(sIneqSor) = v(indexi)-v(indexj) ;                         
        end;   // for I1index = 1:I1elemszam                                   
        for indexi = 1:matrixsize                                              
            sc(indexi) = 0 ;                                                   
        end;    // for indexi = 1:matrixsize                                   
        sc(matrixsize+1) = -1 ;                                                
        for indexi = 1:matrixsize                                              
            slb(indexi) = -10 ;                                                
        end;    // for indexi = 1:matrixsize                                   
        slb(matrixsize+1) = 0 ;                                                
        for indexi = 1:matrixsize+1                                            
            sub(indexi) = 30 ;                                                 
        end;    // for indexi = matrixsize+1:matrixsize+1                      
        for indexi = 1:matrixsize                                              
            sstart(indexi) = v(indexi);                                        
        end;   // for indexi = 1:matrixsize                                    
        sstart(matrixsize+1) = 0;                                              
        for Eqindex = 1:matrixsize+1                                           
            sEq(1,Eqindex) = 0;                                                
        end;   // for Eindex = matrixsize+1:matrixsize+1                       
        sEq(1,1) = 1;                                                          
        sRHSeq(1) = 0 ;                                                        
        sEH = [sEq; sIneq] ;                                                   
        srhs = [sRHSeq; sRHSineq] ;                                            
        [sxopt,slagr,sfopt]=linpro(sc,sEH,srhs,slb,sub,1,sstart);              
        if sfopt <  -10^(-6) 
           R2I1elemszam = 0;
           R2I0elemszam = 0;
           swsum = 0 ;                
               for indexi = 1:matrixsize
                 sw(indexi) = exp(sxopt(indexi));  
                 swsum = swsum + sw(indexi) ;                
               end;  // for indexi = 1:matrixsize        
           for indexi = 1:matrixsize
             R2v(indexi) = log(sw(indexi));
           for indexj = 1:matrixsize
             R2b(indexi,indexj) = log(A(indexi,indexj));
              if real(sw(indexi)/sw(indexj)) - A(indexi,indexj) > 10^(-6) 
                   R2I1elemszam = R2I1elemszam + 1 ;
                   R2I1poziciok(R2I1elemszam,1) = indexi ;
                   R2I1poziciok(R2I1elemszam,2) = indexj ;
              end;
            end;  // for indexj = 1:matrixsize
            end;  // for indexi = 1:matrixsize
            for indexi = 1:matrixsize-1
            for indexj = indexi+1:matrixsize
                 if abs(sw(indexi)/sw(indexj) - A(indexi,indexj)) < 10^(-8)   
                     R2I0elemszam = R2I0elemszam + 1  ;
                     R2I0poziciok(R2I0elemszam,1) = indexi ;
                     R2I0poziciok(R2I0elemszam,2) = indexj ;
                 end;
            end;  // for indexj = 1:matrixsize
            end;  // for indexi = 1:matrixsize
            if R2I0elemszam < matrixsize*(matrixsize-1)/2
                R2IneqSor = 0;     
                for R2I1index = 1:R2I1elemszam
                       R2IneqSor = R2IneqSor + 1 ;
                       R2indexi = R2I1poziciok(R2I1index,1);
                       R2indexj = R2I1poziciok(R2I1index,2);    
                       R2Ineq(R2IneqSor,indexi) = -1 ;
                       R2Ineq(R2IneqSor,indexj) =  1 ;  
                       R2RHSineq(R2IneqSor) = -b(indexi,indexj) ;
                       R2IneqSor = R2IneqSor + 1 ;
                       R2Ineq(R2IneqSor,indexi) =   1 ;
                       R2Ineq(R2IneqSor,indexj) =  -1 ;  
                       R2Ineq(R2IneqSor,matrixsize+R2I1index) =  1 ;      
                       R2RHSineq(R2IneqSor) = R2v(indexi)-R2v(indexj) ;
                end;   // for I1index = 1:I1elemszam
                for indexi = 1:matrixsize
                       R2c(indexi) = 0 ;
                end;    // for indexi = 1:matrixsize
                for indexi = matrixsize+1:matrixsize+R2I1elemszam
                   R2c(indexi) = -1 ;
                end;    // for indexi = matrixsize+1:matrixsize+I1elemszam
                for indexi = 1:matrixsize
                       R2lb(indexi) = -10 ;
                end;    // for indexi = 1:matrixsize
                for indexi = matrixsize+1:matrixsize+R2I1elemszam
                       R2lb(indexi) = 0 ;
                end;    // for indexi = matrixsize+1:matrixsize+I1elemszam
                for indexi = 1:matrixsize+R2I1elemszam
                       R2ub(indexi) = 30 ;
                end;    // for indexi = matrixsize+1:matrixsize+I1elemszam
                for indexi = 1:matrixsize
                      R2start(indexi) = R2v(indexi);
                end;   // for indexi = 1:matrixsize-1
                for indexi = matrixsize+1:matrixsize+R2I1elemszam
                      R2start(indexi) = 0;
                end;   // for indexi = matrixsize+1:matrixsize+I1elemszam
                for R2Eqindex = 1:matrixsize+R2I1elemszam
                       R2Eq(1,R2Eqindex) = 0;
                end;   // for Eindex = matrixsize+1:matrixsize+I1elemszam
                   R2Eq(1,1) = 1;
                   R2RHSeq(1) = 0 ;
                   for R2I0index = 1:R2I0elemszam
                       indexi = R2I0poziciok(R2I0index,1);
                       indexj = R2I0poziciok(R2I0index,2);    
                       R2Eq(R2I0index+1,indexi) =  1 ;
                       R2Eq(R2I0index+1,indexj) =  -1 ;  
                       R2RHSeq(R2I0index+1) = -b(indexj,indexi) ;
                   end;   // for I0index = 1:I0elemszam
                if R2I1elemszam > 0 
                      R2EH = [R2Eq; R2Ineq] ;
                      R2rhs = [R2RHSeq; R2RHSineq] ;
                else 
                      R2EH = R2Eq ;
                      R2rhs = R2RHSeq ;  
                end;
                [R2xopt,R2lagr,R2fopt]=linpro(R2c,R2EH,R2rhs,R2lb,R2ub,1+R2I0elemszam,R2start); 
                R2wsum = 0 ; 
                for indexi = 1:matrixsize
                   R2w(indexi) = exp(R2xopt(indexi));
                   R2wsum = R2wsum + R2w(indexi) ;
                end;         // for indexi = 1:matrixsize
                for indexi = 1:matrixsize
                   R2w(indexi) = R2w(indexi)/R2wsum;
                end;        // for indexi = 1:matrixsize


                disp('The weight vector you entered is strongly inefficient. The efficient weight vector below dominates it strongly and internally.'); 
                text=strcat(['The weight vector you entered is strongly inefficient. The efficient weight vector below dominates it strongly and internally.']);     
// //                for indexi = 1:matrixsize
// //                   disp(strcat(["w",string(indexi),"=",string(R2w(indexi))]));
// //                       tablazat(indexi,1) = strcat(["w",string(indexi),"=",]);          
// //                       tablazat(indexi,2) = string(R2w(indexi)) ;                       
// //                    end;   // for indexi = 1:matrixsize
                disp(R2w) ;
                wout = R2w;
                custom = 2; // 


            else              // R2I0elemszam < matrixsize*(matrixsize-1)/2

                for indexi = 1:matrixsize
                   sw(indexi) = sw(indexi)/swsum;
                end;         // for indexi = 1:matrixsize

                disp('The weight vector you entered is strongly inefficient. The efficient weight vector below dominates it strongly and internally.');
                text=strcat(['The weight vector you entered is strongly inefficient. The efficient weight vector below dominates it strongly and internally.']);     
// //                for indexi = 1:matrixsize
// //                   disp(strcat(["w",string(indexi),"=",string(sw(indexi))]));
// //                   tablazat(indexi,1) = strcat(["w",string(indexi),"=",]);          
// //                   tablazat(indexi,2) = string(sw(indexi)) ;                       
// //                end;   // for indexi = 1:matrixsize
                disp(sw) ;
                wout = sw;
                custom = 2; // 
            end;  // R2I0elemszam < matrixsize*(matrixsize-1)/2

        end; // if sfopt <  -10^(-6) 
    end;  // custom == 1 & I1elemszam == matrixsize*(matrixsize-1)/2 



    if custom ~= 2
       [xopt,lagr,fopt]=linpro(c,EH,rhs,lb,ub,1+I0elemszam,start); 
       
       if custom == 1 & fopt >= -10^(-6)  
              disp('The weight vector you entered is efficient. ');
              text=strcat(['The weight vector you entered is efficient. ']);     
              wout = w;
       end;  // if custom == 1 & fopt >= -10^(-6)
  
       if custom == 0 & fopt >= -10^(-6)  
              disp('The principal right eigenvector is efficient.');
              text=strcat(['The principal right eigenvector is efficient.']);     
              wout = w;
       end;  // if custom == 0 & fopt >= -10^(-6)  

       if custom == 0 & fopt < -10^(-6)  
            dwsum = 0 ; 
             for indexi = 1:matrixsize
                dw(indexi) = exp(xopt(indexi));
                dwsum = dwsum + dw(indexi) ;
             end;   // for indexi = 1:matrixsize
             for indexi = 1:matrixsize
                dw(indexi) = dw(indexi)/dwsum;
             end;   // for indexi = 1:matrixsize

             disp('The principal right eigenvector is inefficient (but weakly efficient). The efficient weight vector below dominates it internally.'); 
             text=strcat(['The principal right eigenvector is inefficient (but weakly efficient). The efficient weight vector below dominates it internally.']);
// //             for indexi = 1:matrixsize
// //              disp(strcat(["w",string(indexi),"=",string(dw(indexi))]));
// //             tablazat(indexi,1) = strcat(["w",string(indexi),"=",]) ;               
// //             tablazat(indexi,2) = string(dw(indexi)) ;                              
// //             end;   // for indexi = 1:matrixsize
             disp(dw) ;
             wout = dw;

       end;  // custom == 0 & fopt < -10^(-6)

       if custom == 1 & fopt < -10^(-6)  
          dwsum = 0 ; 
          for indexi = 1:matrixsize
           dw(indexi) = exp(xopt(indexi));
           dwsum = dwsum + dw(indexi) ;
          end;   // for indexi = 1:matrixsize
          for indexi = 1:matrixsize
           dw(indexi) = dw(indexi)/dwsum;
          end;   // for indexi = 1:matrixsize

          disp('The weight vector you entered is inefficient (but weakly efficient). The efficient weight vector below dominates it internally.');
          text=strcat(['The weight vector you entered is inefficient (but weakly efficient). The efficient weight vector below dominates it internally.']);     
// //          for indexi = 1:matrixsize
// //              disp(strcat(["w",string(indexi),"=",string(dw(indexi))]));
// //              tablazat(indexi,1) = strcat(["w",string(indexi),"=",]) ;               
// //              tablazat(indexi,2) = string(dw(indexi)) ;                              
// //          end;   // for indexi = 1:matrixsize
          disp(dw) ;
          wout = dw;
       end; // custom == 1 & fopt < -10^(-6) 

    end;  // if custom ~= 2 
 
end; // if custom ~= 2   
 

