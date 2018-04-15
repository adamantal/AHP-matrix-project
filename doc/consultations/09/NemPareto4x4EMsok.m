
tablazat =  xlsread('NemPareto4x4EM.xls');   

disp('');
fid = fopen('NemPareto4x4EM.tex','w');
fprintf(fid,'\\documentclass{article} \n');  
fprintf(fid,'\\usepackage[latin1]{inputenc} \n');  
fprintf(fid,'\\usepackage{amssymb, amsmath, amsthm} \n');  
fprintf(fid,'\\usepackage{color} \n');  
fprintf(fid,'\\definecolor{gr}{rgb}{0,0.5,0} \n');  
fprintf(fid,'\\begin{document} \n');  
fprintf(fid,'\\newtheorem{example}{Example}[section] \n');  
fprintf(fid,'\\section*{Appendix} \n');      
for matrixindex=1:591
vk = 0 ;
tk = 0 ;
a12 = tablazat(matrixindex,2) ;  a21 = 1/a12 ;
a13 = tablazat(matrixindex,3) ;  a31 = 1/a13 ;
a14 = tablazat(matrixindex,4) ;  a41 = 1/a14 ;
a23 = tablazat(matrixindex,5) ;  a32 = 1/a23 ;
a24 = tablazat(matrixindex,6) ;  a42 = 1/a24 ;
a34 = tablazat(matrixindex,7) ;  a43 = 1/a34 ;

A = [    1            a12        a13        a14    ; ...
        1/a12          1         a23        a24    ; ...
        1/a13        1/a23        1         a34    ; ...
        1/a14        1/a24      1/a34        1     ] ;

wEM1 = tablazat(matrixindex,11) ;
wEM2 = tablazat(matrixindex,12) ;
wEM3 = tablazat(matrixindex,13) ;
wEM4 = tablazat(matrixindex,14) ;

wEMsum = wEM1 + wEM2 + wEM3 + wEM4 ;
wEM1n = wEM1/wEMsum ; 
wEM2n = wEM2/wEMsum ; 
wEM3n = wEM3/wEMsum ; 
wEM4n = wEM4/wEMsum ; 

dw1 = tablazat(matrixindex,16) ;
dw2 = tablazat(matrixindex,17) ;
dw3 = tablazat(matrixindex,18) ;
dw4 = tablazat(matrixindex,19) ;

dwsum = dw1 + dw2 + dw3 + dw4 ;
dw1n = dw1/dwsum ; 
dw2n = dw2/dwsum ; 
dw3n = dw3/dwsum ; 
dw4n = dw4/dwsum ; 

XEM = [    1       wEM1/wEM2  wEM1/wEM3  wEM1/wEM4 ; ...
       wEM2/wEM1       1      wEM2/wEM3  wEM2/wEM4 ; ...
       wEM3/wEM1   wEM3/wEM2      1      wEM3/wEM4 ; ...
       wEM4/wEM1   wEM4/wEM2  wEM4/wEM3      1     ] ;

Xdw = [    1        dw1/dw2    dw1/dw3     dw1/dw4 ; ...
        dw2/dw1        1       dw2/dw3     dw2/dw4 ; ...
        dw3/dw1     dw3/dw2       1        dw3/dw4 ; ...
        dw4/dw1     dw4/dw2    dw4/dw3        1     ] ;

Elteresek = abs(XEM-Xdw);
Elteresek = Elteresek + Elteresek';

for indexi=1:4
    elteresekszama = 0 ;    
    for indexj=1:4
       if Elteresek(indexi,indexj) > 10^(-4)
          elteresekszama = elteresekszama + 1 ; 
       end;   % if abs(XEM(indexi,indexj)-Xdw(indexi,indexj)) > 10^(-5)
    end;   %    for indexj=1:4
    if elteresekszama > 2
       vk = indexi ; 
       disp([elteresekszama vk]);
   end;  %  if elteresekszama > 2   
end;   %    for indexi=1:4    

XdwAElteresek = abs(A-Xdw);

XdwAElteresek(vk,vk) = 10;
for indexj=1:4
    if XdwAElteresek(vk,indexj) < 10^(-4)
       tk = indexj ; 
    end;   % if Elteresek(vk,indexj) < 10^(-6)
end;   %    for indexj=1:4
        
fprintf(fid,'\\begin{example} \n');  
fprintf(fid,'\\begin{equation*} \n');  
fprintf(fid,'\\mathbf{A} = \n');  
fprintf(fid,'\\begin{pmatrix} \n');  
fprintf(fid,'$\\,\\,$');  
fprintf(fid,' 1 ');  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 

if a12 > 0.9  
fprintf(fid,'%1.0f',a12);   % a_{12}
else
fprintf(fid,' 1/');      
fprintf(fid,'%1.0f',round(1/a12));     
end;
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 

if a13 > 0.9  
fprintf(fid,'%1.0f',a13);   % a_{13}
else
fprintf(fid,' 1/');      
fprintf(fid,'%1.0f',round(1/a13));     
end;
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 

if a14 > 0.9  
fprintf(fid,'%1.0f',a14);   % a_{14}
else
fprintf(fid,' 1/');      
fprintf(fid,'%1.0f',round(1/a14));      
end;
fprintf(fid,' $\\,\\,$ \\\\ \n');   

fprintf(fid,'$\\,\\,$');  
if a12 > 0.9   % a_{21}
   if a12 == 1    
   fprintf(fid,' 1 ');   
   else
   fprintf(fid,' 1/');      
   fprintf(fid,'%1.0f',round(a12));     
   end;
else
fprintf(fid,'%1.0f',round(1/a12));   
end;

fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
fprintf(fid,' 1 ');   % 1
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 

if a23 > 0.9    % a_{23}
fprintf(fid,'%1.0f',a23);   
else
fprintf(fid,' 1/');      
fprintf(fid,'%1.0f',round(1/a23));     
end;

fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 

if a24 > 0.9    % a_{24}
fprintf(fid,'%1.0f',a24);   
else
fprintf(fid,' 1/');      
fprintf(fid,'%1.0f',round(1/a24));     
end;

fprintf(fid,' $\\,\\,$ \\\\ \n');   

fprintf(fid,'$\\,\\,$');  

if a13 > 0.9   % a_{31}
   if a13 == 1    
   fprintf(fid,' 1 ');   
   else
   fprintf(fid,' 1/');      
   fprintf(fid,'%1.0f',round(a13));     
   end;
else
fprintf(fid,'%1.0f',round(1/a13));   
end;

fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 

if a23 > 0.9   % a_{32}
   if a23 == 1    
   fprintf(fid,' 1 ');   
   else
   fprintf(fid,' 1/');      
   fprintf(fid,'%1.0f',round(a23));     
   end;
else
fprintf(fid,'%1.0f',round(1/a23));   
end;

fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
fprintf(fid,' 1 ');   % 1
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 

if a34 > 0.9  % a_{34}
fprintf(fid,'%1.0f',a34);   
else
fprintf(fid,' 1/');      
fprintf(fid,'%1.0f',round(1/a34));      
end;

fprintf(fid,' $\\,\\,$ \\\\ \n');   

fprintf(fid,'$\\,\\,$');  

if a14 > 0.9   % a_{41}
   if a14 == 1    
   fprintf(fid,' 1 ');   
   else
   fprintf(fid,' 1/');      
   fprintf(fid,'%1.0f',round(a14));     
   end;
else
fprintf(fid,'%1.0f',round(1/a14));   
end;


fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 

if a24 > 0.9   % a_{42}
   if a24 == 1    
   fprintf(fid,' 1 ');   
   else
   fprintf(fid,' 1/');      
   fprintf(fid,'%1.0f',round(a24));     
   end;
else
fprintf(fid,'%1.0f',round(1/a24));   
end;

fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 

if a34 > 0.9   % a_{34}
   if a34 == 1    
   fprintf(fid,' 1 ');   
   else
   fprintf(fid,' 1/');      
   fprintf(fid,'%1.0f',round(a34));     
   end;
else
fprintf(fid,'%1.0f',round(1/a34));   
end;

fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
fprintf(fid,' 1  $\\,\\,$ \\\\ \n');   % 1 

fprintf(fid,'\\end{pmatrix}, \n');  
fprintf(fid,'\\qquad \n');  
fprintf(fid,'\\lambda_{\\max} =  \n');  
fprintf(fid,'%6.4f',tablazat(matrixindex,8));  
fprintf(fid,', \n');  
fprintf(fid,'\\qquad \n');  
fprintf(fid,'CR = ');  
fprintf(fid,'%6.4f',tablazat(matrixindex,9));  
fprintf(fid,' \n');  
fprintf(fid,'\\end{equation*} \n');  
fprintf(fid,'\n'); 
fprintf(fid,'\\begin{equation*} \n');  
fprintf(fid,'\\mathbf{w}^{EM} = \n');  
fprintf(fid,'\\begin{pmatrix} \n');  
if vk == 1 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%8.6f',wEM1n);  
if vk == 1 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'\\\\ \n');  
if vk == 2 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%8.6f',wEM2n);  
if vk == 2 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'\\\\ \n');  
if vk == 3 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%8.6f',wEM3n);  
if vk == 3 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'\\\\ \n');  
if vk == 4 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%8.6f',wEM4n);  
if vk == 4 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,' \n');  
fprintf(fid,'\\end{pmatrix}');  
fprintf(fid,'\\end{equation*} \n');  
fprintf(fid,'\\begin{equation*} \n');  
fprintf(fid,'\\left[ \\frac{{w}^{EM}_i}{{w}^{EM}_j} \\right] = \n');  
fprintf(fid,'\\begin{pmatrix} \n');  
fprintf(fid,'$\\,\\,$');  
fprintf(fid,' 1 ');  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 1 || vk == 2 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%6.4f',wEM1/wEM2);   % a_{12}
if vk == 1 || vk == 2 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 1 || vk == 3 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%6.4f',wEM1/wEM3);   % a_{13}
if vk == 1 || vk == 3 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 1 || vk == 4 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%6.4f',wEM1/wEM4);   % a_{14}
if vk == 1 || vk == 4 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ \\\\ \n'); 

fprintf(fid,'$\\,\\,$');  
if vk == 2 || vk == 1 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%6.4f',wEM2/wEM1);   % a_{21}
if vk == 2 || vk == 1 fprintf(fid,'\\color{black} '); end;  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
fprintf(fid,' 1 ');   % 1
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 2 || vk == 3 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%6.4f',wEM2/wEM3);   % a_{23}
if vk == 2 || vk == 3 fprintf(fid,'\\color{black} '); end;  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 2 || vk == 4 fprintf(fid,'\\color{red} '); end;    
fprintf(fid,'%6.4f',wEM2/wEM4);   % a_{24}
if vk == 2 || vk == 4 fprintf(fid,'\\color{black} '); end;  
fprintf(fid,'  $\\,\\,$ \\\\ \n');   

fprintf(fid,'$\\,\\,$');  
if vk == 3 || vk == 1 fprintf(fid,'\\color{red} '); end;  
fprintf(fid,'%6.4f',wEM3/wEM1);  % a_{31}
if vk == 3 || vk == 1 fprintf(fid,'\\color{black} '); end;  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 3 || vk == 2 fprintf(fid,'\\color{red} '); end;  
fprintf(fid,'%6.4f',wEM3/wEM2);   % a_{32}
if vk == 3 || vk == 2 fprintf(fid,'\\color{black} '); end;  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
fprintf(fid,' 1 ');   % 1
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 3 || vk == 4 fprintf(fid,'\\color{red} '); end;  
fprintf(fid,'%6.4f',wEM3/wEM4);   % a_{34}
if vk == 3 || vk == 4 fprintf(fid,'\\color{black} '); end;  
fprintf(fid,' $\\,\\,$ \\\\ \n');   

fprintf(fid,'$\\,\\,$');  
if vk == 4 || vk == 1 fprintf(fid,'\\color{red} '); end;  
fprintf(fid,'%6.4f',wEM4/wEM1);  % a_{41}
if vk == 4 || vk == 1 fprintf(fid,'\\color{black} '); end;  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 4 || vk == 2 fprintf(fid,'\\color{red} '); end;  
fprintf(fid,'%6.4f',wEM4/wEM2);   % a_{42}
if vk == 4 || vk == 2 fprintf(fid,'\\color{black} '); end;  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 4 || vk == 3 fprintf(fid,'\\color{red} '); end;  
fprintf(fid,'%6.4f',wEM4/wEM3);   % a_{43}
if vk == 4 || vk == 3 fprintf(fid,'\\color{black} '); end;  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
fprintf(fid,' 1  $\\,\\,$ \\\\ \n');   % 1 

fprintf(fid,'\\end{pmatrix}, \n');  
fprintf(fid,'\\end{equation*} \n');  
fprintf(fid,'\n'); 
fprintf(fid,'\\begin{equation*} \n');  
fprintf(fid,'\\mathbf{w}^{\\prime} = \n');  
fprintf(fid,'\\begin{pmatrix} \n');  
fprintf(fid,'%8.6f',dw1n);  
fprintf(fid,'\\\\ \n');  
fprintf(fid,'%8.6f',dw2n);  
fprintf(fid,'\\\\ \n');  
fprintf(fid,'%8.6f',dw3n);  
fprintf(fid,'\\\\ \n');  
fprintf(fid,'%8.6f',dw4n);  
fprintf(fid,' \n');  
fprintf(fid,'\\end{pmatrix} = \n');  
if vk == 1 
fprintf(fid,'%8.6f',dw2n/wEM2n);  
else
fprintf(fid,'%8.6f',dw1n/wEM1n);      
end;    
fprintf(fid,'\\cdot \n');  
fprintf(fid,'\\begin{pmatrix} \n');  
if vk == 1 
    fprintf(fid,'\\color{gr} '); 
    fprintf(fid,'%8.6f',dw1n*wEM2n/dw2n);  
else
   fprintf(fid,'%8.6f',wEM1n);  
end;    
if vk == 1 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'\\\\ \n');  
if vk == 2 
    fprintf(fid,'\\color{gr} '); 
    fprintf(fid,'%8.6f',dw2n*wEM1n/dw1n);      
else    
    fprintf(fid,'%8.6f',wEM2n);  
end;        
if vk == 2 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'\\\\ \n');  
if vk == 3 
    fprintf(fid,'\\color{gr} '); 
    fprintf(fid,'%8.6f',dw3n*wEM1n/dw1n);          
else
    fprintf(fid,'%8.6f',wEM3n);  
end;    
if vk == 3 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'\\\\ \n');  
if vk == 4 
    fprintf(fid,'\\color{gr} '); 
    fprintf(fid,'%8.6f',dw4n*wEM1n/dw1n);          
else    
    fprintf(fid,'%8.6f',wEM4n);  
end;    
if vk == 4 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,' \n');  
fprintf(fid,'\\end{pmatrix}, \n');  
fprintf(fid,'\\end{equation*} \n');  
fprintf(fid,'\\begin{equation*} \n');  
fprintf(fid,'\\left[ \\frac{{w}^{\\prime}_i}{{w}^{\\prime}_j} \\right] = \n');  
fprintf(fid,'\\begin{pmatrix} \n');  
fprintf(fid,'$\\,\\,$');  
fprintf(fid,' 1 ');  
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 1 || vk == 2 fprintf(fid,'\\color{gr} '); end;    

if (vk == 1 && tk == 2) || (vk == 2 && tk == 1) 
    fprintf(fid,'\\color{blue} '); 
    if a12 > 0.9   % a_{12}
       fprintf(fid,'%1.0f',round(a12));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a21));     
    end;
else % if (vk == 1 && tk == 2) || (vk == 2 && tk == 1) 
  fprintf(fid,'%6.4f',dw1/dw2);  % a_{12}
end; % if (vk == 1 && tk == 2) || (vk == 2 && tk == 1) 

if vk == 1 || vk == 2 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 1 || vk == 3 fprintf(fid,'\\color{gr} '); end;    

if (vk == 1 && tk == 3) || (vk == 3 && tk == 1) 
    fprintf(fid,'\\color{blue} '); 
    if a13 > 0.9   % a_{13}
       fprintf(fid,'%1.0f',round(a13));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a31));     
    end;
else % if (vk == 1 && tk == 3) || (vk == 3 && tk == 1) 
  fprintf(fid,'%6.4f',dw1/dw3);  % a_{13}
end; % if (vk == 1 && tk == 3) || (vk == 3 && tk == 1) 

if vk == 1 || vk == 3 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 1 || vk == 4 fprintf(fid,'\\color{gr} '); end;    

if (vk == 1 && tk == 4) || (vk == 4 && tk == 1) 
    fprintf(fid,'\\color{blue} '); 
    if a14 > 0.9   % a_{14}
       fprintf(fid,'%1.0f',round(a14));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a41));     
    end;
else % if (vk == 1 && tk == 4) || (vk == 4 && tk == 1) 
  fprintf(fid,'%6.4f',dw1/dw4);  % a_{14}
end; % if (vk == 1 && tk == 4) || (vk == 4 && tk == 1) 

if vk == 1 || vk == 4 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ \\\\ \n'); 

fprintf(fid,'$\\,\\,$');  
if vk == 2 || vk == 1 fprintf(fid,'\\color{gr} '); end;    

if (vk == 2 && tk == 1) || (vk == 1 && tk == 2) 
    fprintf(fid,'\\color{blue} '); 
    if a21 > 0.9   % a_{21}
       fprintf(fid,'%1.0f',round(a21));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a12));     
    end;
else % if (vk == 2 && tk == 1) || (vk == 1 && tk == 2) 
  fprintf(fid,'%6.4f',dw2/dw1);  % a_{21}
end; % if (vk == 2 && tk == 1) || (vk == 1 && tk == 2) 

if vk == 2 || vk == 1 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
fprintf(fid,' 1 ');   % 1
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 2 || vk == 3 fprintf(fid,'\\color{gr} '); end;    

if (vk == 2 && tk == 3) || (vk == 3 && tk == 2) 
    fprintf(fid,'\\color{blue} '); 
    if a23 > 0.9   % a_{23}
       fprintf(fid,'%1.0f',round(a23));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a32));     
    end;
else % if (vk == 2 && tk == 3) || (vk == 3 && tk == 2) 
  fprintf(fid,'%6.4f',dw2/dw3);  % a_{23}
end; % if (vk == 2 && tk == 3) || (vk == 3 && tk == 2) 

if vk == 2 || vk == 3 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 2 || vk == 4 fprintf(fid,'\\color{gr} '); end;    

if (vk == 2 && tk == 4) || (vk == 4 && tk == 2) 
    fprintf(fid,'\\color{blue} '); 
    if a24 > 0.9   % a_{24}
       fprintf(fid,'%1.0f',round(a24));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a42));     
    end;
else % if (vk == 2 && tk == 4) || (vk == 4 && tk == 2) 
  fprintf(fid,'%6.4f',dw2/dw4);  % a_{23}
end; % if (vk == 2 && tk == 4) || (vk == 4 && tk == 2) 

if vk == 2 || vk == 4 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'  $\\,\\,$ \\\\ \n');   

fprintf(fid,'$\\,\\,$');  
if vk == 3 || vk == 1 fprintf(fid,'\\color{gr} '); end;    

if (vk == 3 && tk == 1) || (vk == 1 && tk == 3) 
    fprintf(fid,'\\color{blue} '); 
    if a31 > 0.9   % a_{31}
       fprintf(fid,'%1.0f',round(a31));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a13));     
    end;
else % if (vk == 3 && tk == 1) || (vk == 1 && tk == 3) 
  fprintf(fid,'%6.4f',dw3/dw1);  % a_{31}
end; % if (vk == 3 && tk == 1) || (vk == 1 && tk == 3) 

if vk == 3 || vk == 1 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 3 || vk == 2 fprintf(fid,'\\color{gr} '); end;    

if (vk == 3 && tk == 2) || (vk == 2 && tk == 3) 
    fprintf(fid,'\\color{blue} '); 
    if a32 > 0.9   % a_{32}
       fprintf(fid,'%1.0f',round(a32));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a23));     
    end;
else % if (vk == 3 && tk == 2) || (vk == 2 && tk == 3) 
  fprintf(fid,'%6.4f',dw3/dw2);  % a_{32}
end; % if (vk == 3 && tk == 2) || (vk == 2 && tk == 3) 

if vk == 3 || vk == 2 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
fprintf(fid,' 1 ');   % 1
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 3 || vk == 4 fprintf(fid,'\\color{gr} '); end;    

if (vk == 3 && tk == 4) || (vk == 4 && tk == 3) 
    fprintf(fid,'\\color{blue} '); 
    if a34 > 0.9   % a_{34}
       fprintf(fid,'%1.0f',round(a34));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a43));     
    end;
else % if (vk == 3 && tk == 4) || (vk == 4 && tk == 3) 
  fprintf(fid,'%6.4f',dw3/dw4);  % a_{34}
end; % if (vk == 3 && tk == 4) || (vk == 4 && tk == 3) 
if vk == 3 || vk == 4 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,' $\\,\\,$ \\\\ \n');   

fprintf(fid,'$\\,\\,$');  
if vk == 4 || vk == 1 fprintf(fid,'\\color{gr} '); end;    
if (vk == 4 && tk == 1) || (vk == 1 && tk == 4) 
    fprintf(fid,'\\color{blue} '); 
    if a41 > 0.9   % a_{41}
       fprintf(fid,'%1.0f',round(a41));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a14));     
    end;
else % if (vk == 4 && tk == 1) || (vk == 1 && tk == 4) 
  fprintf(fid,'%6.4f',dw4/dw1);  % a_{41}
end; % if (vk == 4 && tk == 1) || (vk == 1 && tk == 4)    

if vk == 4 || vk == 1 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 4 || vk == 2 fprintf(fid,'\\color{gr} '); end;    

if (vk == 4 && tk == 2) || (vk == 2 && tk == 4)
    fprintf(fid,'\\color{blue} '); 
    if a42 > 0.9   % a_{42}
       fprintf(fid,'%1.0f',round(a42));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a24));     
    end;
else % if (vk == 4 && tk == 2) || (vk == 2 && tk == 4)
  fprintf(fid,'%6.4f',dw4/dw2);   % a_{42}
end; % if (vk == 4 && tk == 2) || (vk == 2 && tk == 4)   

if vk == 4 || vk == 2 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
if vk == 4 || vk == 3 fprintf(fid,'\\color{gr} '); end;    

if (vk == 4 && tk == 3) || (vk == 3 && tk == 4) 
    fprintf(fid,'\\color{blue} '); 
    if a43 > 0.9   % a_{43}
       fprintf(fid,'%1.0f',round(a43));           
    else
       fprintf(fid,' 1/');      
       fprintf(fid,'%1.0f',round(a34));     
    end;
else % if if (vk == 4 && tk == 3) || (vk == 3 && tk == 4) 
  fprintf(fid,'%6.4f',dw4/dw3);   % a_{43}
end;  % if if (vk == 4 && tk == 3) || (vk == 3 && tk == 4)    

if vk == 4 || vk == 3 fprintf(fid,'\\color{black} '); end;    
fprintf(fid,'$\\,\\,$ & $\\,\\,$'); 
fprintf(fid,' 1  $\\,\\,$ \\\\ \n');   % 1 

fprintf(fid,'\\end{pmatrix}, \n');  fprintf(fid,'\\end{equation*} \n');  
fprintf(fid,'\\end{example} \n');  
fprintf(fid,'\\newpage \n');  
end % for matrixindex=1:1

fprintf(fid,'\\end{document} \n');  
fclose(fid);