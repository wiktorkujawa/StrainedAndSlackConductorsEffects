
load danei.txt
 #figure('renderer','painters','pos',[600 300 1200 700])
 
 minimum=min(danei(:,2));
 maksimum=max(danei(:,2));
 
 for i=3:5
   if(minimum>min(danei(:,i)))
   minimum=min(danei(:,i));
   endif
   
   if(maksimum<max(danei(:,i)))
   maksimum=max(danei(:,i));
   endif
 endfor
   
   maksimum=int32(maksimum);
   minimum=int32(minimum);
x=danei(:,1);

l=length(danei(:,1));

subplot(2,2,1);
 plot(x,danei(:,5));
 xlabel(['\bfI[kA]']);
 ylabel(['\bfF[kN]']);
 title('Overal Maximal tensile force');
grid on


subplot(2,2,2);
 plot(x,danei(:,2));
  xlabel(['\bfI[kA]']);
 ylabel(['\bfF[kN]']);
 title('Tensile force caused by swing-out');
grid on

axis([int32(danei(1,1)) int32(danei(l,1)) minimum maksimum])

 subplot(2,2,3);
 plot(x,danei(:,3));
 xlabel(['\bfI[kA]']);
 ylabel(['\bfF[kN]']);
 title('Tensile force caused by drop');
  grid on
  
 axis([int32(danei(1,1)) int32(danei(l,1)) minimum maksimum])
  
 subplot(2,2,4);
 plot(x,danei(:,4));
 xlabel(['\bfI[kA]']);
 ylabel(['\bfF[kN]']);
 title('Tensile force caused by pinch');
grid on


axis([int32(danei(1,1)) int32(danei(l,1)) minimum maksimum])

figure
 plot(x,danei(:,2),x,danei(:,3),x,danei(:,4));
 xlabel(['\bfI[kA]']);
 ylabel(['\bfF[kN]']);
 title('All forces on one chart');
grid on
legend('Swing-out','Drop force','Pinch force',4);

figure
subplot(1,2,1)
plot(x,danei(:,6));
title('Horizontal span diplacement');
xlabel(['\bfI[kA]']);
 ylabel(['\bfbmax[m]']);
 grid on
 
 subplot(1,2,2)
plot(x,danei(:,7));
title('Minimum air clearance');
xlabel(['\bfI[kA]']);
 ylabel(['\bfamin[m]']);
 grid on
