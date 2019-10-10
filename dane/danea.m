
load danea.txt
 #figure('renderer','painters','pos',[600 300 1200 700])
 
 minimum=min(danea(:,2));
 maksimum=max(danea(:,2));
 
 for i=3:5
   if(minimum>min(danea(:,i)))
   minimum=min(danea(:,i));
   endif
   
   if(maksimum<max(danea(:,i)))
   maksimum=max(danea(:,i));
   endif
 endfor
   
   maksimum=int32(maksimum);
   minimum=int32(minimum);
x=danea(:,1);

l=length(danea(:,1));

subplot(2,2,1);
 plot(x,danea(:,5));
 xlabel(['\bfa[m]']);
 ylabel(['\bfF[kN]']);
 title('Overal Maximal tensile force');
grid on


subplot(2,2,2);
 plot(x,danea(:,2));
  xlabel(['\bfa[m]']);
 ylabel(['\bfF[kN]']);
 title('Tensile force caused by swing-out');
grid on

axis([int32(danea(1,1)) int32(danea(l,1)) minimum maksimum])

 subplot(2,2,3);
 plot(x,danea(:,3));
 xlabel(['\bfa[m]']);
 ylabel(['\bfF[kN]']);
 title('Tensile force caused by drop');
  grid on
  
 axis([int32(danea(1,1)) int32(danea(l,1)) minimum maksimum])
  
 subplot(2,2,4);
 plot(x,danea(:,4));
 xlabel(['\bfa[m]']);
 ylabel(['\bfF[kN]']);
 title('Tensile force caused by pinch');
grid on


axis([int32(danea(1,1)) int32(danea(l,1)) minimum maksimum])

figure
 plot(x,danea(:,2),x,danea(:,3),x,danea(:,4));
 xlabel(['\bfa[m]']);
 ylabel(['\bfF[kN]']);
 title('All forces on one chart');
grid on
legend('Swing-out','Drop force','Pinch force',4);

figure
subplot(1,2,1)
plot(x,danea(:,6));
title('Horizontal span diplacement');
xlabel(['\bfa[m]']);
 ylabel(['\bfbmax[m]']);
 grid on
 
 subplot(1,2,2)
plot(x,danea(:,7));
title('Minimum air clearance');
xlabel(['\bfa[m]']);
 ylabel(['\bfamin[m]']);
 grid on
 
 