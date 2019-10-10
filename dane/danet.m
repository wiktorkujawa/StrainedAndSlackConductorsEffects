
load danet.txt
 #figure('renderer','painters','pos',[600 300 1200 700])
 
 minimum=min(danet(:,2));
 maksimum=max(danet(:,2));
 
 for i=3:5
   if(minimum>min(danet(:,i)))
   minimum=min(danet(:,i));
   endif
   
   if(maksimum<max(danet(:,i)))
   maksimum=max(danet(:,i));
   endif
 endfor
   
   maksimum=int32(maksimum);
   minimum=int32(minimum);
x=danet(:,1);

l=length(danet(:,1));

subplot(2,2,1);
 plot(x,danet(:,5));
 xlabel(['\bft[s]']);
 ylabel(['\bfF[kN]']);
 title('Overal Maximal tensile force');
grid on


subplot(2,2,2);
 plot(x,danet(:,2));
  xlabel(['\bft[s]']);
 ylabel(['\bfF[kN]']);
 title('Tensile force caused by swing-out');
grid on

axis([int32(danet(1,1)) int32(danet(l,1)) minimum maksimum])

 subplot(2,2,3);
 plot(x,danet(:,3));
 xlabel(['\bft[s]']);
 ylabel(['\bfF[kN]']);
 title('Tensile force caused by drop');
  grid on
  
 axis([int32(danet(1,1)) int32(danet(l,1)) minimum maksimum])
  
 subplot(2,2,4);
 plot(x,danet(:,4));
 xlabel(['\bft[s]']);
 ylabel(['\bfF[kN]']);
 title('Tensile force caused by pinch');
grid on


axis([int32(danet(1,1)) int32(danet(l,1)) minimum maksimum])

figure
 plot(x,danet(:,2),x,danet(:,3),x,danet(:,4));
 xlabel(['\bft[s]']);
 ylabel(['\bfF[kN]']);
 title('All forces on one chart');
grid on
legend('Swing-out','Drop force','Pinch force',4);

figure
subplot(1,2,1)
plot(x,danet(:,6));
title('Horizontal span displacement');
xlabel(['\bft[s]']);
 ylabel(['\bfbmax[m]']);
 grid on
 
 subplot(1,2,2)
plot(x,danet(:,7));
title('Minimum air clearance');
xlabel(['\bft[s]']);
 ylabel(['\bfamin[m]']);
 grid on
