#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<process.h>
#include<windows.h>
#define oneAndAHalfDivByTenDivByStandardGravity 0.015290519877
#define twoPiXF 314.159265359
int main()
{
    puts("Choose a parameter in terms of calculations will be performed:\n *'a'-center line distance between main conductor mid-points; \n *'i' - three phase initial symmetrical short-circuit current'; \n *'t'-duration of short-circuit current");
    //Mechanical data, temperatures, forces, Young module, spring constant etc.
    double staticTensileForceOfConductorAtLowerTemperature,staticTensileForceOfConductorAtUpperTemperature,E,S,Tk,r,fi1,lowerTemperature,upperTemperature,Timin,Timax,Trmin,Trmax,Emin,Emax,stresmin,stresmax;
    //sag(static)
    double fmin,fmax;
    double i,l,li,lc,SubConductorCrossSection,num,amax,amin,da,a;
    //quantities describing conductor
    double massPerUnitLengthOfOneSubConductor,massOfOneConnection,massOfOneSpacer,d,m,ac,horizontalSpanDisplacement;
    int it,numberOfSubConductors,numberOfSpacers;
    double fiendmin,fiendmax;
    double Xmin,Xmax;
    double maximumSwingOutAngleAtLowerTemperature, maximumSwingOutAngleAtUpperTemperature;
    char wybor;
    double parmin,parmax;
    double row;
    int lineSystem;
    double loadParameterAtLowerTemperature;
    double loadParameterAtUpperTemperature;
    //Searched forces - results
    double maximumShortCircuitTensileForce,shortCircuitTensileForceAtLowerTemperature,shortCircuitTensileForceAtUpperTemperature,dropForceAtLowerTemperature,dropForceAtUpperTemperature,maximumDropForce;
    double elasticExpansionAtLowerTemperature,elasticExpansionAtUpperTemperature,thermalExpansionAtLowerTemperature,thermalExpansionAtUpperTemperature;
    double Nmin,Nmax,k;
    int dropper;
    char conductorMaterialType;
    //dynamic sag
    double materialConstant,dilatationFactorAtLowerTemperature,dilatationFactorAtUpperTemperature,formFactor,fedmin,fedmax;
    double deltamin,deltamax;
    //dropper parametres
    double heightOfTheDropperAtUpperTemperture,widthOfTheDropperAtUpperTemperture;
    double dropperCordLength=0.0;
    char uklad;
    char conductorType;
    double bhmin,bhmax;
    double ls,lk;
    //pinch force parameters
    double pinchForceAtLowerTemperature,pinchForceAtUpperTemperature,maximumPinchForce,shortCircuitCurrentBundleForce;
    double v1kwf,v2,v3,v4min,v4max,Tpi,tal,vemin,vemax;
    double jmin,jmax,pinchStrainFactorAtLowerTemperature,pinchStrainFactorAtUpperTemperature,staticStrainFactorAtLowerTemperature,staticStrainFactorAtUpperTemperature;
    double sinus;
    double etamin,etamax;
    double zetamin, zetamax;
    double minimumAirClearance;
    double lh,lf;
    double maks;
    int zwod=0;
    FILE *plik;

    for(;;)
    {
        etamin=0.0;
        etamax=0.0;
        Tpi=0.0;
        deltamin=-1.0;
        deltamax=-1.0;
        zwod=0;
        parmin=0.0;
        parmax=0.0;
        wybor=getch();
        switch(wybor)
        {
            case 'a':
            {
                printf("Input range of calculation of center line distance between main conductor mid-points[m] in format [amin amax da]: ");
                scanf("%lf %lf %lf",&amin,&amax,&da);
                puts("Arrangement data:");
                printf("If it is two-line single-phase system input 2(else 1): ");
                scanf("%d",&lineSystem);
                printf("Input line length[m]: ");
                scanf("%lf",&l);
                printf("Choose conductor type:\n * slack conductors = l \n * strained conductors = n:\n ");
                printf("Conductor type: ");
                scanf("%c",&conductorType);
                scanf("%c",&conductorType);
                printf("Input number of sub-conductors ");
                scanf("%d",&numberOfSubConductors);
                if(numberOfSubConductors>1)
                {
                    printf("Input center-line distance[m] between sub-conductors: ");
                    scanf("%lf",&ac);
                    printf("Input sub-conductor outer diameter[mm]: ");
                    scanf("%lf",&d);
                }

                if(conductorType=='l')
                {
                    printf("Input extend of one head armature and clamp[m]: ");
                    scanf("%lf",&lh);
                    printf("Input form factor[m]: ");
                    scanf("%lf",&lf);

                    lc=l-2.0*(lh+lf);
                    l=lc;
                }
                if(conductorType=='n')
                {
                    printf("Input insulator chain length[m]: ");
                    scanf("%lf",&li);
                    lc=l-2.0*li;
                }

                printf("Input short-circuit current[kA]: ");
                scanf("%lf",&i);
                printf("Input short-circuit current duration[s]: ");
                scanf("%lf",&Tk);
                printf("Input simulation temperatures[degrees of C] in format[lowerTemperature upperTemperature]: ");
                scanf("%lf %lf",&lowerTemperature,&upperTemperature);
                puts("Mechanical data:");
                printf("Input values of static tensile forces[kN] affecting on conductor in temperature of %0.lf C and temperature of %0.lf C: ",lowerTemperature,upperTemperature);
                scanf("%lf %lf",&staticTensileForceOfConductorAtLowerTemperature,&staticTensileForceOfConductorAtUpperTemperature);
                printf("Input Young module[N/mm2]: ");
                scanf("%lf",&E);
                printf("Input resultant spring constant of both span supports of one span[N/mm]: ");
                scanf("%lf",&S);
                printf("Input mass per unit length of one sub-conductor[kg/m]: ");
                scanf("%lf",&massPerUnitLengthOfOneSubConductor);
                printf("Input number of spacers: ");
                scanf("%d",&numberOfSpacers);
                m=massPerUnitLengthOfOneSubConductor;

                if(numberOfSpacers!=0)
                {
                   ls=0.0;
                   printf("Input mass of one connection[kg]: ");
                   scanf("%lf",&massOfOneConnection);
                   printf("Input mass of one spacer[kg]: ");
                   scanf("%lf",&massOfOneSpacer);
                   m=massPerUnitLengthOfOneSubConductor+((numberOfSpacers-1)*massOfOneConnection+massOfOneSpacer)/(numberOfSubConductors*lc);
                   printf("Input distances[m] between connections[%d] from first insulator to first connection, between next connections and from last connection to last insulator:  ",numberOfSpacers+1);
                   for(it=0;it<=numberOfSpacers;it++)
                   {
                       scanf("%lf",&lk);
                       ls+=lk;
                   }
                   ls=ls/(numberOfSpacers+1);
                   puts("");
                }

                printf("Input cross-section of one sub-coductor[mm2]: ");
                scanf("%lf",&SubConductorCrossSection);
                printf("Choose conductor material: \n * copper['c']\n * Cross-Section ratio Al/Steel>6['a']\n * Cross-secrion ratio Al\Steel<=6['s']\n");
                printf("Chosen material: ");
                scanf("%c",&conductorMaterialType);
                scanf("%c",&conductorMaterialType);

                if(conductorMaterialType=='c')
                {
                    materialConstant=0.088;
                }
                if(conductorMaterialType=='a')
                {
                    materialConstant=0.27;
                }
                if(conductorMaterialType=='s')
                {
                    materialConstant=0.17;
                }

                printf("Dropper in the middle of the span?[1=Yes , 0=No]: ");
                scanf("%d",&dropper);

                if(dropper==1)
                {
                    printf("Choose arrangement: [parallel='r', perpendicular='p']: ");
                    scanf("%c",&uklad);
                    scanf("%c",&uklad);
                    printf("Input height of the dropper at maximum operating temperature of %.0lf C: ",upperTemperature);
                    scanf("%lf",&heightOfTheDropperAtUpperTemperture);
                    printf("Input width of the dropper at maximum operating temperature of %.0lf C: ",upperTemperature );
                    scanf("%lf",&widthOfTheDropperAtUpperTemperture);
                    printf("Input cord length of dropper at temperature of %.0lf C: ",upperTemperature);
                    scanf("%lf",&dropperCordLength);
                }
                printf("Is the current have to flow along half of the main conductor and along the dropper?[1=Yes, 0=No]:");
                scanf("%d",&zwod);
                if(staticTensileForceOfConductorAtLowerTemperature/(numberOfSubConductors*SubConductorCrossSection)*1000.0>50.0)
                {
                    Emin=E;
                }
                else
                    Emin=E*(0.3+0.7*sin(31.4159265359*staticTensileForceOfConductorAtLowerTemperature/(numberOfSubConductors*SubConductorCrossSection)));

                if(staticTensileForceOfConductorAtUpperTemperature/(numberOfSubConductors*SubConductorCrossSection)*1000.0>50.0)
                {
                    Emax=E;
                }
                else
                    Emax=E*(0.3+0.7*sin(31.4159265359*staticTensileForceOfConductorAtUpperTemperature/(numberOfSubConductors*SubConductorCrossSection)));

                Nmin=1.0/(S*l*1000.0)+1.0/(numberOfSubConductors*Emin*SubConductorCrossSection);
                Nmax=1.0/(S*l*1000.0)+1.0/(numberOfSubConductors*Emax*SubConductorCrossSection);




                stresmin=0.0000000040098375*numberOfSubConductors*numberOfSubConductors*m*m*l*l/(staticTensileForceOfConductorAtLowerTemperature*staticTensileForceOfConductorAtLowerTemperature*staticTensileForceOfConductorAtLowerTemperature*Nmin);
                stresmax=0.0000000040098375*numberOfSubConductors*numberOfSubConductors*m*m*l*l/(staticTensileForceOfConductorAtUpperTemperature*staticTensileForceOfConductorAtUpperTemperature*staticTensileForceOfConductorAtUpperTemperature*Nmax);
                fmin=0.00122625*numberOfSubConductors*m*l*l/staticTensileForceOfConductorAtLowerTemperature;
                fmax=0.00122625*numberOfSubConductors*m*l*l/staticTensileForceOfConductorAtUpperTemperature;



                Timin=1.79428058619*sqrt(fmin);
                Timax=1.79428058619*sqrt(fmax);
                pinchForceAtLowerTemperature=0.0;
                pinchForceAtUpperTemperature=0.0;

                if(numberOfSubConductors!=1)
                {
                    if(((1000.0*ac/d<=2.0)&&(ls>=50.0*ac))||((1000.0*ac/d<=2.5)&&(ls>=70.0*ac)))
                    {
                        k=0.0;
                    }

                    else
                    {
                        printf("Input factor for the calculation of the peak short-circuit current(depends on the ratio of R/X): ");
                        scanf("%lf",&k);
                        sinus=sin(3.14159265359/numberOfSubConductors);
                        v1kwf=(ac-d/1000.0)*massPerUnitLengthOfOneSubConductor/(0.2*sinus*sinus)*numberOfSubConductors*numberOfSubConductors*ac/(i*i*(numberOfSubConductors-1));
                        tal=-3.0/twoPiXF/log((k-1.02)/0.98);

                        while(Tpi*Tpi+Tpi/(1.0+twoPiXF*twoPiXF*tal*tal)*((twoPiXF*twoPiXF*tal*tal-1.0)*(tal+sin(2.0*twoPiXF*Tpi)/(2.0*twoPiXF))+tal*cos(2.0*twoPiXF*Tpi)-tal*tal*twoPiXF*exp(-Tpi/tal)*(tal*twoPiXF*exp(-Tpi/tal)+4.0*sin(twoPiXF*Tpi)))-v1kwf<0)
                        {
                            Tpi+=0.00001;
                        }
                        v2=v1kwf/(Tpi*Tpi);
                        v3=d*sqrt(ac*1000.0/d-1.0)/(1000.0*ac*sinus*atan(sqrt(ac*1000.0/d-1.0)));
                        shortCircuitCurrentBundleForce=(numberOfSubConductors-1)*0.2*i*i*ls*v2/(numberOfSubConductors*numberOfSubConductors*ac*v3);

                        staticStrainFactorAtLowerTemperature=1500*staticTensileForceOfConductorAtLowerTemperature*ls*ls*Nmin*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0));
                        pinchStrainFactorAtLowerTemperature=0.375*numberOfSubConductors*shortCircuitCurrentBundleForce*ls*ls*ls*Nmin*sinus*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0)*(ac-d/1000.0));
                        staticStrainFactorAtUpperTemperature=1500*staticTensileForceOfConductorAtUpperTemperature*ls*ls*Nmax*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0));
                        pinchStrainFactorAtUpperTemperature=0.375*numberOfSubConductors*shortCircuitCurrentBundleForce*ls*ls*ls*Nmax*sinus*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0)*(ac-d/1000.0));

                        jmin=sqrt(pinchStrainFactorAtLowerTemperature/(1.0+staticStrainFactorAtLowerTemperature));
                        jmax=sqrt(pinchStrainFactorAtUpperTemperature/(1.0+staticStrainFactorAtUpperTemperature));
                        zetamin=pow(jmin,0.6666);
                        zetamax=pow(jmax,0.6666);


                        if(jmin<1.0)
                        {
                            while(etamin*etamin*etamin+etamin*staticStrainFactorAtLowerTemperature-jmin*jmin*(1.0+staticStrainFactorAtLowerTemperature)*v3*sinus*atan(etamin*(1.0-d/(1000.0*ac))/(1.0-etamin*(1.0-d/(1000.0*ac))))<0.0)
                            {
                                etamin+=0.0001;
                            }

                            v4min=etamin*(ac-d/1000.0)/(ac-etamin*(ac-d/1000.0));
                            vemin=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*i*i/numberOfSubConductors*Nmin*v2*pow(ls*sinus/(etamin*(ac-d/1000.0)),4.0)*etamin*(1.0-atan(sqrt(v4min))/sqrt(v4min))-0.25);
                            pinchForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+vemin/staticStrainFactorAtLowerTemperature*etamin*etamin);

                        }

                        else
                        {
                            while(zetamin*zetamin*zetamin+zetamin*zetamin*staticStrainFactorAtLowerTemperature-jmin*jmin*(1.0+staticStrainFactorAtLowerTemperature)<0.0)
                            {
                                zetamin+=0.0001;
                            }

                            v4min=(ac*1000.0-d)/d;
                            vemin=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*i*i/numberOfSubConductors*Nmin*v2*pow(ls*sinus/(zetamin*(ac-d/1000.0)),4.0)*zetamin*(1.0-atan(sqrt(v4min))/sqrt(v4min))-0.25);
                            pinchForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+vemin*zetamin/staticStrainFactorAtLowerTemperature);
                        }


                        if(jmax<1.0)
                        {
                            while(etamax*etamax*etamax+etamax*staticStrainFactorAtUpperTemperature-jmax*jmax*(1.0+staticStrainFactorAtUpperTemperature)*v3*sinus*atan(etamax*(1.0-d/(1000.0*ac))/(1.0-etamax*(1.0-d/(1000.0*ac))))<0.0)
                            {
                                etamax+=0.0001;
                            }

                            v4max=etamax*(ac-d/1000.0)/(ac-etamax*(ac-d/1000.0));
                            vemax=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*i*i/numberOfSubConductors*Nmax*v2*pow(ls*sinus/(etamax*(ac-d/1000.0)),4.0)*etamax*(1.0-atan(sqrt(v4max))/sqrt(v4max))-0.25);
                            pinchForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+vemax/staticStrainFactorAtUpperTemperature*etamax*etamax);
                        }

                        else
                        {
                            while(zetamax*zetamax*zetamax+zetamax*zetamax*staticStrainFactorAtUpperTemperature-jmax*jmax*(1.0+staticStrainFactorAtUpperTemperature)<0.0)
                            {
                                zetamax+=0.0001;
                            }

                            v4max=(ac*1000.0-d)/d;
                            vemax=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*i*i/numberOfSubConductors*Nmax*v2*pow(ls*sinus/(zetamax*(ac-d/1000.0)),4.0)*zetamax*(1.0-atan(sqrt(v4max))/sqrt(v4max))-0.25);
                            pinchForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+vemax*zetamax/staticStrainFactorAtUpperTemperature);
                        }
                    }
                 }

                plik=fopen("dane/danea.txt","w");
                for(num=amin;num<=amax;num=num+da)
                {
                    r=oneAndAHalfDivByTenDivByStandardGravity*i*i*lc/(numberOfSubConductors*m*num*l);
                    if(zwod==1)
                        r=oneAndAHalfDivByTenDivByStandardGravity*i*i*(lc+dropperCordLength)/(2.0*numberOfSubConductors*m*num*l);
                    if(lineSystem==2)
                    {
                        r=r/0.75;
                    }
                    fi1=atan(r);

                    Trmin=Timin/(pow(1.0+r*r,0.25)*(1.0-fi1*fi1/16.0));
                    Trmax=Timax/(pow(1.0+r*r,0.25)*(1.0-fi1*fi1/16.0));

                    // Calculating end swing-out angle
                    if(Tk/Trmin>0.5)
                    {
                        fiendmin=2*fi1;
                    }
                    else
                        fiendmin=fi1*(1.0-cos(Tk/Trmin*6.28318530718));

                     if(Tk/Trmax>0.5)
                    {
                        fiendmax=2*fi1;
                    }
                    else
                        fiendmax=fi1*(1.0-cos(Tk/Trmax*6.28318530718));
                     //////////////////////////////////////////////

                    //Calculating X parameter
                    if(fiendmin>1.57079632679)
                    {
                        Xmin=1.0-r;
                    }
                    else
                        Xmin=1.0-r*sin(fiendmin);

                    if(fiendmax>1.57079632679)
                    {
                        Xmax=1.0-r;
                    }
                    else
                        Xmax=1.0-r*sin(fiendmax);
                    //////////////////////////////

                    //Calculating maximum swing-out angle
                    if(Xmin<-0.985)
                    {
                        maximumSwingOutAngleAtLowerTemperature=3.14159265359;
                    }
                    else if(Xmin<=0.766)
                    {
                        maximumSwingOutAngleAtLowerTemperature=0.17453292519+acos(Xmin);
                    }
                    else
                        maximumSwingOutAngleAtLowerTemperature=1.25*acos(Xmin);


                    if(Xmax<-0.985)
                    {
                        maximumSwingOutAngleAtUpperTemperature=3.14159265359;
                    }
                    else if(Xmax<=0.766)
                    {
                        maximumSwingOutAngleAtUpperTemperature=0.17453292519+acos(Xmax);
                    }
                    else
                        maximumSwingOutAngleAtUpperTemperature=1.25*acos(Xmax);
                     ///////////////////////////

                    //Calculating load parameter
                    if(Tk<Trmax/4.0)
                    {
                        loadParameterAtUpperTemperature=3.0*(r*sin(fiendmax)+cos(fiendmax)-1.0);
                    }
                    else
                        loadParameterAtUpperTemperature=3.0*(sqrt(1.0+r*r)-1.0);

                    if(Tk<Trmin/4.0)
                    {
                        loadParameterAtLowerTemperature=3.0*(r*sin(fiendmin)+cos(fiendmin)-1.0);
                    }
                    else
                        loadParameterAtLowerTemperature=3.0*(sqrt(1.0+r*r)-1.0);

                 //////////////////////////////////////

                 //Calculating 0-place parameter

                    row=-stresmin*(2.0+loadParameterAtLowerTemperature);
                    if(row>0.0)
                    {
                        while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row>0.0)
                        {
                            parmin=parmin+0.001;
                        }
                    }
                    else
                    {
                        while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row<0.0)
                        {
                            parmin=parmin+0.001;
                        }
                    }

                    row=-stresmax*(2.0+loadParameterAtUpperTemperature);

                    if(row>0.0)
                    {
                        while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax*loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row>0.0)
                        {
                            parmax=parmax+0.001;
                        }
                    }
                    else
                    {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax+loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row<0.0)
                      {
                          parmax=parmax+0.001;
                      }
                    }

                  ////////////////////////////////////////////////////////////////

                  //Calculating Dynamic maximum tensile force

                    shortCircuitTensileForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+loadParameterAtLowerTemperature*parmin);
                    shortCircuitTensileForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+loadParameterAtUpperTemperature*parmax);


                  //Calculating maximum dynamic span
                    elasticExpansionAtUpperTemperature=Nmax*(shortCircuitTensileForceAtUpperTemperature-staticTensileForceOfConductorAtUpperTemperature)*1000.0;
                    elasticExpansionAtLowerTemperature=Nmin*(shortCircuitTensileForceAtLowerTemperature-staticTensileForceOfConductorAtLowerTemperature)*1000.0;


                    if(Tk<Trmax/4.0)
                    {
                        thermalExpansionAtUpperTemperature=materialConstant*i*i/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Tk;
                    }
                    else
                        thermalExpansionAtUpperTemperature=0.25*materialConstant*i*i/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Trmax;


                    if(Tk<Trmin/4.0)
                    {
                        thermalExpansionAtLowerTemperature=materialConstant*i*i/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Tk;
                    }
                    else
                        thermalExpansionAtLowerTemperature=0.25*materialConstant*i*i/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Trmin;



                    dilatationFactorAtUpperTemperature=sqrt(1.0+0.375*l*l/(fmax*fmax)*(elasticExpansionAtUpperTemperature+thermalExpansionAtUpperTemperature));
                    dilatationFactorAtLowerTemperature=sqrt(1.0+0.375*l*l/(fmin*fmin)*(elasticExpansionAtLowerTemperature+thermalExpansionAtLowerTemperature));
                    if(r>=1.8)
                    {
                        formFactor=1.15;
                    }
                    else if(r>0.8)
                    {
                        formFactor=0.97+0.1*r;
                    }
                    else
                        {formFactor=1.05;}

                    fedmax=formFactor*dilatationFactorAtUpperTemperature*fmax;
                    fedmin=formFactor*dilatationFactorAtLowerTemperature*fmin;

                    //////////////////////////////

                    //calculating tensile force caused by swing-out

                    if(uklad=='r')
                    {
                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax+fedmin)*(heightOfTheDropperAtUpperTemperture+fmin+fedmin)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture))
                        {
                            deltamin=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmin*fedmin-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmin*(heightOfTheDropperAtUpperTemperture+fmax)));
                        }

                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax+fedmax)*(heightOfTheDropperAtUpperTemperture+fmax+fedmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture))
                        {
                            deltamax=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmax*fedmax-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmax*(heightOfTheDropperAtUpperTemperture+fmax)));
                        }
                    }

                    if(uklad=='p')
                    {
                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)+fedmin)
                        {
                            deltamin=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmin*fedmin-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmin*sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)))+acos((heightOfTheDropperAtUpperTemperture+fmax)/sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture));
                        }

                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)+fedmax)
                        {
                            deltamax=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmax*fedmax-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmax*sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)))+acos((heightOfTheDropperAtUpperTemperture+fmax)/sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture));
                        }
                    }

                    if(deltamin>=0)
                    {
                        if(deltamin>=fi1)
                        {
                            if(Tk<Trmin/4.0)
                            {
                                loadParameterAtLowerTemperature=3.0*(r*sin(fiendmin)+cos(fiendmin)-1.0);
                            }
                            else
                                loadParameterAtLowerTemperature=3.0*(sqrt(1.0+r*r)-1.0);
                        }
                        else
                        {
                            if(deltamin>fiendmin)
                            {
                                loadParameterAtLowerTemperature=3.0*(r*sin(fiendmin)+cos(fiendmin)-1.0);
                            }
                            else
                            {
                                loadParameterAtLowerTemperature=3.0*(r*sin(deltamin)+cos(deltamin)-1.0);
                            }
                        }

                        parmin=0.0;
                        row=-stresmin*(2.0+loadParameterAtLowerTemperature);
                  if(row>0.0)
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row>0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row<0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }


                     shortCircuitTensileForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+loadParameterAtLowerTemperature*parmin);

                    }


                    if(deltamax>=0)
                    {
                        if(deltamax>=fi1)
                        {
                            if(Tk<Trmax/4.0)
                            {
                                loadParameterAtUpperTemperature=3.0*(r*sin(fiendmax)+cos(fiendmax)-1.0);
                            }
                            else
                                loadParameterAtUpperTemperature=3.0*(sqrt(1.0+r*r)-1.0);
                        }
                        else
                        {
                            if(deltamax>fiendmax)
                            {
                                loadParameterAtUpperTemperature=3.0*(r*sin(fiendmax)+cos(fiendmax)-1.0);
                            }
                            else
                            {
                                loadParameterAtUpperTemperature=3.0*(r*sin(deltamax)+cos(deltamax)-1.0);
                            }
                        }


                     parmax=0.0;
                  row=-stresmax*(2.0+loadParameterAtUpperTemperature);

                  if(row>0.0)
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax+loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row>0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax+loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row<0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }

                     shortCircuitTensileForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+loadParameterAtUpperTemperature*parmax);


                    }


                    //Calculating a force caused by drop force
                    dropForceAtLowerTemperature=0.0;
                    dropForceAtUpperTemperature=0.0;

                    if((r>0.6)&&(maximumSwingOutAngleAtLowerTemperature>=1.221730476))
                    {
                        if(dropper==0)
                        {
                            dropForceAtLowerTemperature=1.2*staticTensileForceOfConductorAtLowerTemperature*sqrt(1.0+2.54647908947*stresmin*maximumSwingOutAngleAtLowerTemperature);
                        }
                        if((dropper==1)&&(deltamin>=1.047197551))
                        {
                            dropForceAtLowerTemperature=1.2*staticTensileForceOfConductorAtLowerTemperature*sqrt(1.0+2.54647908947*stresmin*maximumSwingOutAngleAtLowerTemperature);
                        }
                    }

                    if((r>0.6)&&(maximumSwingOutAngleAtUpperTemperature>=1.221730476))
                    {
                        if(dropper==0)
                        {
                            dropForceAtUpperTemperature=1.2*staticTensileForceOfConductorAtUpperTemperature*sqrt(1.0+2.54647908947*stresmax*maximumSwingOutAngleAtUpperTemperature);
                        }
                        if((dropper==1)&&(deltamax>=1.047197551))
                        {
                            dropForceAtUpperTemperature=1.2*staticTensileForceOfConductorAtUpperTemperature*sqrt(1.0+2.54647908947*stresmax*maximumSwingOutAngleAtUpperTemperature);
                        }
                    }



                    if(conductorType=='l')
                    {
                        if(maximumSwingOutAngleAtLowerTemperature<1.57079632679)
                        {
                            bhmin=fedmin*sin(maximumSwingOutAngleAtLowerTemperature);
                        }
                        else
                            bhmin=fedmin;

                        if(maximumSwingOutAngleAtUpperTemperature<1.57079632679)
                        {
                            bhmax=fedmax*sin(maximumSwingOutAngleAtUpperTemperature);
                        }
                        else
                            bhmax=fedmax;
                    }
                    if(conductorType=='n')
                    {
                        if(dropper==0)
                        {
                            if(maximumSwingOutAngleAtLowerTemperature<fi1)
                            {
                                bhmin=fedmin*sin(maximumSwingOutAngleAtLowerTemperature);
                            }
                            else
                                bhmin=fedmin*sin(fi1);

                            if(maximumSwingOutAngleAtUpperTemperature<fi1)
                            {
                                bhmax=fedmax*sin(maximumSwingOutAngleAtUpperTemperature);
                            }
                            else
                                bhmax=fedmax*sin(fi1);
                        }
                        if(dropper==1)
                        {
                            if(deltamin<maximumSwingOutAngleAtLowerTemperature)
                            {
                                if(maximumSwingOutAngleAtLowerTemperature<fi1)
                                {
                                    bhmin=fedmin*sin(maximumSwingOutAngleAtLowerTemperature);
                                }
                                else
                                    bhmin=fedmin*sin(fi1);
                            }
                            else
                            {
                                if(deltamin<fi1)
                                {
                                    bhmin=fedmin*sin(deltamin);
                                }
                                else
                                    bhmin=fedmin*sin(fi1);
                            }


                            if(deltamax<maximumSwingOutAngleAtUpperTemperature)
                            {
                                if(deltamax<fi1)
                                {
                                    bhmax=fedmax*sin(deltamax);
                                }
                                else
                                    bhmax=fedmax*sin(fi1);
                            }
                            else
                            {

                                if(maximumSwingOutAngleAtUpperTemperature<fi1)
                                {
                                    bhmax=fedmax*sin(maximumSwingOutAngleAtUpperTemperature);
                                }
                                else
                                    bhmax=fedmax*sin(fi1);
                            }

                        }
                    }


                    horizontalSpanDisplacement=(bhmin>bhmax) ? bhmin : bhmax;
                    minimumAirClearance=num-2.0*horizontalSpanDisplacement;

                    if(numberOfSubConductors!=1)
                    {
                        if(k==0.0)
                        {
                            pinchForceAtLowerTemperature=1.1*shortCircuitTensileForceAtLowerTemperature;
                            pinchForceAtUpperTemperature=1.1*shortCircuitTensileForceAtUpperTemperature;
                        }
                        //printf("Maximum Dynamic tensile force caused by pinch force: %.3lf\n",(pinchForceAtLowerTemperature >= pinchForceAtUpperTemperature) ? pinchForceAtLowerTemperature : pinchForceAtUpperTemperature);

                    }

                    maximumShortCircuitTensileForce=(shortCircuitTensileForceAtLowerTemperature>=shortCircuitTensileForceAtUpperTemperature)? shortCircuitTensileForceAtLowerTemperature : shortCircuitTensileForceAtUpperTemperature;
                    maximumDropForce=(dropForceAtLowerTemperature>=dropForceAtUpperTemperature)? dropForceAtLowerTemperature : dropForceAtUpperTemperature;
                    maximumPinchForce=(pinchForceAtLowerTemperature>=pinchForceAtUpperTemperature)? pinchForceAtLowerTemperature : pinchForceAtUpperTemperature;
                    if(conductorType=='l')
                        maximumShortCircuitTensileForce*=1.5;
                    maks=maximumShortCircuitTensileForce;
                    if(maks<maximumDropForce)
                    {
                        maks=maximumDropForce;
                    }
                    if(maks<maximumPinchForce)
                    {
                        maks=maximumPinchForce;
                    }

                if(amin+da>amax)
                {
                    printf("Additional parameters: \n");

                    printf("    r= %.2lf\n",r);
                    printf("    Maximum Swing-Out Angle At Lower Temperature = %lf , Maximum Swing-Out Angle At Upper Temperature = %lf\n",maximumSwingOutAngleAtLowerTemperature,maximumSwingOutAngleAtUpperTemperature);
                    printf("    deltamin = %lf , deltamax = %lf\n",deltamin,deltamax);
                    printf("    parmin: %.3lf , parmax: %.3lf\n",parmin,parmax);
                    printf("    loadParameterAtLowerTemperature: %.2lf ,loadParameterAtUpperTemperature: %.2lf\n",loadParameterAtLowerTemperature,loadParameterAtUpperTemperature);
                    printf("    stresmin: %.2lf , stresmax: %.2lf\n",stresmin,stresmax);
                    printf("    fesmin: %.2lf , fesmax: %.2lf\n",fmin,fmax);
                    printf("    fedmin: %.2lf , fedmax: %.2lf\n",fedmin,fedmax);
                    printf("    Elastic Expansion At Lower Temperature: %.3lf , Elastic Expansion At Upper Temperature: %.3lf\n",elasticExpansionAtLowerTemperature*1000.0,elasticExpansionAtUpperTemperature*1000.0);
                    printf("    Thermal Expansion At Lower Temperature: %.3lf , Thermal Expansion At Upper Temperature: %.3lf\n",thermalExpansionAtLowerTemperature*10000.0,thermalExpansionAtUpperTemperature*10000.0);
                    printf("    Material constant: %.3lf\n",materialConstant);
                    printf("    Dilatation Factor At Lower Temperature: %.3lf , Dilatation Factor At Upper Temperature: %.3lf\n",dilatationFactorAtLowerTemperature,dilatationFactorAtUpperTemperature);
                    printf("    formFactor: %.3lf\n",formFactor);
                puts("Forces:");
                    printf("    Maximum Short-Circuit Tensile Force: %.2lf kN\n",maximumShortCircuitTensileForce);
                    printf("    Maximum Drop Force: %.2lf kN\n",maximumDropForce);
                    printf("    Maximum Pinch Force: %.2lf kN\n",maximumPinchForce);
                    printf("    Overall Maximum Force: %.2lf kN\n",maks);
                puts("Effects of the force");
                    printf("    Horizontal span displacement: %.2lf m\n",horizontalSpanDisplacement);
                    printf("    Minimum air clearance: %.2lf m\n",minimumAirClearance);
                }
                else
                    fprintf(plik,"%lf %lf %lf %lf %lf %lf %lf\n",num,maximumShortCircuitTensileForce,maximumDropForce,maximumPinchForce,maks,horizontalSpanDisplacement,minimumAirClearance);
                }
                fclose(plik);
                if(amin+da<=amax)
                    ShellExecute(NULL, NULL, "dane\\danea.m", NULL, NULL, SW_SHOWNORMAL);
                system("pause");
                system("cls");
                puts("Wybierz parmetr wzgledem, ktorego chcesz wykonac obliczenia:\n *'a'-odstepy miedzy liniami; \n *'i' - prad zwarciowy'; \n *'s'-przekroj przewodu");
            }
            break;

            case 'i':
            {
                printf("Wpisz zakres pradow zwarciowych[m] w formacie [imin imax di]: ");
                scanf("%lf %lf %lf",&amin,&amax,&da);
                puts("Arrangement data:");
                printf("If it is two-line single-phase system input 2(else 1): ");
                scanf("%d",&lineSystem);
                printf("Input line length[m]: ");
                scanf("%lf",&l);
                printf("Wprowadz srednia odleglosc miedzy liniami: ");
                scanf("%lf",&i);
                printf("Choose conductor type:\n * slack conductors = l \n * strained conductors = n:\n ");
                printf("Conductor type: ");
                scanf("%c",&conductorType);
                scanf("%c",&conductorType);
                printf("Input number of sub-conductors ");
                scanf("%d",&numberOfSubConductors);
                if(numberOfSubConductors>1)
                {
                    printf("Input center-line distance[m] between sub-conductors: ");
                    scanf("%lf",&ac);
                    printf("Input sub-conductor outer diameter[mm]: ");
                    scanf("%lf",&d);
                }


                if(conductorType=='l')
                {
                    printf("Input extend of one head armature and clamp[m]: ");
                    scanf("%lf",&lh);
                    printf("Input form factor[m]: ");
                    scanf("%lf",&lf);

                    lc=l-2.0*(lh+lf);
                    l=lc;
                }
                if(conductorType=='n')
                {
                    printf("Input insulator chain length[m]: ");
                    scanf("%lf",&li);
                    lc=l-2.0*li;
                }

                printf("Input short-circuit current duration[s]: ");
                scanf("%lf",&Tk);
                printf("Input simulation temperatures[degrees of C] in format[lowerTemperature upperTemperature]: ");
                scanf("%lf %lf",&lowerTemperature,&upperTemperature);
                puts("Mechanical data:");
                printf("Input values of static tensile forces[kN] affecting on conductor in temperature of %0.lf C and temperature of %0.lf C: ",lowerTemperature,upperTemperature);
                scanf("%lf %lf",&staticTensileForceOfConductorAtLowerTemperature,&staticTensileForceOfConductorAtUpperTemperature);
                printf("Input Young module[N/mm2]: ");
                scanf("%lf",&E);
                printf("Input resultant spring constant of both span supports of one span[N/mm]: ");
                scanf("%lf",&S);
                printf("Input mass per unit length of one sub-conductor[kg/m]: ");
                scanf("%lf",&massPerUnitLengthOfOneSubConductor);
                printf("Input number of spacers: ");
                scanf("%d",&numberOfSpacers);
                m=massPerUnitLengthOfOneSubConductor;

                if(numberOfSpacers!=0)
                {
                   ls=0.0;
                   printf("Input mass of one connection[kg]: ");
                   scanf("%lf",&massOfOneConnection);
                   printf("Input mass of one spacer[kg]: ");
                   scanf("%lf",&massOfOneSpacer);
                   m=massPerUnitLengthOfOneSubConductor+((numberOfSpacers-1)*massOfOneConnection+massOfOneSpacer)/(numberOfSubConductors*lc);
                   printf("Input distances[m] between connections[%d] from first insulator to first connection, between next connections and from last connection to last insulator:  ",numberOfSpacers+1);
                   for(it=0;it<=numberOfSpacers;it++)
                   {
                       scanf("%lf",&lk);
                       ls+=lk;
                   }
                   ls=ls/(numberOfSpacers+1);
                }

                printf("Conductor Cross-section[mm2]: ");
                scanf("%lf",&SubConductorCrossSection);
                printf("Choose conductor material: \n * copper['c']\n * Cross-Section ratio Al/Steel>6['a']\n * Cross-secrion ratio Al\Steel<=6['s']\n");
                printf("Chosen material: ");
                scanf("%c",&conductorMaterialType);
                scanf("%c",&conductorMaterialType);


                if(conductorMaterialType=='c')
                {
                    materialConstant=0.088;
                }
                if(conductorMaterialType=='a')
                {
                    materialConstant=0.27;
                }
                if(conductorMaterialType=='s')
                {
                    materialConstant=0.17;
                }

                printf("Dropper in the middle of the span?[1=Yes , 0=No]: ");
                scanf("%d",&dropper);

                if(dropper==1)
                {
                    printf("Choose arrangement: [parallel='r', perpendicular='p']: ");
                    scanf("%c",&uklad);
                    scanf("%c",&uklad);
                    printf("Input height of the dropper at maximum operating temperature of %.0lf C: ",upperTemperature);
                    scanf("%lf",&heightOfTheDropperAtUpperTemperture);
                    printf("Input width of the dropper at maximum operating temperature of %.0lf C: ",upperTemperature );
                    scanf("%lf",&widthOfTheDropperAtUpperTemperture);
                    printf("Input cord length of dropper at temperature of %.0lf C: ",upperTemperature);
                    scanf("%lf",&dropperCordLength);
                }
                printf("Is the current have to flow along half of the main conductor and along the dropper?[1=Yes, 0=No]:");
                scanf("%d",&zwod);
                if(staticTensileForceOfConductorAtLowerTemperature/(numberOfSubConductors*SubConductorCrossSection)*1000.0>50.0)
                {
                    Emin=E;
                }
                else
                    Emin=E*(0.3+0.7*sin(31.4159265359*staticTensileForceOfConductorAtLowerTemperature/(numberOfSubConductors*SubConductorCrossSection)));

                if(staticTensileForceOfConductorAtUpperTemperature/(numberOfSubConductors*SubConductorCrossSection)*1000.0>50.0)
                {
                    Emax=E;
                }
                else
                    Emax=E*(0.3+0.7*sin(31.4159265359*staticTensileForceOfConductorAtUpperTemperature/(numberOfSubConductors*SubConductorCrossSection)));

                Nmin=1.0/(S*l*1000.0)+1.0/(numberOfSubConductors*Emin*SubConductorCrossSection);
                Nmax=1.0/(S*l*1000.0)+1.0/(numberOfSubConductors*Emax*SubConductorCrossSection);


                stresmin=0.0000000040098375*numberOfSubConductors*numberOfSubConductors*m*m*l*l/(staticTensileForceOfConductorAtLowerTemperature*staticTensileForceOfConductorAtLowerTemperature*staticTensileForceOfConductorAtLowerTemperature*Nmin);
                stresmax=0.0000000040098375*numberOfSubConductors*numberOfSubConductors*m*m*l*l/(staticTensileForceOfConductorAtUpperTemperature*staticTensileForceOfConductorAtUpperTemperature*staticTensileForceOfConductorAtUpperTemperature*Nmax);
                fmin=0.00122625*numberOfSubConductors*m*l*l/staticTensileForceOfConductorAtLowerTemperature;
                fmax=0.00122625*numberOfSubConductors*m*l*l/staticTensileForceOfConductorAtUpperTemperature;



                Timin=1.79428058619*sqrt(fmin);
                Timax=1.79428058619*sqrt(fmax);
                pinchForceAtLowerTemperature=0.0;
                pinchForceAtUpperTemperature=0.0;

                if(((1000.0*ac/d<=2.0)&&(ls>=50.0*ac))||((1000.0*ac/d<=2.5)&&(ls>=70.0*ac)))
                {
                    k=0.0;
                }

                else
                {
                    printf("Input factor for the calculation of the peak short-circuit current(depends on the ratio of R/X): ");
                    scanf("%lf",&k);
                }

                plik=fopen("dane/danei.txt","w");
                for(num=amin;num<=amax;num=num+da)
                {

                    if(k!=0.0)
                    {
                        sinus=sin(3.14159265359/numberOfSubConductors);
                        v1kwf=(ac-d/1000.0)*massPerUnitLengthOfOneSubConductor/(0.2*sinus*sinus)*numberOfSubConductors*numberOfSubConductors*ac/(num*num*(numberOfSubConductors-1));
                        tal=-3.0/twoPiXF/log((k-1.02)/0.98);

                        while(Tpi*Tpi+Tpi/(1.0+twoPiXF*twoPiXF*tal*tal)*((twoPiXF*twoPiXF*tal*tal-1.0)*(tal+sin(2.0*twoPiXF*Tpi)/(2.0*twoPiXF))+tal*cos(2.0*twoPiXF*Tpi)-tal*tal*twoPiXF*exp(-Tpi/tal)*(tal*twoPiXF*exp(-Tpi/tal)+4.0*sin(twoPiXF*Tpi)))-v1kwf<0)
                        {
                            Tpi+=0.00001;
                        }
                        v2=v1kwf/(Tpi*Tpi);
                        v3=d*sqrt(ac*1000.0/d-1.0)/(1000.0*ac*sinus*atan(sqrt(ac*1000.0/d-1.0)));
                        shortCircuitCurrentBundleForce=(numberOfSubConductors-1)*0.2*num*num*ls*v2/(numberOfSubConductors*numberOfSubConductors*ac*v3);



                        staticStrainFactorAtLowerTemperature=1500*staticTensileForceOfConductorAtLowerTemperature*ls*ls*Nmin*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0));
                        pinchStrainFactorAtLowerTemperature=0.375*numberOfSubConductors*shortCircuitCurrentBundleForce*ls*ls*ls*Nmin*sinus*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0)*(ac-d/1000.0));
                        staticStrainFactorAtUpperTemperature=1500*staticTensileForceOfConductorAtUpperTemperature*ls*ls*Nmax*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0));
                        pinchStrainFactorAtUpperTemperature=0.375*numberOfSubConductors*shortCircuitCurrentBundleForce*ls*ls*ls*Nmax*sinus*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0)*(ac-d/1000.0));

                        jmin=sqrt(pinchStrainFactorAtLowerTemperature/(1.0+staticStrainFactorAtLowerTemperature));
                        jmax=sqrt(pinchStrainFactorAtUpperTemperature/(1.0+staticStrainFactorAtUpperTemperature));
                        zetamin=pow(jmin,0.6666);
                        zetamax=pow(jmax,0.6666);


                        if(jmin<1.0)
                        {
                            while(etamin*etamin*etamin+etamin*staticStrainFactorAtLowerTemperature-jmin*jmin*(1.0+staticStrainFactorAtLowerTemperature)*v3*sinus*atan(etamin*(1.0-d/(1000.0*ac))/(1.0-etamin*(1.0-d/(1000.0*ac))))<0.0)
                            {
                                etamin+=0.0001;
                            }

                            v4min=etamin*(ac-d/1000.0)/(ac-etamin*(ac-d/1000.0));
                            vemin=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*num*num/numberOfSubConductors*Nmin*v2*pow(ls*sinus/(etamin*(ac-d/1000.0)),4.0)*etamin*(1.0-atan(sqrt(v4min))/sqrt(v4min))-0.25);
                            pinchForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+vemin/staticStrainFactorAtLowerTemperature*etamin*etamin);

                        }

                        else
                        {
                            while(zetamin*zetamin*zetamin+zetamin*zetamin*staticStrainFactorAtLowerTemperature-jmin*jmin*(1.0+staticStrainFactorAtLowerTemperature)<0.0)
                            {
                                zetamin+=0.0001;
                            }

                        v4min=(ac*1000.0-d)/d;
                        vemin=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*num*num/numberOfSubConductors*Nmin*v2*pow(ls*sinus/(zetamin*(ac-d/1000.0)),4.0)*zetamin*(1.0-atan(sqrt(v4min))/sqrt(v4min))-0.25);
                        pinchForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+vemin*zetamin/staticStrainFactorAtLowerTemperature);
                        }


                     if(jmax<1.0)
                    {
                        while(etamax*etamax*etamax+etamax*staticStrainFactorAtUpperTemperature-jmax*jmax*(1.0+staticStrainFactorAtUpperTemperature)*v3*sinus*atan(etamax*(1.0-d/(1000.0*ac))/(1.0-etamax*(1.0-d/(1000.0*ac))))<0.0)
                        {
                            etamax+=0.0001;
                        }

                        v4max=etamax*(ac-d/1000.0)/(ac-etamax*(ac-d/1000.0));
                        vemax=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*num*num/numberOfSubConductors*Nmax*v2*pow(ls*sinus/(etamax*(ac-d/1000.0)),4.0)*etamax*(1.0-atan(sqrt(v4max))/sqrt(v4max))-0.25);
                        pinchForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+vemax/staticStrainFactorAtUpperTemperature*etamax*etamax);
                    }

                    else
                        {
                            while(zetamax*zetamax*zetamax+zetamax*zetamax*staticStrainFactorAtUpperTemperature-jmax*jmax*(1.0+staticStrainFactorAtUpperTemperature)<0.0)
                            {
                                zetamax+=0.0001;
                            }

                            v4max=(ac*1000.0-d)/d;
                            vemax=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*num*num/numberOfSubConductors*Nmax*v2*pow(ls*sinus/(zetamax*(ac-d/1000.0)),4.0)*zetamax*(1.0-atan(sqrt(v4max))/sqrt(v4max))-0.25);
                            pinchForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+vemax*zetamax/staticStrainFactorAtUpperTemperature);
                        }

                    }

                        r=oneAndAHalfDivByTenDivByStandardGravity*num*num*lc/(numberOfSubConductors*m*i*l);

                        if(zwod==1)
                            r=oneAndAHalfDivByTenDivByStandardGravity*num*num*(lc+dropperCordLength)/(2.0*numberOfSubConductors*m*i*l);

                        if(lineSystem==2)
                        {
                            r=r/0.75;
                        }
                    fi1=atan(r);

                    Trmin=Timin/(pow(1.0+r*r,0.25)*(1.0-fi1*fi1/16.0));
                    Trmax=Timax/(pow(1.0+r*r,0.25)*(1.0-fi1*fi1/16.0));

                    // Calculating end swing-out angle
                    if(Tk/Trmin>0.5)
                    {
                        fiendmin=2*fi1;
                    }
                    else
                        fiendmin=fi1*(1.0-cos(Tk/Trmin*6.28318530718));

                     if(Tk/Trmax>0.5)
                    {
                        fiendmax=2*fi1;
                    }
                    else
                        fiendmax=fi1*(1.0-cos(Tk/Trmax*6.28318530718));
                     //////////////////////////////////////////////

                    //Calculating X parameter
                    if(fiendmin>1.57079632679)
                    {
                        Xmin=1.0-r;
                    }
                    else
                        Xmin=1.0-r*sin(fiendmin);

                    if(fiendmax>1.57079632679)
                    {
                        Xmax=1.0-r;
                    }
                    else
                        Xmax=1.0-r*sin(fiendmax);
                    //////////////////////////////

                    //Calculating maximum swing-out angle
                    if(Xmin<-0.985)
                    {
                        maximumSwingOutAngleAtLowerTemperature=3.14159265359;
                    }
                    else if(Xmin<=0.766)
                    {
                        maximumSwingOutAngleAtLowerTemperature=0.17453292519+acos(Xmin);
                    }
                    else
                        maximumSwingOutAngleAtLowerTemperature=1.25*acos(Xmin);


                    if(Xmax<-0.985)
                    {
                        maximumSwingOutAngleAtUpperTemperature=3.14159265359;
                    }
                    else if(Xmax<=0.766)
                    {
                        maximumSwingOutAngleAtUpperTemperature=0.17453292519+acos(Xmax);
                    }
                    else
                        maximumSwingOutAngleAtUpperTemperature=1.25*acos(Xmax);
                     ///////////////////////////

                    //Calculating load parameter
                    if(Tk<Trmax/4.0)
                    {
                        loadParameterAtUpperTemperature=3.0*(r*sin(fiendmax)+cos(fiendmax)-1.0);
                    }
                    else
                        loadParameterAtUpperTemperature=3.0*(sqrt(1.0+r*r)-1.0);

                    if(Tk<Trmin/4.0)
                    {
                        loadParameterAtLowerTemperature=3.0*(r*sin(fiendmin)+cos(fiendmin)-1.0);
                    }
                    else
                        loadParameterAtLowerTemperature=3.0*(sqrt(1.0+r*r)-1.0);

                 //////////////////////////////////////

                 //Calculating 0-place parameter

                  row=-stresmin*(2.0+loadParameterAtLowerTemperature);
                  if(row>0.0)
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row>0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row<0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }

                  row=-stresmax*(2.0+loadParameterAtUpperTemperature);

                  if(row>0.0)
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax*loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row>0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax+loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row<0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }

                  ////////////////////////////////////////////////////////////////

                  //Calculating Dynamic maximum tensile force

                    shortCircuitTensileForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+loadParameterAtLowerTemperature*parmin);
                    shortCircuitTensileForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+loadParameterAtUpperTemperature*parmax);


                  //Calculating maximum dynamic span
                    elasticExpansionAtUpperTemperature=Nmax*(shortCircuitTensileForceAtUpperTemperature-staticTensileForceOfConductorAtUpperTemperature)*1000.0;
                    elasticExpansionAtLowerTemperature=Nmin*(shortCircuitTensileForceAtLowerTemperature-staticTensileForceOfConductorAtLowerTemperature)*1000.0;


                    if(Tk<Trmax/4.0)
                    {
                        thermalExpansionAtUpperTemperature=materialConstant*num*num/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Tk;
                    }
                    else
                        thermalExpansionAtUpperTemperature=0.25*materialConstant*num*num/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Trmax;


                    if(Tk<Trmin/4.0)
                    {
                        thermalExpansionAtLowerTemperature=materialConstant*num*num/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Tk;
                    }
                    else
                        thermalExpansionAtLowerTemperature=0.25*materialConstant*num*num/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Trmin;



                    dilatationFactorAtUpperTemperature=sqrt(1.0+0.375*l*l/(fmax*fmax)*(elasticExpansionAtUpperTemperature+thermalExpansionAtUpperTemperature));
                    dilatationFactorAtLowerTemperature=sqrt(1.0+0.375*l*l/(fmin*fmin)*(elasticExpansionAtLowerTemperature+thermalExpansionAtLowerTemperature));
                    if(r>=1.8)
                    {
                        formFactor=1.15;
                    }
                    else if(r>0.8)
                    {
                        formFactor=0.97+0.1*r;
                    }
                    else
                        {formFactor=1.05;}

                    fedmax=formFactor*dilatationFactorAtUpperTemperature*fmax;
                    fedmin=formFactor*dilatationFactorAtLowerTemperature*fmin;

                    //////////////////////////////

                    //calculating tensile force caused by swing-out

                    if(uklad=='r')
                    {
                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax+fedmin)*(heightOfTheDropperAtUpperTemperture+fmin+fedmin)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture))
                        {
                            deltamin=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmin*fedmin-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmin*(heightOfTheDropperAtUpperTemperture+fmax)));
                        }

                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax+fedmax)*(heightOfTheDropperAtUpperTemperture+fmax+fedmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture))
                        {
                            deltamax=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmax*fedmax-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmax*(heightOfTheDropperAtUpperTemperture+fmax)));
                        }
                    }

                    if(uklad=='p')
                    {
                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)+fedmin)
                        {
                            deltamin=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmin*fedmin-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmin*sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)))+acos((heightOfTheDropperAtUpperTemperture+fmax)/sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture));
                        }

                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)+fedmax)
                        {
                            deltamax=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmax*fedmax-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmax*sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)))+acos((heightOfTheDropperAtUpperTemperture+fmax)/sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture));
                        }
                    }

                    if(deltamin>=0)
                    {
                        if(deltamin>=fi1)
                        {
                            if(Tk<Trmin/4.0)
                            {
                                loadParameterAtLowerTemperature=3.0*(r*sin(fiendmin)+cos(fiendmin)-1.0);
                            }
                            else
                                loadParameterAtLowerTemperature=3.0*(sqrt(1.0+r*r)-1.0);
                        }
                        else
                        {
                            if(deltamin>fiendmin)
                            {
                                loadParameterAtLowerTemperature=3.0*(r*sin(fiendmin)+cos(fiendmin)-1.0);
                            }
                            else
                            {
                                loadParameterAtLowerTemperature=3.0*(r*sin(deltamin)+cos(deltamin)-1.0);
                            }
                        }

                        parmin=0.0;
                        row=-stresmin*(2.0+loadParameterAtLowerTemperature);
                  if(row>0.0)
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row>0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row<0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }


                     shortCircuitTensileForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+loadParameterAtLowerTemperature*parmin);

                    }


                    if(deltamax>=0)
                    {
                        if(deltamax>=fi1)
                        {
                            if(Tk<Trmax/4.0)
                            {
                                loadParameterAtUpperTemperature=3.0*(r*sin(fiendmax)+cos(fiendmax)-1.0);
                            }
                            else
                                loadParameterAtUpperTemperature=3.0*(sqrt(1.0+r*r)-1.0);
                        }
                        else
                        {
                            if(deltamax>fiendmax)
                            {
                                loadParameterAtUpperTemperature=3.0*(r*sin(fiendmax)+cos(fiendmax)-1.0);
                            }
                            else
                            {
                                loadParameterAtUpperTemperature=3.0*(r*sin(deltamax)+cos(deltamax)-1.0);
                            }
                        }


                     parmax=0.0;
                  row=-stresmax*(2.0+loadParameterAtUpperTemperature);

                  if(row>0.0)
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax+loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row>0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax+loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row<0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }

                     shortCircuitTensileForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+loadParameterAtUpperTemperature*parmax);

                }

                    //Calculating a force caused by drop force
                    dropForceAtLowerTemperature=0.0;
                    dropForceAtUpperTemperature=0.0;

                    if((r>0.6)&&(maximumSwingOutAngleAtLowerTemperature>=1.221730476))
                    {
                        if(dropper==0)
                        {
                            dropForceAtLowerTemperature=1.2*staticTensileForceOfConductorAtLowerTemperature*sqrt(1.0+2.54647908947*stresmin*maximumSwingOutAngleAtLowerTemperature);
                        }
                        if((dropper==1)&&(deltamin>=1.047197551))
                        {
                            dropForceAtLowerTemperature=1.2*staticTensileForceOfConductorAtLowerTemperature*sqrt(1.0+2.54647908947*stresmin*maximumSwingOutAngleAtLowerTemperature);
                        }
                    }

                    if((r>0.6)&&(maximumSwingOutAngleAtUpperTemperature>=1.221730476))
                    {
                        if(dropper==0)
                        {
                            dropForceAtUpperTemperature=1.2*staticTensileForceOfConductorAtUpperTemperature*sqrt(1.0+2.54647908947*stresmax*maximumSwingOutAngleAtUpperTemperature);
                        }
                        if((dropper==1)&&(deltamax>=1.047197551))
                        {
                            dropForceAtUpperTemperature=1.2*staticTensileForceOfConductorAtUpperTemperature*sqrt(1.0+2.54647908947*stresmax*maximumSwingOutAngleAtUpperTemperature);
                        }
                    }



                    if(conductorType=='l')
                    {
                        if(maximumSwingOutAngleAtLowerTemperature<1.57079632679)
                        {
                            bhmin=fedmin*sin(maximumSwingOutAngleAtLowerTemperature);
                        }
                        else
                            bhmin=fedmin;

                        if(maximumSwingOutAngleAtUpperTemperature<1.57079632679)
                        {
                            bhmax=fedmax*sin(maximumSwingOutAngleAtUpperTemperature);
                        }
                        else
                            bhmax=fedmax;
                    }
                    if(conductorType=='n')
                    {
                        if(dropper==0)
                        {
                            if(maximumSwingOutAngleAtLowerTemperature<fi1)
                            {
                                bhmin=fedmin*sin(maximumSwingOutAngleAtLowerTemperature);
                            }
                            else
                                bhmin=fedmin*sin(fi1);

                            if(maximumSwingOutAngleAtUpperTemperature<fi1)
                            {
                                bhmax=fedmax*sin(maximumSwingOutAngleAtUpperTemperature);
                            }
                            else
                                bhmax=fedmax*sin(fi1);
                        }
                        if(dropper==1)
                        {
                            if(deltamin<maximumSwingOutAngleAtLowerTemperature)
                            {
                                if(maximumSwingOutAngleAtLowerTemperature<fi1)
                                {
                                    bhmin=fedmin*sin(maximumSwingOutAngleAtLowerTemperature);
                                }
                                else
                                    bhmin=fedmin*sin(fi1);
                            }
                            else
                            {
                                if(deltamin<fi1)
                                {
                                    bhmin=fedmin*sin(deltamin);
                                }
                                else
                                    bhmin=fedmin*sin(fi1);
                            }


                            if(deltamax<maximumSwingOutAngleAtUpperTemperature)
                            {
                                if(deltamax<fi1)
                                {
                                    bhmax=fedmax*sin(deltamax);
                                }
                                else
                                    bhmax=fedmax*sin(fi1);
                            }
                            else
                            {

                                if(maximumSwingOutAngleAtUpperTemperature<fi1)
                                {
                                    bhmax=fedmax*sin(maximumSwingOutAngleAtUpperTemperature);
                                }
                                else
                                    bhmax=fedmax*sin(fi1);
                            }

                        }
                    }

                    horizontalSpanDisplacement=(bhmin>bhmax) ? bhmin : bhmax;
                    minimumAirClearance=i-2.0*horizontalSpanDisplacement;


                if(numberOfSubConductors!=1)
                {
                    if(k==0.0)
                    {
                        pinchForceAtLowerTemperature=1.1*shortCircuitTensileForceAtLowerTemperature;
                        pinchForceAtUpperTemperature=1.1*shortCircuitTensileForceAtUpperTemperature;
                    }
                }
                    maximumShortCircuitTensileForce=(shortCircuitTensileForceAtLowerTemperature>=shortCircuitTensileForceAtUpperTemperature)? shortCircuitTensileForceAtLowerTemperature : shortCircuitTensileForceAtUpperTemperature;
                    maximumDropForce=(dropForceAtLowerTemperature>=dropForceAtUpperTemperature)? dropForceAtLowerTemperature : dropForceAtUpperTemperature;
                    maximumPinchForce=(pinchForceAtLowerTemperature>=pinchForceAtUpperTemperature)? pinchForceAtLowerTemperature : pinchForceAtUpperTemperature;
                    if(conductorType=='l')
                        maximumShortCircuitTensileForce*=1.5;
                    maks=maximumShortCircuitTensileForce;
                    if(maks<maximumDropForce)
                    {
                        maks=maximumDropForce;
                    }
                    if(maks<maximumPinchForce)
                    {
                        maks=maximumPinchForce;
                    }

                if(amin+da>amax)
                {
                printf("Additional parameters: \n");

                    printf("    r= %.2lf\n",r);
                    printf("    Maximum Swing-Out Angle At Lower Temperature = %lf , Maximum Swing-Out Angle At Upper Temperature = %lf\n",maximumSwingOutAngleAtLowerTemperature,maximumSwingOutAngleAtUpperTemperature);
                    printf("    deltamin = %lf , deltamax = %lf\n",deltamin,deltamax);
                    printf("    parmin: %.3lf , parmax: %.3lf\n",parmin,parmax);
                    printf("    loadParameterAtLowerTemperature: %.2lf ,loadParameterAtUpperTemperature: %.2lf\n",loadParameterAtLowerTemperature,loadParameterAtUpperTemperature);
                    printf("    stresmin: %.2lf , stresmax: %.2lf\n",stresmin,stresmax);
                    printf("    fesmin: %.2lf , fesmax: %.2lf\n",fmin,fmax);
                    printf("    fedmin: %.2lf , fedmax: %.2lf\n",fedmin,fedmax);
                    printf("    Elastic Expansion At Lower Temperature: %.3lf , Elastic Expansion At Upper Temperature: %.3lf\n",elasticExpansionAtLowerTemperature*1000.0,elasticExpansionAtUpperTemperature*1000.0);
                    printf("    Thermal Expansion At Lower Temperature: %.3lf , Thermal Expansion At Upper Temperature: %.3lf\n",thermalExpansionAtLowerTemperature*10000.0,thermalExpansionAtUpperTemperature*10000.0);
                    printf("    Material constant: %.3lf\n",materialConstant);
                    printf("    Dilatation Factor At Lower Temperature: %.3lf , Dilatation Factor At Upper Temperature: %.3lf\n",dilatationFactorAtLowerTemperature,dilatationFactorAtUpperTemperature);
                    printf("    formFactor: %.3lf\n",formFactor);
                puts("Forces:");
                    printf("    Maximum Short-Circuit Tensile Force: %.2lf kN\n",maximumShortCircuitTensileForce);
                    printf("    Maximum Drop Force: %.2lf kN\n",maximumDropForce);
                    printf("    Maximum Pinch Force: %.2lf kN\n",maximumPinchForce);
                    printf("    Overall Maximum Force: %.2lf kN\n",maks);
                puts("Effects of the force");
                    printf("    Horizontal span displacement: %.2lf m\n",horizontalSpanDisplacement);
                    printf("    Minimum air clearance: %.2lf m\n",minimumAirClearance);
                }
                else
                    fprintf(plik,"%lf %lf %lf %lf %lf %lf %lf\n",num,maximumShortCircuitTensileForce,maximumDropForce,maximumPinchForce,maks,horizontalSpanDisplacement,minimumAirClearance);
                }
                fclose(plik);
                if(amin+da<=amax)
                    ShellExecute(NULL, NULL, "dane\\danei.m", NULL, NULL, SW_SHOWNORMAL);

                system("pause");
                system("cls");
                puts("Choose a parameter in terms of calculations will be performed:\n *'a'-center line distance between main conductor mid-points; \n *'i' - three phase initial symmetrical short-circuit current'; \n *'t'-duration of short-circuit current");

            }
            break;

            case 't':
            {
                printf("Input range of short-circuit duration in format[tmin tmax dt] ");
                scanf("%lf %lf %lf",&amin,&amax,&da);
                puts("Arrangement data:");
                printf("If it is two-line single-phase system input 2(else 1): ");
                scanf("%d",&lineSystem);
                printf("Input line length[m]: ");
                scanf("%lf",&l);
                printf("Input centre-line distane between conductors[m]: ");
                scanf("%lf",&a);
                printf("Choose conductor type:\n * slack conductors = l \n * strained conductors = n:\n ");
                printf("Conductor type: ");
                scanf("%c",&conductorType);
                scanf("%c",&conductorType);
                printf("Input number of sub-conductors ");
                scanf("%d",&numberOfSubConductors);
                if(numberOfSubConductors>1)
                {
                    printf("Input center-line distance[m] between sub-conductors: ");
                    scanf("%lf",&ac);
                    printf("Input sub-conductor outer diameter[mm]: ");
                    scanf("%lf",&d);
                }

                if(conductorType=='l')
                {
                    printf("Input extend of one head armature and clamp[m]: ");
                    scanf("%lf",&lh);
                    printf("Input form factor[m]: ");
                    scanf("%lf",&lf);

                    lc=l-2.0*(lh+lf);
                    l=lc;
                }
                if(conductorType=='n')
                {
                    printf("Input insulator chain length[m]: ");
                    scanf("%lf",&li);
                    lc=l-2.0*li;
                }

                printf("Input short-circuit current[kA]: ");
                scanf("%lf",&i);
                printf("Input simulation temperatures[degrees of C] in format[lowerTemperature upperTemperature]: ");
                scanf("%lf %lf",&lowerTemperature,&upperTemperature);
                puts("Mechanical data:");
                printf("Input values of static tensile forces[kN] affecting on conductor in temperature of %0.lf C and temperature of %0.lf C: ",lowerTemperature,upperTemperature);
                scanf("%lf %lf",&staticTensileForceOfConductorAtLowerTemperature,&staticTensileForceOfConductorAtUpperTemperature);
                printf("Input Young module[N/mm2]: ");
                scanf("%lf",&E);
                printf("Input resultant spring constant of both span supports of one span[N/mm]: ");
                scanf("%lf",&S);
                printf("Input mass per unit length of one sub-conductor[kg/m]: ");
                scanf("%lf",&massPerUnitLengthOfOneSubConductor);
                printf("Input number of spacers: ");
                scanf("%d",&numberOfSpacers);
                m=massPerUnitLengthOfOneSubConductor;

                if(numberOfSpacers!=0)
                {
                   ls=0.0;
                   printf("Input mass of one connection[kg]: ");
                   scanf("%lf",&massOfOneConnection);
                   printf("Input mass of one spacer[kg]: ");
                   scanf("%lf",&massOfOneSpacer);
                   m=massPerUnitLengthOfOneSubConductor+((numberOfSpacers-1)*massOfOneConnection+massOfOneSpacer)/(numberOfSubConductors*lc);
                   printf("Input distances[m] between connections[%d] from first insulator to first connection, between next connections and from last connection to last insulator:  ",numberOfSpacers+1);
                   for(it=0;it<=numberOfSpacers;it++)
                   {
                       scanf("%lf",&lk);
                       ls+=lk;
                   }
                   ls=ls/(numberOfSpacers+1);
                }

                printf("Conductor Cross-section[mm2]: ");
                scanf("%lf",&SubConductorCrossSection);
                printf("Choose conductor material: \n * copper['c']\n * Cross-Section ratio Al/Steel>6['a']\n * Cross-section ratio Al\Steel<=6['s']\n");
                printf("Chosen material: ");
                scanf("%c",&conductorMaterialType);
                scanf("%c",&conductorMaterialType);


                if(conductorMaterialType=='c')
                {
                    materialConstant=0.088;
                }
                if(conductorMaterialType=='a')
                {
                    materialConstant=0.27;
                }
                if(conductorMaterialType=='s')
                {
                    materialConstant=0.17;
                }


                printf("Dropper in the middle of the span?[1=Yes , 0=No]: ");
                scanf("%d",&dropper);

                if(dropper==1)
                {
                    printf("Choose arrangement: [parallel='r', perpendicular='p']: ");
                    scanf("%c",&uklad);
                    scanf("%c",&uklad);
                    printf("Input height of the dropper at maximum operating temperature of %.0lf C: ",upperTemperature);
                    scanf("%lf",&heightOfTheDropperAtUpperTemperture);
                    printf("Input width of the dropper at maximum operating temperature of %.0lf C: ",upperTemperature );
                    scanf("%lf",&widthOfTheDropperAtUpperTemperture);
                    printf("Input cord length of dropper at temperature of %.0lf C: ",upperTemperature);
                    scanf("%lf",&dropperCordLength);
                }
                printf("Is the current have to flow along half of the main conductor and along the dropper?[1=Yes, 0=No]:");
                scanf("%d",&zwod);
                if(staticTensileForceOfConductorAtLowerTemperature/(numberOfSubConductors*SubConductorCrossSection)*1000.0>50.0)
                {
                    Emin=E;
                }
                else
                    Emin=E*(0.3+0.7*sin(31.4159265359*staticTensileForceOfConductorAtLowerTemperature/(numberOfSubConductors*SubConductorCrossSection)));

                if(staticTensileForceOfConductorAtUpperTemperature/(numberOfSubConductors*SubConductorCrossSection)*1000.0>50.0)
                {
                    Emax=E;
                }
                else
                    Emax=E*(0.3+0.7*sin(31.4159265359*staticTensileForceOfConductorAtUpperTemperature/(numberOfSubConductors*SubConductorCrossSection)));

                Nmin=1.0/(S*l*1000.0)+1.0/(numberOfSubConductors*Emin*SubConductorCrossSection);
                Nmax=1.0/(S*l*1000.0)+1.0/(numberOfSubConductors*Emax*SubConductorCrossSection);

                stresmin=0.0000000040098375*numberOfSubConductors*numberOfSubConductors*m*m*l*l/(staticTensileForceOfConductorAtLowerTemperature*staticTensileForceOfConductorAtLowerTemperature*staticTensileForceOfConductorAtLowerTemperature*Nmin);
                stresmax=0.0000000040098375*numberOfSubConductors*numberOfSubConductors*m*m*l*l/(staticTensileForceOfConductorAtUpperTemperature*staticTensileForceOfConductorAtUpperTemperature*staticTensileForceOfConductorAtUpperTemperature*Nmax);
                fmin=0.00122625*numberOfSubConductors*m*l*l/staticTensileForceOfConductorAtLowerTemperature;
                fmax=0.00122625*numberOfSubConductors*m*l*l/staticTensileForceOfConductorAtUpperTemperature;

                Timin=1.79428058619*sqrt(fmin);
                Timax=1.79428058619*sqrt(fmax);
                pinchForceAtLowerTemperature=0.0;
                 pinchForceAtUpperTemperature=0.0;

                 if(numberOfSubConductors!=1)
                 {


                if(((1000.0*ac/d<=2.0)&&(ls>=50.0*ac))||((1000.0*ac/d<=2.5)&&(ls>=70.0*ac)))
                {
                    k=0.0;
                }

                else
                {
                    printf("Input factor for the calculation of the peak short-circuit current(depends on the ratio of R/X): ");
                    scanf("%lf",&k);
                    sinus=sin(3.14159265359/numberOfSubConductors);
                    v1kwf=(ac-d/1000.0)*massPerUnitLengthOfOneSubConductor/(0.2*sinus*sinus)*numberOfSubConductors*numberOfSubConductors*ac/(i*i*(numberOfSubConductors-1));
                    tal=-3.0/twoPiXF/log((k-1.02)/0.98);

                    while(Tpi*Tpi+Tpi/(1.0+twoPiXF*twoPiXF*tal*tal)*((twoPiXF*twoPiXF*tal*tal-1.0)*(tal+sin(2.0*twoPiXF*Tpi)/(2.0*twoPiXF))+tal*cos(2.0*twoPiXF*Tpi)-tal*tal*twoPiXF*exp(-Tpi/tal)*(tal*twoPiXF*exp(-Tpi/tal)+4.0*sin(twoPiXF*Tpi)))-v1kwf<0)
                    {
                        Tpi+=0.00001;
                    }
                    v2=v1kwf/(Tpi*Tpi);
                    v3=d*sqrt(ac*1000.0/d-1.0)/(1000.0*ac*sinus*atan(sqrt(ac*1000.0/d-1.0)));
                    shortCircuitCurrentBundleForce=(numberOfSubConductors-1)*0.2*i*i*ls*v2/(numberOfSubConductors*numberOfSubConductors*ac*v3);

                    staticStrainFactorAtLowerTemperature=1500*staticTensileForceOfConductorAtLowerTemperature*ls*ls*Nmin*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0));
                    pinchStrainFactorAtLowerTemperature=0.375*numberOfSubConductors*shortCircuitCurrentBundleForce*ls*ls*ls*Nmin*sinus*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0)*(ac-d/1000.0));
                    staticStrainFactorAtUpperTemperature=1500*staticTensileForceOfConductorAtUpperTemperature*ls*ls*Nmax*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0));
                    pinchStrainFactorAtUpperTemperature=0.375*numberOfSubConductors*shortCircuitCurrentBundleForce*ls*ls*ls*Nmax*sinus*sinus*sinus/((ac-d/1000.0)*(ac-d/1000.0)*(ac-d/1000.0));

                    jmin=sqrt(pinchStrainFactorAtLowerTemperature/(1.0+staticStrainFactorAtLowerTemperature));
                    jmax=sqrt(pinchStrainFactorAtUpperTemperature/(1.0+staticStrainFactorAtUpperTemperature));
                    zetamin=pow(jmin,0.6666);
                    zetamax=pow(jmax,0.6666);


                    if(jmin<1.0)
                    {
                        while(etamin*etamin*etamin+etamin*staticStrainFactorAtLowerTemperature-jmin*jmin*(1.0+staticStrainFactorAtLowerTemperature)*v3*sinus*atan(etamin*(1.0-d/(1000.0*ac))/(1.0-etamin*(1.0-d/(1000.0*ac))))<0.0)
                        {
                            etamin+=0.0001;
                        }

                        v4min=etamin*(ac-d/1000.0)/(ac-etamin*(ac-d/1000.0));
                        vemin=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*i*i/numberOfSubConductors*Nmin*v2*pow(ls*sinus/(etamin*(ac-d/1000.0)),4.0)*etamin*(1.0-atan(sqrt(v4min))/sqrt(v4min))-0.25);
                        pinchForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+vemin/staticStrainFactorAtLowerTemperature*etamin*etamin);

                    }

                    else
                        {
                            while(zetamin*zetamin*zetamin+zetamin*zetamin*staticStrainFactorAtLowerTemperature-jmin*jmin*(1.0+staticStrainFactorAtLowerTemperature)<0.0)
                            {
                                zetamin+=0.0001;
                            }

                        v4min=(ac*1000.0-d)/d;
                        vemin=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*i*i/numberOfSubConductors*Nmin*v2*pow(ls*sinus/(zetamin*(ac-d/1000.0)),4.0)*zetamin*(1.0-atan(sqrt(v4min))/sqrt(v4min))-0.25);
                        pinchForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+vemin*zetamin/staticStrainFactorAtLowerTemperature);
                        }
                     if(jmax<1.0)
                    {
                        while(etamax*etamax*etamax+etamax*staticStrainFactorAtUpperTemperature-jmax*jmax*(1.0+staticStrainFactorAtUpperTemperature)*v3*sinus*atan(etamax*(1.0-d/(1000.0*ac))/(1.0-etamax*(1.0-d/(1000.0*ac))))<0.0)
                        {
                            etamax+=0.0001;
                        }

                        v4max=etamax*(ac-d/1000.0)/(ac-etamax*(ac-d/1000.0));
                        vemax=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*i*i/numberOfSubConductors*Nmax*v2*pow(ls*sinus/(etamax*(ac-d/1000.0)),4.0)*etamax*(1.0-atan(sqrt(v4max))/sqrt(v4max))-0.25);
                        pinchForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+vemax/staticStrainFactorAtUpperTemperature*etamax*etamax);
                    }

                    else
                        {
                            while(zetamax*zetamax*zetamax+zetamax*zetamax*staticStrainFactorAtUpperTemperature-jmax*jmax*(1.0+staticStrainFactorAtUpperTemperature)<0.0)
                            {
                                zetamax+=0.0001;
                            }

                            v4max=(ac*1000.0-d)/d;
                            vemax=0.5+sqrt(1.125*(numberOfSubConductors-1)*0.2*i*i/numberOfSubConductors*Nmax*v2*pow(ls*sinus/(zetamax*(ac-d/1000.0)),4.0)*zetamax*(1.0-atan(sqrt(v4max))/sqrt(v4max))-0.25);
                            pinchForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+vemax*zetamax/staticStrainFactorAtUpperTemperature);
                        }

                    }
                 }
                    r=oneAndAHalfDivByTenDivByStandardGravity*i*i*lc/(numberOfSubConductors*m*a*l);
                        if(zwod==1)
                            r=oneAndAHalfDivByTenDivByStandardGravity*i*i*(lc+dropperCordLength)/(2.0*numberOfSubConductors*m*a*l);

                        if(lineSystem==2)
                        {
                            r=r/0.75;
                        }
                    fi1=atan(r);

                    Trmin=Timin/(pow(1.0+r*r,0.25)*(1.0-fi1*fi1/16.0));
                    Trmax=Timax/(pow(1.0+r*r,0.25)*(1.0-fi1*fi1/16.0));

                plik=fopen("dane/danet.txt","w");
                for(num=amin;num<=amax;num=num+da)
                {
                    // Calculating end swing-out angle
                    if(num/Trmin>0.5)
                    {
                        fiendmin=2*fi1;
                    }
                    else
                        fiendmin=fi1*(1.0-cos(num/Trmin*6.28318530718));

                     if(num/Trmax>0.5)
                    {
                        fiendmax=2*fi1;
                    }
                    else
                        fiendmax=fi1*(1.0-cos(num/Trmax*6.28318530718));
                     //////////////////////////////////////////////

                    //Calculating X parameter
                    if(fiendmin>1.57079632679)
                    {
                        Xmin=1.0-r;
                    }
                    else
                        Xmin=1.0-r*sin(fiendmin);

                    if(fiendmax>1.57079632679)
                    {
                        Xmax=1.0-r;
                    }
                    else
                        Xmax=1.0-r*sin(fiendmax);
                    //////////////////////////////

                    //Calculating maximum swing-out angle
                    if(Xmin<-0.985)
                    {
                        maximumSwingOutAngleAtLowerTemperature=3.14159265359;
                    }
                    else if(Xmin<=0.766)
                    {
                        maximumSwingOutAngleAtLowerTemperature=0.17453292519+acos(Xmin);
                    }
                    else
                        maximumSwingOutAngleAtLowerTemperature=1.25*acos(Xmin);


                    if(Xmax<-0.985)
                    {
                        maximumSwingOutAngleAtUpperTemperature=3.14159265359;
                    }
                    else if(Xmax<=0.766)
                    {
                        maximumSwingOutAngleAtUpperTemperature=0.17453292519+acos(Xmax);
                    }
                    else
                        maximumSwingOutAngleAtUpperTemperature=1.25*acos(Xmax);
                     ///////////////////////////

                    //Calculating load parameter
                    if(num<Trmax/4.0)
                    {
                        loadParameterAtUpperTemperature=3.0*(r*sin(fiendmax)+cos(fiendmax)-1.0);
                    }
                    else
                        loadParameterAtUpperTemperature=3.0*(sqrt(1.0+r*r)-1.0);

                    if(num<Trmin/4.0)
                    {
                        loadParameterAtLowerTemperature=3.0*(r*sin(fiendmin)+cos(fiendmin)-1.0);
                    }
                    else
                        loadParameterAtLowerTemperature=3.0*(sqrt(1.0+r*r)-1.0);

                 //////////////////////////////////////

                 //Calculating 0-place parameter

                  row=-stresmin*(2.0+loadParameterAtLowerTemperature);
                  if(row>0.0)
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row>0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row<0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }

                  row=-stresmax*(2.0+loadParameterAtUpperTemperature);

                  if(row>0.0)
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax*loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row>0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax+loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row<0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }

                  ////////////////////////////////////////////////////////////////

                  //Calculating Dynamic maximum tensile force

                    shortCircuitTensileForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+loadParameterAtLowerTemperature*parmin);
                    shortCircuitTensileForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+loadParameterAtUpperTemperature*parmax);


                  //Calculating maximum dynamic span
                    elasticExpansionAtUpperTemperature=Nmax*(shortCircuitTensileForceAtUpperTemperature-staticTensileForceOfConductorAtUpperTemperature)*1000.0;
                    elasticExpansionAtLowerTemperature=Nmin*(shortCircuitTensileForceAtLowerTemperature-staticTensileForceOfConductorAtLowerTemperature)*1000.0;


                    if(num<Trmax/4.0)
                    {
                        thermalExpansionAtUpperTemperature=materialConstant*i*i/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*num;
                    }
                    else
                        thermalExpansionAtUpperTemperature=0.25*materialConstant*i*i/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Trmax;


                    if(num<Trmin/4.0)
                    {
                        thermalExpansionAtLowerTemperature=materialConstant*i*i/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*num;
                    }
                    else
                        thermalExpansionAtLowerTemperature=0.25*materialConstant*i*i/(numberOfSubConductors*numberOfSubConductors*SubConductorCrossSection*SubConductorCrossSection)*Trmin;



                    dilatationFactorAtUpperTemperature=sqrt(1.0+0.375*l*l/(fmax*fmax)*(elasticExpansionAtUpperTemperature+thermalExpansionAtUpperTemperature));
                    dilatationFactorAtLowerTemperature=sqrt(1.0+0.375*l*l/(fmin*fmin)*(elasticExpansionAtLowerTemperature+thermalExpansionAtLowerTemperature));
                    if(r>=1.8)
                    {
                        formFactor=1.15;
                    }
                    else if(r>0.8)
                    {
                        formFactor=0.97+0.1*r;
                    }
                    else
                        {formFactor=1.05;}

                    fedmax=formFactor*dilatationFactorAtUpperTemperature*fmax;
                    fedmin=formFactor*dilatationFactorAtLowerTemperature*fmin;

                    //////////////////////////////

                    //calculating tensile force caused by swing-out

                    if(uklad=='r')
                    {
                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax+fedmin)*(heightOfTheDropperAtUpperTemperture+fmin+fedmin)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture))
                        {
                            deltamin=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmin*fedmin-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmin*(heightOfTheDropperAtUpperTemperture+fmax)));
                        }

                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax+fedmax)*(heightOfTheDropperAtUpperTemperture+fmax+fedmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture))
                        {
                            deltamax=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmax*fedmax-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmax*(heightOfTheDropperAtUpperTemperture+fmax)));
                        }
                    }

                    if(uklad=='p')
                    {
                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)+fedmin)
                        {
                            deltamin=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmin*fedmin-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmin*sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)))+acos((heightOfTheDropperAtUpperTemperture+fmax)/sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture));
                        }

                        if(dropperCordLength<sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)+fedmax)
                        {
                            deltamax=acos(((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+fedmax*fedmax-dropperCordLength*dropperCordLength+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)/(2.0*fedmax*sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture)))+acos((heightOfTheDropperAtUpperTemperture+fmax)/sqrt((heightOfTheDropperAtUpperTemperture+fmax)*(heightOfTheDropperAtUpperTemperture+fmax)+widthOfTheDropperAtUpperTemperture*widthOfTheDropperAtUpperTemperture));
                        }
                    }

                    if(deltamin>=0)
                    {
                        if(deltamin>=fi1)
                        {
                            if(num<Trmin/4.0)
                            {
                                loadParameterAtLowerTemperature=3.0*(r*sin(fiendmin)+cos(fiendmin)-1.0);
                            }
                            else
                                loadParameterAtLowerTemperature=3.0*(sqrt(1.0+r*r)-1.0);
                        }
                        else
                        {
                            if(deltamin>fiendmin)
                            {
                                loadParameterAtLowerTemperature=3.0*(r*sin(fiendmin)+cos(fiendmin)-1.0);
                            }
                            else
                            {
                                loadParameterAtLowerTemperature=3.0*(r*sin(deltamin)+cos(deltamin)-1.0);
                            }
                        }

                        parmin=0.0;
                        row=-stresmin*(2.0+loadParameterAtLowerTemperature);
                  if(row>0.0)
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row>0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtLowerTemperature*loadParameterAtLowerTemperature*parmin*parmin*parmin+loadParameterAtLowerTemperature*(2.0+stresmin)*parmin*parmin+(1.0+2.0*stresmin)*parmin+row<0.0)
                      {
                          parmin=parmin+0.001;
                      }
                  }


                     shortCircuitTensileForceAtLowerTemperature=staticTensileForceOfConductorAtLowerTemperature*(1.0+loadParameterAtLowerTemperature*parmin);

                    }


                    if(deltamax>=0)
                    {
                        if(deltamax>=fi1)
                        {
                            if(num<Trmax/4.0)
                            {
                                loadParameterAtUpperTemperature=3.0*(r*sin(fiendmax)+cos(fiendmax)-1.0);
                            }
                            else
                                loadParameterAtUpperTemperature=3.0*(sqrt(1.0+r*r)-1.0);
                        }
                        else
                        {
                            if(deltamax>fiendmax)
                            {
                                loadParameterAtUpperTemperature=3.0*(r*sin(fiendmax)+cos(fiendmax)-1.0);
                            }
                            else
                            {
                                loadParameterAtUpperTemperature=3.0*(r*sin(deltamax)+cos(deltamax)-1.0);
                            }
                        }


                     parmax=0.0;
                  row=-stresmax*(2.0+loadParameterAtUpperTemperature);

                  if(row>0.0)
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax+loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row>0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }
                  else
                  {
                      while(loadParameterAtUpperTemperature*loadParameterAtUpperTemperature*parmax*parmax*parmax+loadParameterAtUpperTemperature*(2.0+stresmax)*parmax*parmax+(1.0+2.0*stresmax)*parmax+row<0.0)
                      {
                          parmax=parmax+0.001;
                      }
                  }

                     shortCircuitTensileForceAtUpperTemperature=staticTensileForceOfConductorAtUpperTemperature*(1.0+loadParameterAtUpperTemperature*parmax);


                    }


                    //Calculating a force caused by drop force
                    dropForceAtLowerTemperature=0.0;
                    dropForceAtUpperTemperature=0.0;

                    if((r>0.6)&&(maximumSwingOutAngleAtLowerTemperature>=1.221730476))
                    {
                        if(dropper==0)
                        {
                            dropForceAtLowerTemperature=1.2*staticTensileForceOfConductorAtLowerTemperature*sqrt(1.0+2.54647908947*stresmin*maximumSwingOutAngleAtLowerTemperature);
                        }
                        if((dropper==1)&&(deltamin>=1.047197551))
                        {
                            dropForceAtLowerTemperature=1.2*staticTensileForceOfConductorAtLowerTemperature*sqrt(1.0+2.54647908947*stresmin*maximumSwingOutAngleAtLowerTemperature);
                        }
                    }

                    if((r>0.6)&&(maximumSwingOutAngleAtUpperTemperature>=1.221730476))
                    {
                        if(dropper==0)
                        {
                            dropForceAtUpperTemperature=1.2*staticTensileForceOfConductorAtUpperTemperature*sqrt(1.0+2.54647908947*stresmax*maximumSwingOutAngleAtUpperTemperature);
                        }
                        if((dropper==1)&&(deltamax>=1.047197551))
                        {
                            dropForceAtUpperTemperature=1.2*staticTensileForceOfConductorAtUpperTemperature*sqrt(1.0+2.54647908947*stresmax*maximumSwingOutAngleAtUpperTemperature);
                        }
                    }



                    if(conductorType=='l')
                    {
                        if(maximumSwingOutAngleAtLowerTemperature<1.57079632679)
                        {
                            bhmin=fedmin*sin(maximumSwingOutAngleAtLowerTemperature);
                        }
                        else
                            bhmin=fedmin;

                        if(maximumSwingOutAngleAtUpperTemperature<1.57079632679)
                        {
                            bhmax=fedmax*sin(maximumSwingOutAngleAtUpperTemperature);
                        }
                        else
                            bhmax=fedmax;
                    }
                    if(conductorType=='n')
                    {
                        if(dropper==0)
                        {
                            if(maximumSwingOutAngleAtLowerTemperature<fi1)
                            {
                                bhmin=fedmin*sin(maximumSwingOutAngleAtLowerTemperature);
                            }
                            else
                                bhmin=fedmin*sin(fi1);

                            if(maximumSwingOutAngleAtUpperTemperature<fi1)
                            {
                                bhmax=fedmax*sin(maximumSwingOutAngleAtUpperTemperature);
                            }
                            else
                                bhmax=fedmax*sin(fi1);
                        }
                        if(dropper==1)
                        {
                            if(deltamin<maximumSwingOutAngleAtLowerTemperature)
                            {
                                if(maximumSwingOutAngleAtLowerTemperature<fi1)
                                {
                                    bhmin=fedmin*sin(maximumSwingOutAngleAtLowerTemperature);
                                }
                                else
                                    bhmin=fedmin*sin(fi1);
                            }
                            else
                            {
                                if(deltamin<fi1)
                                {
                                    bhmin=fedmin*sin(deltamin);
                                }
                                else
                                    bhmin=fedmin*sin(fi1);
                            }


                            if(deltamax<maximumSwingOutAngleAtUpperTemperature)
                            {
                                if(deltamax<fi1)
                                {
                                    bhmax=fedmax*sin(deltamax);
                                }
                                else
                                    bhmax=fedmax*sin(fi1);
                            }
                            else
                            {

                                if(maximumSwingOutAngleAtUpperTemperature<fi1)
                                {
                                    bhmax=fedmax*sin(maximumSwingOutAngleAtUpperTemperature);
                                }
                                else
                                    bhmax=fedmax*sin(fi1);
                            }

                        }
                    }


                    horizontalSpanDisplacement=(bhmin>bhmax) ? bhmin : bhmax;
                    minimumAirClearance=a-2.0*horizontalSpanDisplacement;


                    if(numberOfSubConductors!=1)
                    {
                        if(k==0.0)
                        {
                            pinchForceAtLowerTemperature=1.1*shortCircuitTensileForceAtLowerTemperature;
                            pinchForceAtUpperTemperature=1.1*shortCircuitTensileForceAtUpperTemperature;
                        }
                    }

                    maximumShortCircuitTensileForce=(shortCircuitTensileForceAtLowerTemperature>=shortCircuitTensileForceAtUpperTemperature)? shortCircuitTensileForceAtLowerTemperature : shortCircuitTensileForceAtUpperTemperature;
                    maximumDropForce=(dropForceAtLowerTemperature>=dropForceAtUpperTemperature)? dropForceAtLowerTemperature : dropForceAtUpperTemperature;
                    maximumPinchForce=(pinchForceAtLowerTemperature>=pinchForceAtUpperTemperature)? pinchForceAtLowerTemperature : pinchForceAtUpperTemperature;
                    if(conductorType=='l')
                        maximumShortCircuitTensileForce*=1.5;
                    maks=maximumShortCircuitTensileForce;
                    if(maks<maximumDropForce)
                    {
                        maks=maximumDropForce;
                    }
                    if(maks<maximumPinchForce)
                    {
                        maks=maximumPinchForce;
                    }

                if(amin+da>amax)
                {
                printf("Additional parameters: \n");
                    printf("    r= %.2lf\n",r);
                    printf("    Maximum Swing-Out Angle At Lower Temperature = %lf , Maximum Swing-Out Angle At Upper Temperature = %lf\n",maximumSwingOutAngleAtLowerTemperature,maximumSwingOutAngleAtUpperTemperature);
                    printf("    deltamin = %lf , deltamax = %lf\n",deltamin,deltamax);
                    printf("    parmin: %.3lf , parmax: %.3lf\n",parmin,parmax);
                    printf("    loadParameterAtLowerTemperature: %.2lf ,loadParameterAtUpperTemperature: %.2lf\n",loadParameterAtLowerTemperature,loadParameterAtUpperTemperature);
                    printf("    stresmin: %.2lf , stresmax: %.2lf\n",stresmin,stresmax);
                    printf("    fesmin: %.2lf , fesmax: %.2lf\n",fmin,fmax);
                    printf("    fedmin: %.2lf , fedmax: %.2lf\n",fedmin,fedmax);
                    printf("    Elastic Expansion At Lower Temperature: %.3lf , Elastic Expansion At Upper Temperature: %.3lf\n",elasticExpansionAtLowerTemperature*1000.0,elasticExpansionAtUpperTemperature*1000.0);
                    printf("    Thermal Expansion At Lower Temperature: %.3lf , Thermal Expansion At Upper Temperature: %.3lf\n",thermalExpansionAtLowerTemperature*10000.0,thermalExpansionAtUpperTemperature*10000.0);
                    printf("    Material constant: %.3lf\n",materialConstant);
                    printf("    Dilatation Factor At Lower Temperature: %.3lf , Dilatation Factor At Upper Temperature: %.3lf\n",dilatationFactorAtLowerTemperature,dilatationFactorAtUpperTemperature);
                    printf("    formFactor: %.3lf\n",formFactor);
                puts("Forces:");
                    printf("    Maximum Short-Circuit Tensile Force: %.2lf kN\n",maximumShortCircuitTensileForce);
                    printf("    Maximum Drop Force: %.2lf kN\n",maximumDropForce);
                    printf("    Maximum Pinch Force: %.2lf kN\n",maximumPinchForce);
                    printf("    Overall Maximum Force: %.2lf kN\n",maks);
                puts("Effects of the force");
                    printf("    Horizontal span displacement: %.2lf m\n",horizontalSpanDisplacement);
                    printf("    Minimum air clearance: %.2lf m\n",minimumAirClearance);

                }
                else
                    fprintf(plik,"%lf %lf %lf %lf %lf %lf %lf\n",num,maximumShortCircuitTensileForce,maximumDropForce,maximumPinchForce,maks,horizontalSpanDisplacement,minimumAirClearance);
                }
                fclose(plik);
                if(amin+da<=amax)
                    ShellExecute(NULL, NULL, "dane\\danet.m", NULL, NULL, SW_SHOWNORMAL);
                system("pause");
                system("cls");
                puts("Choose a parameter in terms of calculations will be performed:\n *'a'-center line distance between main conductor mid-points; \n *'i' - three phase initial symmetrical short-circuit current'; \n *'t'-duration of short-circuit current");
            }
            break;
        }


    }
    return 0;
}
