#define Nmax_x	10000

#define DD	1
#define DF	2
#define FF	3
#define FD	4
#define DG	5
#define DFD	6
#define DDF	7
#define DFFi	8
#define DFFo	9
#define	LGL	10
#define	DGG	11

#define PONTSZAM 10000

#define KEZD_X	10
#define KEZD_Y	10
#define DOT_X	2
#define DOT_Y	2
//#include "/afs/nd.edu/user1/eravasz/hasznos/includes/Epslib_withoutX.c"

struct line_param{
	double slope,inc,cor;
    double slope_SD,inc_SD;
} plnow;

line_param *fit_line(unsigned long long NNNN,double xxxx[],double yyyy[]){
    unsigned int ijk;
    double szum_x=0,szum_y=0,szum_xx=0,szum_yy=0,szum_xy=0,szum_e2=0,szum_xres2=0;
    
    line_param *plnow=NULL;
    plnow=new line_param;
    //printf("Npoints=%d\n",NNNN);getchar();
    for(ijk=1;ijk<=NNNN;ijk++){
         szum_x+=xxxx[ijk];
         szum_y+=yyyy[ijk];
         szum_xx+=(xxxx[ijk])*(xxxx[ijk]);
         szum_yy+=(yyyy[ijk])*(yyyy[ijk]);
         szum_xy+=(xxxx[ijk])*(yyyy[ijk]);
     }
     //printf("szum_x=%lf,szum_xy=%lf\n",szum_x,szum_y);getchar();
    if(NNNN*szum_xx-szum_x*szum_x!=0){
        plnow->slope=(NNNN*szum_xy-szum_x*szum_y)/(NNNN*szum_xx-szum_x*szum_x);
        plnow->inc=(szum_y-plnow->slope*szum_x)/(double)NNNN;
        //printf("%lf  X  +  %lf\n",p->slope,p->inc);getchar();
        if(NNNN*szum_yy-szum_y*szum_y!=0)
           plnow->cor=(NNNN*szum_xy-szum_x*szum_y)/
                        sqrt((NNNN*szum_xx-szum_x*szum_x)*(NNNN*szum_yy-szum_y*szum_y));
        else {plnow->cor=1.;}
        //printf("Y = %lg  X  +  %lg\n",p->slope,p->inc);//getchar();
         for(ijk=1;ijk<=NNNN;ijk++) {
             szum_e2+= (yyyy[ijk]-plnow->inc-plnow->slope*xxxx[ijk])
                      *(yyyy[ijk]-plnow->inc-plnow->slope*xxxx[ijk]);
             szum_xres2+=(xxxx[ijk]-szum_x/(double)NNNN)*(xxxx[ijk]-szum_x/(double)NNNN);
         }
        plnow->slope_SD=sqrt((szum_e2/(double)(NNNN-2.)) / szum_xres2);
        plnow->inc_SD=sqrt((szum_e2/(double)(NNNN*(NNNN-2.))) * szum_xx / szum_xres2);
        return(plnow);
    }
    else  {printf("infinite slope!\n\n");
        return(NULL);
    }
}

void fit_line_now(unsigned long long NNNN,double xxxx[],double yyyy[]){
    unsigned int ijk;
    double szum_x=0,szum_y=0,szum_xx=0,szum_yy=0,szum_xy=0,szum_e2=0,szum_xres2=0;
    
    //printf("Npoints=%d\n",NNNN);getchar();
    for(ijk=1;ijk<=NNNN;ijk++){
        szum_x+=xxxx[ijk];
        szum_y+=yyyy[ijk];
        szum_xx+=(xxxx[ijk])*(xxxx[ijk]);
        szum_yy+=(yyyy[ijk])*(yyyy[ijk]);
        szum_xy+=(xxxx[ijk])*(yyyy[ijk]);
    }
    //printf("szum_x=%lf,szum_xy=%lf\n",szum_x,szum_y);getchar();
    if(NNNN*szum_xx-szum_x*szum_x!=0){
        plnow.slope=(NNNN*szum_xy-szum_x*szum_y)/(NNNN*szum_xx-szum_x*szum_x);
        plnow.inc=(szum_y-plnow.slope*szum_x)/(double)NNNN;
        //printf("%lf  X  +  %lf\n",p->slope,p->inc);getchar();
        if(NNNN*szum_yy-szum_y*szum_y!=0)
            plnow.cor=(NNNN*szum_xy-szum_x*szum_y)/
            sqrt((NNNN*szum_xx-szum_x*szum_x)*(NNNN*szum_yy-szum_y*szum_y));
        else {plnow.cor=1.;}
        //printf("Y = %lg  X  +  %lg\n",p->slope,p->inc);//getchar();
        for(ijk=1;ijk<=NNNN;ijk++) {
            szum_e2+= (yyyy[ijk]-plnow.inc-plnow.slope*xxxx[ijk])
            *(yyyy[ijk]-plnow.inc-plnow.slope*xxxx[ijk]);
            szum_xres2+=(xxxx[ijk]-szum_x/(double)NNNN)*(xxxx[ijk]-szum_x/(double)NNNN);
        }
        plnow.slope_SD=sqrt((szum_e2/(double)(NNNN-2.)) / szum_xres2);
        plnow.inc_SD=sqrt((szum_e2/(double)(NNNN*(NNNN-2.))) * szum_xx / szum_xres2);
    }
  //  else  {printf("infinite slope!\n\n");}
}

void fit_xy_line(const char pth[],const char fnm[], int N){
    double *x,*y;
    line_param *lp;
    FILE *f;
    char fn[300];
    sprintf(fn,"%s/%s",pth,fnm);
    f=fopen(fn,"r");
    
    x=new double[N+1];
    y=new double[N+1];
    
    for (int i=1; i<=N; i++)
        fscanf(f,"%lg\t%lg",&(x[i]),&(y[i]));
    
    lp=fit_line(N,x,y);
    
    printf("%s:\n\t",fnm);
    printf("y = (%lg +- %lg) * x + (%lg +- %lg)   c=%lg\n",lp->slope,lp->slope_SD, lp->inc,lp->inc_SD,lp->cor);
    
}

void fit_log10_x_log10_y_line(const char pth[],const char fnm[],double x_int[],double y_int[]){
    double *x,*y,a,b;
    FILE *f;
    
    char fn[300];
    sprintf(fn,"%s/%s",pth,fnm);
    f=fopen(fn,"r");
    
    int N=0;
    
    
    int vege=1;
    while (vege!=EOF) {
        vege=fscanf(f,"%lg\t%lg",&a,&b);
        if (vege!=EOF){
            if ((x_int==NULL)&&(y_int==NULL)&&(a>0)&&(b>0)) N++;
            else {
                int ok=1;
                if ((x_int!=NULL)&&((x_int[0]>a)||(x_int[1]<a))&&(a>0)) ok=0;
                if ((y_int!=NULL)&&((y_int[0]>b)||(y_int[1]<b))&&(b>0)) ok=0;
                if (ok)  N++;
            }
        }
    }
    fclose(f);
    x=new double[N+1];
    y=new double[N+1];
    
    f=fopen(fn,"r");
    N=0;
    vege=1;
    while (vege!=EOF) {
        vege=fscanf(f,"%lg\t%lg\n",&a,&b);
        if (vege!=EOF){
            if ((x_int==NULL)&&(y_int==NULL)&&(a>0)&&(b>0)) {
                N++;
                x[N]=log10(a);
                y[N]=log10(b);
                //printf("x=%lg\t y=%lg\n",x[N],y[N]);//getchar();
            }
            else {
                int ok=1;
                if ((x_int!=NULL)&&((x_int[0]>a)||(x_int[1]<a))&&(a>0)) ok=0;
                if ((y_int!=NULL)&&((y_int[0]>b)||(y_int[1]<b))&&(b>0))  ok=0;
                if (ok)  {
                    N++;
                    x[N]=log10(a);
                    y[N]=log10(b);
                    //  printf("x=%lg\t y=%lg\n",x[N],y[N]);//getchar();
                }
            }
        }
    }
    fclose(f);
    fit_line_now(N,x,y);
    
    printf("x=%lg\t y=%lg\n",x[N],y[N]);//getchar();
    printf("%s:\n\tN=%d\n",fnm,N);
    printf("y = (%lg +- %lg) * x + (%lg +- %lg)   c=%lg\n",plnow.slope,plnow.slope_SD, plnow.inc,plnow.inc_SD,plnow.cor);
    
    delete[] x; x=NULL;
    delete[] y; y=NULL;
}


void fit_x_log10_y_line_ignore_smallest_nonzero_y(const char pth[],const char fnm[],double x_int[],double y_int[]){
    double *x,*y,b,sumnow,sumout;
    FILE *f;
   
    char fn[300];
    sprintf(fn,"%s/%s",pth,fnm);
   // printf("%s\n",fn);
    f=fopen(fn,"r");
    
    int N=0;
    double k;
   
    int vege=1;
    vege=fscanf(f,"%lg\t%lg\n",&k,&b);
    double y_min;
    y_min=b;
  //  printf("k=%d\t b=%lg   y_min=%lg\n",k,b,y_min);getchar();
    while (vege!=EOF) {
        // vege=fscanf(f,"%lg\t%lg",&a,&b);
        vege=fscanf(f,"%lg\t%lg\n",&k,&b);
        if (((b<y_min)&&(b>0))) y_min=b;
        //printf("k=%d\t b=%lg   y_min=%lg\n",k,b,y_min);getchar();
        
    }
    fclose(f);
  //  printf("y_min=%lg to ignore!\n",y_min);
    
    f=fopen(fn,"r");
    vege=1;
    while (vege!=EOF) {
       // vege=fscanf(f,"%lg\t%lg",&a,&b);
        vege=fscanf(f,"%lg\t%lg\n",&k,&b);
       // printf("a=%lg\t b=%lg\n",k,b);
        if (vege!=EOF){
            if ((x_int==NULL)&&(y_int==NULL)&&(b>y_min)) N++;
            else {
                int ok=1;
                if ((x_int!=NULL)&&((x_int[0]>k)||(x_int[1]<k))) ok=0;
                if ((y_int!=NULL)&&((y_int[0]>b)||(y_int[1]<b))) ok=0;
                if ((ok)&&(b>y_min))  N++;
            }
        }
    }
    fclose(f);
    x=new double[N+1];
    y=new double[N+1];
    
    f=fopen(fn,"r");
    N=0;sumnow=0;sumout=0;
    vege=1;
    while (vege!=EOF) {
        vege=fscanf(f,"%lg\t%lg\n",&k,&b);
       if (vege!=EOF){
            if ((x_int==NULL)&&(y_int==NULL)&&(b>y_min))  {
                N++;
                x[N]=k;
                y[N]=log10(b);
                //printf("x=%lg\t y=%lg\n",x[N],y[N]);//getchar();
           }
            else {
                int ok=1;
                if ((x_int!=NULL)&&((x_int[0]>k)||(x_int[1]<k))) ok=0;
                if ((y_int!=NULL)&&((y_int[0]>b)||(y_int[1]<b))) ok=0;
                if ((ok)&&(b>y_min))  {
                    N++;
                    x[N]=k;
                    y[N]=log10(b);
                    sumnow+=b;
                   // printf("x=%lg\t y=%lg\n",x[N],y[N]);getchar();
                }
                if (ok==0) {
                    sumout+=b;
                }
            }
        }
    }
    fclose(f);
    fit_line_now(N,x,y);
 
//    printf("x=%lg\t y=%lg\n",x[N],y[N]);//getchar();
    printf("%s:\tN=%d\t",fnm,N);
    //printf("y = (%lg +- %lg) * x + (%lg +- %lg)   c=%lg\n",plnow.slope,plnow.slope_SD, plnow.inc,plnow.inc_SD,plnow.cor);

    //printf("Sum of y in interval = %lg\n",sumnow);

    printf("%lg\t%lg\t%lg\t",-plnow.slope, plnow.inc,plnow.cor);
    
    printf("%lg\n",sumnow);

    
    //    if((plnow.slope<0)&&(plnow.cor<-0.5)&&(N>=10)) printf("%lg\t",-plnow.slope);
//    else printf("0\t");
 //   printf("%lg\t",1.-sumout);
    delete[] x; x=NULL;
    delete[] y; y=NULL;
}

line_param *fit_line_INC_0(unsigned long long NNNN,double xxxx[],double yyyy[])
{unsigned int ijk;
    double szum_x=0,szum_y=0,szum_xx=0,szum_yy=0,szum_xy=0,szum_e2=0,szum_xres2=0;
    line_param *p;
    p=new line_param;
    //printf("Npoints=%d\n",NNNN);getchar();
    for(ijk=1;ijk<=NNNN;ijk++){
        szum_x+=xxxx[ijk];
        szum_y+=yyyy[ijk];
        szum_xx+=(xxxx[ijk])*(xxxx[ijk]);
        szum_yy+=(yyyy[ijk])*(yyyy[ijk]);
        szum_xy+=(xxxx[ijk])*(yyyy[ijk]);
    }
    //printf("szum_x=%lf,szum_xy=%lf\n",szum_x,szum_y);getchar();
    if(NNNN*szum_xx-szum_x*szum_x!=0){
        p->slope=szum_xy/szum_xx;
        p->inc=0.0;
        //printf("%lf  X  +  %lf\n",p->slope,p->inc);getchar();
        for(ijk=1;ijk<=NNNN;ijk++) {
            szum_e2+= (yyyy[ijk]-p->inc-p->slope*xxxx[ijk])
            *(yyyy[ijk]-p->inc-p->slope*xxxx[ijk]);
            szum_xres2+=(xxxx[ijk]-szum_x/(double)NNNN)*(xxxx[ijk]-szum_x/(double)NNNN);
        }
        p->slope_SD=sqrt((szum_e2/(double)(NNNN-2.)) / szum_xres2);
        p->inc_SD=sqrt((szum_e2/(double)(NNNN*(NNNN-2.))) * szum_xx / szum_xres2);
        
        if(NNNN*szum_yy-szum_y*szum_y!=0){
            p->cor=sqrt(1-szum_e2/(szum_yy-szum_y*szum_y/(double)NNNN));
        //            p->cor=(NNNN*szum_xy-szum_x*szum_y)/
        //          sqrt((NNNN*szum_xx-szum_x*szum_x)*(NNNN*szum_yy-szum_y*szum_y));
        }
        else {p->cor=1.;}

        // printf("%lf  X  +  %lf\n",p->slope,p->inc);//getchar();
        return(p);		   
    }
    else  {printf("infinite slope!\n\n");
        return(NULL);
    }
}


line_param *make_binned_PK(char benev[400],unsigned long long PSZ,unsigned int type2)
{
double  y[Nmax_x+2];
double  x[Nmax_x+2];
double  xvon[PONTSZAM+2];
double ertx[PONTSZAM+1];
double erty[PONTSZAM+1];
double xf[PONTSZAM+1];
double yf[PONTSZAM+1];
line_param *p;
FILE *date2;
 char logbinnev[450];
  unsigned long long vege=0,i=1,j,num;
  double a,bla,xmax,x0,x0_real;
  double  nr=0,step=0;
  unsigned long long nri,stepi,pk,PSZf;
  unsigned long long l_nri,l_stepi;
  double atlx=0,atly=0,db=0;

  date2=fopen(benev,"r");
   //printf("%s\n",benev); getchar();
 do{
       if(type2==DD) {vege=fscanf(date2,"%lld%lld",&nri,&stepi);  step=stepi; nr=nri;}
       if(type2==DF) {vege=fscanf(date2,"%lld%lf",&nri,&step); nr=nri;}
       if(type2==FD) {vege=fscanf(date2,"%lf%lld",&nr,&stepi); step=stepi; }
       if(type2==FF) vege=fscanf(date2,"%lf%lf",&nr,&step);
       if(type2==DG) {vege=fscanf(date2,"%lld%lg",&nri,&step); nr=nri;}
	   if(type2==DGG) {vege=fscanf(date2,"%lld%lg%lld",&nri,&step,&pk); nr=nri;}
	   if(type2==DFD) {vege=fscanf(date2,"%lld%lf%lld",&nri,&step,&pk);nr=nri;}
       if(type2==DDF) {vege=fscanf(date2,"%lld%lld%lf",&pk,&nri,&step);nr=nri;}
       if(type2==DFFi) {vege=fscanf(date2,"%lld%lf%lf",&nri,&step,&bla);nr=nri;}
       if(type2==DFFo) {vege=fscanf(date2,"%lld%lf%lf",&nri,&bla,&step);nr=nri;}
       if(type2==LGL) {vege=fscanf(date2,"%lld%lg%lld",&l_nri,&step,&l_stepi);nr=l_nri;}
  y[i]=step;
  x[i]=nr;
  //printf("%d\t%lf\t%lf\n",i,nr,step); getchar();
  i++;
  }while(vege!=EOF);
 fclose(date2);
 xmax=0;  num=i-1; x0=x[num];x0_real=x[1];
 for(i=1;i<=num;i++)
 	{if(x[i]>xmax) xmax=x[i];
	 if((x[i]<x0)&&(x[i]>0.)) x0=x[i];
	 if(x[i]<x0_real) x0_real=x[i];
	}

  if(!PSZ) PSZ=num/10;
  //printf("x0_real=%lf\tx0=%lf\txmax=%lf\n",x0_real,x0,xmax);
  //if(x0==0) getchar();
  a=(log(xmax)-log(x0))/(double)(PSZ);
   //printf("a=%lf\n",a);

  //getchar();
  xvon[0]=x0;
  for(i=1;i<=PSZ;i++)
   {xvon[i]=xvon[i-1]*exp(a);
 //   printf("%lf\n",xvon[i]);
   }

for(i=1;i<=PSZ;i++)
 {atlx=xvon[i-1]*exp(a/2.); 
 //printf("interv: %lf -- %lf -- %lf\n\n",xvon[i-1],atlx,xvon[i]);
  for(j=1;j<=num;j++)
   {if( (xvon[i-1]<=x[j] ) && (x[j]<xvon[i]) )
          {atly+=y[j];
           db++;
	  // printf("j=%ld\ti=%ld\tx[j]=%lf   y[j]=%lf\n",j,i,x[j],y[j]);
	  // printf("\t%lf\t%lf\t%lf\t%lf\n",db,atlx,atly,atly/(xvon[i]-xvon[i-1]));getchar();
          }
   }
  if(db!=0){ ertx[i]=atlx;
             erty[i]=atly*(xmax-x0_real)/((xvon[i]-xvon[i-1])*num);
           }
  else {ertx[i]=0;erty[i]=0;}
atlx=0;
atly=0;
 db=0;
  }

if(type2==DFFi) sprintf(logbinnev,"%s_bin_in",benev);
else if(type2==DFFo) sprintf(logbinnev,"%s_bin_out",benev);
	else sprintf(logbinnev,"%s_bin",benev);
date2=fopen(logbinnev,"w");
PSZf=1;
for(i=1;i<=PSZ;i++)
    if(erty[i]!=0) { fprintf(date2,"%lf\t%lg\n",ertx[i],erty[i]);
					 xf[PSZf]=log(ertx[i]);yf[PSZf]=log(erty[i]);PSZf++;	
	}
fclose(date2);
 p=fit_line(PSZf,xf,yf);
 printf("fit done on %lld points\n",PSZf);//getchar();
 return(p);
}


void make_cumulative_PK(char benev[400],unsigned int type2)
{
double  y[Nmax_x+2];
double  x[Nmax_x+2];
double ertx[PONTSZAM+1];
double erty[PONTSZAM+1];
FILE *date2;
 char logbinnev[700];
unsigned long long vege=0,i,PSZ;
double  nr=0,step=0,bla;
unsigned long long nri,stepi,pk=0;

//  printf("%s\n",benev);getchar();

  for(i=0;i<=Nmax_x;i++)
     {y[i]=0.;x[i]=0.;}
  i=1;
  //printf("%s\n",benev);getchar();
  
  date2=fopen(benev,"r");
  printf("%s\n",benev);
 do{
       if(type2==DD) {vege=fscanf(date2,"%lld%lld",&nri,&stepi);  step=stepi; nr=nri;}
       if(type2==DF) {vege=fscanf(date2,"%lld%lf",&nri,&step); nr=nri;}
       if(type2==FD) {vege=fscanf(date2,"%lf%lld",&nr,&stepi); step=stepi; }
       if(type2==FF) vege=fscanf(date2,"%lf%lf",&nr,&step);
       if(type2==DG) {vege=fscanf(date2,"%lld%lg",&nri,&step); nr=nri;}
       if(type2==DFD) {vege=fscanf(date2,"%lld%lf%lld",&nri,&step,&pk);nr=nri;}
       if(type2==DDF) {vege=fscanf(date2,"%lld%lld%lf",&pk,&nri,&step);nr=nri;}
       if(type2==DFFi) {vege=fscanf(date2,"%lld%lf%lf",&nri,&step,&bla);nr=nri;}
       if(type2==DFFo) {vege=fscanf(date2,"%lld%lf%lf",&nri,&bla,&step);nr=nri;}
       if(type2==LGL) {vege=fscanf(date2,"%lld%lg%lld",&nri,&step,&stepi);nr=nri;}
  y[i]=step;
  x[i]=nr;
  //printf("%d\t%lf\t%lf\n",i,nr,step); getchar();
  i++;
  }while(vege!=EOF);
 fclose(date2);
 
 PSZ=i-2;
 for(i=PSZ+2;i>0;i--)
 	{erty[i]=0.;ertx[i]=0.;}
 for(i=PSZ;i>0;i--)
   	{erty[i]=erty[i+1]+y[i];
	 ertx[i]=x[i];
	}
if(type2==DFFi) sprintf(logbinnev,"%s_cum_in",benev);
else if(type2==DFFo) sprintf(logbinnev,"%s_cum_out",benev);
	else sprintf(logbinnev,"%s_cum",benev);
date2=fopen(logbinnev,"w");
for(i=1;i<=PSZ;i++)
    if(erty[i]!=0)  fprintf(date2,"%lf\t%lg\n",ertx[i],erty[i]);
fclose(date2);
printf("cumul megvan\n");//getchar();
}


double decide_color(unsigned long long hany,double hanymax)
{
printf("%lg\n",hany/(double)hanymax);
return(hany/(double)hanymax);
}


void make_binned_averaging(char benev[400],unsigned long long PSZ,unsigned int type2)
{
  double  y[Nmax_x+2];
  double  x[Nmax_x+2];
  double  xvon[PONTSZAM+2];
  double ertx[PONTSZAM+1];
  double erty[PONTSZAM+1];
  FILE *date2;
   char logbinnev[450];
  unsigned long long vege=0,i=1,j,num;
  double a,xmax,x0;
  double  nr=0,step=0;
  unsigned long long nri,stepi,pk=0;
  double atlx=0,atly=0,db=0;

for(i=1;i<=Nmax_x+1;i++)
 {x[i]=0;y[i]=0;}

i=1;
  date2=fopen(benev,"r");
   printf("%s\n",benev);
 do{
       if(type2==DD) {vege=fscanf(date2,"%lld%lld",&nri,&stepi);      step=stepi; nr=nri;}
       if(type2==DF) {vege=fscanf(date2,"%lld%lf",&nri,&step);       nr=nri;}
       if(type2==FD) {vege=fscanf(date2,"%lf%lld",&nr,&stepi);       step=stepi; }
       if(type2==FF)  vege=fscanf(date2,"%lf%lf",&nr,&step);
       if(type2==DG) {vege=fscanf(date2,"%lld%lg",&nri,&step);       nr=nri;}
       if(type2==DFD){vege=fscanf(date2,"%lld%lf%lld",&nri,&step,&pk);nr=nri;}
       if(type2==DDF){vege=fscanf(date2,"%lld%lld%lf",&pk,&nri,&step);nr=nri;}
  if(type2==DDF)
  	{ y[pk]=step;
  	  x[pk]=nr;
	}
	else{  y[i]=step;
  	       x[i]=nr;
	    }
  //printf("%ld\t%lf\t%lf\n",pk,nr,step); getchar();
  i++;
  }
 while(vege!=EOF);
 fclose(date2);

// printf("%ld\n",i); getchar();
 xmax=0;
 if(type2==DDF) num=pk;
 	else num=i-1;
 x0=10000000;
 for(i=1;i<=num;i++)
 	{if(x[i]>xmax) xmax=x[i];
	 if((x[i]<x0)&&(x[i])) x0=x[i];
	}

  if(!PSZ) PSZ=num/10;

//  printf("i=%d\tnum=%ld\nx0=%ld\txmax=%ld\n",i,num,x0,xmax);
 // printf("a=%lg\n",a);


  a=(log(xmax)-log(x0))/(double)(PSZ);

// printf("a=%lf\n",a);

   //getchar();
  xvon[0]=x0;

  for(i=1;i<=PSZ;i++)
   {xvon[i]=xvon[i-1]*exp(a);
  //  printf("%lf\n",xvon[i]);
   }

for(i=1;i<=PSZ;i++)
 {atlx=xvon[i-1]*exp(a/2.); //printf("interv: %lf -- %lf -- %lf\n\n",xvon[i-1],atlx,xvon[i]);
  for(j=1;j<=num;j++)
      {if( (xvon[i-1]<=x[j] ) && (x[j]<xvon[i]) )
          {//atlx+=x[j];
	   atly+=y[j];
           db++;
	   //printf("j=%ld\ti=%ld\tx[j]=%lf   y[j]=%lf\n",j,i,x[j],y[j]);
	   //printf("\tdb=%lf\t%lf\t%lf\t%lf\n",db,atlx,atly,atly/db);getchar();
          }
      }
  if(db!=0){ ertx[i]=atlx;
             erty[i]=atly/db;
           }
  if(db==0) {ertx[i]=0;
          erty[i]=0;}
  atlx=0;
  atly=0;
  db=0;
 }

 sprintf(logbinnev,"%s_bin",benev);
date2=fopen(logbinnev,"w");
for(i=1;i<=PSZ;i++)
    if(erty[i]!=0)  fprintf(date2,"%lf\t%lg\n",ertx[i],erty[i]);
fclose(date2);

}


