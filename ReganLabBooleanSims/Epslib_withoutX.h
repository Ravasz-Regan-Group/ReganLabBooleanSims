/***************************************************************
 This file is a collection of funtions to draw simple objects
 into an eps file.
 *****************************************************************/
#define GWR 1
#define GYR 2
#define BPR 3
#define BWR 4
#define BYR 5
#define BBR 6
#define GBR 7

#define CARMELO 200
#define CARMELO_GENOTYPE 201

#define MIN_0_MAX_1	1

std::string rgb2hex(int r, int g, int b, bool with_head)
{ std::stringstream ss;
    if (with_head)
        ss << "#";
    ss << std::hex << (r << 16 | g << 8 | b );
    return ss.str();
}

unsigned int COLOR_SCHEME=GWR,VALUE_TRANSFORM_SCHEME=1,EXP_COL_NOW=0;
double CAP_MAX=0;

void EpsSetRgb( FILE *EpsFile, double red, double green, double blue ) {
    fprintf( EpsFile, "%g %g %g srgb\n", red, green, blue); }

void EpsSetRgb_255( FILE *EpsFile, double red, double green, double blue ) {
    fprintf( EpsFile, "%g %g %g srgb\n", red/255., green/255., blue/255.); }

void set_RGB_color_for_value_between_0_and_1(double v,FILE *draw){
    switch (COLOR_SCHEME) {
        case GYR: {if(v<0.5)	EpsSetRgb(draw,2*v,1,0);
        else	        EpsSetRgb(draw,1,1-2*(v-0.5),0);
        }break;
        case GWR: {if(v<0.5)	EpsSetRgb(draw,2*v,1,2*v);
        else		EpsSetRgb(draw,1,1-2*(v-0.5),1-2*(v-0.5));
        }break;
        case BPR: {if(v<0.5)	EpsSetRgb(draw,1-2*v,0,1);
        else		EpsSetRgb(draw,1,0,1-2*(v-0.5));
        }break;
        case BWR: {if(v<0.5)	EpsSetRgb(draw,2*v,2*v,1);
        else		EpsSetRgb(draw,1,1-2*(v-0.5),1-2*(v-0.5));
        }break;
        case BYR: {if(v<0.5)	EpsSetRgb(draw,2*v,2*v,1-2*v); //(0,0,1) -> (1,1,0)
        else		EpsSetRgb(draw,1,1-2*(v-0.5),0); //(1,1,0) -> (1,0,0)
        }break;
        case BBR: {if(v<0.5)	EpsSetRgb(draw,0,0,1-2*v);
        else		EpsSetRgb(draw,2*(v-0.5),0,0);
        }break;
        case GBR: {if(v<0.5)	EpsSetRgb(draw,0,1-2*v,0);
        else		EpsSetRgb(draw,2*(v-0.5),0,0);
        }break;
        default: {if(v<0.5)		EpsSetRgb(draw,1-2*v,2*v,0);
        else		EpsSetRgb(draw,1,1-2*(v-0.5),0);
        }break;
    }
}



/* definition of the point structure */
typedef struct {
    double	x;
    double	y;
} EpsPoint;

typedef struct {
    int red_comp,green_comp,blue_comp;
    char given_name[30];
} RGB_color_ID;

char szinnev[150][30];
unsigned int COL_ORDER_MAX;

void EpsInit( FILE *EpsFile,
             int boundboxlowleftx, int boundboxlowlefty,
             int boundboxuprightx, int boundboxuprighty ) {
    
    fprintf( EpsFile, "%%!PS-Adobe-2.0 EPSF-2.0\n" );
    fprintf( EpsFile, "%%%%BoundingBox: %d %d %d %d\n",
            boundboxlowleftx, boundboxlowlefty,
            boundboxuprightx, boundboxuprighty );
    fprintf( EpsFile, "/cp {closepath} bind def\n" );
    fprintf( EpsFile, "/ef {eofill} bind def\n" );
    fprintf( EpsFile, "/f {fill} bind def\n" );
    fprintf( EpsFile, "/gr {grestore} bind def\n" );
    fprintf( EpsFile, "/gs {gsave} bind def\n" );
    fprintf( EpsFile, "/sa {save} bind def\n" );
    fprintf( EpsFile, "/rs {restore} bind def\n" );
    fprintf( EpsFile, "/l {lineto} bind def\n" );
    fprintf( EpsFile, "/m {moveto} bind def\n" );
    fprintf( EpsFile, "/rm {rmoveto} bind def\n" );
    fprintf( EpsFile, "/n {newpath} bind def\n" );
    fprintf( EpsFile, "/s {stroke} bind def\n" );
    fprintf( EpsFile, "/sh {show} bind def\n" );
    fprintf( EpsFile, "/slc {setlinecap} bind def\n" );
    fprintf( EpsFile, "/slj {setlinejoin} bind def\n" );
    fprintf( EpsFile, "/slw {setlinewidth} bind def\n" );
    fprintf( EpsFile, "/srgb {setrgbcolor} bind def\n" );
    fprintf( EpsFile, "/rot {rotate} bind def\n" );
    fprintf( EpsFile, "/sc {scale} bind def\n" );
    fprintf( EpsFile, "/sd {setdash} bind def\n" );
    fprintf( EpsFile, "/ff {findfont} bind def\n" );
    fprintf( EpsFile, "/sf {setfont} bind def\n" );
    fprintf( EpsFile, "/scf {scalefont} bind def\n" );
    fprintf( EpsFile, "/sw {stringwidth} bind def\n" );
    fprintf( EpsFile, "/tr {translate} bind def\n" );
    fprintf( EpsFile, "/Times-Bold ff 50 scf sf\n" );
}


/* Available eps fonts are: Times-Roman, Times-Italic,
 Times-Bold, Times-Bold-Italic, etc. */

void EpsSetFont( FILE *EpsFile, const char *EpsFontStr, double EpsFontHeight ) {
    fprintf( EpsFile, "/%s ff %g scf sf\n", EpsFontStr, EpsFontHeight ); }

void EpsDrawString( FILE *EpsFile, double deg, double x, double y, const char *TextStr ) {
    fprintf( EpsFile, "newpath\ngsave\n%g %g moveto\n %lf rotate\n (%s) sh s\ngrestore\n",x,y,deg,TextStr);
}

/*
 void EpsDrawString( FILE *EpsFile, double deg, double x, double y, char *TextStr ) {
	
	fprintf( EpsFile, "%.2f rot\n", deg );
	fprintf( EpsFile, "n %g %g m (%s) sh s\n",
 x*cos(deg)+y*sin(deg),
 y*cos(deg)-x*sin(deg),
 TextStr );
 }
 */
void EpsDrawString_Centered( FILE *EpsFile, double x, double y, const char *TextStr ) {
    /* draws text centered horizontally around x and y, pozitioned above the line of y*/
    fprintf( EpsFile, "%g %g moveto\n",x,y);
    fprintf( EpsFile, "(%s) dup stringwidth pop 2 div neg 0 rmoveto show\n",TextStr);
    
}

void EpsDrawString_Centered_Vertical( FILE *EpsFile, double x, double y, const char *TextStr ) {
    /* draws text centered horizontally around x and y, pozitioned above the line of y*/
    fprintf( EpsFile, "%g %g moveto\n",x,y);
    fprintf( EpsFile, "(%s) dup stringwidth pop 2 div neg 0 rmoveto show\n",TextStr);
    
    fprintf( EpsFile, "newpath\ngsave\n%g %g moveto\n",x,y);
    fprintf( EpsFile, "%lf rotate\n (%s) dup stringwidth pop 2 div neg 0 rmoveto show\n",90.,TextStr);
}


void EpsDrawLine( FILE *EpsFile, double x1, double y1, double x2, double y2 ) {
    fprintf( EpsFile, "n %g %g m %g %g l s\n", x1, y1, x2, y2 ); }


void EpsDrawLines( FILE *EpsFile, EpsPoint *P, int npoints ) {
    
    int i;
    
    fprintf( EpsFile, "n %g %g m\n", P[0].x, P[0].y );
    for( i = 1; i < npoints; i++ ) {
        fprintf( EpsFile, "%g %g l\n", P[i].x, P[i].y );
    }
    fprintf( EpsFile, "s\n" );
}


void EpsDrawRectangle( FILE *EpsFile,
                      double lowleftx, double lowlefty,
                      double uprightx, double uprighty ) {
    
    fprintf( EpsFile, "n %g %g m\n", lowleftx, lowlefty );
    fprintf( EpsFile, "%g %g l\n", lowleftx, uprighty );
    fprintf( EpsFile, "%g %g l\n", uprightx, uprighty );
    fprintf( EpsFile, "%g %g l\n", uprightx, lowlefty );
    fprintf( EpsFile, "cp s\n" );
}


void EpsFillRectangle( FILE *EpsFile,
                      double lowleftx, double lowlefty,
                      double uprightx, double uprighty ) {
    
    fprintf( EpsFile, "n %g %g m\n", lowleftx, lowlefty );
    fprintf( EpsFile, "%g %g l\n", lowleftx, uprighty );
    fprintf( EpsFile, "%g %g l\n", uprightx, uprighty );
    fprintf( EpsFile, "%g %g l\n", uprightx, lowlefty );
    fprintf( EpsFile, "cp f s\n" );
}


void EpsFillRectangle_2( FILE *EpsFile,
                        double lowleftx, double lowlefty,
                        double uprightx, double uprighty ) {
    
    fprintf( EpsFile, "n %.2f %.2f m\n", lowleftx, lowlefty );
    fprintf( EpsFile, "%.2f %.2f l\n", lowleftx, uprighty );
    fprintf( EpsFile, "%.2f %.2f l\n", uprightx, uprighty );
    fprintf( EpsFile, "%.2f %.2f l\n", uprightx, lowlefty );
    fprintf( EpsFile, "cp f s\n" );
}



void EpsDrawPoint( FILE *EpsFile, double size, double x, double y ) {
    EpsFillRectangle( EpsFile, x, y, x+size, y+size );
}

void EpsDrawPoint_2( FILE *EpsFile, double size, double x, double y ) {
    EpsFillRectangle_2( EpsFile, x, y, x+size, y+size );
}



void EpsDraw4Poly( FILE *EpsFile,
                  double x1, double y1, double x2, double y2,
                  double x3, double y3, double x4, double y4 ) {
    
    fprintf( EpsFile, "n %g %g m\n", x1, y1 );
    fprintf( EpsFile, "%g %g l\n", x2, y2 );
    fprintf( EpsFile, "%g %g l\n", x3, y3 );
    fprintf( EpsFile, "%g %g l\n", x4, y4 );
    fprintf( EpsFile, "cp s\n" );
}


void EpsFill4Poly( FILE *EpsFile,
                  double x1, double y1, double x2, double y2,
                  double x3, double y3, double x4, double y4 ) {
    
    fprintf( EpsFile, "n %g %g m\n", x1, y1 );
    fprintf( EpsFile, "%g %g l\n", x2, y2 );
    fprintf( EpsFile, "%g %g l\n", x3, y3 );
    fprintf( EpsFile, "%g %g l\n", x4, y4 );
    fprintf( EpsFile, "cp f s\n" );
}

void EpsDrawTriangle( FILE *EpsFile,
                     double x1, double y1, double x2, double y2,double x3, double y3) {
    
    fprintf( EpsFile, "n %g %g m\n", x1, y1 );
    fprintf( EpsFile, "%g %g l\n", x2, y2 );
    fprintf( EpsFile, "%g %g l\n", x3, y3 );
    fprintf( EpsFile, "cp s\n" );
}


void EpsFillTriangle( FILE *EpsFile,
                     double x1, double y1, double x2, double y2,double x3, double y3) {
    
    fprintf( EpsFile, "n %g %g m\n", x1, y1 );
    fprintf( EpsFile, "%g %g l\n", x2, y2 );
    fprintf( EpsFile, "%g %g l\n", x3, y3 );
    fprintf( EpsFile, "cp f s\n" );
}


void EpsDrawPolygon( FILE *EpsFile, EpsPoint *P, int npoints ) {
    
    int i;
    
    fprintf( EpsFile, "n %g %g m\n", P[0].x, P[0].y );
    for( i = 1; i < npoints; i++ ) {
        fprintf( EpsFile, "%g %g l\n", P[i].x, P[i].y );
    }
    fprintf( EpsFile, "cp s\n" );
}


void EpsFillPolygon( FILE *EpsFile, EpsPoint *P, int npoints ) {
    
    int i;
    
    fprintf( EpsFile, "n %g %g m\n", P[0].x, P[0].y );
    for( i = 1; i < npoints; i++ ) {
        fprintf( EpsFile, "%g %g l\n", P[i].x, P[i].y );
    }
    fprintf( EpsFile, "cp f s\n" );
}

void EpsDrawArc( FILE *EpsFile, double centerx, double centery, double rad, double degr1 )
{
    fprintf( EpsFile, "%g %g %g %g %g arc s\n", \
            centery, centerx, centery, rad , degr1);
}

void EpsDrawCircle( FILE *EpsFile, double centerx, double centery, double rad )
{
    fprintf( EpsFile, "n %g %g m %g %g %g 0 360 arc s\n", \
            centerx + rad, centery, centerx, centery, rad );
}

void EpsFillCircle( FILE *EpsFile, double centerx, double centery,
                   double rad ) {
    fprintf( EpsFile, "n %g %g m %g %g %g 0 360 arc f s\n", \
            centerx + rad, centery, centerx, centery, rad );
}

void EpsDrawEllipse( FILE *EpsFile, double centerx, double centery, double a, double b )
{
    fprintf( EpsFile, "gsave newpath \n %g %g translate 1 %g scale \n 0 0 %g 0 360 arc stroke\n stroke grestore\n",
            centerx,centery, b/a, a);
}
void EpsFillEllipse( FILE *EpsFile, double centerx, double centery, double a, double b) {
    fprintf( EpsFile, "gsave newpath \n %g %g translate 1 %g scale \n 0 0 %g 0 360 arc fill stroke\n stroke grestore\n",
            centerx,centery, b/a, a);
}

void EpsFillEllipseArc( FILE *EpsFile, double centerx, double centery, double a, double b,double deg1, double deg2) {
    fprintf( EpsFile, "gsave newpath \n %g %g translate 1 %g scale \n 0 0 %g %g %g arc fill stroke\n stroke grestore\n",
            centerx,centery, b/a, a,deg1,deg2);
}

void EpsDrawEllipseArc( FILE *EpsFile, double centerx, double centery, double a, double b,double deg1, double deg2) {
    fprintf( EpsFile, "gsave newpath \n %g %g translate 1 %g scale \n 0 0 %g %g %g arc stroke\n stroke grestore\n",
            centerx,centery, b/a, a,deg1,deg2);
}


void EpsSetLinewidth( FILE *EpsFile, double linewidth ) {
    fprintf( EpsFile, "%g slw\n", linewidth ); }

void EpsSetLineDash( FILE *EpsFile, const char *dash ) {
    fprintf( EpsFile, "[%s] 0 setdash\n",dash);
}

void setcolors_PAJEK_grow()
{
		  sprintf(szinnev[100], "White");
		  sprintf(szinnev[0], "BrickRed");
		  sprintf(szinnev[1], "Red");
		  sprintf(szinnev[2], "Orange");
		  sprintf(szinnev[3], "Apricot");
		  sprintf(szinnev[4], "Yellow");
    sprintf(szinnev[5], "GreenYellow");
		  sprintf(szinnev[6], "YellowGreen");
		  sprintf(szinnev[7], "ForestGreen");
		  sprintf(szinnev[8], "BlueGreen");
		  sprintf(szinnev[9], "Cyan");
		  sprintf(szinnev[10], "RoyalBlue");
		  sprintf(szinnev[11], "Blue");
		  sprintf(szinnev[12], "Fuchsia");
    sprintf(szinnev[13], "DarkOrchid");
		  sprintf(szinnev[14], "Orchid");
		  sprintf(szinnev[15], "Lavender");
		  sprintf(szinnev[16], "Pink");
		  sprintf(szinnev[17], "Magenta");
		  sprintf(szinnev[18], "RedViolet");
		  sprintf(szinnev[19], "Sepia");
		  sprintf(szinnev[20], "Gray");
		  COL_ORDER_MAX=20;
}

void setcolors_PAJEK(){
    unsigned int i;
    for(i=0;i<150;i++) sprintf(szinnev[i], "Red");
    sprintf(szinnev[100], "White");
		  sprintf(szinnev[0], "Black");
    //                                  sprintf(szinnev[2], "TealBlue");
		  sprintf(szinnev[2], "Blue");
		  sprintf(szinnev[1], "Red");
		  //sprintf(szinnev[2], "Red");
    
    sprintf(szinnev[12], "SeaGreen");
    
		  sprintf(szinnev[16], "Yellow");
		  sprintf(szinnev[8], "MidnightBlue");
		  sprintf(szinnev[17], "LightOrange");
		  sprintf(szinnev[3], "Green");
		  sprintf(szinnev[4], "Cyan");
		  sprintf(szinnev[11], "Magenta");
		  sprintf(szinnev[13], "Pink");
		  sprintf(szinnev[18], "SpringGreen");
		  sprintf(szinnev[27],"GreenYellow");
    sprintf(szinnev[14], "Apricot");
		  sprintf(szinnev[9], "Maroon");
		  sprintf(szinnev[5], "RubineRed");
		  sprintf(szinnev[15], "RedViolet");
		  sprintf(szinnev[20], "Thistle");
		  sprintf(szinnev[10], "Plum");
		  sprintf(szinnev[7], "Violet");
		  sprintf(szinnev[19], "CornflowerBlue");
    sprintf(szinnev[6], "RawSienna");
    sprintf(szinnev[21], "Salmon");
    
    sprintf(szinnev[22], "LimeGreen");
    
		  sprintf(szinnev[53], "Gray");
		  sprintf(szinnev[23], "LightCyan");
		  sprintf(szinnev[25], "LightPurple");
		  sprintf(szinnev[30], "LightGreen");
		  sprintf(szinnev[28], "LSkyBlue");
		  sprintf(szinnev[29], "DarkOrchid");
		  sprintf(szinnev[26], "Melon");
		  sprintf(szinnev[31], "Melon");
		  sprintf(szinnev[32], "LFadedGreen");
		  sprintf(szinnev[33], "WildStrawberry");
    sprintf(szinnev[34], "Fuchsia");
		  sprintf(szinnev[38], "Periwinkle");
		  sprintf(szinnev[36], "JungleGreen");
		  sprintf(szinnev[37], "Cerulean");
		  sprintf(szinnev[50], "SeaGreen");
		  sprintf(szinnev[39], "BrickRed");
		  sprintf(szinnev[40], "Lavender");
		  sprintf(szinnev[41], "Goldenrod");
		  sprintf(szinnev[42], "Sepia");
		  sprintf(szinnev[43], "PineGreen");
    
		  sprintf(szinnev[44], "Tan");
		  sprintf(szinnev[45], "SkyBlue");
		  sprintf(szinnev[46], "OliveGreen");
		  sprintf(szinnev[47], "NavyBlue");
		  sprintf(szinnev[48], "Mulberry");
		  sprintf(szinnev[49], "Turquoise");
    sprintf(szinnev[35], "Tan");
		  sprintf(szinnev[51], "Orange");
		  sprintf(szinnev[52], "LightYellow");
		  sprintf(szinnev[24], "ForestGreen");
    
    //                                  sprintf(szinnev[2], "TealBlue");
    sprintf(szinnev[53], "Blue");
    sprintf(szinnev[54], "Red");
    //sprintf(szinnev[2], "Red");
    
    sprintf(szinnev[55], "SeaGreen");
    
    sprintf(szinnev[56], "Yellow");
    sprintf(szinnev[58], "MidnightBlue");
    sprintf(szinnev[57], "LightOrange");
    sprintf(szinnev[59], "Green");
    sprintf(szinnev[60], "Cyan");
    sprintf(szinnev[61], "Magenta");
    sprintf(szinnev[61], "Pink");
    sprintf(szinnev[63], "SpringGreen");
    sprintf(szinnev[64],"GreenYellow");
    sprintf(szinnev[65], "Apricot");
    sprintf(szinnev[66], "Maroon");
    sprintf(szinnev[67], "RubineRed");
    sprintf(szinnev[68], "RedViolet");
    sprintf(szinnev[69], "Thistle");
    sprintf(szinnev[70], "Plum");
    sprintf(szinnev[71], "Violet");
    sprintf(szinnev[72], "CornflowerBlue");
    sprintf(szinnev[73], "RawSienna");
    sprintf(szinnev[74], "Salmon");
    
    sprintf(szinnev[75], "LimeGreen");
    
    sprintf(szinnev[76], "Gray");
    sprintf(szinnev[77], "LightCyan");
    sprintf(szinnev[78], "LightPurple");
    sprintf(szinnev[79], "LightGreen");
    sprintf(szinnev[80], "LSkyBlue");
    sprintf(szinnev[81], "DarkOrchid");
    sprintf(szinnev[82], "Melon");
    sprintf(szinnev[83], "Melon");
    sprintf(szinnev[84], "LFadedGreen");
    sprintf(szinnev[85], "WildStrawberry");
    sprintf(szinnev[86], "Fuchsia");
    sprintf(szinnev[87], "Periwinkle");
    sprintf(szinnev[88], "JungleGreen");
    sprintf(szinnev[89], "Cerulean");
    sprintf(szinnev[90], "SeaGreen");
    sprintf(szinnev[91], "BrickRed");
    sprintf(szinnev[92], "Lavender");
    sprintf(szinnev[93], "Goldenrod");
    sprintf(szinnev[94], "Sepia");
    sprintf(szinnev[95], "PineGreen");
    
    sprintf(szinnev[96], "Tan");
    sprintf(szinnev[97], "SkyBlue");
    sprintf(szinnev[98], "OliveGreen");
    sprintf(szinnev[99], "NavyBlue");
    sprintf(szinnev[101], "Mahogany");
}



void EpsSetRgb_GREEN_BLACK_RED(FILE *draw,double x){
    double red=0,green=0,blue=0;
    if(x<0) green=1;red=0.5;blue=0.5;
    if(x>1) green=0.5;red=1;blue=0.5;
    if(x<0.5){
        green=1-2*x;
        red=(1-2*x)/2.;
        blue=(1-2*x)/2.;
    }
    if(x==0.50000){red=0;green=0;blue=0;}
    if(x>0.5){
        green=(1-2*(1-x))/2.;
        red=1-2*(1-x);
        blue=(1-2*(1-x))/2.;
    }
    EpsSetRgb(draw,red,green,blue);
}

void EpsSetRgb_GREEN_YELLOW_RED_redheavy(FILE *draw,double x){
    double red,green,blue,v;
    
    if(x<0) green=1;red=0;blue=0;
    if(x>1) green=0;red=1;blue=0;
    
    v=3.*x-1.;
    if(v<0)		EpsSetRgb(draw,1+v,1,0);	 //G-Y-R
    if(v>=0)	EpsSetRgb(draw,1,1-v/2.,0);
}
void EpsSetRgb_GREEN_YELLOW_RED_greenheavy(FILE *draw,double x){
    double red,green,blue,v;
    
    if(x<0) green=1;red=0;blue=0;
    if(x>1) green=0;red=1;blue=0;
    
    v=3.*x-2.;
    if(v<0)		EpsSetRgb(draw,1+v/2.,1,0);	 //G-Y-R
    if(v>=0)	EpsSetRgb(draw,1,1-v,0);
}

void EpsSetRgb_GREEN_YELLOW_RED_balanced(FILE *draw,double x){
    double red,green,blue,v;
    
    if(x<0) green=1;red=0;blue=0;
    if(x>1) green=0;red=1;blue=0;
    
    v=2.*x-1.;
    if(v<0)		EpsSetRgb(draw,1+v,1,0);	 //G-Y-R
    if(v>=0)	EpsSetRgb(draw,1,1-v,0);
}

void EpsSetColor(FILE *EpsFile, char col[50]){
    FILE *colors;
    int vege=2,re = 0,gr =0 ,bl =0 ;
    int nincsmeg, egyenlok,iii;
    char given[50];
    
    colors=fopen("/Users/eravasz/Work/Includes/Eps_Draw/szinek.inf","r");
    nincsmeg=1;
    while((vege!=EOF)&&(nincsmeg))
    {vege=fscanf(colors,"%d%d%d%s",&re,&gr,&bl,given);
        // printf("vege=%d\t%s\t%s\t%d\t%d\t%d\n",vege,col,given,re,gr,bl);//getchar();
        if(vege!=EOF)
        {iii=0;egyenlok=1;
            while((given[iii]!=0)&&egyenlok)
            {if(given[iii]!=col[iii]) egyenlok=0;
                iii++;
            }
            if((egyenlok)&&(col[iii]==0)) nincsmeg=0;
        }
        if(vege==EOF) {printf("%s nincsmeg!",col);//getchar();
        }
    }
    //printf("%s\t%s\t%d\t%d\t%d\n",col,given,re,gr,bl);
    fclose(colors);
    EpsSetRgb(EpsFile,re/255.,gr/255.,bl/255.);
}


void EpsSetColor_PAJEK(FILE *EpsFile, char col[50]){
    FILE *colors;
    int vege=2;
    double re=0.0,gr=0.0,bl = 0.0,inver =0.0;
    int nincsmeg, egyenlok,iii;
    char given[200];
    
    colors=fopen("/Users/eravasz/Dropbox/_CODE/Includes/Eps_Draw/Pajek_szin.inf","r");
    //printf("vege=%d\t%s\t%s\t%d\t%d\t%d\n",vege,col,given,re,gr,bl);   getchar();
    nincsmeg=1;
    while((vege!=EOF)&&(nincsmeg))
    { vege=fscanf(colors,"%s%lf%lf%lf%lf",given,&re,&gr,&bl,&inver);
        //printf("vege=%d\t%s\t%s\t%lf\t%lf\t%lf\n",vege,col,given,re,gr,bl);
        // getchar();
        if(vege!=EOF)
        {iii=0;egyenlok=1;
            while((given[iii]!=0)&&egyenlok)
            {if(given[iii]!=col[iii]) egyenlok=0;
                iii++;
            }
            if((egyenlok)&&(col[iii]==0)) nincsmeg=0;
        }
        if(vege==EOF) {printf("%s Is not available!",col);//getchar();
        }
    }
    //  printf("%s\t%s\t%lf\t%lf\t%lf\n",col,given,re,gr,bl); //getchar();
    fclose(colors);
    
    re=1-re-inver; if(re<0) re=0;
    gr=1-gr-inver;  if(gr<0) gr=0;
    bl=1-bl-inver;  if(bl<0) bl=0;
    //printf("%s\t%s\t%lf\t%lf\t%lf\n",col,given,re,gr,bl); getchar();
    EpsSetRgb(EpsFile,re,gr,bl);
}

void Get_RGB_Color_PAJEK(char col[50],double *r,double *g,double *b){
    FILE *colors;
    int vege=2;
    double re=0.0,gr=0.0,bl = 0.0,inver=0;
    int nincsmeg, egyenlok,iii;
    char given[200];
    
    colors=fopen("/Users/eravasz/Work/Includes/Eps_Draw/Pajek_szin.inf","r");
    //printf("vege=%d\t%s\t%s\t%d\t%d\t%d\n",vege,col,given,re,gr,bl);   getchar();
    nincsmeg=1;
    while((vege!=EOF)&&(nincsmeg))
    { vege=fscanf(colors,"%s%lf%lf%lf%lf",given,&re,&gr,&bl,&inver);
        //printf("vege=%d\t%s\t%s\t%lf\t%lf\t%lf\n",vege,col,given,re,gr,bl);
        // getchar();
        if(vege!=EOF)
        {iii=0;egyenlok=1;
            while((given[iii]!=0)&&egyenlok)
            {if(given[iii]!=col[iii]) egyenlok=0;
                iii++;
            }
            if((egyenlok)&&(col[iii]==0)) nincsmeg=0;
        }
        if(vege==EOF) {printf("%s nincsmeg!",col);getchar();}
    }
    //  printf("%s\t%s\t%lf\t%lf\t%lf\n",col,given,re,gr,bl); //getchar();
    fclose(colors);
    
    re=1-re-inver; if(re<0) re=0;
    gr=1-gr-inver;  if(gr<0) gr=0;
    bl=1-bl-inver;  if(bl<0) bl=0;
    //printf("%s\t%s\t%lf\t%lf\t%lf\n",col,given,re,gr,bl); getchar();
    *r=re;*g=gr;*b=bl;
    //EpsSetRgb(EpsFile,re,gr,bl);
}

/*   int red_comp,green_comp,blue_comp;
 char given_name[30];
 } RGB_color_ID;
 
 
 
 void EpsSetGCColor( FILE *EpsFile, Display* display, GC gc, Colormap colormap,
 char *color_name ) {
 XColor tmpxc, xc;
 
 XLookupColor( display, colormap, color_name, &tmpxc, &xc );
 EpsSetRgb( EpsFile, xc.red/65535.0, xc.green/65535.0, xc.blue/65535.0 ); 
 } 
 */


void EpsEndPage( FILE *EpsFile ) { fprintf( EpsFile, "showpage\n" ); }


void EpsClose( FILE *EpsFile ) { 
    fprintf( EpsFile, "showpage\n" ); 
    fflush( EpsFile );
    fclose( EpsFile );
}


void EpsFlush( FILE *EpsFile ) { fflush( EpsFile ); }


void close_EPSDRAW_file(FILE *draw,char fn[],int conv){
    char fn2[600];
    
    EpsClose(draw);
    
    
    //sprintf(fn2,"echo $PATH\n",fn); system(fn2);
    //getchar();
    
    if(conv==1) {sprintf(fn2,"PATH=${PATH}:/sw/bin\n gm convert -density 300 %s.eps %s.png\n",fn,fn); 
        printf("%s\n",fn2);
        system(fn2);
        printf("\t\t\tdone!\n");
        getchar();
    }
    //	sprintf(fn2,"rm %s.eps\n",fn); system(fn2);
    //getchar();
    sprintf(fn2,"pstopdf %s.eps\n",fn); system(fn2);
    sprintf(fn2,"rm %s.eps\n",fn); system(fn2);
    
    //"convert -verbose -colorspace RGB -resize 800 -interlace none -density 300 -quality 80 $file ${file%%pdf}jpg"
    
}

/******************************************************/
