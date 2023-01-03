
#define E_YES	1
#define E_NO	0
#define CLOSE_E	0.000001

#define CONF_80	1.1826
#define CONF_85	1.4395
#define CONF_90	1.645
#define CONF_95	1.96
#define CONF_99	2.576

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define MAX_LINE 50000

#define MIN_MAX_CENTER	1
#define AVERAGE_CENTER	2
#define DIVIDE_BY_FIRST_GRPOUP_AVE 3

#define EUCLIDES	1
#define CORRELATION	2
#define MI	3
#define MANHATTAN	4
#define CORRELATION_ABS	5
#define ORDERED	0


char GSE_line[MAX_LINE];

//for classes: typedef classname * classname_pt;
// classname * *variable;


typedef double * doublepointer;
typedef doublepointer * triplepointer;
typedef triplepointer * fourpointer;
typedef fourpointer * fivepointer;

typedef int * intpointer;
typedef intpointer * tripleintpointer;

typedef unsigned int * unsignedintpointer;
typedef unsignedintpointer * tripleunsignedintpointer;

typedef unsigned short int * us_shortintpointer;
typedef long int * longintpointer;

typedef unsigned long int * us_longintpointer;

typedef unsigned long long * longnaturalpointer;

typedef longnaturalpointer * double_unsigned_long_int_pointer;
typedef double_unsigned_long_int_pointer * triple_unsigned_long_int_pointer;
 
typedef char * charpointer;
typedef charpointer * stringpointer;
 
typedef FILE * FILEpointer;

class longint_ARRAY{
public:	
	unsigned long int N;
	long int *A;
	
	longint_ARRAY(){
		N=0; A=NULL;
	}
	
	unsigned long int check_for_element(unsigned long int id1){
		unsigned long int i;
		if(N==0) return(0);
		for(i=1;i<=N;i++) if(A[i]==id1) return(i);
		return(0);
	}
	
	void linear_fill(unsigned long int nnow){
		unsigned long int i;
		A=NULL;
		N=nnow;
		A=new long int[N+1];
		for(i=1;i<=N;i++)	A[i]=i;
	}
	unsigned long int add_element(unsigned long int id1){
		 long int *newp;
		unsigned long int i,j;
		
		j=check_for_element(id1);
		if(j==0){
			if(A==NULL){
				N=1;
				A=new long int[2];
				A[1]=id1;
				return(1);
			}
			else{
				newp=A;A=NULL;
				N++;
				A=new long int[N+1];
				for(i=1;i<N;i++)	A[i]=newp[i];
				A[N]=id1;
				//printf("about to delete element list. N=%ld\n",N);getchar();
				delete[] newp;newp=NULL;
				//printf("\t\tdone\n");getchar();
				return(N);
			}
		}
		return(j);
	}
	
	void export_longint_ARRAY(char grname[]){
		FILE *f;
		unsigned long int i;
		
		f=fopen(grname,"w");
		fprintf(f,"%ld\n\n",N);
		for(i=1;i<=N;i++)
			fprintf(f,"%ld\n",A[i]);
		fclose(f);
	}
	void import_longint_ARRAY(char grname[]){
		FILE *f;
		unsigned long int i;
		
		if(N>0){printf("Trying to import into a non-empty array!\n");exit(1);}
		
		f=fopen(grname,"r");
		if(f==NULL) return;
		fscanf(f,"%ld\n\n",&N);
		A=new long int[N+1];
		for(i=1;i<=N;i++)
			fscanf(f,"%ld\n",&(A[i]));
		fclose(f);
	}
	
	void add_longint_ARRAY(longint_ARRAY *to_add){
		unsigned long int i;
		for(i=1;i<=to_add->N;i++)
			add_element(to_add->A[i]);
	}
    void delete_element(unsigned long long int j){
        unsigned long int i;
        
        if(j<N) for(i=j;i<N;i++){ A[i]=A[i+1];}
        N--;
    }
    
    void delete_value(unsigned long int v){
        unsigned long int i,j;
        do {j=check_for_element(v);
            if((j>0)&&(j<=N)){
                for(i=j;i<N;i++){ A[i]=A[i+1];}
                A[N]=0;
                N--;
            }
        }
        while (j>0);
    }
    
	void print(){
        unsigned long int i;
        printf("%ld elements:",N);
		for(i=1;i<=N;i++)
			printf("  A[%ld]=%ld",i,A[i]);
        printf("\n");
    }
	~longint_ARRAY(){
		if(A!=NULL) {delete[] A;A=NULL;}
	}
};


class longlongint_ARRAY_pair{
public:	
	unsigned long long int N;
	char Aname[100],Bname[100];
	long long int *A,*B;
	
	longlongint_ARRAY_pair(){
		N=0; A=NULL;B=NULL;
	}
	
	unsigned long long int check_for_element_A(unsigned long long int id1){
		unsigned long long int i;
		for(i=1;i<=N;i++) if(A[i]==id1) return(i);
		return(0);
	}
	unsigned long long int check_for_element_B(unsigned long long int id1){
		unsigned long long int i;
		for(i=1;i<=N;i++) if(B[i]==id1) return(i);
		return(0);
	}
	unsigned long long int check_for_pair(unsigned long long int id1,unsigned long long int id2){
		unsigned long long int i;
		for(i=1;i<=N;i++) if((A[i]==id1)&&(B[i]==id2)) return(i);
		return(0);
	}
	
	unsigned long long int pair_in_B(unsigned long long int id1){
		unsigned long long int i;
		for(i=1;i<=N;i++) if(A[i]==id1) return(B[i]);
		return(0);
	}
	
	unsigned long long int pair_in_A(unsigned long long int id1){
		unsigned long long int i;
		for(i=1;i<=N;i++) if(B[i]==id1) return(A[i]);
		return(0);
	}
	
	unsigned long long int add_element(unsigned long long int id1,unsigned long long int id2){
		long long int *newpA,*newpB;
		unsigned long int i,j;
		
		j=check_for_pair(id1,id2);
		if(j==0){
			if(A==NULL){
				N=1;
				A=new long long int[2];B=new long long int[2];
				A[1]=id1;B[1]=id2;
				return(1);
			}
			else{
				newpA=A;A=NULL;newpB=B;B=NULL;
				N++;
				A=new long long int[N+1];B=new long long int[N+1];
				for(i=1;i<N;i++)	{
					A[i]=newpA[i];B[i]=newpB[i];
				}
				A[N]=id1;B[N]=id2;
				//printf("about to delete element list. N=%ld\n",N);getchar();
				delete[] newpA;newpA=NULL;
				delete[] newpB;newpB=NULL;
				//printf("\t\tdone\n");getchar();
				return(N);
			}
		}
		return(j);
	}
	
	void delete_element(unsigned long long int j){
		unsigned long long int i;
		
		for(i=j;i<N;i++)	{ A[i]=A[i+1];B[i]=B[i+1];}
		N--;
	}
	void print_all(){
		unsigned int i;
		printf("%lld pairs: ",N);
		for(i=1;i<=N;i++) printf("{%d: %lld, %lld}, ",i,A[i],B[i]);
		printf("\n");
	}
	void add_longlongint_ARRAY_pair(longlongint_ARRAY_pair *to_add){
		unsigned long int i;
		for(i=1;i<=to_add->N;i++)
			add_element(to_add->A[i],to_add->B[i]);
	}
	
	~longlongint_ARRAY_pair(){
		if(N>0){ delete[] A;A=NULL;		delete[] B;B=NULL;}
	}
};

class longint_ARRAY_pair{
public:	
	unsigned  long int N;
	char Aname[100],Bname[100];
	long int *A,*B;
	
	longint_ARRAY_pair(){
		N=0; A=NULL;B=NULL;
        sprintf(Aname,"");sprintf(Bname,"");
	}
    
    longint_ARRAY_pair(longint_ARRAY_pair *a){
        N=a->N;
        A=new long int[N+1];
        B=new long int[N+1];
        for(unsigned  long int i=1;i<=N;i++){
            A[i]=a->A[i];B[i]=a->B[i];
        }
        strcpy(Aname,a->Aname);
        strcpy(Bname,a->Bname);
    }
    
	unsigned long  int check_for_element_A(unsigned  long int id1){
		unsigned  long int i;
		for(i=1;i<=N;i++) if(A[i]==id1) return(i);
		return(0);
	}
	unsigned long  int check_for_element_B(unsigned  long int id1){
		unsigned  long int i;
		for(i=1;i<=N;i++) if(B[i]==id1) return(i);
		return(0);
	}
	unsigned long  int check_for_pair(unsigned  long int id1,unsigned  long int id2){
		unsigned  long int i;
		for(i=1;i<=N;i++) if((A[i]==id1)&&(B[i]==id2)) return(i);
		return(0);
	}
	
	unsigned long  int pair_in_B(unsigned  long int id1){
		unsigned  long int i;
		for(i=1;i<=N;i++) if(A[i]==id1) return(B[i]);
		return(0);
	}
	
	unsigned  long int pair_in_A(unsigned  long int id1){
		unsigned  long int i;
		for(i=1;i<=N;i++) if(B[i]==id1) return(A[i]);
		return(0);
	}
	
	unsigned  long int add_element(unsigned  long int id1,unsigned  long int id2){
		long  int *newpA,*newpB;
		unsigned long int i,j;
		
		j=check_for_pair(id1,id2);
		if(j==0){
			if(A==NULL){
				N=1;
				A=new long  int[2];B=new  long int[2];
				A[1]=id1;B[1]=id2;
				return(1);
			}
			else{
				newpA=A;A=NULL;newpB=B;B=NULL;
				N++;
				A=new long  int[N+1];B=new  long int[N+1];
				for(i=1;i<N;i++)	{
					A[i]=newpA[i];B[i]=newpB[i];
				}
				A[N]=id1;B[N]=id2;
				//printf("about to delete element list. N=%ld\n",N);getchar();
				delete[] newpA;newpA=NULL;
				delete[] newpB;newpB=NULL;
				//printf("\t\tdone\n");getchar();
				return(N);
			}
		}
		return(j);
	}
	
	void delete_element(unsigned int j){
		int i;
		
		for(i=j;i<N;i++)	{ A[i]=A[i+1];B[i]=B[i+1];}
		N--;
	}
	void print_all(){
		unsigned int i;
		printf("%ld pairs: ",N);
		for(i=1;i<=N;i++) printf("{%d: %ld, %ld}, ",i,A[i],B[i]);
		printf("\n");
	}
	void add_longint_ARRAY_pair(longint_ARRAY_pair *to_add){
		unsigned long int i;
		for(i=1;i<=to_add->N;i++)
			add_element(to_add->A[i],to_add->B[i]);
	}
	/*
	~longint_ARRAY_pair(){
		if(N>0){
            N=0;
            delete[] A;A=NULL;
            delete[] B;B=NULL;
        }
	}
    */
    longint_ARRAY *extract_longintarray_from_A(){
        longint_ARRAY *AA;
        AA=new longint_ARRAY();
        for(long int i=1;i<=N;i++)
            AA->add_element(A[i]);
        return(AA);
    }
};

typedef longint_ARRAY_pair * longint_ARRAY_pair_array;


class longlongint_ARRAY{
public:	
	unsigned long long int N;
	long long int *A;
	
	longlongint_ARRAY(){
		N=0; A=NULL;
	}
	
	unsigned long long int check_for_element(unsigned long int id1){
		unsigned long long int i;
		for(i=1;i<=N;i++) if(A[i]==id1) return(i);
		return(0);
	}
	
	unsigned long long int add_element(unsigned long int id1){
		long long int *newp;
		unsigned long int i,j;
		
		j=check_for_element(id1);
		if(j==0){
			if(A==NULL){
				N=1;
				A=new long long int[2];
				A[1]=id1;
				return(1);
			}
			else{
				newp=A;A=NULL;
				N++;
				A=new long long int[N+1];
				for(i=1;i<N;i++)	A[i]=newp[i];
				A[N]=id1;
				//printf("about to delete element list. N=%ld\n",N);getchar();
				delete[] newp;newp=NULL;
				//printf("\t\tdone\n");getchar();
				return(N);
			}
		}
		return(j);
	}
	
	void export_longint_ARRAY(char grname[]){
		FILE *f;
		unsigned long long int i;
		
		f=fopen(grname,"w");
		fprintf(f,"%lld\n\n",N);
		for(i=1;i<=N;i++)
			fprintf(f,"%lld\n",A[i]);
		fclose(f);
	}
	void add_longlongint_ARRAY(longlongint_ARRAY *to_add){
		unsigned long int i;
		for(i=1;i<=to_add->N;i++)
			add_element(to_add->A[i]);
	}
	
	~longlongint_ARRAY(){
		delete[] A;A=NULL;
	}
};

void add_list_to_longint_ARRAY(longint_ARRAY *final,longint_ARRAY *to_add){
	unsigned long int i;
	for(i=1;i<=to_add->N;i++)
		final->add_element(to_add->A[i]);
}

typedef longint_ARRAY * longint_ARRAY_list;
typedef longint_ARRAY_list * longint_ARRAY_matrix;
typedef longlongint_ARRAY * longlongint_ARRAY_list;

longint_ARRAY *common_elements_of_longint_ARRAYs(longint_ARRAY *peak1,longint_ARRAY *peak2){
	longint_ARRAY *gl;
	unsigned long int i;
	
	gl=new longint_ARRAY();
	for(i=1;i<=peak1->N;i++){
		if(peak2->check_for_element(peak1->A[i])) gl->add_element(peak1->A[i]);
	}
	return(gl);
}


void VENN_DIAGRAM(longint_ARRAY *peak1,longint_ARRAY *peak2,longint_ARRAY *peak1_only,longint_ARRAY *peak2_only,longint_ARRAY *both_peaks){
	unsigned int i;
	for(i=1;i<=peak1->N;i++){
		if(peak2->check_for_element(peak1->A[i])) both_peaks->add_element(peak1->A[i]);
		else peak1_only->add_element(peak1->A[i]);
	}
	for(i=1;i<=peak2->N;i++){
		if(peak1->check_for_element(peak2->A[i])) both_peaks->add_element(peak2->A[i]);
		else peak2_only->add_element(peak2->A[i]);
	}
}
unsigned long int Nr_of_common_elements(longint_ARRAY *peak1,longint_ARRAY *peak2){
	unsigned int i;
	unsigned long int n;
	n=0;
	for(i=1;i<=peak1->N;i++){
		if(peak2->check_for_element(peak1->A[i])) n++;
	}
	return(n);
}

unsigned long int Nr_of_union(longint_ARRAY *peak1,longint_ARRAY *peak2){
	unsigned int i;
	unsigned long int n;
	longint_ARRAY *a;
	a=new longint_ARRAY();
	for(i=1;i<=peak1->N;i++) a->add_element(peak1->A[i]);
	for(i=1;i<=peak2->N;i++) a->add_element(peak2->A[i]);
	n=a->N;
	delete a;a=NULL;
	return(n);
}


void VENN_DIAGRAM_in_three(longint_ARRAY *peak1,longint_ARRAY *peak2,longint_ARRAY *peak3,
							longint_ARRAY *peak1_only,longint_ARRAY *peak2_only,longint_ARRAY *peak3_only,
							longint_ARRAY *peak1_peak2_only,longint_ARRAY *peak1_peak3_only,longint_ARRAY *peak2_peak3_only,
							longint_ARRAY *all_three){
	unsigned int i;
	
	for(i=1;i<=peak1->N;i++){
		if(peak2->check_for_element(peak1->A[i])) {
			if(peak3->check_for_element(peak1->A[i])) all_three->add_element(peak1->A[i]);
			else peak1_peak2_only->add_element(peak1->A[i]);
		}
		else{
			if(peak3->check_for_element(peak1->A[i])) peak1_peak3_only->add_element(peak1->A[i]);
			else peak1_only->add_element(peak1->A[i]);
		}
	}
	for(i=1;i<=peak2->N;i++){
		if(peak1->check_for_element(peak2->A[i])) {
			if(peak3->check_for_element(peak2->A[i])) all_three->add_element(peak2->A[i]);
			else peak1_peak2_only->add_element(peak2->A[i]);
		}
		else{
			if(peak3->check_for_element(peak2->A[i])) peak2_peak3_only->add_element(peak2->A[i]);
			else peak2_only->add_element(peak2->A[i]);
		}
	}
	for(i=1;i<=peak3->N;i++){
		if(peak1->check_for_element(peak3->A[i])) {
			if(peak2->check_for_element(peak3->A[i])) all_three->add_element(peak3->A[i]);
			else peak1_peak3_only->add_element(peak3->A[i]);
		}
		else{
			if(peak2->check_for_element(peak3->A[i])) peak2_peak3_only->add_element(peak3->A[i]);
			else peak3_only->add_element(peak3->A[i]);
		}
	}
}

void VENN_DIAGRAM_in_three_firstlist(longint_ARRAY_pair_array *peaks){
    unsigned int i;
	
	for(i=1;i<=peaks[0]->N;i++){
		if(peaks[1]->check_for_element_A(peaks[0]->A[i])) {
			if(peaks[2]->check_for_element_A(peaks[0]->A[i])) peaks[9]->add_element(peaks[0]->A[i],peaks[0]->B[i]);
			else peaks[6]->add_element(peaks[0]->A[i],peaks[0]->B[i]);
		}
		else{
			if(peaks[2]->check_for_element_A(peaks[0]->A[i])) peaks[7]->add_element(peaks[0]->A[i],peaks[0]->B[i]);
			else peaks[3]->add_element(peaks[0]->A[i],peaks[0]->B[i]);
		}
	}
	for(i=1;i<=peaks[1]->N;i++){
		if(peaks[0]->check_for_element_A(peaks[1]->A[i])) {
			if(peaks[2]->check_for_element_A(peaks[1]->A[i])) peaks[9]->add_element(peaks[1]->A[i],peaks[1]->B[i]);
			else peaks[6]->add_element(peaks[1]->A[i],peaks[1]->B[i]);
		}
		else{
			if(peaks[2]->check_for_element_A(peaks[1]->A[i])) peaks[8]->add_element(peaks[1]->A[i],peaks[1]->B[i]);
			else peaks[4]->add_element(peaks[1]->A[i],peaks[1]->B[i]);
		}
	}
	for(i=1;i<=peaks[2]->N;i++){
		if(peaks[0]->check_for_element_A(peaks[2]->A[i])) {
			if(peaks[1]->check_for_element_A(peaks[2]->A[i])) peaks[9]->add_element(peaks[2]->A[i],peaks[2]->B[i]);
			else peaks[7]->add_element(peaks[2]->A[i],peaks[2]->B[i]);
		}
		else{
			if(peaks[1]->check_for_element_A(peaks[2]->A[i])) peaks[8]->add_element(peaks[2]->A[i],peaks[2]->B[i]);
			else peaks[5]->add_element(peaks[2]->A[i],peaks[2]->B[i]);
		}
	}
}

double get_LCB(double av,double sz,int CONF_inter){
	double lcb,xc = 0.0;
	
	if(CONF_inter==90) xc=CONF_90*sz;
	else {
		if(CONF_inter==95) xc=CONF_95*sz;
		else {
			if(CONF_inter==90) xc=CONF_99*sz;
			else{ printf("Need contant for CONF %d\n",CONF_inter);getchar();}
		}
	}
	
	if(!((av-xc<=0)&&(av+xc>=0))){
		if(av<0) lcb=av+xc;
		else lcb=av-xc;
	}
	else lcb=0;
	return(lcb);
}

struct poz2d{
		double x,y;
};

typedef poz2d * poz2d_list;

struct line{
		double a,b,c,m;
		unsigned short int vertical;	
};

double minimum(unsigned long int N, double *x){
	unsigned long int i;
	double m;
	m=x[1];
	for(i=2;i<=N;i++) 
		if(x[i]<m) m=x[i];
	return(m);
}
double maximum(unsigned long int N, double *x){
	unsigned long int i;
	double m;
	m=x[1];
	for(i=2;i<=N;i++) 
		if(x[i]>m) m=x[i];
	return(m);
}

unsigned int maximum(unsigned long int N, unsigned int  *x){
	unsigned long int i;
	double m;
	m=x[1];
	for(i=2;i<=N;i++)
		if(x[i]>m) m=x[i];
	return(m);
}

double maximum(double a,double b){
	if(a>b) return(a);
	else return(b);
}
/*
double factorial(unsigned short int n){
	double f=1;
	int i;
	if(n<=1) return(1);
	for(i=2;i<=n;i++) f*=i;
	return(f);
}
*/

double Euclidean_distance(unsigned long int N, double x[], double y[]){
	unsigned long int i,ok;
	double dist=0;
	
	dist=0;ok=0;
	for(i=1;i<=N;i++)
		dist+=(x[i]-y[i])*(x[i]-y[i]);
	return(sqrt(dist));
}

double correlation(unsigned long int N, double *x, double *y){
	unsigned long int i;
	double av_egy=0,av_ket=0,av_egyket=0,s_egy=0,s_ket=0;
	
	if(N<=1) return(0);
	av_egy=0;av_ket=0;av_egyket=0;s_egy=0;s_ket=0;
	for(i=1;i<=N;i++){
			av_egy+=x[i];
			av_ket+=y[i];
			av_egyket+=x[i]*y[i];
			s_egy+=x[i]*x[i];
			s_ket+=y[i]*y[i];
		}
	av_egy/=(double)N;		
	av_ket/=(double)N;		
	av_egyket/=(double)N;		
	s_egy/=(double)N;		
	s_ket/=(double)N;	
	if((s_egy>0.00000000001)&&(s_ket>0.00000000001)){
		s_egy=sqrt(s_egy-av_egy*av_egy);
		s_ket=sqrt(s_ket-av_ket*av_ket);
		return((av_egyket-av_egy*av_ket)/(s_egy*s_ket));
	}
	else return(0.);
}


double factorial(int n) {
	//Returns the value n! as a double number.  
	
	static int ntop=4; 
	static float a[33]={1.0,1.0,2.0,6.0,24.0}; //Fill in table only as required. 
	int j; 
	if (n < 0) {printf("Negative factorial in routine factortial\n");getchar();} 
	if (n > 32) return exp(gammln(n+1.0)); //Larger value than size of table is required.
				//Actually, this big a value is going to overflow on many computers, 
				//but no harm in trying. 
	while (ntop<n) { ///Fill in table up to desired value. 
		j=ntop++; 
		a[ntop]=a[j]*ntop; 
	} 
	return a[n]; 
}

double Gamma_half_int(unsigned short int n){ // input TWICE the half integer needed!!!
	double f=2./M_2_SQRTPI;
	int i;
	if(n<1) return(1);
	if(n==1) return(f);
	if(n/2==n/2.) return(factorial(n/2-1));
	for(i=1;i<n;i+=2) f*=i/2.;
	return(f);
}

void hpsort(int n, double ra[])
// sorts an array ra[1,n] into ascending order (small to large)
{
	int i,ir,j,l;
	double rra;

	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
		} else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra;
	}
}

void hpsort(unsigned long long int n, double ra[])
// sorts an array ra[1,n] into ascending order (small to large)
{
	unsigned long long i,ir,j,l;
	double rra;
    
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
		} else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra;
	}
}

void hpsort_coshuffle_indices(unsigned long int n, double ra[],unsigned long int ind_ra[])
// sorts an array ra[1,n] into ascending order (small to large)
// original indices of the ra members will end up in ind_ra (it does not matter what it contains)
{
	unsigned long int i,ir,j,l,ind_rra;
	double rra;
	
	for(i=1;i<=n;i++) ind_ra[i]=i; 
	
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	//printf("l=%d\tir=%d\n",l,ir);getchar();
	for (;;) {
		if (l > 1) {
			rra=ra[--l];	ind_rra=ind_ra[l];
			//printf("rra=%lf from %d\n",rra,l);getchar();
		} else {
			rra=ra[ir];		ind_rra=ind_ra[ir];
			ra[ir]=ra[1];	ind_ra[ir]=ind_ra[1];
			//printf("rra=%lf\tra[%d] becomes ra[1]=%lf\n",rra,ir,ra[1]);getchar();
			if (--ir == 1) {
				ra[1]=rra;  ind_ra[1]=ind_rra;
				//printf("\tra[1] becomes rra=%lf\n",rra);getchar();
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			//printf("i=%d\tj=%d\n",i,j);getchar();
			if (j < ir && ra[j] < ra[j+1]) j++;
			//printf("\tj=%d\n",j);getchar();
			if (rra < ra[j]) {
				ra[i]=ra[j];   ind_ra[i]=ind_ra[j];
				//printf("\t\tra[%d] becomes ra[j]=%lf\n",i,j,ra[j]);getchar();
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra; ind_ra[i]=ind_rra;
		//printf("\t\tra[%d] becomes rra=%lf\n",i,rra);getchar();
	}
}

void hpsort_coshuffle_indices(int n, unsigned long long ra[],unsigned long long ind_ra[])
// sorts an array ra[1,n] into ascending order (small to large)
// original indices of the ra members will end up in ind_ra (it does not matter what it contains)
{
	unsigned long long i,ir,j,l,ind_rra;
	unsigned long long rra;
	
	for(i=1;i<=n;i++) ind_ra[i]=i; 
	
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	//printf("l=%d\tir=%d\n",l,ir);getchar();
	for (;;) {
		if (l > 1) {
			rra=ra[--l];	ind_rra=ind_ra[l];
			//printf("rra=%lf from %d\n",rra,l);getchar();
		} else {
			rra=ra[ir];		ind_rra=ind_ra[ir];
			ra[ir]=ra[1];	ind_ra[ir]=ind_ra[1];
			//printf("rra=%lf\tra[%d] becomes ra[1]=%lf\n",rra,ir,ra[1]);getchar();
			if (--ir == 1) {
				ra[1]=rra;  ind_ra[1]=ind_rra;
				//printf("\tra[1] becomes rra=%lf\n",rra);getchar();
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			//printf("i=%d\tj=%d\n",i,j);getchar();
			if (j < ir && ra[j] < ra[j+1]) j++;
			//printf("\tj=%d\n",j);getchar();
			if (rra < ra[j]) {
				ra[i]=ra[j];   ind_ra[i]=ind_ra[j];
				//printf("\t\tra[%d] becomes ra[j]=%lf\n",i,j,ra[j]);getchar();
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra; ind_ra[i]=ind_rra;
		//printf("\t\tra[%d] becomes rra=%lf\n",i,rra);getchar();
	}
}

void hpsort_coshuffle_indices(unsigned int n, unsigned int ra[],unsigned int ind_ra[])
// sorts an array ra[1,n] into ascending order (small to large)
// original indices of the ra members will end up in ind_ra (it does not matter what it contains)
{
	unsigned int i,ir,j,l,ind_rra;
	unsigned int rra;
	
	for(i=1;i<=n;i++) ind_ra[i]=i; 
	
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	//printf("l=%d\tir=%d\n",l,ir);getchar();
	for (;;) {
		if (l > 1) {
			rra=ra[--l];	ind_rra=ind_ra[l];
			//printf("rra=%lf from %d\n",rra,l);getchar();
		} else {
			rra=ra[ir];		ind_rra=ind_ra[ir];
			ra[ir]=ra[1];	ind_ra[ir]=ind_ra[1];
			//printf("rra=%lf\tra[%d] becomes ra[1]=%lf\n",rra,ir,ra[1]);getchar();
			if (--ir == 1) {
				ra[1]=rra;  ind_ra[1]=ind_rra;
				//printf("\tra[1] becomes rra=%lf\n",rra);getchar();
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			//printf("i=%d\tj=%d\n",i,j);getchar();
			if (j < ir && ra[j] < ra[j+1]) j++;
			//printf("\tj=%d\n",j);getchar();
			if (rra < ra[j]) {
				ra[i]=ra[j];   ind_ra[i]=ind_ra[j];
				//printf("\t\tra[%d] becomes ra[j]=%lf\n",i,j,ra[j]);getchar();
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra; ind_ra[i]=ind_rra;
		//printf("\t\tra[%d] becomes rra=%lf\n",i,rra);getchar();
	}
}

void hpsort_coshuffle_indices(unsigned long long n, unsigned long long ra[],unsigned long long ind_ra[])
// sorts an array ra[1,n] into ascending order (small to large)
// original indices of the ra members will end up in ind_ra (it does not matter what it contains)
{
	unsigned long long i,ir,j,l,ind_rra;
	unsigned long long rra;
	
	for(i=1;i<=n;i++) ind_ra[i]=i; 
	
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	//printf("l=%d\tir=%d\n",l,ir);getchar();
	for (;;) {
		if (l > 1) {
			rra=ra[--l];	ind_rra=ind_ra[l];
			//printf("rra=%lf from %d\n",rra,l);getchar();
		} else {
			rra=ra[ir];		ind_rra=ind_ra[ir];
			ra[ir]=ra[1];	ind_ra[ir]=ind_ra[1];
			//printf("rra=%lf\tra[%d] becomes ra[1]=%lf\n",rra,ir,ra[1]);getchar();
			if (--ir == 1) {
				ra[1]=rra;  ind_ra[1]=ind_rra;
				//printf("\tra[1] becomes rra=%lf\n",rra);getchar();
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			//printf("i=%d\tj=%d\n",i,j);getchar();
			if (j < ir && ra[j] < ra[j+1]) j++;
			//printf("\tj=%d\n",j);getchar();
			if (rra < ra[j]) {
				ra[i]=ra[j];   ind_ra[i]=ind_ra[j];
				//printf("\t\tra[%d] becomes ra[j]=%lf\n",i,j,ra[j]);getchar();
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra; ind_ra[i]=ind_rra;
		//printf("\t\tra[%d] becomes rra=%lf\n",i,rra);getchar();
	}
}

void hpsort_coshuffle_indices(unsigned long long n, double ra[],unsigned long long ind_ra[])
// sorts an array ra[1,n] into ascending order (small to large)
// original indices of the ra members will end up in ind_ra (it does not matter what it contains)
{
	unsigned long long i,ir,j,l,ind_rra;
	double rra;
	
	for(i=1;i<=n;i++) ind_ra[i]=i; 
	
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	//printf("l=%d\tir=%d\n",l,ir);getchar();
	for (;;) {
		if (l > 1) {
			rra=ra[--l];	ind_rra=ind_ra[l];
			//printf("rra=%lf from %d\n",rra,l);getchar();
		} else {
			rra=ra[ir];		ind_rra=ind_ra[ir];
			ra[ir]=ra[1];	ind_ra[ir]=ind_ra[1];
			//printf("rra=%lf\tra[%d] becomes ra[1]=%lf\n",rra,ir,ra[1]);getchar();
			if (--ir == 1) {
				ra[1]=rra;  ind_ra[1]=ind_rra;
				//printf("\tra[1] becomes rra=%lf\n",rra);getchar();
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			//printf("i=%d\tj=%d\n",i,j);getchar();
			if (j < ir && ra[j] < ra[j+1]) j++;
			//printf("\tj=%d\n",j);getchar();
			if (rra < ra[j]) {
				ra[i]=ra[j];   ind_ra[i]=ind_ra[j];
				//printf("\t\tra[%d] becomes ra[j]=%lf\n",i,j,ra[j]);getchar();
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra; ind_ra[i]=ind_rra;
		//printf("\t\tra[%d] becomes rra=%lf\n",i,rra);getchar();
	}
}

void hpsort_coshuffle(unsigned long long n, unsigned long long ra[],unsigned long long ind_ra[])
// sorts an array ra[1,n] into ascending order (small to large)
// ind_ra[] will be shuffled together with ra

{
	unsigned long long i,ir,j,l,ind_rra;
	unsigned long long rra;
	
	
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	//printf("l=%d\tir=%d\n",l,ir);getchar();
	for (;;) {
		if (l > 1) {
			rra=ra[--l];	ind_rra=ind_ra[l];
			//printf("rra=%lf from %d\n",rra,l);getchar();
		} else {
			rra=ra[ir];		ind_rra=ind_ra[ir];
			ra[ir]=ra[1];	ind_ra[ir]=ind_ra[1];
			//printf("rra=%lf\tra[%d] becomes ra[1]=%lf\n",rra,ir,ra[1]);getchar();
			if (--ir == 1) {
				ra[1]=rra;  ind_ra[1]=ind_rra;
				//printf("\tra[1] becomes rra=%lf\n",rra);getchar();
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			//printf("i=%d\tj=%d\n",i,j);getchar();
			if (j < ir && ra[j] < ra[j+1]) j++;
			//printf("\tj=%d\n",j);getchar();
			if (rra < ra[j]) {
				ra[i]=ra[j];   ind_ra[i]=ind_ra[j];
				//printf("\t\tra[%d] becomes ra[j]=%lf\n",i,j,ra[j]);getchar();
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra; ind_ra[i]=ind_rra;
		//printf("\t\tra[%d] becomes rra=%lf\n",i,rra);getchar();
	}
}

void hpsort_coshuffle(long long unsigned n, double ra[], long long unsigned ind_ra[])
// sorts an array ra[1,n] into ascending order (small to large)
// ind_ra[] will be shuffled together with ra

{
	unsigned long long i,ir,j,l,ind_rra;
	double rra;
	
	
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	//printf("l=%d\tir=%d\n",l,ir);getchar();
	for (;;) {
		if (l > 1) {
			rra=ra[--l];	ind_rra=ind_ra[l];
			//printf("rra=%lf from %d\n",rra,l);getchar();
		} else {
			rra=ra[ir];		ind_rra=ind_ra[ir];
			ra[ir]=ra[1];	ind_ra[ir]=ind_ra[1];
			//printf("rra=%lf\tra[%d] becomes ra[1]=%lf\n",rra,ir,ra[1]);getchar();
			if (--ir == 1) {
				ra[1]=rra;  ind_ra[1]=ind_rra;
				//printf("\tra[1] becomes rra=%lf\n",rra);getchar();
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			//printf("i=%d\tj=%d\n",i,j);getchar();
			if (j < ir && ra[j] < ra[j+1]) j++;
			//printf("\tj=%d\n",j);getchar();
			if (rra < ra[j]) {
				ra[i]=ra[j];   ind_ra[i]=ind_ra[j];
				//printf("\t\tra[%d] becomes ra[j]=%lf\n",i,j,ra[j]);getchar();
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra; ind_ra[i]=ind_rra;
		//printf("\t\tra[%d] becomes rra=%lf\n",i,rra);getchar();
	}
}

void hpsort(int n, unsigned long long ra[])
// sorts an array ra[1,n] into ascending order
{
	int i,ir,j,l;
	unsigned long long rra;

	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
		} else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra;
	}
}

void hpsort(longint_ARRAY *ra_in)
// sorts an array ra[1,n] into ascending order
{
	unsigned long int i,ir,j,l;
	unsigned long int rra,n;
	
	n=ra_in->N;
	
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra_in->A[--l];
		} else {
			rra=ra_in->A[ir];
			ra_in->A[ir]=ra_in->A[1];
			if (--ir == 1) {
				ra_in->A[1]=rra;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra_in->A[j] < ra_in->A[j+1]) j++;
			if (rra < ra_in->A[j]) {
				ra_in->A[i]=ra_in->A[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra_in->A[i]=rra;
	}
}

int length(char *s)
{ int hatso=0;
	while(s[hatso])
		hatso++;
	return(hatso);
}

void boxin_coordinates(unsigned int N, poz2d *coord,double x1,double y1,double x2,double y2,short int KEEP_ASPECT){
	unsigned int i;
	int min_x,max_x,min_y,max_y,sc_x,sc_y;
	
	if((N==0)||(x2-x1<=0)||(y2-y1<=0)) {
		printf("Bad boxin_coordinates call\n");
		return;
	}
	min_x=max_x=coord[1].x;
	min_y=max_y=coord[1].y;
	for(i=2;i<=N;i++){
			if(coord[i].x<min_x) min_x=coord[i].x;
			if(coord[i].x>max_x) max_x=coord[i].x;
			if(coord[i].y<min_y) min_y=coord[i].y;
			if(coord[i].y>max_y) max_y=coord[i].y;
	}
	sc_x=max_x-min_x;
	sc_y=max_y-min_y;
	
	for(i=1;i<=N;i++){
		coord[i].x-=min_x;
		coord[i].y-=min_y;
	}
	
	if(!KEEP_ASPECT){
		if(sc_x!=0) for(i=1;i<=N;i++)
			coord[i].x=coord[i].x/sc_x;
		if(sc_y!=0) for(i=1;i<=N;i++)	
			coord[i].y=coord[i].y/sc_y;
		
		for(i=1;i<=N;i++){
			coord[i].x=coord[i].x*(x2-x1)+x1;
			coord[i].y=coord[i].y*(y2-y1)+y1;
		}
		return;
	}
	else{
		if((sc_y)&&((sc_y/sc_x)<(y2-y1)/(x2-x1))){
			for(i=1;i<=N;i++){ //scale by x axis down to 1 the up to new window
				coord[i].x=coord[i].x/sc_x;
				coord[i].y=coord[i].y/sc_x;
			}
			for(i=1;i<=N;i++){
				coord[i].x=coord[i].x*(x2-x1)+x1;
				coord[i].y=coord[i].y*(x2-x1)+y1;
				coord[i].y+=((y2-y1)-(x2-x1))/2.;
			}
		}
		if((sc_x)&&((sc_y/sc_x)>(y2-y1)/(x2-x1))){
			for(i=1;i<=N;i++){ //scale by y axis down to 1 the up to new window
				coord[i].x=coord[i].x/sc_y;
				coord[i].y=coord[i].y/sc_y;
			}
			for(i=1;i<=N;i++){
				coord[i].x=coord[i].x*(y2-y1)+x1;
				coord[i].y=coord[i].y*(y2-y1)+y1;
				coord[i].y+=((x2-x1)-(y2-y1))/2.;
			}
		}
		
	}	
}

unsigned long int separate_line_smallchunks(char line[],char a[][500],char bc,int MAX){
	unsigned long int i,j,n,start,k;
	n=0;
	char *bp;
	char break_char[4];
	
	bp=new char[strlen(line)+1];
	for(i=0;i<=strlen(line);i++) bp[i]=0;
	
	for(i=0;i<=3;i++) break_char[i]=0;
	break_char[0]=bc;
	
	i=0;
	while(break_char[i]!=0){
		for(j=0;j<strlen(line);j++) if(line[j]==break_char[i]) {bp[j]=1;}
		i++;
	}
	
	i=0;
	while(break_char[i]==0){
		i++;
	}
	
	start=i; j=1;
	while((line[i]!=0)&&(i<strlen(line))){
		//printf("start of %d th is %d (i=%d, line[i]=%c %d)\n",j,start,i,line[i],line[i]);
		while((bp[i]==0)&&(i<strlen(line))&&(i-start<5000)){
			a[j][i-start]=line[i];
			i++;
		}
		a[j][i-start]=0;
		if(i-start>=5000) {printf("separate_line error: piece %ld does not fit in length 5000\n%s\n",j,line);}
		for(k=i-start+1;k<5000;k++) a[j][k]=0;
		//printf("Break at i=%d\n",i);//getchar();
		j++; 
		if(j>=MAX-1) {printf("separate_line_smallchunks: Incoming array too small for string fragments!\n%s\n",line);exit(1);}
		i++;
		start=i;
	}
	//printf("returning n=%d\n",j-1);getchar();
	delete[] bp;bp=NULL;
	return(j-1);
}


unsigned long int separate_line_tiny(char line[],char a[][50],char bc,int MAX){
    unsigned long int i,j,n,start,k;
    n=0;
    char *bp;
    char break_char[4];
    
    bp=new char[strlen(line)+1];
    for(i=0;i<=strlen(line);i++) bp[i]=0;
    
    for(i=0;i<=3;i++) break_char[i]=0;
    break_char[0]=bc;
    
    i=0;
    while(break_char[i]!=0){
        for(j=0;j<strlen(line);j++) if(line[j]==break_char[i]) {bp[j]=1;}
        i++;
    }
    
    i=0;
    while(break_char[i]==0){
        i++;
    }
    
    start=i; j=1;
    while((line[i]!=0)&&(i<strlen(line))){
        //printf("start of %d th is %d (i=%d, line[i]=%c %d)\n",j,start,i,line[i],line[i]);
        while((bp[i]==0)&&(i<strlen(line))&&(i-start<50)){
            a[j][i-start]=line[i];
            i++;
        }
        a[j][i-start]=0;
        if(i-start>=50) {printf("separate_line error: piece %ld does not fit in length 5000\n%s\n",j,line);}
        for(k=i-start+1;k<50;k++) a[j][k]=0;
        //printf("Break at i=%d\n",i);//getchar();
        j++;
        if(j>=MAX-1) {printf("separate_line_tiny: Incoming array too small for string fragments!\n%s\n",line);exit(1);}
        i++;
        start=i;
    }
    //printf("returning n=%d\n",j-1);getchar();
    delete[] bp;bp=NULL;
    return(j-1);
}

unsigned long int separate_large_line(char line[],char a[][30000],char bc,int MAX){
	unsigned long int i,j,n,start,k;
	n=0;
	char *bp;
	char break_char[4];
	
	bp=new char[strlen(line)+1];
	for(i=0;i<=strlen(line);i++) bp[i]=0;
	
	for(i=0;i<=3;i++) break_char[i]=0;
	break_char[0]=bc;
	
	i=0;
	while(break_char[i]!=0){
		for(j=0;j<strlen(line);j++) if(line[j]==break_char[i]) {bp[j]=1;}
		i++;
	}
	
	i=0;
	while(break_char[i]==0){
		i++;
	}
	
	start=i; j=1;
	while((line[i]!=0)&&(i<strlen(line))){
		//printf("start of %d th is %d (i=%d, line[i]=%c %d)\n",j,start,i,line[i],line[i]);
		while((bp[i]==0)&&(i<strlen(line))&&(i-start<30000)){
			a[j][i-start]=line[i];
			i++;
		}
		a[j][i-start]=0;
		if(i-start>=30000) {printf("separate_line error: piece %ld does not fit in length 30000\n%s\n",j,line);}
		for(k=i-start+1;k<30000;k++) a[j][k]=0;
		//printf("Break at i=%d\n",i);//getchar();
		j++; 
		if(j>=MAX-1) {printf("separate_large_line: Incoming array too small for string fragments!\n%s\n",line);exit(1);}
		i++;
		start=i;
	}
	//printf("returning n=%d\n",j-1);getchar();
	delete[] bp;bp=NULL;
	return(j-1);
}


unsigned long int separate_line(char line[],char a[][5000],char bc,int MAX){
	unsigned long int i,j,n,start,k;
	n=0;
	char *bp;
	char break_char[4];
	
	bp=new char[strlen(line)+1];
	for(i=0;i<=strlen(line);i++) bp[i]=0;
	
	for(i=0;i<=3;i++) break_char[i]=0;
	break_char[0]=bc;
	
	i=0;
	while(break_char[i]!=0){
		for(j=0;j<strlen(line);j++) if(line[j]==break_char[i]) {bp[j]=1;}
			i++;
	}
	
	i=0;
	while(break_char[i]==0){
		i++;
	}
	
	start=i; j=1;
	while((line[i]!=0)&&(i<strlen(line))){
		//printf("start of %d th is %d (i=%d, line[i]=%c %d)\n",j,start,i,line[i],line[i]);
		while((bp[i]==0)&&(i<strlen(line))&&(i-start<5000)){
			a[j][i-start]=line[i];
			i++;
		}
		a[j][i-start]=0;
		if(i-start>=5000) {
			printf("separate_line error: piece %ld does not fit in length 5000\n%s\n",j,line);
			getchar();
		}
		for(k=i-start+1;k<5000;k++) a[j][k]=0;
		//printf("Break at i=%d\n",i);//getchar();
		j++; 
		if(j>=MAX-1) {printf("separate_line: Incoming array too small for string fragments!\n%s\n",line);exit(1);}
		i++;
		start=i;
	}
	//printf("returning n=%d\n",j-1);getchar();
	delete[] bp;bp=NULL;
	return(j-1);
}

unsigned long int separate_line_2(char line[],char a[][7000],char bc,int MAX){
	unsigned long int i,j,n,start,k;
	n=0;
	char *bp;
	char break_char[4];
	
	bp=new char[strlen(line)+1];
	for(i=0;i<=strlen(line);i++) bp[i]=0;
	
	for(i=0;i<=3;i++) break_char[i]=0;
	break_char[0]=bc;
	
	i=0;
	while(break_char[i]!=0){
		for(j=0;j<strlen(line);j++) if(line[j]==break_char[i]) {bp[j]=1;}
        i++;
	}
	
	i=0;
	while(break_char[i]==0){
		i++;
	}
	
	start=i; j=1;
	while((line[i]!=0)&&(i<strlen(line))){
		//printf("start of %d th is %d (i=%d, line[i]=%c %d)\n",j,start,i,line[i],line[i]);
		while((bp[i]==0)&&(i<strlen(line))&&(i-start<7000)){
			a[j][i-start]=line[i];
			i++;
		}
		a[j][i-start]=0;
		if(i-start>=7000) {
			printf("separate_line error: piece %ld does not fit in length 7000\n%s\n",j,line);
			getchar();
		}
		for(k=i-start+1;k<7000;k++) a[j][k]=0;
		//printf("Break at i=%d\n",i);//getchar();
		j++;
		if(j>=MAX-1) {printf("separate_line: Incoming array too small for string fragments!\n%s\n",line);exit(1);}
		i++;
		start=i;
	}
	//printf("returning n=%d\n",j-1);getchar();
	delete[] bp;bp=NULL;
	return(j-1);
}

unsigned long int separate_line_by_string(char line[],char a[][5000],const char bc[],int MAX){
	unsigned long int i,j,n,k,ok;

	n=1;

	i=0;j=0;
	strcpy(a[1],"");
	
	while(line[i]!=0){
		if(line[i]!=bc[0]) {a[n][j]=line[i];j++;}
		else{ 
			ok=1;
			for(k=1;k<strlen(bc);k++)
				if(line[i+k]!=bc[k]) ok=0;
			if(ok) { a[n][j]=0;
					 n++; 
					 if(n>=MAX-1) {printf("Incoming array too small for string fragments!\n%s\n",line);
						 exit(1);}
                     strcpy(a[n], "");
                     i+=strlen(bc);
					 j=0;
			}
			a[n][j]=line[i];j++;
			if(j>=5000){ printf("separate_line_by_string error: piece %ld does not fit in length 5000\n%s\n",j, line);exit(1);}
		}
		i++;
	}
	a[n][j]=0;
	return(n);
}

unsigned long int separate_line_by_string_2(char line[],char a[][7000],const char bc[],int MAX){
	unsigned long int i,j,n,k,ok;
    
	n=1;
    
	i=0;j=0;
	strcpy(a[1],"");
	
	while(line[i]!=0){
		if(line[i]!=bc[0]) {a[n][j]=line[i];j++;}
		else{
			ok=1;
			for(k=1;k<strlen(bc);k++)
				if(line[i+k]!=bc[k]) ok=0;
			if(ok) { a[n][j]=0;
                n++;
                if(n>=MAX-1) {printf("Incoming array too small for string fragments!\n%s\n",line);
                    exit(1);}
                strcpy(a[n], "");
                i+=strlen(bc);
                j=0;
			}
			a[n][j]=line[i];j++;
			if(j>=7000){ printf("separate_line_by_string error: piece %ld does not fit in length 7000\n%s\n",j, line);exit(1);}
		}
		i++;
	}
	a[n][j]=0;
	return(n);
}

char *strcasestr (char *haystack, char *needle)
{
	char *p, *startn = 0, *np = 0;
	
	for (p = haystack; *p; p++) {
		if (np) {
			if (toupper(*p) == toupper(*np)) {
				if (!*++np)
					return startn;
			} else
				np = 0;
		} else if (toupper(*p) == toupper(*needle)) {
			np = needle + 1;
			startn = p;
		}
	}
	
	return 0;
}

char *export_string(char st[]){
	char *nst;
	unsigned long i;
	
	nst=new char[strlen(st)+2];
	for(i=0;i<=strlen(st);i++)
		if(((st[i]<=32)||(st[i]=='/')||(st[i]=='(')||(st[i]==')')||(st[i]=='%')||(st[i]=='.')||(st[i]==39))&&(st[i]!=0)) nst[i]=95;
		else nst[i]=st[i];
	nst[strlen(st)+1]=0;
	i=strlen(nst);
	return(nst);
}

char *export_GB_string(char st[]){
	char *nst;
	int i;
	
	nst=new char[strlen(st)+2];
	for(i=0;i<=strlen(st);i++)
		if((st[i]<=32)&&(st[i]!=0)) nst[i]=219;
		else nst[i]=st[i];
	nst[strlen(st)+1]=0;
	return(nst);
}

char *nice_string(char st[]){
	char *nst;
	int i;
	
	nst=new char[strlen(st)+2];
	for(i=0;i<=strlen(st);i++)
		if((st[i]==95)||(st[i]==219)) nst[i]=32;
		else nst[i]=st[i];
	nst[strlen(st)+1]=0;
	return(nst);
}
char *nodollar(char st[]){
	char *nst;
	int i;
	
	nst=new char[strlen(st)+2];
	for(i=0;i<=strlen(st);i++)
		if(st[i]=='$') nst[i]='-';
		else nst[i]=st[i];
	nst[strlen(st)+1]=0;
	return(nst);
}

char *no_vertical_line(char st[]){
	char *nst;
	int i;
	
	nst=new char[strlen(st)+2];
	for(i=0;i<=strlen(st);i++)
		if(st[i]=='$') nst[i]='-';
		else nst[i]=st[i];
	nst[strlen(st)+1]=0;
	return(nst);
}

unsigned int get_hist_box(unsigned int BOX,double min,double max,double x){
	if(x<min) {printf("Bad call of get_hist_box: x < min (%lg < %lg)\n",x,min);exit(1);}
	if(x>max) {printf("Bad call of get_hist_box: x > max (%lg > %lg)\n",x,max);exit(1);}
	if(min>=max) {printf("Bad call of get_hist_box: min >= max (%lg > %lg)\n",min,max);exit(1);}
	return((unsigned int)((x-min)*BOX/(max-min)));
}

unsigned int get_middle_of_box(unsigned int BOX,double min,double max,unsigned int i){
	if(min>=max) {printf("Bad call of get_hist_box: min >= max (%lg > %lg)\n",min,max);exit(1);}
	return( min + (max-min)*(i-0.5)/(double)BOX);
}

/*char *nice_GB_string(char st[]){
	char *nst;
	int i;
	
	nst=new char[strlen(st)+2];
	for(i=0;i<=strlen(st);i++)
		if(st[i]==219) nst[i]=32;
		else nst[i]=st[i];
	nst[strlen(st)+1]=0;
	return(nst);
}
*/

void sample_average_sigma(unsigned long int N, double *x,double *aver,double *szig){
	unsigned long int i;
	double av=0,sz=0;
	for(i=1;i<=N;i++) {av+=x[i];sz+=x[i]*x[i];} 
	av/=(double)N;
	sz=sqrt((sz- av*av*(double)N)/(double)(N - 1));
	*aver=av;
	*szig=sz;
}

void sample_average_variance(unsigned long int N, double *x,double *aver,double *szig){
	unsigned long int i;
	double av=0,sz=0;
	for(i=1;i<=N;i++) {av+=x[i];sz+=x[i]*x[i];} 
	av/=(double)N;
	sz=(sz- av*av*(double)N)/(double)(N - 1);
	*aver=av;
	*szig=sz;
}

double sigma_from_average_and_sum_of_squares(double av,double sumsq,unsigned long int n){
	double sz;	
    if(n<=0) return(0);
    if(n==1) return(0);
	//printf("\t\t\t\tvar=%lg\n",(sumsq-n*av*av)/(double)(n-1));
    if(sumsq<n*av*av) {
        getchar();
        return(0);
    }
	sz=sqrt((sumsq-n*av*av)/(double)(n-1));
	return(sz);
}


double get_general_median(unsigned long int N_list,unsigned long int k,double *x){
	double *xx, mednow;;
	unsigned long int i,id1,n_sm=0,n_eq=0,n_larg=0;
	
	//		S is a list of numerical values 
	//		k is the rank of the kth smallest element in S. 
	//		procedure select(k,S) 
	
	if(N_list<1) return(0); 
	if(N_list==1) return (x[0]); 
	id1=rand_int(0,N_list-1);
	
	//		choose an element a randomly from S; 
	n_sm=0;n_eq=0;n_larg=0;
	for(i=0;i<N_list;i++){
		if(x[i]<x[id1]) n_sm++;
		else if(x[i]>x[id1]) n_larg++;
		else n_eq++;
	}
	///		let S1, S2, S3 be the sequences of 
	//		elements in S respectively less than, 
	//		equal to, and greater than a; 
	if(n_sm>=k) {
		xx=new double[n_sm];
		n_sm=0;
		for(i=0;i<N_list;i++) if(x[i]<x[id1]) {xx[n_sm]=x[i];n_sm++;}
		mednow=get_general_median(n_sm,k,xx);
		delete[] xx;xx=NULL;
		return (mednow);
	}
	//if (|S1| >= k) then 
	//	return select(k,S1); 
	else {
		if(n_sm+n_eq>=k) return(x[id1]);
		//if (|S1|+|S2| >= k) then 
		//	return a; 
		else{
			xx=new double[n_larg];
			n_larg=0;
			for(i=0;i<N_list;i++) if(x[i]>x[id1]) {xx[n_larg]=x[i];n_larg++;}
			mednow=get_general_median(n_larg,k-n_sm-n_eq,xx);
			delete[] xx;xx=NULL;
			return (mednow);
		} 
		//		return select(k-|S1|-|S2|, S3); 
	}
}

double get_median(unsigned long int N_list,double *x){
	double mednow;
	if(N_list==1) return(x[0]);
	mednow=get_general_median(N_list,(unsigned long int)(N_list/2),x);
	return(mednow);
}

void get_Quartiles(unsigned long int N_list,double *x,double *q1,double *q2,double *q3){

	switch (N_list-4*(unsigned long int)(N_list/4)) {
		case 0:{
			*q1=get_general_median(N_list,(unsigned long int)(N_list/4),x);
			*q2=get_general_median(N_list,(unsigned long int)(N_list/2),x);
			*q3=get_general_median(N_list,3*(unsigned long int)(N_list/4),x);
		}
			break;
		case 1:{
			*q1=get_general_median(N_list,(unsigned long int)(N_list/4),x);
			*q2=0.5*(get_general_median(N_list,(unsigned long int)(N_list/2),x)+get_general_median(N_list,(unsigned long int)(N_list/2)+1,x));
			*q3=get_general_median(N_list,3*(unsigned long int)(N_list/4)+1,x);
		}
			break;
		case 2:{
			*q1=0.5*(get_general_median(N_list,(unsigned long int)(N_list/4),x)+get_general_median(N_list,(unsigned long int)(N_list/4)+1,x));
			*q2=get_general_median(N_list,(unsigned long int)(N_list/2),x);
			*q3=0.5*(get_general_median(N_list,3*(unsigned long int)(N_list/4),x)+get_general_median(N_list,3*(unsigned long int)(N_list/4)+1,x));
		}
			break;
		case 3:{
			*q1=get_general_median(N_list,(unsigned long int)(N_list/4)+1,x);
			*q2=0.5*(get_general_median(N_list,(unsigned long int)(N_list/2),x)+get_general_median(N_list,(unsigned long int)(N_list/2)+1,x));
			*q3=get_general_median(N_list,3*(unsigned long int)(N_list/4)+1,x);
		}
			break;
	}
}

struct Ttest_results {
	double p,t,av_a,av_b,var_a,var_b,sig_ab;
};


double get_Mann_Whitney_Wilcoxon_rank_sum_U(unsigned int nx,double x[],unsigned int ny,double y[]){
	double *xy,*rank_equals_divide;
	unsigned int i,j,st;
	unsigned long int *orig_index;
	double R1,R2,U1,U2;
	
	xy=new double[nx+ny+1];
	for(i=1;i<=nx;i++) xy[i]=x[i];
	for(i=1;i<=ny;i++) xy[nx+i]=y[i];
	orig_index=new unsigned long int[nx+ny+1];
	rank_equals_divide=new double[nx+ny+1];
	for(i=1;i<=nx+ny;i++) rank_equals_divide[i]=i;
	
	hpsort_coshuffle_indices(nx+ny, xy,orig_index);
	
	for(i=1;i<nx+ny;i++) {
		j=0;
		while((i+j+1<=nx+ny)&&(xy[i+j+1]-xy[i+j])<0.000000000000000000000001) j++;
		for(st=i;st<=i+j;st++) rank_equals_divide[st]=i+j/2.;
		i+=j;
		//if(j>0) printf("Averaged rank of %d equan values to %lg\n",j,i+j/2.);
	}
	R1=0;R2=0;
	for(i=1;i<=nx+ny;i++) {
		if(orig_index[i]<=nx) R1+=rank_equals_divide[i]; 
		else R2+=rank_equals_divide[i];
	}
	
	U1=R1-(nx*(nx+1))/2.;
	U2=R2-(ny*(ny+1))/2.;
	
	delete[] xy;xy=NULL;	
	delete[] orig_index;orig_index=NULL;	
	delete[] rank_equals_divide;rank_equals_divide=NULL;	
	
	if(U1<U2) return(U1); 
	else return(-U2);
}


double P_Cumulative_for_Normal_Distrib(double x){ // Z-score 
	double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
	
    // Save the sign of x
    int sign = 1;
    if (x < 0) sign = -1;
    x = fabs(x)/sqrt(2.0);
	
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    return 0.5*(1.0 + sign*y);
}

char bla1[4000],bla2[1000],good[5000];
char desc[25000];

short int get_line_and_check_size(FILE *f){
    if(fgets(GSE_line,MAX_LINE,f)!=NULL){
        if(strlen(GSE_line)>MAX_LINE-3) {printf("Raize GSE_line size, I have one %lu long (%s)!\n",strlen(GSE_line),GSE_line);
            getchar();
        }
        return(1);
    }
    else return(0);
}

short int get_line_and_check(FILE *f){
    
    if(fgets(GSE_line,MAX_LINE,f)!=NULL){
        if(strlen(GSE_line)>MAX_LINE-3) {printf("Raize GSE_line size, I have one %lu long (%s)!\n",strlen(GSE_line),GSE_line); exit(1);}
        
        strcpy(bla1,"");strcpy(bla2,"");strcpy(good,"");
        sscanf(GSE_line,"%s %s %s",bla1,bla2,good);
        if(strlen(bla1)>=4000-3) {printf("Raize bla1 size, I have one %lu long (%s)!\n",strlen(bla1),bla1); exit(1);}
        if(strlen(bla2)>=1000-3) {printf("Raize bla2 size, I have one %lu long (%s)!\n",strlen(bla2),bla2); exit(1);}
        if(strlen(good)>=5000-3) {printf("Raize good size, I have one %lu long (%s)!\n",strlen(good),good); exit(1);}
        return(1);
    }
    else return(0);
}

std::vector<std::string> Extract_String_Between_Delimiters_On_Same_Level(const std::string &original_string,
                                                               char opening_delimiter, // Should be a single char, like '(' or '{'
                                                               char closing_delimiter){ // Should be a single char, like '(' or '{'
    std::vector<std::string> word;
    auto b=original_string.begin(), e=original_string.end(), p=b;
    int nesting_level=0;
    while (b != e){
        if (*b == closing_delimiter){
            if (nesting_level > 0 && --nesting_level == 0){
                word.push_back(std::string(p, b));
            }
        }
        
        if (*b++ == opening_delimiter){ // anything before this preceeds a top-level parantneses
            if (nesting_level++ == 0){
                // after last
                p=b; // marks start of one top-level bracket
            }
        }
    }
    return (word);
}


std::string getLastWord( std::string str ){
    std::string token ;
    
    // read the first token from an input string stream constructed with the string in reverse
    std::istringstream( { str.rbegin(), str.rend() } ) >> token ; // #include <sstream>
    
    // return the reverse of the token that was read
    return { token.rbegin(), token.rend() } ;
}

/* Linear equation solution by Gauss-Jordan elimination. 
 a[0..n-1][0..n-1] is the input matrix. b[0..n-1] is input
 containing the right-hand side vectors. On output a is
 replaced by its matrix inverse, and b is replaced by the 
 corresponding set of solution vectors 
 */
/*
#define NRANSI
#include "nrutil.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv,temp;
	
	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
				for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
					}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
#undef SWAP
#undef NRANSI
*/

