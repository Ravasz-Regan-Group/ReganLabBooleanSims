#define NOBOWTIE_IT_MAX 1000
#define FUNCTIONAL_BOOLEAN_GATE 1

unsigned int DRAW_EXPORT;

double p_0=0.25;

double beta;
// logic of network links is OPPOSITE from sgraph! So: INCOMING links are stored under a node-s own list to make gate
// updates easyer. k_in and k_out are inveresely labeled because of this. notice update_Gate();
char attractor_szinnev[15][30];
unsigned int attractor_szinnev_MAX;
void set_attractor_szinnev(){
	sprintf(attractor_szinnev[1],"Black");
	sprintf(attractor_szinnev[2],"Blue");
	sprintf(attractor_szinnev[3],"OliveGreen");
	sprintf(attractor_szinnev[4],"Mulberry");
	sprintf(attractor_szinnev[5],"Brown");
	sprintf(attractor_szinnev[6],"Gray");
	sprintf(attractor_szinnev[7],"Cerulean");
	sprintf(attractor_szinnev[8],"Lavender");
	sprintf(attractor_szinnev[9],"Tan");
	sprintf(attractor_szinnev[10],"Pink");
	sprintf(attractor_szinnev[11],"RedOrange");
	attractor_szinnev_MAX=11;
}

double maximum_distance_between_p_Binary_string_pairs(unsigned int p){
    
    if (p % 2 == 0) {
       // printf("p=%d returning %lg\n",p,p/((double)2*(p-1)));
        return(p/((double)2*(p-1)));
    }
    else {
       // printf("p=%d returning %lg\n",p,((double)p/2.+0.5)/((double)p));
        return(((double)p/2.+0.5)/((double)p));
    }
    
    //    p=2  akkor   a max distance N
    //    p=3 kijon a 2/3 N
    //    p=4 akkor is 2/3 N
    //    p=5,6   3/5N
    //    p=7,8   4/7N
    //    p=9,10  5/9N
}



class Boolean_Gate{
public:	
	unsigned int K;
	unsigned long long N;
	unsigned short int *O;
	
	Boolean_Gate(unsigned int k,unsigned int F){
		K=k;
		N=(unsigned long long)(pow(2.,(double)K));
		generate_Rnd_Boolean_Gate(p_0);
		if(F==FUNCTIONAL_BOOLEAN_GATE) 
			guarantee_function_for_all_inputs();
	}
    
    Boolean_Gate(unsigned int k){
		K=k;
		N=(unsigned long long)(pow(2.,(double)K));
        O=new unsigned short int[N];
	}
	
	Boolean_Gate(Boolean_Gate *g_to_copy){
		unsigned long int i;
		K=g_to_copy->K;
		N=(unsigned long long)(pow(2,(double)K));
		O=new unsigned short int[N];
		for(i=0;i<N;i++) {
			O[i]=g_to_copy->O[i];
		}
	}
    
    unsigned long int get_n0(){
        unsigned long int n0=0;
        for(unsigned long int i=0;i<N;i++) if(O[i]==0) n0++;
        return(n0);
    }
    
    void knockout(){
        for(unsigned long int i=0;i<N;i++) O[i]=0;
    }
    
    void overexpress(){
        for(unsigned long int i=0;i<N;i++) O[i]=1;
    }
	
	void generate_Rnd_Boolean_Gate(double p){
		unsigned long long i;
		double r;
		O=new unsigned short int[N];
		for(i=0;i<N;i++) {
			r=rand_real(0,1);
			if(r<=p_0) O[i]=0;else O[i]=1;
		}
	}
	
    
	void add_one_rnd_input_as_it_was_0_before(){
		unsigned int K_new;
		unsigned long long N_new;
		unsigned short int *O_new;
		unsigned long long i;
		double r;
	
		K_new=K+1;
		N_new=(unsigned long long)(pow(2.,(double)K_new));
		O_new=new unsigned short int[N_new];
		for(i=0;i<N;i++) {
			O_new[i]=O[i];
		}
		for(i=N;i<2*N;i++){
			r=rand_real(0,1);
			if(r<=p_0) O_new[i]=0;else O_new[i]=1;
		}
		

		K=K_new;
		N=N_new;
		delete[] O;O=NULL;
		O=O_new;
		guarantee_function_for_all_inputs();
	}
    
    short int check_overall_sign_of_input(unsigned int inp1){
		unsigned int *sin_now;
		unsigned long long i,ones,zeros;
        
        ones=0;zeros=0;
        sin_now=new unsigned int[K+1];
        for(i=0;i<N;i++) {
			decimal_to_any_base(i,2,K,sin_now); // sin_now is the input which generates the ith output
//            if (K<inp1) {
//                for (int ii=0; ii<K; ii++) printf("%d",sin_now[ii]);
//                printf(" K=%d i=%lld O[i]=%d  inp1=%d\n",K,i,O[i],inp1);
//               // getchar();
//            }
            if(sin_now[K-inp1]==0) zeros+=O[i]; // if the inp1th input is 0, sum up output
            else ones+=O[i];
		}
        delete[] sin_now;sin_now=NULL;
        if(ones>zeros) return(1);
        else if(ones<zeros) return(-1);
        return(0);
        
	}
	
	unsigned int gate_value_for_input(unsigned int sin[]){
		unsigned long long k_1,k;
		
		k_1=0;
		k=0;
		for(k=1;k<=K;k++){
			if(sin[k]==1) k_1+=(unsigned long long)(pow(2.,(double)(K-k)));
		}
		// if(BRN->N>10) printf("k1=%lld\t returning %d\n",k_1,BRN->Gate[node_ID]->O[k_1]);
		return(O[k_1]);
	}
	
	unsigned long long gate_value_INDEX_for_input(unsigned int sin[]){
		unsigned long long k_1,k;

		k_1=0;
		k=0;
		//printf("\nstate: ");
		for(k=1;k<=K;k++){
			//printf("%d",sin[k]);
			if(sin[k]==1) k_1+=(unsigned long long)(pow(2.,(double)(K-k)));
		}
		//printf("\n");
		// if(BRN->N>10) printf("k1=%lld\t returning %d\n",k_1,BRN->Gate[node_ID]->O[k_1]);
		return(k_1);
	}
	double Entropy_of_gate(){
        double p = 0.0;
        for(unsigned long long i=0;i<N;i++)
            if(O[i]==0) p++;
        p/=(double)N;
        return(-p*log(p)-(1.-p)*log(1-p));
    }
    
	longint_ARRAY *test_function_for_all_inputs_return_bad_ones(){
		unsigned int nodifference,k,*s_n,*sin_0,*sin_1;
		unsigned long long k_0,k_1,i,l;
        longint_ARRAY *badin; badin=new longint_ARRAY();
		if(N==1) {
			if(O[0]==O[1]) {
               // printf("\t\tSingle input 1 does not affect gate\n");
                badin->add_element(1);
                return(badin);
                //getchar();
			}
            else return(badin);
		}
        
        s_n=new unsigned int[K];
        sin_0=new unsigned int[K+1];
        sin_1=new unsigned int[K+1];
        for(k=1;k<=K;k++){
            //printf("\tinput in position k=%d\n",k);
            nodifference=1;
            for(i=0;i<N/2;i++){
                //printf("i=%lld  sn=",i);
                decimal_to_any_base(i,2,K-1,s_n);
                //for(l=0;l<K-1;l++) printf("%d",s_n[l]);
                //printf("\n");
                for(l=1;l<=K;l++) {
                    if(l<k) {sin_0[l]=s_n[K-l-1];   
                             sin_1[l]=s_n[K-l-1]; 
                        //printf("l=%d sn[K-l-1]=%d\n",l,s_n[K-l-1]);
                    }
                    if(l==k){sin_0[l]=0;sin_1[l]=1;//printf("l=%d =k \n",l);
                    }
                    if(l>k){ sin_0[l]=s_n[K-l]; 
                             sin_1[l]=s_n[K-l];
                        //printf("l=%d sn[K-l]=%d\n",l,s_n[K-l]);
                    }
                }
                k_0=gate_value_INDEX_for_input(sin_0);
                k_1=gate_value_INDEX_for_input(sin_1);
                if(O[k_0]!=O[k_1]) { nodifference=0;i=N;}
                //printf("\t\t output k_0=%lld -> O=%d vs k_1=%lld ->O=%d   NODIFF=%d\n",k_0,O[k_0],k_1,O[k_1],nodifference);
            }
            
            if(nodifference==1){
               // printf("\t\tInput %d does not affect gate\n",k);getchar();
                badin->add_element(k);
            }
            //getchar();
        }
        delete[] s_n;s_n=NULL;
        delete[] sin_0;sin_0=NULL;
        delete[] sin_1;sin_1=NULL;
        return(badin);
	}
	
    unsigned int test_function_for_input(unsigned int k){
		unsigned int nodifference,*s_n,*sin_0,*sin_1;
		unsigned long long k_0,k_1,i,l;
		
		if(N==1) {
			if(O[0]==O[1]) {
                return(0);
			}
            else return(1);
		}
        
        s_n=new unsigned int[K];
        sin_0=new unsigned int[K+1];
        sin_1=new unsigned int[K+1];
        nodifference=1;
        for(i=0;i<N/2;i++){
            //printf("i=%lld  sn=",i);
            decimal_to_any_base(i,2,K-1,s_n);
            //for(l=0;l<K-1;l++) printf("%d",s_n[l]);
            //printf("\n");
            for(l=1;l<=K;l++) {
                if(l<k) {sin_0[l]=s_n[K-l-1];   
                    sin_1[l]=s_n[K-l-1]; 
                    //printf("l=%d sn[K-l-1]=%d\n",l,s_n[K-l-1]);
                }
                if(l==k){sin_0[l]=0;sin_1[l]=1;//printf("l=%d =k \n",l);
                }
                if(l>k){ sin_0[l]=s_n[K-l]; 
                    sin_1[l]=s_n[K-l];
                    //printf("l=%d sn[K-l]=%d\n",l,s_n[K-l]);
                }
            }
            k_0=gate_value_INDEX_for_input(sin_0);
            k_1=gate_value_INDEX_for_input(sin_1);
            if(O[k_0]!=O[k_1]) { nodifference=0;i=N;}
            //printf("\t\t output k_0=%lld -> O=%d vs k_1=%lld ->O=%d   NODIFF=%d\n",k_0,O[k_0],k_1,O[k_1],nodifference);
        }
        delete[] s_n;s_n=NULL;
        delete[] sin_0;sin_0=NULL;
        delete[] sin_1;sin_1=NULL;
        
        if(nodifference==1){
            return(0);
        }
        else return(1);
	}

    void randomize_function_of_input(unsigned int k,int preserve_output_for_k_value){
		unsigned int bit_difference,*s_n,*sin_0,*sin_1;
        unsigned long long r;
		unsigned long long k_0,k_1,i,l;
		
		if(N==1) {
			if(O[0]==O[1]) {
                r=rand_int(0,1);
                O[0]=r;O[1]=r;
                return;
			}
            else {
                r=rand_int(0,1);
                O[0]=r;O[1]=1-r;
                return;
            }
		}
        
        s_n=new unsigned int[K];
        sin_0=new unsigned int[K+1];
        sin_1=new unsigned int[K+1];
        bit_difference=0;
        
        // cound and erase all difference
        for(i=0;i<N/2;i++){
            //printf("i=%lld  sn=",i);
            decimal_to_any_base(i,2,K-1,s_n);
            //for(l=0;l<K-1;l++) printf("%d",s_n[l]);
            //printf("\n");
            for(l=1;l<=K;l++) {
                if(l<k) {sin_0[l]=s_n[K-l-1];
                         sin_1[l]=s_n[K-l-1];
                    //printf("l=%d sn[K-l-1]=%d\n",l,s_n[K-l-1]);
                }
                if(l==k){sin_0[l]=0;sin_1[l]=1;//printf("l=%d =k \n",l);
                }
                if(l>k){ sin_0[l]=s_n[K-l];
                         sin_1[l]=s_n[K-l];
                    //printf("l=%d sn[K-l]=%d\n",l,s_n[K-l]);
                }
            }
            k_0=gate_value_INDEX_for_input(sin_0);
            k_1=gate_value_INDEX_for_input(sin_1);
            if(O[k_0]!=O[k_1]) {
               // printf("\t\t erase: output k_0=%lld -> O=%d vs k_1=%lld ->O=%d \n",k_0,O[k_0],k_1,O[k_1]);
                bit_difference++;
                if(preserve_output_for_k_value==0) O[k_1]=O[k_0];
                else O[k_0]=O[k_1];
               // printf("\t\t\t goes to output k_0=%lld -> O=%d vs k_1=%lld ->O=%d \n",k_0,O[k_0],k_1,O[k_1]);
            }
        }
        
        // pick an equal number of random places to add a difference
        
        longint_ARRAY *place_diff;
        place_diff=new longint_ARRAY();
        while (bit_difference) {
            r=rand_int(0,N/2-1);
            if(place_diff->check_for_element(r)==0){place_diff->add_element(r);bit_difference--;}
        }
        
        // add difference in random spots
        
        for (int kk=1; kk<=place_diff->N; kk++) {
            i=place_diff->A[kk];
            //printf("Inserting diff. at i=%d\n",i);
            decimal_to_any_base(i,2,K-1,s_n);
            //for(l=0;l<K-1;l++) printf("%d",s_n[l]);
            //printf("\n");
            for(l=1;l<=K;l++) {
                if(l<k) {sin_0[l]=s_n[K-l-1];
                    sin_1[l]=s_n[K-l-1];
                    //printf("l=%d sn[K-l-1]=%d\n",l,s_n[K-l-1]);
                }
                if(l==k){sin_0[l]=0;sin_1[l]=1;//printf("l=%d =k \n",l);
                }
                if(l>k){ sin_0[l]=s_n[K-l];
                    sin_1[l]=s_n[K-l];
                    //printf("l=%d sn[K-l]=%d\n",l,s_n[K-l]);
                }
            }
            k_0=gate_value_INDEX_for_input(sin_0);
            k_1=gate_value_INDEX_for_input(sin_1);
            //printf("\t\t new random: output k_0=%lld -> O=%d vs k_1=%lld ->O=%d \n",k_0,O[k_0],k_1,O[k_1]);
            if(preserve_output_for_k_value==0) O[k_1]=1-O[k_0];
            else O[k_0]=1-O[k_1];
            //printf("\t\t\t goes to output k_0=%lld -> O=%d vs k_1=%lld ->O=%d \n",k_0,O[k_0],k_1,O[k_1]);
        }
        
        delete place_diff;place_diff=NULL;
        delete[] s_n;s_n=NULL;
	}

    void guarantee_function_for_all_inputs(){
		unsigned int allok,nodifference,k,*s_n,*sin_0,*sin_1;
		unsigned long long k_0,k_1,i,l;
		
		if(N==1) {
			if(O[0]==O[1]) {
				k_1=rand_int(0, 1);
				O[k_1]=1-O[k_1];
			}
			return;
		}
		do{ allok=1;
			//printf("Guaranteeing function for gate: ");
			//for(i=0;i<N;i++) printf("%d ",O[i]);
		//	printf("\n");

			s_n=new unsigned int[K];
			sin_0=new unsigned int[K+1];
			sin_1=new unsigned int[K+1];
			for(k=1;k<=K;k++){
				//printf("\tinput in position k=%d\n",k);
				nodifference=1;
				for(i=0;i<N/2;i++){
					//printf("i=%lld  sn=",i);
					decimal_to_any_base(i,2,K-1,s_n);
					//for(l=0;l<K-1;l++) printf("%d",s_n[l]);
					//printf("\n");
					for(l=1;l<=K;l++) {
						if(l<k) {sin_0[l]=s_n[K-l-1];   
						         sin_1[l]=s_n[K-l-1]; 
									//printf("l=%d sn[K-l-1]=%d\n",l,s_n[K-l-1]);
						}
						if(l==k){sin_0[l]=0;sin_1[l]=1;//printf("l=%d =k \n",l);
						}
						if(l>k){ sin_0[l]=s_n[K-l]; 
								 sin_1[l]=s_n[K-l];
									//printf("l=%d sn[K-l]=%d\n",l,s_n[K-l]);
						}
					}
					k_0=gate_value_INDEX_for_input(sin_0);
					k_1=gate_value_INDEX_for_input(sin_1);
					if(O[k_0]!=O[k_1]) { nodifference=0;i=N;}
					//printf("\t\t output k_0=%lld -> O=%d vs k_1=%lld ->O=%d   NODIFF=%d\n",k_0,O[k_0],k_1,O[k_1],nodifference);
				}
				
				if(nodifference==1){
					i=rand_int(0, N/2);
					decimal_to_any_base(i,2,K-1,s_n);
					for(l=1;l<=K;l++) {
						if(l<k) {sin_0[l]=s_n[K-l-1];   sin_1[l]=s_n[K-l-1]; //printf("l=%d sn[K-l-1]=%d\n",l,s_n[K-l-1]);
						}
						if(l==k){sin_0[l]=0;sin_1[l]=1;//printf("l=%d =k \n",l);
						}
						if(l>k){ sin_0[l]=s_n[K-l]; sin_1[l]=s_n[K-l];//printf("l=%d sn[K-l]=%d\n",l,s_n[K-l]);
						}
					}
					k_0=gate_value_INDEX_for_input(sin_0);
					k_1=gate_value_INDEX_for_input(sin_1);
					l=rand_int(0,1);
					if(l==0) {O[k_0]=1-O[k_0]; //printf("\tFlipping  position %lld from %d to %d\n",k_0,1-O[k_0],O[k_0]);
					}
					else {O[k_1]=1-O[k_1];//printf("\tFlipping  position %lld from %d to %d\n",k_1,1-O[k_1],O[k_1]);
					}
					allok=0;
					//getchar();
				}
				//getchar();
			}
			delete[] s_n;s_n=NULL;
            delete[] sin_0;sin_0=NULL;
            delete[] sin_1;sin_1=NULL;
		}
		while(allok==0);	
	}
    
   	~Boolean_Gate(){
		delete[] O;O=NULL;	
	}
};

typedef Boolean_Gate *Boolean_Gate_list;

class Boolean_RegNet_METADATA{
public :
	unsigned int N_NODE;
    char MOD_Assignment_name[600];
	charpointer *Node_name;
	longint_ARRAY_list *GL;
    unsigned int *Module_ID,Module_NR;
    unsigned int *BIG_ID;
    double *x,*y,*r,*g,*b;
	// GOAL: all genes assoc by ENTREZ. Start from species, go to all homologues.
	
	Boolean_RegNet_METADATA(unsigned int n){
		unsigned int i;
		N_NODE=n;
        Module_NR=0;
		GL=new longint_ARRAY_list[N_NODE+1];
		Node_name=new charpointer[N_NODE+1];
		x=new double[N_NODE+1];
        y=new double[N_NODE+1];

        for(i=1;i<=N_NODE;i++){
			GL[i]=new longint_ARRAY();
			Node_name[i]=new char[50];
		}
        Module_ID=NULL;
        BIG_ID=NULL;
	}
	
	~Boolean_RegNet_METADATA(){
		for(int i=1;i<=N_NODE;i++){
            delete GL[i];GL[i]=NULL;
            delete[] Node_name[i];Node_name[i]=NULL;
        }
        delete[] GL;GL=NULL;
        delete[] Node_name;Node_name=NULL;
        delete[] Module_ID;Module_ID=NULL;
        delete[] x;x=NULL;
        delete[] y;y=NULL;
        if (BIG_ID!=NULL) {
             delete[] BIG_ID;BIG_ID=NULL;
        }
 	}
    
	void list_node_names(){
        unsigned int i;
        if(Module_ID==NULL)
            for(i=1;i<=N_NODE;i++)
                printf("%d   %s\n",i,Node_name[i]);
        else for(i=1;i<=N_NODE;i++)
            printf("%d   %s -> M-%d\n",i,Node_name[i],Module_ID[i]);
        printf("\n");
    }
    
    void print_positions(){
        unsigned int i;
        for(i=1;i<=N_NODE;i++)
            printf("%d   %s  %lf  %lf\n",i,Node_name[i],x[i],y[i]);
        printf("\n");
    }
    unsigned int get_node_ID(const char nm[]){
        unsigned int i;
        for(i=1;i<=N_NODE;i++)
            if(!strcmp(Node_name[i],nm)) return(i);
        return(0);
    }
    unsigned int get_node_ID_BIGID(unsigned long int ibig){
        unsigned int i;
        for(i=1;i<=N_NODE;i++)
            if(BIG_ID[i]==ibig) return(i);
        return(0);
    }
    

};

class Boolean_RegNet{
public:	
	char path[600];
	char name[600];
	unsigned int N,K,Evolved,NoBowTie,RND_seed;
	sgraph *BN;
	Boolean_Gate_list *Gate;
	Boolean_RegNet_METADATA *MDAT;
	longint_ARRAY_list *Draw_in_Modules;
    longint_ARRAY_pair *Module_Order;
    longint_ARRAY *active_inputs,*active_inputs_noGF,*silenced_inputs;
    std::vector<std::string> Module_Names;
    
	Boolean_RegNet(unsigned int n,unsigned int k,unsigned int rnds){
		unsigned int i;
		N=n;K=k;MDAT=NULL;
		
		RND_seed=rnds;
		sprintf(path, "Rnd_NoBowTie_BRN");
		strcpy(name,"");sprintf(name,"Rnd_BRN_%d_%d_seed-%d",N,K,RND_seed);
		if(!import_Boolean_RegNet()){
            printf("generating Rnd_NoBowTie_BRN network N=%d, K=%d, Seed=%d\n",N,K,RND_seed);
            generate_Rnd_Boolean_RegNet(RND_seed);
            Gate=new Boolean_Gate_list[N+1];
            for(i=1;i<=N;i++) {
                Gate[i]=new Boolean_Gate((unsigned int)(BN->node[i].k_out),FUNCTIONAL_BOOLEAN_GATE);
                //printf("Gate [%d] = ",i);
                //for(k=0;k<Gate[i]->N;k++)
                //	printf("%d ",Gate[i]->O[k]);
                //printf("\n");
            }
            Evolved=0;
            export_Boolean_RegNet();
            
            //evolve_Rnd_Boolean_RegNet();
            //sprintf(name,"");sprintf(name,"%s_%d_%d",nm,N,K);
            //export_Boolean_RegNet(1);
        }
	}
	
	//Boolean_RegNet(char nm[],unsigned int n,unsigned int RND_seed){
	//	sprintf(name,"");sprintf(name,"%s",nm);
	//	N=n;MDAT=NULL;
		//import_Boolean_RegNet();
	//}
	
	Boolean_RegNet(const char pn[],const char nm[],unsigned int n){
		// no-bowtie not guaranteed, empty shell
		unsigned int i;
		strcpy(path,"");sprintf(path,"%s",pn);
		strcpy(name,"");sprintf(name,"%s",nm);
		N=n;MDAT=NULL;
		K=0;
		NoBowTie=0;
		BN=new sgraph(N,name);
		Gate=new Boolean_Gate_list[N+1];
		for(i=1;i<=N;i++) Gate[i]=NULL;
	}
    
    Boolean_RegNet(Boolean_RegNet *A,const char subset[]){
        // no-bowtie not guaranteed, empty shell
        unsigned int i,mid,j;
        FILE *f;
        char fn[600],bla[100];
        
        strcpy(path,"");sprintf(path,"%s/%s_SUB_%s",A->path,A->name,subset);
        strcpy(name,"");sprintf(name,"%s_SUB_%s",A->name,subset);
        sprintf(fn,"mkdir %s/%s\n",MODEL_Directory_system,path);
        system(fn);
        
        sprintf(fn,"%s/%s/%s_%s_%s.txt",MODEL_Directory,A->path,A->name,A->MDAT->MOD_Assignment_name,subset);
        printf("%s\n",fn);//getchar();
        f=fopen(fn,"r");
        if(f==NULL) {printf("Module assignment %s not found\n%s\t",A->MDAT->MOD_Assignment_name,fn);exit(1);}
        longint_ARRAY *exclude;
        exclude=new longint_ARRAY();
        
        int vege=1;
        while (vege!=EOF) {
            vege=fscanf(f,"%s %d %d",bla,&mid,&j);
            unsigned long int nid=A->MDAT->get_node_ID(bla);
            if(nid>0){
                exclude->add_element(nid); // A id to be excluded from subset
            }
        }
        fclose(f);
        
        N=A->N-(unsigned int)exclude->N;
        
        K=0;
        NoBowTie=0;
        BN=new sgraph(N,name);
        Gate=new Boolean_Gate_list[N+1];
        for(i=1;i<=N;i++) Gate[i]=NULL;
       
        MDAT=new Boolean_RegNet_METADATA(N);
        MDAT->Module_ID=new unsigned int[N+1];
        MDAT->BIG_ID=new unsigned int[N+1];
        strcpy(MDAT->MOD_Assignment_name,A->MDAT->MOD_Assignment_name);
        
        int sub_index=0;
        for(unsigned int iA=1;iA<=A->N;iA++)
            if(exclude->check_for_element(iA)==0){
                sub_index++;
                MDAT->BIG_ID[sub_index]=iA;
                strcpy(MDAT->Node_name[sub_index],A->MDAT->Node_name[iA]);
               // if(MDAT->GL[iA]!=NULL)
               //     MDAT->GL[sub_index]->add_longint_ARRAY(MDAT->GL[iA]);
                MDAT->Module_ID[sub_index]=A->MDAT->Module_ID[iA];
                if(A->MDAT->Module_ID[iA] > MDAT->Module_NR)
                        MDAT->Module_NR=A->MDAT->Module_ID[iA];
                MDAT->x[sub_index]=A->MDAT->x[iA];
                MDAT->y[sub_index]=A->MDAT->y[iA];
        }
        
        Module_Order=new longint_ARRAY_pair();
        for(unsigned long int k=1;k<=A->Module_Order->N;k++){
            if (A->Module_Order->A[k]==0){
                unsigned long int nidok = MDAT->get_node_ID_BIGID(A->Module_Order->B[k]);
                if (nidok != 0) Module_Order->add_element(0,nidok);
            }
            else{   Module_Order->add_element(A->Module_Order->A[k],0);
                    // some modules may have 0 nodes in the subsetted network
            }
        }
        
        Draw_in_Modules=new longint_ARRAY_list[MDAT->Module_NR+1];
        for (int i=0; i<=MDAT->Module_NR; i++) {
            Draw_in_Modules[i]=new longint_ARRAY();
            for(unsigned long int k=1;k<=A->Draw_in_Modules[i]->N;k++) {
                unsigned long int nidok = MDAT->get_node_ID_BIGID(A->Draw_in_Modules[i]->A[k]);
                if (nidok != 0)
                    Draw_in_Modules[i]->add_element(nidok);
            }
        }
    }
	
    Boolean_RegNet(Boolean_RegNet *A,const char mutantnode[],int onoff){
        // no-bowtie not guaranteed, empty shell
        char fn[600];
        
        if (onoff==0) {
            strcpy(path,"");sprintf(path,"%s/%s_KO_%s",A->path,A->name,mutantnode);
            strcpy(name,"");sprintf(name,"%s_KO_%s",A->name,mutantnode);
        }
        sprintf(fn,"mkdir %s/%s\n",MODEL_Directory_system,path);
        system(fn);
        
        
        N=A->N;
        
        K=0;
        NoBowTie=0;
        BN=new sgraph(N,name);
        Gate=new Boolean_Gate_list[N+1];
        for(unsigned int i=1;i<=N;i++) Gate[i]=NULL;
        
        MDAT=new Boolean_RegNet_METADATA(N);
        MDAT->Module_ID=new unsigned int[N+1];
        MDAT->BIG_ID=new unsigned int[N+1];
        strcpy(MDAT->MOD_Assignment_name,A->MDAT->MOD_Assignment_name);
        
        for(unsigned int iA=1;iA<=A->N;iA++){
                MDAT->BIG_ID[iA]=iA;
                strcpy(MDAT->Node_name[iA],A->MDAT->Node_name[iA]);
                // if(MDAT->GL[iA]!=NULL)
                //     MDAT->GL[sub_index]->add_longint_ARRAY(MDAT->GL[iA]);
                MDAT->Module_ID[iA]=A->MDAT->Module_ID[iA];
                if(A->MDAT->Module_ID[iA] > MDAT->Module_NR)
                    MDAT->Module_NR=A->MDAT->Module_ID[iA];
                MDAT->x[iA]=A->MDAT->x[iA];
                MDAT->y[iA]=A->MDAT->y[iA];
            }
        
        Module_Order=new longint_ARRAY_pair();
        for(unsigned long int k=1;k<=A->Module_Order->N;k++){
            if (A->Module_Order->A[k]==0){
                unsigned long int nidok = MDAT->get_node_ID_BIGID(A->Module_Order->B[k]);
                if (nidok != 0) Module_Order->add_element(0,nidok);
            }
            else{   Module_Order->add_element(A->Module_Order->A[k],0);
                // some modules may have 0 nodes in the subsetted network
            }
        }
        
        Draw_in_Modules=new longint_ARRAY_list[MDAT->Module_NR+1];
        for (int i=0; i<=MDAT->Module_NR; i++) {
            Draw_in_Modules[i]=new longint_ARRAY();
            for(unsigned long int k=1;k<=A->Draw_in_Modules[i]->N;k++) {
                unsigned long int nidok = MDAT->get_node_ID_BIGID(A->Draw_in_Modules[i]->A[k]);
                if (nidok != 0)
                    Draw_in_Modules[i]->add_element(nidok);
            }
        }
    }

    int *generate_best_canalyzing_value_for_single_input(unsigned long int BRN_index,unsigned long int inputnode){
        int *can,*outside;
        snode::slinklist::slink *sc2;
        
        can=new int[N+1];
        outside=new int[N+1];
        for(int i=0;i<=N;i++) {can[i]=-1;outside[i]=0;} // not a neighbor

        sc2=BN->node[BRN_index].Li.first_link_pt;
        while(sc2!=NULL){
            if(sc2->sneighbor == inputnode){
                outside[0]++;
                outside[sc2->sneighbor]=outside[0];
            }
            else can[sc2->sneighbor]=-2; // non-canalyzed neighbor
            sc2=sc2->next(sc2);
        }
        unsigned int Ncomb=(unsigned int)pow(2,outside[0]);
        unsigned short int *s_outmodlinks;
        
        double *H_left_for_input; H_left_for_input=new double[Ncomb+1];
        short int *i_OK_input; i_OK_input=new short int[Ncomb+1];

        for(unsigned long long int i=0;i<Ncomb;i++){
            s_outmodlinks=new unsigned short int[Ncomb+1];
            decimal_to_any_base(i,2,outside[0], s_outmodlinks);
          //  printf("i=%lld\t Canalyzing values set to:\n",i);
            for(int k=1;k<=N;k++)
                if(outside[k]>0){
                    can[k]=s_outmodlinks[outside[k]-1];
          //          printf("\t k=%d (%s) value = %d\n",k,BRN->MDAT->Node_name[k],can[k]);
                }
            longint_ARRAY *gl;
            gl=new longint_ARRAY();
            Boolean_Gate *G_test=get_PARTIAL_gate(BRN_index,can,gl);
            
        //    printf("Partial gate returned is:\n");
        //    BRN->print_PARTIAL_gate(BRN_index,gl,G_test);
            
            H_left_for_input[i]=G_test->Entropy_of_gate();
            if(G_test->test_function_for_all_inputs_return_bad_ones()->N>0){
                i_OK_input[i]=0;
           //     printf("\t\t\t NO GO!\n");
            }
            else {i_OK_input[i]=1;
            }
            // need to check how many of G_test's inputs are real,
            // order all combinations by largest real inputs
            // if more than 1 winner: select gate with largest H. if more than 1 winner: random choice.
            delete[] s_outmodlinks; s_outmodlinks=NULL;
            delete gl;gl=NULL;
        }
        
        double hmax=0; long long int i_choice=-1;
        for(unsigned long long int i=0;i<Ncomb;i++){
            if(i_OK_input[i]==1) if(H_left_for_input[i]>hmax) {hmax=H_left_for_input[i];i_choice=i;}
        }
        if (i_choice==-1) {
          //  printf("There are no way of setting the out-of-module inputs and still guaranteeing function of all internal links!\n");
          //  BRN->print_gate(BRN_index);
            //getchar();
            for(unsigned long long int i=0;i<Ncomb;i++)
                if(H_left_for_input[i]>hmax) {hmax=H_left_for_input[i];i_choice=i;}
            if(i_choice==-1) {
              //  printf("There is no regulation left for this node! locked in with canal. inputs all 0.\n");
                i_choice=0;
            }
        }
        
        s_outmodlinks=new unsigned short int[Ncomb+1];
        decimal_to_any_base(i_choice,2,outside[0], s_outmodlinks);
       // printf("Final choice: %lld:\n",i_choice);
        for(int k=1;k<=N;k++)
            if(outside[k]>0){
                can[k]=s_outmodlinks[outside[k]-1];
               // printf("\t k=%d (%s) value = %d\n",k,BRN->MDAT->Node_name[k],can[k]);
            }
        delete[] s_outmodlinks; s_outmodlinks=NULL;
        
        s_outmodlinks=new unsigned short int[Ncomb+1];
        decimal_to_any_base(i_choice,2,outside[0], s_outmodlinks);
        for(int k=1;k<=N;k++)
            if(outside[k]>0){
                can[k]=s_outmodlinks[outside[k]-1];
               // printf("\t k=%d (%s) canalyzing value settled at = %d\n",k,BRN->MDAT->Node_name[k],can[k]);
            }
    
        delete[] s_outmodlinks;s_outmodlinks=NULL;
        delete[] outside;outside=NULL;
        delete[] H_left_for_input;H_left_for_input=NULL;
        delete[] i_OK_input;i_OK_input=NULL;
        return(can);
    }
    
    Boolean_RegNet *mutant_network(const char how[],int howmany){
        Boolean_RegNet *mut;
        snode::slinklist::slink *sc2;
        char namenow[300];
        int *canalyze;
        longint_ARRAY *gl;
        
        
        strcpy(namenow,"");sprintf(namenow,"%s_%s_mutant-%d",name,how,howmany);
        mut = new Boolean_RegNet(path,namenow,N);
        
        // make shell with identical metadata
        mut->N=N;
        mut->K=0;
        mut->NoBowTie=0;
        mut->BN=new sgraph(N,namenow);
        mut->Gate=new Boolean_Gate_list[N+1];
        for(unsigned int i=1;i<=mut->N;i++) mut->Gate[i]=NULL;
        
        mut->MDAT=new Boolean_RegNet_METADATA(N);
        mut->MDAT->Module_ID=new unsigned int[N+1];
        mut->MDAT->BIG_ID=new unsigned int[N+1];
        strcpy(mut->MDAT->MOD_Assignment_name,MDAT->MOD_Assignment_name);
        
        for(unsigned int iA=1;iA<=N;iA++){
            mut->MDAT->BIG_ID[iA]=iA;
            strcpy(mut->MDAT->Node_name[iA],MDAT->Node_name[iA]);
            mut->MDAT->Module_ID[iA]=MDAT->Module_ID[iA];
            if(MDAT->Module_ID[iA] > mut->MDAT->Module_NR) mut->MDAT->Module_NR=MDAT->Module_ID[iA];
            mut->MDAT->x[iA]=MDAT->x[iA];
            mut->MDAT->y[iA]=MDAT->y[iA];
        }
        
        mut->Module_Order=new longint_ARRAY_pair();
        for(unsigned long int k=1;k<=Module_Order->N;k++){
            if (Module_Order->A[k]==0){
                unsigned long int nidok = mut->MDAT->get_node_ID_BIGID(Module_Order->B[k]);
                if (nidok != 0) mut->Module_Order->add_element(0,nidok);
            }
            else{   mut->Module_Order->add_element(Module_Order->A[k],0);
                // some modules may have 0 nodes in the subsetted network
            }
        }
        
        mut->Draw_in_Modules=new longint_ARRAY_list[mut->MDAT->Module_NR+1];
        for (int i=0; i<=mut->MDAT->Module_NR; i++) {
            mut->Draw_in_Modules[i]=new longint_ARRAY();
            for(unsigned long int k=1;k<=Draw_in_Modules[i]->N;k++) {
                unsigned long int nidok = mut->MDAT->get_node_ID_BIGID(Draw_in_Modules[i]->A[k]);
                if (nidok != 0)
                    mut->Draw_in_Modules[i]->add_element(nidok);
            }
        }
        
        // populate shell with mutant copy
        
        Boolean_Gate *G_test;
        if(!strcmp(how,"remove_links")){
            long long unsigned int L_total = BN->get_stored_link_number();
            long long unsigned int counter = 0;
            longint_ARRAY *targetlink_index;
            targetlink_index = new longint_ARRAY();
            for (unsigned int i=1; i<=howmany; i++) {
                targetlink_index -> add_element(rand_int(1, L_total));
            }
            
            longint_ARRAY_pair *targetnodes;
            targetnodes = new longint_ARRAY_pair();
            
            int linknr=0;
            for(int j=1;j<=N;j++){
                sc2=BN->node[j].Li.first_link_pt;
                linknr=0;
                while(sc2!=NULL){
                    counter++;
                    if (targetlink_index->check_for_element(counter)) {
                        targetnodes->add_element(j,sc2->sneighbor);
                    }
                    sc2=sc2->next(sc2);
                }
            }
            
            for(int j=1;j<=N;j++){
                if(targetnodes->check_for_element_A(j)==0){
                    mut->Gate[j]=new Boolean_Gate(Gate[j]);
                    sc2=BN->node[j].Li.first_link_pt;
                    gl=new longint_ARRAY();
                    while(sc2!=NULL){
                        gl->add_element(sc2->sneighbor);
                        sc2=sc2->next(sc2);
                    }
                    for(unsigned long int l=0;l<gl->N;l++){
                        mut->BN->add_link(j,gl->A[gl->N-l]);
                        mut->K++;
                    }
                    delete gl;gl=NULL;
                   // print_gate(j); mut->print_gate(j);getchar();
                }
            }
            for(unsigned int alter=1;alter<=targetnodes->N;alter++){
                canalyze=generate_best_canalyzing_value_for_single_input(targetnodes->A[alter],targetnodes->B[alter]);
                gl=new longint_ARRAY();
                G_test=get_PARTIAL_gate(targetnodes->A[alter],canalyze,gl);
                mut->Gate[targetnodes->A[alter]]=G_test;
                
                for(unsigned long int l=0;l<gl->N;l++){
                    mut->BN->add_link(targetnodes->A[alter],gl->A[gl->N-l]);
                    mut->K++;
                }
                if(gl->N==0)  mut->BN->add_link(targetnodes->A[alter],targetnodes->A[alter]);
              //  print_gate(targetnode); mut->print_gate(targetnode);getchar();
               
                delete[] canalyze;canalyze=NULL;
                delete gl;gl=NULL;
            }
        }
        
        if(!strcmp(how,"gate_oneoutput_flips")){
            longint_ARRAY *targetnodes;
            targetnodes = new longint_ARRAY();
            for (unsigned int i=1; i<=howmany; i++) {
                targetnodes -> add_element(rand_int(1, N));
            }
            
            for(int j=1;j<=N;j++){
                mut->Gate[j]=new Boolean_Gate(Gate[j]);
                sc2=BN->node[j].Li.first_link_pt;
                gl=new longint_ARRAY();
                while(sc2!=NULL){
                    gl->add_element(sc2->sneighbor);
                    sc2=sc2->next(sc2);
                }
                for(unsigned long int l=0;l<gl->N;l++){
                    mut->BN->add_link(j,gl->A[gl->N-l]);
                    mut->K++;
                }
                delete gl;gl=NULL;
                
                if(targetnodes->check_for_element(j)){
                    unsigned long int p_to_alter = rand_int(0, mut->Gate[j]->N-1);
                    mut->Gate[j]->O[p_to_alter] = 1 - Gate[j]->O[p_to_alter];
                }
            }
        }
        
        if(!strcmp(how,"lock_nodes_ONOFF")){
            longint_ARRAY *targetnodes;
            targetnodes = new longint_ARRAY();
            for (unsigned int i=1; i<=howmany; i++) {
                targetnodes -> add_element(rand_int(1, N));
            }
            
            for(int j=1;j<=N;j++){
                mut->Gate[j]=new Boolean_Gate(Gate[j]);
                sc2=BN->node[j].Li.first_link_pt;
                gl=new longint_ARRAY();
                while(sc2!=NULL){
                    gl->add_element(sc2->sneighbor);
                    sc2=sc2->next(sc2);
                }
                for(unsigned long int l=0;l<gl->N;l++){
                    mut->BN->add_link(j,gl->A[gl->N-l]);
                    mut->K++;
                }
                delete gl;gl=NULL;
                
                if(targetnodes->check_for_element(j)){
                    double p_to_alter = rand_real(0, 1);
                    for(unsigned long int k=0;k<mut->Gate[j]->N;k++)
                        if (p_to_alter<0.5) mut->Gate[j]->O[k] = 0; else mut->Gate[j]->O[k] = 1;
                }
            }
        }
        return mut;
    }

  /*
    Boolean_RegNet(Boolean_RegNet *A){
        // no-bowtie not guaranteed, empty shell
        char fn[600];
        
        Boolean_Gate_list *Gate;
        Boolean_RegNet_METADATA *MDAT;
        longint_ARRAY_list *Draw_in_Modules;
        longint_ARRAY_pair *Module_Order;
        longint_ARRAY *active_inputs,*active_inputs_noGF,*silenced_inputs;
        std::vector<std::string> Module_Names;
        
        
        strcpy(path,"");sprintf(path,"%s",A->path);
        strcpy(name,"");sprintf(name,"%s",A->name);
        
        sprintf(fn,"mkdir %s/%s\n",MODEL_Directory_system,path);
        system(fn);
        
        N=A->N; K=A->K; NoBowTie=A->NoBowTie;
        Evolved=A->Evolved; RND_seed=A->RND_seed;
        
        BN=new sgraph(A->BN);
        Gate=new Boolean_Gate_list[N+1];
        for(unsigned int i=1;i<=N;i++)
            Gate[i]=new Boolean_Gate(A->Gate[i]);
        
        MDAT=new Boolean_RegNet_METADATA(N);
        MDAT->Module_ID=new unsigned int[N+1];
        MDAT->BIG_ID=new unsigned int[N+1];
        strcpy(MDAT->MOD_Assignment_name,A->MDAT->MOD_Assignment_name);
        
        for(unsigned int iA=1;iA<=A->N;iA++){
            MDAT->BIG_ID[iA]=iA;
            strcpy(MDAT->Node_name[iA],A->MDAT->Node_name[iA]);
            MDAT->Module_ID[iA]=A->MDAT->Module_ID[iA];
            if(A->MDAT->Module_ID[iA] > MDAT->Module_NR)
                MDAT->Module_NR=A->MDAT->Module_ID[iA];
            MDAT->x[iA]=A->MDAT->x[iA];
            MDAT->y[iA]=A->MDAT->y[iA];
        }
        
        Module_Order=new longint_ARRAY_pair();
        for(unsigned long int k=1;k<=A->Module_Order->N;k++){
            if (A->Module_Order->A[k]==0){
                unsigned long int nidok = MDAT->get_node_ID_BIGID(A->Module_Order->B[k]);
                if (nidok != 0) Module_Order->add_element(0,nidok);
            }
            else{   Module_Order->add_element(A->Module_Order->A[k],0);
                // some modules may have 0 nodes in the subsetted network
            }
        }
        
        Draw_in_Modules=new longint_ARRAY_list[MDAT->Module_NR+1];
        for (int i=1; i<=MDAT->Module_NR; i++) {
            Draw_in_Modules[i]=new longint_ARRAY();
            for(unsigned long int k=1;k<=A->Draw_in_Modules[i]->N;k++) {
                unsigned long int nidok = MDAT->get_node_ID_BIGID(A->Draw_in_Modules[i]->A[k]);
                if (nidok != 0)
                    Draw_in_Modules[i]->add_element(nidok);
            }
        }
    }
  */
    
    Boolean_RegNet(const char modelname[]){
		char fn[600];
		FILE *f;
		unsigned int i,j,k,b;
		char * pch;
        char attr_list[500][500];
        unsigned long int Nf,kg;
        
        K=0;
		strcpy(path,"");sprintf(path,"%s",modelname);
		strcpy(name,"");sprintf(name,"%s",modelname);
		
		sprintf(fn,"ls %s/%s/%s/%s_Fine/ > %s/%s/%s_molecule_list.txt",MODEL_Directory_system, path,path,path,MODEL_Directory_system, path,name);
        system(fn);
        strcpy(fn, "");
 		sprintf(fn,"%s/%s/%s_molecule_list.txt",MODEL_Directory, path,name);
       
        f=fopen(fn,"r");
        if(f==NULL) {printf("Network %s not found\n",fn);exit(1);} 
		N=0;
        while (fgets(GSE_line,MAX_LINE,f)!=NULL) {
           // printf("%s --> ",GSE_line);
            pch = strstr (GSE_line,".csv");
          //  printf("pch = %s\n",pch);
            strncpy (pch,"",4);
            if(GSE_line[0]>0) N++;   
       //      printf("%s\n",GSE_line);
        }
        printf("Number of nodes %u\n",N);
        BN=new sgraph(N,name);
        Gate=new Boolean_Gate_list[N+1];
        MDAT=new Boolean_RegNet_METADATA(N);
        for(i=1;i<=N;i++) Gate[i]=NULL;
        fclose(f);
        
        f=fopen(fn,"r");
        i=0;
        while (fgets(GSE_line,MAX_LINE,f)!=NULL) {
            pch = strstr (GSE_line,".csv");
           // printf("pch = %s\n",pch);
            strncpy (pch,"",4);
            if(GSE_line[0]>0) i++;  
            strcpy(MDAT->Node_name[i],GSE_line);
        }
        fclose(f);
       // getchar();
        pch=NULL;
       for(i=1;i<=N;i++){
            sprintf(fn,"%s/%s/%s/%s_Fine/%s.csv",MODEL_Directory, path,path,path,MDAT->Node_name[i]);
            f=fopen(fn,"r");
            if(f==NULL) {printf("Gate .csv file for node %s not found\n",MDAT->Node_name[i]);exit(1);} 
            fgets(GSE_line,MAX_LINE,f);
           Nf=separate_line_smallchunks(GSE_line,attr_list,DEFAULT_DELIMITER,500);
           attr_list[Nf][strlen(attr_list[Nf])-1]=0;
           if(attr_list[Nf][strlen(attr_list[Nf])-1] == '\r')  attr_list[Nf][strlen(attr_list[Nf])-1]=0;
           if(attr_list[Nf][strlen(attr_list[Nf])-1] == '\n')  attr_list[Nf][strlen(attr_list[Nf])-1]=0;
           
           if(strcmp(attr_list[Nf],MDAT->Node_name[i])) {
               printf("Error in gate of node %s: output node %s in csv file does not match file name!\n",MDAT->Node_name[i], attr_list[Nf]); getchar(); exit(1);
           }
           Nf--; // the last one is the node itself
           for(k=1;k<=Nf;k++){
                j=MDAT->get_node_ID(attr_list[Nf-k+1]);
                if(j>0) {
                    BN->add_link(i,j);K++;
                    //printf("%u %s <- %u %s\n",i,MDAT->Node_name[i],j,MDAT->Node_name[j]);
                }
                else {printf("Node name %s not recorded, error in Boolean model definition for node %s\n",attr_list[Nf-k+1],
                             MDAT->Node_name[i]);exit(1);}
            }
            // printf("pch = %s\n",pch);
            Gate[i]=new Boolean_Gate((unsigned int)Nf,0);
			
            for(kg=0;kg<Gate[i]->N;kg++){
                unsigned short int *anybase;
                anybase=new unsigned short int[Nf+1];
                decimal_to_any_base(kg,2,(unsigned int)Nf,anybase);
                for(k=1;k<=Nf;k++) {
                    fscanf(f,"%d ",&b);
                   // printf("%d ",anybase[Nf-k]);
                    if(anybase[Nf-k]!=b) {
                        printf("Boolean rules for gate %s are wrong: line %ld\n",MDAT->Node_name[i],kg);
                        exit(1);
                    }
                    //  printf("needs to be %d ; is %d ",kg/(int)pow(2,Nf-k), b);getchar();
                }
                delete []anybase; anybase=NULL;
				fscanf(f,"%d ",&b);
              //  printf("\n");
                //printf("Gate[%u] pozition %lu will become %d\n",i,kg,b);
                Gate[i]->O[kg]=b;
            }
           fclose(f);
           //printf("______________________________\n");getchar();
        }
        /*
        for(i=1;i<=N;i++){
            sc2=BN->node[i].Li.first_link_pt;	
            printf("i=%u: ",i);
             while(sc2!=NULL){
                printf("%lld ",sc2->sneighbor);
                sc2=sc2->next(sc2);//getchar();
            }
             printf("\n");
            //getchar();
         } 
        */    
		NoBowTie=0;
        export_Boolean_RegNet();
         MDAT->list_node_names();
       // getchar();
        active_inputs = new longint_ARRAY();
        active_inputs_noGF = new longint_ARRAY();
        silenced_inputs = new longint_ARRAY();
	}
	
	
	void export_Boolean_RegNet(){
		char fn[600];
		FILE *f;
		unsigned int i;
		unsigned long long k;
		snode::slinklist::slink *sc2;
		
		//if(!Evolved) 
			sprintf(fn,"%s/%s/%s.inf",MODEL_Directory,path,name);
//		else sprintf(fn,"%s/_NoBowTie_BRN/%s.inf",MODEL_Directory,path,name);
		//printf("About to overwrite file %s\n",fn);getchar();
        
        f=fopen(fn,"w");
        printf("%s\n",fn);// getchar();
        
		for(i=1;i<=N;i++){
			fprintf(f, "%d\t%lld\t%lld\t%lld\n",i,BN->node[i].ID,BN->node[i].k_in,BN->node[i].k_out);
          //  printf( "%d\t%lld\t%lld\t%lld\n",i,BN->node[i].ID,BN->node[i].k_in,BN->node[i].k_out);
            sc2=BN->node[i].Li.first_link_pt;
			while(sc2!=NULL){
				fprintf(f,"\t%lld",sc2->sneighbor);
         //       printf("\t%lld",sc2->sneighbor);
                sc2=sc2->next(sc2);//getchar();
			}
			fprintf(f,"\n\t");
           // printf("\n\t");
            for(k=0;k<Gate[i]->N;k++){
              //  printf("%d ",Gate[i]->O[k]);
				fprintf(f,"%d ",Gate[i]->O[k]);
            }
			fprintf(f,"\n");
          //  printf("\n");

		}
		fclose(f);
       // printf("Exported %s\n",fn);
       // getchar();
    }
	
	int import_Boolean_RegNet(){
		char fn[600];
		FILE *f;
		unsigned int i,k;
		unsigned long long incom,ko,j,*inc_temp;
		
		sprintf(fn,"%s/%s/%s.inf",MODEL_Directory,path,name);
        printf("Looking for %s\n",fn);//getchar();
        
		f=fopen(fn,"r");
		if(f==NULL) {return(0);}
		K=0;
        BN=new sgraph(N,name);
        Gate=new Boolean_Gate_list[N+1];
        for(i=1;i<=N;i++){
			fscanf(f, "%d\t%lld\t%lld\t%lld\n",&k,&(BN->node[i].ID),&j,&ko);
			//printf("Node %d, basin %lld, incoming: %lld\n",i,BN->node[i].ID,ko);//getchar();
			inc_temp=new unsigned long long[ko+1];
			for(k=1;k<=ko;k++){
				fscanf(f,"\t%lld",&incom);
				inc_temp[k]=incom;
				//printf("Link from %lld\n",incom);
			}
			for(k=1;k<=ko;k++){ BN->add_link(i,inc_temp[ko-k+1]);K++;}
            
            delete[] inc_temp;inc_temp=NULL;
			Gate[i]=new Boolean_Gate((unsigned int)ko,0);
			for(k=0;k<Gate[i]->N;k++)
				fscanf(f,"%hd ",&(Gate[i]->O[k]));
		}
		fclose(f);
		return(1);
	}
    void export_Boolean_RegNet_Module_ID(){
		char fn[600];
		FILE *f;
		unsigned int i;
		unsigned long long k;
		snode::slinklist::slink *sc2;
		
		//if(!Evolved)
        sprintf(fn,"%s/%s/%s_M-ID.inf",MODEL_Directory,path,name);
        //		else sprintf(fn,"%s/_NoBowTie_BRN/%s.inf",MODEL_Directory,path,name);
		f=fopen(fn,"w");
		for(i=1;i<=N;i++){
			fprintf(f, "%d\t%d\t%lld\t%lld\t%lld\n",i,MDAT->Module_ID[i],BN->node[i].ID,BN->node[i].k_in,BN->node[i].k_out);
			sc2=BN->node[i].Li.first_link_pt;
			while(sc2!=NULL){
				fprintf(f,"\t%lld",sc2->sneighbor);
				sc2=sc2->next(sc2);//getchar();
			}
			fprintf(f,"\n\t");
			for(k=0;k<Gate[i]->N;k++)
				fprintf(f,"%d ",Gate[i]->O[k]);
			fprintf(f,"\n");
			
		}
		fclose(f);
       // printf("Exported %s\n",fn);
       // getchar();
	}
    
    void export_Boolean_RegNet_David(){
        char fn[600];
        FILE *f;
        unsigned int i;
        unsigned long long k;
        snode::slinklist::slink *sc2;
        
        //if(!Evolved)
        sprintf(fn,"%s/%s/%s_David.inf",MODEL_Directory,path,name);
        //		else sprintf(fn,"%s/_NoBowTie_BRN/%s.inf",MODEL_Directory,path,name);
        f=fopen(fn,"w");
        for(i=1;i<=N;i++){
            sc2=BN->node[i].Li.first_link_pt;
            while(sc2!=NULL){
                fprintf(f,"%lld ",sc2->sneighbor);
                sc2=sc2->next(sc2);//getchar();
            }
            fprintf(f,"\n");
            if((MDAT!=NULL)&&(MDAT->Module_ID!=NULL))
                fprintf(f, "%d\n",MDAT->Module_ID[i]);
            else fprintf(f, "0\n");
            for(k=0;k<Gate[i]->N;k++)
                fprintf(f,"%d ",Gate[i]->O[k]);
            fprintf(f,"\n");
        }
        fclose(f);
    }
	
	int import_Boolean_RegNet_Module_ID(){
		char fn[600];
		FILE *f;
		unsigned int i,k;
		unsigned long long incom,ko,j,*inc_temp;
		
		if(!Evolved) sprintf(fn,"%s/%s/%s_M-ID.inf",MODEL_Directory,path,name);
        printf("Looking for %s\n",fn);//getchar();
        
		f=fopen(fn,"r");
		if(f==NULL) {return(0);}
		K=0;
        MDAT=new Boolean_RegNet_METADATA(N);
        MDAT->Module_ID=new unsigned int[N+1];
		for(i=1;i<=N;i++){
			fscanf(f, "%d\t%d\t%lld\t%lld\t%lld\n",&k,&(MDAT->Module_ID[i]),&(BN->node[i].ID),&j,&ko);
			printf("Node %d, basin %lld, incoming: %lld\n",i,BN->node[i].ID,ko);//getchar();
			inc_temp=new unsigned long long[ko+1];
			for(k=1;k<=ko;k++){
				fscanf(f,"\t%lld",&incom);
				inc_temp[k]=incom;
				printf("Link from %lld\n",incom);
			}
			for(k=1;k<=ko;k++){ BN->add_link(i,inc_temp[ko-k+1]);K++;}
            delete[] inc_temp;inc_temp=NULL;
			Gate[i]=new Boolean_Gate((unsigned int)ko,0);
			for(k=0;k<Gate[i]->N;k++)
				fscanf(f,"%hd ",&(Gate[i]->O[k]));
		}
		fclose(f);
		return(1);
	}
    
    unsigned int *import_Boolean_RegNet_Module_IDs_but_not_the_gates(int Module_randomizer,int rnd_now){
		char fn[600];
		FILE *f;
		unsigned int i,k;
        unsigned long long j,jj,ko;
        unsigned short int bla;
		unsigned int *MIDs;
		
        switch (Module_randomizer) {
            case 1:
                sprintf(fn,"%s/%s/%s_RMOD_MoveOne%d__Uncoupled_M-ID.inf",MODEL_Directory,path,name,rnd_now);
                break;
            case 2:
                sprintf(fn,"%s/%s/%s_RMOD_SwapOne%d__Uncoupled_M-ID.inf",MODEL_Directory,path,name,rnd_now);
                break;
            case 3:
                sprintf(fn,"%s/%s/%s_RMOD_RndPart%d__Uncoupled_M-ID.inf",MODEL_Directory,path,name,rnd_now);
                break;
                
            default:
                break;
        }
        printf("Looking for %s\n",fn);//getchar();
        
		f=fopen(fn,"r");
		if(f==NULL) {return(NULL);}
		MIDs=new unsigned int[N+1];
        MIDs[0]=0;
		for(i=1;i<=N;i++){
			fscanf(f, "%d\t%d\t%lld\t%lld\t%lld\n",&k,&(MIDs[i]),&j,&jj,&ko);
            if (MIDs[i]>MIDs[0]) MIDs[0]=MIDs[i];
           // printf("%d\t%d\t%lld\t%lld\t%lld\n",k,MIDs[i],j,jj,ko);
           // if (MIDs[i]>0)  printf ("\t\t\t%d - %d \n",i,MIDs[i]);

			for(k=1;k<=ko;k++){
				fscanf(f,"\t%lld",&j);
			}
            unsigned long long Ng=(unsigned long long)(pow(2.,(double)ko));
            for(k=0;k<Ng;k++){
				fscanf(f,"%hd ",&bla);
               // printf("\t\t%hd",bla);
            }
           // getchar();
            
		}
		fclose(f);
		return(MIDs);
	}

    double get_overall_gate_bias(){
        double p0now=0,n_all=0;
        for(int i=1;i<=N;i++){
            p0now+=Gate[i]->get_n0();
            n_all+=Gate[i]->N;
        }
        return(p0now/(double)n_all);
    }
    
	void export_Boolean_RegNet_with_names_to_Cytoskape(){
        unsigned long long i;
		unsigned int k,ok;
		char fn[600];
		FILE *f;
		snode::slinklist::slink *sc2;
        short int linksign;
        
        sprintf(fn,"mkdir %s/%s/Cytoskape\n",
				MODEL_Directory_system,path);
		system(fn);
        
        sprintf(fn,"%s/%s/Cytoskape/%s_Boolean_w_Names.txt",
				MODEL_Directory,path,name);
		f=fopen(fn,"w");
		fprintf(f,"start_Node\tregulates\tend_Node\n");
		for(i=1;i<=N;i++){
			sc2=BN->node[i].Li.first_link_pt;	
            k=0;
            while(sc2!=NULL){
                k++;
                ok=Gate[i]->test_function_for_input(k);
                if(ok) {
                    linksign=Gate[i]->check_overall_sign_of_input(k);
                   // printf("%lld: ",i);
                   // printf("%s <-",MDAT->Node_name[i]);
                  //  printf("%lld: ",sc2->sneighbor);
                   //if(sc2->sneighbor>0) printf("%s sign=", MDAT->Node_name[sc2->sneighbor]);
                  // else getchar();
                  //  printf("%d\n",linksign);
                    switch (linksign) {
                        case 1:
                            fprintf(f,"%s\t<--\t%s\n",MDAT->Node_name[i],MDAT->Node_name[sc2->sneighbor]);
                            break;
                        case -1:
                            fprintf(f,"%s\t|--\t%s\n",MDAT->Node_name[i],MDAT->Node_name[sc2->sneighbor]);
                            break;
                        case 0:
                            fprintf(f,"%s\t*--\t%s\n",MDAT->Node_name[i],MDAT->Node_name[sc2->sneighbor]);
                            break;
                        default:
                            fprintf(f,"%s\t*--\t%s\n",MDAT->Node_name[i],MDAT->Node_name[sc2->sneighbor]);
                            break;
                    }
                }
				sc2=sc2->next(sc2);
			}
		}
		fclose(f);
    }	
    
	void generate_Rnd_Boolean_RegNet(unsigned int rs){
		unsigned long long i,a,b;
        if(K<=N) {printf("We need at least as many total links as the number of nodes!\n");exit(1);}

		init_random(rs);
		NoBowTie=0;
		BN=NULL;
		while (NoBowTie==0) {
			//if(BN==NULL) {delete BN;BN=NULL;}
			BN=new sgraph(N,name);
			for(i=1;i<=K;i++){
				a=rand_int(1,N);
				b=rand_int(1,N);
				if(!(BN->add_link(a,b))) i--;
			}
			NoBowTie=No_BowTie_Evolution();
		}
	}
	
	int No_BowTie_Evolution(){
		unsigned int itnow,i,ok;
		unsigned long long a,link_rev;
		longlongint_ARRAY_pair *sparelinks;
		longlongint_ARRAY *leafnodes,*inputnodes;
		snode::slinklist::slink *sc2;
		
		itnow=0;
		while(itnow<=NOBOWTIE_IT_MAX){
			sparelinks=new longlongint_ARRAY_pair();
			leafnodes=new longlongint_ARRAY();
			inputnodes=new longlongint_ARRAY();
			for(i=1;i<=N;i++){
				if((!((BN->node[i].k_out==1)&&(BN->node[i].Li.first_link_pt->sneighbor!=i)))&&
				   (!((BN->node[i].k_out==2)&&(BN->node[i].Li.get_dirlink(i))))){ // more than 1 outlink
					sc2=BN->node[i].Li.first_link_pt;	
					while(sc2!=NULL){
						if(!(BN->node[sc2->sneighbor].k_in==1)) {
							sparelinks->add_element(i,sc2->sneighbor);
							//printf("\t\t%d -> %lld marked as spare!\n",i,sc2->sneighbor);
						}
						sc2=sc2->next(sc2);//getchar();
					}
				}
				if((BN->node[i].k_in==0)||((BN->node[i].k_in==1)&&(BN->node[i].Li.first_link_pt!=NULL)&&(BN->node[i].Li.first_link_pt->sneighbor==i))) 
					leafnodes->add_element(i);
				if((BN->node[i].k_out==0)||((BN->node[i].k_out==1)&&(BN->node[i].Li.first_link_pt->sneighbor==i))) 
					inputnodes->add_element(i);
			}
			if(sparelinks->N<leafnodes->N+inputnodes->N) {printf("Not enough spare links!!!\n");return(0);}
			if((leafnodes->N>0)||(inputnodes->N>0)) itnow++;
				else itnow=NOBOWTIE_IT_MAX+2;
			
		//	printf("%lld spare links, %lld leaf nodes, %lld input nodes\n",sparelinks->N,leafnodes->N,inputnodes->N);
			//sparelinks->print_all();
			for(i=1;i<=leafnodes->N;i++){
				link_rev=rand_int(1,sparelinks->N);
				//printf("\t\tlink_rev_index=%lld\n",link_rev);
				BN->delete_link(sparelinks->A[link_rev],sparelinks->B[link_rev]);
				do{ ok=1;
					a=rand_int(1,N);
					if(!BN->add_link(a,leafnodes->A[i])) ok=0; 
					//else printf("Rewiring %lld -> %lld to %lld -> %lld\n",sparelinks->A[link_rev],sparelinks->B[link_rev],a,leafnodes->A[i]);
				}
				while(ok==0);
				sparelinks->delete_element(link_rev);
				//sparelinks->print_all();
				//printf("\n");
			}
			for(i=1;i<=inputnodes->N;i++){
				link_rev=rand_int(1,sparelinks->N);
				//printf("\t\tlink_rev_index=%lld\n",link_rev);
				BN->delete_link(sparelinks->A[link_rev],sparelinks->B[link_rev]);
				do{ ok=1;
					a=rand_int(1,N);
					if(!BN->add_link(inputnodes->A[i],a)) ok=0; 
					//else printf("Rewiring %lld -> %lld to %lld -> %lld\n",sparelinks->A[link_rev],sparelinks->B[link_rev],inputnodes->A[i],a);
				}
				while(ok==0);
				sparelinks->delete_element(link_rev);
				//sparelinks->print_all();
				//printf("\n");
			}		
			//getchar();
			delete sparelinks;sparelinks=NULL;
			delete leafnodes;leafnodes=NULL;
			delete inputnodes;inputnodes=NULL;
		}
		ok=0;
		while(ok==0){
			ok=test_modularity_in_burning();
			//printf("%d\t Modularity_test_failed\n\n\n",ok);
		}
		if(ok==1) return(1);
		else return(0);
	}
	
    longint_ARRAY *Find_Strongly_Connected_Component(){
        unsignedintpointer *Reachable;
        unsigned int i,j,k,ok;
        longint_ARRAY *SCC;
        
        SCC=new longint_ARRAY();
        Reachable=new unsignedintpointer[N+1];
        for(i=0;i<=N;i++){
            Reachable[i]=new unsigned int[N+1];
            for(j=0;j<=N;j++) Reachable[i][j]=0;
            Reachable[i][i]=1;
        }
        ok=0;
        while(ok==0) {
            ok=1;
            for(i=1;i<=N;i++){
                for(j=1;j<=N;j++)
                    if(Reachable[i][j]==0){
                        for(k=1;k<=N;k++)
                            if((Reachable[i][k]==1)&&(BN->node[k].Li.get_dirlink(j)))
                            {Reachable[i][j]=1; ok=0;}
                    }
            }
        }
        for(i=2;i<=N;i++)
            for(j=1;j<i;j++)
                if (Reachable[i][j]==Reachable[j][i]==1) {
                    SCC->add_element(i);
                    SCC->add_element(j);
                }
        return(SCC);
    }
    
	int test_modularity_in_burning(){
		unsignedintpointer *Reachable;
		unsigned int i,j,k,ok;
		longlongint_ARRAY_pair *sparelinks;
		snode::slinklist::slink *sc2;
		unsigned long long a,b,link_rev;
		
		Reachable=new unsignedintpointer[N+1];
		for(i=0;i<=N;i++){
			Reachable[i]=new unsigned int[N+1];
			for(j=0;j<=N;j++) Reachable[i][j]=0;
			Reachable[i][i]=1;
		}
		ok=0;
		while(ok==0) {
			ok=1;
			for(i=1;i<=N;i++){
				for(j=1;j<=N;j++)
					if(Reachable[i][j]==0){
						for(k=1;k<=N;k++)
							if((Reachable[i][k]==1)&&(BN->node[k].Li.get_dirlink(j))) 
								{Reachable[i][j]=1; ok=0;}
					}
			}
		}
		//for(i=1;i<=N;i++)
		//	for(j=1;j<=N;j++) printf("R[%d][%d]=%d\n",i,j,Reachable[i][j]);
		ok=1;
		for(i=1;i<=N;i++)
			for(j=1;j<=N;j++) if(Reachable[i][j]==0) ok=0;
		if(ok) {
			for(i=0;i<=N;i++) {delete[] Reachable[i];Reachable[i]=NULL;}
			delete[] Reachable;Reachable=NULL;
			return(1);
		}
		
		sparelinks=new longlongint_ARRAY_pair();
		for(i=1;i<=N;i++){
			if((!((BN->node[i].k_out==1)&&(BN->node[i].Li.first_link_pt->sneighbor!=i)))&&
			   (!((BN->node[i].k_out==2)&&(BN->node[i].Li.get_dirlink(i))))){ // more than 1 outlink
				sc2=BN->node[i].Li.first_link_pt;	
				while(sc2!=NULL){
					if(!(BN->node[sc2->sneighbor].k_in==1)) {
						sparelinks->add_element(i,sc2->sneighbor);
						//printf("\t\t%d -> %lld marked as spare!\n",i,sc2->sneighbor);
					}
					sc2=sc2->next(sc2);//getchar();
				}
			}
		}
		if(sparelinks->N<1) { //printf("No more spare links!!!\n");
			return(-1);
		}
		
		link_rev=rand_int(1,sparelinks->N);
		//printf("\t\tlink_rev_index=%lld\n",link_rev);
		BN->delete_link(sparelinks->A[link_rev],sparelinks->B[link_rev]);
		do{ ok=1;
			do{
				a=rand_int(1,N);
				b=rand_int(1,N);
			}
			while(Reachable[a][b]==1);
			if(!BN->add_link(a,b)) ok=0; 
			//else printf("Rewiring %lld -> %lld to %lld -> %lld\n",sparelinks->A[link_rev],sparelinks->B[link_rev],a,b);
		}
		while(ok==0);
		sparelinks->delete_element(link_rev);
		//sparelinks->print_all();
		//printf("\n");
		delete sparelinks;sparelinks=NULL;
		
		for(i=0;i<=N;i++) {delete[] Reachable[i];Reachable[i]=NULL;}
		delete[] Reachable;Reachable=NULL;
		return(0);
	}
	
    void export_two_way_reachability_matrix(){
		unsignedintpointer *Reachable;
		unsigned int i,j,k,ok;
		char fn[600];
        FILE *f;
        
		Reachable=new unsignedintpointer[N+1];
		for(i=0;i<=N;i++){
			Reachable[i]=new unsigned int[N+1];
			for(j=0;j<=N;j++) Reachable[i][j]=0;
			Reachable[i][i]=1;
		}
		ok=0;
		while(ok==0) {
			ok=1;
			for(i=1;i<=N;i++){
				for(j=1;j<=N;j++)
					if(Reachable[i][j]==0){
						for(k=1;k<=N;k++)
							if((Reachable[i][k]==1)&&(BN->node[k].Li.get_dirlink(j))) 
                            {Reachable[i][j]=1; ok=0;}
					}
			}
		}
		//for(i=1;i<=N;i++)
		//	for(j=1;j<=N;j++) printf("R[%d][%d]=%d\n",i,j,Reachable[i][j]);
		 
        sprintf(fn,"%s/%s/%s_Two_way_paths.dat",
				MODEL_Directory,path,name);
		f=fopen(fn,"w");
        if(MDAT!=NULL) {
            fprintf(f,"NM/NM");
            for(i=1;i<=N;i++) fprintf(f,"\t%s",MDAT->Node_name[i]);
            fprintf(f,"\n");
        }
        for(i=1;i<=N;i++){
            if(MDAT!=NULL) fprintf(f,"%s\t",MDAT->Node_name[i]);
			for(j=1;j<=N;j++) fprintf(f,"%d\t",Reachable[i][j]);
            fprintf(f,"\n");
        }
		fclose(f);
        if(MDAT!=NULL){
            sprintf(fn,"%s/%s/%s_Two_way_links.dat",
				MODEL_Directory,path,name);
            f=fopen(fn,"w");
            for(i=1;i<=N;i++){
               for(j=i;j<=N;j++) if((Reachable[i][j]==1)&&(Reachable[j][i]==1)) 
                   fprintf(f,"%s - %s\n",MDAT->Node_name[i],MDAT->Node_name[j]);
            }
            fclose(f);
        }
	}
    
    int load_Module_IDs_and_Drawing_Order(){
        char fn[400];
        char attr_list[600][5000];
        int ok = 0;
        FILE *f;
        
        strcpy(MDAT->MOD_Assignment_name,Model_type);
        Module_Order=new longint_ARRAY_pair();
        MDAT->Module_NR=0;
        Module_Names.push_back("");
        MDAT->Module_ID=new unsigned int[N+1];
        for(int i=0; i<=N; i++) MDAT->Module_ID[i]=0;
       
        sprintf(fn,"%s/%s/%s_ModelMapping.txt",MODEL_Directory,path,name);
        f=fopen(fn,"r");
        if(f==NULL) {
            printf("Module membership of nodes and their drawing order not found\n Please create file %s\t",fn); getchar();
            exit(1);
        }
        
        while (fgets(GSE_line,MAX_LINE,f)!=NULL) {
            char bla1[600], bla3[600];
            sscanf(GSE_line, "%s %s", bla1,bla3);
           // printf("w1 = %s\tw2=%s\n",bla1,bla3); getchar();
            std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(GSE_line, '(',')');
            if (word.size() > 1) {
                printf("Line %s not correctly formatted; codes needs a single list of node-names enclosed in ( ... ) brackets with no space on either side, before / after the first /last nodename. Nodes sgould be separatd by comma and space: e.g.: (ORC, Cdc6, Cdt1, Pre_RC, geminin) \n",GSE_line);
                getchar(); exit (1);
            }
            char bla2[600]; sprintf(bla2,"%s",word[0].c_str());
            unsigned long int Nf=separate_line_by_string(bla2,attr_list,", ",200);
            if (Nf > 0) {
                MDAT->Module_NR++;
                Module_Names.push_back(bla3);
            }
            else {printf("Line %s not correctly formatted; codes needs a single list of node-names enclosed in ( ... ) brackets with no space on either side, before / after the first /last nodename. Nodes sgould be separatd by comma and space: e.g.: (ORC, Cdc6, Cdt1, Pre_RC, geminin) \n",GSE_line);
                getchar(); exit (1);
            }
        }
        fclose(f);
        
        Draw_in_Modules=new longint_ARRAY_list[MDAT->Module_NR+1];
        for (int i=0; i<=MDAT->Module_NR; i++) Draw_in_Modules[i]=new longint_ARRAY();
        unsigned int mnr=0;
        f=fopen(fn,"r");
        while (fgets(GSE_line,MAX_LINE,f)!=NULL) {
            std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(GSE_line, '(',')');
            char bla2[600]; sprintf(bla2,"%s",word[0].c_str());
            unsigned long int Nf=separate_line_by_string(bla2,attr_list,", ",200);
            if (Nf > 0){
                mnr++;
                Module_Order->add_element(mnr,0);
                for (unsigned int i=1; i<=Nf; i++) {
                    unsigned long int nid=MDAT->get_node_ID(attr_list[i]);
                    if(nid>0){
                        MDAT->Module_ID[nid]=mnr;
                        Draw_in_Modules[mnr]->add_element(nid);
                    }
                    else { printf("Node name %s from list of module nodes is not in the Boolean Model!\n",attr_list[i]);
                        getchar(); exit (1);
                    }
                }
            }
        }
        fclose(f);

        for(int i=1; i<=N; i++)
            if (MDAT->Module_ID[i]==0) {
                printf("Boolean node %s is missing from all modules and will NOT be shown on time course!\n",MDAT->Node_name[i]);
                getchar();
            }
        return ok;
    }
    
    int load_Module_Drawing_order(const char mod_ass_name[]){
        char fn[400],bla[100];
        int nid,mid,o,m;
        int vege;
        FILE *f;
        
        m=0;
    
        Module_Order=new longint_ARRAY_pair();
        sprintf(fn,"%s/%s/%s_%s_ORDER.txt",MODEL_Directory,path,name,mod_ass_name);
        f=fopen(fn,"r");
        if(f==NULL) {
            printf("test! fn=%s\n",fn);getchar();
            printf("Module order for drawing %s not found\n%s\t",mod_ass_name,fn);
            exit(1);}
        
        vege=fscanf(f,"%s %s %s",bla,bla,bla);
        while (vege!=EOF) {
            vege=fscanf(f,"%d %s %d",&o,bla,&mid);
            nid=MDAT->get_node_ID(bla);
            if((nid>0)&&(mid>0)){
                printf("Please use module names that DO NOT match a particular node name!\n");
                exit(1);
            }
            if(nid>0){
                Module_Order->add_element(0,nid);
                // printf("%d %s in module %d\n",nid, MDAT->Node_name[nid],MDAT->Module_ID[nid]);
            }
            else {Module_Order->add_element(mid,0);}
            if(m<mid) m=mid;
        }
        fclose(f);
        return(m);
    }
    
    void load_Module_IDs(const char mod_ass_name[],int Mod_NR){
        char fn[600],bla[100];
        int nid,mid;
        int vege;
        FILE *f;
        
        Draw_in_Modules=new longint_ARRAY_list[Mod_NR+1];
        for (int i=0; i<=Mod_NR; i++) {
            Draw_in_Modules[i]=new longint_ARRAY();
        }
        
        sprintf(fn,"%s/%s/%s_%s.MOD",MODEL_Directory,path,name,mod_ass_name);
        f=fopen(fn,"r");
        if(f==NULL) {printf("Module assignment %s not found\n%s\t",mod_ass_name,fn);exit(1);}
        MDAT->Module_ID=new unsigned int[N+1];
        for(int i=0; i<=N; i++) MDAT->Module_ID[i]=0;
        strcpy(MDAT->MOD_Assignment_name,mod_ass_name);
        vege=1;
        while (vege!=EOF) {
            vege=fscanf(f,"%s %d",bla,&mid);
            nid=MDAT->get_node_ID(bla);
            if(nid>0){  MDAT->Module_ID[nid]=mid;
                        Draw_in_Modules[mid]->add_element(nid);
                        if(MDAT->Module_NR<mid) MDAT->Module_NR=mid;
                       // printf("%d %s in module %d\n",nid, MDAT->Node_name[nid],MDAT->Module_ID[nid]);
            }
        }
        fclose(f);
    }
    
    void load_Module_IDs_old(const char mod_ass_name[]){
        char fn[600],bla[100];
        int nid,mid;
        int vege;
        FILE *f;
        
    
        sprintf(fn,"%s/%s/%s_%s.MOD",MODEL_Directory,path,name,mod_ass_name);
        f=fopen(fn,"r");
        if(f==NULL) {printf("Module assignment %s not found\n%s\t",mod_ass_name,fn);exit(1);}
        MDAT->Module_ID=new unsigned int[N+1];
        strcpy(MDAT->MOD_Assignment_name,mod_ass_name);
        vege=1;
        while (vege!=EOF) {
            vege=fscanf(f,"%s %d",bla,&mid);
            nid=MDAT->get_node_ID(bla);
            if(nid>0){  MDAT->Module_ID[nid]=mid;
                if(MDAT->Module_NR<mid) MDAT->Module_NR=mid;
                // printf("%d %s in module %d\n",nid, MDAT->Node_name[nid],MDAT->Module_ID[nid]);
            }
        }
        fclose(f);
    }
    
    void print_gate(unsigned long int ind1){
        snode::slinklist::slink *sc2;
        unsigned int *sin_now;
        
        printf("Gate of node %s\n",MDAT->Node_name[ind1]);
        
        sc2=BN->node[ind1].Li.first_link_pt;
        while(sc2!=NULL){
            printf("%s\t\t",MDAT->Node_name[sc2->sneighbor]);
            sc2=sc2->next(sc2);
        }
        printf(">> %s\t\t",MDAT->Node_name[ind1]);
        
        printf("\n");
        sin_now=new unsigned int[N+1];
        for(unsigned long long k=0;k<Gate[ind1]->N;k++){
            decimal_to_any_base(k,2,Gate[ind1]->K,sin_now);
            for(unsigned long long l=1;l<=Gate[ind1]->K;l++)
                printf("%d\t\t",sin_now[Gate[ind1]->K-l]);
            printf("\t\t | %d\n",Gate[ind1]->O[k]);
        }
        delete[] sin_now;sin_now=NULL;
    }
    
    //void print_gates_to_File(const char fn[]){
    void print_gates_to_File(){
        FILE *f;
        snode::slinklist::slink *sc2;
        unsigned int *sin_now;
        char fn[400];
        sprintf(fn,"%s/%s/%s_GATES.txt",MODEL_Directory,path,name);
        f=fopen(fn,"w");
        for(long int ind1=1;ind1<=N;ind1++){
            fprintf(f,"Gate of node %s\n",MDAT->Node_name[ind1]);
            sc2=BN->node[ind1].Li.first_link_pt;
            while(sc2!=NULL){
                fprintf(f,"%s\t\t",MDAT->Node_name[sc2->sneighbor]);
                sc2=sc2->next(sc2);
            }
            fprintf(f,">> %s\t\t",MDAT->Node_name[ind1]);
            
            fprintf(f,"\n");
            sin_now=new unsigned int[N+1];
            for(unsigned long long k=0;k<Gate[ind1]->N;k++){
                decimal_to_any_base(k,2,Gate[ind1]->K,sin_now);
                for(unsigned long long l=1;l<=Gate[ind1]->K;l++)
                    fprintf(f,"%d\t\t",sin_now[Gate[ind1]->K-l]);
                fprintf(f,"\t\t | %d\n",Gate[ind1]->O[k]);
            }
            delete[] sin_now;sin_now=NULL;
        }
        fclose(f);

    }

    //sprintf(fn, "%s/%s/%s_GATES.txt",MODEL_Directory,PD_A->Coupled->DIRname,PD_A->Coupled->SAMPLEname);
    //PD_A->Coupled->BRN->print_gates_to_File(fn);
    //sprintf(fn, "%s/%s",MODEL_Directory,PD_A->Coupled->DIRname);
    
    void print_gates_to_tt_folder(const char fn[]){
        FILE *f;
        snode::slinklist::slink *sc2;
        unsigned int *sin_now;
        char fn2[600];
        
        for(long int ind1=1;ind1<=N;ind1++){
            sprintf(fn2,"%s/%s/%s_Fine/%s.csv", fn,fn,fn,MDAT->Node_name[ind1]);
            f=fopen(fn2,"w");
            //fprintf(f,"Gate of node %s\n",MDAT->Node_name[ind1]);
            sc2=BN->node[ind1].Li.first_link_pt;
            while(sc2!=NULL){
                fprintf(f,"%s ",MDAT->Node_name[sc2->sneighbor]);
                sc2=sc2->next(sc2);
            }
            fprintf(f,"%s\n",MDAT->Node_name[ind1]);
            sin_now=new unsigned int[N+1];
            for(unsigned long long k=0;k<Gate[ind1]->N;k++){
                decimal_to_any_base(k,2,Gate[ind1]->K,sin_now);
                for(unsigned long long l=1;l<=Gate[ind1]->K;l++)
                    fprintf(f,"%d\t",sin_now[Gate[ind1]->K-l]);
                fprintf(f,"%d\n",Gate[ind1]->O[k]);
            }
            delete[] sin_now;sin_now=NULL;
            fclose(f);
        }
    }
    
    void print_PARTIAL_gate(unsigned long int ind1,longint_ARRAY *gl,Boolean_Gate *G_p){
        unsigned int *sin_now;
        
//        printf("Partial Gate of node %s\n",MDAT->Node_name[ind1]);
//        
//        for(long int i=1;i<=gl->N;i++)
//            printf("%s\t",MDAT->Node_name[gl->A[i]]);
//        printf("%s\t",MDAT->Node_name[ind1]);
//        printf("\n");
        sin_now=new unsigned int[N+1];
        for(unsigned long long k=0;k<G_p->N;k++){
            decimal_to_any_base(k,2,G_p->K,sin_now);
            //for(unsigned long long l=1;l<=G_p->K;l++)
            //    printf("%d\t\t",sin_now[G_p->K-l]);
            //printf("\t\t | %d\n",G_p->O[k]);
        }
    }
    
    Boolean_Gate *get_PARTIAL_gate(unsigned long int ind1, int canalyze[],longint_ARRAY *gl){
        Boolean_Gate *G_test;
        unsigned int *s_n,*sin_full,*full_poz,*cut_poz;
        snode::slinklist::slink *sc2;
        
        full_poz=new unsigned int[N+1];
        cut_poz=new unsigned int[N+1];
        for(int i=1;i<=N;i++) {full_poz[i]=0;cut_poz[i]=0;}
        
        int pf=0;int pc=0;
        sc2=BN->node[ind1].Li.first_link_pt;
        while(sc2!=NULL){
            pf++; full_poz[sc2->sneighbor]=pf;
            //printf("Full poz %s =%d\n",MDAT->Node_name[sc2->sneighbor],full_poz[sc2->sneighbor]);
            if(canalyze[sc2->sneighbor]==-2) {
                gl->add_element(sc2->sneighbor);
                pc++; cut_poz[sc2->sneighbor]=pc;
              //  printf("\t Cut poz %s =%d\n",MDAT->Node_name[sc2->sneighbor],cut_poz[sc2->sneighbor]);
            }
            sc2=sc2->next(sc2);
        }
        
        
        if(gl->N>0){
            G_test=new Boolean_Gate((unsigned int)gl->N);
        
            s_n=new unsigned int[G_test->K+1];
            sin_full=new unsigned int[Gate[ind1]->K+1];
            
            for(unsigned long long k=0;k<G_test->N;k++){
                decimal_to_any_base(k,2,G_test->K,s_n);
                
                //if(prnow){ printf("Figuring out index %lld, string: ",k);
                //    for(unsigned long long l=1;l<=G_test->K;l++)
                //        printf("%d\t\t",s_n[G_test->K-l]);
                //    printf("\n");
                //}
                
                for(int i=1;i<=N;i++){
                    if(canalyze[i]!=-1) { // i is an input of j
                        if(canalyze[i]==-2) { // free input
                            //sin_full[Gate[ind1]->K-full_poz[i]+1]=s_n[cut_poz[i]-1];
                            sin_full[full_poz[i]]=s_n[G_test->K-cut_poz[i]];
                        }
                        else{    // canalyzed input
                            //sin_full[Gate[ind1]->K-full_poz[i]+1]=canalyze[i];
                            sin_full[full_poz[i]]=canalyze[i];
                        }
                    }
                }
                
                //if(prnow){
                //    printf("Turned into full string: ");
                //    for(unsigned long long l=1;l<=Gate[ind1]->K;l++)
                //        printf("%d\t\t",sin_full[l]);
                //    printf("\n");
                //}
                unsigned long long int k_0=Gate[ind1]->gate_value_INDEX_for_input(sin_full);
               // if(prnow) printf("Index in full gate: %lld\n",k_0);
                G_test->O[k]=Gate[ind1]->O[k_0];
               // if(prnow) getchar();
            }
            delete[] s_n;s_n=NULL;
            delete[] sin_full;sin_full=NULL;
        }
        else {
            G_test=new Boolean_Gate(1);
            G_test->O[0]=0;G_test->O[1]=0;
        }
        
        delete[] full_poz;full_poz=NULL;
        delete[] cut_poz;cut_poz=NULL;
        return(G_test);
    }
    
    Boolean_Gate *get_FULL_gate(unsigned long int ind1, longint_ARRAY *gl){
        Boolean_Gate *G_test;
        snode::slinklist::slink *sc2;
        
        sc2=BN->node[ind1].Li.first_link_pt;
        while(sc2!=NULL){
            gl->add_element(sc2->sneighbor);
            sc2=sc2->next(sc2);
        }
    
        if(gl->N>0){
            G_test=new Boolean_Gate((unsigned int)gl->N);
            for(unsigned long long k=0;k<G_test->N;k++){
                G_test->O[k]=Gate[ind1]->O[k];
                // if(prnow) getchar();
            }
        }
        else {
            G_test=new Boolean_Gate(1);
            G_test->O[0]=0;G_test->O[1]=0;
        }
        return(G_test);
    }
    
    void test_gate(unsigned long int ind1, int canalyze[]){
        Boolean_Gate *G_test;
        longint_ARRAY *gl;
        
       // print_gate(ind1);
        gl=new longint_ARRAY();
        G_test=get_PARTIAL_gate(ind1,canalyze,gl);
      //  print_PARTIAL_gate(ind1,gl,G_test);
    }
    
    void identify_input_nodes(){
        longint_ARRAY *badin;
        snode::slinklist::slink *sc2;
       
        active_inputs = new longint_ARRAY();
        active_inputs_noGF = new longint_ARRAY();
        

        for (unsigned int i=1; i<=N; i++) {
            badin=Gate[i]->test_function_for_all_inputs_return_bad_ones();
            if (badin->N>0) {
                if ((badin->N==1) && (BN->node[i].Li.first_link_pt->sneighbor == i)) {
                      // silenced input!
                    silenced_inputs->add_element(i);
                }
                else {
                    printf("The following input node(s) do not affect the outcome of the node %s Gate; please correct:\n",MDAT->Node_name[i]);
                    sc2=BN->node[i].Li.first_link_pt;
                    unsigned int k=1;
                    while(sc2!=NULL){
                        if (badin->check_for_element(k))
                            printf("\t Input %s\n",MDAT->Node_name[sc2->sneighbor]);
                        sc2=sc2->next(sc2);
                        k++;
                    }
                }
            }
            else {
                sc2=BN->node[i].Li.first_link_pt;
                if((sc2->next(sc2)==NULL) && (sc2->sneighbor == i)) {
                    active_inputs->add_element(i);
                    if (strcmp( MDAT->Node_name[i], GF_high_name) != 0 ){
                        active_inputs_noGF -> add_element(i);
                      //  printf("--- %s\n",MDAT->Node_name[i]);
                    }
                }
                if (strcmp( MDAT->Node_name[i], GF_low_name) == 0 )
                    active_inputs->add_element(i);
                if (strcmp( MDAT->Node_name[i], CD_low_name) == 0 )
                    active_inputs->add_element(i);
            }
            delete badin; badin = NULL;
        }
    }
    
    void input_node_report(){
        for(unsigned int i=1;i<=active_inputs->N;i++){
            printf("Active input node: %s\t",MDAT->Node_name[active_inputs->A[i]]);
            if (active_inputs_noGF->check_for_element(i)==0) printf("GF input!\n");
            else printf("\n");
        }
        printf("\n");
        for(unsigned int i=1;i<=silenced_inputs->N;i++){
            printf("Silenced input: %s\n",MDAT->Node_name[silenced_inputs->A[i]]);
        }
       // getchar();
    }
    
	~Boolean_RegNet(){
		delete BN;BN=NULL;
        for(int i=1;i<=N;i++) {delete Gate[i];Gate[i]=NULL;}
        delete[] Gate;Gate=NULL;
        //if(MDAT!=NULL) {delete MDAT; MDAT=NULL;}
	}
};

typedef Boolean_RegNet * Boolean_RegNet_array;

class Boolean_Dynamics{
public:	
	Boolean_RegNet *BRN;
	unsigned short int *s;
	longlongint_ARRAY_list *Attractor_Cycle;
	sgraph *State_G,*Attractor_G;
	unsigned long int Attractor_NR;
	
	double beta,error_p,TM_sync,TM_firsterror;
	doublepointer *MARKOV_TM,*P_EQ_trans;
	double *P_EQ_state,*Dp;
	
	Boolean_Dynamics(Boolean_RegNet *b){
		unsigned int i;
		BRN=b;
		s=new unsigned short int[BRN->N+1];
		for(i=1;i<=BRN->N;i++) s[i]=0;
		State_G=NULL;
		Attractor_NR=0;
		Attractor_Cycle=NULL;
	}
	
	int update_Gate(unsigned long long node_ID){
		snode::slinklist::slink *sc2;
		unsigned long long k_1;
		unsigned int k;
		k_1=0;
		sc2=BRN->BN->node[node_ID].Li.first_link_pt;	
		// if(BRN->N>10) printf("Updating node %lld  inputs:",node_ID);
		k=0;
		while(sc2!=NULL){
			k++;
			if(s[sc2->sneighbor]==1) k_1+=(unsigned long long)(pow(2.,(double)BRN->BN->node[node_ID].k_out-k));
			// if(BRN->N>10) printf("%lld:%d ",sc2->sneighbor,s[sc2->sneighbor]);
			sc2=sc2->next(sc2);//getchar();
		}
		// if(BRN->N>10) printf("k1=%lld\t returning %d\n",k_1,BRN->Gate[node_ID]->O[k_1]);
		if(k>0) return(BRN->Gate[node_ID]->O[k_1]);
		else return(s[node_ID]);
	}
    
    
	void calculate_noisy_landscape(double bet){
		unsigned long long Nst,i,j;
		unsigned int *state_i,*state_j,*s_n,k,g,Gi,ne;
		//gsl_permutation * p;
		double prob,pozmin;
		gsl_matrix *A;
		//gsl_vector *x;//,*b;
		
		beta=bet;
		error_p=1./(1.+exp(2*beta));
		Nst=(unsigned long long)pow(2.,(double)BRN->N);
		//printf("Nst=%lld\n",Nst);
		
		s_n=new unsigned int[BRN->N+1];
		state_i=new unsigned int[BRN->N+1];
		state_j=new unsigned int[BRN->N+1];
		
		MARKOV_TM=new doublepointer[Nst+1];
		for(i=0;i<=Nst;i++) {
			MARKOV_TM[i]=new double[Nst+1];
			for(j=0;j<=Nst;j++) MARKOV_TM[i][j]=0;
		}
		
		for(i=1;i<=Nst;i++){ // on all possible system states
			decimal_to_any_base(i-1,2,BRN->N,s_n);//s_n takes the base 2 code of state i-1
			for(k=0;k<BRN->N;k++)
				state_i[BRN->N-k]=s_n[k]; // state_i is set as starting state
			set_state(i);
			for(j=1;j<=Nst;j++){ // all possible states j to transition to
				decimal_to_any_base(j-1,2,BRN->N,s_n); // s_n encodes these states
				for(k=0;k<BRN->N;k++)
					state_j[BRN->N-k]=s_n[k]; // state_j is set as target state
				prob=1.;ne=0;
				for(g=1;g<=BRN->N;g++){ // all nodes need to go from i to j -> prob of this is prob.
					Gi=update_Gate(g); // Gi=what node g would become deterministically
					if(Gi==state_j[g]) prob=prob*(1.-error_p); // if that is the state in j, non-error prob
					else {prob=prob*error_p;ne++;} // if not, error only prob
				}
				MARKOV_TM[i][j]=prob;
				if(ne==0) TM_sync=prob;
				if(ne==1) TM_firsterror=prob;
				//printf("%lld %lld Markov_p=%lg\n",i,j,MARKOV_TM[i][j]);
			}
			//MARKOV_TM[i][i]-=1.;
			//getchar();
		}
		
		A	= gsl_matrix_calloc (Nst+1, Nst);
		for(i=1;i<=Nst;i++){
			for(j=1;j<=Nst;j++) gsl_matrix_set(A, j-1, i-1, MARKOV_TM[i][j]);
			gsl_matrix_set(A, i-1, i-1, MARKOV_TM[i][i]-1.);
		}
		for(i=1;i<=Nst;i++) gsl_matrix_set(A, Nst, i-1, 1.);
		
		gsl_vector *b = gsl_vector_alloc (Nst+1);
		
		for(i=0;i<Nst;i++) gsl_vector_set(b,i,0);
		gsl_vector_set(b,Nst,1);
		
		gsl_vector *x = gsl_vector_calloc (Nst);
		
		gsl_vector *A_tau = gsl_vector_calloc(Nst);
		gsl_vector *residual = gsl_vector_calloc(Nst+1);
		
		//gsl_matrix_fprintf (stdout, A, "%.8lf");
		//printf("\n\n");
		gsl_linalg_QR_decomp(A, A_tau);
				
		gsl_linalg_QR_lssolve(A, A_tau, b, x,residual);
		gsl_vector_free(A_tau);
		
		//printf ("x = \n");
		//gsl_vector_fprintf (stdout, x, "%g");
		//printf ("\n\nresidual = \n");
		//gsl_vector_fprintf (stdout, residual, "%g");
		
		//getchar();
		P_EQ_state	= new double [Nst+1];
		pozmin=1;
		for(i=1;i<=Nst;i++) {
			P_EQ_state[i]=gsl_vector_get(x,i-1);
			if((P_EQ_state[i]>0) && (P_EQ_state[i]<pozmin)) pozmin=P_EQ_state[i];
			//printf("P_EQ_state[%lld]=%lg\t\t%lg\n",i,P_EQ_state[i],gsl_vector_get(x,i-1));
		}
		for(i=1;i<=Nst;i++) if(P_EQ_state[i]<=0) P_EQ_state[i]=pozmin/10.;
		//getchar();
		P_EQ_trans=new doublepointer[Nst+1];
		for(i=0;i<=Nst;i++) {
			P_EQ_trans[i]=new double[Nst+1];
			for(j=0;j<=Nst;j++) P_EQ_trans[i][j]=P_EQ_state[i]*MARKOV_TM[i][j]-P_EQ_state[j]*MARKOV_TM[j][i];
		}
		gsl_matrix_free(A);
		gsl_vector_free(b);
		gsl_vector_free(x);
		gsl_vector_free(residual);
        
        delete[] state_i;state_i=NULL;
        delete[] state_j;state_j=NULL;
        delete[] s_n;s_n=NULL;
        
        weighted_network_of_Attractors(0);
    }
	
    void export_noisy_landscape_edgelist(double bet,int order){
		unsigned long long Nst,i,stid;
		unsigned int *state_j,k,g,*s_n;
		double p0,p1,p2;
		char fn[600];
		FILE *f;
		snode::slinklist::slink *sc2;

		
        beta=bet;
		error_p=1./(1.+exp(2*beta));
		Nst=(unsigned long long)pow(2.,(double)BRN->N);
		printf("Nst=%lld\n",Nst);
		p0=-log(pow(1.-error_p,(double)BRN->N));
        p1=-log(pow(1.-error_p,(double)BRN->N-1)*error_p);
        p2=-log(2*pow(1.-error_p,(double)BRN->N-2)*error_p*error_p);
        
      //  p0=pow(1.-error_p,(double)BRN->N);
      //  p1=pow(1.-error_p,(double)BRN->N-1)*error_p;
      //  p2=2*pow(1.-error_p,(double)BRN->N-2)*error_p*error_p;

        
		if(!BRN->Evolved) sprintf(fn,"mkdir %s/%s/Approx_Edgelist\n",MODEL_Directory_system,BRN->path);
        
		system(fn);
		
		//synchronous links
		sprintf(fn,"%s/%s/Approx_Edgelist/%s_W-TM_%d.txt",
				MODEL_Directory,BRN->path,BRN->name,order);
		f=fopen(fn,"w");
		for(i=1;i<=State_G->N;i++){
			sc2=State_G->node[i].Li.first_link_pt;
			while(sc2!=NULL){
				fprintf(f,"%lld\t%lld\t%lg\n",sc2->sneighbor,i,p0);
				sc2=sc2->next(sc2);
			}
		}
                
		s_n=new unsigned int[BRN->N+1];
		state_j=new unsigned int[BRN->N+1];
		
        //noisy transition links
		
		for(i=1;i<=Nst;i++){ // on all possible system states
			set_state(i);
			synchronously_update_all_Gates();
            for(g=1;g<=BRN->N;g++){
                for(k=0;k<BRN->N;k++) state_j[BRN->N-k]=s[k];
                state_j[g-1]=1-s[g-1];   // flipped one gate, prob p*(1-p)^{N-1}
                stid=anybase_to_decimal(state_j,BRN->N,2);
                fprintf(f,"%lld\t%lld\t%lg\n",i,stid+1,p1);
                if(order>1){
                    for(int g2=g+1;g2<=BRN->N;g2++){
                        for(k=0;k<BRN->N;k++) state_j[BRN->N-k]=s[k];
                        state_j[g-1]=1-s[g-1];   // flipped two gates, first permut  prob p*p*(1-p)^{N-2}
                        state_j[g2-1]=1-s[g2-1];   // flipped two gates, second permut prob p*p*(1-p)^{N-2}
                        stid=anybase_to_decimal(state_j,BRN->N,2);
                        fprintf(f,"%lld\t%lld\t%lg\n",i,stid+1,p2);
                    }
                }
           }
        }
        fclose(f);
	}
    
	void weighted_network_of_Attractors(int writeout){
		unsigned long long i,j;
		double lw;
		char fn[600];
		sprintf(fn,"%s_NETW_OF_ATTRACTORS",BRN->name);
		Attractor_G=new sgraph(Attractor_NR,fn);
		for(i=1;i<=Attractor_NR;i++){
			Attractor_G->node[i].node_props=new node_attr(1.0);
			Attractor_G->node[i].node_props->nw=0.;
			for(j=1;j<=Attractor_NR;j++)
				//if(j!=i) 
					Attractor_G->add_link(i,j,0.);
		}
		for(i=1;i<=State_G->N;i++){
			Attractor_G->node[State_G->node[i].ID].node_props->nw+=P_EQ_state[i];
			//printf("State %lld: Prob. of attractor %lld increased by %lg   to %lg\n",
			//	   i,State_G->node[i].ID,P_EQ_state[i],Attractor_G->node[State_G->node[i].ID].node_props->nw);
			//getchar();
		}
		for(i=1;i<=State_G->N;i++){
			for(j=1;j<=State_G->N;j++)
				//if(State_G->node[i].ID!=State_G->node[j].ID) 
				{
					lw=Attractor_G->node[State_G->node[i].ID].get_link_w(State_G->node[j].ID);
					//printf("Attractor Link weight  %lld -> %lld going from %lg   to  ",
					//	   State_G->node[i].ID,State_G->node[j].ID,lw);
					lw+=P_EQ_state[i]*MARKOV_TM[i][j]/Attractor_G->node[State_G->node[i].ID].node_props->nw;
					//printf("%lg\n",lw);
					//getchar();
					Attractor_G->node[State_G->node[i].ID].set_link_w(State_G->node[j].ID,lw);
				}
		}
		if(writeout){
            strcpy(fn,"");
            sprintf(fn,"%s/%s",
				MODEL_Directory,BRN->path);
            write_weighted_network(Attractor_G,fn);
            for(i=1;i<=Attractor_NR;i++){
                printf("Weight of attractor %lld = %lg\n",i,Attractor_G->node[i].node_props->nw);
            }
        }
        
	}
	
	void synchronously_update_all_Gates(){
		unsigned short int *s_n;
		unsigned int i;
		
		s_n=new unsigned short int[BRN->N+1];
		for(i=1;i<=BRN->N;i++) s_n[i]=update_Gate(i);
		for(i=1;i<=BRN->N;i++) s[i]=s_n[i];
		delete[] s_n;s_n=NULL;
	}
    
    short int synchronously_update_all_Gates_with_noise(double p_error){
        unsigned short int *s_n;
        unsigned int i;
        short int ok;
        double p;
        
        s_n=new unsigned short int[BRN->N+1];
        for(i=1;i<=BRN->N;i++) s_n[i]=update_Gate(i);
        ok=0;
        for(i=1;i<=BRN->N;i++) {
            p=rand_real(0,1);
            if(p<=p_error) {
                if(s_n[i]==1) s_n[i]=0; else s_n[i]=1;
            }
            if(s_n[i]!=s[i]) {s[i]=s_n[i];ok=1;}
        }
        delete[] s_n;s_n=NULL;
        return(ok);
    }
    
	void Asynchronously_update_all_Gates(){
		unsigned int i,l;
		double *x;
		
		x=new double[BRN->N+1];
		for(l=1;l<=BRN->N;l++) x[l]=l;
		randomize(x,BRN->N);
		
        for(l=1;l<=BRN->N;l++){
			i=(unsigned int)(x[l]);
			s[i]=update_Gate(i);
		}
		delete[] x;x=NULL;
	}
	
    void Asynchronously_update_all_Gates(longint_ARRAY_pair *Hits, double p_hits[]){
        unsigned int i,l;
        double *x;
        
        x=new double[BRN->N+1];
        for(l=1;l<=BRN->N;l++) x[l]=l;
        randomize(x,BRN->N);  // order goes from 1 to N
        
        longint_ARRAY *skip;
        skip=new longint_ARRAY();
        double r;
        if(Hits!=NULL){
            r=rand_real(0,1);
            for(unsigned long int l=1;l<=Hits->N;l++){
                if(r<=p_hits[l]) {
                    s[Hits->A[l]]=Hits->B[l];  // updated node Hits->A[l] by force at the START; otheriwse in random order
                    skip->add_element(Hits->A[l]);  // updated alreay, skip during random order!
                }
            }
        }
        for(l=1;l<=BRN->N;l++){
            i=(unsigned int)(x[l]);  // order goes from 1 to N
            if(Hits!=NULL){
                if(skip->check_for_element(i)==0) s[i]=update_Gate(i);
            }
            else s[i]=update_Gate(i);    // index order goes from 0 to N-1 but gate index is 1 to N
        }
        //        s[BRN->N]=0;
        
        delete skip; skip = NULL;
        delete[] x;x=NULL;
    }
    
    void Asynchronously_update_all_Gates(int preserve_order,longint_ARRAY_pair *Hits, double p_hits[]){
        unsigned int i,l,a;
        double *x;
        
        x=new double[BRN->N+1];
        for(l=1;l<=BRN->N;l++) x[l]=l;
        randomize(x,BRN->N);  // order goes from 1 to N
        
        if (preserve_order==1) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                //   for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                //   if (s[FoxM1-1]==48 + 1)     { iii++; a=x[O_FoxM1]; x[O_FoxM1]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                
                }
        }
        if (preserve_order==2) {    // Replication alone, symmetric
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int O_Rep=0;
                
                int iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
            }
        }
        if (preserve_order==3) {    // 4N_DNA
            if (!strcmp(Model_type,"Life_and_Death")) {
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                
                int iii=0,O_4N=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==4) {    // Replic + 4N_DNA
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                
                int O_Rep=0,O_4N=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==5) {    // U_Kinetochores
            if (!strcmp(Model_type,"Life_and_Death")) {
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                
                int O_uK=0;
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==6) {    // A_Kinetochores
            if (!strcmp(Model_type,"Life_and_Death")) {
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                
                int O_AK=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==7) {    // both kinetochores
            if (!strcmp(Model_type,"Life_and_Death")) {
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                
                int O_AK=0,O_uK=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==8) {    // FoxM1
            if (!strcmp(Model_type,"Life_and_Death")) {
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                
                int O_FoxM1=0;
                int iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==9) {    // Ect2
            if (!strcmp(Model_type,"Life_and_Death")) {
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                
                int O_cytokin=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
            }
        }
        if (preserve_order==10) {    // 4N_DNA + Ect2
            if (!strcmp(Model_type,"Life_and_Death")) {
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                
                int O_cytokin=0,O_4N=0;
                
                int iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==11) {    // Plk1_H
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int O_Plk_H=0;
                int iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==12) {    // Cdc20
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int O_Cdc20=0;
                int iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==13) {    // Pre-RC
            if (!strcmp(Model_type,"Life_and_Death")) {
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int O_PreRC=0;
                int iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
            }
        }
        if (preserve_order==14) {    // CyclinE
            if (!strcmp(Model_type,"Life_and_Death")) {
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                int O_CycE=0;
                int iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==15) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int O_CycB=0;
                int iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
            }
        }
        // 16: without 4N_DNA
        if (preserve_order==16) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                //   for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                //   if (s[FoxM1-1]==48 + 1)     { iii++; a=x[O_FoxM1]; x[O_FoxM1]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 17: without Replicaiton
        if (preserve_order==17) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 18: without U_Kinetochores
        if (preserve_order==18) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 19: without A_Kinetochores
        if (preserve_order==19) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 20: without FoxM1
        if (preserve_order==20) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 21: without Ect2
        if (preserve_order==21) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                //   for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                //   if (s[FoxM1-1]==48 + 1)     { iii++; a=x[O_FoxM1]; x[O_FoxM1]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                
                //  if ((s[uK-1]==48 + 1)||(s[AK-1]==48 + 1)) {
                //                    printf("Order: %s=%d %s=%d %s=%d %s=%d %s=%d %s=%d %s=%d .... %s=%d %s=%d %s=%d %s=%d %s=%d %s=%d %s=%d\n",
                //                           BRN->MDAT->Node_name[(unsigned int)x[1]],s[(unsigned int)x[1]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[2]],s[(unsigned int)x[2]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[3]],s[(unsigned int)x[3]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[4]],s[(unsigned int)x[4]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[5]],s[(unsigned int)x[5]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[6]],s[(unsigned int)x[6]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[7]],s[(unsigned int)x[7]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[BRN->N-6]],s[(unsigned int)x[BRN->N-6]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[BRN->N-5]],s[(unsigned int)x[BRN->N-5]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[BRN->N-4]],s[(unsigned int)x[BRN->N-4]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[BRN->N-3]],s[(unsigned int)x[BRN->N-3]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[BRN->N-2]],s[(unsigned int)x[BRN->N-2]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[BRN->N-1]],s[(unsigned int)x[BRN->N-1]-1]-48,
                //                           BRN->MDAT->Node_name[(unsigned int)x[BRN->N  ]],s[(unsigned int)x[BRN->N  ]-1]-48);
                //               //     after=1;
                // getchar();
                //  }
                //  else after=0;
            }
        }
        // 22: without Plk1_H
        if (preserve_order==22) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 23: without Cdc20
        if (preserve_order==23) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                //   for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                //   if (s[FoxM1-1]==48 + 1)     { iii++; a=x[O_FoxM1]; x[O_FoxM1]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 24: without Pre-RC
        if (preserve_order==24) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 25: without CyclinE
        if (preserve_order==25) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 26: without CyclinB
        if (preserve_order==26) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 27: no U_Kinetochores and no Replication -- BEST
        if (preserve_order==27) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 28: no U_Kinetochores and no Replication and no Ect2
        if (preserve_order==28) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_FoxM1=0,O_AK=0,O_4N=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 29: no U_Kinetochores at start
        if (preserve_order==29) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
            }
        }
        // 30: no U_Kinetochores at end
        if (preserve_order==30) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 31: no Replication at start
        if (preserve_order==31) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                iii=0;
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48 + 1)          { iii++; a=x[O_Rep]; x[O_Rep]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        // 32: no Replication at end
        if (preserve_order==32) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Repl)  O_Rep=l;}
                if (s[Repl-1]==48)          { iii++; a=x[O_Rep]; x[O_Rep]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                
                
                iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48 + 1)      { iii++; a=x[O_uK];       x[O_uK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        //33: no U_Kinetochores at end; no Replication
        if (preserve_order==33) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,
                O_uK=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==uK)  O_uK=l;}
                if (s[uK-1]==48)            { iii++; a=x[O_uK]; x[O_uK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48)         { iii++; a=x[O_Cdc20]; x[O_Cdc20]=x[iii];  x[iii] = a;}
                
                
                
                iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        if (preserve_order==34) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int f4N=BRN->MDAT->get_node_ID("4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Ect2");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre-RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_FoxM1=0,O_cytokin=0,O_AK=0,O_4N=0,O_Cdc20=0,O_Plk_H=0,O_CycB=0,O_PreRC=0,O_CycE=0;
                
                int iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==PreRC)  O_PreRC=l;}
                if (s[PreRC-1]==48 + 1)          { iii++; a=x[O_PreRC]; x[O_PreRC]=x[iii];  x[iii] = a; }
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48 + 1)        { iii++; a=x[O_AK]; x[O_AK]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48 + 1)     { iii++; a=x[O_Plk_H]; x[O_Plk_H]=x[iii];  x[iii] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cycB)  O_CycB=l;}
                if (s[cycB-1]==48 + 1)      { iii++; a=x[O_CycB]; x[O_CycB]=x[iii];  x[iii] = a;}
                
                iii=0;
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==f4N)  O_4N=l;}
                //   if (s[O_4N-1]==48 + 1)
                {iii++; a=x[O_4N];       x[O_4N]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin)  O_cytokin=l;}
                if (s[cytokin-1]==48)     { iii++; a=x[O_cytokin];  x[O_cytokin]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==AK)  O_AK=l;}
                if (s[AK-1]==48)          { iii++; a=x[O_AK];       x[O_AK]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==CycE)  O_CycE=l;}
                if (s[CycE-1]==48 + 1)       { iii++; a=x[O_CycE];    x[O_CycE]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==FoxM1)  O_FoxM1=l;}
                if (s[FoxM1-1]==48 + 1)       { iii++; a=x[O_FoxM1];    x[O_FoxM1]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Cdc20)  O_Cdc20=l;}
                if (s[Cdc20-1]==48 + 1)   { iii++; a=x[O_Cdc20];    x[O_Cdc20]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==Plk_H)  O_Plk_H=l;}
                if (s[Plk_H-1]==48)       { iii++; a=x[O_Plk_H];    x[O_Plk_H]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
            }
        }
        
        longint_ARRAY *skip;
        skip=new longint_ARRAY();
        double r;
        if(Hits!=NULL){
            r=rand_real(0,1);
            for(unsigned long int l=1;l<=Hits->N;l++){
                if(r<=p_hits[l]) {
                    s[Hits->A[l]]=Hits->B[l];  // updated node Hits->A[l] by force at the START; otheriwse in random order
                    skip->add_element(Hits->A[l]);  // updated alreay, skip during random order!
                }
            }
        }
        for(l=1;l<=BRN->N;l++){
            i=(unsigned int)(x[l]);  // order goes from 1 to N
            if(Hits!=NULL){
                if(skip->check_for_element(i)==0) s[i]=update_Gate(i);
            }
            else s[i]=update_Gate(i);    // index order goes from 0 to N-1 but gate index is 1 to N
        }
//        s[BRN->N]=0;
        
        delete skip; skip = NULL;
        delete[] x;x=NULL;
    }
    
	void set_state(unsigned long long stid){
		unsigned int *s_n,i;
		
		s_n=new unsigned int[BRN->N+1];
		decimal_to_any_base(stid-1,2,BRN->N,s_n);
		for(i=0;i<BRN->N;i++)
			s[BRN->N-i]=s_n[i];
		delete[] s_n;s_n=NULL;
	}
    
    void set_state(const char s_n[]){
		unsigned int i;
		for(i=0;i<BRN->N;i++)
			s[i+1]=s_n[i]-48;
        s[BRN->N+1]=0;
	}
    
	void print_state(){
		unsigned int i;
		for(i=1;i<=BRN->N;i++){ printf("%d ",s[i]);
			if(BRN->MDAT!=NULL) printf("(%s)  ",BRN->MDAT->Node_name[i]);
		}
		printf("\n");
	}
	
    char *print_string_for_state(unsigned long long stid){
        char *a;
        unsigned long long stid_old;
        stid_old=get_state_now();
        set_state(stid);
        a=new char[BRN->N+1];
        for(int i=1;i<=BRN->N;i++) a[i-1]=s[i]+48;
        a[BRN->N]=0;
        set_state(stid_old);
        return(a);
    }
    
	unsigned long long get_state_now(){
		unsigned int *s_n,i;
		unsigned long long stid;
		s_n=new unsigned int[BRN->N+1];
		
		for(i=0;i<BRN->N;i++) s_n[i]=s[BRN->N-i];
		stid=anybase_to_decimal(s_n,BRN->N,2);
		delete[] s_n;s_n=NULL;
		return(stid+1);
	}
	
	unsigned int get_Hamming_distance_of_states(unsigned long long s1,unsigned long long s2){
		unsigned int *s_1,*s_2,i,h;
		// no setting of state in the meantime!
		s_1=new unsigned int[BRN->N+1];
		decimal_to_any_base(s1-1,2,BRN->N,s_1);
		s_2=new unsigned int[BRN->N+1];
		decimal_to_any_base(s2-1,2,BRN->N,s_2);
		h=0;
		for(i=0;i<BRN->N;i++)
			if(s_1[i]!=s_2[i]) h++;
		return(h);
	}
	unsigned int get_Hamming_distance_of_states(unsigned  int s_1[],unsigned int s_2[]){
		unsigned int i,h;
		// no setting of state in the meantime!
		h=0;
		for(i=0;i<BRN->N;i++)
			if(s_1[i]!=s_2[i]) h++;
		return(h);
	}
    
    unsigned int get_Hamming_distance_of_states(const char s_1[],const char s_2[]){
		unsigned int i,h;
		// no setting of state in the meantime!
		h=0;
		for(i=0;i<BRN->N;i++)
			if(s_1[i]!=s_2[i]) h++;
		return(h);
	}

	unsigned long long next_state_synchronous(){
		synchronously_update_all_Gates();
		return(get_state_now());	
	}
	
	unsigned long long next_state_Asynchronous(){
		Asynchronously_update_all_Gates();
		return(get_state_now());	
	}
	
	void mark_Attractor_State_from(unsigned long long a_0){
		unsigned long long j,j_old;
		
		set_state(a_0);
		//printf("Setting state to %lld\n",a_0); print_state();
		j=a_0;
		do{ State_G->node[j].node_props->basin=0;
			j_old=j;
			j=next_state_synchronous();
			//printf("j=%lld\n",j); print_state();getchar();
		}
		while(j!=a_0);
	}
	
	void print_Attractor_Cycle(unsigned int aid){
		unsigned long long i;
		printf("Attractor Cycle %d with %lld states:\n",aid,Attractor_Cycle[aid]->N);
		for(i=1;i<=Attractor_Cycle[aid]->N;i++){
			 printf(" %lld - ",Attractor_Cycle[aid]->A[i]); 
			 set_state(Attractor_Cycle[aid]->A[i]);
			 print_state();
		}
	}
    
    unsigned int largest_H_from_all_attractors(){
        unsigned int  max_H_all=0,min_H,h;
        
        for (unsigned long int i=1; i<=pow(2,BRN->N); i++) {
            min_H=BRN->N;
            for(int j=1;j<=Attractor_NR;j++){
                for(int l=1;l<=Attractor_Cycle[j]->N;l++){
                    h=get_Hamming_distance_of_states(i,Attractor_Cycle[j]->A[l]);
                    if (h<min_H) min_H=h;
                }
            }
//            printf("State %ld best match min_H=%d\n",i, min_H);
            if (min_H>max_H_all) {max_H_all=min_H;
                //set_state(i);print_state();
            }
        }
        //printf("largest_H_from_all_attractors H_max_all=%d\n",max_H_all);getchar();
        return(max_H_all);
    }
    
    double normalized_distance_between_attractor_pair(unsigned int aid1, unsigned int aid2){
        double h=0;
        for(int i=1;i<=Attractor_Cycle[aid1]->N;i++){
            for(int j=1;j<=Attractor_Cycle[aid2]->N;j++){
                h+=get_Hamming_distance_of_states(Attractor_Cycle[aid1]->A[i],Attractor_Cycle[aid2]->A[j]);
            }
        }
        h/=(double)(Attractor_Cycle[aid1]->N*Attractor_Cycle[aid2]->N);
        h/=(double)BRN->N;
        return(h);
    }
    
    double normalized_distance_between_attractors(){
        double h=0;
        if(Attractor_NR==1) return(0);
        for(int aid1=2;aid1<=Attractor_NR;aid1++)
            for(int aid2=1;aid2<aid1;aid2++)
                h+=normalized_distance_between_attractor_pair(aid1,aid2);
        h/=(double)((Attractor_NR/2.)*(Attractor_NR-1.));
        printf("\t\t Module - N_attr=%lu\t H_av=%lg  ",Attractor_NR,h);
        h/=maximum_distance_between_p_Binary_string_pairs((int)Attractor_NR);
        printf(" ideal Hav = %lg  => nornalized H_av = %lg\n",maximum_distance_between_p_Binary_string_pairs((int)Attractor_NR),h);//getchar();
        return(h);
    }
    
    void read_attractors(){
        char fn[600];
        FILE *f;
        long long int n,bb;
        
        sprintf(fn,"%s/%s/%s_Attractors.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"r");
        if(f==NULL){
            generate_full_state_graph_synchronously(0);
            
            f=fopen(fn,"w");
            fprintf(f,"%lu\n",Attractor_NR);
 
            for(int aid=1;aid<=Attractor_NR;aid++){
                fprintf(f,"%lld\t",Attractor_Cycle[aid]->N);
                for(int i=1;i<=Attractor_Cycle[aid]->N;i++){
                    fprintf(f,"%lld ",Attractor_Cycle[aid]->A[i]);
                }
                fprintf(f,"\n");
            }
        }
        else{
            fscanf(f,"%lu\n",&Attractor_NR);
            Attractor_Cycle=new longlongint_ARRAY_list[Attractor_NR+1];
            for(int i=1;i<=Attractor_NR;i++) Attractor_Cycle[i]=NULL;
            
            for(int aid=1;aid<=Attractor_NR;aid++){
                Attractor_Cycle[aid]=new longlongint_ARRAY();
                fscanf(f,"%lld\t",&n);
                for(int i=1;i<=n;i++){
                    fscanf(f,"%lld ",&bb);
                    Attractor_Cycle[aid]->add_element(bb);
                }
            }
        }
		fclose(f);
    }
    
    char *get_attractor_state_string(unsigned int  a,unsigned int  aindex){
        char *ss;
        set_state(Attractor_Cycle[a]->A[aindex]);
        
        ss=new char[BRN->N+1];
        for(int i=1;i<=BRN->N;i++) ss[i-1]=s[i]+48;
        ss[BRN->N]=0;
        return(ss);
    }
	
	longlongint_ARRAY *get_Attractor_Cycle_from(unsigned long long a_0){
		unsigned long long j,j_old;
		longlongint_ARRAY *ac;
		
		ac=new longlongint_ARRAY();
		set_state(a_0);
		//printf("Setting state to %lld\n",a_0); print_state();
		j=a_0;
		do{ ac->add_element(j);
			j_old=j;
			j=next_state_synchronous();
			//printf("j=%lld\n",j); print_state();getchar();
		}
		while(j!=a_0);
		return(ac);
	}
	
	longlongint_ARRAY *list_Attractor_State_from(unsigned long long a_0){
		unsigned long long j,j_old;
		longlongint_ARRAY *asl;
		asl=new longlongint_ARRAY();
		asl->add_element(a_0);
		set_state(a_0);
		//print_state();
		j=a_0;
		do{ j_old=j;
			j=next_state_synchronous();
			asl->add_element(j);
		}
		while(j!=a_0);
		return(asl);
	}
	
	longlongint_ARRAY *list_timetrace_to_first_Attractor_State(unsigned long long a_0){
		unsigned long long j,j_old;
		longlongint_ARRAY *asl;
		
		printf("starting timetrace of %lld\n",a_0);
		asl=new longlongint_ARRAY();
		asl->add_element(a_0);
		set_state(a_0);
		print_state();
		j=a_0;
		do{ j_old=j;
			j=next_state_synchronous();
			asl->add_element(j);
			print_state();
			//printf("\t\t%lld basin %lld\n",j,State_G->node[j].node_props->basin);getchar();
		}
		while(State_G->node[j].node_props->basin>0);
		return(asl);
	}
	
	void set_State_weights(){
		longlongint_ARRAY *asl;
		unsigned long long i,j,k;
		
		for(i=1;i<=State_G->N;i++){
			asl=list_timetrace_to_first_Attractor_State(i);
			for(j=2;j<=asl->N;j++) State_G->node[asl->A[j]].node_props->nw+=(j-1);
			for(k=1;k<=State_G->N;k++)
				if((k!=asl->A[j])&&(State_G->node[k].ID==State_G->node[i].ID)&&(State_G->node[j].node_props->basin==0))
					State_G->node[k].node_props->nw+=(asl->N-1); 
					// this is just the inflowing tree, 
					// the attractor state reached is added when considered alone
			//delete asl;asl=NULL;
		}
	}
    
    
    void export_STT_as_state_list_pairs(){
        FILE *net;
        unsigned long long i;
        snode::slinklist::slink *sc2;
   
        char fn[400];
        sprintf(fn,"%s/%s/Cytoskape/%s_STG.txt",MODEL_Directory,BRN->path,BRN->name);
        net=fopen(fn,"w");

        for(i=1;i<=State_G->N;i++){
            sc2=State_G->node[i].Li.first_link_pt;//printf("%d ",sc2->sneighbor->ID);getchar();
            while(sc2!=NULL){
                fprintf(net,"%s\t%s\n",print_string_for_state(i),print_string_for_state(sc2->sneighbor));
                sc2=sc2->next(sc2);//getchar();
            }
        }
        fclose(net);
        
        sprintf(fn,"%s/%s/Cytoskape/%s_ATTRACTORS.txt",MODEL_Directory,BRN->path,BRN->name);
        net=fopen(fn,"w");
        for(i=1;i<=Attractor_NR;i++){
            fprintf(net,"Attractor Cycle %llu with %lld states:\n",i,Attractor_Cycle[i]->N);
            for(int ai=1;ai<=Attractor_Cycle[i]->N;ai++)
                 fprintf(net,"%s\n",print_string_for_state(Attractor_Cycle[i]->A[ai]));
        }
        fclose(net);

        sprintf(fn,"%s/%s/Cytoskape/%s_STATE_Probs.txt",MODEL_Directory,BRN->path,BRN->name);
        net=fopen(fn,"w");
        for(i=1;i<=State_G->N;i++){
            fprintf(net,"%s\t%g\n",print_string_for_state(i),P_EQ_state[i]);
        }
        fclose(net);
    }

	void generate_full_state_graph_synchronously(int mark_and_weights){
		unsigned long long i,j,j_old,k;
		longlongint_ARRAY_pair *remap;
		
		State_G=new sgraph((unsigned long long)pow(2.,(double)BRN->N));
		for(i=1;i<=State_G->N;i++) State_G->node[i].ID=0;
		Attractor_NR=0;
		for(i=1;i<=State_G->N;i++){
			State_G->node[i].node_props=new node_attr(1.0);
			State_G->node[i].node_props->basin=1;
		}
		remap=new longlongint_ARRAY_pair();
		int n1=0,n2=0;
        
        for(i=1;i<=State_G->N;i++){
			if(State_G->node[i].ID==0){
				Attractor_NR++;
				set_state(i);
				// if(BRN->N>10){
                //print_state(); printf("Basin %ld\n",Attractor_NR);//getchar();
                //}
				j=i;
				do{ State_G->node[j].ID=Attractor_NR;
					j_old=j;
					j=next_state_synchronous();
					State_G->add_link(j,j_old);
					// if(BRN->N>10) {print_state();getchar();}
                  //  print_state(); printf("Basin %lld\n",State_G->node[j].ID);//getchar();
                    
				}
				while(State_G->node[j].ID==0);
				
				if(State_G->node[j].ID!=Attractor_NR){
					
                    //ind1=remap->check_for_element_A(State_G->node[j].ID);
					//remap->add_element(i,remap->B[ind1]);
                    for(k=1;k<=State_G->N;k++) {
                        if(State_G->node[k].ID==Attractor_NR) State_G->node[k].ID=State_G->node[j].ID;
                    
                    }

                   // printf("\t\tBasin %ld mapped back into %lld\n",Attractor_NR,State_G->node[j].ID);
                    
					Attractor_NR--;
					//printf("\n\ni=%lld\tTimecourse run into attractor marked from %lld\t remapped to %lld\n",i,
					//	   State_G->node[j].ID,remap->B[ind1]);//getchar();
				}
                
				else {
                   // printf("\n\nTimecourse found new attractor %ld\n",Attractor_NR);getchar();
					mark_Attractor_State_from(j); // marking = basin is set to 0
					//remap->add_element(i, Attractor_NR);
					// if(BRN->N>10) getchar();
				}
            
			}
            if(State_G->node[i].ID==1) n1++;
            if(State_G->node[i].ID==2) n2++;
		}
        //printf("n1=%d\tn2=%d\n",n1,n2);getchar();
		Attractor_Cycle=new longlongint_ARRAY_list[Attractor_NR+1];
		for(i=1;i<=Attractor_NR;i++) Attractor_Cycle[i]=NULL;
		
		for(k=1;k<=State_G->N;k++){
			//ind1=remap->check_for_element_A(State_G->node[k].ID);
			//if(ind1>0) {
			//	State_G->node[k].ID=remap->B[ind1];
			//	if((State_G->node[k].node_props->basin==0)&&(Attractor_Cycle[State_G->node[k].ID]==NULL))
			//		Attractor_Cycle[State_G->node[k].ID]=get_Attractor_Cycle_from(k);
				//printf("State %lld belongs to attractor %lld\n",k,remap->B[ind1]);
			//}
            
            if((State_G->node[k].node_props->basin==0)&&(Attractor_Cycle[State_G->node[k].ID]==NULL))
                Attractor_Cycle[State_G->node[k].ID]=get_Attractor_Cycle_from(k);
			
		}
		delete remap;remap=NULL;
		
		if(mark_and_weights) set_State_weights();
	}
	
	void generate_full_state_graph_Asynchronously(unsigned int MAXL_N,unsigned int MAXL_ST){
		unsigned long long i,j,j_old,st,stst;
		
		State_G=new sgraph((unsigned long long)pow(2.,(double)BRN->N));
		for(i=1;i<=State_G->N;i++) State_G->node[i].ID=0;
		Attractor_NR=0;
		for(i=1;i<=State_G->N;i++){
			State_G->node[i].node_props=new node_attr(1.0);
			State_G->node[i].node_props->basin=1;
			State_G->node[i].node_props->nw=0.;
		}
		for(stst=1;stst<=MAXL_ST;stst++){
			for(i=1;i<=State_G->N;i++){
				set_state(i);
				j=i;
				for(st=1;st<=MAXL_N*BRN->N;st++){
					j_old=j;
					j=next_state_Asynchronous();
					State_G->add_link_and_count(j,j_old);
					State_G->node[j].node_props->nw+=1.;
					if(j==j_old) st=MAXL_N*BRN->N+1;
				}
			}
		}
	}

	longlongint_ARRAY_pair *get_ID_and_size_of_functional_attractors(unsigned long int n){
		double *Bs;
		unsigned long long *aid,i;
		longlongint_ARRAY_pair *fa;
		Bs=new double[Attractor_NR+1];
		aid=new unsigned long long[Attractor_NR+1];
		for(i=0;i<=Attractor_NR;i++) Bs[i]=0;
		for(i=1;i<=State_G->N;i++) Bs[State_G->node[i].ID]++;
		hpsort_coshuffle_indices(Attractor_NR,Bs, aid);
		fa=new longlongint_ARRAY_pair();
		for(i=1;i<=n;i++){ 
			if(i<=Attractor_NR) 
				fa->add_element(aid[Attractor_NR-i+1],(unsigned long long)(Bs[Attractor_NR-i+1]));
			else fa->add_element(i-Attractor_NR,0);
		}
		delete[] Bs;Bs=NULL;
		delete[] aid;aid=NULL;
		return(fa);
	}
	
	double get_Coverage_fitness(longlongint_ARRAY_pair *fa){
		double c;
		unsigned long long i;
		c=0;
		for(i=1;i<=fa->N;i++) c+=fa->B[i];
		c/=pow(2.,(double)BRN->N);
		return(c);
	}
	
	double get_Basin_Similarity_fitness(longlongint_ARRAY_pair *fa){
		double c,s,x;
		unsigned long long i;
		c=0;
		for(i=1;i<=fa->N;i++) c+=fa->B[i];
		s=0;
		for(i=1;i<=fa->N;i++) {
			x=fa->B[i]/c;
			if(x>0) s+=-x*log(x);
		}
		if(fa->N>1) s/=log(fa->N);
		return(s);
	}
	
    double get_normalized_Basin_Size_entropy(){
		double c,s,x,s0,*Bs;
		unsigned long long i;
		
        Bs=new double[Attractor_NR+1];
		for(i=0;i<=Attractor_NR;i++) Bs[i]=0;
		for(i=1;i<=State_G->N;i++) Bs[State_G->node[i].ID]++;
		
        c=0;
		for(i=1;i<=Attractor_NR;i++) c+=Bs[i];
		s=0;
		for(i=1;i<=Attractor_NR;i++) {
			x=Bs[i]/c;
			if(x>0) s+=-x*log(x);
		}
		if(Attractor_NR>1) s/=log(Attractor_NR);
		
        x=1./(double)Attractor_NR;
        s0=-x*log(x);
        s0*=Attractor_NR;
        if(Attractor_NR>1) s0/=log(Attractor_NR);
        return(s/s0);
       
	}
    double get_normalized_occupation_probability_entropy(){
		double s,x,s0;
		unsigned long long i;
		
        if (Attractor_NR==1) return(0);
        
        s=0;
        for(i=1;i<=Attractor_NR;i++) {
			x=Attractor_G->node[i].node_props->nw;
			if(x>0) s+=-x*log(x);
		}
		if(Attractor_NR>1) s/=log(Attractor_NR);
		
        x=1./(double)Attractor_NR;
        s0=-x*log(x);
        s0*=Attractor_NR;
        if(Attractor_NR>1) s0/=log(Attractor_NR);
        return(s/s0);
	}

	
	double get_One_flip_from_attractor_fitness(longlongint_ARRAY_pair *fa){
		double p_attr,p_av,p;
		unsigned long long i,j,stid;
		unsigned int n,k,l,*s_n,*s_p;
		
		p_av=0;
		s_n=new unsigned int[BRN->N+1];
		s_p=new unsigned int[BRN->N+1];
		
		for(i=1;i<=fa->N;i++)
			if(fa->B[i]>0){
				p_attr=0;n=0;
				for(j=1;j<=State_G->N;j++){
					if((State_G->node[j].ID==fa->A[i])&&(State_G->node[j].node_props->basin==0)){
						p=0;n++;
						decimal_to_any_base(j,2,BRN->N,s_n);
						for(k=0;k<BRN->N;k++){ // all node flips from attractor state j
							for(l=0;l<BRN->N;l++) s_p[l]=s_n[l];
							s_p[k]=1-s_p[k];	
							stid=anybase_to_decimal(s_p,BRN->N,2);
							if(State_G->node[stid].ID==fa->A[i]) p++;
						}
						p/=(double)BRN->N; // fraction of flips that stay from A state j
					//	printf("Attractor index %lld Attr ID %lld state %lld   p=%lg\n",i,fa->A[i],j,p);getchar();
						p_attr+=p; 
					}
				}
				p_attr/=(double)(n); // average fraction for all attractor states in fa->A[i]
				p_av+=p_attr;
				//printf("Attractor index %lld Attr ID %lld state average  p_attr=%lg\n",i,fa->A[i],p_attr);getchar();
		}
		p_av/=(double)(fa->N); //average over all funcitonal attractors
	//	printf("Funcitonal Attractors average  p_av=%lg\n",p_av);getchar();
		delete[] s_n;s_n=NULL;
		delete[] s_p;s_p=NULL;
		return(p_av);
	}
/*
    void print_all_oneflips_from_attractor(int AID){
        unsigned int *s_n,*s_p;
		unsigned long long stid;
		longlongint_ARRAY *gl;
        
        s_n=new unsigned int[BRN->N+1];
		s_p=new unsigned int[BRN->N+1];
       	
        for(unsigned long int j=1;j<=State_G->N;j++){
            if(State_G->node[j].ID==AID){
                decimal_to_any_base(j,2,BRN->N,s_n);
                for(int k=0;k<BRN->N;k++){ // all node flips from attractor state j
                    printf("Flipping node %d %s\n",k+1,BRN->MDAT->Node_name[k+1]);
                    for(int l=0;l<BRN->N;l++) s_p[l]=s_n[l];
                    s_p[k]=1-s_p[k];
                    stid=anybase_to_decimal(s_p,BRN->N,2);
                    gl=list_timetrace_to_first_Attractor_State(stid);
                    if(State_G->node[stid].ID==AID) printf("Same attractor maintained\n");
                    else {
                        printf("Flipping to Attractor %llu\n",State_G->node[stid].ID);
                        print_Attractor_Cycle((unsigned int)State_G->node[stid].ID);
                    }
                    getchar();
                }
            }
        }
        delete[] s_n;s_n=NULL;
		delete[] s_p;s_p=NULL;
    }
	*/
    
    void print_all_oneflips_from_attractor(int AID){
        unsigned int *s_n,*s_p;
		unsigned long long stid;
		longlongint_ARRAY *gl;
        
        s_n=new unsigned int[BRN->N+1];
		s_p=new unsigned int[BRN->N+1];
       	
        for(unsigned long int j=1;j<=State_G->N;j++){
            if((State_G->node[j].ID==AID)&&(State_G->node[j].node_props->basin==0)){
                printf("\n\n\nFrom attractor state:\n");
                set_state(j);
                print_state();
                for(int k=1;k<=BRN->N;k++){ // all node flips from attractor state j
                    set_state(j);
                    printf("Flipping node %d %s\n",k,BRN->MDAT->Node_name[k]);
                    s[k]=1-s[k];
                    print_state();
                    stid=get_state_now();
                    gl=list_timetrace_to_first_Attractor_State(stid);
                    if(State_G->node[stid].ID==AID) printf("Same attractor maintained\n");
                    else {
                        printf("Flipping to Attractor %llu\n",State_G->node[stid].ID);
                        print_Attractor_Cycle((unsigned int)State_G->node[stid].ID);
                    }
                    getchar();
                }
            }
        }
        delete[] s_n;s_n=NULL;
		delete[] s_p;s_p=NULL;
    }
    
    void print_all_TWOflips_from_attractor(int AID){
        unsigned int *s_n,*s_p;
		unsigned long long stid;
		longlongint_ARRAY *gl;
        
        s_n=new unsigned int[BRN->N+1];
		s_p=new unsigned int[BRN->N+1];
       	
        for(unsigned long int j=1;j<=State_G->N;j++){
            if((State_G->node[j].ID==AID)&&(State_G->node[j].node_props->basin==0)){
                printf("\n\n\nFrom attractor state:\n");
                set_state(j);
                print_state();
                for(int k=2;k<=BRN->N;k++){ // all node flips from attractor state j
                    for(int k2=1;k2<k;k2++){
                        printf("Flipping nodes %d %s and %d %s\n",k,BRN->MDAT->Node_name[k],k2,BRN->MDAT->Node_name[k2]);
                        s[k]=1-s[k];s[k2]=1-s[k2];
                        print_state();
                        stid=get_state_now();
                        gl=list_timetrace_to_first_Attractor_State(stid);
                        if(State_G->node[stid].ID==AID) printf("Same attractor maintained\n");
                        else {
                            printf("Flipping to Attractor %llu\n",State_G->node[stid].ID);
                            print_Attractor_Cycle((unsigned int)State_G->node[stid].ID);
                            getchar();
                        }
                        
                    }
                }
            }
        }
        delete[] s_n;s_n=NULL;
		delete[] s_p;s_p=NULL;
    }

    
	double get_Double_flip_from_attractor_fitness(longlongint_ARRAY_pair *fa){
		double p_attr,p_av,p;
		unsigned long long i,j,stid;
		unsigned int n,k,kk,l,*s_n,*s_p;
		
		p_av=0;
		s_n=new unsigned int[BRN->N+1];
		s_p=new unsigned int[BRN->N+1];
		
		for(i=1;i<=fa->N;i++)
			if(fa->B[i]>0){
				p_attr=0;n=0;
				for(j=1;j<=State_G->N;j++){
					if((State_G->node[j].ID==fa->A[i])&&(State_G->node[j].node_props->basin==0)){
						p=0;n++;
						decimal_to_any_base(j,2,BRN->N,s_n);
						for(k=1;k<BRN->N;k++)// all double flips from attractor state j
							for(kk=0;kk<k;kk++){ 
								for(l=0;l<BRN->N;l++) s_p[l]=s_n[l];
								s_p[k]=1-s_p[k];s_p[kk]=1-s_p[kk];		
								stid=anybase_to_decimal(s_p,BRN->N,2);
								if(State_G->node[stid].ID==fa->A[i]) p++;
						}
						p/=(double)(BRN->N*(BRN->N-1.)/2.); // fraction of flips that stay from A state j
						p_attr+=p; 
					}
				}
				p_attr/=(double)(n); // average fraction for all attractor states in fa->A[i]
				p_av+=p_attr;
			}
		p_av/=(double)(fa->N); //average over all funcitonal attractors
		delete[] s_n;s_n=NULL;
		delete[] s_p;s_p=NULL;
		return(p_av);
	}
	
	void export_state_Graph_to_Cytoskape_no_noise(){
		unsigned long long i;
		
		char fn[600];
		FILE *f;
		snode::slinklist::slink *sc2;
        
		if(!BRN->Evolved) sprintf(fn,"mkdir %s/%s/Cytoskape\n",
								  MODEL_Directory_system,BRN->path);
        
		system(fn);
		
		//links
		sprintf(fn,"%s/%s/Cytoskape/%s_LINKS.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"start_State\tis_followed_by\tend_State\tFlux\n");
		for(i=1;i<=State_G->N;i++){
			sc2=State_G->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				fprintf(f,"%lld\tis_followed_by\t%lld\n",sc2->sneighbor,i);
				sc2=sc2->next(sc2);
			}
		}
		fclose(f);
        
        sprintf(fn,"%s/%s/Cytoskape/%s_Boolean_Net_LINKS.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"start_Node\tregulates\tend_Node\n");
		for(i=1;i<=BRN->N;i++){
			sc2=BRN->BN->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				if(BRN->MDAT==NULL)
                    fprintf(f,"%lld\tis_controlled_by\t%lld\n",i,sc2->sneighbor);
                else fprintf(f,"%s\tis_controlled_by\t%s\n",BRN->MDAT->Node_name[i],BRN->MDAT->Node_name[sc2->sneighbor]);
				sc2=sc2->next(sc2);
			}
		}
		fclose(f);
        
        sprintf(fn,"%s/%s/Cytoskape/%s_Attractors.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"Attractor_ID\n");
		for(i=1;i<=State_G->N;i++) fprintf(f,"%lld = %lld\n",i,State_G->node[i].ID);
		fclose(f);	
		
        
        for(int a=1;a<=Attractor_NR;a++) {
            //D->print_Attractor_Cycle(a);
            for(int as=1;as<=Attractor_Cycle[a]->N;as++){
                sprintf(fn,"%s/%s/Cytoskape/%s_Boolean_Net_Attractor_%d_S-%d.txt",
                        MODEL_Directory,BRN->path,BRN->name,a,as);
                f=fopen(fn,"w");
                fprintf(f,"A-%d__%d\n",a,as);
                set_state(Attractor_Cycle[a]->A[as]);
               // printf("%s\n",s);
                for(i=1;i<=BRN->N;i++) fprintf(f,"%s = %d\n",BRN->MDAT->Node_name[i],s[i]);
                fclose(f);
            }
        }
    }
    
	void export_state_Graph_to_Cytoskape(){
		unsigned long long i,j;
		unsigned int k;
		char fn[600];
		FILE *f;
		snode::slinklist::slink *sc2;

		if(!BRN->Evolved) sprintf(fn,"mkdir %s/%s/Cytoskape\n",
								  MODEL_Directory_system,BRN->path);

		system(fn);
		
		//links
		sprintf(fn,"%s/%s/Cytoskape/%s_LINKS.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"start_State\tis_followed_by\tend_State\tFlux\n");
		for(i=1;i<=State_G->N;i++){
			sc2=State_G->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				fprintf(f,"%lld\tis_followed_by\t%lld\t%lg\n",sc2->sneighbor,i,MARKOV_TM[sc2->sneighbor][i]);
				sc2=sc2->next(sc2);
			}
		}
		fclose(f);
		
		sprintf(fn,"%s/%s/Cytoskape/%s_MATRIX.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"start_State\tis_followed_by\tend_State\tFLOW\tSYNC_UPD_LINK\n");
		for(i=1;i<=State_G->N;i++){
			for(j=1;j<=State_G->N;j++)
				if((i!=j)&&(MARKOV_TM[j][i]>0.0001)){
					if(State_G->node[i].get_directed_link(j) != NULL)
						 fprintf(f,"%lld\tis_followed_by\t%lld\t%lg\tSYNC\n",j,i,MARKOV_TM[j][i]);
					else fprintf(f,"%lld\tis_followed_by\t%lld\t%lg\tNOT_SYNC\n",j,i,MARKOV_TM[j][i]);
				}
			sc2=State_G->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				fprintf(f,"%lld\tis_followed_by\t%lld\t%lg\tSYNC\n",sc2->sneighbor,i,MARKOV_TM[sc2->sneighbor][i]);
				sc2=sc2->next(sc2);
			}
		}
		fclose(f);

		
		sprintf(fn,"%s/%s/Cytoskape/%s_Energy_of_state_space.noa",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"Energy\n");
		for(i=1;i<=State_G->N;i++){
			if(P_EQ_state[i]>0.00000000000000000000000000001)
				fprintf(f,"%lld = %lg\n",i,-log(P_EQ_state[i])/beta);
			else fprintf(f,"%lld = %lg\n",i,-log(0.00000000000000000000000000001)/beta);
            //printf("Node %lld -  p=%lg\t e=%lg\t (beta=%lg)\n",i,P_EQ_state[i],-log(P_EQ_state[i])/beta,beta);
        }
        fclose(f);
		//Boolean network links

		sprintf(fn,"%s/%s/Cytoskape/%s_Boolean_Net_LINKS.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"start_Node\tregulates\tend_Node\n");
		for(i=1;i<=BRN->N;i++){
			sc2=BRN->BN->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				if(BRN->MDAT==NULL)
                     fprintf(f,"%lld\tis_controlled_by\t%lld\n",i,sc2->sneighbor);
                else fprintf(f,"%s\tis_controlled_by\t%s\n",BRN->MDAT->Node_name[i],BRN->MDAT->Node_name[sc2->sneighbor]);
				sc2=sc2->next(sc2);
			}
		}
		fclose(f);
		
		sprintf(fn,"%s/%s/Cytoskape/%s_Attractors.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"Attractor_ID\n");
		for(i=1;i<=State_G->N;i++) fprintf(f,"%lld = %lld\n",i,State_G->node[i].ID);
		fclose(f);	
		sprintf(fn,"%s/%s/Cytoskape/%s_Weight.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"Weight\n");
		for(i=1;i<=State_G->N;i++) fprintf(f,"%lld = %lg\n",i,State_G->node[i].node_props->nw);
		fclose(f);	
		sprintf(fn,"%s/%s/Cytoskape/%s_Not_covered.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"Fraction_of_states_flowing_through\n");
		for(i=1;i<=State_G->N;i++) fprintf(f,"%lld = %lg\n",i,(double)(State_G->N-State_G->node[i].node_props->nw)/(double)State_G->N);
		fclose(f);	
		
		sprintf(fn,"%s/%s/Cytoskape/%s_Expression.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"Expression\n");
		for(i=1;i<=State_G->N;i++){
			fprintf(f,"%lld = \"",i);
			set_state(i);
			for(k=1;k<=BRN->N;k++) fprintf(f,"%d",s[k]);
			fprintf(f,"\"\n");
		}
		fclose(f);	
		
	}

	void export_ASYNC_state_Graph_to_Cytoskape(){
		unsigned long long i;
		unsigned int k;
		char fn[600];
		FILE *f;
		snode::slinklist::slink *sc2;
		
		if(!BRN->Evolved) sprintf(fn,"mkdir %s/%s/Cytoskape\n",
								  MODEL_Directory_system,BRN->path);
		
		system(fn);
		
		//links
		sprintf(fn,"%s/%s/Cytoskape/ASYNC_%s_LINKS.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"start_State\tis_followed_by\tend_State\tWeight\n");
		for(i=1;i<=State_G->N;i++){
			sc2=State_G->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				fprintf(f,"%lld\tis_followed_by\t%lld\t%lg\n",sc2->sneighbor,i,sc2->link_props->lw);
				sc2=sc2->next(sc2);
			}
		}
		fclose(f);
		//Boolean network links
		
		sprintf(fn,"%s/%s/Cytoskape/%s_Boolean_Net_LINKS.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"start_Node\tregulates\tend_Node\n");
		for(i=1;i<=BRN->N;i++){
			sc2=BRN->BN->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				fprintf(f,"%lld\tregulates\t%lld\n",i,sc2->sneighbor);
				sc2=sc2->next(sc2);
			}
		}
		fclose(f);
		
		sprintf(fn,"%s/%s/Cytoskape/ASYNC_%s_Weight.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"Weight\n");
		for(i=1;i<=State_G->N;i++) fprintf(f,"%lld = %lg\n",i,State_G->node[i].node_props->nw);
		fclose(f);	
		
		sprintf(fn,"%s/%s/Cytoskape/ASYNC_%s_Expression.txt",
				MODEL_Directory,BRN->path,BRN->name);
		f=fopen(fn,"w");
		fprintf(f,"Expression\n");
		for(i=1;i<=State_G->N;i++){
			fprintf(f,"%lld = \"",i);
			set_state(i);
			for(k=1;k<=BRN->N;k++) fprintf(f,"%d",s[k]);
			fprintf(f,"\"\n");
		}
		fclose(f);	
		
	}
	
    double get_link_P_Zero_above_H_errors(const char s_1[],const char s_2[], int *ST_order,double H_ERR){
        // procedure needs to restore current state of system
        char *s_remember,*s_target;
        double f;
        int h1;
        s_remember=new char[BRN->N+1];
        for(int i=0;i<=BRN->N;i++) s_remember[i]=s[i];
        s_remember[BRN->N]=0;
        
        set_state(s_1);
        synchronously_update_all_Gates();
        
        s_target=new char[BRN->N+1];
        for(int i=0;i<=BRN->N;i++) s_target[i]=s[i];
        s_target[BRN->N]=0;
        h1=get_Hamming_distance_of_states(s_2,s_target);
        *ST_order=h1;
        if(h1<=H_ERR) f=pow(error_p,(double)h1)*pow(1.-error_p,(double)(BRN->N-h1));
        else f=0.;
        set_state(s_remember);
        delete[] s_remember;s_remember=NULL;
        delete[] s_target;s_target=NULL;
        return(f);
    }
    
    double get_link_P(const char s_1[],const char s_2[], int *ST_order){
        // procedure needs to restore current state of system
        char *s_remember,*s_target;
        double f;
        int h1;
        s_remember=new char[BRN->N+1];
        for(int i=0;i<=BRN->N;i++) s_remember[i]=s[i];
        s_remember[BRN->N]=0;
        
        set_state(s_1);
        synchronously_update_all_Gates();
        
        s_target=new char[BRN->N+1];
        for(int i=0;i<=BRN->N;i++) s_target[i]=s[i];
        s_target[BRN->N]=0;
        h1=get_Hamming_distance_of_states(s_2,s_target);
        *ST_order=h1;
        f=pow(error_p,(double)h1)*pow(1.-error_p,(double)(BRN->N-h1));
        set_state(s_remember);
        delete[] s_remember;s_remember=NULL;
        delete[] s_target;s_target=NULL;
        return(f);
    }
    
    void clean_landscape_calc(){
        unsigned long long Nst=(unsigned long long)pow(2.,(double)BRN->N);
        if(MARKOV_TM!=NULL){
            for(unsigned long long i=0;i<=Nst;i++) {delete[] MARKOV_TM[i]; MARKOV_TM[i]=NULL;}
            delete[] MARKOV_TM; MARKOV_TM=NULL;
        }
        if(P_EQ_trans!=NULL){
            for(unsigned long long i=0;i<=Nst;i++) {delete[] P_EQ_trans[i]; P_EQ_trans[i]=NULL;}
            delete[] P_EQ_trans; P_EQ_trans=NULL;
            delete[] P_EQ_state; P_EQ_state=NULL;
        }
    }
    
    
	~Boolean_Dynamics(){
		//delete BRN;BRN=NULL;
		delete[] s;s=NULL;
		if(State_G!=NULL) {delete State_G;State_G=NULL;}
        //if(Attractor_G!=NULL) {delete Attractor_G;Attractor_G=NULL;}
	}
};

typedef Boolean_Dynamics *Boolean_Dynamics_list;

void generate_Attractor_number_histogram(unsigned int N, unsigned int L,unsigned int STAT){
	Boolean_RegNet *A;
	Boolean_Dynamics *D;
	unsigned int *hist,i;
	unsigned long int *rnd_id;
	double *av_Basin_simil,*flip_comb;
	FILE *f;
	char fn[600];
	longlongint_ARRAY_pair *fa;
	double p1,p2,*overall_fitness,bsfit;
	
	hist=new unsigned int[N*N+1];
	av_Basin_simil=new double[N*N+1];
	flip_comb=new double[N*N+1];
	overall_fitness=new double[STAT+1];
	for(i=0;i<=N*N;i++) {hist[i]=0;flip_comb[i]=0;av_Basin_simil[i]=0;}
	for(i=1;i<=STAT;i++){
		A=new Boolean_RegNet(N,L,i);
		D=new Boolean_Dynamics(A);
		D->generate_full_state_graph_synchronously(0);
		//if(D->Attractor_NR==2) D->export_state_Graph_to_Cytoskape();
		if(D->Attractor_NR>N*N) {printf("More attractors than N*N=%u \t%ld!!!\n",N*N,D->Attractor_NR);}
		else {
			hist[D->Attractor_NR]++;
			//if(i/100==i/100.) 
			//printf("seed=%d has %ld attractors\n",i,D->Attractor_NR);
			//for(k=1;k<=D->Attractor_NR;k++)
			//	printf("\t %d: %lld     ",k,D->Attractor_Cycle[k]->N);
			//printf("\n\n"); 
			fa=  D->get_ID_and_size_of_functional_attractors(D->Attractor_NR);
			p1=  D->get_One_flip_from_attractor_fitness(fa);
			p2=  D->get_Double_flip_from_attractor_fitness(fa);
			
			flip_comb[D->Attractor_NR]+=p1+p2*p2;
			bsfit=D->get_Basin_Similarity_fitness(fa);
			av_Basin_simil[D->Attractor_NR]+=bsfit;
			if((D->Attractor_NR==2)&&(bsfit>0.95)&&(D->Attractor_Cycle[1]->N==1)&&(D->Attractor_Cycle[2]->N==1)) {//
				overall_fitness[i]=p1+p2*p2;
				printf("i=%u \t flip fitness = %lg\n",i,p1+p2*p2);
			}
			else overall_fitness[i]=0;
		}
		//delete D;D=NULL;A=NULL;
	}
	rnd_id=new unsigned long int[STAT+1];
	hpsort_coshuffle_indices(STAT,overall_fitness, rnd_id);
	
	printf("TOP 20 for 2 point attractors ~ equal size, out of first rnd %d\n",STAT);
	for(i=1;i<=20;i++)
		printf("Rank %d\t RND_ID = %ld   fitness = %lg\n",i,rnd_id[STAT-i+1],overall_fitness[STAT-i+1]);
	
	delete[] rnd_id;rnd_id=NULL;	
	delete[] overall_fitness;overall_fitness=NULL;	
	sprintf(fn,"%s/Rnd_NoBowTie_BRN/N-%d_L-%d__Attractor_nr_hist_STAT-%d.dat",
			MODEL_Directory,N,L,STAT);
	f=fopen(fn,"w");
	for(i=0;i<=N*N;i++) if(hist[i]>0) fprintf(f,"%d\t%d\t%lg\t%lg\n",i,hist[i],av_Basin_simil[i]/(double)hist[i],flip_comb[i]/(double)hist[i]);
	fclose(f);
	exit(1);
	
}

void generate_histograms_and_Scatter_of_robustness(unsigned int N, unsigned int L,unsigned int STAT, unsigned int n_funct){
	Boolean_RegNet *A;
	Boolean_Dynamics *D;
	unsigned int *hist_BC,*hist_BS,*hist_FC,i,BOX=100;
	FILE *f,*f_2;
	char fn[600];
	longlongint_ARRAY_pair *fa;
	double p1,p2,cover,fc,bs;
	
	hist_BC=new unsigned int[BOX+1];
	hist_BS=new unsigned int[BOX+1];
	hist_FC=new unsigned int[BOX+1];
	for(i=0;i<=N*N;i++) {hist_BC[i]=0;hist_BS[i]=0;hist_FC[i]=0;}
	
	sprintf(fn,"%s/Rnd_NoBowTie_BRN/Rnd_NoBowTie_BRN_%d_%d__Roubustness_scatter_Nfunc-%d_STAT-%d.dat",
			MODEL_Directory,N,L,n_funct,STAT);
	f=fopen(fn,"w");
	sprintf(fn,"%s/Rnd_NoBowTie_BRN/Rnd_NoBowTie_BRN_%d_%d__Roubustness_scatter_ATTRNR-%d_STAT-%d.dat",
			MODEL_Directory,N,L,n_funct,STAT);
	f_2=fopen(fn,"w");
	for(i=1;i<=STAT;i++){
		A=new Boolean_RegNet(N,L,i);
		D=new Boolean_Dynamics(A);
		D->generate_full_state_graph_synchronously(0);
		//printf("seed=%d has %ld attractors\n",i,D->Attractor_NR);
		//for(k=1;k<=D->Attractor_NR;k++)
		//	printf("\t %d: %lld     ",k,D->Attractor_Cycle[k]->N);
		//printf("\n\n");
		
		fa=  D->get_ID_and_size_of_functional_attractors(n_funct);
		cover=D->get_Coverage_fitness(fa);
		p1=  D->get_One_flip_from_attractor_fitness(fa);
		p2=  D->get_Double_flip_from_attractor_fitness(fa);
		bs=  D->get_Basin_Similarity_fitness(fa);
		fc= (p1+p2*p2);

		hist_BC[(unsigned int)(cover*BOX)]++;
		hist_BS[(unsigned int)(bs*BOX)]++;
		hist_FC[(unsigned int)(fc*BOX/2.)]++;
		fprintf(f,"%lg\t%lg\t%lg\n",cover,bs,fc);
		if(D->Attractor_NR==2) {
			fprintf(f_2,"%lg\t%lg\t%lg\n",cover,bs,fc);
		}
		if(i/100==i/100.) printf("stat %d\n",i);
		if((cover>0.999)&&(bs>0.9)&&(fc>1)){
			printf("NR_attr=%lu   %d basin similarity %lg    flip_resistance %lg\n",D->Attractor_NR, i,bs,fc);
		}
	}
	fclose(f);	
	fclose(f_2);
	sprintf(fn,"%s/Rnd_NoBowTie_BRN/Rnd_NoBowTie_BRN_%d_%d__Robust-Coverage_distrib_Nfunc-%d_STAT-%d.dat",
			MODEL_Directory,N,L,n_funct,STAT);
	f=fopen(fn,"w");
	for(i=0;i<=BOX;i++) fprintf(f,"%d\t%d\n",i,hist_BC[i]);
	fclose(f);	
	sprintf(fn,"%s/Rnd_NoBowTie_BRN/Rnd_NoBowTie_BRN_%d_%d__Robust-BSimilarity_distrib_Nfunc-%d_STAT-%d.dat",
			MODEL_Directory,N,L,n_funct,STAT);
	f=fopen(fn,"w");
	for(i=0;i<=BOX;i++) fprintf(f,"%d\t%d\n",i,hist_BS[i]);
	fclose(f);	
	sprintf(fn,"%s/Rnd_NoBowTie_BRN/Rnd_NoBowTie_BRN_%d_%d__Robust-StateFlip_distrib_Nfunc-%d_STAT-%d.dat",
			MODEL_Directory,N,L,n_funct,STAT);
	f=fopen(fn,"w");
	for(i=0;i<=BOX;i++)  fprintf(f,"%d\t%d\n",i,hist_FC[i]);
	fclose(f);	
}


void Bio_Boolean_Network(const char bbn[],int calc_landscape){
	Boolean_RegNet *A;
	Boolean_Dynamics *D;
	unsigned int i;
	double cover,p1,p2,bs,fc;
	longlongint_ARRAY_pair *fa;
	
	beta=4.;
//	A=new Boolean_RegNet(bbn,bbn,"HS");
	A=new Boolean_RegNet(bbn);
	A->export_Boolean_RegNet();
	D=new Boolean_Dynamics(A);
    D->generate_full_state_graph_synchronously(0);
	D->export_state_Graph_to_Cytoskape_no_noise(); 
    for(i=1;i<=D->Attractor_NR;i++) D->print_Attractor_Cycle(i);

    
    if(calc_landscape){
        if(A->N<=11){
            D->calculate_noisy_landscape(beta);
            D->BRN->export_Boolean_RegNet();
            D->export_state_Graph_to_Cytoskape();
        }
    }
    D->export_noisy_landscape_edgelist(beta,2);
    exit(1);
    
	fa=  D->get_ID_and_size_of_functional_attractors(2);
	cover=D->get_Coverage_fitness(fa);
	p1=  D->get_One_flip_from_attractor_fitness(fa);
	p2=  D->get_Double_flip_from_attractor_fitness(fa);
	bs=  D->get_Basin_Similarity_fitness(fa);
	fc= (p1+p2*p2);
	
	printf("N=%d  K=%d  Cover = %lg\n Basin_size_entropy=%lg\n Flip_fitness=%lg\n",A->N,A->K,cover,bs,fc);
}


void Boolean_RegNet_Modeling(){
	//char pn[600];
	set_attractor_szinnev();
	
//	Bio_Boolean_Network("Cell_cycle_Faure");
//	Bio_Boolean_Network("Cell_Cycle_Faure_noCycD");
//	Bio_Boolean_Network("Apop_minimal_noCasp3_selfreg");
//	Bio_Boolean_Network("Apop_minimal");
//	Bio_Boolean_Network("Tip_stalk_MARY");
//	Bio_Boolean_Network("Apop_Casp3_overrides_XIAP");
//	Bio_Boolean_Network("Helikar_Fibroblast_GPCRs");
//	Bio_Boolean_Network("Tip_stalk_1Cell");
//	Bio_Boolean_Network("Tip_stalk_2Cell");
    
    //Bio_Boolean_Network("G0_G1");
	//Bio_Boolean_Network("S_M",0);
	
    //Bio_Boolean_Network("CellCycle_Erzso_new",0);
	//Bio_Boolean_Network("CellCycle_Restriction_SW",1);
	//Bio_Boolean_Network("CellCycle_Mitosis_SW",1);
    
    Bio_Boolean_Network("CC_DNA_Damage",1);
	
    //	Bio_Boolean_Network("_EC_DM_network");
	
//	generate_Attractor_number_histogram(5,12,5000);
//	generate_histograms_and_Scatter_of_robustness(5,12,15000,2);
//	study_composites_SYNC();
//	generate_histograms_and_Scatter_of_robustness(10,36,15000,2);
}


void generate_Basin_entropy_dependent_rnd_atreactor_transition_pics(unsigned int N, unsigned int L,unsigned int A_NR,unsigned int BOX,unsigned int STAT_box,double beta){
	Boolean_RegNet *A;
	Boolean_Dynamics *D=NULL;
	//FILE *f,*f_2;
	//char fn[600];
    int *n_box,i;
	double *x;
	doublepointer *W_av,*W_s;
    triplepointer *TR_av,*TR_s;
    unsigned long int *orig_id,*ord_a;
	
	//sprintf(fn,"%s/Rnd_NoBowTie_BRN/Rnd_NoBowTie_BRN_%d_%d__Roubustness_scatter_Nfunc-%d_STAT-%d.dat",MODEL_Directory,N,L,n_funct,STAT);
//	f=fopen(fn,"w");
    n_box=new int[BOX+1];
    W_av=new doublepointer[BOX+1];W_s=new doublepointer[BOX+1];
    TR_av=new triplepointer[BOX+1];TR_s=new triplepointer[BOX+1];
    for (i=0;i<=BOX; i++) {
        n_box[i]=0;
        W_av[i]=new double[A_NR+1];W_s[i]=new double[A_NR+1];
        TR_av[i]=new doublepointer[BOX+1];TR_s[i]=new doublepointer[BOX+1];
        for (int a=1;a<=A_NR;a++) {
            W_av[i][a]=0;W_s[i][a]=0;
            TR_av[i][a]=new double[A_NR+1];TR_s[i][a]=new double[A_NR+1];
            for (int b=1;b<=A_NR;b++) {
                TR_av[i][a][b]=0;TR_s[i][a][b]=0;
            }
        }
    }
    x=new double[BOX+1]; orig_id=new unsigned long int[BOX+1];ord_a=new unsigned long int[BOX+1];
        

    snode::slinklist::slink *sc2;
	    
    i=0;int boxes_not_full=0;
    while(!boxes_not_full) {
        i++;
		A=new Boolean_RegNet(N,L,i);
		D=new Boolean_Dynamics(A);
		D->generate_full_state_graph_synchronously(0);
		double H=D->get_normalized_Basin_Size_entropy();
        
        if ((D->Attractor_NR==A_NR)&&(n_box[(int)(H*BOX)]<STAT_box)) {
            D->calculate_noisy_landscape(beta);
            n_box[(int)(H*BOX)]++;
          //  if((int)(H*BOX)==10)
            //    {printf("\n%d H=%lg   n[%d]=%d\n",i,H,(int)(H*BOX),n_box[(int)(H*BOX)]);getchar();}
            for (int a=1; a<=D->Attractor_NR; a++) {
                x[a]=D->Attractor_G->node[a].node_props->nw;
               // if((int)(H*BOX)==2)
                 //   printf("\t W[%d]=%lg  links: ",a, D->Attractor_G->node[a].node_props->nw);
                sc2=D->Attractor_G->node[a].Li.first_link_pt;
                while(sc2!=NULL){
                   // if((int)(H*BOX)==2)
                     //   printf("%d->%llu %lg  ",a,sc2->sneighbor,sc2->link_props->lw);
                    sc2=sc2->next(sc2);
                }
             //   printf("\n");
            }
            hpsort_coshuffle_indices(D->Attractor_NR,x,orig_id);
            for (int a=1; a<=D->Attractor_NR; a++) ord_a[orig_id[a]]=a;
            
            for (int a=1; a<=D->Attractor_NR; a++){
                W_av[(int)(H*BOX)][a]+=x[a]; W_s[(int)(H*BOX)][a]+=x[a]*x[a];
                sc2=D->Attractor_G->node[orig_id[a]].Li.first_link_pt;
                while(sc2!=NULL){
                    TR_av[(int)(H*BOX)][a][ord_a[sc2->sneighbor]]+=sc2->link_props->lw;
                    TR_s[(int)(H*BOX)][a][ord_a[sc2->sneighbor]]+=sc2->link_props->lw*sc2->link_props->lw;
                    sc2=sc2->next(sc2);
                }
            }
           /* if((int)(H*BOX)==2){
                printf("stats n= %d: ",n_box[(int)(H*BOX)]);
                for (int a=1; a<=D->Attractor_NR; a++){
                    printf("\n Weight %d:  %lg ; %lg\n",a,W_av[(int)(H*BOX)][a],W_s[(int)(H*BOX)][a]);
                    for (int b=1; b<=D->Attractor_NR; b++)
                        printf("\t ->%d: %lg ; %lg    ", b, TR_av[(int)(H*BOX)][a][b],TR_s[(int)(H*BOX)][a][b]);
                }
                getchar();
            }
            */
        }
       // else printf("i=%d N_attr=%lu\n",i,D->Attractor_NR);
        boxes_not_full=1;
        for (int j=1;j<=BOX; j++) if(n_box[j]<STAT_box) boxes_not_full=0;
        //if(i>100) boxes_not_full=1;
	}
    
    for (int j=1;j<=BOX; j++){
        printf("\n\nBOX=%d stats n= %d: ",j,n_box[j]);
        for (int a=1; a<=D->Attractor_NR; a++){
            W_av[j][a]/=(double)STAT_box;
            W_s[j][a]=sigma_from_average_and_sum_of_squares(W_av[j][a],W_s[j][a],STAT_box);
            printf("\n Weight %d:  %lg ; %lg\n",a,W_av[j][a],W_s[j][a]);
            for (int b=1; b<=D->Attractor_NR; b++){
                TR_av[j][a][b]/=(double)STAT_box;
                TR_s[j][a][b]=sigma_from_average_and_sum_of_squares(TR_av[j][a][b],TR_s[j][a][b],STAT_box);
                printf("\t ->%d: %lg ; %lg    ", b, TR_av[j][a][b],TR_s[j][a][b]);
            }
            //getchar();
        }
    }
    
	//fclose(f);
}

void basin_occupancy_entropy_vs_barrier_breach(unsigned int N, unsigned int L,unsigned int A_NR,unsigned int BOX,int STAT,int STAT_start,int RND_start,double beta){
	Boolean_RegNet *A;
	Boolean_Dynamics *D;
	int i;
	snode::slinklist::slink *sc2;
    
    i=RND_start-1;int ok=STAT_start;
    while(ok<STAT) {
        i++;
		A=new Boolean_RegNet(N,L,i);
		D=new Boolean_Dynamics(A);
		D->generate_full_state_graph_synchronously(0);
		
        if (D->Attractor_NR==A_NR) {
            ok++;
            if(STAT_start<=ok){
                D->calculate_noisy_landscape(beta);
                double H=D->get_normalized_occupation_probability_entropy();
                
                double lwmax=0;
                for (int a=1; a<=D->Attractor_NR; a++){
                    sc2=D->Attractor_G->node[a].Li.first_link_pt;
                    while(sc2!=NULL){
                        if (sc2->sneighbor!=a) {
                            if (sc2->link_props->lw>lwmax) {
                                lwmax=sc2->link_props->lw;
                            }
                            //printf("\t\t%lg\t%lg\t%d\n",H,sc2->link_props->lw,ok);getchar();
                        }
                        sc2=sc2->next(sc2);
                    }
                }
                printf("%lg\t%lg\t%d\t\t%d\n",H,lwmax,ok,i);
                D->clean_landscape_calc();
            }
        }
        // else printf("i=%d N_attr=%lu\n",i,D->Attractor_NR);
        delete D;D=NULL;
	}

}


