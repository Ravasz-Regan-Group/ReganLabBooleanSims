
#ifndef _Boolean_modeling_sampling_2_h
#define _Boolean_modeling_sampling_2_h

typedef std::map <std::string, unsigned long long> State_index_MAP_type;

class neighbor_SYNC_const_noise{
public:
    char *s;
    
    neighbor_SYNC_const_noise(int n,const char sin[]){
        unsigned int i;
        s=new char[n+1];
        for(i=0;i<n;i++)s[i]=sin[i];
        s[n]=0;
    }
    ~neighbor_SYNC_const_noise(){
        delete[] s;s=NULL;
    }
};

typedef neighbor_SYNC_const_noise *neighbor_SYNC_const_noise_list;


struct Network_microstate_SYNC_const_noise{
public:
    char *s;   
    unsigned long long visits;
    neighbor_SYNC_const_noise_list * NextDoor;
};

typedef Network_microstate_SYNC_const_noise *Trajectory;

typedef std::map <std::string, Network_microstate_SYNC_const_noise> Network_microstate_SYNC_const_noise_MAP_type;


class Attractor_basin{
public:
    unsigned long int Basin_ID,N_a,Rarest_visit_count;
    unsigned long long N_ms,W_Total,N_MAX,*Attractor;
    State_index_MAP_type States_in_basin;
    double x,y;
    
//    sparse_hash_map<std::string, unsigned long long> States_in_basin; 
    //Rarest_states;
   // Trajectory *Attractor,*Basin_Stored;
    Trajectory *Basin_Stored;
    
    Attractor_basin(unsigned long int b,unsigned long long nm){
        Basin_ID=b;
        N_a=0;
        W_Total=0;
        N_ms=0;
        N_MAX=nm;
        Basin_Stored=new Trajectory[N_MAX+1];
        Rarest_visit_count=0;
    }
    
    ~Attractor_basin(){
        delete[] Attractor;Attractor=NULL;
        for(unsigned long long int i=1;i<=N_ms;i++) {
            delete Basin_Stored[i];Basin_Stored[i]=NULL;
        }
        delete[] Basin_Stored;Basin_Stored=NULL;
        States_in_basin.clear();
    }
    
    unsigned short int check_if_state_is_in_Basin(const char s_now[]){
        State_index_MAP_type::const_iterator It;
        
        It=States_in_basin.find(s_now);
        if(It != States_in_basin.end()) return(1);
        else return(0);
    }
    unsigned short int check_if_state_is_an_Attractor_State(const char s_now[]){
        unsigned int i;
        for(i=1;i<=N_a;i++) {
            //printf("%s__%s\n",s_now,Basin_Stored[Attractor[i]]->s);
            if(!strcmp(s_now,Basin_Stored[Attractor[i]]->s)) {
               // printf("\t ret: %d\n",i);
                return(i);
            }
        }
       //  printf("\t ret: 0\n");
        return(0);
    }
    
    unsigned long long get_visits_to_state__Basin(const char s_now[]){
        State_index_MAP_type::const_iterator It;
        
        It=States_in_basin.find(s_now);
        if(It != States_in_basin.end()) {
            return(Basin_Stored[It->second]->visits);
        }
        return(0);
    }

    
    unsigned short int check_if_state_is_an_Attractor_State(unsigned long long BS_id){
        unsigned int i;
        for(i=1;i<=N_a;i++) if(BS_id==Attractor[i]) return(1);
        return(0);
    }

    
    unsigned long long get_Microstate_and_mark_visit(const char s_now[]){
        State_index_MAP_type::const_iterator It;
        unsigned long long sid;
        
        It=States_in_basin.find(s_now);
        if(It != States_in_basin.end()) {
            sid=It->second;
            Basin_Stored[sid]->visits++;
            W_Total++;
//            printf("sid = %lld\t snow= %s (found %s)\t visits= %d\n",sid,s_now, );
            return(sid);
        }
        else {
            W_Total++;
            //printf("State %s was supposed to be in basin of %ld\n",s_now,Basin_ID);exit(1);
        }
        return (0);
    }
    unsigned long long get_or_increase_visitcount(){
        unsigned long long i;
        unsigned long long N_rare;
        
        N_rare=0;
        for(i=1;i<=N_ms;i++){
            if((Basin_Stored[i]->visits==Rarest_visit_count)&&(!check_if_state_is_an_Attractor_State(Basin_Stored[i]->s))) N_rare++;
        }
        if(N_rare>0) return(N_rare);
        
        //printf_attractor_full_basin_HashMap();
        while (N_rare==0) {
            //Rarest_states.clear();
            N_rare=0; //printf("2 Basin %ld N_rare=%llu\n",Basin_ID,N_rare);
            //            Rarest_states.set_deleted_key("______");
            //printf("Cleared pool with %lu visits\n",Rarest_visit_count);getchar();
            Rarest_visit_count++;
            for(i=1;i<=N_ms;i++){
                if(Basin_Stored[i]->visits==Rarest_visit_count){
                    if(!check_if_state_is_an_Attractor_State(Basin_Stored[i]->s)){
                        //Rarest_states[Basin_Stored[i]->s]=i;
                        N_rare++; //printf("3 Basin %ld N_rare=%llu\n",Basin_ID,N_rare);
                    }
                }
            }
            //printf_attractor_full_basin_HashMap();
            //getchar();
        }
        printf("Lowest visitcount pool for attractor %lu is %lu \n",Basin_ID,Rarest_visit_count);
        //getchar();
        return (N_rare);
    }
    
    void mark_visit(const char s_now[],unsigned long int ct){
        State_index_MAP_type::const_iterator It;
        
        It=States_in_basin.find(s_now);
        W_Total+=ct;
        if(It != States_in_basin.end()) {
            Basin_Stored[It->second]->visits+=ct;
        }
    }
    
    void printf_attractor_basin_info(){
     //   unsigned long int i;
        
        printf("Basin %ld with an attractor cycle length %ld\n\tW=%lld\tStates recorded in basin: %lld \n\tRarest  visitation count %lu\n\tAttractor:\n",
               Basin_ID,N_a,W_Total,N_ms,Rarest_visit_count);
      //  for(i=1;i<=N_a;i++) {
      //      printf("%lu %llu  ",i,Attractor[i]);
      //      printf("%s\n",Basin_Stored[Attractor[i]]->s);
      //  }
    //    printf("\ndone\n");
        //getchar();
    }
    
    void printf_attractor_basin_short_info(){
        
        printf("Basin %ld with an attractor cycle length %ld\n\tW=%lld\tStates recorded in basin: %lld \n\tRarest  visitation count %lu\n",
               Basin_ID,N_a,W_Total,N_ms,Rarest_visit_count);
        //getchar();
    }
    
    /*
    void printf_attractor_full_basin_HashMap(int NN){
        unsigned int i;
        State_index_MAP_type::const_iterator It;
        unsigned int j;
        
        printf("\nBasin %ld with an attractor cycle length %ld\n\tW=%lld\tStates recorded in basin: %lld (%lu)\n\tAttractor: ",
               Basin_ID,N_a,W_Total,N_ms,get_size_of_Basin_HashMap());
        for(i=1;i<=N_a;i++) printf("%s  ",Attractor[i]->s);
        printf("\n");
        
        for (It = States_in_basin.begin(); It != States_in_basin.end(); ++It) {
            if(Basin_Stored[It->second]->visits>0){
                printf("\tin basin: %s  v=%llu \n",
                       Basin_Stored[It->second]->s,Basin_Stored[It->second]->visits);
                for(j=0;j<=NN;j++){
                    printf("\t\t%d %s\n",j,
                           Basin_Stored[It->second]->NextDoor[j]->s);
                }
            }
        }
        //getchar();
    }
    */
    void printf_attractor_full_basin_HashMap(){
        State_index_MAP_type::const_iterator It;
        
      /*  printf("\nBasin %ld with an attractor cycle length %ld\n\tW=%lld\tStates recorded in basin: %lld (%llu)\n\tAttractor: ",
               Basin_ID,N_a,W_Total,N_ms,get_size_of_Basin_HashMap());
        
        for(i=1;i<=N_a;i++) printf("%s  ",Attractor[i]->s);
        printf("\n");
        */
        for (It = States_in_basin.begin(); It != States_in_basin.end(); ++It) {
         //   if(Basin_Stored[It->second]->visits>0){
                printf("\tin basin: %s  v=%llu \n",
                       Basin_Stored[It->second]->s,Basin_Stored[It->second]->visits);
           // }
        }
        //getchar();
    }
    
     void printf_rarest_states_HashMap(){
         State_index_MAP_type::const_iterator It;
     
         printf("\nBasin %ld rarest states pool:\n",Basin_ID);
     
         for (It = States_in_basin.begin(); It != States_in_basin.end(); ++It) {
             if(Basin_Stored[It->second]->visits==Rarest_visit_count)
                 printf("\tin basin: %s  v=%llu \n",
                        Basin_Stored[It->second]->s,Basin_Stored[It->second]->visits);
         }
         getchar();
     }
    
    /*
     unsigned long long get_size_of_rarest_states_HashMap(){
     unsigned long long nr;
     State_index_MAP_type::const_iterator It;
     nr=0;
     printf("Attractor %ld\n",Basin_ID);//getchar();
     
     // if(Rarest_states.empty()) {
     //   printf("Ratest_states_is_somehow_empty!\n");
     //  return (0);
     //}
     for (It = Rarest_states.begin(); It != Rarest_states.end(); ++It) {
     nr++;
     printf("\t\t\t%lld %s %lld - v-%llu\n",nr,Basin_Stored[It->second]->s,It->second,Basin_Stored[It->second]->visits);
     }
     printf("Rarest states has size %lld\n",nr);//getchar();
     return (nr);
     }
     */
    unsigned long long get_RND_Element_of_rarest_states_HashMap(){
        unsigned long long nr,rid,i,N_rare;
        
        N_rare=get_or_increase_visitcount();
        rid=rand_int(1,N_rare);
        //printf("Random choice of %llu out of %llu\n",rid,N_rare);
        nr=0;
        for(i=1;i<=N_ms;i++){
            if(Basin_Stored[i]->visits==Rarest_visit_count) {
                if(!check_if_state_is_an_Attractor_State(Basin_Stored[i]->s)) nr++;
            }
            if(rid==nr) return(i);
        }
        return (0);
    }
    unsigned long long get_size_of_Basin_HashMap(){
        unsigned long long nr;
        State_index_MAP_type::const_iterator It;
        nr=0;
        for (It = States_in_basin.begin(); It != States_in_basin.end(); ++It) {
            nr++;
        }
        return (nr);
    }
};

double get_energy(double p_state,double p_error){
    
    return(-2*(log(p_state)/(log(1-p_error)-log(p_error))));
}

typedef Attractor_basin *Attractor_landscape;

class Boolean_Dynamics_TRACKER{
public:	
	Boolean_RegNet *BRN;
	char *s;
	unsigned long int Attractor_NR;
    Attractor_landscape * Attractor_valleys;
    doublepointer  *Trans_Basins;
	unsigned long long N_MAX,W_total_steps;
    unsigned long int *Node_frozen_in;
    unsigned int *color;
    double *Node_average;
    char DIRname[200],SAMPLEname[200];
    //char *modular_chopped;
   // State_index_MAP_type Transition_track;
    sgraph_array *Attr_Transitions;
    
    
	Boolean_Dynamics_TRACKER(Boolean_RegNet *b,unsigned long long nm){
		BRN=b;
		s=new char[BRN->N+1];
        s[BRN->N]=0;
		Attractor_NR=0;
        Attractor_valleys = NULL;
        Trans_Basins=NULL;
        N_MAX=nm;
        W_total_steps=0;
        strcpy(DIRname,b->name);
        strcpy(SAMPLEname,b->name);
        Node_average=NULL;
        Node_frozen_in=NULL;
        color=NULL;
        Attr_Transitions=NULL;
    }
    
    Boolean_Dynamics_TRACKER(Boolean_RegNet *b,unsigned long long nm,const char dirnm[],const char spnm[]){
		
		BRN=b;
		s=new char[BRN->N+1];
        s[BRN->N]=0;
		Attractor_NR=0;
        Attractor_valleys = NULL;
        Trans_Basins=NULL;
        N_MAX=nm;
        W_total_steps=0;
        strcpy(DIRname,dirnm);
        strcpy(SAMPLEname,spnm);
        Node_average=NULL;
        Node_frozen_in=NULL;
        color=NULL;
    }
    unsigned long int get_NODE_ID_from_node_label_line(){
		unsigned long int ind1;
		char attr_list[50][5000];
        unsigned long int Nf;
        char * pch;
        
        pch = strstr (GSE_line,"label=");
        Nf=separate_line(pch,attr_list,'"',50);
        ind1=BRN->MDAT->get_node_ID(attr_list[1]);
		return(ind1);
	}
    
    void read_Node_Cytoskape_positions(){
		char fn[300];
		FILE *f;
		unsigned long int s,Nf;
		char attr_list[50][5000];
		double *x,*y;
		sprintf(fn,"%s/%s/POZ_%s.xgmml",
				MODEL_Directory,DIRname,SAMPLEname);
		f=NULL;
		while (f==NULL) {
			f=fopen(fn, "r");
			if(f==NULL){ printf("Please use exported Cytoskape link and attribute files to generate\n\t %s\n THEN hit Enter\n",fn);
				getchar();
			}
		}
		x=new double[BRN->N+1];
		y=new double[BRN->N+1];
        int *ok;
        
        ok=new int[BRN->N+1];
        for (int i=1; i<=BRN->N; i++) ok[i]=0;
        
		while (get_line_and_check(f)) {
			if(strstr(GSE_line,"<node ")!=NULL){
               // printf("found <node \n");getchar();
                s=get_NODE_ID_from_node_label_line();
               // printf("found %ld\n",s);getchar();
                
				get_line_and_check(f);
				while(strstr(GSE_line,"<graphics ")==NULL){
					get_line_and_check(f);
				}
                char * pch;
                pch = strstr (GSE_line,"x=");
                Nf=separate_line(pch,attr_list,'"',50);
				x[s]=atof(attr_list[1]);
               
                pch = strstr (GSE_line,"y=");
                Nf=separate_line(pch,attr_list,'"',50);
                y[s]=atof(attr_list[1]);
                ok[s]=1;
				//h=atof(attr_list[4]);
				//if(s>0) printf("s=%ld\t %s\t x=%lg  y=%lg\n",s,BRN->MDAT->Node_name[s],x[s],y[s]);//getchar();
				//if(h_min>h) h_min=h;
				//if(h_max<h) h_max=h;
			}
		}
		fclose(f);
        ok[0]=1;
        for (int i=1; i<=BRN->N; i++)
            if(ok[i]==0){
                printf("Node %s has no coordinates\n",BRN->MDAT->Node_name[i]);
                ok[0]=0;
            }
        if (ok[0]==0)  getchar();
        
        BRN->MDAT->x=x;
        BRN->MDAT->y=y;
	}
    
	~Boolean_Dynamics_TRACKER(){
		unsigned int i;
        
        //delete BRN;BRN=NULL;
		delete[] s;s=NULL;
        if(Trans_Basins!=NULL) {
            for(i=0;i<=Attractor_NR;i++) {delete[] Trans_Basins[i];Trans_Basins[i]=NULL;}
            delete[] Trans_Basins;Trans_Basins=NULL;
        }
        if(Attractor_valleys!=NULL){
            for(i=1;i<=Attractor_NR;i++) {
                for(unsigned long long int j=1;j<=Attractor_valleys[i]->N_ms;j++) {
                    delete[] Attractor_valleys[i]->Basin_Stored[j]->s;Attractor_valleys[i]->Basin_Stored[j]->s=NULL;
                }
                delete Attractor_valleys[i];Attractor_valleys[i]=NULL;
            }
            delete[] Attractor_valleys;Attractor_valleys=NULL;
        }
        if(Node_frozen_in!=NULL) {delete[] Node_frozen_in;Node_frozen_in=NULL;}
        if(Node_average!=NULL) {delete[] Node_average;Node_average=NULL;}
        if(color!=NULL) {delete[] color;color=NULL;}
        if(Attr_Transitions!=NULL) {
            for (unsigned int i=0; i<=BRN->N; i++) {
                if(Attr_Transitions[i]!=NULL){delete Attr_Transitions[i];Attr_Transitions[i]=NULL;}
            delete[] Attr_Transitions;Attr_Transitions=NULL;
            }
        }
    }
    
	char update_Gate(unsigned long long node_ID){
		snode::slinklist::slink *sc2;
		unsigned long long k_1;
		unsigned int k;
        char c;
		k_1=0;
		sc2=BRN->BN->node[node_ID].Li.first_link_pt;	
		// if(BRN->N>10) printf("Updating node %lld  inputs:",node_ID);
		k=0;
		while(sc2!=NULL){
			k++;
			if(s[sc2->sneighbor-1]=='1')
                k_1+=(unsigned long long)(pow(2.,(double)BRN->BN->node[node_ID].k_out-k));
			// if(BRN->N>10) printf("%lld:%d ",sc2->sneighbor,s[sc2->sneighbor-1]);
			sc2=sc2->next(sc2);//getchar();
		}
		// if(BRN->N>10) printf("k1=%lld\t returning %d\n",k_1,BRN->Gate[node_ID]->O[k_1]);
		if(k>0) {
            c=BRN->Gate[node_ID]->O[k_1]+48;
            return(c);    
        }
		else return(s[node_ID-1]);
	}
    
    void Asynchronously_update_all_Gates(){
        unsigned int i,l;
        double *x;
        
        x=new double[BRN->N+1];
        for(l=1;l<=BRN->N;l++) x[l]=l;
        randomize(x,BRN->N);  // order goes from 1 to N
        for(l=1;l<=BRN->N;l++){
            i=(unsigned int)(x[l]);  // order goes from 1 to N
            s[i-1]=update_Gate(i);    // index order goes from 0 to N-1 but gate index is 1 to N
        }
        s[BRN->N]=0;
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
                    s[Hits->A[l]-1]=48+Hits->B[l];  // updated node Hits->A[l] by force at the START; otheriwse in random order
                    skip->add_element(Hits->A[l]);  // updated alreay, skip during random order!
                }
            }
        }
        for(l=1;l<=BRN->N;l++){
            i=(unsigned int)(x[l]);  // order goes from 1 to N
            if(Hits!=NULL){
                if(skip->check_for_element(i)==0) s[i-1]=update_Gate(i);
            }
            else s[i-1]=update_Gate(i);    // index order goes from 0 to N-1 but gate index is 1 to N
        }
        s[BRN->N]=0;
        delete skip; skip = NULL;
        delete[] x;x=NULL;
    }
    
    void Asynchronously_update_all_Gates(int preserve_order,longint_ARRAY_pair *Hits, double p_hits[]){
        unsigned int i,l,a;
        double *x;
        
        x=new double[BRN->N+1];
        for(l=1;l<=BRN->N;l++) x[l]=l;
        randomize(x,BRN->N);  // order goes from 1 to N
        
    
        // 30: no U_Kinetochores at end
        if (preserve_order==30) {    // best so far version of preserve order!
            if (!strcmp(Model_type,"Life_and_Death")) {
                int Repl=BRN->MDAT->get_node_ID("Replication"); if (Repl==0) getchar();
                int cycA=BRN->MDAT->get_node_ID("CyclinA");     if (cycA==0) getchar();
                int f4N=BRN->MDAT->get_node_ID("f4N_DNA");       if (f4N==0) getchar();
                int FoxM1=BRN->MDAT->get_node_ID("FoxM1");      if (FoxM1==0) getchar();
                int cytokin_s=BRN->MDAT->get_node_ID("Ect2");     if (cytokin_s==0) getchar();
                int cytokin=BRN->MDAT->get_node_ID("Cytokinesis");     if (cytokin==0) getchar();
                int AK=BRN->MDAT->get_node_ID("A_Kinetochores");if (AK==0) getchar();
                int Cdh1=BRN->MDAT->get_node_ID("Cdh1");        if (Cdh1==0) getchar();
                int uK=BRN->MDAT->get_node_ID("U_Kinetochores");if (uK==0) getchar();
                int Cdc20=BRN->MDAT->get_node_ID("Cdc20");      if (Cdc20==0) getchar();
                int PreRC=BRN->MDAT->get_node_ID("Pre_RC");     if (PreRC==0) getchar();
                int cycB=BRN->MDAT->get_node_ID("CyclinB");     if (cycB==0) getchar();
                int Plk_H=BRN->MDAT->get_node_ID("Plk1_H");     if (Plk_H==0) getchar();
                int CycE=BRN->MDAT->get_node_ID("CyclinE");     if (CycE==0) getchar();
                
                int O_Rep=0,O_FoxM1=0,O_cytokin_s=0,O_cytokin=0,O_AK=0,O_4N=0,
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
                
                for(l=1;l<=BRN->N;l++){ if (x[l]==cytokin_s)  O_cytokin_s=l;}
                if (s[cytokin_s-1]==48)     { iii++; a=x[O_cytokin_s];  x[O_cytokin_s]=x[BRN->N-iii+1];  x[BRN->N-iii+1] = a;}
               
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
                    s[Hits->A[l]-1]=48+Hits->B[l];  // updated node Hits->A[l] by force at the START; otheriwse in random order
                    skip->add_element(Hits->A[l]);  // updated alreay, skip during random order!
                }
            }
        }
        for(l=1;l<=BRN->N;l++){
            i=(unsigned int)(x[l]);  // order goes from 1 to N
            if(Hits!=NULL){
                if(skip->check_for_element(i)==0) s[i-1]=update_Gate(i);
            }
            else s[i-1]=update_Gate(i);    // index order goes from 0 to N-1 but gate index is 1 to N
        }
        s[BRN->N]=0;
        delete skip; skip = NULL;
        delete[] x;x=NULL;
    }
    
    longint_ARRAY *get_slow_node_IDs(){
        long int idnow;
        longint_ARRAY *l;
        l=new longint_ARRAY();
        
       idnow=BRN->MDAT->get_node_ID("FoxO3"); if (idnow==0) printf("This model does not have FoxO3 node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("FoxO1"); if (idnow==0) printf("This model does not have FoxO1 node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("YAP"); if (idnow==0) printf("This model does not have YAP node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("WT1"); if (idnow==0) printf("This model does not have WT1 node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("TAZ"); if (idnow==0) printf("This model does not have TAZ node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("NfkB"); if (idnow==0) printf("This model does not have NfkB node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("c_Myb"); if (idnow==0) printf("This model does not have c_Myb node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("N_bcatenin"); if (idnow==0) printf("This model does not have N_bcatenin node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("SNAI1"); if (idnow==0) printf("This model does not have SNAI1 node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("LEF1"); if (idnow==0) printf("This model does not have LEF1 node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("SNAI2"); if (idnow==0) printf("This model does not have SNAI2 node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("Twist"); if (idnow==0) printf("This model does not have Twist node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("ZEB1_H"); if (idnow==0) printf("This model does not have ZEB1_H node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("ZEB1"); if (idnow==0) printf("This model does not have ZEB1 node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("N_bcatenin_H"); if (idnow==0) printf("This model does not have N_bcatenin_H node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("pRB"); if (idnow==0) printf("This model does not have pRB node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("Myc"); if (idnow==0) printf("This model does not have Myc node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("E2F1"); if (idnow==0) printf("This model does not have E2F1 node\n"); else l->add_element(idnow);
       idnow=BRN->MDAT->get_node_ID("FoxM1"); if (idnow==0) printf("This model does not have FoxM1 node\n"); else l->add_element(idnow);
   
     //   idnow=BRN->MDAT->get_node_ID("Replication"); if (idnow==0) printf("This model does not have Replication node\n"); else l->add_element(idnow);
        idnow=BRN->MDAT->get_node_ID("f4N_DNA"); if (idnow==0) printf("This model does not have f4N_DNA node\n"); else l->add_element(idnow);
     //   idnow=BRN->MDAT->get_node_ID("A_Kinetochores"); if (idnow==0) printf("This model does not have A_Kinetochores node\n"); else l->add_element(idnow);
     //   idnow=BRN->MDAT->get_node_ID("Cytokinesis"); if (idnow==0) printf("This model does not have Cytokinesis node\n"); else l->add_element(idnow);
        return (l);
    }
    
    void Asynchronously_update_all_Gates_Twoscale(longint_ARRAY_pair *Hits, double p_hits[]){
        unsigned int i,l, xp=0;
        int *x_fast,*x_slow;
        
        longint_ARRAY *slow_IDs;
        slow_IDs=get_slow_node_IDs();
        
        x_fast=new int[BRN->N-slow_IDs->N+1];

        for(l=1;l<=BRN->N;l++)
            if (slow_IDs->check_for_element(l)==0) {xp++; x_fast[xp]=l;} // adding the fast indexes in order
        randomize(x_fast,BRN->N-slow_IDs->N);  // order goes from 1 to N
       
        x_slow=new int[slow_IDs->N+1];
        for(l=1;l<=slow_IDs->N;l++) x_slow[l]=slow_IDs->A[l];  // adding the slow indexes in order
        randomize(x_slow,slow_IDs->N);  // order goes from 1 to N
       
        longint_ARRAY *skip;
        skip=new longint_ARRAY();
        double r;
        if(Hits!=NULL){
            r=rand_real(0,1);
            for(unsigned long int l=1;l<=Hits->N;l++){
                if(r<=p_hits[l]) {
                    s[Hits->A[l]-1]=48+Hits->B[l];  // updated node Hits->A[l] by force at the START; otheriwse in random order
                    skip->add_element(Hits->A[l]);  // updated alreay, skip during random order!
                }
            }
        }
        for(l=1;l<=BRN->N-slow_IDs->N;l++){
            i=(unsigned int)(x_fast[l]);
            if(Hits!=NULL){
                if(skip->check_for_element(i)==0) s[i-1]=update_Gate(i);
            }
            else s[i-1]=update_Gate(i);    // index order goes from 0 to N-1 but gate index is 1 to N
        }
        
        for(l=1;l<=slow_IDs->N;l++){
            i=(unsigned int)(x_slow[l]);
            if(rand_real(0,1)<=1/3.){
        //    if((rand_real(0,1)<=1/3.) || (s[i-1]==49)){
                if(Hits!=NULL){
                    if(skip->check_for_element(i)==0) s[i-1]=update_Gate(i);
                }
                else s[i-1]=update_Gate(i);    // index order goes from 0 to N-1 but gate index is 1 to N
            }
        }
        
        s[BRN->N]=0;
        delete skip; skip = NULL;
        delete[] x_fast;x_fast=NULL;
        delete[] x_slow;x_slow=NULL;
    }
    
	void synchronously_update_all_Gates(){
		char *s_n;
		unsigned int i;
		
		s_n=new char[BRN->N+1];
        for(i=0;i<BRN->N;i++)
            s_n[i]=update_Gate(i+1);
		for(i=0;i<BRN->N;i++) s[i]=s_n[i];
        s[BRN->N]=0;
		delete[] s_n;s_n=NULL;
	}
    
    void synchronously_update_all_Gates_print(){
        char *s_n;
        unsigned int i;
        
        s_n=new char[BRN->N+1];
        for(i=0;i<BRN->N;i++) {
            if(s[i]!=update_Gate(i+1)){
                BRN->print_gate(i+1);
                printf("%c -> %c\n",s[i],update_Gate(i+1));
                getchar();
            }
            s_n[i]=update_Gate(i+1);
        }
        for(i=0;i<BRN->N;i++) s[i]=s_n[i];
        s[BRN->N]=0;
        delete[] s_n;s_n=NULL;
    }
    
	short int synchronously_update_all_Gates_with_noise(double p_error){
		char *s_n;
		unsigned int i;
        short int ok;
		double p;
        
		s_n=new char[BRN->N+1];
		for(i=0;i<BRN->N;i++) s_n[i]=update_Gate(i+1);
		ok=0;
        for(i=0;i<BRN->N;i++){ 
            p=rand_real(0,1);
            if(p<=p_error) {
                if(s_n[i]=='1') s_n[i]='0'; else s_n[i]='1';
            }
            if(s_n[i]!=s[i]) {s[i]=s_n[i];ok=1;}
        }
        s[BRN->N]=0;
		delete[] s_n;s_n=NULL;
        return(ok);
	}
    short int synchronously_update_all_Gates_with_noise(double p_error[]){
		char *s_n;
		unsigned int i;
        short int ok;
		double p;
        
		s_n=new char[BRN->N+1];
		for(i=0;i<BRN->N;i++) s_n[i]=update_Gate(i+1);
		ok=0;
        for(i=0;i<BRN->N;i++){ 
            p=rand_real(0,1);
            if(p<=p_error[i+1]) {
                if(s_n[i]=='1') s_n[i]='0'; else s_n[i]='1';
            }
            if(s_n[i]!=s[i]) {s[i]=s_n[i];ok=1;}
        }
        s[BRN->N]=0;
		delete[] s_n;s_n=NULL;
        return(ok);
	}

	void set_state(const char s_n[]){
		unsigned int i;
		for(i=0;i<BRN->N;i++)
			s[i]=s_n[i];
        s[BRN->N]=0;
	}

    
    void set_random_network_state(){
        unsigned int i;
        double r;
		for(i=0;i<BRN->N;i++){
            r=rand_real(0,1);
        	if(r<=0.50000) s[BRN->N-i-1]='0';else s[BRN->N-i-1]='1';
         //   printf("i=%d\t r=%lg\ts[%d]=%c\n",i,r,BRN->N-i-1,s[BRN->N-i-1]);
        }
        s[BRN->N]=0;
      //  printf("Set random state s=%s\t s[0]='%c'",s,s[0]);getchar();
    }
    
    void set_random_network_state(longint_ARRAY *gin,unsigned long long int binary_init){
        unsigned int i;
        double r;
        char *s_in;
        
        s_in=decimal_to_string(binary_init,2,(int)gin->N);
        for(i=0;i<BRN->N;i++){
            r=rand_real(0,1);
            if(r<=0.50000) s[BRN->N-i-1]='0';else s[BRN->N-i-1]='1';
            //   printf("i=%d\t r=%lg\ts[%d]=%c\n",i,r,BRN->N-i-1,s[BRN->N-i-1]);
        }
        for(i=1;i<gin->N;i++){
            s[BRN->N-gin->A[i]-1]=s_in[i-1];
        }
        s[BRN->N]=0;
        //  printf("Set random state s=%s\t s[0]='%c'",s,s[0]);getchar();
    }
    
    void set_random_network_state(longint_ARRAY_pair *gl_s){
        unsigned int i;
        unsigned long int ginit;
        double r;
		for(i=0;i<BRN->N;i++){
            ginit=gl_s->check_for_element_A(BRN->N-i);
            if(ginit>0) s[BRN->N-i-1]=(unsigned short int)(gl_s->B[ginit])+48;
            else{
                r=rand_real(0,1);
                if(r<=0.50000) s[BRN->N-i-1]='0';else s[BRN->N-i-1]='1';
            }
        }
        s[BRN->N]=0;
    }
    
	void print_state(){
		unsigned int i;
		for(i=0;i<BRN->N;i++){ 
            printf("%c ",s[i]);
			if(BRN->MDAT!=NULL) printf("(%s)  ",BRN->MDAT->Node_name[i+1]);
		}
		printf("   %s\n",s);
	}
	
    void set_node_state(int n_id,int st_now){
            s[n_id-1]=48+st_now;
    }

    
    char *get_attractor_state_string(unsigned int  a,unsigned int  aindex){
        char *ss;
        ss=new char[BRN->N+1];
        for (int i=0; i<BRN->N; i++) {
            ss[i]=Attractor_valleys[a]->Basin_Stored[Attractor_valleys[a]->Attractor[aindex]]->s[i];
        }
        ss[BRN->N]=0;
        return(ss);
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
    
    unsigned long int check_if_state_has_Basin_ID(const char s_now[]){
        unsigned long int i;
        for(i=1;i<=Attractor_NR;i++){
            if(Attractor_valleys[i]->check_if_state_is_in_Basin(s_now)>0) return(i);
        }
        return(0);
    }
    
    unsigned long int check_if_state_was_sampled_enough(const char s_now[],unsigned long int vMin){
        unsigned long int i,w;
        for(i=1;i<=Attractor_NR;i++){
            if(Attractor_valleys[i]->check_if_state_is_in_Basin(s_now)>0){
                w=Attractor_valleys[i]->get_visits_to_state__Basin(s_now);
                if(w>=vMin) return(i);
                else return(0);
            }
        }
        return(0);
    }
    
    unsigned long int get_weight(const char s_now[]){
        unsigned long int i,w;
        for(i=1;i<=Attractor_NR;i++){
            if(Attractor_valleys[i]->check_if_state_is_in_Basin(s_now)>0){
                w=Attractor_valleys[i]->get_visits_to_state__Basin(s_now);
                 return(w);
            }
        }
        return(0);
    }

    int check_if_state_is_an_Attractor_State(const char s_now[]){
        for(int i=1;i<=Attractor_NR;i++){
            if(Attractor_valleys[i]->check_if_state_is_an_Attractor_State(s_now)>0) return(1);
        }
        return(0);
    }
    int check_if_state_is_an_Attractor_State_return_A(const char s_now[]){
        for(int i=1;i<=Attractor_NR;i++){
            if(Attractor_valleys[i]->check_if_state_is_an_Attractor_State(s_now)>0) return(i);
        }
        return(0);
    }
    
    double get_energy_of_state(const char s_now[],double NORM_factor,double p_error){
        unsigned long int i,w;
        double e;
        for(i=1;i<=Attractor_NR;i++){
            if(Attractor_valleys[i]->check_if_state_is_in_Basin(s_now)>0){
                w=Attractor_valleys[i]->get_visits_to_state__Basin(s_now);
                if(w>0) e=get_energy(w/NORM_factor,p_error);
                else  e=get_energy(0.1/NORM_factor,p_error);
                return(e);
            }
        }
        e=get_energy(0.1/NORM_factor,p_error);
        return(e);
    }
    
    Network_microstate_SYNC_const_noise new_MicroState(const char s_now[]){
        Network_microstate_SYNC_const_noise nms;
        unsigned int i;
        char *s_remember;
        
        nms.s=new char[BRN->N+1];
        sprintf(nms.s,"%s",s_now);
        nms.visits=0;
        nms.NextDoor=new neighbor_SYNC_const_noise_list[BRN->N+1];
        
        s_remember=new char[BRN->N+1];
        for(i=0;i<BRN->N;i++) s_remember[i]=s[i];
        s_remember[BRN->N]=0;
        
        //printf("before remember:");print_state();
        //printf("\tstate to remember: %s\n",s_remember);
        
        set_state(nms.s);
        synchronously_update_all_Gates();
        nms.NextDoor[0]=new neighbor_SYNC_const_noise(BRN->N,s);
        for(i=1;i<=BRN->N;i++){
            if(s[i]=='1') s[i]='0'; else s[i]='1'; // flip the ith node of the next state
            nms.NextDoor[i]=new neighbor_SYNC_const_noise(BRN->N,s);
            if(s[i]=='1') s[i]='0'; else s[i]='1'; // flip the ith node of the next state back
        }
        
        set_state(s_remember);
        delete[] s_remember;
        return(nms);    
    }
    
    Network_microstate_SYNC_const_noise new_MicroState_no_neighbors(const char s_now[]){
        Network_microstate_SYNC_const_noise nms;
        
        nms.s=new char[BRN->N+1];
        sprintf(nms.s,"%s",s_now);
        nms.visits=0;
        nms.NextDoor=NULL;
        return(nms);
    }

    
    unsigned long long add_trace_to_basin(unsigned long int b,Network_microstate_SYNC_const_noise_MAP_type tracedown,
                            const char s_now[]){
        Network_microstate_SYNC_const_noise_MAP_type::iterator It,Itsnow;
        State_index_MAP_type::const_iterator It_bs;
        unsigned long long sid;
        
        sid=0;
        
        if(s_now!=NULL) {
            Itsnow=tracedown.find(s_now);
            It_bs=Attractor_valleys[b]->States_in_basin.find(s_now);
            if(It_bs!=Attractor_valleys[b]->States_in_basin.end()) sid=It_bs->second;    
        } 
        for (It = tracedown.begin(); It != tracedown.end(); ++It) {
            It_bs=Attractor_valleys[b]->States_in_basin.find(It->second.s);
            if(It_bs==Attractor_valleys[b]->States_in_basin.end()){      
                if((s_now!=NULL)&&(It==Itsnow)) {
                    sid=add_to_Basin_Stored(b,It->second,1);
                }
                else add_to_Basin_Stored(b,It->second,0);
            }
            else {
                delete[] It->second.s;
            }
        }
        return(sid);
    }
    
    void add_to_latest_Attractor_cycle(unsigned long long sid){
        unsigned long long *Told;
        unsigned long int i;
        
        if(Attractor_valleys[Attractor_NR]->N_a==0){
            Attractor_valleys[Attractor_NR]->N_a=1;
            Attractor_valleys[Attractor_NR]->Attractor =new unsigned long long[Attractor_valleys[Attractor_NR]->N_a+1];
        }
        else{
            Told=Attractor_valleys[Attractor_NR]->Attractor;Attractor_valleys[Attractor_NR]->Attractor=NULL;
            Attractor_valleys[Attractor_NR]->N_a++;
            Attractor_valleys[Attractor_NR]->Attractor=new unsigned long long[Attractor_valleys[Attractor_NR]->N_a+1];
            for(i=1;i<Attractor_valleys[Attractor_NR]->N_a;i++)	
                Attractor_valleys[Attractor_NR]->Attractor[i]=Told[i];
            delete[] Told;Told=NULL;
        }
        
        //Attractor_valleys[Attractor_NR]->Attractor[Attractor_valleys[Attractor_NR]->N_a]=new Network_microstate_SYNC_const_noise;
        Attractor_valleys[Attractor_NR]->Attractor[Attractor_valleys[Attractor_NR]->N_a]=sid;
        
       // printf("Added state %s to attractor %ld into poz %ld\n",
         //      Attractor_valleys[Attractor_NR]->Basin_Stored[Attractor_valleys[Attractor_NR]->Attractor[Attractor_valleys[Attractor_NR]->N_a]]->s,
        //       Attractor_NR,Attractor_valleys[Attractor_NR]->N_a);
        //getchar();
    }
    
    
    unsigned long long add_to_Basin_Stored(unsigned long int b,Network_microstate_SYNC_const_noise MSnow,unsigned int visited){
        unsigned long long sind;
        double p,r;
        
        if(Attractor_valleys[b]->N_ms<N_MAX){
          //  printf("1: adding %s to basin %ld\n",MSnow.s,b);
            Attractor_valleys[b]->N_ms++;
            Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->N_ms]=new Network_microstate_SYNC_const_noise;
            *(Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->N_ms])=MSnow;
            //printf("Added state %s to attractor %ld\n",MSnow.s,Attractor_NR);
            Attractor_valleys[b]->States_in_basin[MSnow.s] = Attractor_valleys[b]->N_ms;
            return(Attractor_valleys[b]->N_ms);
        }
        else{
            if(visited==2){ // reserved for attractor states
                sind=Attractor_valleys[b]->get_RND_Element_of_rarest_states_HashMap();
                if(sind==0){printf("\n\nTROUBLE\n");exit(1);}
                
                Attractor_valleys[b]->States_in_basin.erase(Attractor_valleys[b]->Basin_Stored[sind]->s);
               
                 delete[] Attractor_valleys[b]->Basin_Stored[sind]->s;
                
                *(Attractor_valleys[b]->Basin_Stored[sind])=MSnow;
                Attractor_valleys[b]->States_in_basin[MSnow.s] = sind;
                return(sind);
            }
            else{
           // if(b!=1){printf("State number in attractor %lu will not increase!\n",b);getchar();}
                if(visited==0){
                    if (Attractor_valleys[b]->Rarest_visit_count>0) return(0);
                    else {
                       // printf("Before:\n");
                       // Attractor_valleys[b]->printf_attractor_full_basin_HashMap();
                        
                        //printf("2: adding %s to basin %ld\n",MSnow.s,b);
                        
                        sind=Attractor_valleys[b]->get_RND_Element_of_rarest_states_HashMap();
                        if(sind==0){printf("\n\nTROUBLE\n");exit(1);}
                        //printf("State to remove id %llu, %s\n",sind,Attractor_valleys[b]->Basin_Stored[sind]->s);
                        
                        Attractor_valleys[b]->States_in_basin.erase(Attractor_valleys[b]->Basin_Stored[sind]->s);
                        
                        //printf("After removing %llu %s from rearest list:\n",sind,Attractor_valleys[b]->Basin_Stored[sind]->s);
                       // Attractor_valleys[b]->printf_rarest_states_HashMap();
                        
                       // printf("Now inserting %s in its place:\n",MSnow.s);
                        
                        delete[] Attractor_valleys[b]->Basin_Stored[sind]->s;
                        
                        *(Attractor_valleys[b]->Basin_Stored[sind])=MSnow;
                        Attractor_valleys[b]->States_in_basin[MSnow.s] = sind;
                        
                       // printf("After:\n");
                       // Attractor_valleys[b]->printf_attractor_full_basin_HashMap();
                        
                       // Attractor_valleys[b]->printf_rarest_states_HashMap();
                        //getchar();
                        return(sind);
                    }
                }
                
                else{
                    p=pow(0.5,(double)Attractor_valleys[b]->Rarest_visit_count);
                    r=rand_real(0,1);
                    if(r<=p){
                        //printf("3: adding %s to basin %ld\n",MSnow.s,b);
                        //printf("Attemping to replace an old visited state with probab: %lg\n",p);
                        
                        sind=Attractor_valleys[b]->get_RND_Element_of_rarest_states_HashMap();
                        if(sind==0){printf("\n\nTROUBLE\n");exit(1);}
                        printf("State to remove id %llu, %s\n",sind,Attractor_valleys[b]->Basin_Stored[sind]->s);
                        
                        // Attractor_valleys[b]->printf_rarest_states_HashMap();
                        
                        Attractor_valleys[b]->States_in_basin.erase(Attractor_valleys[b]->Basin_Stored[sind]->s);

                        delete[] Attractor_valleys[b]->Basin_Stored[sind]->s;
                        
                        *(Attractor_valleys[b]->Basin_Stored[sind])=MSnow;
                        Attractor_valleys[b]->States_in_basin[MSnow.s] = sind;
                       // printf("State replaced\n");
                        //Attractor_valleys[b]->printf_rarest_states_HashMap();
                        //getchar();
                        return(sind);
                    }
                 else return(0);
                }
                }
        }
        return(0);
    }
    
    void print_cycling_components_of_attractor(unsigned long long b){
        unsigned int i,j,pr;
        char ok;
 
        if(!strcmp(BRN->name,"Tip_stalk_MARY_Retina")){
            printf("%s\t\t%c\n",BRN->MDAT->Node_name[2],Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[1]]->s[2-1]);
            printf("%s\t\t%c\n",BRN->MDAT->Node_name[10],Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[1]]->s[10-1]);
        }
        if(Attractor_valleys[b]->N_a>1){
            for(i=1;i<=BRN->N;i++){
                pr=0;
                ok=Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[1]]->s[i-1];
                for(j=2;j<=Attractor_valleys[b]->N_a;j++)
                    if(ok!=Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[j]]->s[i-1]) pr=1;
                if(pr){
                   for(j=1;j<=Attractor_valleys[b]->N_a;j++)
                        printf("%c ",Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[j]]->s[i-1]);
                   
                    printf("\t\t%s\n",BRN->MDAT->Node_name[i]);
                  
                }
            }   
        }
       
    }
    void create_new_Attractor_from_trace(Network_microstate_SYNC_const_noise_MAP_type tracedown){
        Attractor_landscape *old_AL;
        doublepointer *old_trans;
        unsigned long int i,b,j;
       // char *sneg;
        unsigned long long sid;
        State_index_MAP_type Atrace;
        State_index_MAP_type::const_iterator It,It_tr; 
        //Network_microstate_SYNC_const_noise_MAP_type::const_iterator It;
        Network_microstate_SYNC_const_noise MSnow;
        //Network_microstate_SYNC_const_noise_MAP_type Atrace;
        char *s_remember;
        
        s_remember=new char[BRN->N+1];
        for(i=0;i<=BRN->N;i++) s_remember[i]=s[i];
        s_remember[BRN->N]=0;
        
       // sneg=new char[BRN->N+1];
      //  for(i=0;i<=BRN->N;i++) sneg[i]='_';
      //  sneg[BRN->N]=0;
        
        if(Attractor_NR==0){
            Attractor_NR=1;
            Attractor_valleys=new Attractor_landscape[Attractor_NR+1];
            Trans_Basins= new doublepointer[Attractor_NR+1];
            for(j=0;j<=Attractor_NR;j++) {
                Trans_Basins[j]=new double[Attractor_NR+1];
                for(i=0;i<=Attractor_NR;i++) Trans_Basins[j][i]=0;
            }
        }
        else{
            old_AL=Attractor_valleys;Attractor_valleys=NULL;
            old_trans=Trans_Basins;Trans_Basins=NULL;
            Attractor_NR++;
            Attractor_valleys=new Attractor_landscape[Attractor_NR+1];
            for(i=1;i<Attractor_NR;i++)	Attractor_valleys[i]=old_AL[i];
           // delete[] old_AL;old_AL=NULL;
            
            Trans_Basins= new doublepointer[Attractor_NR+1];
            for(j=0;j<Attractor_NR;j++) {
                Trans_Basins[j]=new double[Attractor_NR+1];
                for(i=0;i<Attractor_NR;i++) Trans_Basins[j][i]=old_trans[j][i];
            }
            Trans_Basins[Attractor_NR]=new double[Attractor_NR+1];
            for(i=0;i<Attractor_NR;i++) {
                Trans_Basins[Attractor_NR][i]=0;   
                Trans_Basins[i][Attractor_NR]=0;   
            }
            for(i=0;i<Attractor_NR;i++) {delete[] old_trans[i];old_trans[i]=NULL;}
            delete[] old_trans;old_trans=NULL;
        }
        Attractor_valleys[Attractor_NR]=new Attractor_basin(Attractor_NR,N_MAX);
        
        Attractor_valleys[Attractor_NR]->Basin_ID=Attractor_NR;
        
        sid=add_trace_to_basin(Attractor_NR,tracedown,s);
        if(sid==0){
            printf("Index of an attractor state %s came back 0!\n",s);
            Attractor_valleys[Attractor_NR]->printf_attractor_full_basin_HashMap();
            getchar();
        }
        add_to_latest_Attractor_cycle(sid);
        Atrace[s]=sid;
        b=0;
        do{ //print_state();
            synchronously_update_all_Gates();
            if(Atrace.find(s) == Atrace.end()) {
                It_tr=Attractor_valleys[Attractor_NR]->States_in_basin.find(s);
                if(It_tr == Attractor_valleys[Attractor_NR]->States_in_basin.end()){
                    printf("Could not find attractor state \n%s\n in saved trace!\n",s);
                    MSnow=new_MicroState_no_neighbors(s);
                    sid=add_to_Basin_Stored(Attractor_NR,MSnow,2); 
                    
                    add_to_latest_Attractor_cycle(sid);
                    Atrace[s]=sid;
                    //getchar();
                }
                else  {add_to_latest_Attractor_cycle(It_tr->second);
                       Atrace[s]=It_tr->second;
                }
                //print_Network_microstate_SYNC_const_noise_HashMap(Atrace);
            }
            else b=1;
        }
        while(b==0);
        if(Attractor_valleys[Attractor_NR]->N_a>=N_MAX) N_MAX=Attractor_valleys[Attractor_NR]->N_a+2;
        
      //  Attractor_valleys[Attractor_NR]->printf_attractor_basin_short_info();
        Attractor_valleys[Attractor_NR]->printf_attractor_basin_info();
        //if(Attractor_valleys[Attractor_NR]->N_a>1) {
        //    print_cycling_components_of_attractor(Attractor_NR);
        //    getchar();
        //}
        // print_cycling_components_of_attractor(Attractor_NR);
        //Export_attractor_as_Cytoskape_node_attributes(Attractor_NR);
       // export_Attractor_onto_Boolean_network__PDF(Attractor_NR);
       // print_state();
        
        //getchar();
        
        //if(Attractor_NR/20==Attractor_NR/20.) {
        //    print_nodes_that_do_not_freeze_into_steady_state_in_any_attractor();
        //    getchar();
        //}
        set_state(s_remember);
        
        Atrace.clear();
        delete[] s_remember;s_remember=NULL;
    }
    
    
    double get_link_flux(const char s_1[],const char s_2[], int *ST_order,double p_error){
        // procedure needs to restore current state of system
        char *s_remember,*s_target;
        unsigned long int a1,a2;
        unsigned int h1,h2;
        unsigned long long w1,w2;
        double f;
        
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
        
        a1=check_if_state_has_Basin_ID(s_1);
        if(a1>0) w1=Attractor_valleys[a1]->get_visits_to_state__Basin(s_1); else w1=0;
        a2=check_if_state_has_Basin_ID(s_target);
        if(a2>0) w2=Attractor_valleys[a2]->get_visits_to_state__Basin(s_target); else w2=0;

        set_state(s_2);
        synchronously_update_all_Gates();
        h2=get_Hamming_distance_of_states(s_1,s);
        
        f= w1*pow(p_error,(double)h1)*pow(1.-p_error,(double)(BRN->N-h1)) - w2*pow(p_error,(double)h2)*pow(1.-p_error,(double)(BRN->N-h2));
        set_state(s_remember);
        
        delete[] s_remember;s_remember=NULL;
        delete[] s_target;s_target=NULL;
        return(f);
    }
 
    double get_link_flux_Zero_above_H_errors(const char s_1[],const char s_2[], int *ST_order,double p_error,double H_ERR){
        // procedure needs to restore current state of system
        char *s_remember,*s_target;
        unsigned long int a1,a2;
        unsigned int h1,h2;
        unsigned long long w1,w2;
        double f;
        
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
        if(h1<=H_ERR){
            a1=check_if_state_has_Basin_ID(s_1);
            if(a1>0) w1=Attractor_valleys[a1]->get_visits_to_state__Basin(s_1); else w1=0;
            a2=check_if_state_has_Basin_ID(s_target);
            if(a2>0) w2=Attractor_valleys[a2]->get_visits_to_state__Basin(s_target); else w2=0;
            
            set_state(s_2);
            synchronously_update_all_Gates();
            h2=get_Hamming_distance_of_states(s_1,s);
            
            f= w1*pow(p_error,(double)h1)*pow(1.-p_error,(double)(BRN->N-h1)) - w2*pow(p_error,(double)h2)*pow(1.-p_error,(double)(BRN->N-h2));
        }
        else f=0.;
        set_state(s_remember);
        delete[] s_remember;s_remember=NULL;
        delete[] s_target;s_target=NULL;
        return(f);
    }
    double get_link_P(const char s_1[],const char s_2[], int *ST_order,double p_error){
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
        f=pow(p_error,(double)h1)*pow(1.-p_error,(double)(BRN->N-h1));
        set_state(s_remember);
        delete[] s_remember;s_remember=NULL;
        delete[] s_target;s_target=NULL;
        return(f);
    }
    
    double get_link_P_Zero_above_H_errors(const char s_1[],const char s_2[], int *ST_order,double p_error,double H_ERR){
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
        if(h1<=H_ERR) f=pow(p_error,(double)h1)*pow(1.-p_error,(double)(BRN->N-h1));
        else f=0.;
        set_state(s_remember);
        delete[] s_remember;s_remember=NULL;
        delete[] s_target;s_target=NULL;
        return(f);
    }
    
    void print_nodes_that_do_not_freeze_into_steady_state_in_any_attractor(){
        unsigned int i,j,pr;
        unsigned long long b;
        char ok;
        Node_average= new double[BRN->N+1];
        Node_frozen_in=new unsigned long int[BRN->N+1];
        for(i=1;i<=BRN->N;i++) {
            Node_average[i]=0;
            Node_frozen_in[i]=0;
        }
        for(b=1;b<=Attractor_NR;b++){
            if(Attractor_valleys[b]->N_a>1){
                for(i=1;i<=BRN->N;i++){
                    pr=0;
                    ok=Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[1]]->s[i-1];
                    for(j=2;j<=Attractor_valleys[b]->N_a;j++){
                        if(ok!=Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[j]]->s[i-1]) pr=1;
                        Node_average[i]+=(Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[j]]->s[i-1]-48)/(double)Attractor_valleys[b]->N_a;
                    }
                    if(!pr) Node_frozen_in[i]++;
                }   
            }
            else for(i=1;i<=BRN->N;i++) {
                Node_frozen_in[i]++;
                Node_average[i]+=(Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[1]]->s[i-1]-48);
            }
        }
        for(i=1;i<=BRN->N;i++)
            printf("%s\t%ld\t%ld\t%lg\n",BRN->MDAT->Node_name[i],Node_frozen_in[i],Attractor_NR,Node_average[i]);
        delete[] Node_frozen_in; Node_frozen_in=NULL;
        delete[] Node_average; Node_average=NULL;
    }
    
   
    void print_ALL_attractor_info(){
        unsigned long int i,j;
        
        printf("\nTransition counts matrix:\n");
        for(j=1;j<=Attractor_NR;j++) {
            for(i=1;i<=Attractor_NR;i++) 
                printf("%lu -> %lu : %lg\n",j,i,Trans_Basins[j][i]);
        }
       // getchar();
        for(i=1;i<=Attractor_NR;i++){
            Attractor_valleys[i]->printf_attractor_basin_info();
            //Attractor_valleys[i]->printf_attractor_full_basin_HashMap();
           //  getchar();
        }
         
        //if(Attractor_NR>0) getchar();
    }
    
    double normalized_distance_between_attractor_pair(unsigned int aid1, unsigned int aid2){
        double h=0;
        for(int i=1;i<=Attractor_valleys[aid1]->N_a;i++){
            for(int j=1;j<=Attractor_valleys[aid2]->N_a;j++){
                h+=get_Hamming_distance_of_states(Attractor_valleys[aid1]->Basin_Stored[Attractor_valleys[aid1]->Attractor[i]]->s,Attractor_valleys[aid2]->Basin_Stored[Attractor_valleys[aid2]->Attractor[j]]->s);
            }
        }
        h/=(double)(Attractor_valleys[aid1]->N_a * Attractor_valleys[aid2]->N_a);
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
        printf("\t\t Coupled - N_attr=%lu\t H_av=%lg  ",Attractor_NR,h);
        h/=maximum_distance_between_p_Binary_string_pairs((int)Attractor_NR);
       // printf("N_attr=%lu\t max H= %lg\n",Attractor_NR,maximum_distance_between_p_Binary_string_pairs(Attractor_NR));getchar();
        printf(" ideal Hav = %lg  => nornalized H_av = %lg\n",maximum_distance_between_p_Binary_string_pairs((int)Attractor_NR),h);//getchar();
        return(h);
    }
    
    
    void Export_ALL_Attractors_onto_Boolean_Network(int Cytinfo){
        unsigned long int i;
        
        for(i=1;i<=Attractor_NR;i++){
            export_Attractor_onto_Boolean_network__PDF(i);
            if(Cytinfo) Export_attractor_as_Cytoskape_node_attributes(i);
        }
        Export_ALL_attractor_as_Cytoskape_node_attributes();
    }
    
    void Export_ALL_Attractors_onto_Boolean_Network_ALIVE(int Cytinfo){
        unsigned long int i;
        
        for(i=1;i<=Attractor_NR;i++){
            int id_dead=BRN->MDAT->get_node_ID("CAD");
            if (Attractor_valleys[i]->Basin_Stored[Attractor_valleys[i]->Attractor[1]]->s[id_dead-1]==48) {
                export_Attractor_onto_Boolean_network__PDF(i);
                if(Cytinfo) Export_attractor_as_Cytoskape_node_attributes(i);
            }
          }
    }
    
    void Export_ALL_Attractors_onto_Boolean_Network_ALIVE_or_cycle(int Cytinfo){
        unsigned long int i;
        
        for(i=1;i<=Attractor_NR;i++){
            if (Attractor_valleys[i]->N_a!=1)
                export_Attractor_onto_Boolean_network__PDF(i);
            if(Cytinfo) Export_attractor_as_Cytoskape_node_attributes(i);
            else {
                int id_dead=BRN->MDAT->get_node_ID("CAD");
                if (Attractor_valleys[i]->Basin_Stored[Attractor_valleys[i]->Attractor[1]]->s[id_dead-1]==48) {
                    export_Attractor_onto_Boolean_network__PDF(i);
                    if(Cytinfo) Export_attractor_as_Cytoskape_node_attributes(i);
                }
            }
        }
    }
    
    void Export_attractor_as_Cytoskape_node_attributes(unsigned long int b){
        unsigned long int i,k;
        FILE *f;
        char fn[300];
        
        
		for(i=1;i<=Attractor_valleys[b]->N_a;i++) {
            sprintf(fn,"%s/%s/Cytoskape/%s_Attr_%ld_C-%lu.txt", MODEL_Directory,DIRname,BRN->name,b,i);
            
            f=fopen(fn,"w");
            fprintf(f,"Attractor_%lu_C-%lu\tAttr_%lu_C-%lu\n",b,i,b,i);
            for(k=1;k<=BRN->N;k++){
                fprintf(f,"%s\t%c\n",BRN->MDAT->Node_name[k],Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s[k-1]);
            }
            fclose(f);
        }
    }
   
    void Export_ALL_attractor_as_Cytoskape_node_attributes(){
        unsigned long int i,k;
        FILE *f;
        char fn[300];
        
        
        sprintf(fn,"%s/%s/%s_Attr_ALL.txt", MODEL_Directory,DIRname, BRN->name);
        f=fopen(fn,"w");
        fprintf(f,"Node");
        for(unsigned long int b=1;b<=Attractor_NR;b++){
            for(i=1;i<=Attractor_valleys[b]->N_a;i++)
                fprintf(f,"\tAttr_%ld_C-%lu", b,i);
        }
        fprintf(f,"\n");
        
        for(k=1;k<=BRN->N;k++){
            fprintf(f,"%s",BRN->MDAT->Node_name[k]);
            for(unsigned long int b=1;b<=Attractor_NR;b++){
                for(i=1;i<=Attractor_valleys[b]->N_a;i++) {
                    fprintf(f,"\t%c",Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s[k-1]);
                }
            }
            fprintf(f,"\n");
        }
        fclose(f);
    }
    
    
    void export_GML(){
        FILE *f;
        char fn[300];
        std::string bla;
        snode::slinklist::slink *sc2;
        
        sprintf(fn,"%s/%s/%s_AutoGenerated_Hier_Graph.gml", MODEL_Directory,DIRname,BRN->name);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        fprintf(f, "Creator \"yFiles\"\nVersion \"2.17\"\ngraph\n[\n\thierarchic\t1\n\tlabel\t\"\"\n\tdirected\t1\n");
        
        for(int m=1;m<BRN->Module_Names.size();m++){
            fprintf(f, "\tnode\n\t[\n\t\tid\t%d\n\t\tlabel\t\"%s\"\n\t\tgraphics\n\t\t[\n\t\t\thasFill\t0\n\t\t]\n\t\tLabelGraphics\n\t\t[\n\t\tfontSize\t14\n\t\tfontStyle\t\"bold\"\n\t\tanchor\t\"tr\"\n\t\t\tfill\t    \"%s\"\n\t\t]\n\t\tisGroup\t1\n\t]\n",m,BRN->Module_Names[m].c_str(),
                    rgb2hex((int)(255*BRN->MDAT->r[BRN->Draw_in_Modules[m]->A[1]]),
                            (int)(255*BRN->MDAT->g[BRN->Draw_in_Modules[m]->A[1]]),
                            (int)(255*BRN->MDAT->b[BRN->Draw_in_Modules[m]->A[1]]),
                            true).c_str());
        }
       // rgb2hex(int r, int g, int b, bool with_head)
        
        for(int i=1;i<=BRN->N;i++){
            fprintf(f, "\tnode\n\t[\n\t\tid\t%lu\n\t\tlabel\t\"%s\"\n\t\tgraphics\n\t\t[\n\t\tx\t%lg\n\t\ty\t%lg\n\t\t\tw\t%lf\n\t\t\th\t30.0\n\t\t\tfill\t\"%s\"\n\t\t]\n\t\tLabelGraphics\n\t\t[\n\t\tfontSize\t12\n\t\t]\n\t\tgid\t%d\n\t]\n",BRN->Module_Names.size()+i,BRN->MDAT->Node_name[i],BRN->MDAT->x[i],BRN->MDAT->y[i],
                    maximum(30,30*strlen(BRN->MDAT->Node_name[i])/4.),
                    rgb2hex((int)(255*BRN->MDAT->r[i]),(int)(255*BRN->MDAT->g[i]),(int)(255*BRN->MDAT->b[i]),true).c_str(),
                    BRN->MDAT->Module_ID[i]);
        }
        for(int s=1;s<=BRN->N;s++){
            sc2=BRN->BN->node[s].Li.first_link_pt;
            int k=0;
            while(sc2!=NULL){
                k++;
                int linksign=BRN->Gate[s]->check_overall_sign_of_input(k);
                switch (linksign) {
                    case 1: bla = "standard"; break;
                    case -1: bla = "t_shape"; break;
                    case 0: bla = "diamond"; break;
                    default: bla = "diamond"; break;
                }
                fprintf(f, "\tedge\n\t[\n\t\tsource\t%llu\n\t\ttarget\t%lu\n\t\tgraphics\n\t\t[\n\t\t\tfill\t\"#000000\"\n\t\ttargetArrow\t\"%s\"\n\t\t] \n\t]",
                        BRN->Module_Names.size()+sc2->sneighbor,
                        BRN->Module_Names.size()+s,bla.c_str()
                );
                sc2=sc2->next(sc2);
            }
        }
        fprintf(f, "]\n");
        fclose(f);
    }
    
    void export_subgraph_to_GML(sgraph *gcc, const char *subdir){
        FILE *f;
        char fn[300];
        std::string bla;
        snode::slinklist::slink *sc2,*sc_og;
        
        sprintf(fn,"%s/%s/_EXP/%s/%s.gml", MODEL_Directory,BRN->path,subdir,gcc->gname);
        f=fopen(fn,"w");
        fprintf(f, "Creator \"yFiles\"\nVersion \"2.17\"\ngraph\n[\n\tlabel\t\"\"\n\tdirected\t1\n");
        
                       
        for(unsigned long int i=1;i<=gcc->N;i++){
           // printf("gcc node %ld has node ID %ld\t",i,gcc->node[i].ID);// getchar();
           // printf("%s\n",PD->BRN->MDAT->Node_name[gcc->node[i].ID]); getchar();
            fprintf(f, "\tnode\n\t[\n\t\tid\t%lu\n\t\tlabel\t\"%s\"\n\t\tgraphics\n\t\t[\n\t\tx\t%lg\n\t\ty\t%lg\n\t\t\tw\t%lf\n\t\t\th\t30.0\n\t\t\tfill\t\"%s\"\n\t\t]\n\t\tLabelGraphics\n\t\t[\n\t\tfontSize\t12\n\t\t]\n\n\t]\n",i,
                    BRN->MDAT->Node_name[gcc->node[i].ID],
                    BRN->MDAT->x[gcc->node[i].ID],BRN->MDAT->y[gcc->node[i].ID],
                    maximum(30,30*strlen(BRN->MDAT->Node_name[gcc->node[i].ID])/4.),
                    rgb2hex((int)(255*BRN->MDAT->r[gcc->node[i].ID]),(int)(255*BRN->MDAT->g[gcc->node[i].ID]),(int)(255*BRN->MDAT->b[gcc->node[i].ID]),true).c_str());
        }
        for(unsigned long int s=1;s<=gcc->N;s++){
            sc2=gcc->node[s].Li.first_link_pt;
            while(sc2!=NULL){
                int k=0;
                sc_og=BRN->BN->node[gcc->node[s].ID].Li.first_link_pt;
                int linksign=0;
                while(sc_og!=NULL){
                    k++;
                    if (sc_og->sneighbor == gcc->node[sc2->sneighbor].ID)
                        linksign=BRN->Gate[gcc->node[s].ID]->check_overall_sign_of_input(k);
                    sc_og=sc_og->next(sc_og);
                }
                switch (linksign) {
                    case 1: bla = "standard"; break;
                    case -1: bla = "t_shape"; break;
                    case 0: bla = "diamond"; break;
                    default: bla = "diamond"; break;
                }
                fprintf(f, "\tedge\n\t[\n\t\tsource\t%llu\n\t\ttarget\t%lu\n\t\tgraphics\n\t\t[\n\t\t\tfill\t\"#000000\"\n\t\ttargetArrow\t\"%s\"\n\t\t] \n\t]",
                        sc2->sneighbor,s,bla.c_str());
                sc2=sc2->next(sc2);
            }
        }
        fprintf(f, "]\n");
        fclose(f);
    }
    
    void export_Attractor_onto_Boolean_network__PDF(unsigned long int b){
        unsigned long int i;
        FILE *f;
        char fn[300];
        double XMAX,YMAX,XMIN, YMIN,e_min,e_max,ch_min,ch_max,RMAX=15;
        unsigned long int s;
        snode::slinklist::slink *sc2;
        
        sprintf(fn,"mkdir %s/%s/Cytoskape\n",MODEL_Directory_system,DIRname);
        system(fn);
        
        for(i=1;i<=Attractor_valleys[b]->N_a;i++) {
            sprintf(fn,"%s/%s/Cytoskape/%s_Attr_%ld_C-%lu.eps", MODEL_Directory,DIRname,BRN->name,b,i);
            printf("%s\n",fn);
            
            f=fopen(fn,"w");

            COLOR_SCHEME=GYR;
            setcolors_PAJEK();
            
            XMAX=maximum(BRN->N, BRN->MDAT->x);
            XMIN=minimum(BRN->N, BRN->MDAT->x);
            YMAX=maximum(BRN->N, BRN->MDAT->y);
            YMIN=minimum(BRN->N, BRN->MDAT->y);
            //printf("XMIN=%lg\tXMAX=%lg\nYMIN=%lg\tYMAX=%lg\n",XMIN,XMAX,YMIN,YMAX);
            //getchar();
            EpsInit(f,-25,-25,XMAX-XMIN+25,YMAX-YMIN+25);
            e_min=2000000;e_max=0;
            ch_min=2000000;ch_max=-2000000;
            
            EpsSetRgb(f,0.2,0.2,0.2);
            EpsSetLinewidth(f, 1);
            EpsSetFont(f,"Times-Roman",RMAX);
            int linksign;
            for(s=1;s<=BRN->N;s++){
                int k=0;
                sc2=BRN->BN->node[s].Li.first_link_pt;
                while(sc2!=NULL){
                    k++;
                   // printf("Link %s (%d) -> %s (%llu) (input # %d)\n",BRN->MDAT->Node_name[s],s, BRN->MDAT->Node_name[sc2->sneighbor],sc2->sneighbor,k);
                    linksign=BRN->Gate[s]->check_overall_sign_of_input(k);
                    switch (linksign) {
                        case 1:  EpsSetRgb(f,0,0,0.8);break;
                        case -1: EpsSetRgb(f,0.8,0,0);break;
                        default: EpsSetRgb(f,0.2,0.2,0.2); break;
                    }
                    EpsDrawLine(f,BRN->MDAT->x[s]-XMIN,YMAX-BRN->MDAT->y[s],BRN->MDAT->x[sc2->sneighbor]-XMIN,YMAX-BRN->MDAT->y[sc2->sneighbor]);
                    EpsFillCircle(f,BRN->MDAT->x[s]-XMIN - (BRN->MDAT->x[s] - BRN->MDAT->x[sc2->sneighbor])/4.,
                                    YMAX-BRN->MDAT->y[s] + (BRN->MDAT->y[s] - BRN->MDAT->y[sc2->sneighbor])/4.,3);
                    sc2=sc2->next(sc2);
                }
            }
            
            //all nodes with energy color, size and attracror border
            
            EpsSetLinewidth(f,3);
            
            for(s=1;s<=BRN->N;s++){
                switch(Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s[s-1]){
              //      case 48:  EpsSetRgb(f,0,1,0);break; //0
              //      case 49:  EpsSetRgb(f,1,0,0);break; //1
                      case 48:  EpsSetRgb(f,0.5,0.5,1);break; //0
                      case 49:  EpsSetRgb(f,1,0.5,0);break; //1
                    default:  EpsSetRgb(f,1,1,1);break;
                }
                //printf("node %s\t color %d\n",BRN->MDAT->Node_name[s],Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s[s-1]);
                
                EpsFillCircle(f,BRN->MDAT->x[s]-XMIN,YMAX-BRN->MDAT->y[s],RMAX);
                
                if (BRN->MDAT->Module_ID!=NULL) {
                    switch(BRN->MDAT->Module_ID[s]){
                        case 0:  EpsSetRgb(f,0,0,0);break;
                        case 1:  EpsSetRgb(f,0.8,0.8,0);break; // dark yellow
                        case 2: EpsSetRgb(f,0.8,0,0.8);break;  // violet
                        case 3: EpsSetRgb(f,0.5,0.5,0.5);break;  // gray
                        case 4: EpsSetRgb(f,0.8,0,0);break;     // dark red
                        case 5: EpsSetRgb(f,0,0.8,0.8);break;
                        case 6: EpsSetRgb(f,0,0,0.8);break; // blue
                        case 7: EpsSetRgb(f,0.5,0.5,0.8);break; // 
                        default: EpsSetRgb(f,0,0,0);break;
                    }
                }
                EpsDrawCircle(f,BRN->MDAT->x[s]-XMIN,YMAX-BRN->MDAT->y[s],RMAX);
                EpsSetRgb(f,0,0,0);
                EpsDrawString_Centered(f, BRN->MDAT->x[s]-XMIN,YMAX-BRN->MDAT->y[s], BRN->MDAT->Node_name[s]);
            }
            
            sprintf(fn,"%s/%s/Cytoskape/%s_Attr_%ld_C-%lu", MODEL_Directory_system,DIRname,BRN->name,b,i);
            close_EPSDRAW_file(f,fn,2);
        }
      //  printf("\ndone\n");
    }
    
    void export_Attractor_onto_Boolean_network__PDF(unsigned long int b, char fn_med[]){
        unsigned long int i;
        FILE *f;
        char fn[600];
        double XMAX,YMAX,XMIN, YMIN,e_min,e_max,ch_min,ch_max,RMAX=15;
        unsigned long int s;
        snode::slinklist::slink *sc2;
        
        for(i=1;i<=Attractor_valleys[b]->N_a;i++) {
            sprintf(fn,"%s/%s/Cytoskape/%s/%s_Attr_%ld_C-%lu.eps", MODEL_Directory,DIRname,fn_med,BRN->name,b,i);
          //  printf("%s\n",fn);
            
            f=fopen(fn,"w");
            
            COLOR_SCHEME=GYR;
            setcolors_PAJEK();
            
            XMAX=maximum(BRN->N, BRN->MDAT->x);
            XMIN=minimum(BRN->N, BRN->MDAT->x);
            YMAX=maximum(BRN->N, BRN->MDAT->y);
            YMIN=minimum(BRN->N, BRN->MDAT->y);
            //printf("XMIN=%lg\tXMAX=%lg\nYMIN=%lg\tYMAX=%lg\n",XMIN,XMAX,YMIN,YMAX);
            //getchar();
            EpsInit(f,-25,-25,XMAX-XMIN+25,YMAX-YMIN+25);
            e_min=2000000;e_max=0;
            ch_min=2000000;ch_max=-2000000;
            
            EpsSetRgb(f,0.2,0.2,0.2);
            EpsSetLinewidth(f, 1);
            EpsSetFont(f,"Times-Roman",RMAX);
            int linksign;
            for(s=1;s<=BRN->N;s++){
                int k=0;
                sc2=BRN->BN->node[s].Li.first_link_pt;
                while(sc2!=NULL){
                    k++;
                    // printf("Link %s (%d) -> %s (%llu) (input # %d)\n",BRN->MDAT->Node_name[s],s, BRN->MDAT->Node_name[sc2->sneighbor],sc2->sneighbor,k);
                    linksign=BRN->Gate[s]->check_overall_sign_of_input(k);
                    switch (linksign) {
                        case 1:  EpsSetRgb(f,0,0,0.8);break;
                        case -1: EpsSetRgb(f,0.8,0,0);break;
                        default: EpsSetRgb(f,0.2,0.2,0.2); break;
                    }
                    EpsDrawLine(f,BRN->MDAT->x[s]-XMIN,YMAX-BRN->MDAT->y[s],BRN->MDAT->x[sc2->sneighbor]-XMIN,YMAX-BRN->MDAT->y[sc2->sneighbor]);
                    EpsFillCircle(f,BRN->MDAT->x[s]-XMIN - (BRN->MDAT->x[s] - BRN->MDAT->x[sc2->sneighbor])/4.,
                                  YMAX-BRN->MDAT->y[s] + (BRN->MDAT->y[s] - BRN->MDAT->y[sc2->sneighbor])/4.,3);
                    sc2=sc2->next(sc2);
                }
            }
            
            //all nodes with energy color, size and attracror border
            
            EpsSetLinewidth(f,3);
            
            for(s=1;s<=BRN->N;s++){
                switch(Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s[s-1]){
                        //      case 48:  EpsSetRgb(f,0,1,0);break; //0
                        //      case 49:  EpsSetRgb(f,1,0,0);break; //1
                    case 48:  EpsSetRgb(f,0.5,0.5,1);break; //0
                    case 49:  EpsSetRgb(f,1,0.5,0);break; //1
                    default:  EpsSetRgb(f,1,1,1);break;
                }
                //printf("node %s\t color %d\n",BRN->MDAT->Node_name[s],Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s[s-1]);
                
                EpsFillCircle(f,BRN->MDAT->x[s]-XMIN,YMAX-BRN->MDAT->y[s],RMAX);
                
                if (BRN->MDAT->Module_ID!=NULL) {
                    switch(BRN->MDAT->Module_ID[s]){
                        case 0:  EpsSetRgb(f,0,0,0);break;
                        case 1:  EpsSetRgb(f,0.8,0.8,0);break; // dark yellow
                        case 2: EpsSetRgb(f,0.8,0,0.8);break;  // violet
                        case 3: EpsSetRgb(f,0.5,0.5,0.5);break;  // gray
                        case 4: EpsSetRgb(f,0.8,0,0);break;     // dark red
                        case 5: EpsSetRgb(f,0,0.8,0.8);break;
                        case 6: EpsSetRgb(f,0,0,0.8);break; // blue
                        case 7: EpsSetRgb(f,0.5,0.5,0.8);break; //
                        default: EpsSetRgb(f,0,0,0);break;
                    }
                }
                EpsDrawCircle(f,BRN->MDAT->x[s]-XMIN,YMAX-BRN->MDAT->y[s],RMAX);
                EpsSetRgb(f,0,0,0);
                EpsDrawString_Centered(f, BRN->MDAT->x[s]-XMIN,YMAX-BRN->MDAT->y[s], BRN->MDAT->Node_name[s]);
            }
            
            //sprintf(fn,"%s/%s/Cytoskape/%s/%s_Attr_%ld_C-%lu.eps", MODEL_Directory,DIRname,fn_med,BRN->name,b,i);
            sprintf(fn,"%s/%s/Cytoskape/%s/%s_Attr_%ld_C-%lu", MODEL_Directory_system,DIRname,fn_med,BRN->name,b,i);
            close_EPSDRAW_file(f,fn,2);
        }
        //  printf("\ndone\n");
    }
    
    void print_Network_microstate_SYNC_const_noise_HashMap(Network_microstate_SYNC_const_noise_MAP_type tracedown){
        Network_microstate_SYNC_const_noise_MAP_type::const_iterator It;
        unsigned int n,j;
        n=0;
        for (It = tracedown.begin(); It != tracedown.end(); ++It) {
            n++;
            printf("%d %s  v=%llu \n",n,
                   It->second.s,It->second.visits);
            //printf("\tstrcmp key (in struct %s) with current state %s is %d\n",It->second.s,s,strcmp(It->second.s,s));
            for(j=0;j<=BRN->N;j++){
                if(j>0) printf("\t");
               // printf("\t%d %s\n",j,It->second.NextDoor[j]->s);
            }
        }
        //getchar();
    }
    
    
    void print_microstate(Network_microstate_SYNC_const_noise MsSnow){
        printf("%s  v=%llu \n",MsSnow.s,MsSnow.visits);
    }
    
    unsigned long int find_attractor_basin_for_state(const char s_now[],int addforsure){
		//unsigned long long i,j,j_old,k,ind1;
		unsigned long int b,i;
        Network_microstate_SYNC_const_noise_MAP_type tracedown;
        Network_microstate_SYNC_const_noise MSnow;
        Network_microstate_SYNC_const_noise_MAP_type::iterator It;
        char *s_remember;
        
        b=check_if_state_has_Basin_ID(s_now);
        if(b>0) return(b);
        
        MSnow=new_MicroState_no_neighbors(s_now);
        tracedown[MSnow.s]=MSnow;
        
        //print_Network_microstate_SYNC_const_noise_HashMap(tracedown);
        
        s_remember=new char[BRN->N+1];
        for(i=0;i<=BRN->N;i++) s_remember[i]=s[i];
        
        set_state(s_now);
        
        do{ //print_ALL_attractor_info();
            
            synchronously_update_all_Gates();
            //print_Network_microstate_SYNC_const_noise_HashMap(tracedown);
            b=check_if_state_has_Basin_ID(s);
            if(b>0){
                if(addforsure) add_trace_to_basin(b,tracedown,s_now);
                else add_trace_to_basin(b,tracedown,NULL);
            }
            else {
                //It=tracedown.find(s);
                if(tracedown.find(s) != tracedown.end()) {
                   // printf("Starting the creation of a new attractor %ld!\n",Attractor_NR+1);//getchar();
                    create_new_Attractor_from_trace(tracedown);
                    b=Attractor_NR;
                }              
                else {
                    MSnow=new_MicroState_no_neighbors(s);
                    tracedown[MSnow.s]=MSnow;
                    //print_Network_microstate_SYNC_const_noise_HashMap(tracedown);
                }
            }
        }
        while(b==0);
        set_state(s_remember);
        //empty_map_of_Network_microstate_SYNC_const_noise(tracedown); 
      
        delete[] s_remember;s_remember=NULL;
        
       // tracedown.clear();
        //printf("Attractor basin ");
        return(b);
    }
    
    unsigned long int find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(const char s_now[]){
		//unsigned long long i,j,j_old,k,ind1;
		unsigned long int b,i;
        Network_microstate_SYNC_const_noise_MAP_type tracedown;
        Network_microstate_SYNC_const_noise MSnow;
        Network_microstate_SYNC_const_noise_MAP_type::iterator It;
        char *s_remember;
        
        b=check_if_state_has_Basin_ID(s_now);
        if(b>0) {//printf("Returing b=%ld\n",b);
            return(b);
        }
        
        MSnow=new_MicroState_no_neighbors(s_now);
        tracedown[MSnow.s]=MSnow;
        
        //print_Network_microstate_SYNC_const_noise_HashMap(tracedown);
        
        s_remember=new char[BRN->N+1];
        for(i=0;i<=BRN->N;i++) s_remember[i]=s[i];
        
        set_state(s_now);
        
        do{ //print_ALL_attractor_info();
            
            synchronously_update_all_Gates();
            //print_Network_microstate_SYNC_const_noise_HashMap(tracedown);
            b=check_if_state_has_Basin_ID(s);
            if(b>0){
                set_state(s_remember);
                delete[] s_remember;s_remember=NULL;
                
                for (It = tracedown.begin(); It != tracedown.end(); ++It) {
                    delete[] It->second.s;It->second.s=NULL;
                }
                tracedown.clear();
                return(b);
            }
            else {
                //It=tracedown.find(s);
                if(tracedown.find(s) != tracedown.end()) {
                    // printf("Starting the creation of a new attractor %ld!\n",Attractor_NR+1);//getchar();
                    create_new_Attractor_from_trace(tracedown);
                    b=Attractor_NR;
                    printf("Returing new attractor b=%ld\n",b);
                    set_state(s_remember);
                    return(b);
                }
                else {
                    MSnow=new_MicroState_no_neighbors(s);
                    tracedown[MSnow.s]=MSnow;
                    //print_Network_microstate_SYNC_const_noise_HashMap(tracedown);
                }
            }
        }
        while(b==0);
        set_state(s_remember);
        //empty_map_of_Network_microstate_SYNC_const_noise(tracedown);
        delete[] MSnow.s;MSnow.s=NULL;
        
        delete[] s_remember;s_remember=NULL;
        
        //printf("Returing b=%ld\n",b);
        return(b);
    }
    
    unsigned long int find_position_state_hits_attractor_cycle_without_altering_tracking(const char s_now[],int C_AID){
		unsigned long int b,i;
        Network_microstate_SYNC_const_noise_MAP_type tracedown;
        Network_microstate_SYNC_const_noise MSnow;
        Network_microstate_SYNC_const_noise_MAP_type::iterator It;
        char *s_remember;
       
        b=Attractor_valleys[C_AID]->check_if_state_is_an_Attractor_State(s_now);
        if(b>0){ return(b);}
        
        MSnow=new_MicroState_no_neighbors(s_now);
        tracedown[MSnow.s]=MSnow;
        
        //print_Network_microstate_SYNC_const_noise_HashMap(tracedown);
        
        s_remember=new char[BRN->N+1];
        for(i=0;i<=BRN->N;i++) s_remember[i]=s[i];
        
        set_state(s_now);
       // printf("\t\t\t\t%s\n",s_now);
        do{ synchronously_update_all_Gates();
          //  printf("\t\t\t\t%s\n",s);
            b=Attractor_valleys[C_AID]->check_if_state_is_an_Attractor_State(s);
            
            if(b>0){
                delete[] s_remember;s_remember=NULL;
                delete[] MSnow.s;MSnow.s=NULL;
                
                for (It = tracedown.begin(); It != tracedown.end(); ++It) {
                    if(It->second.s!=NULL) {
                       // delete[] It->second.s;It->second.s=NULL;
                        It->second.NextDoor=NULL;
                    }
                }
                tracedown.clear();
                return(b);
            }
            else {
                //It=tracedown.find(s);
                if(tracedown.find(s) != tracedown.end()) {
                    printf("State \n_%s_-->down to_%s_\nis not in the marked attractor C_AID=%d !!!\n",s_now,s,C_AID);
                    Attractor_valleys[C_AID]->printf_attractor_basin_info();
                   // return(0);
                    getchar();
                }
                else {
                    delete[] MSnow.s;MSnow.s=NULL;
                    MSnow=new_MicroState_no_neighbors(s);
                    tracedown[MSnow.s]=MSnow;
                }
            }
        }
        while(b==0);
        set_state(s_remember);
        
        delete[] s_remember;s_remember=NULL;
        delete[] MSnow.s;MSnow.s=NULL;
        
        for (It = tracedown.begin(); It != tracedown.end(); ++It) {
            if(It->second.s!=NULL) {
                delete[] It->second.s;It->second.s=NULL;
                It->second.NextDoor=NULL;
            }
        }
        tracedown.clear();
        return(b);
    }
    
    double *w_attractor_share_of_cycle_states(int C_AID){
        double *w;
        unsigned long int b,nt;
        
        w=new double[Attractor_valleys[C_AID]->N_a+1];
        nt=0;
        for (int i=0; i<=Attractor_valleys[C_AID]->N_a; i++) w[i]=0;
        
        for (long int ms=1; ms<=Attractor_valleys[C_AID]->N_ms; ms++) {
            //if(Attractor_valleys[C_AID]->N_a>1){
            
            b=find_position_state_hits_attractor_cycle_without_altering_tracking(Attractor_valleys[C_AID]->Basin_Stored[ms]->s,C_AID);
            
            w[b]+=Attractor_valleys[C_AID]->Basin_Stored[ms]->visits;
           // }
           // else {w[1]+=Attractor_valleys[C_AID]->Basin_Stored[ms]->visits;}
            nt+=Attractor_valleys[C_AID]->Basin_Stored[ms]->visits;
        }
        for (int i=1; i<=Attractor_valleys[C_AID]->N_a; i++) {
            w[i]/=(double)nt;
            //if(Attractor_valleys[C_AID]->N_a>1) {printf("Share of Cycle state %d = %lg\n",i,w[i]);getchar();}
        }
        w[0]=nt;
       // printf("W_0=%lg\n",w[0]);//getchar();
        return(w);
    }
    
    
    longint_ARRAY *LIST_attractor_portion_landing_on_cycle_state(int C_AID,int cycle_state_index){
        longint_ARRAY *P;
        unsigned long int b;
        
        
        P=new longint_ARRAY();
        for (long int ms=1; ms<=Attractor_valleys[C_AID]->N_ms; ms++) {
            //if(Attractor_valleys[C_AID]->N_a>1){
            b=find_position_state_hits_attractor_cycle_without_altering_tracking(Attractor_valleys[C_AID]->Basin_Stored[ms]->s,C_AID);
            if(b==cycle_state_index) P->add_element(ms);
        }
        return(P);
    }

    
    int check_if_current_state_is_direct_neighbour_of_Network_microstate_SYNC_const_noise(Network_microstate_SYNC_const_noise *MSnow){
        int i;
        for(i=0;i<=BRN->N;i++)
            if(!strcmp(s,MSnow->NextDoor[i]->s)) return(i);
        return(-1);
    }
    
    
    void record_time_trace_from_current_state_SYNC_const_noise(unsigned long int N_trace,double p_error){
        unsigned long int b,ii,j,c,inplace;
        Network_microstate_SYNC_const_noise *MSnow;
        double alpha;
        unsigned long long sid;
        short int newst;
        //int MSneigh;
        char *double_state;
        State_index_MAP_type::iterator IT_tr;
        
        double_state=new char[BRN->N*2+2];
        alpha=1.-p_error*BRN->N;
        // trans basin probab is NOT normalized to timesteps in this routine, so multiple timetraces may be piled on top of each other
        MSnow=NULL;
       // int DELnow=0;
        
        for(ii=1;ii<=N_trace;ii++){
           // printf("Going from %s: ",s); 
            b=find_attractor_basin_for_state(s,1);
            sid=Attractor_valleys[b]->get_Microstate_and_mark_visit(s);
            W_total_steps++;
          
                MSnow=new Network_microstate_SYNC_const_noise;
                *MSnow=new_MicroState(s);  // this one lists the neighbors!!!
           //     DELnow=1;
         //   }
            //printf("ii=%lu: Attractor state for %s is %ld, stored id=%lld\n",ii,s,b,sid);
            //getchar();
           // print_microstate(*MSnow);       
            //getchar();
            
            inplace=1;
            do{ 
                strcpy(double_state,"");
                sprintf(double_state, "%s",s);
               
                newst=synchronously_update_all_Gates_with_noise(p_error);
                sprintf(double_state, "%s %s",double_state,s);
                //IT_tr=Transition_track.find(double_state);
                //if(IT_tr==Transition_track.end())                    
                //    Transition_track[double_state]=1;
                //else IT_tr->second++;
                if(newst==0) {
                    inplace++;
                    ii++;
                    //printf("Staying put\n");//getchar();
                }
               // else printf("Going: %s  -> %s\n",MSnow->s,s);
            }
            while((newst==0)&&(ii<N_trace));
            
            //MSneigh=check_if_current_state_is_direct_neighbour_of_Network_microstate_SYNC_const_noise(MSnow);
            //if(MSneigh>-1) MSnow->NextDoor[MSneigh]->out_to+=inplace;
            
            if(inplace>1) {
                Attractor_valleys[b]->mark_visit(MSnow->s,inplace-1);
                W_total_steps+=inplace-1;
            }
            c=find_attractor_basin_for_state(MSnow->NextDoor[0]->s,0);
            
            Trans_Basins[b][c]+=inplace*alpha;
            for(j=1;j<=BRN->N;j++) {
                c=find_attractor_basin_for_state(MSnow->NextDoor[j]->s,0);
                Trans_Basins[b][c]+=inplace*p_error; 
                //printf("\t\t\tj=%lu Adding %lg to Trans_M[%lu][%lu]\n",j,inplace*p_error,b,c);
            }
            
            //print_ALL_attractor_info();
           // printf("%ld steps done (%ld attractors now)\n",ii,Attractor_NR); 
            //getchar();
            delete[] MSnow->s;MSnow->s=NULL;
            for(int i=0;i<=BRN->N;i++) {
                delete[] MSnow->NextDoor[i]->s;MSnow->NextDoor[i]->s=NULL;
                delete MSnow->NextDoor[i];MSnow->NextDoor[i]=NULL;
            }
            delete[] MSnow->NextDoor;MSnow->NextDoor=NULL;
            delete MSnow;MSnow=NULL;
        }
        delete[] double_state;double_state=NULL;
    }
    
    void record_time_trace_from_current_state_SYNC_node_specific_noise(unsigned long int N_trace,double p_error[]){
        unsigned long int b,ii,j,c,inplace;
        Network_microstate_SYNC_const_noise *MSnow;
        double alpha;
        unsigned long long sid;
        short int newst;
        //int MSneigh;
        char *double_state;
        State_index_MAP_type::iterator IT_tr;
        
        double_state=new char[BRN->N*2+2];
        alpha=1.;
        for(j=1;j<=BRN->N;j++)
            alpha-=p_error[j];
        // trans basin probab is NOT normalized to timesteps in this routine, so multiple timetraces may be piled on top of each other
        MSnow=NULL;
        
        for(ii=1;ii<=N_trace;ii++){
            // printf("Going from %s: ",s); 
            b=find_attractor_basin_for_state(s,1);
            sid=Attractor_valleys[b]->get_Microstate_and_mark_visit(s);
            W_total_steps++;
            //if(sid>0) {
            //    MSnow=Attractor_valleys[b]->Basin_Stored[sid];
            //    DELnow=0;
           // }
           // else {
            //    if(DELnow==1){
                    delete[] MSnow->s;MSnow->s=NULL;
                    for(int i=0;i<=BRN->N;i++) {
                        delete[] MSnow->NextDoor[i]->s;MSnow->NextDoor[i]->s=NULL;
                        delete MSnow->NextDoor[i];MSnow->NextDoor[i]=NULL;
                    }
                    delete[] MSnow->NextDoor;MSnow->NextDoor=NULL;
                    delete MSnow;MSnow=NULL;
             //   }
                MSnow=new Network_microstate_SYNC_const_noise;
                *MSnow=new_MicroState(s);
              //  DELnow=1;
           // }
            // printf("ii=%lu: Attractor state for %s is %ld, stored id=%lld\n",ii,s,b,sid); 
            //getchar();
            // print_microstate(*MSnow);       
            //getchar();
            
            inplace=1;
            do{ 
                strcpy(double_state,"");
                sprintf(double_state, "%s",s);
                
                newst=synchronously_update_all_Gates_with_noise(p_error);
                sprintf(double_state, "%s %s",double_state,s);
                //IT_tr=Transition_track.find(double_state);
                //if(IT_tr==Transition_track.end())                    
                //    Transition_track[double_state]=1;
                //else IT_tr->second++;
                if(newst==0) {
                    inplace++;
                    ii++;
                    //printf("Staying put\n");//getchar();
                }
                // else printf("Going: %s  -> %s\n",MSnow->s,s);
            }
            while((newst==0)&&(ii<N_trace));
            
           // MSneigh=check_if_current_state_is_direct_neighbour_of_Network_microstate_SYNC_const_noise(MSnow);
            //if(MSneigh>-1) MSnow->NextDoor[MSneigh]->out_to+=inplace;
            
            if(inplace>1) {
                Attractor_valleys[b]->mark_visit(MSnow->s,inplace-1);
                W_total_steps+=inplace-1;
            }
            c=find_attractor_basin_for_state(MSnow->NextDoor[0]->s,0);
            
            Trans_Basins[b][c]+=inplace*alpha;
            for(j=1;j<=BRN->N;j++) {
                c=find_attractor_basin_for_state(MSnow->NextDoor[j]->s,0);
                Trans_Basins[b][c]+=inplace*p_error[j]; 
                //printf("\t\t\tj=%lu Adding %lg to Trans_M[%lu][%lu]\n",j,inplace*p_error,b,c);
            }
            
            //print_ALL_attractor_info();
            // printf("%ld steps done (%ld attractors now)\n",ii,Attractor_NR); 
            //getchar();
        }
        delete[] double_state;double_state=NULL;
    }
    
    
    void record_timetraces_from_collection_of_random_inputs(unsigned long int N_trace,unsigned int N_RND,double p_error){
        unsigned int ii,r_done;
        r_done=import_state_of_tracking(N_trace,p_error);
        for(ii=r_done+1;ii<=N_RND;ii++){
            init_random(ii);
            set_random_network_state();
            record_time_trace_from_current_state_SYNC_const_noise(N_trace,p_error);
            export_state_of_tracking(N_trace,ii,p_error);//
            printf("Random seed %d\n",ii);
            //if(ii/5==ii/5.) 
                //exit(1);
        }
        //export_landscape(N_trace,N_RND,p_error);
    }
    
    void record_timetraces_from_collection_of_random_inputs_per_inputBOX(unsigned long int N_trace,unsigned int N_RND,double p_error, longint_ARRAY *gl_fixedIN){
        unsigned int ii,r_done;
        r_done=import_state_of_tracking(N_trace,p_error);
        if (r_done<N_RND*pow(2,gl_fixedIN->N)) {
        // track phase space in the rest of input-boxes
            for (int box=0; box<pow(2,gl_fixedIN->N); box++) {
                for(ii=1;ii<=N_RND;ii++){
                    init_random(ii+box*N_RND);
                    set_random_network_state(gl_fixedIN,box);
                    record_time_trace_from_current_state_SYNC_const_noise(N_trace,p_error);
                   // printf("Box %d of %ld -- seed %d -- total random seed %d\n",box,(long int)pow(2, gl_fixedIN->N),ii,ii+box*N_RND);
                }
                export_state_of_tracking(N_trace,box*N_RND,p_error);//
                printf("Box %d of %ld -- seed %d \n",box,(long int)pow(2, gl_fixedIN->N),ii);
            }
         //   export_landscape(N_trace,N_RND,p_error);
        }
    }
    
    void record_timetraces_from_collection_of_random_inputs(unsigned long int N_trace,unsigned int N_RND,double p_error[]){
        unsigned int ii,r_done;
        r_done=import_state_of_tracking(N_trace);
        for(ii=r_done+1;ii<=N_RND;ii++){
            init_random(ii);
            set_random_network_state();
            record_time_trace_from_current_state_SYNC_node_specific_noise(N_trace,p_error);
            export_state_of_tracking(N_trace,ii,p_error);//
            printf("Random seed %d\n",ii);
           // if(ii/5==ii/5.) exit(1);
        }
        //export_landscape_p_error_b(N_trace,N_RND,p_error[0]);
    }
    
    void export_state_of_tracking(unsigned long int N_trace,unsigned int N_RND,double p_error){
        unsigned long int i,j,b;
        FILE *f;
        char fn[300];
        
        sprintf(fn,"%s/%s/%s_Landscape__Ntr-%lu_Nmax-%llu_p-%.3lf",MODEL_Directory,DIRname,SAMPLEname,N_trace,N_MAX,p_error);
        if(BRN->silenced_inputs->N>0)
            for (unsigned k=1; k<=BRN->silenced_inputs->N; k++) {
                sprintf(fn,"%s_%s",fn,BRN->MDAT->Node_name[BRN->silenced_inputs->A[k]]);
            }
        sprintf(fn,"%s.dat",fn);

        f=fopen(fn,"w");
        fprintf(f,"%ld\tN_RND %u\n",Attractor_NR,N_RND);
        for(j=1;j<=Attractor_NR;j++) {
            for(i=1;i<=Attractor_NR;i++) fprintf(f,"%lu %lu %lg\n",j,i,Trans_Basins[j][i]);
        }
        fprintf(f,"\n");
        for(b=1;b<=Attractor_NR;b++){
            fprintf(f,"%lu\t%lu\t%lu\n\t",Attractor_valleys[b]->Basin_ID,Attractor_valleys[b]->N_a,Attractor_valleys[b]->Rarest_visit_count);
              //  printf("A-%lu\t%lu\t%lu\n\t",Attractor_valleys[b]->Basin_ID,Attractor_valleys[b]->N_a,Attractor_valleys[b]->Rarest_visit_count);
            fprintf(f,"%llu\t%llu\t%llu\n",Attractor_valleys[b]->N_ms,Attractor_valleys[b]->W_Total,Attractor_valleys[b]->N_MAX);
              //  printf("B-%llu\t%llu\t%llu\n",Attractor_valleys[b]->N_ms,Attractor_valleys[b]->W_Total,Attractor_valleys[b]->N_MAX);
            for(i=1;i<=Attractor_valleys[b]->N_a;i++){
                fprintf(f,"\t%lu %s   %llu\n",i,Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s,
                                          Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->visits);
                   // printf("C-\t%lu %s   %llu\n",i,Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s,
                    //    Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->visits);
            }
            for(i=1;i<=Attractor_valleys[b]->N_ms;i++){
                if(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)==0){
                    fprintf(f,"\t\t%s   %llu\n",Attractor_valleys[b]->Basin_Stored[i]->s,
                            Attractor_valleys[b]->Basin_Stored[i]->visits);
                       // printf("D-\ti=%d\t%s   %llu\n",i,Attractor_valleys[b]->Basin_Stored[i]->s,
                       //     Attractor_valleys[b]->Basin_Stored[i]->visits);
                   // getchar();
                }
            }
        }
        fclose(f);//getchar();
    }
    
    void export_state_of_tracking(unsigned long int N_trace,unsigned int N_RND,double p_error[]){
        unsigned long int i,j,b;
        FILE *f;
        char fn[300];
        
        
        sprintf(fn,"%s/%s/%s_Landscape__Ntr-%lu_Nmax-%llu_pnode.dat",MODEL_Directory,DIRname,SAMPLEname,N_trace,N_MAX);
        sprintf(fn,"%s/%s/%s_Landscape__Ntr-%lu_Nmax-%llu_pnode", MODEL_Directory,DIRname,SAMPLEname,N_trace,N_MAX);
        if(BRN->silenced_inputs->N>0)
            for (unsigned k=1; k<=BRN->silenced_inputs->N; k++) {
                sprintf(fn,"%s_%s",fn,BRN->MDAT->Node_name[BRN->silenced_inputs->A[k]]);
            }
        sprintf(fn,"%s.dat",fn);

        f=fopen(fn,"w");
        fprintf(f,"%ld\tN_RND %u\n",Attractor_NR,N_RND);
        for(j=1;j<=Attractor_NR;j++) {
            for(i=1;i<=Attractor_NR;i++) fprintf(f,"%lu %lu %lg\n",j,i,Trans_Basins[j][i]);
        }
        fprintf(f,"\n");
        for(i=1;i<=BRN->N;i++) fprintf(f,"%ld\t%lg\n",i,p_error[i]);
        fprintf(f,"\n");
        for(b=1;b<=Attractor_NR;b++){
            fprintf(f,"%lu\t%lu\t%lu\n\t",Attractor_valleys[b]->Basin_ID,Attractor_valleys[b]->N_a,Attractor_valleys[b]->Rarest_visit_count);
            fprintf(f,"%llu\t%llu\t%llu\n",Attractor_valleys[b]->N_ms,Attractor_valleys[b]->W_Total,Attractor_valleys[b]->N_MAX);
            for(i=1;i<=Attractor_valleys[b]->N_a;i++){
                fprintf(f,"\t%lu %s   %llu\n",i,Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s,
                        Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->visits);
            }
            for(i=1;i<=Attractor_valleys[b]->N_ms;i++){
                if(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)==0)
                    fprintf(f,"\t\t%s   %llu\n",Attractor_valleys[b]->Basin_Stored[i]->s,
                            Attractor_valleys[b]->Basin_Stored[i]->visits);
            }
        }
        fclose(f);
    }
    
    unsigned int import_state_of_tracking(unsigned long int N_trace,double p_error){
        unsigned long int i,j,b,na,N_a_in_attr,rv;
        FILE *f;
        char fn[300],bla[300],*s_now;
        unsigned int r_done;
        unsigned long long int Ms,v,nmax;
        Network_microstate_SYNC_const_noise MSnow;
        
        s_now=new char[BRN->N+1];
        
        sprintf(fn,"%s/%s/%s_Landscape__Ntr-%lu_Nmax-%llu_p-%.3lf",MODEL_Directory,DIRname,SAMPLEname,N_trace,N_MAX,p_error);
        if(BRN->silenced_inputs->N>0)
            for (unsigned k=1; k<=BRN->silenced_inputs->N; k++) {
                sprintf(fn,"%s_%s",fn,BRN->MDAT->Node_name[BRN->silenced_inputs->A[k]]);
            }
        sprintf(fn,"%s.dat",fn);
        
        printf("Importing %s\n",fn);//getchar();
        f=fopen(fn,"r");
        if(f==NULL) return(0);
        fscanf(f,"%lu %s %u\n",&na,bla,&r_done);
        Attractor_NR=na;
        Attractor_valleys=new Attractor_landscape[Attractor_NR+1];
            
        Trans_Basins= new doublepointer[Attractor_NR+1];
        for(j=0;j<=Attractor_NR;j++) {
            Trans_Basins[j]=new double[Attractor_NR+1];
            for(i=0;i<=Attractor_NR;i++) Trans_Basins[j][i]=0;
        }
        for(j=1;j<=Attractor_NR;j++) {
            for(i=1;i<=Attractor_NR;i++) fscanf(f,"%lu %lu %lg\n",&b,&b,&(Trans_Basins[j][i]));
        }
            
        
        for(b=1;b<=na;b++){
            fscanf(f,"%lu\t%lu\t%lu\n\t",&i,&N_a_in_attr,&rv);
            fscanf(f,"%llu\t%llu\t%llu\n",&Ms,&v,&nmax);
            
            Attractor_valleys[b]=new Attractor_basin(b,nmax);
            Attractor_valleys[b]->Basin_ID=b;
            Attractor_valleys[b]->N_a=N_a_in_attr;
            Attractor_valleys[b]->Attractor=new unsigned long long[Attractor_valleys[b]->N_a+1];
            Attractor_valleys[b]->N_ms=Ms;
            Attractor_valleys[b]->W_Total=v;
            Attractor_valleys[b]->Rarest_visit_count=rv;
            W_total_steps+=Attractor_valleys[b]->W_Total;
            
            for(i=1;i<=Attractor_valleys[b]->N_a;i++)	{
                fscanf(f,"\t%lu %s   %llu\n",&j,s_now,&v);
                s_now[BRN->N]=0;
                MSnow=new_MicroState_no_neighbors(s_now);
               
                Attractor_valleys[b]->Attractor[i]=i;
                Attractor_valleys[b]->Basin_Stored[i]=new Network_microstate_SYNC_const_noise;
                *(Attractor_valleys[b]->Basin_Stored[i])=MSnow;
                Attractor_valleys[b]->Basin_Stored[i]->visits=v;
                Attractor_valleys[b]->States_in_basin[MSnow.s] = i;
            }
            for(i=Attractor_valleys[b]->N_a+1;i<=Attractor_valleys[b]->N_ms;i++){
                fscanf(f,"\t%s   %llu\n",s_now,&v);
                s_now[BRN->N]=0;
                MSnow=new_MicroState_no_neighbors(s_now);
                
                Attractor_valleys[b]->Basin_Stored[i]=new Network_microstate_SYNC_const_noise;
                *(Attractor_valleys[b]->Basin_Stored[i])=MSnow;
                Attractor_valleys[b]->Basin_Stored[i]->visits=v;
                Attractor_valleys[b]->States_in_basin[MSnow.s] = i;
            } 
        }
        fclose(f);
        delete[] s_now;
        //print_ALL_attractor_info();//getchar();
        return r_done;
    }
    
    unsigned int import_state_of_tracking(unsigned long int N_trace){
        unsigned long int i,j,b,na,N_a_in_attr,rv;
        FILE *f;
        char fn[300],bla[300],*s_now;
        unsigned int r_done;
        unsigned long long int Ms,v,nmax;
        Network_microstate_SYNC_const_noise MSnow;
        double x;
        
        s_now=new char[BRN->N+1];
        
        sprintf(fn,"%s/%s/%s_Landscape__Ntr-%lu_Nmax-%llu_pnode",MODEL_Directory,DIRname,SAMPLEname,N_trace,N_MAX);
        if(BRN->silenced_inputs->N>0)
            for (unsigned k=1; k<=BRN->silenced_inputs->N; k++) {
                sprintf(fn,"%s_%s",fn,BRN->MDAT->Node_name[BRN->silenced_inputs->A[k]]);
            }
        sprintf(fn,"%s.dat",fn);
        
        f=fopen(fn,"r");
        if(f==NULL) return(0);
        fscanf(f,"%lu %s %u\n",&na,bla,&r_done);
        Attractor_NR=na;
        Attractor_valleys=new Attractor_landscape[Attractor_NR+1];
        
        Trans_Basins= new doublepointer[Attractor_NR+1];
        for(j=0;j<=Attractor_NR;j++) {
            Trans_Basins[j]=new double[Attractor_NR+1];
            for(i=0;i<=Attractor_NR;i++) Trans_Basins[j][i]=0;
        }
        for(j=1;j<=Attractor_NR;j++) {
            for(i=1;i<=Attractor_NR;i++) fscanf(f,"%lu %lu %lg\n",&b,&b,&(Trans_Basins[j][i]));
        }
        for(i=1;i<=BRN->N;i++) {
            fscanf(f,"%ld\t%lg\n",&b,&x);
        }
        
        for(b=1;b<=na;b++){
            fscanf(f,"%lu\t%lu\t%lu\n\t",&i,&N_a_in_attr,&rv);
            fscanf(f,"%llu\t%llu\t%llu\n",&Ms,&v,&nmax);
            
            Attractor_valleys[b]=new Attractor_basin(b,nmax);
            Attractor_valleys[b]->Basin_ID=b;
            Attractor_valleys[b]->N_a=N_a_in_attr;
            Attractor_valleys[b]->Attractor=new unsigned long long[Attractor_valleys[b]->N_a+1];
            Attractor_valleys[b]->N_ms=Ms;
            Attractor_valleys[b]->W_Total=v;
            Attractor_valleys[b]->Rarest_visit_count=rv;
            W_total_steps+=Attractor_valleys[b]->W_Total;
            
            for(i=1;i<=Attractor_valleys[b]->N_a;i++)	{
                fscanf(f,"\t%lu %s   %llu\n",&j,s_now,&v);
                s_now[BRN->N]=0;
                MSnow=new_MicroState_no_neighbors(s_now);
                
                Attractor_valleys[b]->Attractor[i]=i;
                Attractor_valleys[b]->Basin_Stored[i]=new Network_microstate_SYNC_const_noise;
                *(Attractor_valleys[b]->Basin_Stored[i])=MSnow;
                Attractor_valleys[b]->Basin_Stored[i]->visits=v;
                Attractor_valleys[b]->States_in_basin[MSnow.s] = i;
            }
            for(i=Attractor_valleys[b]->N_a+1;i<=Attractor_valleys[b]->N_ms;i++){
                fscanf(f,"\t%s   %llu\n",s_now,&v);
                s_now[BRN->N]=0;
                MSnow=new_MicroState_no_neighbors(s_now);
                
                Attractor_valleys[b]->Basin_Stored[i]=new Network_microstate_SYNC_const_noise;
                *(Attractor_valleys[b]->Basin_Stored[i])=MSnow;
                Attractor_valleys[b]->Basin_Stored[i]->visits=v;
                Attractor_valleys[b]->States_in_basin[MSnow.s] = i;
            } 
        }
        fclose(f);
        return r_done;
    }
    
    void export_landscape(unsigned long int N_trace,unsigned int N_RND,double p_error){
        //N_MAX
        unsigned long int i,j,k,b,a;
       // State_index_MAP_type::iterator It_tr;
        FILE *f;
        char fn[300];//,*s1,*s2;
        //int ok1,ok2;
        unsigned int h,V_MIN;
        
        V_MIN=2;
        sprintf(fn,"mkdir %s/%s/Cytoskape/Landscape\n",MODEL_Directory_system,DIRname);
        
        system(fn);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_p-%.3lf.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX,p_error);
        f=fopen(fn,"w");
        fprintf(f,"State_st\tTo\tState_end\tState_Transition_p\tLinkP\n");
        for(b=1;b<=Attractor_NR;b++)
            for(a=1;a<=Attractor_NR;a++)
                for(i=1;i<=Attractor_valleys[b]->N_ms;i++){
                    if((Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                        for(j=1;j<=Attractor_valleys[a]->N_ms;j++)
                            if((Attractor_valleys[a]->Basin_Stored[j]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(j)>0)){
                                set_state(Attractor_valleys[b]->Basin_Stored[i]->s);
                                synchronously_update_all_Gates();   
                                h=get_Hamming_distance_of_states(s,Attractor_valleys[a]->Basin_Stored[j]->s);
                                if(h<=2){
                                    if(Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)
                                        fprintf(f,"%s\t->\t%s\t%lg\t%lg\n",Attractor_valleys[b]->Basin_Stored[i]->s,Attractor_valleys[a]->Basin_Stored[j]->s,-log(Attractor_valleys[b]->Basin_Stored[i]->visits*pow(p_error, (double)h)*pow(1-p_error, (double)BRN->N-h)/(double)((N_trace+1)*N_RND)),-log(pow(p_error, (double)h)*pow(1-p_error, (double)BRN->N-h)/(double)((N_trace+1)*N_RND)));
                                    else  fprintf(f,"%s\t->\t%s\t%lg\t%lg\n",Attractor_valleys[b]->Basin_Stored[i]->s,Attractor_valleys[a]->Basin_Stored[j]->s,-log(0.01*pow(p_error, (double)h)*pow(1-p_error, (double)BRN->N-h)/(double)((N_trace+1)*N_RND)),-log(pow(p_error, (double)h)*pow(1-p_error, (double)BRN->N-h)/(double)((N_trace+1)*N_RND)));
                                }
                            }
                    }
                }
        fclose(f);
        
         sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_p-%.3lf__BasinID.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX,p_error);
        f=fopen(fn,"w");
        fprintf(f,"Basin_ID\n");
        for(b=1;b<=Attractor_NR;b++){
            for(i=1;i<=Attractor_valleys[b]->N_ms;i++)
                if((Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                    fprintf(f,"%s = %lu\n",Attractor_valleys[b]->Basin_Stored[i]->s,
                                        Attractor_valleys[b]->Basin_ID);
            }
        }
        fclose(f);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_p-%.3lf__Probability.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX,p_error);
        f=fopen(fn,"w");
        fprintf(f,"Energy\n");
        for(b=1;b<=Attractor_NR;b++){
            for(i=1;i<=Attractor_valleys[b]->N_ms;i++)
                if((Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                    fprintf(f,"%s = %lg\n",Attractor_valleys[b]->Basin_Stored[i]->s,
                        -log(Attractor_valleys[b]->Basin_Stored[i]->visits/(double)((N_trace+1)*N_RND)));
            }
        }
        fclose(f);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_p-%.3lf__Attractor_State.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX,p_error);
        f=fopen(fn,"w");
        fprintf(f,"Attractor_state\n");
        for(b=1;b<=Attractor_NR;b++){
            for(i=1;i<=Attractor_valleys[b]->N_a;i++){
                fprintf(f,"%s = %lu\n",Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s,
                        Attractor_valleys[b]->Basin_ID);
           }
        }
        fclose(f);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_p-%.3lf__GenesON.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX,p_error);
        f=fopen(fn,"w");
        fprintf(f,"Genes_ON\n");
        for(b=1;b<=Attractor_NR;b++){
            for(i=1;i<=Attractor_valleys[b]->N_ms;i++)
                if((Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                    fprintf(f,"%s = ",Attractor_valleys[b]->Basin_Stored[i]->s);
                    for(k=1;k<=BRN->N;k++)
                        if(Attractor_valleys[b]->Basin_Stored[i]->s[k-1]!=48) fprintf(f,"%s ",BRN->MDAT->Node_name[k]);
                    fprintf(f,"\n");
                }
        }
        fclose(f);

    }
    void export_landscape_p_error_b(unsigned long int N_trace,unsigned int N_RND,double p_error_b){
        //N_MAX
        unsigned long int i,j,k,b,a;
        // State_index_MAP_type::iterator It_tr;
        FILE *f;
        char fn[300];//,*s1,*s2;
        //int ok1,ok2;
        unsigned int h,V_MIN;
        
        V_MIN=2;
        sprintf(fn,"mkdir %s/%s/Cytoskape/Landscape\n",         MODEL_Directory_system,DIRname);
        system(fn);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_pnode.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX);
        f=fopen(fn,"w");
        fprintf(f,"State_st\tTo\tState_end\tState_Transition_p\tLinkP\n");
        for(b=1;b<=Attractor_NR;b++)
            for(a=1;a<=Attractor_NR;a++)
                for(i=1;i<=Attractor_valleys[b]->N_ms;i++){
                    if((Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                        for(j=1;j<=Attractor_valleys[a]->N_ms;j++)
                            if((Attractor_valleys[a]->Basin_Stored[j]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(j)>0)){
                                set_state(Attractor_valleys[b]->Basin_Stored[i]->s);
                                synchronously_update_all_Gates();   
                                h=get_Hamming_distance_of_states(s,Attractor_valleys[a]->Basin_Stored[j]->s);
                                if(h<=2){
                                    if(Attractor_valleys[b]->Basin_Stored[i]->visits>=1)
                                        fprintf(f,"%s\t->\t%s\t%lg\t%lg\n",Attractor_valleys[b]->Basin_Stored[i]->s,Attractor_valleys[a]->Basin_Stored[j]->s,-log(Attractor_valleys[b]->Basin_Stored[i]->visits*pow(p_error_b, (double)h)*pow((double)(1-p_error_b), (double)BRN->N-h)/(double)((N_trace+1)*N_RND)),-log(pow(p_error_b, (double)h)*pow((double)(1-p_error_b), (double)BRN->N-h)/(double)((N_trace+1)*N_RND)));
                                    else  fprintf(f,"%s\t->\t%s\t%lg\t%lg\n",Attractor_valleys[b]->Basin_Stored[i]->s,Attractor_valleys[a]->Basin_Stored[j]->s,-log(0.01*pow(p_error_b, (double)h)*pow((double)(1-p_error_b), (double)BRN->N-h)/(double)((N_trace+1)*N_RND)),-log(pow(p_error_b, (double)h)*pow((double)(1-p_error_b), (double)BRN->N-h)/(double)((N_trace+1)*N_RND)));
                                }
                            }
                    }
                }
        fclose(f);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_pnode__BasinID.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX);
        f=fopen(fn,"w");
        fprintf(f,"Basin_ID\n");
        for(b=1;b<=Attractor_NR;b++){
            for(i=1;i<=Attractor_valleys[b]->N_ms;i++)
                if((Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                    fprintf(f,"%s = %lu\n",Attractor_valleys[b]->Basin_Stored[i]->s,
                            Attractor_valleys[b]->Basin_ID);
                }
        }
        fclose(f);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_pnode__Probability.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX);
        f=fopen(fn,"w");
        fprintf(f,"Energy\n");
        for(b=1;b<=Attractor_NR;b++){
            for(i=1;i<=Attractor_valleys[b]->N_ms;i++)
                if((Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                    fprintf(f,"%s = %lg\n",Attractor_valleys[b]->Basin_Stored[i]->s,
                            -log(Attractor_valleys[b]->Basin_Stored[i]->visits/(double)((N_trace+1)*N_RND)));
                }
        }
        fclose(f);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_pnode__Attractor_State.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX);
        f=fopen(fn,"w");
        fprintf(f,"Attractor_state\n");
        for(b=1;b<=Attractor_NR;b++){
            for(i=1;i<=Attractor_valleys[b]->N_a;i++){
                fprintf(f,"%s = %lu\n",Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s,
                        Attractor_valleys[b]->Basin_ID);
            }
        }
        fclose(f);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_Nmax-%llu_pnode__GenesON.txt\n",
                MODEL_Directory,DIRname,BRN->name,N_trace,N_RND,N_MAX);
        f=fopen(fn,"w");
        fprintf(f,"Genes_ON\n");
        for(b=1;b<=Attractor_NR;b++){
            for(i=1;i<=Attractor_valleys[b]->N_ms;i++)
                if((Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)||(Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                    fprintf(f,"%s = ",Attractor_valleys[b]->Basin_Stored[i]->s);
                    for(k=1;k<=BRN->N;k++)
                        if(Attractor_valleys[b]->Basin_Stored[i]->s[k-1]!=48) fprintf(f,"%s ",BRN->MDAT->Node_name[k]);
                    fprintf(f,"\n");
                }
        }
        fclose(f);
        
    }
   
    unsigned short int cascading_change( unsigned long int a,  unsigned long int inp, unsigned long int outp){
        if (Attractor_valleys[a]->N_a==1) return(0);
        for(unsigned long int k=1;k<=Attractor_valleys[a]->N_a;k++)
            if(   ((k<Attractor_valleys[a]->N_a) &&
                  (Attractor_valleys[a]->Basin_Stored[Attractor_valleys[a]->Attractor[k]]->s[inp-1] != Attractor_valleys[a]->Basin_Stored[Attractor_valleys[a]->Attractor[k+1]]->s[inp-1]))
               || ((k==Attractor_valleys[a]->N_a) &&
                   (Attractor_valleys[a]->Basin_Stored[Attractor_valleys[a]->Attractor[k]]->s[inp-1] != Attractor_valleys[a]->Basin_Stored[Attractor_valleys[a]->Attractor[1]]->s[inp-1]))){
                if(   ((k+1<Attractor_valleys[a]->N_a) &&
                      (Attractor_valleys[a]->Basin_Stored[Attractor_valleys[a]->Attractor[k+1]]->s[outp-1] != Attractor_valleys[a]->Basin_Stored[Attractor_valleys[a]->Attractor[k+2]]->s[outp-1]))
                   || ((k+1==Attractor_valleys[a]->N_a) &&
                       (Attractor_valleys[a]->Basin_Stored[Attractor_valleys[a]->Attractor[k+1]]->s[outp-1] != Attractor_valleys[a]->Basin_Stored[Attractor_valleys[a]->Attractor[1]]->s[outp-1])))
                    return(1);
            }
        return(0);
    }
    
    
    Boolean_Dynamics_TRACKER *mutant_dynamics(const char how[],int howmany){
        Boolean_Dynamics_TRACKER *mut;
        Boolean_RegNet *b_mut_net;
        b_mut_net = BRN->mutant_network(how,howmany);
        mut = new Boolean_Dynamics_TRACKER(b_mut_net,N_MAX);
        return(mut);
    }
};

typedef Boolean_Dynamics_TRACKER * Boolean_Dynamics_TRACKER_list;

void compare_tracking_to_exact_solution(Boolean_RegNet *A){
	Boolean_Dynamics_TRACKER *D;
    Boolean_Dynamics *D_exact;
    double beta,error_p;
    unsigned int i,j;
    unsigned long long k,exact_index;
    FILE *f;
    char fn[300];
    
    beta=2.;
    error_p=1./(1.+exp(2.*beta));
    //error_p=0.05;
    D=new Boolean_Dynamics_TRACKER(A,10);
    D->record_timetraces_from_collection_of_random_inputs(50,1000,error_p);
    D->print_ALL_attractor_info();
    printf("D->W_total_steps = %llu (should be %d)\n",D->W_total_steps,500*100);
   // getchar();
    
    D_exact=new Boolean_Dynamics(A);
    D_exact->generate_full_state_graph_synchronously(0);
	D_exact->calculate_noisy_landscape(beta);
	D_exact->weighted_network_of_Attractors(1);
    
    sprintf(fn,"%s/%s/__Scatter_exact_vs_sample.dat",
            MODEL_Directory,A->path);
    
    f=fopen(fn,"w");
    for(i=1;i<=D->Attractor_NR;i++){
        printf("Exact weight of attractor %d %lg\n\tsampled weight=%lg\n",i,
               D_exact->Attractor_G->node[i].node_props->nw, 
               D->Attractor_valleys[i]->W_Total/((double)D->W_total_steps));
        for(k=1;k<= D->Attractor_valleys[i]->N_ms;k++){
            printf("Sampled state %s ",D->Attractor_valleys[i]->Basin_Stored[k]->s);
            for(j=1;j<=A->N;j++) D_exact->s[j]=D->Attractor_valleys[i]->Basin_Stored[k]->s[j-1]-48;
            //printf("Matched state  ");D_exact->print_state();
            exact_index=D_exact->get_state_now();
            printf("ID %llu\n ",exact_index);
            fprintf(f,"%lg\t%lg\t%lg\t%d\n",D_exact->P_EQ_state[exact_index],
                    fabs(D_exact->P_EQ_state[exact_index]-D->Attractor_valleys[i]->Basin_Stored[k]->visits/((double)D->W_total_steps))/D_exact->P_EQ_state[exact_index],
                    D->Attractor_valleys[i]->Basin_Stored[k]->visits/((double)D->W_total_steps),i);
            printf("%lg\t%lg\t%d\n",D_exact->P_EQ_state[exact_index],
                    D->Attractor_valleys[i]->Basin_Stored[k]->visits/((double)D->W_total_steps),i);
           // getchar();
            
        }
    }
    fclose(f);
}
#endif
