//
//  Dynamical_Mod_Paired_Sampling.h
//  Erzso_Code_Central
//
//  Created by Erzsó Ravasz Regan on 1/29/13.
//  Copyright (c) 2013 Beth Israel Deaconess Medical Center. All rights reserved.
//

#ifndef Erzso_Code_Central_Dynamical_Mod_Paired_Sampling_h
#define Erzso_Code_Central_Dynamical_Mod_Paired_Sampling_h

/*
class Module_Attractor_Combination{
public:
    short int N_mod;
    unsigned int *Module_Attractor_IDs;
    longint_ARRAY *Uncoupled_Attractors_under_combination;

    Module_Attractor_Combination(int n){
        DynMod[i]->N_mod=n;
        DynMod[i]->Module_Attractor_IDs=new unsigned int[n+1];
        for(int k=1;k<=n;k++) Module_Attractor_IDs[i]=0;
        DynMod[i]->Uncoupled_Attractors_under_combination=new longint_ARRAY();
    }
};
*/
typedef std::map <std::string, longint_ARRAY> Module_Attactor_ID_Combination_MAP_type;


class Module_Landscape_in_Coupled_Attr{
public:
    int N_module_Phenotypes;
    doublepointer *Phenotype_Transition;
    double *Phenotype_Stability;
    
    Module_Landscape_in_Coupled_Attr(int n){
        N_module_Phenotypes=n;
        Phenotype_Stability=new double[N_module_Phenotypes+1];;
        Phenotype_Transition=new doublepointer[N_module_Phenotypes+1];
        for (int i=0; i<=N_module_Phenotypes; i++) {
            Phenotype_Stability[i]=0;
            Phenotype_Transition[i]=new double[N_module_Phenotypes+1];
            for (int j=0; j<=N_module_Phenotypes; j++) Phenotype_Transition[i][j]=0;
        }
    }
    ~Module_Landscape_in_Coupled_Attr(){
        for (int i=0; i<=N_module_Phenotypes; i++){
            delete[] Phenotype_Transition[i];Phenotype_Transition[i]=NULL;
        }
        delete[] Phenotype_Transition;Phenotype_Transition=NULL;
        delete[] Phenotype_Stability;Phenotype_Stability=NULL;
    }
};

typedef Module_Landscape_in_Coupled_Attr *Module_Landscape_in_Coupled_Attr_array;

class Module_Dynamical_Modularity_Measure{
 public:
    int N_module_Phenotypes,N_Coupled_Attr_nr;
    double *Module_Phenotype_Score;
    double Module_Score;
    Module_Landscape_in_Coupled_Attr_array *Land_CH;
    
    Module_Dynamical_Modularity_Measure(int n,int CA_nr){
        N_module_Phenotypes=n;
        N_Coupled_Attr_nr=CA_nr;
        Module_Score=0;
        Module_Phenotype_Score=new double[N_module_Phenotypes+1];
        for (int i=0; i<=N_module_Phenotypes; i++) {
            Module_Phenotype_Score[i]=0;
        }
        Land_CH=new Module_Landscape_in_Coupled_Attr_array[N_Coupled_Attr_nr+1];
        for (int i=1; i<=N_Coupled_Attr_nr; i++) {
            Land_CH[i]=new Module_Landscape_in_Coupled_Attr(N_module_Phenotypes);
        }
    }
    ~Module_Dynamical_Modularity_Measure(){
        //if(Module_Phenotype_Score!=NULL) {delete[] Module_Phenotype_Score; Module_Phenotype_Score=NULL;}
        for (int i=1; i<=N_Coupled_Attr_nr; i++) {delete Land_CH[i];Land_CH[i]=NULL;}
        delete[] Land_CH;Land_CH=NULL;
   }
};

typedef Module_Dynamical_Modularity_Measure *Module_Dynamical_Modularity_Measure_array;

class Dynamical_Modularity_Measure{
public:
    int Module_NR,N_nodes_in_Mod,C_Attr_NR;
    double *Attractor_modularity_measure,Switch_stability_measure,Module_coordination_measure,Module_Attractor_Distinction_measure,*Module_Barrier_Quality;
    Module_Dynamical_Modularity_Measure_array *M_DynMod;
   
    Dynamical_Modularity_Measure(int n1,int n2,int n3){
        Module_NR=n1;
        N_nodes_in_Mod=n2;
        C_Attr_NR=n3;
        Attractor_modularity_measure=new double[C_Attr_NR+1];
        for (int i=0; i<=C_Attr_NR; i++) {
            Attractor_modularity_measure[i]=0;
        }
        Attractor_modularity_measure[0]=1;
        
        M_DynMod=new Module_Dynamical_Modularity_Measure_array[Module_NR+1];
    }
    ~Dynamical_Modularity_Measure(){
        delete[] Attractor_modularity_measure;Attractor_modularity_measure=NULL;
       // for (int i=1; i<=Module_NR; i++) { delete M_DynMod[i];M_DynMod[i]=NULL;}
        delete[] M_DynMod;M_DynMod=NULL;
    }
};

class Module_Attractors_in_a_row{
public:
    unsigned int N_UAS;
    unsigned int *MA_Membership,*MAttractor_index;

    Module_Attractors_in_a_row( unsigned int n){
        N_UAS=n;
        MA_Membership=new unsigned int[N_UAS+1];
        MAttractor_index=new unsigned int[N_UAS+1];
        for (int i=0;i<=N_UAS; i++) {
            MA_Membership[i]=0;MAttractor_index[i]=0;
        }
    }
    ~Module_Attractors_in_a_row(){
        delete[] MA_Membership;MA_Membership=NULL;
        delete[] MAttractor_index;MAttractor_index=NULL;
    }
    void print_Module_Attractors_in_a_row(){
        printf("Total of %d module-attractor columns:\n",N_UAS);
        for (int i=1; i<=N_UAS; i++)  printf("%d\t",MA_Membership[i]);
        printf("\n");
        for (int i=1; i<=N_UAS; i++)  printf("%d\t",MAttractor_index[i]);
        printf("\n");
    }
};

class MODULAR_Boolean_Dynamics_PAIRED_TRACKER{
public:
	Boolean_RegNet *BRN,*BRN_U;
    Boolean_RegNet_array *BRN_Module;
	unsigned long long N_MAX;
    Boolean_Dynamics_TRACKER *Unconnected,*Coupled;
    Boolean_Dynamics_list *BRN_Module_Dynamics;
    
    Boolean_Dynamics_TRACKER_list *BRN_Module_Dynamics_TR;
    
    unsigned int *Module_ID,Module_NR,*Module_size,*ID_in_Module;
    double *x,*y;
    Module_Attactor_ID_Combination_MAP_type DynMod;
    Dynamical_Modularity_Measure *DynMod_Measure;
    char Bigname[390];
    
    intpointer *Unconnected_attractor_match_of_Composite_attractor;
    intpointer *Rank_in_matching_Unconnected_attractor_for_Composite_attractor;
    doublepointer *Unconnected_attractor_minHamming_of_Composite_attractor;
    int *CHOSEN_Unconnected_attractor_match_of_Composite_attractor;
    
	MODULAR_Boolean_Dynamics_PAIRED_TRACKER(Boolean_RegNet *b,unsigned long long nm){
        Unconnected_attractor_match_of_Composite_attractor=NULL;
        Rank_in_matching_Unconnected_attractor_for_Composite_attractor=NULL;
        Unconnected_attractor_minHamming_of_Composite_attractor=NULL;
        CHOSEN_Unconnected_attractor_match_of_Composite_attractor=NULL;
  
        sprintf(Bigname,"%s",b->name);
        BRN=b;
		N_MAX=nm;
        Module_ID=new unsigned int[BRN->N+1];
        ID_in_Module=new unsigned int[BRN->N+1];
        Module_NR=0;
        for(int i=1;i<=BRN->N;i++) {
            Module_ID[i]=BRN->MDAT->Module_ID[i];
            if(Module_ID[i]>Module_NR) Module_NR=Module_ID[i];
        }
        
        BRN_Module=NULL;  x=NULL;y=NULL;DynMod_Measure=NULL;
        generate_Unconnected_RNB(0);
        
        Unconnected =new Boolean_Dynamics_TRACKER(BRN_U ,N_MAX,Bigname,BRN_U->name);
        Coupled     =new Boolean_Dynamics_TRACKER(BRN   ,N_MAX);
        
        BRN_Module_Dynamics=NULL;
        BRN_Module_Dynamics_TR=NULL;
        
    }
    
    MODULAR_Boolean_Dynamics_PAIRED_TRACKER(Boolean_RegNet *b,unsigned long long nm,int rnd_version,int Keepgate,int Module_randomizer,int rand_all,longint_ARRAY *gl_keep_asis){
        Unconnected_attractor_match_of_Composite_attractor=NULL;
        Rank_in_matching_Unconnected_attractor_for_Composite_attractor=NULL;
        Unconnected_attractor_minHamming_of_Composite_attractor=NULL;
        CHOSEN_Unconnected_attractor_match_of_Composite_attractor=NULL;
  
		N_MAX=nm;
        sprintf(Bigname,"%s",b->name);
        BRN_Module=NULL;x=NULL;y=NULL; DynMod_Measure=NULL;
        
        if(rand_all==1){
            char fn[100];
            sprintf(fn,"%s_RndNetwork_GF_%d_Coupled",Bigname,rnd_version);
            BRN=new Boolean_RegNet(Bigname,fn,b->N);
            BRN->MDAT=new Boolean_RegNet_METADATA(BRN->N);
            BRN->MDAT=b->MDAT;
            
            if(! BRN->import_Boolean_RegNet()){
                
                // making a random network that is one smaller than cell cycle, and tacking on a self-looped GF node (at the same node-id as in cell cycle, with 2 random out-links to the rest of the network).
                Boolean_RegNet *b_prelim;
                
                // count links in and out of protected nodes
                int khold=0;
                for(int i=1;i<=gl_keep_asis->N;i++){
                    khold+=b->BN->node[gl_keep_asis->A[i]].k;
                }
                b_prelim=new Boolean_RegNet(b->N-(unsigned int)gl_keep_asis->N,b->K-khold,rnd_version);
                //b_prelim->export_Boolean_RegNet();
                //getchar();
                
               // unsigned int gf_id=b->MDAT->get_node_ID("GF");
               // if (gf_id==0) { printf("Problem: gf_id=0\n");exit(1);}
                
                snode::slinklist::slink *sc2;
                longint_ARRAY *gl;
                int shift_i;
                shift_i=0;
                for(int i=1;i<=b->N;i++){
                    if(gl_keep_asis->check_for_element(i)) shift_i++;
                    else{
                        sc2=b_prelim->BN->node[i-shift_i].Li.first_link_pt;
                        gl=new longint_ARRAY();
                        while(sc2!=NULL){
                            int shift=0;
                            for(int ii=1;ii<=gl_keep_asis->N;ii++)
                                if (sc2->sneighbor>=gl_keep_asis->A[ii]-shift) shift++;
                            gl->add_element(sc2->sneighbor+shift);
                            BRN->K++;
                            sc2=sc2->next(sc2);
                        }
                        
                        BRN->Gate[i]=new Boolean_Gate(b_prelim->Gate[i-shift_i]);
                       // printf("%d + %d  = %d - %llu\n",i,shift,i+shift,BRN->Gate[i+shift]->N);getchar();
                        for(unsigned long int j=1;j<=gl->N;j++)
                            BRN->BN->add_link(i,gl->A[gl->N-j+1]);
                        delete gl;gl=NULL;
                    }
                }
                
                delete b_prelim;b_prelim=NULL;
                
              //  for(int i=1;i<=BRN->N;i++)
               //     if(BRN->Gate[i]!=NULL) printf("%d - %llu\n",i,BRN->Gate[i]->N); else printf("%d - no gate yet!\n",i);
               // getchar();
                
                // links coming in to gate (outgoing on net, so  in node's list)

                for(int ii=1;ii<=gl_keep_asis->N;ii++){
                    sc2=b->BN->node[gl_keep_asis->A[ii]].Li.first_link_pt;
                    gl=new longint_ARRAY();
                    while(sc2!=NULL){
                        gl->add_element(sc2->sneighbor);
                        BRN->K++;
                        sc2=sc2->next(sc2);
                    }
                    BRN->Gate[gl_keep_asis->A[ii]]=new Boolean_Gate(b->Gate[gl_keep_asis->A[ii]]);
                   
                    // printf("%ld - %llu\n",gl_keep_asis->A[ii],BRN->Gate[gl_keep_asis->A[ii]]->N);getchar();
                    for(unsigned long int j=1;j<=gl->N;j++){
                        if (gl->A[gl->N-j+1]==gl_keep_asis->A[ii]) {
                            BRN->BN->add_link(gl_keep_asis->A[ii],gl_keep_asis->A[ii]);
                        }
                        else{
                            int ok=0;
                            while (ok<1) {
                                long long int a=rand_int(1,BRN->N);
                                if (!gl_keep_asis->check_for_element(a)) {
                                    BRN->BN->add_link(gl_keep_asis->A[ii],a);
                                    ok++;
                                }
                            }
                        }
                    }
                    delete gl;gl=NULL;
                   // BRN->print_gate(gl_keep_asis->A[ii]);getchar();
                    
                    for(unsigned long int j=1;j<=b->N;j++){
                        if ((j!=gl_keep_asis->A[ii])&&(b->BN->linked(j,gl_keep_asis->A[ii]))) {
                            int ok=0;
                            while (ok<1) {
                                long long int a=rand_int(1,BRN->N);
                                if((!gl_keep_asis->check_for_element(a))
                                    &&(!BRN->BN->linked(a,gl_keep_asis->A[ii]))) {
                                    BRN->BN->add_link(a,gl_keep_asis->A[ii]);
                                    ok++;
                                    BRN->Gate[a]->add_one_rnd_input_as_it_was_0_before();
                                   // BRN->print_gate(a);getchar();
                                }
                            }
                           
                        }
                    }
                }
                // links goint out of node (incoming on net, not in node's list)
                BRN->export_Boolean_RegNet();
            }
            
            Module_ID=new unsigned int[b->N+1];
            ID_in_Module=new unsigned int[b->N+1];
            Module_NR=0;
            for(int i=1;i<=b->N;i++) {
                Module_ID[i]=b->MDAT->Module_ID[i];
                if(Module_ID[i]>Module_NR) Module_NR=Module_ID[i];
            }
        
            
           // for(int j=1;j<=BRN->N;j++) {
           //     printf("Node j=%d %s\n",j,BRN->MDAT->Node_name[j]);
           //     BRN->print_gate(j);
             //
            //}
           // getchar();

            generate_Unconnected_RNB(rnd_version);
            
          //  BRN_U->export_Boolean_RegNet_Module_ID();
          //  BRN->export_Boolean_RegNet_Module_ID();
            
            Unconnected =new Boolean_Dynamics_TRACKER(BRN_U ,N_MAX,Bigname,BRN_U->name);
            //  printf("Unconnected: DIR %s   name %s    Sample %s\n", Unconnected->DIRname,Unconnected->BRN->name,Unconnected->SAMPLEname);
            Coupled     =new Boolean_Dynamics_TRACKER(BRN   ,N_MAX,Bigname,BRN->name);
            // printf("Coupled: DIR %s   name %s    Sample %s\n",Coupled->DIRname,Coupled->BRN->name,Coupled->SAMPLEname);getchar();
        }
        else {
            if(Module_randomizer==0){
                Module_ID=new unsigned int[b->N+1];
                ID_in_Module=new unsigned int[b->N+1];
                Module_NR=0;
                for(int i=1;i<=b->N;i++) {
                    Module_ID[i]=b->MDAT->Module_ID[i];
                    if(Module_ID[i]>Module_NR) Module_NR=Module_ID[i];
                }
                BRN=b;
                
                generate_Unconnected_RNB(0);
                Boolean_RegNet *b1;
                b1=randomly_shuffle_inter_module_links(b,rnd_version,Keepgate);
                BRN=b1;
                
                
                Unconnected =new Boolean_Dynamics_TRACKER(BRN_U ,N_MAX,Bigname,BRN_U->name);
              //  printf("Unconnected: DIR %s   name %s    Sample %s\n", Unconnected->DIRname,Unconnected->BRN->name,Unconnected->SAMPLEname);
                Coupled     =new Boolean_Dynamics_TRACKER(BRN   ,N_MAX,Bigname,BRN->name);
               // printf("Coupled: DIR %s   name %s    Sample %s\n",Coupled->DIRname,Coupled->BRN->name,Coupled->SAMPLEname);getchar();
            }
     
            else{
                // make sure name of coupled is the same as original network
                // name of Unconnected needs the random tag, and rules of decoupling are needded, now that the modules are randomized!!!
                
                //printf("starting with:\n");
                //b->MDAT->list_node_names();

                BRN=randomize_module_IDs(b,rnd_version,Module_randomizer,gl_keep_asis);
                // printf("ending with:\n");
               // BRN->MDAT->list_node_names();
                //getchar();
                
                
               /* int idgf=BRN->MDAT->get_node_ID("GF");
                if(BRN->MDAT->Module_ID[idgf]==0){
                    BRN->MDAT->Module_ID[idgf]=BRN->MDAT->Module_NR+1;
                    BRN->MDAT->Module_NR++;
                }
                */
                
               // printf("Module number of BRN now = %d\n",BRN->MDAT->Module_NR);getchar();
                
                Module_ID=new unsigned int[BRN->N+1];
                ID_in_Module=new unsigned int[BRN->N+1];
                Module_NR=0;
                for(int i=1;i<=BRN->N;i++) {
                    Module_ID[i]=BRN->MDAT->Module_ID[i];
                    if(Module_ID[i]>Module_NR) Module_NR=Module_ID[i];
                }
                generate_Unconnected_RNB_randomized(Module_randomizer,rnd_version);
                
               // printf("Number of nodes assigned to modules=%d\n",get_NR_of_module_assigned_nodes());getchar();
           
                Unconnected =new Boolean_Dynamics_TRACKER(BRN_U ,N_MAX,Bigname,BRN_U->name);
               // printf("Unconnected: DIR %s   name %s    Sample %s\n", Unconnected->DIRname,Unconnected->BRN->name,Unconnected->SAMPLEname);
                Coupled     =new Boolean_Dynamics_TRACKER(BRN   ,N_MAX,Bigname,b->name);
                //printf("Coupled: DIR %s   name %s    Sample %s\n",Coupled->DIRname,Coupled->BRN->name,Coupled->SAMPLEname);getchar();
                
                 printf("End: Module number  = %d\n",Module_NR);//getchar();
            }
        }
        
        BRN_Module_Dynamics=NULL;
        BRN_Module_Dynamics_TR=NULL;
        
       // printf("Unconnected:\n");
       // Unconnected->BRN->MDAT->list_node_names();
       // printf("Coupled:\n");
       // Coupled->BRN->MDAT->list_node_names();
       // getchar();
        
    }
    
	~MODULAR_Boolean_Dynamics_PAIRED_TRACKER(){
		if(Module_ID!=NULL) delete[] Module_ID;Module_ID=NULL;
        if(Module_size!=NULL) delete[] Module_size;Module_size=NULL;
        if(ID_in_Module!=NULL) delete[] ID_in_Module;ID_in_Module=NULL;
        //if(BRN_U!=NULL) delete BRN_U;BRN_U=NULL;
        //delete BRN;BRN=NULL;
       
        if(BRN_Module!=NULL){
            for(int i=1;i<=Module_NR;i++){
              //  if(BRN_Module[i]->MDAT!=NULL) {delete BRN_Module[i]->MDAT;  BRN_Module[i]->MDAT=NULL;}
              // delete BRN_Module[i];  BRN_Module[i]=NULL;
            }
            delete[] BRN_Module; BRN_Module=NULL;
        }
    
        delete Unconnected;Unconnected=NULL;
        
        if(Unconnected_attractor_match_of_Composite_attractor!=NULL){
            for(int i=1;i<=Coupled->Attractor_NR;i++){
                    if(Unconnected_attractor_match_of_Composite_attractor[i]!=NULL) delete[] Unconnected_attractor_match_of_Composite_attractor[i];Unconnected_attractor_match_of_Composite_attractor[i]=NULL;
                    if(Rank_in_matching_Unconnected_attractor_for_Composite_attractor[i]!=NULL) delete[] Rank_in_matching_Unconnected_attractor_for_Composite_attractor[i];Rank_in_matching_Unconnected_attractor_for_Composite_attractor[i]=NULL;
                    if(Unconnected_attractor_minHamming_of_Composite_attractor[i]!=NULL) delete[] Unconnected_attractor_minHamming_of_Composite_attractor[i];Unconnected_attractor_minHamming_of_Composite_attractor[i]=NULL;
            }
            delete[] Unconnected_attractor_match_of_Composite_attractor;Unconnected_attractor_match_of_Composite_attractor=NULL;
            delete[] Rank_in_matching_Unconnected_attractor_for_Composite_attractor;Unconnected_attractor_match_of_Composite_attractor=NULL;
            delete[] CHOSEN_Unconnected_attractor_match_of_Composite_attractor;Unconnected_attractor_match_of_Composite_attractor=NULL;
            delete[] Unconnected_attractor_minHamming_of_Composite_attractor;Unconnected_attractor_match_of_Composite_attractor=NULL;
        }
        delete Coupled;Coupled=NULL;
     
        
        if(DynMod_Measure!=NULL) delete DynMod_Measure;DynMod_Measure=NULL;
        if(x!=NULL) delete[] x;x=NULL;
        if(y!=NULL) delete[] y;y=NULL;
        if(BRN_Module_Dynamics!=NULL){
            delete[] BRN_Module_Dynamics;BRN_Module_Dynamics=NULL;
        }
    }
    
    Boolean_RegNet *randomize_module_IDs(Boolean_RegNet *b,int rnd_version,int Module_randomizer,longint_ARRAY *gl_keep_asis){
       
        Boolean_RegNet *br1;//,*br2;
        
        init_random(rnd_version);
        
        br1=new Boolean_RegNet(b->path);
        br1->MDAT->Module_ID=b->import_Boolean_RegNet_Module_IDs_but_not_the_gates(Module_randomizer,rnd_version);
        // printf("Number of nodes assigned to modules=%d\n",get_NR_of_module_assigned_nodes(br1));getchar();
        
        
        if (br1->MDAT->Module_ID==NULL) {
            br1->load_Module_IDs_old(b->MDAT->MOD_Assignment_name);
        
            switch(Module_randomizer){
                case 0: {// test case, leave module IDs unchanged
                    sprintf(br1->name,"%s_RMOD_test",b->name);
                } break;
                case 1: {
                    sprintf(br1->name,"%s_RMOD_MoveOne%d",b->name,rnd_version);
                    int moved=0;
                    do{
                        long int kmove=rand_int(1,br1->N);
                        if((br1->MDAT->Module_ID[kmove]!=0)&&(!gl_keep_asis->check_for_element(kmove))){
                            int destination=0;
                            do{
                                int kmodto=(int)rand_int(1,br1->MDAT->Module_NR);
                                if(kmodto!=br1->MDAT->Module_ID[kmove]){
                                    destination=1;
                                    br1->MDAT->Module_ID[kmove]=kmodto;
                                    moved=1;
                                }
                            }
                            while(destination==0);
                        }
                    }
                    while(moved==0);
                } break;
                case 2: {
                    sprintf(br1->name,"%s_RMOD_SwapOne%d",b->name,rnd_version);
                    int moved=0;
                    do{
                        long int kmove1=rand_int(1,br1->N);
                        long int kmove2=rand_int(1,br1->N);
                        if((!gl_keep_asis->check_for_element(kmove1))&&(!gl_keep_asis->check_for_element(kmove2))&&
                           (br1->MDAT->Module_ID[kmove1]!=br1->MDAT->Module_ID[kmove2])&&(br1->MDAT->Module_ID[kmove1]>0)&&(br1->MDAT->Module_ID[kmove2]>0)){
                                int m=br1->MDAT->Module_ID[kmove1];
                                br1->MDAT->Module_ID[kmove1]=br1->MDAT->Module_ID[kmove2];
                                br1->MDAT->Module_ID[kmove2]=m;
                                moved=1;
                        }
                    }
                    while(moved==0);
                } break;
                case 3: {
                    sprintf(br1->name,"%s_RMOD_RndPart%d",b->name,rnd_version);
                   // int *modcounts,*modnow;
                  //  modcounts=new int[br1->MDAT->Module_NR+1];
                  //  modnow=new int[br1->MDAT->Module_NR+1];
                   // for(int iii=0;iii<=br1->MDAT->Module_NR;iii++){
                   //     modcounts[iii]=0;modnow[iii]=0;
                   // }
                   // for(int iii=1;iii<=br1->N;iii++){
                   //     modcounts[br1->MDAT->Module_ID[iii]]++;
                   // }
                    
                    longint_ARRAY *glrnd=new longint_ARRAY();
                    for(int iii=1;iii<=br1->N;iii++){
                        if((br1->MDAT->Module_ID[iii]!=0)&&(!gl_keep_asis->check_for_element(iii))){
                            glrnd->add_element(iii);
                        }
                    }
                    int *xind=new int[glrnd->N+1];
                    for(int iii=1;iii<=glrnd->N;iii++) xind[iii]=br1->MDAT->Module_ID[glrnd->A[iii]];
                    randomize(xind,glrnd->N);
                    
                    for(int iii=1;iii<=glrnd->N;iii++){
                        br1->MDAT->Module_ID[glrnd->A[iii]]=xind[iii];
                        printf("node %d %s -- new module ID = %d\n",iii,br1->MDAT->Node_name[iii],br1->MDAT->Module_ID[iii]);
                    }
                    // getchar();
                } break;
            }
        }
        else {
            br1->MDAT->Module_NR=br1->MDAT->Module_ID[0];
            switch(Module_randomizer){
                case 0: {sprintf(br1->name,"%s_RMOD_test",b->name);} break;
                case 1: {sprintf(br1->name,"%s_RMOD_MoveOne%d",b->name,rnd_version);} break;
                case 2: { sprintf(br1->name,"%s_RMOD_SwapOne%d",b->name,rnd_version);} break;
                case 3: { sprintf(br1->name,"%s_RMOD_RndPart%d",b->name,rnd_version);} break;
            }
        }

        return(br1);
    }
    
    Boolean_RegNet *randomly_shuffle_inter_module_links(Boolean_RegNet *b,int rnd_version,int Keepgate){
        Boolean_RegNet *br1;//,*br2;
        
        char fn[300];
        snode::slinklist::slink *sc2;
        int *canalyze;
        longint_ARRAY *gl;
        unsigned long int newn;
        init_random(rnd_version);
        
        if(Keepgate) sprintf(fn,"%s_R%d",b->name,rnd_version);
        else sprintf(fn,"%s_RG%d",b->name,rnd_version);
        
        br1=new Boolean_RegNet(b->path,fn,b->N);
        br1->MDAT=b->MDAT;
        if(br1->import_Boolean_RegNet()) {
            printf("Imported %s from file!\n",br1->name);
          //  for(int i=1;i<=br1->N;i++) {
          //      printf("Imported gate:\n");
          //      br1->print_gate(i);
           //     printf("vs. Cell Cycle gate:\n");
           //     b->print_gate(i);
             //   getchar();
           // }
            return(br1);
        }
        
        printf("Generating random intermodule links, version %d\n",rnd_version);
      //  getchar();
		for(int j=1;j<=b->N;j++){
     //       printf("Node j=%d %s\n",j,b->MDAT->Node_name[j]);
           // canalyze=read_canalyzing_values(j,b);
            canalyze=generate_best_canalyzing_value_choice(j,BRN);
            
          //  for(int l=1;l<=b->N;l++)
          //      if(canalyze[l]!=-1)
           //         printf("\t <- %d %s : %d\n",l,b->MDAT->Node_name[l],canalyze[l]);
            
            gl=new longint_ARRAY();
            
            sc2=b->BN->node[j].Li.first_link_pt;
            while(sc2!=NULL){
                if(Module_ID[sc2->sneighbor]==Module_ID[j]){
                    if(canalyze[sc2->sneighbor]>-1) {printf("Error: in canalizing list of node %s:\nLink %d <- %lld is internal to module %d, but input %s is set to canalysing value %d\n",b->MDAT->Node_name[j],j,sc2->sneighbor, Module_ID[j],b->MDAT->Node_name[sc2->sneighbor],canalyze[sc2->sneighbor]);//exit(1);
                        getchar();
                    }
                    else {canalyze[sc2->sneighbor]=-2;
                        gl->add_element(sc2->sneighbor);
                    }
                }
                else {
                    if(Module_ID[sc2->sneighbor]==0) gl->add_element(sc2->sneighbor);
                    else{
                        do{ newn=random_node_from_module(Module_ID[sc2->sneighbor],b);
                            if(gl->check_for_element(newn)) newn=0;
                        }
                        while(newn==0);
                        gl->add_element(newn);
                       // printf("Rewired link %lld %s -> %d %s to %ld %s -> %d %s\n",sc2->sneighbor,b->MDAT->Node_name[sc2->sneighbor],j,b->MDAT->Node_name[j],newn,b->MDAT->Node_name[newn],j,b->MDAT->Node_name[j]);
                       
                    }
                    if(canalyze[sc2->sneighbor]==-1) {printf("Error: in canalizing list of node %s:\nLink %d <- %lld is external to module %d, but the canalysing value of input %s (mod %d)is not set\n",b->MDAT->Node_name[j],j,sc2->sneighbor, Module_ID[j],b->MDAT->Node_name[sc2->sneighbor],Module_ID[sc2->sneighbor]);
                        getchar();
                      //  exit(1);
                    }
                }
                sc2=sc2->next(sc2);
            }
           // getchar();
            
            Boolean_Gate *G_test=new Boolean_Gate(b->Gate[j]);
            for(unsigned long int l=0;l<gl->N;l++){
                br1->BN->add_link(j,gl->A[gl->N-l]);
                br1->K++;
            }
            
            br1->Gate[j]=G_test;
            delete[] canalyze;canalyze=NULL;
            delete gl;gl=NULL;
        }
        
        if (!Keepgate) {
            randomize_gate_for_intermodule_inputs(br1);
        }
        br1->export_Boolean_RegNet();
       
        //for(int i=1;i<=br1->N;i++) {
        //    printf("Exported gate:\n");
        //    br1->print_gate(i);
        //    printf("vs. Cell Cycle gate:\n");
        //    b->print_gate(i);
        //    getchar();
        //}
        return(br1);
    }
    
    void randomize_gate_for_intermodule_inputs(Boolean_RegNet *b){
        snode::slinklist::slink *sc2;
        int *canalyze;
        
        for(int j=1;j<=b->N;j++){
           // printf("Node j=%d %s\n",j,b->MDAT->Node_name[j]);
            //canalyze=read_canalyzing_values(j,b);
            canalyze=generate_best_canalyzing_value_choice(j,b);
           // for(int l=1;l<=b->N;l++)
            //    if(canalyze[l]!=-1)
            //        printf("\t <- %d %s : %d\n",l,b->MDAT->Node_name[l],canalyze[l]);
            
            int k=0;
            sc2=b->BN->node[j].Li.first_link_pt;
            while(sc2!=NULL){
                k++;
                if((Module_ID[sc2->sneighbor]!=Module_ID[j])&&(Module_ID[sc2->sneighbor]!=0)){
                  //  printf("Before randomizing gate for input %s\n",b->MDAT->Node_name[sc2->sneighbor]);
                  //  b->print_gate(j);
            
                    b->Gate[j]->randomize_function_of_input(k,canalyze[sc2->sneighbor]);
                    
                 //   printf("AFTER randomizing gate for input %s\n",b->MDAT->Node_name[sc2->sneighbor]);
                 //   b->print_gate(j);
                 //   getchar();
                }
                sc2=sc2->next(sc2);
            }
        }
    }
    
    unsigned long int random_node_from_module(int mid,Boolean_RegNet *b){
        unsigned long ninmodule=0;
        for(int i=1;i<=b->N;i++){
            if (b->MDAT->Module_ID[i]==mid) {
                ninmodule++;
            }
        }
        unsigned long id1=rand_int(1,ninmodule);
        ninmodule=0;
        for(int i=1;i<=b->N;i++){
            if (b->MDAT->Module_ID[i]==mid) {
                ninmodule++;
                if(ninmodule==id1) return(i);
            }
        }
        printf("Something's worng, should not have gone this far!\n");exit(1);
    }

    unsigned long int get_NODE_ID_from_node_label_line(){
        unsigned long int ind1;
        char attr_list[50][5000];
        unsigned long int Nf;
        char * pch;
        
        pch = strstr (GSE_line,"label=");
        Nf=separate_line(pch,attr_list,'"',50);
        //printf("found <node for %s\n",attr_list[2]);getchar();
        ind1=BRN->MDAT->get_node_ID(attr_list[2]);
        return(ind1);
    }
    
    unsigned long int get_NODE_ID_from_node_label_line_gDetails(){
        unsigned long int ind1;
        char attr_list[50][5000];
        unsigned long int Nf;
        
        Nf=separate_line(GSE_line,attr_list,'"',50);
        //printf("Found node %s\n",attr_list[2]); getchar();
        ind1=BRN->MDAT->get_node_ID(attr_list[2]);
        return(ind1);
    }
    
    void read_Node_positions_from_DMMS_gDetails(){
        char fn[300], bla[100];
        FILE *f;
        unsigned long int s;
        double *x,*y,*r,*g,*b;
       // double gg;
        sprintf(fn,"%s/%s/%s_Fine_gDetails.txt",
                MODEL_Directory,Bigname,Bigname);
        f=NULL;
        while (f==NULL) {
            f=fopen(fn, "r");
            if(f==NULL){ printf("Please use dynamicalverify -g %s.dmms to export current node coordinates \n\t THEN hit Enter\n",Bigname);
                getchar();
            }
        }
        x=new double[BRN->N+1];
        y=new double[BRN->N+1];
        r=new double[BRN->N+1];
        g=new double[BRN->N+1];
        b=new double[BRN->N+1];
        int *ok;
        
        ok=new int[BRN->N+1];
        for (int i=1; i<=BRN->N; i++) ok[i]=0;
        
        while (get_line_and_check(f)) {
            if(strstr(GSE_line,"    ( ")!=NULL){
                s=get_NODE_ID_from_node_label_line_gDetails();
                // printf("found %ld\n",s);getchar();
                
                get_line_and_check(f);
                sscanf(GSE_line,"%s%s%lg%lg%lg",bla,bla,&(r[s]),&(g[s]),&(b[s]));
                get_line_and_check(f); // ,
                
                get_line_and_check(f); //  [ 280.1    // x-coord
                sscanf(GSE_line,"%s%lg",bla,&(x[s]));
                get_line_and_check(f); //  , 100.3    // y-coord
                sscanf(GSE_line,"%s%lg",bla,&(y[s]));
                get_line_and_check(f); // ]
                get_line_and_check(f); // )
                get_line_and_check(f); // ,
                ok[s]=1;
            }
        }
        fclose(f);
        ok[0]=1;
        for (int i=1; i<=BRN->N; i++)
            if(ok[i]==0){
                printf("Node %s has no coordinates\n",BRN->MDAT->Node_name[i]);
                ok[0]=0;
            }
           // else printf("%s: %lg , %lg\n",BRN->MDAT->Node_name[i],x[i],y[i]);
        if (ok[0]==0)  getchar();
        
        
        BRN->MDAT->x=x;
        BRN->MDAT->y=y;
        BRN->MDAT->r=r;
        BRN->MDAT->g=g;
        BRN->MDAT->b=b;
    }
    
    void read_Node_Cytoskape_positions(){
        char fn[300];
        FILE *f;
        unsigned long int s,Nf;
        char attr_list[50][5000];
        double *x,*y;
       	sprintf(fn,"%s/%s/POZ_%s.xgmml",
                MODEL_Directory,Bigname,Bigname);
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
                s=get_NODE_ID_from_node_label_line();
               // printf("found %ld\n",s);getchar();
                
                get_line_and_check(f);
                while(strstr(GSE_line,"<graphics ")==NULL){
                    get_line_and_check(f);
                }
                char * pch;
                pch = strstr (GSE_line,"x=");
                Nf=separate_line(pch,attr_list,'"',50);
                x[s]=atof(attr_list[2]);
                
                pch = strstr (GSE_line,"y=");
                Nf=separate_line(pch,attr_list,'"',50);
                y[s]=atof(attr_list[2]);
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
    
    void transfer_Node_Cytoskape_positions(){
        
        BRN->MDAT->x=x;
        BRN->MDAT->y=y;
        BRN_U->MDAT->x=x;
        BRN_U->MDAT->y=y;
        
        for(int j=1;j<=BRN->N;j++){
            if(Module_ID[j]>0) {
                BRN_Module[Module_ID[j]]->MDAT->x[ID_in_Module[j]]=BRN->MDAT->x[j];
                BRN_Module[Module_ID[j]]->MDAT->y[ID_in_Module[j]]=BRN->MDAT->y[j];
            }
        }
    }
    
    
    /*
    int *read_canalyzing_values(int BRN_index,Boolean_RegNet *b){
        FILE *f;
        char fn[300];
        int *can,icanal,c;
        
        can=new int[b->N+1];
        for(int i=0;i<=b->N;i++) can[i]=-1;
       
        sprintf(fn,"%s/%s/tt_decomp/BIOCOMP_%s.csv",MODEL_Directory,Bigname,b->MDAT->Node_name[BRN_index]);
        
        f=fopen(fn,"r");
        if(f==NULL) { //printf("No Decomposition rule for node %s\n",b->MDAT->Node_name[BRN_index]);
            return(can);
        }
		char bla[100];
        
        while (fgets(GSE_line,MAX_LINE,f)!=NULL) {
            sscanf(GSE_line,"%s %d",bla,&c);
            icanal=b->MDAT->get_node_ID(bla);
            can[icanal]=c;
        }
        fclose(f);
        return(can);
    }
    */
    
    int *generate_best_canalyzing_value_choice(int BRN_index,Boolean_RegNet *b){
        int *can,*outside;
        snode::slinklist::slink *sc2;
        
        can=new int[b->N+1];
        outside=new int[b->N+1];
        for(int i=0;i<=b->N;i++) {can[i]=-1;outside[i]=0;} // not a neighbor
        
        sc2=BRN->BN->node[BRN_index].Li.first_link_pt;
        while(sc2!=NULL){
            if(Module_ID[sc2->sneighbor]!=Module_ID[BRN_index]){
                outside[0]++;
                outside[sc2->sneighbor]=outside[0];
               // printf("Link %s to %s comes from outside:\n", BRN->MDAT->Node_name[sc2->sneighbor],BRN->MDAT->Node_name[BRN_index]);
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
            for(int k=1;k<=b->N;k++)
                if(outside[k]>0){
                    can[k]=s_outmodlinks[outside[k]-1];
          //          printf("\t k=%d (%s) value = %d\n",k,BRN->MDAT->Node_name[k],can[k]);
                }
            longint_ARRAY *gl;
            gl=new longint_ARRAY();
            Boolean_Gate *G_test=BRN->get_PARTIAL_gate(BRN_index,can,gl);
            
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
            //    printf("There is no regulation left for this node! locked in with canal. inputs all 0.\n");
                i_choice=0;
            }
        }
        
        
        s_outmodlinks=new unsigned short int[Ncomb+1];
        decimal_to_any_base(i_choice,2,outside[0], s_outmodlinks);
       // printf("Final choice: %lld:\n",i_choice);
        for(int k=1;k<=b->N;k++)
            if(outside[k]>0){
                can[k]=s_outmodlinks[outside[k]-1];
               // printf("\t k=%d (%s) value = %d\n",k,BRN->MDAT->Node_name[k],can[k]);
            }
        delete[] s_outmodlinks; s_outmodlinks=NULL;
        
        s_outmodlinks=new unsigned short int[Ncomb+1];
        decimal_to_any_base(i_choice,2,outside[0], s_outmodlinks);
        for(int k=1;k<=b->N;k++)
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
    
   
    
    int *generate_ENV_canalyzing_values(int BRN_index,Boolean_RegNet *b,char g_s[]){
        int *can,*outside;
        snode::slinklist::slink *sc2;
        
        can=new int[b->N+1];
        outside=new int[b->N+1];
        for(int i=0;i<=b->N;i++) {can[i]=-1;outside[i]=0;} // not a neighbor
        
        sc2=BRN->BN->node[BRN_index].Li.first_link_pt;
        while(sc2!=NULL){
            if(Module_ID[sc2->sneighbor]!=Module_ID[BRN_index]){
                outside[0]++;
                outside[sc2->sneighbor]=outside[0];
                // printf("Link %s to %s comes from outside:\n", BRN->MDAT->Node_name[sc2->sneighbor],BRN->MDAT->Node_name[BRN_index]);
            }
            else can[sc2->sneighbor]=-2; // non-canalyzed neighbor
            sc2=sc2->next(sc2);
        }
        
        for(int k=1;k<=b->N;k++)
            if(outside[k]>0){
                can[k]=g_s[k-1]-48;
                printf("\t k=%d (%s) value = %d\n",k,BRN->MDAT->Node_name[k],can[k]);
            }
      //  longint_ARRAY *gl;
      //  gl=new longint_ARRAY();
      //  Boolean_Gate *G_test=BRN->get_PARTIAL_gate(BRN_index,can,gl);
        
       // printf("Partial gate returned is:\n");
      //  BRN->print_PARTIAL_gate(BRN_index,gl,G_test);
      //  getchar();
        
        delete[] outside;outside=NULL;
        return(can);
    }
    
    Boolean_RegNet *Module_Network_with_global_state_as_environment(int m_id,char s_glob[],char env_name[]){
        Boolean_RegNet *br_env;
        char fn[300];
        int *canalyze;
        longint_ARRAY *gl;
        snode::slinklist::slink *sc2;
        
        sprintf(fn,"%s__M-%d_env_%s",Bigname,m_id,env_name);
        br_env=new Boolean_RegNet(BRN->path,fn,Module_size[m_id]);
        br_env->MDAT=new Boolean_RegNet_METADATA(br_env->N);

        for(int j=1;j<=BRN->N;j++){
            if(Module_ID[j]==m_id) {
                strcpy(br_env->MDAT->Node_name[ID_in_Module[j]],BRN->MDAT->Node_name[j]);
                br_env->MDAT->x[ID_in_Module[j]]=BRN->MDAT->x[j];
                br_env->MDAT->y[ID_in_Module[j]]=BRN->MDAT->y[j];
            }
        }
        
        for(int j=1;j<=BRN->N;j++)
            if(Module_ID[j]==m_id){
                 printf("Node j=%d %s\n",j,BRN->MDAT->Node_name[j]);
                  BRN->print_gate(j);
            
                canalyze=generate_ENV_canalyzing_values(j,BRN,s_glob);
            
//                for(int l=1;l<=BRN->N;l++)
//                  if(canalyze[l]!=-1)
//                       printf("\t <- %d %s : %d\n",l,BRN->MDAT->Node_name[l],canalyze[l]);
//
                sc2=BRN->BN->node[j].Li.first_link_pt;
                while(sc2!=NULL){
                    if(Module_ID[sc2->sneighbor]==Module_ID[j]){
                        if(canalyze[sc2->sneighbor]>-1) {printf("Error: in canalizing list of node %s:\nLink %d <- %lld is internal to module %d, but input %s is set to canalysing value %d\n",BRN->MDAT->Node_name[j],j,sc2->sneighbor, Module_ID[j],BRN->MDAT->Node_name[sc2->sneighbor],canalyze[sc2->sneighbor]);//exit(1);
                            getchar();
                        }
                        else canalyze[sc2->sneighbor]=-2;
                    }
                    else {if(canalyze[sc2->sneighbor]==-1) {printf("Error: in canalizing list of node %s:\nLink %d <- %lld is external to module %d, but the canalysing value of input %s (mod %d)is not set\n",BRN->MDAT->Node_name[j],j,sc2->sneighbor, Module_ID[j],BRN->MDAT->Node_name[sc2->sneighbor],Module_ID[sc2->sneighbor]);//exit(1);
                        getchar();
                    }
                    }
                    sc2=sc2->next(sc2);
                }
                
                gl=new longint_ARRAY();
                
                Boolean_Gate *G_test=BRN->get_PARTIAL_gate(j,canalyze,gl);
                
                 BRN->print_PARTIAL_gate(j,gl,G_test);//getchar();
                
                for(unsigned long int l=0;l<gl->N;l++){
                    br_env->BN->add_link(ID_in_Module[j],
                                                           ID_in_Module[gl->A[gl->N-l]]);
                    //printf("Adding link in module %d into %s (module id = %d) from %s (module id =  %d)\n",
                    //        Module_ID[j],br_env[Module_ID[j]]->MDAT->Node_name[ID_in_Module[j]],             ID_in_Module[j],
                    //                     br_env[Module_ID[j]]->MDAT->Node_name[ID_in_Module[gl->A[gl->N-l]]],ID_in_Module[gl->A[gl->N-l]]);
                    br_env->K++;
                }
                if(gl->N==0) {
                    br_env->BN->add_link(ID_in_Module[j],ID_in_Module[j]);
                        //printf("Adding link in module %d into %s (module id = %d) from %s (module id =  %d)\n",
                        //       Module_ID[j],br_env[Module_ID[j]]->MDAT->Node_name[ID_in_Module[j]],ID_in_Module[j],
                        //                    br_env[Module_ID[j]]->MDAT->Node_name[ID_in_Module[j]],ID_in_Module[j]);
                        //snode::slinklist::slink *sc2;
                        //sc2=br_env[Module_ID[j]]->BN->node[ID_in_Module[j]].Li.first_link_pt;
                        //if (sc2==NULL) {
                        //    printf("\t\t Link was not added!\n");getchar();
                        //}
                        //else printf("\t\t%lld\n",sc2->sneighbor);
                }
                
               br_env->Gate[ID_in_Module[j]]=new Boolean_Gate(G_test);
                
                //    br_env[Module_ID[j]]->print_gate(ID_in_Module[j]);
                // printf("\nbr_env[%d]->Gate[%d]->K=%d",Module_ID[j],ID_in_Module[j],br_env[Module_ID[j]]->Gate[ID_in_Module[j]]->K);getchar();
                //}
                
                
                delete[] canalyze;canalyze=NULL;
                delete gl;gl=NULL;
        }
        return(br_env);
    }
    
    void generate_Unconnected_RNB(int newmodules){
        char fn[300];
        snode::slinklist::slink *sc2;
        int *canalyze;
        longint_ARRAY *gl;
        
        if (newmodules) sprintf(fn,"%s_RndNetwork%d_Uncoupled",Bigname,newmodules);
        else sprintf(fn,"%s__Uncoupled",Bigname);
        BRN_U=new Boolean_RegNet(BRN->path,fn,BRN->N);
        BRN_U->MDAT=BRN->MDAT;
        
        Module_size=new unsigned int[Module_NR+1];
        for(int i=0;i<=Module_NR;i++) Module_size[i]=0;
        for(int i=1;i<=BRN->N;i++) {
            Module_size[Module_ID[i]]++;
            ID_in_Module[i]=Module_size[Module_ID[i]];
        }

        BRN_Module=new Boolean_RegNet_array[Module_NR+1];

        for(int i=1;i<=Module_NR;i++) {
            if(newmodules) sprintf(fn,"%s_RndNetwork%d_Uncoupled__M-%d",Bigname,newmodules,i);
                else  //sprintf(fn,"%s__M-%d",Bigname,i);
                      sprintf(fn,"MODULE_%d_%s",i,BRN->Module_Names[i].c_str());
            char fn3[300]; sprintf(fn3,"%s/%s",BRN->path,fn);
            BRN_Module[i]=new Boolean_RegNet(fn3,fn,Module_size[i]);
            BRN_Module[i]->MDAT=new Boolean_RegNet_METADATA(BRN_Module[i]->N);
        }

 
        for(int j=1;j<=BRN->N;j++){
            printf("j=%d  Mod ID=%d ",j,Module_ID[j]);
            printf("ID_in_Mod=%d\n",ID_in_Module[j]);
            
            if(Module_ID[j]>0) {
                strcpy(BRN_Module[Module_ID[j]]->MDAT->Node_name[ID_in_Module[j]],BRN->MDAT->Node_name[j]);
                BRN_Module[Module_ID[j]]->MDAT->x[ID_in_Module[j]]=BRN->MDAT->x[j];
                BRN_Module[Module_ID[j]]->MDAT->y[ID_in_Module[j]]=BRN->MDAT->y[j];
            }
        }
        

		for(int j=1;j<=BRN->N;j++){
           // printf("Node j=%d %s\n",j,BRN->MDAT->Node_name[j]);
           // BRN->print_gate(j);
       
            //canalyze=read_canalyzing_values(j,BRN);
            
            canalyze=generate_best_canalyzing_value_choice(j,BRN);
            
            // for(int l=1;l<=BRN->N;l++)
           //     if(canalyze[l]!=-1)
           //         printf("\t <- %d %s : %d\n",l,BRN->MDAT->Node_name[l],canalyze[l]);
           
            sc2=BRN->BN->node[j].Li.first_link_pt;
            while(sc2!=NULL){
                if (sc2->sneighbor==0) {
                    printf("Node %d has a neighbor 0 for some reason!\n",j);getchar();
                }
                if(Module_ID[sc2->sneighbor]==Module_ID[j]){
                    if(canalyze[sc2->sneighbor]>-1) {printf("Error: in canalizing list of node %s:\nLink %d <- %lld is internal to module %d, but input %s is set to canalysing value %d\n",BRN->MDAT->Node_name[j],j,sc2->sneighbor, Module_ID[j],BRN->MDAT->Node_name[sc2->sneighbor],canalyze[sc2->sneighbor]);//exit(1);
                        getchar();
                    }
                    else canalyze[sc2->sneighbor]=-2;
                }
                else {if(canalyze[sc2->sneighbor]==-1) {printf("Error: in canalizing list of node %s:\nLink %d <- %lld is external to module %d, but the canalysing value of input %s (mod %d)is not set\n",BRN->MDAT->Node_name[j],j,sc2->sneighbor, Module_ID[j],BRN->MDAT->Node_name[sc2->sneighbor],Module_ID[sc2->sneighbor]);//exit(1);
                        getchar();
                    }
                }
                sc2=sc2->next(sc2);
            }
            

            gl=new longint_ARRAY();

            Boolean_Gate *G_test=BRN->get_PARTIAL_gate(j,canalyze,gl);

           // BRN->print_PARTIAL_gate(j,gl,G_test);getchar();

           for(unsigned long int l=0;l<gl->N;l++){
                BRN_U->BN->add_link(j,gl->A[gl->N-l]);
                BRN_U->K++;
                if(Module_ID[j]>0) {
                    BRN_Module[Module_ID[j]]->BN->add_link(ID_in_Module[j],
                                                           ID_in_Module[gl->A[gl->N-l]]);
                    //printf("Adding link in module %d into %s (module id = %d) from %s (module id =  %d)\n",
                    //        Module_ID[j],BRN_Module[Module_ID[j]]->MDAT->Node_name[ID_in_Module[j]],             ID_in_Module[j],
                    //                     BRN_Module[Module_ID[j]]->MDAT->Node_name[ID_in_Module[gl->A[gl->N-l]]],ID_in_Module[gl->A[gl->N-l]]);
                    BRN_Module[Module_ID[j]]->K++;
                }
            }
            if(gl->N==0) {
                BRN_U->BN->add_link(j,j);
                if(Module_ID[j]>0) {
                    BRN_Module[Module_ID[j]]->BN->add_link(ID_in_Module[j],ID_in_Module[j]);
                    //printf("Adding link in module %d into %s (module id = %d) from %s (module id =  %d)\n",
                    //       Module_ID[j],BRN_Module[Module_ID[j]]->MDAT->Node_name[ID_in_Module[j]],ID_in_Module[j],
                    //                    BRN_Module[Module_ID[j]]->MDAT->Node_name[ID_in_Module[j]],ID_in_Module[j]);
                    //snode::slinklist::slink *sc2;
                    //sc2=BRN_Module[Module_ID[j]]->BN->node[ID_in_Module[j]].Li.first_link_pt;
                    //if (sc2==NULL) {
                    //    printf("\t\t Link was not added!\n");getchar();
                    //}
                    //else printf("\t\t%lld\n",sc2->sneighbor);
                }
            }
            
            BRN_U->Gate[j]=G_test;
            if(Module_ID[j]>0) BRN_Module[Module_ID[j]]->Gate[ID_in_Module[j]]=new Boolean_Gate(G_test);
         
            //  BRN_U->print_gate(j);
           // if(Module_ID[j]>0) {
            //    BRN_Module[Module_ID[j]]->print_gate(ID_in_Module[j]);
               // printf("\nBRN_Module[%d]->Gate[%d]->K=%d",Module_ID[j],ID_in_Module[j],BRN_Module[Module_ID[j]]->Gate[ID_in_Module[j]]->K);getchar();
            //}
 

            delete[] canalyze;canalyze=NULL;
            delete gl;gl=NULL;
        }
    }
   
    void generate_Unconnected_RNB_randomized(int Module_randomizer,int rnd_version){
        char fn[300];
        snode::slinklist::slink *sc2;
        int *canalyze;
        longint_ARRAY *gl;
        
        sprintf(fn,"%s__Uncoupled",BRN->name);
       // printf("Trying to read in %s\n",fn);getchar();
        BRN_U=new Boolean_RegNet(BRN->path,fn,BRN->N);
        BRN_U->MDAT=BRN->MDAT;
        
              
        Module_size=new unsigned int[Module_NR+1];
        for(int i=0;i<=Module_NR;i++) Module_size[i]=0;
        for(int i=1;i<=BRN->N;i++) {
            Module_size[Module_ID[i]]++;
            ID_in_Module[i]=Module_size[Module_ID[i]];
        }
        
        BRN_Module=new Boolean_RegNet_array[Module_NR+1];
        for(int i=1;i<=Module_NR;i++) {
            sprintf(fn,"%s/MODULE_%d_%s",BRN_U->name,i,BRN->Module_Names[i].c_str());
            //printf("RNd module name = %s\n",fn);getchar();
            BRN_Module[i]=new Boolean_RegNet(BRN->path,fn,Module_size[i]);
            BRN_Module[i]->MDAT=new Boolean_RegNet_METADATA(BRN_Module[i]->N);
        }
        for(int j=1;j<=BRN->N;j++){
           // printf("j=%d %s Mod ID=%d\n ",j,BRN->MDAT->Node_name[j],Module_ID[j]);
           // printf("ID_in_Mod=%d\n",ID_in_Module[j]);
            
            if(Module_ID[j]>0) strcpy(BRN_Module[Module_ID[j]]->MDAT->Node_name[ID_in_Module[j]],BRN->MDAT->Node_name[j]);
        }
        
		for(int j=1;j<=BRN->N;j++){
           // printf("\n\n\nNode j=%d %s\n",j,BRN->MDAT->Node_name[j]);
           // BRN->print_gate(j);
            canalyze=generate_best_canalyzing_value_choice(j,BRN);
          //  for(int l=1;l<=BRN->N;l++)
          //      if(canalyze[l]!=-1)
          //          printf("\t <- %d %s : %d\n",l,BRN->MDAT->Node_name[l],canalyze[l]);
           // getchar();
            
            sc2=BRN->BN->node[j].Li.first_link_pt;
            while(sc2!=NULL){
                if(Module_ID[sc2->sneighbor]==Module_ID[j]){
                    if(canalyze[sc2->sneighbor]>-1) {printf("Error: in canalizing list of node %s:\nLink %d <- %lld is internal to module %d, but input %s is set to canalysing value %d\n",BRN->MDAT->Node_name[j],j,sc2->sneighbor, Module_ID[j],BRN->MDAT->Node_name[sc2->sneighbor],canalyze[sc2->sneighbor]);exit(1);}
                    else canalyze[sc2->sneighbor]=-2;
                }
                else {if(canalyze[sc2->sneighbor]==-1) {printf("Error: in canalizing list of node %s:\nLink %d <- %lld is external to module %d, but the canalysing value of input %s (mod %d)is not set\n",BRN->MDAT->Node_name[j],j,sc2->sneighbor, Module_ID[j],BRN->MDAT->Node_name[sc2->sneighbor],Module_ID[sc2->sneighbor]);exit(1);};
                }
                sc2=sc2->next(sc2);
            }
            
            gl=new longint_ARRAY();
            Boolean_Gate *G_test=BRN->get_PARTIAL_gate(j,canalyze,gl);
            //BRN->print_PARTIAL_gate(j,gl,G_test);
            
            for(unsigned long int l=0;l<gl->N;l++){
                BRN_U->BN->add_link(j,gl->A[gl->N-l]);
                BRN_U->K++;
                if(Module_ID[j]>0) {
                    BRN_Module[Module_ID[j]]->BN->add_link(ID_in_Module[j],
                                                           ID_in_Module[gl->A[gl->N-l]]);
                    BRN_Module[Module_ID[j]]->K++;
                }
            }
            if(gl->N==0) {
                BRN_U->BN->add_link(j,j);
                if(Module_ID[j]>0) BRN_Module[Module_ID[j]]->BN->add_link(ID_in_Module[j],ID_in_Module[j]);
            }
            BRN_U->Gate[j]=G_test;
            if(Module_ID[j]>0) BRN_Module[Module_ID[j]]->Gate[ID_in_Module[j]]=G_test;
            //  BRN_U->print_gate(j);
            //if(Module_ID[j]>0) BRN_Module[Module_ID[j]]->print_gate(ID_in_Module[j]);
            //if(BRN_U->Gate[j]->K!=BRN->Gate[j]->K)
            //if(Module_ID[j]>0) getchar();
        }
        BRN_U->export_Boolean_RegNet_Module_ID();
    }
    
    unsigned int get_Hamming_distance_of_states(const char s_1[],const char s_2[]){
		unsigned int i,h;
		// no setting of state in the meantime!
		h=0;
		for(i=0;i<BRN->N;i++)
			if(s_1[i]!=s_2[i]) h++;
		return(h);
	}
    
    unsigned int get_OVERLAP_of_states_IN_MODULE(int ModuleID,const char s_1[],int a,int aindex){
		unsigned int i,h,n;
        char *s_11,*s_22 = nullptr;
        
        //printf("This module has %d nodes\n",BRN_Module[ModuleID]->N);
		// no setting of state in the meantime!
		h=0;n=0;
        s_11=get_module_substate(s_1,ModuleID);
        if(BRN_Module_Dynamics!=NULL) s_22=BRN_Module_Dynamics[ModuleID]->get_attractor_state_string(a,aindex);
        else if(BRN_Module_Dynamics_TR!=NULL) s_22=BRN_Module_Dynamics_TR[ModuleID]->get_attractor_state_string(a,aindex);
        
        for(i=0;i<BRN_Module[ModuleID]->N;i++) {
            if(s_11[i]!=s_22[i]) h++;
            //printf("\t\t\tmodule state now: i=%d  %c ---  attr: i=%d  %c\n",i,s_11[i],i,s_22[i]);
        }
        delete[] s_11;s_11=NULL;
        delete[] s_22;s_22=NULL;
		return(BRN_Module[ModuleID]->N-h);
	}
    
    longint_ARRAY_pair *get_CLOSEST_Module_attractor_states(int Module_ID,const char s_1[],double *o_max){
        longint_ARRAY_pair *BestMatches;
        int o;

        BestMatches=NULL;
        (*o_max)=-1;
        if (BRN_Module_Dynamics!=NULL) {
            for(int b=1;b<=BRN_Module_Dynamics[Module_ID]->Attractor_NR;b++) {
                for (int i=1; i<=BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[b]->N; i++) {
                    o=get_OVERLAP_of_states_IN_MODULE(Module_ID,s_1,b,i);
                   // printf("Module %d A %d i=%d\t o=%d\n",Module_ID,b,i,o);//getchar();
                    if(o>(*o_max)){
                        if(BestMatches!=NULL) {
                           // printf("Cleaned BestMatches\n");
                            delete BestMatches;BestMatches=NULL;
                        }
                        BestMatches=new longint_ARRAY_pair();
                        BestMatches->add_element(b,i);
                        (*o_max)=o;
                        //printf("\tAdded %d(%d) o_max=%lg\n",b,i,o_max);
                    }
                    else {if(o==(*o_max)) {
                         BestMatches->add_element(b,i);
                        // printf("\tAdded %d(%d) as well\n",b,i);
                        }
                    }
                }
            }
        }
        else { if (BRN_Module_Dynamics_TR!=NULL) {
                     for(int b=1;b<=BRN_Module_Dynamics_TR[Module_ID]->Attractor_NR;b++) {
                        for (int i=1; i<=BRN_Module_Dynamics_TR[Module_ID]->Attractor_valleys[b]->N_a; i++) {
                            o=get_OVERLAP_of_states_IN_MODULE(Module_ID,s_1,b,i);
                            //printf("Module %d A %d i=%d\t o=%d\n",Module_ID,b,i,o);getchar();
                            if(o>(*o_max)){
                                if(BestMatches!=NULL) {
                                    // printf("Cleaned BestMatches\n");
                                    delete BestMatches;BestMatches=NULL;
                                }
                                BestMatches=new longint_ARRAY_pair();
                                BestMatches->add_element(b,i);
                                (*o_max)=o;
                                //printf("\tAdded %d(%d) o_max=%lg\n",b,i,o_max);
                            }
                            else {if(o==(*o_max)) {
                                BestMatches->add_element(b,i);
                                // printf("\tAdded %d(%d) as well\n",b,i);
                            }
                            }
                        }
                    }
                }
                else {printf("Module dynamics is not defined!\n");exit(1);}
        }
        (*o_max)/=(double)BRN_Module[Module_ID]->N;
        //printf("returing best overlap o=%lg\n",*o_max);
        //getchar();
        
        return(BestMatches);
    }

    unsigned int get_Module_attractor_states_by_module_attractor_membership_of_substate(int Module_ID,const char s_1[],double *o_max){
        char *s_11;
        if (BRN_Module_Dynamics!=NULL) {
            for(int b=1;b<=BRN_Module_Dynamics[Module_ID]->Attractor_NR;b++) {
                s_11=get_module_substate(s_1,Module_ID);
                BRN_Module_Dynamics[Module_ID]->set_state(s_11);
                unsigned long long substate_id=BRN_Module_Dynamics[Module_ID]->get_state_now();
                return((int)BRN_Module_Dynamics[Module_ID]->State_G->node[substate_id].ID);
            }
        }
        else {
            if (BRN_Module_Dynamics_TR!=NULL) {
                for(int b=1;b<=BRN_Module_Dynamics_TR[Module_ID]->Attractor_NR;b++) {
                    s_11=get_module_substate(s_1,Module_ID);
                    unsigned long a= BRN_Module_Dynamics_TR[Module_ID]->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(s_11);
                    return((unsigned int)a);
            }
        }
        else {printf("Module dynamics is not defined!\n");exit(1);}
        }
        return(0);
    }

    void print_attractor_matching(){
        for(int b=1;b<=Coupled->Attractor_NR;b++){
            printf("Coupled %d: chosen best Unconnected match: %d (h=%f)\n",b,CHOSEN_Unconnected_attractor_match_of_Composite_attractor[b],
                Unconnected_attractor_minHamming_of_Composite_attractor[b][0]);
            for(int i=1;i<=Unconnected_attractor_match_of_Composite_attractor[b][0];i++)
            printf("\tCoupled attractor %d loopstate %d maps best onto old %d  <h>=%lg\t Rank=%d\n",
                   b, i, Unconnected_attractor_match_of_Composite_attractor[b][i],
                   Unconnected_attractor_minHamming_of_Composite_attractor[b][i],
                   Rank_in_matching_Unconnected_attractor_for_Composite_attractor[b][i]);
        }
    }
    /*
    void Quantitate_Dynamical_Modularity(){
        short int *hmin_to_U;
        
        
        // initialize
        DynMod_Measure=new Dynamical_Modularity_Measure(Module_NR,get_NR_of_module_assigned_nodes(),
                                                        (int)(Coupled->Attractor_NR));
        for (int i=1; i<=Module_NR; i++){
            DynMod_Measure->M_DynMod[i]=new Module_Dynamical_Modularity_Measure(nr_mod_phen_from_DynMod(i),(int)Coupled->Attractor_NR);
         // printf("Module %d has %d phenotypes\n",i,DynMod_Measure->M_DynMod[i]->N_module_Phenotypes);
          // getchar();
        }

        //printf("a: Module 1 has %d phenotypes\n",DynMod_Measure->M_DynMod[1]->N_module_Phenotypes);
       // printf("a: Module 2 has %d phenotypes\n",DynMod_Measure->M_DynMod[2]->N_module_Phenotypes);
        
        // Calculate attractor modularity measure
        
        for(int b=1;b<=Coupled->Attractor_NR;b++){
            hmin_to_U=new short int[Unconnected->Attractor_NR+1];
            for(int i=1;i<=Unconnected->Attractor_NR;i++) hmin_to_U[i]=BRN->N+1;
            // find the closes Coupled cycle member to each unconnected attractor.
            for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
                if(Unconnected_attractor_minHamming_of_Composite_attractor[b][i]<hmin_to_U[Unconnected_attractor_match_of_Composite_attractor[b][i]]){
                            hmin_to_U[Unconnected_attractor_match_of_Composite_attractor[b][i]]=Unconnected_attractor_minHamming_of_Composite_attractor[b][i];
                    }
            }
            int nc=0;
            for(int i=1;i<=Unconnected->Attractor_NR;i++)
                if (hmin_to_U[i]<BRN->N) {
                    DynMod_Measure->Attractor_modularity_measure[b]+=hmin_to_U[i];
                    nc++;
                 //   printf("Unconn attr %d - Shift goes up by %d\n (%d members)\n",i, hmin_to_U[i], nc );
                    //getchar();
                }
            DynMod_Measure->Attractor_modularity_measure[b]/=(double)nc;
            DynMod_Measure->Attractor_modularity_measure[0]*=(1.-DynMod_Measure->Attractor_modularity_measure[b]/(double)DynMod_Measure->N_nodes_in_Mod);
          //  printf("AM(%d) = %lg  (N=%d)\n",b,1.-DynMod_Measure->Attractor_modularity_measure[b]/(double)DynMod_Measure->N_nodes_in_Mod,DynMod_Measure->N_nodes_in_Mod);
          //  printf("Coupled %d: attractor shift %lg \n",b,DynMod_Measure->Attractor_modularity_measure[b]);
            
            delete[] hmin_to_U;hmin_to_U=NULL;
           // printf("b: Module 1 has %d phenotypes\n",DynMod_Measure->M_DynMod[1]->N_module_Phenotypes);
           // printf("b: Module 2 has %d phenotypes\n",DynMod_Measure->M_DynMod[2]->N_module_Phenotypes);
           // getchar();
        }
        
        
     //   printf("Attractor modularity measure (AMM) = %lg \n",DynMod_Measure->Attractor_modularity_measure[0]);
      //  getchar();
        
        // Calculate switch stability measure
        double *ma_freq,*ma_freq_a,xx,*w_share ,*w_share_a;
        
        DynMod_Measure->Switch_stability_measure=0;
        for (int k=1;k<=Module_NR; k++) {
        //    printf("c: Module %d has %d phenotypes\n",k,DynMod_Measure->M_DynMod[k]->N_module_Phenotypes);
            // calc. the new weighted network for a module, in each coupled attractor.
            for(int b=1;b<=Coupled->Attractor_NR;b++){
         //       printf("\t Coupled Attr %d\n",b);
                
                w_share=Coupled->w_attractor_share_of_cycle_states(b);
                
                ma_freq=Module_Attr_spread_in_Coupled(b,k,DynMod_Measure->M_DynMod[k]->N_module_Phenotypes,w_share);
              //  xx=w_share[0]/(double)Coupled->W_total_steps;
                xx=1.;
                // ma_freq[1]=0.2 means module k phenotype 1 takes up 20% of the coupled attractor basin b
                // w_share cyles through the coupled attractor, storing w_attractor_share_of_cycle_states.
                for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++){
                    DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]=ma_freq[i];
        //            printf("\t\tPhenotype %d: new stability: %lg\n",i,DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]);
                }
                // transitions within the attractor basin itself: these are along the cycle.
                for (int j=2; j<=Coupled->Attractor_valleys[b]->N_a; j++){
                    DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Transition
                            [Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[b][j-1],k)]
                            [Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[b][j],k)]
                        +=Coupled->Trans_Basins[b][b]*w_share[j-1]*xx/w_share[0];
                  //  printf("\t\t\t\t Step %d->%d adds %lg (T_%lg * wsh %lg * xx_%lg / Wsh %lg) to %d-%d transition\n",j-1,j,Coupled->Trans_Basins[b][b]*w_share[j-1]*xx/w_share[0],Coupled->Trans_Basins[b][b],w_share[j-1],xx,w_share[0],Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[b][j-1],k),Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[b][j],k));
                }
                DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Transition
                        [Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[b][Coupled->Attractor_valleys[b]->N_a],k)]
                        [Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[b][1],k)]
                    +=Coupled->Trans_Basins[b][b]*w_share[Coupled->Attractor_valleys[b]->N_a]*xx/w_share[0];
               
               // printf("\t\t\t\t Step %ld->%d adds %lg (T_%lg * wsh %lg * xx_%lg / Wsh %lg) to %d-%d transition\n",Coupled->Attractor_valleys[b]->N_a,1,Coupled->Trans_Basins[b][b]*w_share[Coupled->Attractor_valleys[b]->N_a]*xx/w_share[0],Coupled->Trans_Basins[b][b],w_share[Coupled->Attractor_valleys[b]->N_a],xx,w_share[0],                       Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[b][1],k),Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[b][Coupled->Attractor_valleys[b]->N_a],k));
                
                // transitions our of the basin of b, these are overall transitions
               
               // for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++){
               //     for (int j=1; j<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; j++){
               //         printf("\t\t\t T_within[%d][%d]= %lg\n",i,j,DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Transition[i][j]);
               //     }
               // }
                
                for(int a=1;a<=Coupled->Attractor_NR;a++)
                    if(a!=b){
                        w_share_a=Coupled->w_attractor_share_of_cycle_states(a);
                        ma_freq_a=Module_Attr_spread_in_Coupled(a,k,DynMod_Measure->M_DynMod[k]->N_module_Phenotypes,w_share_a);
                        for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++)
                            for (int j=1; j<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; j++){
                                DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Transition[i][j]
                                    +=Coupled->Trans_Basins[b][a]*ma_freq[i]*xx*ma_freq_a[j]/w_share[0];
                                
                              //  printf("\t\t\t\t Basin transition %d->%d adds %lg (T_%lg * f_%lg * xx_%lg * f_a %lg / Wsh %lg) to %d-%d transition\n",b,a,Coupled->Trans_Basins[b][a]*ma_freq[i]*xx*ma_freq_a[j]/w_share[0],Coupled->Trans_Basins[b][a],ma_freq[i],xx,ma_freq_a[j],w_share[0],i,j);
                            }
                    }
                
              //  for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++){
               //     for (int j=1; j<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; j++){
               //         printf("\t\t\t T_overall[%d][%d]= %lg\n",i,j,DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Transition[i][j]);
               //     }
              //  }
               // getchar();
                delete[] ma_freq;ma_freq=NULL;
                delete[] ma_freq_a;ma_freq_a=NULL;
                delete[] w_share;w_share=NULL;
                delete[] w_share_a;w_share_a=NULL;
            }
            
            // calc. dynamical modularity scores for each module, as in grant (more or less. + cycles...)
            
            DynMod_Measure->M_DynMod[k]->Module_Score=1;
            
            for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++){
                double maxw_of_phen=0;
                double xx,y;
                y=Unconnected_module_phenotype_weight_scaled(k,i);
          //      printf("\t Module %d phenotype %d has weight %lg in the uncoupled system:\n",k,i,y);
                
                for(int b=1;b<=Coupled->Attractor_NR;b++){
                    xx=Coupled->Attractor_valleys[b]->W_Total/(double)Coupled->W_total_steps;
       //             printf("\t\t In coupled attr %d (weight %lg): it's stability is %lg\n",b,xx,DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]);
        //            printf("\t\t\t Conclusion to maximize: %lg (max= %lg):\n",xx*DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]/y,maxw_of_phen);
                    if(maxw_of_phen<xx*DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]/y)
                       maxw_of_phen=xx*DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]/y;
                }
                DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i]=maxw_of_phen;
         //       printf("\t\tPhenotype %d:  Dynmod phenotype score = %lg\n",i,DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i]);
                DynMod_Measure->M_DynMod[k]->Module_Score*=DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i];
                // getchar();
            }
            DynMod_Measure->Switch_stability_measure+=DynMod_Measure->M_DynMod[k]->Module_Score;
            
        //    printf("\t Module %d Dynmod score = %lg\n",k,DynMod_Measure->M_DynMod[k]->Module_Score);
        //    getchar();
            
     
             // calc. the overall change
            member_U=new longint_ARRAY_list[DynMod_Measure->M_DynMod[k]->N_module_Phenotypes+1];
            weight_C=new doublepointer[DynMod_Measure->M_DynMod[k]->N_module_Phenotypes+1];
            
            for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++) {
                member_U[i]=new longint_ARRAY();
                weight_C[i]=new double[Coupled->Attractor_NR+1];
            }
            
            for(int b=1;b<=Unconnected->Attractor_NR;b++){
                int mphen=Module_Attr_in_Uncoupled(b,k);
                member_U[mphen]->add_element(b);
            }
            for(int b=1;b<=Coupled->Attractor_NR;b++){
                double *ma_freq;
                ma_freq=Module_Attr_spread_in_Coupled(b,k,DynMod_Measure->M_DynMod[k]->N_module_Phenotypes);
                    // ma_freq[1]=0.2 means module k phenotype 1 takes up 20% of the coupled attractor cycle b
                for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++) 
                    weight_C[i][b]=ma_freq[i];
            }
            
            printf("\nModule %d",k);
            for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++){
                printf("\n\tPhenotype %d\n\t\tUncoupled: ",i);
                for (int l=1; l<=member_U[i]->N; l++) {
                    printf("Attr-%ld   ",member_U[i]->A[l]);
                }
                printf("\n\t\tCoupled: ");
                for(int b=1;b<=Coupled->Attractor_NR;b++) {
                    printf("Attr %d--%lg   ",b,weight_C[i][b]);
                }
            }
    
            double wu=0;double wc=0;
            
            for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++){
                wu=0; wc=0;
                for (int l=1; l<=member_U[i]->N; l++) {
                    wu+=Unconnected->Attractor_valleys[member_U[i]->A[l]]->W_Total/((double)(Unconnected->W_total_steps));
                    printf("\n\t\t\t l=%d U - Valley %ld w=%lld - %lg",l,member_U[i]->A[l],Unconnected->Attractor_valleys[member_U[i]->A[l]]->W_Total,Unconnected->Attractor_valleys[member_U[i]->A[l]]->W_Total/((double)(Unconnected->W_total_steps)));
                }
                for(int b=1;b<=Coupled->Attractor_NR;b++) {
                    wc+=weight_C[i][b]*Coupled->Attractor_valleys[b]->W_Total/((double)Coupled->W_total_steps);
                }
                DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i][i]=wc/wu;
                printf("\n\tPhenotype %d\n\t\tPhenotype_score_diagonal = %lg ",i,DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i][i]);
                
                for (int j=1; j<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; j++){
                    wu=0; wc=0;
                    for (int l=1; l<=member_U[i]->N; l++) {
                        for (int l2=1; l2<=member_U[j]->N; l2++) {
                            wu+=Unconnected->Trans_Basins[member_U[i]->A[l]][member_U[j]->A[l2]]/((double)Unconnected->W_total_steps);
                        }
                    }
                    for(int b=1;b<=Coupled->Attractor_NR;b++)
                        for(int b2=1;b2<=Coupled->Attractor_NR;b2++) {
                            wc+=weight_C[i][b]*weight_C[j][b2]*Coupled->Trans_Basins[b][b2]/((double)Coupled->W_total_steps);
                    }
                    if(wu>0)
                        DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i][j]=wc/wu;
                    else DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i][j]=0;
                    printf("n\t\tPhenotype_score with %d = %lg ",j,DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i][j]);
                }
            }
     
           // getchar();
        }
     //   printf("\t Switch stability measure (SSM) = %lg\n",DynMod_Measure->Switch_stability_measure);
    //    getchar();
        
        // calculate Module coordination measure
    //    DynMod_Measure->Module_coordination_measure=(Unconnected->Attractor_NR - Coupled->Attractor_NR)/(double)Unconnected->Attractor_NR;
    //    printf("\t Module coordination measure (MCM) = %lg\n",DynMod_Measure->Module_coordination_measure);
        
       
        printf("SUMMARY:\n\tAttractor modularity measure (AMM) = %lg \n",DynMod_Measure->Attractor_modularity_measure[0]);
        printf("\t Switch stability measure (SSM) = %lg\n",DynMod_Measure->Switch_stability_measure);
        printf("\t Module coordination measure (MCM) = %lg\n",DynMod_Measure->Module_coordination_measure);
        
    }
    
    */
    
    
    tripleunsignedintpointer *old_and_new_attractors_match_matrix(int Module_ID,Module_Attractors_in_a_row *UA_lineup){
        tripleunsignedintpointer *New_ROW_Old_COL_MATCH;
        
        int o,maxo;
        
        New_ROW_Old_COL_MATCH =new tripleunsignedintpointer[Coupled->Attractor_NR+1];
        for(int b=1;b<=Coupled->Attractor_NR;b++){
            New_ROW_Old_COL_MATCH[b]=new unsignedintpointer[Coupled->Attractor_valleys[b]->N_a+1];
            for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
                New_ROW_Old_COL_MATCH[b][i]=new unsigned int[UA_lineup->N_UAS+1];
                for(int k=1;k<=UA_lineup->N_UAS;k++){
                    
                    BRN_Module_Dynamics[Module_ID]->set_state(BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[UA_lineup->MA_Membership[k]]->A[UA_lineup->MAttractor_index[k]]);
                    o=get_OVERLAP_of_states_IN_MODULE(Module_ID,Coupled->Attractor_valleys[b]->Basin_Stored[Coupled->Attractor_valleys[b]->Attractor[i]]->s,
                                                      UA_lineup->MA_Membership[k],UA_lineup->MAttractor_index[k]);
                    New_ROW_Old_COL_MATCH[b][i][k]=o;
                   // printf("Overlap of C-%d state %d with M-%d A-%d state %d = %d\n",b,i,Module_ID,UA_lineup->MA_Membership[k],UA_lineup->MAttractor_index[k],New_ROW_Old_COL_MATCH[b][i][k]);//getchar();
                }
            }
        }
        
        for(int b=1;b<=Coupled->Attractor_NR;b++){
            for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){ // for each row: meaning the points along the bth coupled attractor
                maxo=maximum(UA_lineup->N_UAS,New_ROW_Old_COL_MATCH[b][i]);
                for(int k=1;k<=UA_lineup->N_UAS;k++){
                    if(New_ROW_Old_COL_MATCH[b][i][k]<maxo) New_ROW_Old_COL_MATCH[b][i][k]=0;
                    else {
                        printf("\t  Max in row: C-%d state %d with M-%d A-%d state %d = %d\n",b,i,Module_ID,UA_lineup->MA_Membership[k],UA_lineup->MAttractor_index[k],New_ROW_Old_COL_MATCH[b][i][k]);//getchar();
                    }
                }
            }
        }
        //getchar();
        return(New_ROW_Old_COL_MATCH);
    }
    
    tripleunsignedintpointer *old_and_Uncoupled_attractors_match_matrix(int Module_ID,Module_Attractors_in_a_row *UA_lineup){
        tripleunsignedintpointer *Uncoupled_ROW_Old_COL_MATCH;
        
        int o,maxo;
        
        Uncoupled_ROW_Old_COL_MATCH =new tripleunsignedintpointer[Unconnected->Attractor_NR+1];
        for(int b=1;b<=Unconnected->Attractor_NR;b++){
            Uncoupled_ROW_Old_COL_MATCH[b]=new unsignedintpointer[Unconnected->Attractor_valleys[b]->N_a+1];
            for(int i=1;i<=Unconnected->Attractor_valleys[b]->N_a;i++){
                Uncoupled_ROW_Old_COL_MATCH[b][i]=new unsigned int[UA_lineup->N_UAS+1];
                for(int k=1;k<=UA_lineup->N_UAS;k++){
                    BRN_Module_Dynamics[Module_ID]->set_state(BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[UA_lineup->MA_Membership[k]]->A[UA_lineup->MAttractor_index[k]]);
                    o=get_OVERLAP_of_states_IN_MODULE(Module_ID,Unconnected->Attractor_valleys[b]->Basin_Stored[Unconnected->Attractor_valleys[b]->Attractor[i]]->s,
                                                      UA_lineup->MA_Membership[k],UA_lineup->MAttractor_index[k]);
                    Uncoupled_ROW_Old_COL_MATCH[b][i][k]=o;
                    //printf("Overlap of C-%d state %d with M-%d A-%d state %d = %d\n",b,i,Module_ID,UA_lineup->MA_Membership[k],UA_lineup->MAttractor_index[k],Uncoupled_ROW_Old_COL_MATCH[b][i][k]);//getchar();
                }
            }
        }
        
        for(int b=1;b<=Unconnected->Attractor_NR;b++){
            for(int i=1;i<=Unconnected->Attractor_valleys[b]->N_a;i++){ // for each row: meaning the points along the bth coupled attractor
                maxo=maximum(UA_lineup->N_UAS,Uncoupled_ROW_Old_COL_MATCH[b][i]);
                for(int k=1;k<=UA_lineup->N_UAS;k++){
                    if(Uncoupled_ROW_Old_COL_MATCH[b][i][k]<maxo) Uncoupled_ROW_Old_COL_MATCH[b][i][k]=0;
                    //else {
                      //  printf("\t  Max in row: C-%d state %d with M-%d A-%d state %d = %d\n",b,i,Module_ID,UA_lineup->MA_Membership[k],UA_lineup->MAttractor_index[k],Uncoupled_ROW_Old_COL_MATCH[b][i][k]);//getchar();
                    //}
                }
            }
        }
        //getchar();
        return(Uncoupled_ROW_Old_COL_MATCH);
    }
    
    tripleunsignedintpointer *old_and_new_attractors_closest_flyby_matrix(Module_Attractors_in_a_row *UA_lineup,tripleunsignedintpointer *New_ROW_Old_COL_MATCH){
        tripleunsignedintpointer *New_ROW_Old_ROW_FLyBY_MATCH;
        int maxo;
        unsigned int *x;
        
        printf("c. Separate matrix into horizontal strips: separate bunch of rows for the states of each coupled attractor:\nc.1 -> keep nonzero MAXIMUM in each COLUMN within the strip:\n");
        
        New_ROW_Old_ROW_FLyBY_MATCH =new tripleunsignedintpointer[Coupled->Attractor_NR+1];
        for(int b=1;b<=Coupled->Attractor_NR;b++){
            New_ROW_Old_ROW_FLyBY_MATCH[b]=new unsignedintpointer[Coupled->Attractor_valleys[b]->N_a+1];
            for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
                New_ROW_Old_ROW_FLyBY_MATCH[b][i]=new unsigned int[UA_lineup->N_UAS+1];
                for(int k=1;k<=UA_lineup->N_UAS;k++){
                    New_ROW_Old_ROW_FLyBY_MATCH[b][i][k]=New_ROW_Old_COL_MATCH[b][i][k];
                }
            }
        }
        
        for(int b=1;b<=Coupled->Attractor_NR;b++){
            printf("     Coupled attractor C-%d:\n",b);
            for(int k=1;k<=UA_lineup->N_UAS;k++){
                // for each column: meaning the points along uncoupled attractors, all in a row
                x=new unsigned int[Coupled->Attractor_valleys[b]->N_a+1];
                for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
                    x[i]=New_ROW_Old_ROW_FLyBY_MATCH[b][i][k];
                }
                maxo=maximum(Coupled->Attractor_valleys[b]->N_a,x);
                if (maxo>0) {printf("\t  - col %d (part of Module A-%d):\n",k,UA_lineup->MA_Membership[k]);}
                for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
                    if(New_ROW_Old_ROW_FLyBY_MATCH[b][i][k]<maxo) New_ROW_Old_ROW_FLyBY_MATCH[b][i][k]=0;
                    else if(New_ROW_Old_ROW_FLyBY_MATCH[b][i][k]>0){
                        printf("\t\t\tC-%d state %d with A-%d state %d = %d\n",b,i,UA_lineup->MA_Membership[k],UA_lineup->MAttractor_index[k],New_ROW_Old_ROW_FLyBY_MATCH[b][i][k]);//getchar();
                    }
                }
            }
        }
        
        return(New_ROW_Old_ROW_FLyBY_MATCH);
    }
   
   Module_Attractors_in_a_row *get_UA_lineup(int Module_ID){
       Module_Attractors_in_a_row *UA_lineup;
       
       unsigned int n=0;
       for(int b=1;b<=BRN_Module_Dynamics[Module_ID]->Attractor_NR;b++) {
            n+=BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[b]->N;
       }
       UA_lineup=new Module_Attractors_in_a_row(n);
       n=0;
       for(int b=1;b<=BRN_Module_Dynamics[Module_ID]->Attractor_NR;b++) {
           for (int i=1; i<=BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[b]->N; i++) {
               n++;
               UA_lineup->MA_Membership[n]=b;
               UA_lineup->MAttractor_index[n]=i;
           }
       }
       return(UA_lineup);
   }
    
    void delete_Overlap_matrix(tripleunsignedintpointer *ooo){
        for(int b=1;b<=Coupled->Attractor_NR;b++){
            for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){ delete[] ooo[b][i];ooo[b][i]=NULL;}
            delete[] ooo[b];ooo[b]=NULL;
        }
        delete[] ooo;ooo=NULL;
    }
    
    double Calculate_cycle_to_cycle_NORM_overlap(int Module_ID,
                                                 int a, // module attractor ID
                                                 int b, // Coupled Attractor basin ID
                                                 longint_ARRAY *Coupled_cycle_on_flyby, //ordered list if Coupled attractor states on closest flyby to Module attractor a
                                                 int Couple_states_are_a_loop,
                                                 tripleunsignedintpointer *New_ROW_Old_ROW_FLyBY_MATCH,
                                                 Module_Attractors_in_a_row *UA_lineup){
        double weighted_overlap=0;
        longint_ARRAY_list *Module_cycle_index_match,*Coupled_cycle_k_match;
        unsigned long int delta;
        unsigned int *Module_Cycle_transition_Mimicked;
        
        Module_cycle_index_match=new longint_ARRAY_list[Coupled_cycle_on_flyby->N+1];
        Coupled_cycle_k_match=new longint_ARRAY_list[Coupled_cycle_on_flyby->N+1];
        
        Module_Cycle_transition_Mimicked= new unsigned int[BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N+1];
        for (int i=0; i<=BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N; i++) Module_Cycle_transition_Mimicked[i]=(int)BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N;
        
        for (int i=1; i<=Coupled_cycle_on_flyby->N; i++) {
            Coupled_cycle_k_match[i]=new longint_ARRAY();;
            Module_cycle_index_match[i]=new longint_ARRAY();
            for(int k=1;k<=UA_lineup->N_UAS;k++){
                if((UA_lineup->MA_Membership[k]==a) && (New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[i]][k]>0)) {
                    Module_cycle_index_match[i]->add_element(UA_lineup->MAttractor_index[k]); // matching Moduleattractor-state index of the ith flyby coupled state
                    Coupled_cycle_k_match[i]->add_element(k);  // matching column of the ith flyby coupled state
                    // this is where douplicate matches get dropped:
                     // if a coupled cycle state is equally close to 2 module cycle states (and these are in the New_ROW_Old_ROW_FLyBY_MATCH matrix) then the condition above is true twice.
                    // if this happens, the second choice gets auto-overwritten in Module_cycle_index_match and Coupled_cycle_k_match
                }
            }
            if( Coupled_cycle_k_match[i]==0){
                printf("Problem this should never happen!\n");getchar();
            }
           // printf("\t\tCoupled cycle on flyby %d = %ld\n",i,Coupled_cycle_on_flyby->A[i]);
        }
        
        for (int i=1; i<Coupled_cycle_on_flyby->N; i++) {
            double w=0;
            for (int l2=1; l2<=Module_cycle_index_match[i+1]->N;l2++) {
                for (int l1=1; l1<=Module_cycle_index_match[i]->N;l1++) {
                    if (Module_cycle_index_match[i+1]->A[l2]>=Module_cycle_index_match[i]->A[l1]) {
                        delta=Module_cycle_index_match[i+1]->A[l2]-Module_cycle_index_match[i]->A[l1];
                        for (unsigned long int msid=Module_cycle_index_match[i]->A[l1]; msid<Module_cycle_index_match[i+1]->A[l2]; msid++) {
                            if (Module_Cycle_transition_Mimicked[msid]>delta) Module_Cycle_transition_Mimicked[msid]=(int)delta;
                        }
                    }
                    else {
                        delta=Module_cycle_index_match[i+1]->A[l2]-Module_cycle_index_match[i]->A[l1]+(int)BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N;
                        for (unsigned long int msid=Module_cycle_index_match[i]->A[l1]; msid<BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N; msid++) {
                            if (Module_Cycle_transition_Mimicked[msid]>delta) Module_Cycle_transition_Mimicked[msid]=(int)delta;
                        }
                        if (Module_Cycle_transition_Mimicked[BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N]>delta)
                            Module_Cycle_transition_Mimicked[BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N]=(int)delta;
                        for (unsigned long int msid=1; msid<Module_cycle_index_match[i+1]->A[l2]; msid++) {
                            if (Module_Cycle_transition_Mimicked[msid]>delta) Module_Cycle_transition_Mimicked[msid]=(int)delta;
                        }
                    }
                    double wnow;
                    if (BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N>1) {
                        double oo=( New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[i]][Coupled_cycle_k_match[i]->A[l1]]
                                   +New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[i+1]][Coupled_cycle_k_match[i+1]->A[l2]])
                        /((double)(2.*BRN_Module[Module_ID]->N));
                        
                       // printf("oo_start=%lg\n",oo);//getchar();
                        wnow= oo * exp(-std::abs(delta-1.)/(double)(BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N-1.));
                    }
                    else {
                        wnow=( New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[i]][Coupled_cycle_k_match[i]->A[l1]]
                              +New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[i+1]][Coupled_cycle_k_match[i+1]->A[l2]])
                                    /((double)(2.*BRN_Module[Module_ID]->N));
                       // printf("wnow_start=%lg\n",wnow);//getchar();
                    }
                    printf("\t\ti:%d->%d (subset %d %d) CA transition %ld (o-%d) -> %ld (o-%d): delta=%lu   => w = %lg   ",
                           i,i+1,l1,l2,Coupled_cycle_on_flyby->A[i],
                           New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[i]][Coupled_cycle_k_match[i]->A[l1]],
                           Coupled_cycle_on_flyby->A[i+1],
                           New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[i+1]][Coupled_cycle_k_match[i+1]->A[l2]],
                           delta,wnow);
                  
                    printf("wnow=%lg\n",wnow);//getchar();
                    if (wnow<0) wnow=0.;
                    w+=wnow;
                }
            }
            w/=(double)((Module_cycle_index_match[i]->N)*(Module_cycle_index_match[i+1]->N));
            weighted_overlap+=w;
           // printf("Weighted overlap = %lg\n",w);
        }
        
        // leftover last to first
        if ((Couple_states_are_a_loop)||(Coupled_cycle_on_flyby->N==1)){
            double w=0;
            for (int l2=1; l2<=Module_cycle_index_match[1]->N;l2++) {
                for (int l1=1; l1<=Module_cycle_index_match[Coupled_cycle_on_flyby->N]->N;l1++) {
                    if (Module_cycle_index_match[1]->A[l2]>=Module_cycle_index_match[Coupled_cycle_on_flyby->N]->A[l1]) {
                        delta=Module_cycle_index_match[1]->A[l2]-Module_cycle_index_match[Coupled_cycle_on_flyby->N]->A[l1];
                        for (unsigned long int msid=Module_cycle_index_match[Coupled_cycle_on_flyby->N]->A[l1]; msid<Module_cycle_index_match[1]->A[l2]; msid++) {
                            if (Module_Cycle_transition_Mimicked[msid]>delta) Module_Cycle_transition_Mimicked[msid]=(int)delta;
                        }
                    }
                    else {
                        delta=Module_cycle_index_match[1]->A[l2]-Module_cycle_index_match[Coupled_cycle_on_flyby->N]->A[l1]+(int)BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N;
                        for (unsigned long int msid=Module_cycle_index_match[Coupled_cycle_on_flyby->N]->A[l1]; msid<BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N; msid++) {
                            if (Module_Cycle_transition_Mimicked[msid]>delta) Module_Cycle_transition_Mimicked[msid]=(int)delta;
                        }
                        if (Module_Cycle_transition_Mimicked[BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N]>delta)
                            Module_Cycle_transition_Mimicked[BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N]=(int)delta;
                        for (unsigned long int msid=1; msid<Module_cycle_index_match[1]->A[l2]; msid++) {
                            if (Module_Cycle_transition_Mimicked[msid]>delta) Module_Cycle_transition_Mimicked[msid]=(int)delta;
                        }
                    }
                    double wnow;
                    if (BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N>1) {
                        double oo=( New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[Coupled_cycle_on_flyby->N]][Coupled_cycle_k_match[Coupled_cycle_on_flyby->N]->A[l1]]
                                   + New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[1]][Coupled_cycle_k_match[1]->A[l2]])
                        /((double)(2.*BRN_Module[Module_ID]->N));
                        wnow=oo * exp(-std::abs(delta-1.)/(double)(BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N-1.));
                    }
                    else{
                        wnow=( New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[Coupled_cycle_on_flyby->N]][Coupled_cycle_k_match[Coupled_cycle_on_flyby->N]->A[l1]]
                                      +New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[1]][Coupled_cycle_k_match[1]->A[l2]])
                        /((double)(2.*BRN_Module[Module_ID]->N));
                    }
                    printf("\t\ti:%lu->%d (subset %d %d) CA transition %ld (o-%d) -> %ld (o-%d): delta=%lu   => w = %lg   ",Coupled_cycle_on_flyby->N,1,l1,l2,Coupled_cycle_on_flyby->A[Coupled_cycle_on_flyby->N],
                         New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[Coupled_cycle_on_flyby->N]][Coupled_cycle_k_match[Coupled_cycle_on_flyby->N]->A[l1]],
                           Coupled_cycle_on_flyby->A[1],
                           New_ROW_Old_ROW_FLyBY_MATCH[b][Coupled_cycle_on_flyby->A[1]][Coupled_cycle_k_match[1]->A[l2]],
                           delta,wnow);
                    
                    if (wnow<0) wnow=0.;
                    w+=wnow;
                }
            }
            w/=(double)((Module_cycle_index_match[Coupled_cycle_on_flyby->N]->N)*(Module_cycle_index_match[1]->N));
            weighted_overlap+=w;
            printf("Weighted overlap = %lg\n",w);
        }
       
        if ((Couple_states_are_a_loop)||(Coupled_cycle_on_flyby->N==1))
             weighted_overlap/=(double)Coupled_cycle_on_flyby->N;
        else weighted_overlap/=(double)(Coupled_cycle_on_flyby->N-1);
        
        printf("\t\t Average weigthed overlap = %lg\n",weighted_overlap);
        
        if (BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N>1) {
            int penalty=0;
            for (int i=1; i<=BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N; i++){
                if(Module_Cycle_transition_Mimicked[i]>penalty) penalty=Module_Cycle_transition_Mimicked[i];
             //   if (i<BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N) {
             //       printf("score of module cycle step %d -> %d is %d\n",i,i+1,Module_Cycle_transition_Mimicked[i]);
             //   }
             //   else printf("score of module cycle step %d -> %d is %d\n",i,1,Module_Cycle_transition_Mimicked[i]);
            }
          //  printf("\t\t worst jump = %d\n",penalty);
          //  printf("\t\t penalty for worst jump over Module cycle steps: %lg\n",(BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N-penalty)/(double(BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N-1.)));
            weighted_overlap*=(BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N-penalty)
                              /(double(BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N-1.));
          //  if((BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N-penalty)/(double(BRN_Module_Dynamics[Module_ID]->Attractor_Cycle[a]->N-1.)) > 0) getchar();
        }
        printf("\t\t\t Final result for cycle comparisons = %lg\n",weighted_overlap);
        
        
        delete[] Module_cycle_index_match;Module_cycle_index_match=NULL;
        delete[] Coupled_cycle_k_match;Coupled_cycle_k_match=NULL;
        return(weighted_overlap);
    }
    
    longint_ARRAY_list *break_Coupled_cycle_on_flyby_pieces_into_chunks(int b,longint_ARRAY *Coupled_cycle_on_flyby,
                                                                        longint_ARRAY *Coupled_cycle_NOT_closer_to_other_Module_attractor,int *n){
        longint_ARRAY_list *Coupled_cycle_on_flyby_pieces;
        
        int talk=0;
       // printf("Coupled length = %ld\n",Coupled->Attractor_valleys[b]->N_a);
        if(talk){
            printf("RED list:  ");        Coupled_cycle_on_flyby->print();
            printf("Gray list:  ");        Coupled_cycle_NOT_closer_to_other_Module_attractor->print();
        }
        int stgap=1;
        (*n)=0; // starting on gap, no pieces.
        for (int i=1; i<=Coupled->Attractor_valleys[b]->N_a;i++) {
            if (Coupled_cycle_NOT_closer_to_other_Module_attractor->check_for_element(i)==0) {
                if(stgap==0) {stgap=1;
                   //printf("\t\trestarging a gap\n");
                }
            }
            else {
                if(stgap==1) {stgap=0;(*n)++;
                    // ending gap, starging a piece, n goes up.
                }
            }
        }
        
        Coupled_cycle_on_flyby_pieces=new longint_ARRAY_list[(*n)+1];
        stgap=1;
        int n2=0; // starting on gap, no pieces.
        for (int i=1; i<=Coupled->Attractor_valleys[b]->N_a;i++) {
            if(talk) printf("%d\n",i);
            if (Coupled_cycle_NOT_closer_to_other_Module_attractor->check_for_element(i)==0) {
                if(stgap==0) {stgap=1;
                    printf("\t\t restarging a gap\n");
                    if(Coupled_cycle_on_flyby_pieces[n2]->N==0) {
                       n2--;
                       if(talk) printf("\t\t previous piece had no elements, removing it, n2 back to %d \n",n2);
                    }
                }
            }
            else {
                if(stgap==1) {
                    stgap=0;n2++;
                    if(talk) printf("\t\t ending gap, starging a piece, n goes up to %d\n",n2); // ending gap, starging a piece, n goes up.
                    
                    Coupled_cycle_on_flyby_pieces[n2]=new longint_ARRAY();
                    if (Coupled_cycle_on_flyby->check_for_element(i)>0) {
                        Coupled_cycle_on_flyby_pieces[n2]->add_element(i);
                        if(talk) printf("---Piece %d accruing %d\n",n2,i);
                    }
                    else {if(talk) printf("----------Piece %d: state %d not a red link\n",n2,i);}
                }
                else {
                    if (Coupled_cycle_on_flyby->check_for_element(i)>0) {
                        Coupled_cycle_on_flyby_pieces[n2]->add_element(i);
                        if(talk)printf("---Piece %d accruing %d\n",n2,i);
                    }
                    else {if(talk) printf("----------Piece %d: state %d not a red link\n",n2,i);}
                }
            }
        }

        
        if((Coupled_cycle_on_flyby_pieces[n2]->N>0)&&(Coupled_cycle_NOT_closer_to_other_Module_attractor->check_for_element(1)>0)&&(Coupled_cycle_NOT_closer_to_other_Module_attractor->check_for_element(Coupled->Attractor_valleys[b]->N_a)>0)) {
            // the first and last piece are joined here.
            if(n2>1){
                Coupled_cycle_on_flyby_pieces[n2]->add_longint_ARRAY(Coupled_cycle_on_flyby_pieces[1]);
                delete Coupled_cycle_on_flyby_pieces[1];
                Coupled_cycle_on_flyby_pieces[1]=Coupled_cycle_on_flyby_pieces[n2];
                n2--;
                if(talk) printf("\t\t uniting last and first, n goes down to %d\n",n2); // ending gap, starging a piece, n goes up.
            }
        }
        
        if(Coupled_cycle_on_flyby_pieces[n2]->N==0) {
            n2--;
            if(talk) printf("\t\t last piece had no elements, removing it, n2 back to %d \n",n2);
        }
        
        (*n)=n2;
        if(talk){
            printf("\t\tpieces=%d\n",(*n));
            for (int i=1; i<=(*n); i++) {
                printf("\t\t piece %d:  ",i); Coupled_cycle_on_flyby_pieces[i]->print();
            }
            getchar();
        }
        return(Coupled_cycle_on_flyby_pieces);
    }
    

    double get_best_Normalized_overlap_for_Coupled_Attractor_in_Module(int Module_ID,
                                                                       int b,
                                                                       tripleunsignedintpointer *New_ROW_Old_ROW_FLyBY_MATCH,
                                                                       tripleunsignedintpointer *New_ROW_Old_COL_MATCH,
                                                                       Module_Attractors_in_a_row *UA_lineup,
                                                                       int *nonzero_score_for_Module_Attr){
        longint_ARRAY *Coupled_cycle_on_flyby,*Coupled_cycle_NOT_closer_to_other_Module_attractor;
        unsigned int *ModuleAttractor_touched_by_flyby;
        
        ModuleAttractor_touched_by_flyby=new unsigned int [BRN_Module_Dynamics[Module_ID]->Attractor_NR+1];
        
        for(int a=0;a<=BRN_Module_Dynamics[Module_ID]->Attractor_NR;a++)   ModuleAttractor_touched_by_flyby[a]=0;
        
        for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
            for(int k=1;k<=UA_lineup->N_UAS;k++){
                if(New_ROW_Old_ROW_FLyBY_MATCH[b][i][k]>0) {
                    if(ModuleAttractor_touched_by_flyby[UA_lineup->MA_Membership[k]]==0) {
                        ModuleAttractor_touched_by_flyby[UA_lineup->MA_Membership[k]]=1;
                        ModuleAttractor_touched_by_flyby[0]++;
                    }
                }
            }
        }
        
        double av_FlyBy_Overlap=0;
        for(int a=1;a<=BRN_Module_Dynamics[Module_ID]->Attractor_NR;a++){
            if (ModuleAttractor_touched_by_flyby[a]==1) {
                printf("\t- flyby module A-%d\n\t\t",a);
                Coupled_cycle_on_flyby=new longint_ARRAY();
                Coupled_cycle_NOT_closer_to_other_Module_attractor=new longint_ARRAY();
                
                for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
                    for(int k=1;k<=UA_lineup->N_UAS;k++){
                        if(UA_lineup->MA_Membership[k]==a){
                           if (New_ROW_Old_ROW_FLyBY_MATCH[b][i][k]>0) {
                               Coupled_cycle_on_flyby->add_element(i);
                               Coupled_cycle_NOT_closer_to_other_Module_attractor->add_element(i);
                           }
                           else{ if(New_ROW_Old_COL_MATCH[b][i][k]>0)
                                Coupled_cycle_NOT_closer_to_other_Module_attractor->add_element(i);
                           }
                        }
                    }
                }
                Coupled_cycle_on_flyby->print(); //getchar();
                
                if (Coupled_cycle_NOT_closer_to_other_Module_attractor->N==Coupled->Attractor_valleys[b]->N_a) {
                    // this means the global attractor stays in the same module attractor the whole type: counts as CLOSED loop:
                    double F_Overlap=Calculate_cycle_to_cycle_NORM_overlap(Module_ID,a,b,Coupled_cycle_on_flyby,1,New_ROW_Old_ROW_FLyBY_MATCH,UA_lineup);
                    if(F_Overlap>0) nonzero_score_for_Module_Attr[a]=1; else nonzero_score_for_Module_Attr[a]=0;
                    av_FlyBy_Overlap+=F_Overlap;
                }
                else{
                    int nflyby_pieces=0;
                    longint_ARRAY_list *Coupled_cycle_on_flyby_pieces;
                    Coupled_cycle_on_flyby_pieces=break_Coupled_cycle_on_flyby_pieces_into_chunks(b,Coupled_cycle_on_flyby,Coupled_cycle_NOT_closer_to_other_Module_attractor,&nflyby_pieces);
                    double F_Overlap=0;
                    for (int i=1;i<=nflyby_pieces; i++) {
                        double w=Calculate_cycle_to_cycle_NORM_overlap(Module_ID,a,b,Coupled_cycle_on_flyby_pieces[i],0,New_ROW_Old_ROW_FLyBY_MATCH,UA_lineup);
                        F_Overlap+=w;
                    }
                    F_Overlap/=(double)nflyby_pieces;
                    if(F_Overlap>0) nonzero_score_for_Module_Attr[a]=1; else nonzero_score_for_Module_Attr[a]=0;
                    av_FlyBy_Overlap+=F_Overlap;
                    for (int i=1;i<=nflyby_pieces; i++) {delete Coupled_cycle_on_flyby_pieces[i];Coupled_cycle_on_flyby_pieces[i]=NULL;}
                    delete[] Coupled_cycle_on_flyby_pieces;Coupled_cycle_on_flyby_pieces=NULL;
                }
                delete Coupled_cycle_on_flyby;Coupled_cycle_on_flyby=NULL;
                delete Coupled_cycle_NOT_closer_to_other_Module_attractor;Coupled_cycle_NOT_closer_to_other_Module_attractor=NULL;
                printf("Sum of flybys so far = %lg\n",av_FlyBy_Overlap);
            }
            else nonzero_score_for_Module_Attr[a]=0;
        }
        if (ModuleAttractor_touched_by_flyby[0]>0)
            av_FlyBy_Overlap/=(double)ModuleAttractor_touched_by_flyby[0];
        else av_FlyBy_Overlap=0;
        
        printf("  C-%d: Average of flybys = %lg\n\n",b,av_FlyBy_Overlap);//getchar();
        return(av_FlyBy_Overlap);
    }
    
    void print_attractor_comparison(){
        char *sm;
        if (BRN_Module_Dynamics!=NULL) {
            for (int k=1;k<=Module_NR; k++){
                printf("\n\nModule %d:\n",k);
                for(int a=1;a<=BRN_Module_Dynamics[k]->BRN->N;a++)
                    //BRN_Module_Dynamics[k]->BRN->print_gate(a);
                    printf("%s\t",BRN_Module_Dynamics[k]->BRN->MDAT->Node_name[a]);
                printf("\n");
                for(int a=1;a<=BRN_Module_Dynamics[k]->Attractor_NR;a++) {
                    printf("\tM-%d A-%d (cycle length %llu):\n",k,a,BRN_Module_Dynamics[k]->Attractor_Cycle[a]->N);
                    for(int i=1;i<=BRN_Module_Dynamics[k]->Attractor_Cycle[a]->N;i++){
                        sm=BRN_Module_Dynamics[k]->get_attractor_state_string(a,i);
                        printf("\t\t %s\n",sm);
                        delete[] sm;sm=NULL;
                    }
                }
            }
        }
        else
            if (BRN_Module_Dynamics_TR!=NULL) {
                for (int k=1;k<=Module_NR; k++){
                    printf("\n\nModule %d:\n",k);
                    for(int a=1;a<=BRN_Module_Dynamics_TR[k]->BRN->N;a++)
                        //BRN_Module_Dynamics[k]->BRN->print_gate(a);
                        printf("%s\t",BRN_Module_Dynamics_TR[k]->BRN->MDAT->Node_name[a]);
                    printf("\n");
                    for(int a=1;a<=BRN_Module_Dynamics_TR[k]->Attractor_NR;a++) {
                        printf("\tM-%d A-%d (cycle length %lu):\n",k,a,BRN_Module_Dynamics_TR[k]->Attractor_valleys[a]->N_a);
                        for(int i=1;i<=BRN_Module_Dynamics_TR[k]->Attractor_valleys[a]->N_a;i++){
                            sm=BRN_Module_Dynamics_TR[k]->get_attractor_state_string(a,i);
                            printf("\t\t %s\n",sm);
                            delete[] sm;sm=NULL;
                        }
                    }
                }
            }

        printf("\n\nCoupled networks:\n");
        for(int b=1;b<=Coupled->Attractor_NR;b++){
           printf("Attractor Cycle %d with %lu states:\n",b,Coupled->Attractor_valleys[b]->N_a);
            for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
                printf(" %d - ",i);
                for (int k=0;k<=Module_NR; k++){
                    sm=get_module_substate(Coupled->Attractor_valleys[b]->Basin_Stored[Coupled->Attractor_valleys[b]->Attractor[i]]->s,k);
                    printf("%s | ",sm);
                    delete[] sm;sm=NULL;
                }
                printf("\n");
            }
        }
        //getchar();
    }
    
    double *Module_Attr_spread_in_Coupled_DAVID(int C_AID,
                                                tripleunsignedintpointer *New_ROW_Old_COL_MATCH,
                                                Module_Attractors_in_a_row *UA_lineup,
                                                int *Nonzero_score_for_Module_attractor,
                                                int n_phen,
                                                double w_share[]){
        double *w;
        
        // printf("STARTING\n");getchar();
        w=new double[n_phen+1];
        for (int i=1; i<=n_phen; i++) w[i]=0;
        
        for(int i=1;i<=Coupled->Attractor_valleys[C_AID]->N_a;i++){
            //      printf("Coupled state %d gets share %lg of coupled valley %d\n",i,w_share[i],C_AID);
            for (int k=1;k<=UA_lineup->N_UAS; k++) {
                if ((Nonzero_score_for_Module_attractor[UA_lineup->MA_Membership[k]]==1)&&( New_ROW_Old_COL_MATCH[C_AID][i][k]>0)) {
                    w[UA_lineup->MA_Membership[k]]+=w_share[i];
                   // if(Coupled->Attractor_valleys[C_AID]->N_a>1) printf("coupled i=%d: W for Module attr %d going up by %lg\n",i,UA_lineup->MA_Membership[k],w_share[i]);
                }
            }
        }
      //  for (int i=1; i<=n_phen; i++)
        //    printf("\t\t\t Within Coupled attr %d, total fraction of time spent in module phenotype %d, f[%d]=%lg\n",C_AID,i,i,w[i]);
        return(w);
    }
    
    double Unconnected_module_phenotype_weight_scaled_DAVID(tripleunsignedintpointer *Uncoupled_ROW_Old_COL_MATCH,
                                                            Module_Attractors_in_a_row *UA_lineup,
                                                            int phen_ID){
        double w=0,wt=0;
        
        for(int b=1;b<=Unconnected->Attractor_NR;b++){
            wt+=Unconnected->Attractor_valleys[b]->W_Total;
   
            int ok=0;
            for(int i=1;i<=Unconnected->Attractor_valleys[b]->N_a;i++){ // all states in uncoupled attr b
                for (int k=1;k<=UA_lineup->N_UAS; k++) { // all linedup states of module attractors for which Uncoupled_ROW_Old_COL_MATCH and UA_lineup were computed
                    if ( (Uncoupled_ROW_Old_COL_MATCH[b][i][k]>0)&&(UA_lineup->MA_Membership[k]==phen_ID)) {
                        ok=1;  // module attractor k falls into is part of the uncoupled combo  b AND the attractor number k falls into = phen_ID
                    }
                }
            }
            if(ok==1)  w+=Unconnected->Attractor_valleys[b]->W_Total;
        }
        return(w/(double)wt);
    }
    
    FILE *Quantitate_Dynamical_Modularity_DAVID(int rnd_version, unsigned long int N_trace,unsigned int N_RND,double p_error, int talk){
        Module_Attractors_in_a_row *UA_lineup;
        tripleunsignedintpointer *New_ROW_Old_COL_MATCH,*New_ROW_Old_ROW_FLyBY_MATCH,*Uncoupled_ROW_Old_COL_MATCH;
        
        BRN_Module_Dynamics=new Boolean_Dynamics_list[Module_NR+1];
        
        for (int k=1;k<=Module_NR; k++) {
            BRN_Module_Dynamics[k]=new Boolean_Dynamics(BRN_Module[k]);
            BRN_Module_Dynamics[k]->read_attractors();
        }
        
        print_attractor_comparison();// getchar();
        
        DynMod_Measure=new Dynamical_Modularity_Measure(Module_NR,get_NR_of_module_assigned_nodes(),
                                                        (int)(Coupled->Attractor_NR));
        for (int i=1; i<=Module_NR; i++){
            DynMod_Measure->M_DynMod[i]=new Module_Dynamical_Modularity_Measure((int)BRN_Module_Dynamics[i]->Attractor_NR,(int)Coupled->Attractor_NR);
      //      printf("Module %d has %d phenotypes\n",i,DynMod_Measure->M_DynMod[i]->N_module_Phenotypes);
        }
        
        DynMod_Measure->Attractor_modularity_measure[0]=1.;
        DynMod_Measure->Switch_stability_measure=1;
        double *ma_freq,*w_share;
        
        for (int k=1;k<=Module_NR; k++) {
           if(talk) printf("\n\nMODULE %d RESULTS:\nLAW 1\n\na. Calculating First Matrix:\n\t - columns: attractor states of Module %d, one after the other\n\t- rows: coupled attractor states, one after the other\n\nb. Cutting First Matrix to keep only maxima in each ROW:\n",k,k);
            
            UA_lineup=get_UA_lineup(k);
            
            // all states of all coupled attractors: best matching module attractor state
            New_ROW_Old_COL_MATCH=old_and_new_attractors_match_matrix(k,UA_lineup);
            
            // all module attractors: couple attractor states on closest flyby
            New_ROW_Old_ROW_FLyBY_MATCH=old_and_new_attractors_closest_flyby_matrix(UA_lineup,New_ROW_Old_COL_MATCH);
            
            // all states of all Unconnected attractors: best matching module attractor state
            Uncoupled_ROW_Old_COL_MATCH=old_and_Uncoupled_attractors_match_matrix(k,UA_lineup);
            
            intpointer *nonzero_score_for_Module_Attr;
            nonzero_score_for_Module_Attr=new intpointer[Coupled->Attractor_NR+1];
            
            // LAW 1

            if(talk) printf("c.2 -> for each coupled attractor, compare “closest fly-by” to module attractor (cycle vs cycle comparison):\n");
            double Obest_CA=1;
            
         //   double H_max=BRN_Module_Dynamics[k]->largest_H_from_all_attractors();
           // printf("\t\t Module %d -> H_max = %lg\n",k,H_max);//getchar();
          //  if(BRN_Module_Dynamics[k]->BRN->N>1) H_max=H_max/(double)(BRN_Module_Dynamics[k]->BRN->N);
         //   printf("\t\t Module %d -> H_max = %lg\n",k,H_max);//getchar();
            
            for(int b=1;b<=Coupled->Attractor_NR;b++){
                nonzero_score_for_Module_Attr[b]=new int[DynMod_Measure->M_DynMod[k]->N_module_Phenotypes+1];
                if(talk) printf("     Coupled attractor C-%d:\n",b);
                // averages over best flyby close to uncoupled attractors this global cycle touches
                double ooo=get_best_Normalized_overlap_for_Coupled_Attractor_in_Module(k,b,New_ROW_Old_ROW_FLyBY_MATCH,New_ROW_Old_COL_MATCH,UA_lineup,nonzero_score_for_Module_Attr[b]);
              
                ooo=2*(ooo-0.5);
                if (ooo<0) ooo=0;
                Obest_CA=Obest_CA*ooo;
                //printf("     ooo=%lg   Obest_CA=%lg:\n",ooo,Obest_CA);getchar();
            }
            
           // DynMod_Measure->Attractor_modularity_measure[k]=pow(Obest_CA,(double)(1./(double)DynMod_Measure->M_DynMod[k]->N_module_Phenotypes));
            
            
            DynMod_Measure->Attractor_modularity_measure[k]=pow(Obest_CA,(double)(1./(double)Coupled->Attractor_NR));
            
            //old
            DynMod_Measure->Attractor_modularity_measure[k]=pow(Obest_CA,(double)(1./(double)Coupled->Attractor_NR));
            DynMod_Measure->Attractor_modularity_measure[0]*=DynMod_Measure->Attractor_modularity_measure[k];
            
            if(talk) printf("Attractor modularity measure for Module %d  (%d phenotypes)(AMM[%d]) = %lg \n",k,DynMod_Measure->M_DynMod[k]->N_module_Phenotypes,k,DynMod_Measure->Attractor_modularity_measure[k]);
           
            // getchar();
       
            
            // LAW 2

           
            if(talk) printf("\n\n LAW 2\n\nc: Module %d has %d phenotypes\n",k,DynMod_Measure->M_DynMod[k]->N_module_Phenotypes);
          //  getchar();
            // calc. the new weighted network for a module, in each coupled attractor.
            
            for(int b=1;b<=Coupled->Attractor_NR;b++){
               if(talk)  printf("\t Coupled Attr %d\n",b);
                
                w_share=Coupled->w_attractor_share_of_cycle_states(b);  // for every state in C. attr. b, the sum of weights that flows in
                double xx=w_share[0]/(double)Coupled->W_total_steps; // total prob of basin b
                printf("\t\t Probability of the system being anywhere in the entire attractor basin %d, W[%d] =%lg\n",b,b,xx);//getchar();
                
                if(talk) { for (int iiii=1; iiii<=Coupled->Attractor_valleys[b]->N_a; iiii++)
                    printf("\t\t In Coupled attr %d, state %d: total probl: w[%d][%d] = %lg\n",b, iiii,b, iiii,w_share[iiii]*xx);
                }
                
                ma_freq=Module_Attr_spread_in_Coupled_DAVID(b,New_ROW_Old_COL_MATCH,UA_lineup,
                                                            nonzero_score_for_Module_Attr[b],
                                                            DynMod_Measure->M_DynMod[k]->N_module_Phenotypes,w_share);
                // xx=1.;
                // ma_freq[1]=0.2 means module k phenotype 1 takes up 20% of the coupled attractor basin b
                // w_share cyles through the coupled attractor, storing w_attractor_share_of_cycle_states.
                for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++){
                    DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]=ma_freq[i]*xx; //normalized by weight of global attr b
                    if(talk) printf("\t\tPhenotype %d: stability score = f[%d]*W[%d] = %lg\n",i,i,b,DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]);//getchar();
                }
                delete[] ma_freq;ma_freq=NULL;
                delete[] w_share;w_share=NULL;
            }
            
            // calc. dynamical modularity scores for each module
            
            DynMod_Measure->M_DynMod[k]->Module_Score=1.;
            
            for (int i=1; i<=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes; i++){
               double maxw_of_phen=0,y=0;
                y=Unconnected_module_phenotype_weight_scaled_DAVID(Uncoupled_ROW_Old_COL_MATCH, UA_lineup,i);

                if(talk) printf("\n\t Module %d phenotype %d has total weight %lg in the UNCONNECTED system:\n",k,i,y); //getchar();
                
                for(int b=1;b<=Coupled->Attractor_NR;b++){
                    // if the overall time spent in module k phenotype i decreased in coupled landscape, normalize by time spent in uncoupled.
                    // otherwise leave it as 1. This prevents stabilization of a weak module attractor from compensating for the destabilization of a good one.
                    // also, this keeps the measure between 0 and 1: ideal case =  phenotype is at least as "relevant" as alone overall.
                    
                    //if (DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]<y) {
                    //    DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]/=y;
                   // }
                   // else DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i]=1.;
                    maxw_of_phen+=DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i];
//                    if(maxw_of_phen<DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i])
  //                      maxw_of_phen=DynMod_Measure->M_DynMod[k]->Land_CH[b]->Phenotype_Stability[i];

                }
                // if the overall time spent in module k phenotype i decreased in coupled landscape, normalize by time spent in uncoupled.
                // otherwise leave it as 1. This prevents stabilization of a weak module attractor from compensating for the destabilization of a good one.
                // also, this keeps the measure between 0 and 1: ideal case =  phenotype is at least as "relevant" as alone overall.
                if(talk) printf("\t Module %d phenotype %d has total weight %lg in the COUPLED system:\n",k,i,maxw_of_phen); //getchar();
                
                if(maxw_of_phen<y) DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i]=maxw_of_phen/y;
                else DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i]=1.;
                
                if(talk) printf("\t\tPhenotype %d: Phenotype Stabilty[%d:%d] = %lg\n",i,k,i,DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i]);
                DynMod_Measure->M_DynMod[k]->Module_Score*=DynMod_Measure->M_DynMod[k]->Module_Phenotype_Score[i];
               // getchar();
            }
            DynMod_Measure->M_DynMod[k]->Module_Score=pow(DynMod_Measure->M_DynMod[k]->Module_Score,1./(double)DynMod_Measure->M_DynMod[k]->N_module_Phenotypes);
           
            DynMod_Measure->Switch_stability_measure*=DynMod_Measure->M_DynMod[k]->Module_Score;
           
            if(talk) printf("\n\t Module %d Switch Stability SS[%d] = %lg\n",k,k,DynMod_Measure->M_DynMod[k]->Module_Score);
            //getchar();
             
            
            delete_Overlap_matrix(New_ROW_Old_COL_MATCH);
            delete_Overlap_matrix(New_ROW_Old_ROW_FLyBY_MATCH);
            delete UA_lineup;
        }

    
        DynMod_Measure->Attractor_modularity_measure[0]=pow(DynMod_Measure->Attractor_modularity_measure[0],1./(double)Module_NR);
        if(talk) printf("\nFINAL Attractor Modularity Measure AMM = %lg\n",DynMod_Measure->Attractor_modularity_measure[0]);

       
        DynMod_Measure->Switch_stability_measure=pow(DynMod_Measure->Switch_stability_measure,1./(double)Module_NR);
        
        if(talk) printf("\nFINAL Switch Stability Measure SSM = %lg\n",DynMod_Measure->Switch_stability_measure);
        
        //LAW 3

        if(talk) printf("\n\n LAW 3\n\n");
        
        unsigned int Mod_comb_nr=1;
        for(int k=1;k<=Module_NR; k++)
            Mod_comb_nr*=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes;
        
        DynMod_Measure->Module_coordination_measure=((double)Mod_comb_nr - (double)Coupled->Attractor_NR)/(double)Mod_comb_nr;
        if(talk){
            printf("\t Number of module phenotype combinations, Mod_comb_nr=%d\n",Mod_comb_nr);
            printf("\t Coordination measure (how much this got reduced) MCM =%lg\n",DynMod_Measure->Module_coordination_measure);
        }
        
     //   printf("Mod_comb_nr=%d\t\t",Mod_comb_nr);
      //  printf(" DynMod_Measure->Module_coordination_measure=%lg\n",DynMod_Measure->Module_coordination_measure);
      //  getchar();
        
        if (N_trace>0) {
            char fn[300];
            FILE *f;
            sprintf(fn,"%s/%s/DynMod_RND/DM_%s__Ntr-%lu_Nr-%u_p-%.3lf.txt",
                    MODEL_Directory,Coupled->DIRname,BRN->name,N_trace,N_RND,p_error);
           // printf("%s\n", fn);
            
            double hprod=1;double h=0;
            double mcm_modules=1;
            DynMod_Measure->Module_Barrier_Quality=new double[Module_NR+2];
            
            for (int m=1; m<=Module_NR; m++) {
                Boolean_Dynamics *D;
                
                D=new Boolean_Dynamics(BRN_Module[m]);
                D->generate_full_state_graph_synchronously(0);
                h=D->normalized_distance_between_attractors();
                if(talk) printf("\t\tModule %d Mod_h_av=%lg  mcm =%lg \n",m,h,(pow(2,BRN_Module[m]->N) - (double)D->Attractor_NR)/(double)pow(2,BRN_Module[m]->N));//getchar();
                DynMod_Measure->Module_Barrier_Quality[m]=h;
                hprod=hprod*h;
                if (BRN_Module[m]->N>1) {
                    mcm_modules*=(pow(2,BRN_Module[m]->N) - (double)D->Attractor_NR)/(double)pow(2,BRN_Module[m]->N);
                }
                else {
                    if ((double)D->Attractor_NR==1) {
                        mcm_modules=0;
                    }
                }
 
               // delete D;D=NULL;
            }
            hprod=pow(hprod,1./(double)Module_NR);
            mcm_modules=pow(mcm_modules,1./(double)Module_NR);
            if(talk) printf("\t\t\tModule Geom Av: H=%lg\n\t\t\tModule coord Geom Av =%lg\n",hprod,mcm_modules);//getchar();
           
            h=Coupled->normalized_distance_between_attractors();
            if(talk) printf("\t\tH_coupled=%lg\n",h);
            hprod=hprod*h;
            if(talk) printf("2-level Attractor Difference Measure: H = Mod_h_av* H_coupled = %lg\n",hprod);
            
            DynMod_Measure->Module_Attractor_Distinction_measure=hprod;
           // getchar();
            
            
            // Module_Barrier_Quality:
            /*
             
            double beta_now=0.5 *log(1./p_error-1.);
            
            double B_prod=0.;
            for (int m=1; m<=Module_NR; m++) {
                Boolean_Dynamics *D;
                
                D=new Boolean_Dynamics(BRN_Module[m]);
                D->generate_full_state_graph_synchronously(0);
                if (D->Attractor_NR>1) {
                    D->calculate_noisy_landscape(beta_now);
                    double Bprod_mod=0.;
                    snode::slinklist::slink *sc2;
                    for (int a=1; a<=D->Attractor_NR; a++){
                        sc2=D->Attractor_G->node[a].Li.first_link_pt;
                        while(sc2!=NULL){
                            if (sc2->sneighbor!=a) {
                                if(sc2->link_props->lw>Bprod_mod) Bprod_mod=sc2->link_props->lw;
                            }
                            sc2=sc2->next(sc2);
                        }
                    }
                   // Bprod_mod=pow(Bprod_mod,1./(double)(D->Attractor_NR*(D->Attractor_NR-1)));
                    
                    D->clean_landscape_calc();
                    if(talk) printf("\t\t\tModule %d, BARRIER Geom Av: H=%lg\n",m,Bprod_mod);
                    
                    DynMod_Measure->Module_Barrier_Quality[m]=-log(Bprod_mod);
                    B_prod+=-log(Bprod_mod);
                    // delete D;D=NULL;
                }
                else {
                        DynMod_Measure->Module_Barrier_Quality[m]=0;
                }
            }
            //B_prod=pow(B_prod,1./(double)Module_NR);
            B_prod/=(double)Module_NR;
            
            // average of transition prob. across weakest barrier
            if(talk) printf("\t\t\tModule worst BARRIER  Av: H=%lg\n",B_prod);
            
            double Bc=0.;
            for (int a=1; a<=Coupled->Attractor_NR; a++)
                for (int b=1; b<=Coupled->Attractor_NR; b++)
                    if (a!=b) {
                        double x=Coupled->Trans_Basins[a][b]/(double)Coupled->Attractor_valleys[b]->W_Total;
                        if(x>Bc) Bc=x;
                    }
            //Bc=pow(Bc,1./(double)(Coupled->Attractor_NR*(Coupled->Attractor_NR-1)));
            if(talk) printf("\t\t\tCoupled worst BARRIER  : H=%lg\n",Bc);
            
            DynMod_Measure->Module_Barrier_Quality[0]=-log(Bc);
            
            //B_prod=sqrt((1.-B_prod)*(1.-Bc)); // geom mean of (1-words global)*(1- avmodule)
            
            //DynMod_Measure->Module_Barrier_Quality[0]=B_prod;
            //if(talk) printf("\t\t\tModule_Barrier_Quality  : %lg\n",DynMod_Measure->Module_Barrier_Quality);
            //getchar();
            */
            
            
            f=fopen(fn,"w");
            fprintf(f,"%lg\t%lg\t%lg\t\t%lg\t%lg\t%lg\t\t",DynMod_Measure->Attractor_modularity_measure[0],
                                               DynMod_Measure->Switch_stability_measure,
                                               DynMod_Measure->Module_Attractor_Distinction_measure*DynMod_Measure->Module_coordination_measure,
                                               DynMod_Measure->Module_Attractor_Distinction_measure,
                                               DynMod_Measure->Module_coordination_measure,
                                               mcm_modules);
            
           // for(int k=1;k<=Module_NR; k++) fprintf(f,"%lg\t",DynMod_Measure->Module_Barrier_Quality[k]);
            
            int mpprod=1;
            for(int k=1;k<=Module_NR; k++) {
               fprintf(f,"%d\t",DynMod_Measure->M_DynMod[k]->N_module_Phenotypes);
                mpprod*=DynMod_Measure->M_DynMod[k]->N_module_Phenotypes;
            }
            fprintf(f,"%d\t",mpprod);
            
            fprintf(f,"%d\t",rnd_version);
           // fclose(f);
            
            if(talk) printf("FINAL Module coordination & and attractor distinction measure (MCM) = %lg \n",DynMod_Measure->Module_coordination_measure*DynMod_Measure->Module_Attractor_Distinction_measure);
            return(f);
        }
        else return(NULL);
        
        if(talk) {
            printf("\n\n\n SUMMARY:\n");
            printf("Attractor modularity measure (AMM) = %lg \n",DynMod_Measure->Attractor_modularity_measure[0]);
            //getchar();
            
            printf("Switch stability measure (SSM) = %lg \n",DynMod_Measure->Switch_stability_measure);
            printf("Module coordination & and attractor distinction measure (MCM) = %lg \n",DynMod_Measure->Module_coordination_measure*DynMod_Measure->Module_Attractor_Distinction_measure);
            printf("\tModule coordination measure (MCM) = %lg \n",DynMod_Measure->Module_coordination_measure);
            printf("\t2-level Attractor distinction measure (MQM) = %lg \n",DynMod_Measure->Module_Attractor_Distinction_measure);
            //printf(" Barrier quality Measure (BQM) = %lg \n",DynMod_Measure->Module_Barrier_Quality[0]);
        }
      //  getchar();
    }
    
    double Unconnected_module_phenotype_weight_scaled(int Mod_ID,int phen_ID){
        double w=0,wt=0;
        for(int b=1;b<=Unconnected->Attractor_NR;b++){
            wt+=Unconnected->Attractor_valleys[b]->W_Total;
            if(Module_Attr_in_Uncoupled(b,Mod_ID)==phen_ID)
                w+=Unconnected->Attractor_valleys[b]->W_Total;
        }
        return(w/(double)wt);
    }
    
    double *Module_Attr_spread_in_Coupled(int C_AID,int Mod_ID,int n_phen,double w_share[]){
        double *w;
        
       // printf("STARTING\n");getchar();
        w=new double[n_phen+1];
        for (int i=1; i<=n_phen; i++) w[i]=0;
        
        if (Coupled->Attractor_valleys[C_AID]->N_a>1) {
            for(int i=1;i<=Coupled->Attractor_valleys[C_AID]->N_a;i++){
          //      printf("Coupled state %d gets share %lg of coupled valley %d\n",i,w_share[i],C_AID);
                w[Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[C_AID][i],Mod_ID)]+=
                        w_share[i];
          //      printf("\t For module %d - share goes to phenotype %d \n",Mod_ID,Module_Attr_in_Uncoupled(Unconnected_attractor_match_of_Composite_attractor[C_AID][i],Mod_ID));
               // getchar();
            }
        }
        else {
           // w_share=new double[2];
           // w_share[0]=Coupled->Attractor_valleys[C_AID]->W_Total;
            //w_share[1]=1.;
            w[Module_Attr_in_Uncoupled(CHOSEN_Unconnected_attractor_match_of_Composite_attractor[C_AID],Mod_ID)]+=1;
        //    printf("Coupled state 1 is the entire coupled valley %d\n",C_AID);
        //    printf("\t For module %d - share goes to phenotype %d \n",Mod_ID,Module_Attr_in_Uncoupled(CHOSEN_Unconnected_attractor_match_of_Composite_attractor[C_AID],Mod_ID));
           //getchar();
        }
        return(w);
    }
    
    
    
    void match_old_and_new_attractors(){
        unsigned int h,h_min,n_mapped_onto_a,a_match,onem = 0,stone = 0;
        unsigned long long *order_b;
        double *order_h;
        
        Unconnected_attractor_match_of_Composite_attractor=new intpointer[Coupled->Attractor_NR+1];
        Rank_in_matching_Unconnected_attractor_for_Composite_attractor=new intpointer[Coupled->Attractor_NR+1];
        Unconnected_attractor_minHamming_of_Composite_attractor=new doublepointer[Coupled->Attractor_NR+1];
        CHOSEN_Unconnected_attractor_match_of_Composite_attractor=new int[Coupled->Attractor_NR+1];
        
        for(int b=1;b<=Coupled->Attractor_NR;b++){
            Unconnected_attractor_match_of_Composite_attractor[b]=new int[Coupled->Attractor_valleys[b]->N_a+1];
            Rank_in_matching_Unconnected_attractor_for_Composite_attractor[b]=new int[Coupled->Attractor_valleys[b]->N_a+1];
            Unconnected_attractor_minHamming_of_Composite_attractor[b]=new double[Coupled->Attractor_valleys[b]->N_a+1];
            for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
                Unconnected_attractor_match_of_Composite_attractor[b][i]=0;
                Rank_in_matching_Unconnected_attractor_for_Composite_attractor[b][i]=0;
                Unconnected_attractor_minHamming_of_Composite_attractor[b][i]=0;
            }
            Unconnected_attractor_minHamming_of_Composite_attractor[b][0]=BRN->N+1;
            Unconnected_attractor_match_of_Composite_attractor[b][0]=(int)Coupled->Attractor_valleys[b]->N_a;
            CHOSEN_Unconnected_attractor_match_of_Composite_attractor[b]=0;
        }
        
        
        for(int b=1;b<=Coupled->Attractor_NR;b++){
            for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++){
                h_min=BRN->N;a_match=0;
                for(int a=1;a<=Unconnected->Attractor_NR;a++){
                    for(int j=1;j<=Unconnected->Attractor_valleys[a]->N_a;j++){
                        h=get_Hamming_distance_of_states(Coupled->Attractor_valleys[b]->Basin_Stored[Coupled->Attractor_valleys[b]->Attractor[i]]->s,Unconnected->Attractor_valleys[a]->Basin_Stored[Unconnected->Attractor_valleys[a]->Attractor[j]]->s);
                        if(h<h_min) {h_min=h;a_match=a;}
                    }
                }
                Unconnected_attractor_match_of_Composite_attractor[b][i]=a_match;
                Unconnected_attractor_minHamming_of_Composite_attractor[b][i]=h_min;
                if(h_min<Unconnected_attractor_minHamming_of_Composite_attractor[b][0]){
                    CHOSEN_Unconnected_attractor_match_of_Composite_attractor[b]=a_match;
                    Unconnected_attractor_minHamming_of_Composite_attractor[b][0]=h_min;
                }
            }
        }
    
        for(int a=1;a<=Unconnected->Attractor_NR;a++){
            n_mapped_onto_a=0;
            for(int b=1;b<=Coupled->Attractor_NR;b++)
                 for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++)
                     if(Unconnected_attractor_match_of_Composite_attractor[b][i]==a) {n_mapped_onto_a++;onem=b;stone=i;}
            if(n_mapped_onto_a==1) Rank_in_matching_Unconnected_attractor_for_Composite_attractor[onem][stone]=1;
            else{
                order_h=new double[n_mapped_onto_a+1];
                order_b=new unsigned long long[n_mapped_onto_a+1];
                n_mapped_onto_a=0;
                for(int b=1;b<=Coupled->Attractor_NR;b++)
                    for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++)
                        if(Unconnected_attractor_match_of_Composite_attractor[b][i]==a) {
                            n_mapped_onto_a++;
                            order_h[n_mapped_onto_a]=Unconnected_attractor_minHamming_of_Composite_attractor[b][i];
                            order_b[n_mapped_onto_a]=n_mapped_onto_a;
                    }
                hpsort_coshuffle((unsigned long long)n_mapped_onto_a, order_h, order_b);
                
                int nn=0;
                for(int b=1;b<=Coupled->Attractor_NR;b++)
                    for(int i=1;i<=Coupled->Attractor_valleys[b]->N_a;i++)
                        if(Unconnected_attractor_match_of_Composite_attractor[b][i]==a) {
                            nn++;
                            for(int k=1;k<=n_mapped_onto_a;k++)
                                if(order_b[k]==nn)
                                    Rank_in_matching_Unconnected_attractor_for_Composite_attractor[b][i]=k;
                        }
                delete[] order_b;order_b=NULL;
                delete[] order_h;order_h=NULL;
            }
        }
    }
    
    char *get_module_substate(const char sfull[],int Module_k){
        int ninmodule=0;
        char *smod;
        for(int i=1;i<=BRN->N;i++){
            if (BRN->MDAT->Module_ID[i]==Module_k) {
                ninmodule++;
            }
        }
        smod=new char[ninmodule+1];
        ninmodule=0;
        for(int i=1;i<=BRN->N;i++){
            if (BRN->MDAT->Module_ID[i]==Module_k) {
                smod[ninmodule]=sfull[i-1];
                ninmodule++;
            }
        }
        
        smod[ninmodule]=0;
        return(smod);
    }
    
    void SET_module_substate_Coupled(char smod[],int Module_k){
        int ninmodule=0;
        
        for(int i=1;i<=BRN->N;i++){
            if (BRN->MDAT->Module_ID[i]==Module_k) {
                Coupled->s[i-1]=smod[ninmodule];
                ninmodule++;
            }
        }
        Coupled->s[BRN->N]=0;
    }
    
    void SET_module_substate_Unconnected(char smod[],int Module_k){
        int ninmodule=0;
        
        for(int i=1;i<=BRN->N;i++){
            if (BRN->MDAT->Module_ID[i]==Module_k) {
                Unconnected->s[i-1]=smod[ninmodule];
                ninmodule++;
            }
        }
        Coupled->s[BRN->N]=0;
    }
    
    void print_module_substate_Unconnected(int Module_k){
        
        for(int i=1;i<=BRN->N;i++){
            if (BRN->MDAT->Module_ID[i]==Module_k) {
                printf("%d",Unconnected->s[i-1]-48);
            }
        }
    }

    int get_NR_of_module_assigned_nodes(){
        int ninmodule=0;
        for(int i=1;i<=BRN->N;i++){
            if (BRN->MDAT->Module_ID[i]>0) {
                ninmodule++;
            }
        }
        return(ninmodule);
    }
    
    int get_NR_of_module_assigned_nodes(Boolean_RegNet *A){
        int ninmodule=0;
        if (A->MDAT->Module_ID==NULL) return(0);
        for(int i=1;i<=A->N;i++){
            if (A->MDAT->Module_ID[i]>0) {
                ninmodule++;
            }
        }
        return(ninmodule);
    }

    void group_Uncoupled_Attractors_by_Module_Combinations(){
        intpointer *colored_Uncoupled;
        State_index_MAP_type::iterator It_tr;
        Module_Attactor_ID_Combination_MAP_type::iterator It_tr_M;
        char *snow,*s_M;
        unsigned long int N_comb;
        longint_ARRAY *gl;
  
        colored_Uncoupled=new intpointer[Module_NR+1];
        N_comb=1;
        for(int k=0;k<=Module_NR;k++){
            colored_Uncoupled[k]=new int[Unconnected->Attractor_NR+1];
            for(int i=0;i<=Unconnected->Attractor_NR;i++){
                colored_Uncoupled[k][i]=0;
            }
        }
        
        State_index_MAP_type tracedown;
        
        for(int k=1;k<=Module_NR;k++){
            int m_a=0;
            for(int i=1;i<=Unconnected->Attractor_NR;i++){
                int newa=1;
               // printf("For unconnected attr %d\n with %ld states:\n",i,Unconnected->Attractor_valleys[i]->N_a);
                for(int j=1;j<=Unconnected->Attractor_valleys[i]->N_a;j++){
                  //  printf("\t state %d: %s  \n\t",j,Unconnected->Attractor_valleys[i]->Basin_Stored[Unconnected->Attractor_valleys[i]->Attractor[j]]->s);
                    snow=get_module_substate(Unconnected->Attractor_valleys[i]->Basin_Stored[Unconnected->Attractor_valleys[i]->Attractor[j]]->s,k);
                  //  printf("module %d substate is %s\n",k,snow);
                    if(tracedown.find(snow) != tracedown.end()) {
                        It_tr=tracedown.find(snow);
                        colored_Uncoupled[k][i]=(int)It_tr->second;
                       // printf("\t\t\tmodule %d substate is %s - module attractor colored_Uncoupled[%d][%d]= %d\n",k,snow,k,i,colored_Uncoupled[k][i]);
                       // printf("\n\n\na: colored_Uncoupled[1][5] = %d\n",colored_Uncoupled[1][5]);
                       // getchar();
                    }
                    else {
                        if(newa==1){
                            newa=0;
                            m_a++;
                           // printf("\n\n\nb1: colored_Uncoupled[1][5] = %d\n",colored_Uncoupled[1][5]);
                           //  printf("\t\t\t About to set: module %d substate is %s - module attractor colored_Uncoupled[%d][%d]= %d\n",k,snow,k,i,m_a);
                          //  getchar();
                            colored_Uncoupled[k][i]=m_a;
                          //   printf("\t\t\tmodule %d substate is %s - module attractor colored_Uncoupled[%d][%d]= %d\n",k,snow,k,i,colored_Uncoupled[k][i]);
                            
                           // getchar();
                        }
                        tracedown[snow]=m_a;
                    }
                   // printf("\n\n\nc: colored_Uncoupled[1][5] = %d\n",colored_Uncoupled[1][5]);
                   // getchar();
                    delete[] snow;snow=NULL;
                }
            }
            N_comb*=m_a;
            tracedown.clear();
          //  printf("\n\n\ne: colored_Uncoupled[1][5] = %d\n",colored_Uncoupled[1][5]);
          //  getchar();
        }
        
       
       // for(int i=1;i<=Unconnected->Attractor_NR;i++)
       //      for(int k=1;k<=Module_NR;k++)
       //          printf("colored_Uncoupled[%d][%d] = %d\n",k,i,colored_Uncoupled[k][i]);
      //  getchar();
        
        s_M=new char[Module_NR+1];
        for(int i=1;i<=Unconnected->Attractor_NR;i++){
            for(int k=1;k<=Module_NR;k++) {
                s_M[k-1]=colored_Uncoupled[k][i]+48;
               // printf("\t\t s_M[%d-1]=%d\n",k,colored_Uncoupled[k][i]);
            }
            s_M[Module_NR]=0;
           // printf("Unconn attr %d\t s_M = %s\n",i,s_M);
            if(DynMod.find(s_M) != DynMod.end()) {
                It_tr_M=DynMod.find(s_M);
                It_tr_M->second.add_element(i);
               //  printf("\t adding %d\n",i);
            }
            else {
                gl=new longint_ARRAY();
                DynMod[s_M]=*gl;
                It_tr_M=DynMod.find(s_M);
                It_tr_M->second.add_element(i);
               // printf("\t adding %d\n",i);
            }
         //   getchar();
        }
        delete[] s_M;s_M=NULL;
        
        for(int k=0;k<=Module_NR;k++){ delete[] colored_Uncoupled[k];colored_Uncoupled[k]=NULL;}
        delete[] colored_Uncoupled;colored_Uncoupled=NULL;
    }
    
    void printf_DynMod(){
        Module_Attactor_ID_Combination_MAP_type::const_iterator It_tr;
        
        for (It_tr = DynMod.begin(); It_tr != DynMod.end(); ++It_tr) {
            printf("\t Module Combination %s covers %ld Uncopupled Attractors:\n",
                   It_tr->first.c_str(),It_tr->second.N);
            for(int l=1;l<=It_tr->second.N;l++) printf("    %ld",It_tr->second.A[l]);
            printf("\n");
        }
    }
   
    int nr_mod_phen_from_DynMod(int Mod_ID){
        Module_Attactor_ID_Combination_MAP_type::const_iterator It_tr;
        
        int n=0;
        for (It_tr = DynMod.begin(); It_tr != DynMod.end(); ++It_tr) {
            if((int)(It_tr->first.c_str()[Mod_ID-1]-48)>n) n=(int)(It_tr->first.c_str()[Mod_ID-1]-48);
            //printf("\t Module Combination %s: mod %d id is c:%c - d:%d\t\t n=%d\n",It_tr->first.c_str(),Mod_ID,It_tr->first.c_str()[Mod_ID-1],It_tr->first.c_str()[Mod_ID-1]-48,n);
        }
        return(n);
    }
    
    int Module_Attr_in_Uncoupled(int UC_A,int Mod_ID){
        Module_Attactor_ID_Combination_MAP_type::iterator It_tr;
        
        for (It_tr = DynMod.begin(); It_tr != DynMod.end(); ++It_tr) {
            if(It_tr->second.check_for_element(UC_A))
                return(It_tr->first.c_str()[Mod_ID-1]-48);
        }
        printf("Error!\n Uncoupled attractor %d not found in any DynMod combination list!\n",UC_A);getchar();
        return(0);
    }
    
    int N_DynMod(){
        Module_Attactor_ID_Combination_MAP_type::const_iterator It_tr;
        int n=0;
        for (It_tr = DynMod.begin(); It_tr != DynMod.end(); ++It_tr) {
            n++;
        }
        return(n);
    }
    
    unsigned int Valley_index(unsigned long int b){
        Module_Attactor_ID_Combination_MAP_type::const_iterator It_tr;
        int n;
        n=0;
        for (It_tr = DynMod.begin(); It_tr != DynMod.end(); ++It_tr) {
            n++;
            for(int l=1;l<=It_tr->second.N;l++)
                if(It_tr->second.A[l]==b) return(n);
        }
        printf("ERROR: Unconnected attractor %ld is not found on ANY combination list!!!\n",b);
        exit(1);
        return(0);
    }
    
    
    longint_ARRAY *get_non_module_node_list(){
        longint_ARRAY *gl;
        gl=new longint_ARRAY();
        for(int i=1;i<=BRN->N;i++)
            if(Module_ID[i]==0) gl->add_element(i);
        return(gl);
    }
    
    void export_paired_landscape(unsigned long int N_trace,unsigned int N_RND,double p_error,int V_MIN,double minflux_fraction_of_p){
        unsigned long int i,b,a,a1,a2;
        unsigned long long w2;
        // State_index_MAP_type::iterator It_tr;
        FILE *f,*f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8;
        FILEpointer *f_Extras;
        char fn[300];//,*s1,*s2;
        //int ok1,ok2;
        int ST_order1,ST_order2;
        double lf1,lf2;
    
        sprintf(fn,"mkdir %s/%s/Cytoskape/Landscape_PAIRED\n", MODEL_Directory_system,BRN->name);
        
        system(fn);
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf.txt\n",
                     MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error);
        f=fopen(fn,"w");
        fprintf(f,"State_st\tTo\tState_end\tPaired_Weight\tST_order_Uncoupled\tLinkFlux_Uncoupled\tST_order_Coupled\tLinkFlux_Coupled\n");
        
        sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf__Uncoupled_BasinID.noa\n",
                MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error);
        f1=fopen(fn,"w");
        sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf__Coupled_BasinID_Match.noa\n",
                MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error);
        f2=fopen(fn,"w");
        sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf__Coupled_BasinID.noa\n",
                MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error);
        f5=fopen(fn,"w");
       
        sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf__Uncoupled_Energy.noa\n",
                MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error);
        f3=fopen(fn,"w");
        sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf__Coupled_Energy.noa\n",
                MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error);
        f4=fopen(fn,"w");
        sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf__ATTR_Uncoupled.noa\n",
                MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error);
        f6=fopen(fn,"w");
        sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf__ATTR_Coupled.noa\n",
                MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error);
        f7=fopen(fn,"w");
        sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf__GenesON.noa\n",
                MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error);
        f8=fopen(fn,"w");
        
        fprintf(f1,"Uncoupled_Basin_ID\n");
        //fprintf(f2,"Coupled_Basin_ID_Matched_rank\n");
        fprintf(f5,"Coupled_Basin_ID\n");
        fprintf(f3,"Uncoupled_Energy\n");
        fprintf(f4,"Coupled_Energy\n");
        fprintf(f6,"ATTR_Uncoupled\n");
        fprintf(f7,"ATTR_Coupled\n");
        fprintf(f8,"GenesON\n");
        
        longint_ARRAY *gl;
        gl=get_non_module_node_list();
        f_Extras=new FILEpointer[gl->N+1];
        for(int i=1;i<=gl->N;i++){
            sprintf(fn,"%s/%s/Cytoskape/Landscape_PAIRED/%s_State_Transition_Graph__Ntr-%lu_Nr-%u_p-%.3lf__%s.noa\n",
                    MODEL_Directory,BRN->name,BRN->name,N_trace,N_RND,p_error,BRN->MDAT->Node_name[gl->A[i]]);
            f_Extras[i]=fopen(fn,"w");
            fprintf(f_Extras[i],"%s\n",BRN->MDAT->Node_name[gl->A[i]]);
        }
        
        for(b=1;b<=Unconnected->Attractor_NR;b++){
            for(i=1;i<=Unconnected->Attractor_valleys[b]->N_ms;i++){
                if(  (Unconnected->Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)
                   ||(Unconnected->Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                    
                    
                    Unconnected->set_state(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                    Unconnected->synchronously_update_all_Gates();
                    // taking care of all synchronous links out of a heavy enough Uncoupled state
                    
                    lf1=Unconnected->get_link_flux(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->s,&ST_order1,p_error);
                    lf2=Coupled->get_link_flux(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->s,&ST_order2,p_error);
                    lf1=lf1/(double)((N_trace+1)*N_RND);lf2=lf2/(double)((N_trace+1)*N_RND);
                    
                    fprintf(f,"%s\t->\t%s\t%lg\t%d\t%lg\t%d\t%lg\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->s,(lf1+lf2)/2.,ST_order1,lf1,ST_order2,lf2);
                    
                    
                    a=Coupled->check_if_state_has_Basin_ID(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                    if(a==0) a=Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                    fprintf(f1,"%s = %ld\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,b);
                    fprintf(f5,"%s = %ld\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,a);
                    //fprintf(f2,"%s = %d_%d\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,
                    //        Unconnected_attractor_match_of_Composite_attractor[a],
                    //        Rank_in_matching_Unconnected_attractor_for_Composite_attractor[a]);
                    fprintf(f3,"%s = %lg\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->get_energy_of_state(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,(double)((N_trace+1)*N_RND),p_error));
                    fprintf(f4,"%s = %lg\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->get_energy_of_state(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,(double)((N_trace+1)*N_RND),p_error));
                    fprintf(f6,"%s = %d\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->check_if_state_is_an_Attractor_State(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s));
                    fprintf(f7,"%s = %d\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->check_if_state_is_an_Attractor_State(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s));
                    
                    fprintf(f8,"%s = ",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                    for(int k=1;k<=BRN->N;k++) if(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s[k-1]!=48) fprintf(f8,"%s ",BRN->MDAT->Node_name[k]);
                    fprintf(f8,"\n");
                    for(int k=1;k<=gl->N;k++){
                        if(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s[gl->A[k]-1]!=48)
                             fprintf(f_Extras[k],"%s = 1\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                        else fprintf(f_Extras[k],"%s = 0\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                    }
                    
                    a=Coupled->check_if_state_has_Basin_ID(Unconnected->s);
                    if(a==0) a=Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Unconnected->s);
                    fprintf(f1,"%s = %ld\n",Unconnected->s,b);
                    fprintf(f5,"%s = %ld\n",Unconnected->s,a);
                    //fprintf(f2,"%s = %d_%d\n",Unconnected->s,
                    //        Unconnected_attractor_match_of_Composite_attractor[a],
                   //         Rank_in_matching_Unconnected_attractor_for_Composite_attractor[a]);
                    fprintf(f3,"%s = %lg\n",Unconnected->s,Unconnected->get_energy_of_state(Unconnected->s,(double)((N_trace+1)*N_RND),p_error));
                    fprintf(f4,"%s = %lg\n",Unconnected->s,Coupled->get_energy_of_state(Unconnected->s,(double)((N_trace+1)*N_RND),p_error));
                    fprintf(f6,"%s = %d\n",Unconnected->s,Unconnected->check_if_state_is_an_Attractor_State(Unconnected->s));
                    fprintf(f7,"%s = %d\n",Unconnected->s,Coupled->check_if_state_is_an_Attractor_State(Unconnected->s));
                    fprintf(f8,"%s = ",Unconnected->s);
                    for(int k=1;k<=BRN->N;k++) if(Unconnected->s[k-1]!=48) fprintf(f8,"%s ",BRN->MDAT->Node_name[k]);
                    fprintf(f8,"\n");
                    for(int k=1;k<=gl->N;k++){
                        if(Unconnected->s[gl->A[k]-1]!=48)
                            fprintf(f_Extras[k],"%s = 1\n",Unconnected->s);
                        else fprintf(f_Extras[k],"%s = 0\n",Unconnected->s);
                    }
                    // taking care of all high-enough flux OTHER links out of this state
                    
                    for(int ipoz=1;ipoz<=BRN->N;ipoz++){ // iterates through the positions in the end of sync link state, flips one and uses the new state as alt end to weak link. Then flips position back.
                        if(Unconnected->s[ipoz]=='1') Unconnected->s[ipoz]='0'; else Unconnected->s[ipoz]='1'; // flip the ith node of the next state
                       
                        a1=Unconnected->check_if_state_was_sampled_enough(Unconnected->s,V_MIN);
                        a2=Coupled->check_if_state_was_sampled_enough(Unconnected->s,V_MIN);
                        if((a1>0)||(a2>0)){
                            lf1=Unconnected->get_link_flux(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->s,&ST_order1,p_error);
                            lf2=Coupled->get_link_flux(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->s,&ST_order2,p_error);
                            lf1=lf1/(double)((N_trace+1)*N_RND);lf2=lf2/(double)((N_trace+1)*N_RND);
                            if((lf1>=p_error*minflux_fraction_of_p)||(lf2>=p_error*minflux_fraction_of_p)){
                                fprintf(f,"%s\t->\t%s\t%lg\t%d\t%lg\t%d\t%lg\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->s,(lf1+lf2)/2.,ST_order1,lf1,ST_order2,lf2);
                            
                                if(a1==0) a1=Unconnected->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Unconnected->s);
                                if(a2==0) a2=Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Unconnected->s);
                                fprintf(f1,"%s = %ld\n",Unconnected->s,a1);
                                fprintf(f5,"%s = %ld\n",Unconnected->s,a2);
                                //fprintf(f2,"%s = %d_%d\n",Unconnected->s,
                                //        Unconnected_attractor_match_of_Composite_attractor[a2],
                                //        Rank_in_matching_Unconnected_attractor_for_Composite_attractor[a2]);
                                fprintf(f3,"%s = %lg\n",Unconnected->s,Unconnected->get_energy_of_state(Unconnected->s,(double)((N_trace+1)*N_RND),p_error));
                                fprintf(f4,"%s = %lg\n",Unconnected->s,Coupled->get_energy_of_state(Unconnected->s,(double)((N_trace+1)*N_RND),p_error));
                                fprintf(f6,"%s = %d\n",Unconnected->s,Unconnected->check_if_state_is_an_Attractor_State(Unconnected->s));
                                fprintf(f7,"%s = %d\n",Unconnected->s,Coupled->check_if_state_is_an_Attractor_State(Unconnected->s));
                                fprintf(f8,"%s = ",Unconnected->s);
                                for(int k=1;k<=BRN->N;k++) if(Unconnected->s[k-1]!=48) fprintf(f8,"%s ",BRN->MDAT->Node_name[k]);
                                fprintf(f8,"\n");
                                for(int k=1;k<=gl->N;k++){
                                    if(Unconnected->s[gl->A[k]-1]!=48)
                                        fprintf(f_Extras[k],"%s = 1\n",Unconnected->s);
                                    else fprintf(f_Extras[k],"%s = 0\n",Unconnected->s);
                                }
                            }
                        }
                        if(Unconnected->s[ipoz]=='1') Unconnected->s[ipoz]='0'; else Unconnected->s[ipoz]='1'; // flip the ith node of the next state back
                    }
                }
                else {
                    a=Coupled->check_if_state_has_Basin_ID(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                    if(a>0) {
                        w2=Coupled->Attractor_valleys[a]->get_visits_to_state__Basin(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                        if(w2>=V_MIN){
                            Coupled->set_state(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                            Coupled->synchronously_update_all_Gates();
                            // taking care of all synchronous links out of a heavy enough Coupled state (weakly sampled in Unconnected)
                    
                            lf1=Unconnected->get_link_flux(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,&ST_order1,p_error);
                            lf2=Coupled->get_link_flux(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,&ST_order2,p_error);
                            lf1=lf1/(double)((N_trace+1)*N_RND);lf2=lf2/(double)((N_trace+1)*N_RND);
                            fprintf(f,"%s\t->\t%s\t%lg\t%d\t%lg\t%d\t%lg\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,(lf1+lf2)/2.,ST_order1,lf1,ST_order2,lf2);
                            
                            fprintf(f1,"%s = %ld\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,b);
                            fprintf(f5,"%s = %ld\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,a);
                           // fprintf(f2,"%s = %d_%d\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,
                            //        Unconnected_attractor_match_of_Composite_attractor[a],
                            //        Rank_in_matching_Unconnected_attractor_for_Composite_attractor[a]);
                            fprintf(f3,"%s = %lg\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->get_energy_of_state(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,(double)((N_trace+1)*N_RND),p_error));
                            fprintf(f4,"%s = %lg\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->get_energy_of_state(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,(double)((N_trace+1)*N_RND),p_error));
                            fprintf(f6,"%s = %d\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->check_if_state_is_an_Attractor_State(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s));
                            fprintf(f7,"%s = %d\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->check_if_state_is_an_Attractor_State(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s));
                            fprintf(f8,"%s = ",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                            for(int k=1;k<=BRN->N;k++) if(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s[k-1]!=48) fprintf(f8,"%s ",BRN->MDAT->Node_name[k]);
                            fprintf(f8,"\n");
                            for(int k=1;k<=gl->N;k++){
                                if(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s[gl->A[k]-1]!=48)
                                    fprintf(f_Extras[k],"%s = 1\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                                else fprintf(f_Extras[k],"%s = 0\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s);
                            }
                            
                            unsigned long int b2;
                            b2=Unconnected->check_if_state_has_Basin_ID(Coupled->s);
                            if(b2==0) b2=Unconnected->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Coupled->s);
                            a1=Coupled->check_if_state_has_Basin_ID(Coupled->s);
                            if(a1==0) a1=Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Coupled->s);
                            fprintf(f1,"%s = %ld\n",Coupled->s,b2);
                            fprintf(f5,"%s = %ld\n",Coupled->s,a1);
                           // fprintf(f2,"%s = %d_%d\n",Coupled->s,
                            //        Unconnected_attractor_match_of_Composite_attractor[a1],
                            //        Rank_in_matching_Unconnected_attractor_for_Composite_attractor[a1]);
                            fprintf(f3,"%s = %lg\n",Coupled->s,Unconnected->get_energy_of_state(Coupled->s,(double)((N_trace+1)*N_RND),p_error));
                            fprintf(f4,"%s = %lg\n",Coupled->s,Coupled->get_energy_of_state(Coupled->s,(double)((N_trace+1)*N_RND),p_error));
                            fprintf(f6,"%s = %d\n",Coupled->s,Unconnected->check_if_state_is_an_Attractor_State(Coupled->s));
                            fprintf(f7,"%s = %d\n",Coupled->s,Coupled->check_if_state_is_an_Attractor_State(Coupled->s));
                            fprintf(f8,"%s = ",Coupled->s);
                            for(int k=1;k<=BRN->N;k++) if(Coupled->s[k-1]!=48) fprintf(f8,"%s ",BRN->MDAT->Node_name[k]);
                            fprintf(f8,"\n");
                            for(int k=1;k<=gl->N;k++){
                                if(Coupled->s[gl->A[k]-1]!=48)
                                    fprintf(f_Extras[k],"%s = 1\n",Coupled->s);
                                else fprintf(f_Extras[k],"%s = 0\n",Coupled->s);
                            }
                            
                            // taking care of all high-enough flux OTHER links out of this state
                            
                            for(int ipoz=1;ipoz<=BRN->N;ipoz++){ // iterates through the positions in the end of sync link state, flips one and uses the new state as alt end to weak link. Then flips position back.
                                if(Coupled->s[ipoz]=='1') Coupled->s[ipoz]='0'; else Coupled->s[ipoz]='1'; // flip the ith node of the next state
                                a1=Unconnected->check_if_state_was_sampled_enough(Coupled->s,V_MIN);
                                a2=Coupled->check_if_state_was_sampled_enough(Coupled->s,V_MIN);
                                if((a1>0)||(a2>0)){
                                    lf1=Unconnected->get_link_flux(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,&ST_order1,p_error);
                                    lf2=Coupled->get_link_flux(Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,&ST_order2,p_error);
                                    lf1=lf1/(double)((N_trace+1)*N_RND);lf2=lf2/(double)((N_trace+1)*N_RND);
                                    if((lf1>=p_error*minflux_fraction_of_p)||(lf2>=p_error*minflux_fraction_of_p)){
                                        fprintf(f,"%s\t->\t%s\t%lg\t%d\t%lg\t%d\t%lg\n",Unconnected->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,(lf1+lf2)/2.,ST_order1,lf1,ST_order2,lf2);
                                        
                                        if(a1==0) a1=Unconnected->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Coupled->s);
                                        if(a2==0) a2=Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Coupled->s);
                                        fprintf(f1,"%s = %ld\n",Coupled->s,a1);
                                        fprintf(f5,"%s = %ld\n",Coupled->s,a2);
                                       // fprintf(f2,"%s = %d_%d\n",Coupled->s,
                                       //         Unconnected_attractor_match_of_Composite_attractor[a2],
                                       //         Rank_in_matching_Unconnected_attractor_for_Composite_attractor[a2]);
                                        fprintf(f3,"%s = %lg\n",Coupled->s,Unconnected->get_energy_of_state(Coupled->s,(double)((N_trace+1)*N_RND),p_error));
                                        fprintf(f4,"%s = %lg\n",Coupled->s,Coupled->get_energy_of_state(Coupled->s,(double)((N_trace+1)*N_RND),p_error));
                                        fprintf(f6,"%s = %d\n",Coupled->s,Unconnected->check_if_state_is_an_Attractor_State(Coupled->s));
                                        fprintf(f7,"%s = %d\n",Coupled->s,Coupled->check_if_state_is_an_Attractor_State(Coupled->s));
                                        fprintf(f8,"%s = ",Coupled->s);
                                        for(int k=1;k<=BRN->N;k++) if(Coupled->s[k-1]!=48) fprintf(f8,"%s ",BRN->MDAT->Node_name[k]);
                                        fprintf(f8,"\n");
                                        for(int k=1;k<=gl->N;k++){
                                            if(Coupled->s[gl->A[k]-1]!=48)
                                                fprintf(f_Extras[k],"%s = 1\n",Coupled->s);
                                            else fprintf(f_Extras[k],"%s = 0\n",Coupled->s);
                                        }
                                    }
                                }
                                if(Coupled->s[ipoz]=='1') Coupled->s[ipoz]='0'; else Coupled->s[ipoz]='1'; // flip the ith node of the next state back
                            }
                        }
                    }
                }
            }
        }
    
        for(b=1;b<=Coupled->Attractor_NR;b++){
            for(i=1;i<=Coupled->Attractor_valleys[b]->N_ms;i++){
                if(  (Coupled->Attractor_valleys[b]->Basin_Stored[i]->visits>=V_MIN)
                   ||(Coupled->Attractor_valleys[b]->check_if_state_is_an_Attractor_State(i)>0)){
                    
                    a=Unconnected->check_if_state_has_Basin_ID(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s);
                    if(a==0){
                        Coupled->set_state(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s);
                        Coupled->synchronously_update_all_Gates();
                        // taking care of all synchronous links out of a heavy enough Coupled state, assuming they were not sampled by Uncoupled
                        
                        lf1=Unconnected->get_link_flux(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,&ST_order1,p_error);
                        lf2=Coupled->get_link_flux(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,&ST_order2,p_error);
                        
                        lf1=lf1/(double)((N_trace+1)*N_RND);lf2=lf2/(double)((N_trace+1)*N_RND);
                        fprintf(f,"%s\t->\t%s\t%lg\t%d\t%lg\t%d\t%lg\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,(lf1+lf2)/2.,ST_order1,lf1,ST_order2,lf2);
                        
                        a=Unconnected->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s);
                        if(a==0) a=Unconnected->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s);
                        
                        fprintf(f1,"%s = %ld\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,a);
                        fprintf(f5,"%s = %ld\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,b);
                       // fprintf(f2,"%s = %d_%d\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,
                        //        Unconnected_attractor_match_of_Composite_attractor[b],
                        //        Rank_in_matching_Unconnected_attractor_for_Composite_attractor[b]);
                        fprintf(f3,"%s = %lg\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->get_energy_of_state(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,(double)((N_trace+1)*N_RND),p_error));
                        fprintf(f4,"%s = %lg\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->get_energy_of_state(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,(double)((N_trace+1)*N_RND),p_error));
                        fprintf(f6,"%s = %d\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Unconnected->check_if_state_is_an_Attractor_State(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s));
                        fprintf(f7,"%s = %d\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->check_if_state_is_an_Attractor_State(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s));
                        fprintf(f8,"%s = ",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s);
                        for(int k=1;k<=BRN->N;k++) if(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s[k-1]!=48) fprintf(f8,"%s ",BRN->MDAT->Node_name[k]);
                        fprintf(f8,"\n");
                        for(int k=1;k<=gl->N;k++){
                            if(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s[gl->A[k]-1]!=48)
                                fprintf(f_Extras[k],"%s = 1\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s);
                            else fprintf(f_Extras[k],"%s = 0\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s);
                        }
                        
                        a=Unconnected->check_if_state_has_Basin_ID(Coupled->s);
                        if(a==0){ a=Unconnected->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Coupled->s);
                            //printf("Retured %ld\n",a);getchar();
                        }
                        fprintf(f1,"%s = %ld\n",Coupled->s,a);
                        fprintf(f5,"%s = %ld\n",Coupled->s,b);
                        //fprintf(f2,"%s = %d_%d\n",Coupled->s,
                        //        Unconnected_attractor_match_of_Composite_attractor[b],
                       //         Rank_in_matching_Unconnected_attractor_for_Composite_attractor[b]);
                        fprintf(f3,"%s = %lg\n",Coupled->s,Unconnected->get_energy_of_state(Coupled->s,(double)((N_trace+1)*N_RND),p_error));
                        fprintf(f4,"%s = %lg\n",Coupled->s,Coupled->get_energy_of_state(Coupled->s,(double)((N_trace+1)*N_RND),p_error));
                        fprintf(f6,"%s = %d\n",Coupled->s,Unconnected->check_if_state_is_an_Attractor_State(Coupled->s));
                        fprintf(f7,"%s = %d\n",Coupled->s,Coupled->check_if_state_is_an_Attractor_State(Coupled->s));
                        fprintf(f8,"%s = ",Coupled->s);
                        for(int k=1;k<=BRN->N;k++) if(Coupled->s[k-1]!=48) fprintf(f8,"%s ",BRN->MDAT->Node_name[k]);
                        fprintf(f8,"\n");
                        for(int k=1;k<=gl->N;k++){
                            if(Coupled->s[gl->A[k]-1]!=48)
                                fprintf(f_Extras[k],"%s = 1\n",Coupled->s);
                            else fprintf(f_Extras[k],"%s = 0\n",Coupled->s);
                        }
                        
                        
                        // taking care of all high-enough flux OTHER links out of this state
                        
                        for(int ipoz=1;ipoz<=BRN->N;ipoz++){ // iterates through the positions in the end of sync link state, flips one and uses the new state as alt end to weak link. Then flips position back.
                            if(Coupled->s[ipoz]=='1') Coupled->s[ipoz]='0'; else Coupled->s[ipoz]='1'; // flip the ith node of the next state
                            a1=Unconnected->check_if_state_was_sampled_enough(Coupled->s,V_MIN);
                            a2=Coupled->check_if_state_was_sampled_enough(Coupled->s,V_MIN);
                            if((a1>0)||(a2>0)){
                                lf1=Unconnected->get_link_flux(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,&ST_order1,p_error);
                                lf2=Coupled->get_link_flux(Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,&ST_order2,p_error);
                                lf1=lf1/(double)((N_trace+1)*N_RND);lf2=lf2/(double)((N_trace+1)*N_RND);
                                if((lf1>=p_error*minflux_fraction_of_p)||(lf2>=p_error*minflux_fraction_of_p)){
                                    fprintf(f,"%s\t->\t%s\t%lg\t%d\t%lg\t%d\t%lg\n",Coupled->Attractor_valleys[b]->Basin_Stored[i]->s,Coupled->s,(lf1+lf2)/2.,ST_order1,lf1,ST_order2,lf2);
                                    
                                    if(a1==0) a1=Unconnected->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Coupled->s);
                                    if(a2==0) a2=Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(Coupled->s);
                                    fprintf(f1,"%s = %ld\n",Coupled->s,a1);
                                    fprintf(f5,"%s = %ld\n",Coupled->s,a2);
                                    //fprintf(f2,"%s = %d_%d\n",Coupled->s,
                                    //        Unconnected_attractor_match_of_Composite_attractor[a2],
                                    //        Rank_in_matching_Unconnected_attractor_for_Composite_attractor[a2]);
                                    fprintf(f3,"%s = %lg\n",Coupled->s,Unconnected->get_energy_of_state(Coupled->s,(double)((N_trace+1)*N_RND),p_error));
                                    fprintf(f4,"%s = %lg\n",Coupled->s,Coupled->get_energy_of_state(Coupled->s,(double)((N_trace+1)*N_RND),p_error));
                                    fprintf(f6,"%s = %d\n",Coupled->s,Unconnected->check_if_state_is_an_Attractor_State(Coupled->s));
                                    fprintf(f7,"%s = %d\n",Coupled->s,Coupled->check_if_state_is_an_Attractor_State(Coupled->s));
                                    fprintf(f8,"%s = ",Coupled->s);
                                    for(int k=1;k<=BRN->N;k++) if(Coupled->s[k-1]!=48) fprintf(f8,"%s ",BRN->MDAT->Node_name[k]);
                                    fprintf(f8,"\n");
                                    for(int k=1;k<=gl->N;k++){
                                        if(Coupled->s[gl->A[k]-1]!=48)
                                            fprintf(f_Extras[k],"%s = 1\n",Coupled->s);
                                        else fprintf(f_Extras[k],"%s = 0\n",Coupled->s);
                                    }
                                }
                            }
                            if(Coupled->s[ipoz]=='1') Coupled->s[ipoz]='0'; else Coupled->s[ipoz]='1'; // flip the ith node of the next state back
                        }
                    }
                }
            }
        }
        fclose(f);
        fclose(f1);
        fclose(f2);
        fclose(f3);
        fclose(f4);
        fclose(f5);
        fclose(f6);
        fclose(f7);
        fclose(f8);
        for(int i=1;i<=gl->N;i++) fclose(f_Extras[i]);
    }
    
};


#endif
