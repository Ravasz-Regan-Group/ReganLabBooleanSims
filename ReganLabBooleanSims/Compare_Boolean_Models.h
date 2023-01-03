//
//  Compare_Boolean_Models.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 11/21/17.
//  Copyright © 2017 Regan_Group. All rights reserved.
//

#define AUTO    -1
#define OFF    0
#define ON      1
#define COMB    10
#define BOXN    10


#include "Model_is_Awesome_checklist.h"


unsigned int N_trace=5;
unsigned int N_RND=100;
double p_error=0.02;
unsigned int N_MAX=1000;


void build_both_Models_and_sample_state_space(const char modelname_REF[],
                                              const char modelname_MUT[],
                                              const char modules_name[]){
    Boolean_RegNet *A,*B;
    unsigned int N_trace,N_RND,N_MAX;
    MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,*PD_B;
    
    longint_ARRAY *gl_keep_asis;
    
    N_trace=5;N_RND=150;
    
    p_error=0.02;N_MAX=1000;
    
    A=new Boolean_RegNet(modelname_REF);
    B=new Boolean_RegNet(modelname_MUT);
    
    A->load_Module_IDs_and_Drawing_Order();
// int m1=A->load_Module_Drawing_order(modules_name);
//    A->load_Module_IDs(modules_name,m1);
    B->load_Module_IDs_and_Drawing_Order();
//    int m2=B->load_Module_Drawing_order(modules_name);
//    B->load_Module_IDs(modules_name,m2);
//
    gl_keep_asis=get_input_node_names(A);
    
    PD_A=new MODULAR_Boolean_Dynamics_PAIRED_TRACKER(A,N_MAX);
    PD_A->read_Node_Cytoskape_positions();
    PD_A->Coupled->record_timetraces_from_collection_of_random_inputs_per_inputBOX(N_trace,N_RND,p_error,gl_keep_asis);
  
    PD_B=new MODULAR_Boolean_Dynamics_PAIRED_TRACKER(B,N_MAX);
    PD_B->read_Node_Cytoskape_positions();
    PD_B->Coupled->record_timetraces_from_collection_of_random_inputs_per_inputBOX(N_trace,N_RND,p_error,gl_keep_asis);

    
    for(int i=1;i<=PD_A->Coupled->Attractor_NR;i++){
        write_link_for_All_environmental_changes_from_ONE_attractor(gl_keep_asis,i,PD_A);
    }
    
    for(int i=1;i<=PD_B->Coupled->Attractor_NR;i++){
        write_link_for_All_environmental_changes_from_ONE_attractor(gl_keep_asis,i,PD_B);
    }
    
    PD_A->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis->N)),p_error);
    PD_B->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis->N)),p_error);
    
    //  PD->Coupled->Export_ALL_Attractors_onto_Boolean_Network(0);
    //PD->Coupled->Export_ALL_Attractors_onto_Boolean_Network_ALIVE(0);
    //  PD->Coupled->Export_ALL_Attractors_onto_Boolean_Network_ALIVE_or_cycle(0);
    
   
    longint_ARRAY *l_sim_B;
    l_sim_B= Find_Similar_Attractors__Life_and_Death(PD_B,gl_keep_asis,1);
    delete l_sim_B;l_sim_B=NULL;
    
    Export_Attractor_Life_and_Death_table(PD_B,gl_keep_asis,1);
}


int *generate_best_canalyzing_value_choice(unsigned int A_SUB_index,Boolean_RegNet *A_SUB,Boolean_RegNet *A){
    int *can,*outside;
    snode::slinklist::slink *sc2;
    
    can=new int[A->N+1];
    outside=new int[A->N+1];
    for(int i=0;i<=A->N;i++) {can[i]=-1;outside[i]=0;} // not a neighbor
    
    unsigned int BRN_index=A_SUB->MDAT->BIG_ID[A_SUB_index];
    
    sc2=A->BN->node[BRN_index].Li.first_link_pt;
    while(sc2!=NULL){
        unsigned int subI_neigh= A_SUB->MDAT->get_node_ID_BIGID(sc2->sneighbor);
        if (subI_neigh == 0){
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
        // printf("i=%lld\t Canalyzing values set to:\n",i);
        for(int k=1;k<=A->N;k++)
            if(outside[k]>0){
                can[k]=s_outmodlinks[outside[k]-1];
                // printf("\t k=%d (%s) value = %d\n",k,BRN->MDAT->Node_name[k],can[k]);
            }
        longint_ARRAY *gl;
        gl=new longint_ARRAY();
        Boolean_Gate *G_test=A->get_PARTIAL_gate(BRN_index,can,gl);
        
       // printf("Partial gate returned is:\n");
       // A->print_PARTIAL_gate(BRN_index,gl,G_test);
        
        H_left_for_input[i]=G_test->Entropy_of_gate();
        if(G_test->test_function_for_all_inputs_return_bad_ones()->N>0){
            i_OK_input[i]=0;
           //printf("\t\t\t NO GO!\n");
        }
        else {i_OK_input[i]=1; // printf("\t\t\t All IN OK!\n");
        }
        // need to check how many of G_test's inputs are real,
        // order all combinations by largest real inputs
        // if more than 1 winner: select gate with largest H. if more than 1 winner: random choice.
        delete[] s_outmodlinks; s_outmodlinks=NULL;
        delete gl;gl=NULL;
    }
    double hmax=0;unsigned long long int i_choice=-1;
    for(unsigned long long int i=0;i<Ncomb;i++){
        if(i_OK_input[i]==1) if(H_left_for_input[i]>hmax) {hmax=H_left_for_input[i];i_choice=i;}
    }
    if (i_choice==-1) {
         printf("There are no way of setting the out-of-module inputs and still guaranteeing function of all internal links!\n");
        for(unsigned long long int i=0;i<Ncomb;i++)
            if(H_left_for_input[i]>hmax) {hmax=H_left_for_input[i];i_choice=i;}
        if(i_choice==-1) {
            printf("There is no regulation left for this node; locked in with all canal. inputs set to 0.\n");
            printf("\t k=%d (%s)\n",A_SUB_index,A_SUB->MDAT->Node_name[A_SUB_index]);
            i_choice=0;
           // getchar();
        }
    }
    
    s_outmodlinks=new unsigned short int[Ncomb+1];
    decimal_to_any_base(i_choice,2,outside[0], s_outmodlinks);
  //  printf("Final choice: %lld:\n",i_choice);
    for(int k=1;k<=A->N;k++)
        if(outside[k]>0){
            can[k]=s_outmodlinks[outside[k]-1];
            // printf("\t k=%d (%s) value = %d\n",k,BRN->MDAT->Node_name[k],can[k]);
        }
    delete[] s_outmodlinks; s_outmodlinks=NULL;
    
    s_outmodlinks=new unsigned short int[Ncomb+1];
    decimal_to_any_base(i_choice,2,outside[0], s_outmodlinks);
    for(int k=1;k<=A->N;k++)
        if(outside[k]>0){
            can[k]=s_outmodlinks[outside[k]-1];
            printf("\t k=%d (%s) canalyzing value settled at = %d\n",k,A->MDAT->Node_name[k],can[k]);
        }
    
    delete[] s_outmodlinks;s_outmodlinks=NULL;
    delete[] H_left_for_input;H_left_for_input=NULL;
    delete[] i_OK_input;i_OK_input=NULL;
    
    for(int k=1;k<=A->N;k++)
        if(outside[k]>0){
            if (k==A->MDAT->get_node_ID("SAHF")) {
                can[k]=0;
            }
        }
    delete[] outside;outside=NULL;
    return(can);
}

int *generate_LOCK_canalyzing_value(unsigned int A_SUB_index,Boolean_RegNet *A_SUB,Boolean_RegNet *A,int can_value){
    int *can,*outside;
    snode::slinklist::slink *sc2;
    
    can=new int[A->N+1];
    outside=new int[A->N+1];
    for(int i=0;i<=A->N;i++) {can[i]=-1;outside[i]=0;} // not a neighbor
    
    unsigned int BRN_index=A_SUB->MDAT->BIG_ID[A_SUB_index];
    
    sc2=A->BN->node[BRN_index].Li.first_link_pt;
    while(sc2!=NULL){
        unsigned int subI_neigh= A_SUB->MDAT->get_node_ID_BIGID(sc2->sneighbor);
        if (subI_neigh == 0){
            outside[0]++;
            outside[sc2->sneighbor]=outside[0];
            // printf("Link %s to %s comes from outside:\n", BRN->MDAT->Node_name[sc2->sneighbor],BRN->MDAT->Node_name[BRN_index]);
        }
        else can[sc2->sneighbor]=-2; // non-canalyzed neighbor
        sc2=sc2->next(sc2);
    }
    
    for(int k=1;k<=A->N;k++)
        if(outside[k]>0){
            can[k]=can_value;
        }
    longint_ARRAY *gl;
    gl=new longint_ARRAY();
    Boolean_Gate *G_test=A->get_PARTIAL_gate(BRN_index,can,gl);
    
    if(G_test->test_function_for_all_inputs_return_bad_ones()==0){
        printf("\t\t\t Partial gate in mutant network NOT OK:\n");
        A->print_PARTIAL_gate(BRN_index,gl,G_test);
    }
    else{ printf("\t\t\t Partial gate in mutant network OK -- all inputs functional!\n");
       }
    delete gl;gl=NULL;
    delete[] outside;outside=NULL;
    return(can);
}

Boolean_RegNet *generate_SUBSET_RNB(Boolean_RegNet *A, const char subset[],int subset_mode){
  //  char fn[300];
    snode::slinklist::slink *sc2;
    int *canalyze;
    longint_ARRAY *gl;
    Boolean_RegNet *A_SUB;
    
    A_SUB=new Boolean_RegNet(A,subset);
    
    for(int jsub=1;jsub<=A_SUB->N;jsub++){
        unsigned int jfull=A->MDAT->get_node_ID(A_SUB->MDAT->Node_name[jsub]);
        printf("Node_sub j=%d %s\t Node_full jf=%d %s\n",
               jsub,A_SUB->MDAT->Node_name[jsub],
               jfull,A->MDAT->Node_name[jfull]);
       // A->print_gate(jfull);
        
        switch (subset_mode) {
            case -1:
                canalyze=generate_best_canalyzing_value_choice(jsub,A_SUB,A);
                break;
            case 0:
                canalyze=generate_LOCK_canalyzing_value(jsub,A_SUB,A,0);
                break;
            case 1:
                canalyze=generate_LOCK_canalyzing_value(jsub,A_SUB,A,1);
                break;
            default:
                canalyze=generate_best_canalyzing_value_choice(jsub,A_SUB,A);
                break;
        }
        
        for(int l=1;l<=A->N;l++)
            if(canalyze[l]!=-1)
                 printf("\t <- %d %s : %d\n",l,A->MDAT->Node_name[l],canalyze[l]);
        
        
        sc2=A->BN->node[jfull].Li.first_link_pt;
        while(sc2!=NULL){
            if (sc2->sneighbor==0) {
                printf("Node %d has a neighbor 0 for some reason!\n",jfull);getchar();
            }
            unsigned int subI_neigh= A_SUB->MDAT->get_node_ID_BIGID(sc2->sneighbor);
            if (subI_neigh != 0){
                if(canalyze[sc2->sneighbor]>-1) {
                    printf("Error: in canalizing list of node %s:\nLink %d <- %lld is internal to module %d, but input %s is set to canalysing value %d\n",A->MDAT->Node_name[jfull],jfull,sc2->sneighbor, A->MDAT->Module_ID[jfull],A->MDAT->Node_name[sc2->sneighbor],canalyze[sc2->sneighbor]);//exit(1);
                    getchar();
                }
                else canalyze[sc2->sneighbor]=-2;
            }
            else {
                if(canalyze[sc2->sneighbor]==-1) {
                    printf("Error: in canalizing list of node %s:\nLink %d <- %lld is external to module %d, but the canalysing value of input %s (mod %d)is not set\n",A->MDAT->Node_name[jfull],jfull,sc2->sneighbor, A->MDAT->Module_ID[jfull],A->MDAT->Node_name[sc2->sneighbor],A->MDAT->Module_ID[sc2->sneighbor]);//exit(1);
                    getchar();
                }
            }
            sc2=sc2->next(sc2);
        }
        
        
        gl=new longint_ARRAY();
        
        Boolean_Gate *G_test=A->get_PARTIAL_gate(jfull,canalyze,gl);
       // A->print_PARTIAL_gate(jfull,gl,G_test);
       // getchar();
        
        for(unsigned long int l=0;l<gl->N;l++){
            unsigned int subI_neigh= A_SUB->MDAT->get_node_ID_BIGID(gl->A[gl->N-l]);
            A_SUB->BN->add_link(jsub,subI_neigh);
            A_SUB->K++;
        }
        if(gl->N==0) {
            A_SUB->BN->add_link(jsub,jsub);
        }
        
        A_SUB->Gate[jsub]=G_test;
        //A_SUB->print_gate(jsub); getchar();
        delete[] canalyze;canalyze=NULL;
        delete gl;gl=NULL;
    }
    
  //  getchar();
    return (A_SUB);
}

Boolean_RegNet *generate_MUTANT_RNB(Boolean_RegNet *A, const char mutantnode[],int setMutant){
    int *canalyze=NULL;
    longint_ARRAY *gl;
    Boolean_RegNet *A_SUB;
    
    A_SUB=new Boolean_RegNet(A,mutantnode,setMutant);
    
    unsigned int jmut=A->MDAT->get_node_ID(mutantnode);
    
    for(int jsub=1;jsub<=A_SUB->N;jsub++){
        unsigned int jfull=A->MDAT->get_node_ID(A_SUB->MDAT->Node_name[jsub]);
        printf("Node_sub j=%d %s\t Node_full jf=%d %s\n",
               jsub,A_SUB->MDAT->Node_name[jsub],
               jfull,A->MDAT->Node_name[jfull]);
        // A->print_gate(jfull);
        
        gl=new longint_ARRAY();
        Boolean_Gate *G_test=A->get_FULL_gate(jfull,gl);
        // A->print_PARTIAL_gate(jfull,gl,G_test);
        // getchar();
        
        for(unsigned long int l=0;l<gl->N;l++){
            unsigned int subI_neigh= A_SUB->MDAT->get_node_ID_BIGID(gl->A[gl->N-l]);
            A_SUB->BN->add_link(jsub,subI_neigh);
            A_SUB->K++;
        }
        if(gl->N==0) {
            A_SUB->BN->add_link(jsub,jsub);
        }
        
        A_SUB->Gate[jsub]=G_test;
        if(jfull==jmut) {
            if(setMutant==0) A_SUB->Gate[jsub]->knockout();
            if(setMutant==1) A_SUB->Gate[jsub]->overexpress();
            //A_SUB->print_gate(jsub); getchar();
        }
        delete[] canalyze;canalyze=NULL;
        delete gl;gl=NULL;
    }
    return (A_SUB);
}

void run_cell_cycle_entry_exit(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_B,
                               longint_ARRAY *gl_keep_asis_A,
                               longint_ARRAY *gl_keep_asis_B,int cycles_norm){
    longint_ARRAY *l;
    int na;
    
    
    l=get_Healthy_G0_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            int pause = 0;
            na=0;
            while((na<=cycles_norm)&&(pause<=100)){
                pause++;
                na= count_cell_cycles_following_OneNode_pulse_experiment_MOD ("GF_High",(int)(l->A[j]),1,PD_A,pause);
            }
            run_OneNode_pulse_experiment_MOD("GF_High",(int)(l->A[j]),1,PD_A,0.000000001,pause-1);
            run_OneNode_pulse_experiment_MOD("GF_High",(int)(l->A[j]),1,PD_A,0.000000001,pause);
        }
    }
    delete l;l=NULL;
    
    l=get_Healthy_G0_attractors(PD_B);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            int pause = 0;
            na=0;
            while((na<=cycles_norm)&&(pause<=100)){
                pause++;
                na= count_cell_cycles_following_OneNode_pulse_experiment_MOD ("GF_High",(int)(l->A[j]),1,PD_B,pause);
            }
            run_OneNode_pulse_experiment_MOD("GF_High",(int)(l->A[j]),1,PD_B,0.000000001,pause-1);
            run_OneNode_pulse_experiment_MOD("GF_High",(int)(l->A[j]),1,PD_B,0.000000001,pause);
        }
    }
    delete l;l=NULL;
}

void run_apoptosis_from_G0(const char nodename[], MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_B,
                               longint_ARRAY *gl_keep_asis_A,
                               longint_ARRAY *gl_keep_asis_B){
    longint_ARRAY *l;
    
    l=get_Healthy_G0_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            int dead=chech_time_of_death_following_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_A);
            if(dead>0){
                run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_A,0.000000001,dead-1);
                run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_A,0.000000001,dead);
                //run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_A,0.000000001,dead+1);
            }
            else {printf("Full model: cell does not die out of G0 in respose to %s!",nodename);getchar();}
        }
    }
    delete l;l=NULL;
    
    l=get_Healthy_G0_attractors(PD_B);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            int dead=chech_time_of_death_following_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_B);
            if(dead>0){
               run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_B,0.000000001,dead-1);
                run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_B,0.000000001,dead);
               // run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_B,0.000000001,dead+1);
            }
            else {printf("Small model: cell does not die out of G0 in respose to %s!",nodename);getchar();}
        }
    }
    delete l;l=NULL;
}

/*
void run_cell_cycle_arrest(const char nodename[],MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_B,
                               longint_ARRAY *gl_keep_asis_A,
                               longint_ARRAY *gl_keep_asis_B){
    longint_ARRAY *l;
    FILE *f;
    Attractors_on_MAP *AMAP;
    double X_BOX=20.,Y_BOX=20.;
    Boolean_Dynamics_TRACKER *D;
    int ll_now;
    char bla[300];
    
    l=get_Healthy_CellCycle_attractors(PD_A,gl_keep_asis_A);
    unsigned long int pulse_ID = PD_A->BRN->MDAT->get_node_ID(nodename);
    ll_now=0;
    for(unsigned int ll=1;ll<=gl_keep_asis_A->N;ll++){
        if (pulse_ID == gl_keep_asis_A->A[ll]) {
            ll_now=ll;
        }
    }
    if (ll_now==0) {
        printf("Pulse of %s does now work, %s is not one of the pre-set imput nodes!\n",nodename,nodename);
        return;
    }
    unsigned long int cytokin = PD_A->BRN->MDAT->get_node_ID("4N_DNA");
    if (cytokin == 0) cytokin = PD_A->BRN->MDAT->get_node_ID("f4N_DNA");
    unsigned long int GF_id = PD_A->BRN->MDAT->get_node_ID("GF");
    if(cytokin<=0) return;
        
    D=PD_A->Coupled;
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            int pmax=0,pmin=-1;
            printf("j=%d l->A[j]=%ld\t",j, l->A[j]);
            if(D->Attr_Transitions[ll_now]!=NULL){
                snode::slinklist::slink *sc2;
                    sc2=D->Attr_Transitions[ll_now]->node[l->A[j]].Li.first_link_pt;
                    while(sc2!=NULL){
                        printf("\tlink from attractor %ld to %llu :ID = %d  label = %d (pulseid= %ld)  w= %lg\n",
                               l->A[j],sc2->sneighbor,
                               sc2->link_props->ID,
                               sc2->link_props->distinc_link_labels,pulse_ID,
                               sc2->link_props->lw);//getchar();
                        if(sc2->link_props->distinc_link_labels == pulse_ID){
                            if(sc2->link_props->lw > pmax) pmax=(int)sc2->link_props->lw;
                            if(pmin<0) pmin=(int)sc2->link_props->lw;
                            if(sc2->link_props->lw < pmin) pmin=(int)sc2->link_props->lw;
                        }
                        sc2=sc2->next(sc2);
                    }
            }
            else {printf("PD_A->Coupled->Attr_Transitions[%d] does not exist!!!\n",ll_now);getchar();}
           
            char fn[400];
            sprintf(fn,"%s/%s/_EXP/%s_CellCycle_dept_Response_of_ATTR-%ld.eps",MODEL_Directory,PD_A->BRN->path,nodename,l->A[j]);
            printf("%s\npmax = %d",fn,pmax);
            if (pmax==0) pmax=20;
            f=fopen(fn,"w");
            EpsInit(f,-5,-5,X_BOX*(pmax+2),Y_BOX*(D->Attractor_valleys[l->A[j]]->N_a+2));
            EpsSetFont(f,"Times-Roman",5);
            int ind_G1st=0;
            for(unsigned int i=2;i<=D->Attractor_valleys[l->A[j]]->N_a;i++){
                if ((D->Attractor_valleys[l->A[j]]->Basin_Stored[D->Attractor_valleys[l->A[j]]->Attractor[i-1]]->s[cytokin-1]==49)&&
                    (D->Attractor_valleys[l->A[j]]->Basin_Stored[D->Attractor_valleys[l->A[j]]->Attractor[i]]->s[cytokin-1]==48))
                    ind_G1st=i;
            }
            if ((D->Attractor_valleys[l->A[j]]->Basin_Stored[D->Attractor_valleys[l->A[j]]->Attractor[D->Attractor_valleys[l->A[j]]->N_a]]->s[cytokin-1]==49)&&
                (D->Attractor_valleys[l->A[j]]->Basin_Stored[D->Attractor_valleys[l->A[j]]->Attractor[1]]->s[cytokin-1]==48))
                ind_G1st=1;
          //  if(pmax==0) pmax=D->Attractor_valleys[l->A[j]]->N_a+2;
           
            double y;
            for (unsigned int pp=1; pp<=pmax; pp++) {
               // if(((pp==pmin)&&(pmin>1))||(pp==pmax))
               //     run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_A,0.000000001,pp);
                
                sprintf(bla, "%d",pp);
                EpsDrawString(f,0, pp*X_BOX,(D->Attractor_valleys[l->A[j]]->N_a+1)*Y_BOX, bla);
                for(int i=ind_G1st;i<=D->Attractor_valleys[l->A[j]]->N_a;i++){
                    unsigned long int bback=attractor_after_pulse(PD_A,l->A[j],i,pulse_ID,GF_id,pp);
                    AMAP = new Attractors_on_MAP(PD_A,bback);
                    y=D->Attractor_valleys[l->A[j]]->N_a - (i-ind_G1st);
                    draw_attractor_symbol(f, AMAP,1,pp*X_BOX,y*Y_BOX,30,30,0.2);
                    delete AMAP; AMAP=NULL;
                    sprintf(bla, "%d",i-ind_G1st);
                    EpsDrawString(f,0, 0 ,y*Y_BOX, bla);
                }
                for(int i=1;i<ind_G1st;i++){
                    unsigned long int bback=attractor_after_pulse(PD_A,l->A[j],i,pulse_ID,GF_id,pp);
                    AMAP = new Attractors_on_MAP(PD_A,bback);
                    y=D->Attractor_valleys[l->A[j]]->N_a - (D->Attractor_valleys[l->A[j]]->N_a + (i-ind_G1st));
                    draw_attractor_symbol(f, AMAP,1,pp*X_BOX,y*Y_BOX,30,30,0.2);
                    delete AMAP; AMAP=NULL;
                    sprintf(bla, "%d",(int)D->Attractor_valleys[l->A[j]]->N_a + (i-ind_G1st));
                    EpsDrawString(f,0, 0 ,y*Y_BOX, bla);
                }
            }
            run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_A,0.000000001,100);
            
            sprintf(fn,"%s/%s/_EXP/%s_CellCycle_dept_Response_of_ATTR-%ld",MODEL_Directory_system,PD_A->BRN->path,nodename,l->A[j]);
            close_EPSDRAW_file(f,fn,2);
        }
    }
    delete l;l=NULL;
}
 */

void run_cell_cycle_arrest(const char nodename[],MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                           longint_ARRAY *gl_keep_asis_A){
    longint_ARRAY *l;
    FILE *f;
    Attractors_on_MAP *AMAP;
    double X_BOX=20.,Y_BOX=20.;
    Boolean_Dynamics_TRACKER *D;
    int ll_now;
    char bla[300];
    
    l=get_Healthy_CellCycle_attractors(PD_A);
    unsigned long int pulse_ID = PD_A->BRN->MDAT->get_node_ID(nodename);
    ll_now=0;
    for(unsigned int ll=1;ll<=gl_keep_asis_A->N;ll++){
        if (pulse_ID == gl_keep_asis_A->A[ll]) {
            ll_now=ll;
        }
    }
    if (ll_now==0) {
        printf("Pulse of %s does now work, %s is not one of the pre-set imput nodes!\n",nodename,nodename);
        return;
    }
    unsigned long int cytokin = PD_A->BRN->MDAT->get_node_ID("4N_DNA");
    if (cytokin == 0) cytokin = PD_A->BRN->MDAT->get_node_ID("f4N_DNA");

    unsigned long int GF_id = PD_A->BRN->MDAT->get_node_ID(GF_low_name);
    unsigned long int CD_id = PD_A->BRN->MDAT->get_node_ID(CD_low_name);
    if(cytokin<=0) return;
    
    D=PD_A->Coupled;
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            int pmax=0,pmin=-1;
            printf("j=%d l->A[j]=%ld\t",j, l->A[j]);
            if(D->Attr_Transitions[ll_now]!=NULL){
                snode::slinklist::slink *sc2;
                sc2=D->Attr_Transitions[ll_now]->node[l->A[j]].Li.first_link_pt;
                while(sc2!=NULL){
                    printf("\tlink from attractor %ld to %llu :ID = %d  label = %d (pulseid= %ld)  w= %lg\n",
                           l->A[j],sc2->sneighbor,
                           sc2->link_props->ID,
                           sc2->link_props->distinc_link_labels,pulse_ID,
                           sc2->link_props->lw);//getchar();
                    if(sc2->link_props->distinc_link_labels == pulse_ID){
                        if(sc2->link_props->lw > pmax) pmax=(int)sc2->link_props->lw;
                        if(pmin<0) pmin=(int)sc2->link_props->lw;
                        if(sc2->link_props->lw < pmin) pmin=(int)sc2->link_props->lw;
                    }
                    sc2=sc2->next(sc2);
                }
            }
            else {printf("PD_A->Coupled->Attr_Transitions[%d] does not exist!!!\n",ll_now);getchar();}
            
            char fn[400];
            sprintf(fn,"%s/%s/_EXP/%s_CellCycle_dept_Response_of_ATTR-%ld.eps",MODEL_Directory,PD_A->BRN->path,nodename,l->A[j]);
            printf("%s\npmax = %d",fn,pmax);
            if (pmax==0) pmax=20;
            f=fopen(fn,"w");
            EpsInit(f,-5,-5,X_BOX*(pmax+2),Y_BOX*(D->Attractor_valleys[l->A[j]]->N_a+2));
            EpsSetFont(f,"Times-Roman",5);
            int ind_G1st=0;
            for(unsigned int i=2;i<=D->Attractor_valleys[l->A[j]]->N_a;i++){
                if ((D->Attractor_valleys[l->A[j]]->Basin_Stored[D->Attractor_valleys[l->A[j]]->Attractor[i-1]]->s[cytokin-1]==49)&&
                    (D->Attractor_valleys[l->A[j]]->Basin_Stored[D->Attractor_valleys[l->A[j]]->Attractor[i]]->s[cytokin-1]==48))
                    ind_G1st=i;
            }
            if ((D->Attractor_valleys[l->A[j]]->Basin_Stored[D->Attractor_valleys[l->A[j]]->Attractor[D->Attractor_valleys[l->A[j]]->N_a]]->s[cytokin-1]==49)&&
                (D->Attractor_valleys[l->A[j]]->Basin_Stored[D->Attractor_valleys[l->A[j]]->Attractor[1]]->s[cytokin-1]==48))
                ind_G1st=1;
            //  if(pmax==0) pmax=D->Attractor_valleys[l->A[j]]->N_a+2;
            
            double y;
            for (unsigned int pp=1; pp<=pmax; pp++) {
                // if(((pp==pmin)&&(pmin>1))||(pp==pmax))
                //     run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_A,0.000000001,pp);
                
                sprintf(bla, "%d",pp);
                EpsDrawString(f,0, pp*X_BOX,(D->Attractor_valleys[l->A[j]]->N_a+1)*Y_BOX, bla);
                for(int i=ind_G1st;i<=D->Attractor_valleys[l->A[j]]->N_a;i++){
                    unsigned long int bback=attractor_after_pulse(PD_A,l->A[j],i,pulse_ID,GF_id,CD_id,pp);
                    AMAP = new Attractors_on_MAP(PD_A,bback);
                    y=D->Attractor_valleys[l->A[j]]->N_a - (i-ind_G1st);
                    draw_attractor_symbol(f, AMAP,1,pp*X_BOX,y*Y_BOX,30,30,0.2);
                    delete AMAP; AMAP=NULL;
                    sprintf(bla, "%d",i-ind_G1st);
                    EpsDrawString(f,0, 0 ,y*Y_BOX, bla);
                }
                for(int i=1;i<ind_G1st;i++){
                    unsigned long int bback=attractor_after_pulse(PD_A,l->A[j],i,pulse_ID,GF_id,CD_id,pp);
                    AMAP = new Attractors_on_MAP(PD_A,bback);
                    y=D->Attractor_valleys[l->A[j]]->N_a - (D->Attractor_valleys[l->A[j]]->N_a + (i-ind_G1st));
                    draw_attractor_symbol(f, AMAP,1,pp*X_BOX,y*Y_BOX,30,30,0.2);
                    delete AMAP; AMAP=NULL;
                    sprintf(bla, "%d",(int)D->Attractor_valleys[l->A[j]]->N_a + (i-ind_G1st));
                    EpsDrawString(f,0, 0 ,y*Y_BOX, bla);
                }
            }
            run_OneNode_pulse_experiment_MOD(nodename,(int)(l->A[j]),1,PD_A,0.000000001,100);
            
            sprintf(fn,"%s/%s/_EXP/%s_CellCycle_dept_Response_of_ATTR-%ld",MODEL_Directory_system,PD_A->BRN->path,nodename,l->A[j]);
            close_EPSDRAW_file(f,fn,2);
        }
    }
    delete l;l=NULL;
}


unsigned long int states_match(Boolean_Dynamics_TRACKER *A,Boolean_Dynamics_TRACKER *A_sub){
    unsigned long int O=0;
    for (unsigned long int i=1;i<=A_sub->BRN->N;i++){
        if(A_sub->s[i-1]==A->s[A_sub->BRN->MDAT->BIG_ID[i]-1]) O++;
    }
    return(O);
}

void compare_attractor_graphs(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                              MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_B,
                              longint_ARRAY *gl_keep_asis_A,
                              longint_ARRAY *gl_keep_asis_B){
    Phenotype_now P1,P2;
    
    if(gl_keep_asis_B->N!=gl_keep_asis_A->N){
        printf("Number of environment-setting variables need to match!\n"); exit(1);
    }
    for(unsigned long int l=1;l<=gl_keep_asis_B->N;l++)
        if(PD_B->BRN->MDAT->BIG_ID[gl_keep_asis_B->A[l]]!= gl_keep_asis_A->A[l]){
            printf("Order of environment-setting variables need to match!\n\t Mismatch at poz %ld: %s (index %ld) of full network vs.  %s (index %ld, bigID %d) of subset\n",l,PD_A->Coupled->BRN->MDAT->Node_name[gl_keep_asis_A->A[l]],gl_keep_asis_A->A[l],PD_B->Coupled->BRN->MDAT->Node_name[gl_keep_asis_B->A[l]],gl_keep_asis_B->A[l],PD_B->BRN->MDAT->BIG_ID[gl_keep_asis_B->A[l]]); exit(1);
        }
    
    for(unsigned long int ia=1;ia<=PD_A->Coupled->Attractor_NR;ia++){
        P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD_A,ia,1);
        for (unsigned int ll=1; ll<=gl_keep_asis_A->N; ll++) {
            PD_A->Coupled->Attr_Transitions[ll]->node[ia].ID = 0;
            PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw =0;
            for(unsigned long int ib=1;ib<=PD_B->Coupled->Attractor_NR;ib++){
                if (PD_A->Coupled->Attr_Transitions[1]->node[ia].node_props->basin == PD_B->Coupled->Attr_Transitions[1]->node[ib].node_props->basin) {
                    P2=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD_B,ib,1);
                    if(identical_phenotype(P1,P2)){
                         for (unsigned int l2=1; l2<=gl_keep_asis_B->N; l2++)
                             PD_B->Coupled->Attr_Transitions[l2]->node[ib].ID = ia;
                        // mapped to main network attractor with identical phenotype (may not be unique)
                        if ((PD_A->Coupled->Attractor_valleys[ia]->N_a == 1) && (PD_B->Coupled->Attractor_valleys[ib]->N_a == 1)) {
                            PD_A->Coupled->set_state(PD_A->Coupled->Attractor_valleys[ia]->Basin_Stored[PD_A->Coupled->Attractor_valleys[ia]->Attractor[1]]->s);
                            PD_B->Coupled->set_state(PD_B->Coupled->Attractor_valleys[ib]->Basin_Stored[PD_B->Coupled->Attractor_valleys[ib]->Attractor[1]]->s);
                            unsigned long int nm=states_match(PD_A->Coupled,PD_B->Coupled);
                            for (unsigned int l2=1; l2<=gl_keep_asis_B->N; l2++)
                                PD_B->Coupled->Attr_Transitions[l2]->node[ib].node_props->nw =nm;
                                // weight matches number of identical nodes
                            if (nm>PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw){
                                PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw =nm;
                                PD_A->Coupled->Attr_Transitions[ll]->node[ia].ID = ib;
                            }
                        }
                        else {
                            if(PD_A->Coupled->Attractor_valleys[ia]->N_a == PD_B->Coupled->Attractor_valleys[ib]->N_a){
                                for (unsigned int l2=1; l2<=gl_keep_asis_B->N; l2++)
                                    PD_B->Coupled->Attr_Transitions[l2]->node[ib].node_props->nw = - 1.;
                                if ((PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw < 1 )){
                                        PD_A->Coupled->Attr_Transitions[ll]->node[ia].ID = ib;
                                        PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw = - 1;
                                }
                            }
                            else    {
                                for (unsigned int l2=1; l2<=gl_keep_asis_B->N; l2++)
                                    PD_B->Coupled->Attr_Transitions[l2]->node[ib].node_props->nw = - 100.;
                                if  ((PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw < 1 ) &&
                                    (PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw > - 1 )){
                                        PD_A->Coupled->Attr_Transitions[ll]->node[ia].ID = ib;
                                        PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw = -100;
                                }
                                else    {
                                    for (unsigned int l2=1; l2<=gl_keep_asis_B->N; l2++)
                                        PD_B->Coupled->Attr_Transitions[l2]->node[ib].node_props->nw -=  100.;
                                    if  (PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw < - 99 ){
                                         PD_A->Coupled->Attr_Transitions[ll]->node[ia].ID = -1;
                                         PD_A->Coupled->Attr_Transitions[ll]->node[ia].node_props->nw -= 100;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void run_paper_figs(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                    MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_B,
                    longint_ARRAY *gl_keep_asis_A,
                    longint_ARRAY *gl_keep_asis_B){
    
//    — Fig. 3,4 - check for pre-commitment (time course series, output before and after pre-commitment)
 
  //  run_cell_cycle_entry_exit(PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B,3);
  //  run_cell_cycle_entry_exit(PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B,1);
 //   run_cell_cycle_entry_exit(PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B,0);

 //   run_apoptosis_from_G0("Trail",PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
 //   run_apoptosis_from_G0("UV",PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
 //   run_apoptosis_from_G0("GF",PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
 //   run_apoptosis_from_G0("ROX_ext",PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
    
 
   // run_cell_cycle_arrest("GF_High",PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
 //   run_cell_cycle_arrest("Trail",PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
 
  //  run_cell_cycle_arrest("Gamma",PD_A,gl_keep_asis_A);
  //  run_cell_cycle_arrest("Gamma",PD_B,gl_keep_asis_B);
    
   // run_transition_from_cycle(20,"GF_High",PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);

    run_OneNode_pulse_experiment_MOD("GF_High",20,1,PD_A,0.000000001,130);

   run_OneNode_pulse_experiment_MOD("GF_High",16,1,PD_B,0.000000001,130);
//    run_OneNode_pulse_experiment_MOD("GF_High",35,1,PD_B,0.000000001,130);
//    run_OneNode_pulse_experiment_MOD("GF_High",36,1,PD_B,0.000000001,130);
//    run_OneNode_pulse_experiment_MOD("GF_High",35,1,PD_B,0.000000001,130);
//    run_OneNode_pulse_experiment_MOD("GF_High",30,1,PD_B,0.000000001,130);
//    run_OneNode_pulse_experiment_MOD("GF_High",31,1,PD_B,0.000000001,130);
  
    //  run_OneNode_pulse_experiment_MOD("GF_High",77,1,PD_B,0.000000001,130);
 //  run_OneNode_pulse_experiment_MOD("GF_High",16,1,PD_B,0.000000001,130);
  //  run_OneNode_pulse_experiment_MOD("GF_High",26,1,PD_B,0.000000001,100);
  //  run_OneNode_pulse_experiment_MOD("GF_High",22,1,PD_B,0.000000001,100);

    
//    run_transition_from_cycle(18,"GF_High",PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
//   run_OneNode_pulse_experiment_MOD("Gamma",15,1,PD_A,0.000000001,100);
//    run_OneNode_pulse_experiment_MOD("Gamma",28,1,PD_A,0.000000001,100);
//   run_OneNode_pulse_experiment_MOD("Gamma",21,1,PD_A,0.000000001,120);
 //  run_OneNode_pulse_experiment_MOD("Gamma",54,1,PD_A,0.000000001,120);
 //   run_OneNode_pulse_experiment_MOD("Gamma",51,1,PD_A,0.000000001,120);
  //  run_OneNode_pulse_experiment_MOD("Trail",14,1,PD_A,0.000000001,120);
   //   run_OneNode_pulse_experiment_MOD("Trail",15,1,PD_A,0.000000001,120);

  //   run_OneNode_pulse_experiment_MOD("GF_High",16,1,PD_B,0.000000001,120);
//     run_OneNode_pulse_experiment_MOD("GF_High",35,1,PD_B,0.000000001,120);
 //   run_OneNode_pulse_experiment_MOD("GF_High",37,1,PD_B,0.000000001,120);

//    run_OneNode_pulse_experiment_MOD("GF_High",38,1,PD_A,0.000000001,120);
 //   run_OneNode_pulse_experiment_MOD("Trail",14,1,PD_B,0.000000001,120);
 //   run_OneNode_pulse_experiment_MOD("Trail",18,1,PD_B,0.000000001,120);
 //   run_OneNode_pulse_experiment_MOD("Trail",19,1,PD_B,0.000000001,120);
 //   run_OneNode_pulse_experiment_MOD("Trail",22,1,PD_B,0.000000001,120);
}

MODULAR_Boolean_Dynamics_PAIRED_TRACKER *build_and_sample(Boolean_RegNet *A, longint_ARRAY *gl_keep_asis_A, const char modules_name[],int MODULES_ATTR,const char DIRName[],int subm){
    MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A;
    
    
    PD_A=new MODULAR_Boolean_Dynamics_PAIRED_TRACKER(A,N_MAX);
    strcpy(PD_A->Coupled->DIRname,DIRName);
    
    if(subm==0) PD_A->read_Node_Cytoskape_positions();
    PD_A->Coupled->record_timetraces_from_collection_of_random_inputs_per_inputBOX(N_trace,N_RND,p_error,gl_keep_asis_A);
    
    
    for(int i=1;i<=PD_A->Coupled->Attractor_NR;i++)
        write_link_for_All_environmental_changes_from_ONE_attractor(gl_keep_asis_A,i,PD_A);
    
    PD_A->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis_A->N)),p_error);
    PD_A->BRN->export_Boolean_RegNet_with_names_to_Cytoskape();
    char fn[400];
    //sprintf(fn, "%s/%s/%s_GATES.txt",MODEL_Directory,PD_A->Coupled->DIRname,PD_A->Coupled->SAMPLEname);
    PD_A->Coupled->BRN->print_gates_to_File();
    
    sprintf(fn, "mkdir %s/%s/%s\n",MODEL_Directory_system,PD_A->Coupled->DIRname,PD_A->Coupled->SAMPLEname);
    system(fn);
    sprintf(fn, "mkdir %s/%s/%s/tt\n",MODEL_Directory_system,PD_A->Coupled->DIRname,PD_A->Coupled->SAMPLEname);
    system(fn);
    sprintf(fn, "%s/%s/%s",MODEL_Directory,PD_A->Coupled->DIRname,PD_A->Coupled->SAMPLEname);
    PD_A->Coupled->BRN->print_gates_to_tt_folder(fn);
    
    longint_ARRAY *yztu_inputs_A;
    yztu_inputs_A=get_input_nodes_for_5input_drawing(A);
    for (int i=1; i<=gl_keep_asis_A->N; i++) {
        make_3D_map_of_steady_states_new(PD_A,yztu_inputs_A,gl_keep_asis_A,i,0);
    }
    PD_A->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis_A->N)),p_error);
    run_OnePulse_analysis(gl_keep_asis_A,PD_A);
    PD_A->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis_A->N)),p_error);
    
    PD_A->BRN->export_Boolean_RegNet_with_names_to_Cytoskape();
    longint_ARRAY *l_sim_A;
    l_sim_A= Find_Similar_Attractors__Life_and_Death(PD_A,gl_keep_asis_A,1);
    delete l_sim_A;l_sim_A=NULL;
    Export_Attractor_Life_and_Death_table(PD_A,gl_keep_asis_A,1);
    
    longint_ARRAY *l;
    l=get_Healthy_G0_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
          //  run_OneNode_pulse_experiment_MOD("GF_High",(int)l->A[j],1,PD_A,0.000000001,120);
            //          run_OneNode_pulse_experiment_MOD("Gamma",l->A[j],1,PD_A,0.000000001,120);
            //          run_OneNode_pulse_experiment_MOD("GF",l->A[j],1,PD_A,0.000000001,120);
            //          run_OneNode_pulse_experiment_MOD("Trail",l->A[j],1,PD_A,0.000000001,120);
        }
    }
    delete l; l=NULL;
    
    if(MODULES_ATTR==1){
            set_up_module_MetaData(PD_A);
        
//            work_with_individual_Modules(1,PD_A,0.0000001); // R-SW
//            work_with_individual_Modules(2,PD_A,0.0000001); // PH-SW
//            work_with_individual_Modules(9,PD_A,0.0000001); // ORC
//            work_with_individual_Modules(4,PD_A,0.0000001); // Apoptosis
           //  work_with_individual_Modules(5,PD_A,0.0000001); // AKT
//
    
 //           work_with_individual_Modules(8,PD_A,0.0000001); // Replication, 4N
        
        //   work_with_individual_Modules(3,PD_A,0.0000001); // DNA damage core
          //  work_with_individual_Modules(13,PD_A,0.0000001); // DNA repair stuff
  //          work_with_individual_Modules(6,PD_A,0.0000001); // Mitochondria
  //          work_with_individual_Modules(7,PD_A,0.0000001);// CH senesence
           // work_with_individual_Modules(8,PD_A,0.0000001); // Replication, 4N
    //        work_with_individual_Modules(3,PD_A,0.0000001);  // DNA damage
    //        work_with_individual_Modules(11,PD_A,0.0000001); // p53 kill
        
        // exit(1);
    }
    return (PD_A);
}

void run_cell_cycle_Delayed_KO(const char KOname[],
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                               longint_ARRAY *gl_keep_asis_A){
    longint_ARRAY *l;
    
    l=get_Healthy_CellCycle_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_OneNode_Delayed_KO_experiment_MOD(KOname,(int)l->A[j],1,PD_A,0);
        }
    }
}

void run_cell_cycle_Delayed_OE(const char KOname[],
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                               longint_ARRAY *gl_keep_asis_A){
    longint_ARRAY *l;
    
    l=get_Healthy_CellCycle_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_OneNode_Delayed_OE_experiment_MOD(KOname,(int)l->A[j],1,PD_A,0);
        }
    }
}

void run_cell_cycle_Delayed_HITS(longint_ARRAY_pair *Hits,
                              MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                              longint_ARRAY *gl_keep_asis_A){
    longint_ARRAY *l;
    
    l=get_Healthy_CellCycle_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_Delayed_Hits_experiment_MOD(Hits,(int)l->A[j],1,PD_A,0);
        }
    }
}


void run_stochastic_cell_cycle(const char InputNode[], double p_input, longint_ARRAY_pair *Hits,
                                 MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                 longint_ARRAY *gl_keep_asis_A, int T_max_cc){
    longint_ARRAY *l;
    Cell_cycle_stats_in_run STnow;
    char fn[400],fn2[400];
    FILE *f;
    unsigned long int G1,R,G2,M,T;
    
    sprintf(fn2, "Hits");
    if((Hits != NULL) && (Hits->N>0)){
        for (unsigned long int l=1; l<=Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD_A->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
        }
    }
    else sprintf(fn2,"%s_WILDTYPE",fn2);
    
    unsigned long int id_input=PD_A->BRN->MDAT->get_node_ID(InputNode);
    l=get_Healthy_CellCycle_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            Cell_cycle_stats_in_run Statsnow;
            Statsnow.N_cycles=0;
            Statsnow.N_endoredupl_T=0;
            Statsnow.N_endoredupl_G2=0;
            Statsnow.N_endoredupl_Aneup=0;
            Statsnow.T_live=0;
            Statsnow.Time_G1=0; G1=0;
            Statsnow.Time_R=0; R=0;
            Statsnow.Time_G2=0; G2=0;
            Statsnow.Time_M=0; M=0;
            Statsnow.Time_Telo=0; T=0;
            for (int t=0;t<=MAX_ARREST; t++){
                Statsnow.Dist_G1[t]=0;
                Statsnow.Dist_G2[t]=0;
                Statsnow.Dist_Mitosis[t]=0;
                Statsnow.Dist_Replication[t]=0;
                Statsnow.Dist_Telophase[t]=0;
            }
            
            run_stochastic_Input_experiment_MOD(id_input,p_input,Hits,(int)l->A[j],1,PD_A,0);
           
            unsigned long long int t_live_total =0; int deaths =0;
            while (t_live_total<T_max_cc) {
                STnow=count_cell_cycles_stochastic_input(InputNode,p_input,Hits,(int)l->A[j],1,PD_A,0,T_max_cc);
                t_live_total+=STnow.T_live;
                if (STnow.T_live<T_max_cc) {
                    deaths++;
                }
                
                Statsnow.N_cycles+=STnow.N_cycles;
                Statsnow.N_endoredupl_T+=STnow.N_endoredupl_T;
                Statsnow.N_endoredupl_G2+=STnow.N_endoredupl_G2;
                Statsnow.N_endoredupl_Aneup+=STnow.N_endoredupl_Aneup;
                for (int t=0;t<=MAX_ARREST; t++){
                    Statsnow.Dist_G1[t]+=STnow.Dist_G1[t];
                        if(STnow.Dist_G1[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_G2[t]+=STnow.Dist_G2[t];
                        if(STnow.Dist_G2[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_Mitosis[t]+=STnow.Dist_Mitosis[t];
                        if(STnow.Dist_Mitosis[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_Replication[t]+=STnow.Dist_Replication[t]; if(STnow.Dist_Replication[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_Telophase[t]+=STnow.Dist_Telophase[t];
                        if(STnow.Dist_Telophase[t]>0) Statsnow.MAX_INTERVAL=t;
                 }
                Statsnow.T_live+=STnow.T_live;
                Statsnow.Time_G1+=STnow.Time_G1;
                Statsnow.Time_R+=STnow.Time_R;
                Statsnow.Time_G2+=STnow.Time_G2;
                Statsnow.Time_M+=STnow.Time_M;
                Statsnow.Time_Telo+=STnow.Time_Telo;
                
            }
            for (int t=0;t<=MAX_ARREST; t++){
                G1+=Statsnow.Dist_G1[t];
                R+=Statsnow.Dist_Replication[t];
                G2+=Statsnow.Dist_G2[t];
                M+=Statsnow.Dist_Mitosis[t];
                T+=Statsnow.Dist_Telophase[t];
            }
            
            sprintf(fn,"%s/%s/_EXP/%s/%s_Stochastic_%s-%.2lf__A%d_As-%d__STATS_%d-CC_%d-Dead.txt",MODEL_Directory,PD_A->BRN->path,fn2,fn2,InputNode,p_input,(int)l->A[j],1,Statsnow.N_cycles,deaths);
           // printf("%s\n",fn);
            printf("%lg\t%d\t%d\t%lld\n",p_input,Statsnow.N_cycles,deaths,Statsnow.T_live);
            
            f=fopen(fn,"w");
            for (int t=0;t<=Statsnow.MAX_INTERVAL; t++) {
                fprintf(f, "%d\t%lg\t%lg\t%lg\t%lg\t%lg\n",t,Statsnow.Dist_G1[t]/(double)G1,Statsnow.Dist_Replication[t]/(double)R,Statsnow.Dist_G2[t]/(double)G2,Statsnow.Dist_Mitosis[t]/(double)M,Statsnow.Dist_Telophase[t]/(double)T);
            }
            fclose(f);
        }
    }
}

void run_stochastic_cell_cycle(const char InputNode[], double p_input,
                               const char InputNode2[], double p_input2,
                               longint_ARRAY_pair *Hits,
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                               longint_ARRAY *gl_keep_asis_A, int T_max_cc){
    longint_ARRAY *l;
    Cell_cycle_stats_in_run STnow;
    char fn[400],fn2[400];
    FILE *f;
    unsigned long int G1,R,G2,M,T;
    
    sprintf(fn2, "Hits");
    if((Hits != NULL) && (Hits->N>0)){
        for (unsigned long int l=1; l<=Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD_A->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
        }
    }
    else sprintf(fn2,"%s_WILDTYPE",fn2);
    
    unsigned long int id_input=PD_A->BRN->MDAT->get_node_ID(InputNode);
    unsigned long int id_input2=PD_A->BRN->MDAT->get_node_ID(InputNode2);
    l=get_Healthy_CellCycle_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            Cell_cycle_stats_in_run Statsnow;
            Statsnow.N_cycles=0;
            Statsnow.N_endoredupl_T=0;
            Statsnow.N_endoredupl_G2=0;
            Statsnow.N_endoredupl_Aneup=0;
            Statsnow.T_live=0;
            Statsnow.Time_G1=0;
            Statsnow.Time_R=0;
            Statsnow.Time_G2=0;
            Statsnow.Time_M=0;
            Statsnow.Time_Telo=0;
            G1=0;R=0;G2=0;M=0;T=0;
            for (int t=0;t<=MAX_ARREST; t++){
                Statsnow.Dist_G1[t]=0;
                Statsnow.Dist_G2[t]=0;
                Statsnow.Dist_Mitosis[t]=0;
                Statsnow.Dist_Replication[t]=0;
                Statsnow.Dist_Telophase[t]=0;
            }
            
            run_stochastic_Input_experiment_MOD(id_input,p_input,id_input2,p_input2,Hits,(int)l->A[j],1,PD_A,0);
            
            int t_live_total =0; int deaths =0;
            while (t_live_total<T_max_cc) {
                STnow=count_cell_cycles_stochastic_input(InputNode,p_input,InputNode2,p_input2,Hits,(int)l->A[j],1,PD_A,0,T_max_cc);
                t_live_total+=STnow.T_live;
                if (STnow.T_live<T_max_cc) {
                    deaths++;
                }
                
                Statsnow.N_cycles+=STnow.N_cycles;
                Statsnow.N_endoredupl_T+=STnow.N_endoredupl_T;
                Statsnow.N_endoredupl_G2+=STnow.N_endoredupl_G2;
                Statsnow.N_endoredupl_Aneup+=STnow.N_endoredupl_Aneup;
                for (int t=0;t<=MAX_ARREST; t++){
                    Statsnow.Dist_G1[t]+=STnow.Dist_G1[t];
                    if(STnow.Dist_G1[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_G2[t]+=STnow.Dist_G2[t];
                    if(STnow.Dist_G2[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_Mitosis[t]+=STnow.Dist_Mitosis[t];
                    if(STnow.Dist_Mitosis[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_Replication[t]+=STnow.Dist_Replication[t];
                    if(STnow.Dist_Replication[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_Telophase[t]+=STnow.Dist_Telophase[t];
                    if(STnow.Dist_Telophase[t]>0) Statsnow.MAX_INTERVAL=t;
                }
                Statsnow.T_live+=STnow.T_live;
                Statsnow.Time_G1+=STnow.Time_G1;
                Statsnow.Time_R+=STnow.Time_R;
                Statsnow.Time_G2+=STnow.Time_G2;
                Statsnow.Time_M+=STnow.Time_M;
                Statsnow.Time_Telo+=STnow.Time_Telo;
           }
            
            for (int t=0;t<=MAX_ARREST; t++){
                G1+=Statsnow.Dist_G1[t];
                R+=Statsnow.Dist_Replication[t];
                G2+=Statsnow.Dist_G2[t];
                M+=Statsnow.Dist_Mitosis[t];
                T+=Statsnow.Dist_Telophase[t];
            }
            sprintf(fn,"%s/%s/_EXP/%s/%s_Stochastic_%s-%.2lf_%s-%.2lf__A%d_As-%d__STATS_%d-CC_%d-Dead.txt",MODEL_Directory,PD_A->BRN->path,fn2,fn2,InputNode,p_input,InputNode2,p_input2,(int)l->A[j],1,Statsnow.N_cycles,deaths);
            // printf("%s\n",fn);
            printf("%lg\t%lg\t%d\t%d\t%lld\n",p_input,p_input2,Statsnow.N_cycles,deaths,Statsnow.T_live);
            
            f=fopen(fn,"w");
            for (int t=0;t<=Statsnow.MAX_INTERVAL; t++) {
                fprintf(f, "%d\t%lg\t%lg\t%lg\t%lg\t%lg\n",t,Statsnow.Dist_G1[t]/(double)G1,Statsnow.Dist_Replication[t]/(double)R,Statsnow.Dist_G2[t]/(double)G2,Statsnow.Dist_Mitosis[t]/(double)M,Statsnow.Dist_Telophase[t]/(double)T);
            }
            fclose(f);
        }
    }
}



void run_stochastic_cell_cycle_SWEEP(const char InputNode[], longint_ARRAY_pair *Hits,
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                               longint_ARRAY *gl_keep_asis_A, int T_max_cc){
    longint_ARRAY *l;
    Cell_cycle_stats_in_run STnow;
    char fn[400],fn2[400];
    FILE *f;
    
    unsigned long long int maxG1=0,maxR=0,maxG2=0,maxM=0,maxT=0;
    unsigned long int *G1,*R,*G2,*M,*T;
    
    sprintf(fn2, "Hits");
    if((Hits != NULL) && (Hits->N>0)){
        for (unsigned long int l=1; l<=Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD_A->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
        }
    }
    else sprintf(fn2,"%s_WILDTYPE",fn2);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD_A->BRN->path,fn2);
    system(fn);
    
    unsigned long int id_input=PD_A->BRN->MDAT->get_node_ID(InputNode);
    Cell_cycle_stats_in_run *Statsnow;
    
    l=get_Healthy_CellCycle_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            Statsnow=new Cell_cycle_stats_in_run[BOXN+2];
            G1=new unsigned long int[BOXN+2];
            R=new unsigned long int[BOXN+2];
            G2=new unsigned long int[BOXN+2];
            M=new unsigned long int[BOXN+2];
            T=new unsigned long int[BOXN+2];
            unsigned int box=0;
            for(double p_input=0;p_input<=1.; p_input+=1./(double)BOXN){
                box++;
                Statsnow[box].N_cycles=0;
                Statsnow[box].N_endoredupl_T=0;
                Statsnow[box].N_endoredupl_G2=0;
                Statsnow[box].N_endoredupl_Aneup=0;
                Statsnow[box].T_live=0;
                Statsnow[box].N_deaths=0;
                Statsnow[box].Time_G1=0;
                Statsnow[box].Time_R=0;
                Statsnow[box].Time_G2=0;
                Statsnow[box].Time_M=0;
                Statsnow[box].Time_Telo=0;
                G1[box]=0;R[box]=0;G2[box]=0; M[box]=0;T[box]=0;
                for (int t=0;t<=MAX_ARREST; t++){
                    Statsnow[box].Dist_G1[t]=0;
                    Statsnow[box].Dist_G2[t]=0;
                    Statsnow[box].Dist_Mitosis[t]=0;
                    Statsnow[box].Dist_Replication[t]=0;
                    Statsnow[box].Dist_Telophase[t]=0;
                }
                
                run_stochastic_Input_experiment_MOD(id_input,p_input,Hits,(int)l->A[j],1,PD_A,0);
                
                unsigned long long int t_live_total =0;
                while (t_live_total<T_max_cc) {
                    STnow=count_cell_cycles_stochastic_input(InputNode,p_input,Hits,(int)l->A[j],1,PD_A,0,T_max_cc);
                    
                    t_live_total+=STnow.T_live;
                    if (STnow.T_live<T_max_cc) {
                        STnow.N_deaths++;
                    }
                    Statsnow[box].N_cycles+=STnow.N_cycles;
                    Statsnow[box].N_endoredupl_T+=STnow.N_endoredupl_T;
                    Statsnow[box].N_endoredupl_G2+=STnow.N_endoredupl_G2;
                    Statsnow[box].N_endoredupl_Aneup+=STnow.N_endoredupl_Aneup;
                    Statsnow[box].T_live+=STnow.T_live;
                    Statsnow[box].N_deaths+=STnow.N_deaths;
                    Statsnow[box].Time_G1+=STnow.Time_G1;
                    Statsnow[box].Time_R+=STnow.Time_R;
                    Statsnow[box].Time_G2+=STnow.Time_G2;
                    Statsnow[box].Time_M+=STnow.Time_M;
                    Statsnow[box].Time_Telo+=STnow.Time_Telo;
                    
                    for (int t=0;t<=MAX_ARREST; t++){
                        Statsnow[box].Dist_G1[t]+=STnow.Dist_G1[t];
                        Statsnow[box].Dist_G2[t]+=STnow.Dist_G2[t];
                        Statsnow[box].Dist_Mitosis[t]+=STnow.Dist_Mitosis[t];
                        Statsnow[box].Dist_Replication[t]+=STnow.Dist_Replication[t];
                        Statsnow[box].Dist_Telophase[t]+=STnow.Dist_Telophase[t];
                        
                        if(STnow.Dist_G1[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_G2[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_Mitosis[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_Replication[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_Telophase[t]>0) Statsnow[box].MAX_INTERVAL=t;
                       
                        if((STnow.Dist_G1[t]>0)&&(t>maxG1)) maxG1=t;
                        if((STnow.Dist_Replication[t]>0)&&(t>maxR)) maxR=t;
                        if((STnow.Dist_G2[t]>0)&&(t>maxG2)) maxG2=t;
                        if((STnow.Dist_Mitosis[t]>0)&&(t>maxM)) maxM=t;
                        if((STnow.Dist_Telophase[t]>0)&&(t>maxT)) maxT=t;

                    }
                }
            }
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;
                for (int t=0;t<=MAX_ARREST; t++){
                    G1[box]+=Statsnow[box].Dist_G1[t];
                    R[box]+=Statsnow[box].Dist_Replication[t];
                    G2[box]+=Statsnow[box].Dist_G2[t];
                    M[box]+=Statsnow[box].Dist_Mitosis[t];
                    T[box]+=Statsnow[box].Dist_Telophase[t];
                }
            }
            sprintf(fn,"%s/%s/_EXP/%s/%s_Stoch-%s__A%d_As-%d__DISTRIB_G1-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,(int)l->A[j],1,BOXN);
            printf("\t\t%s\n",fn);
            f=fopen(fn,"w");
            box=0;
            for (int t=0;t<=maxG1; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_G1[t]/(double)G1[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_Stoch-%s__A%d_As-%d__DISTRIB_Replication-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxR; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_Replication[t]/(double)R[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_Stoch-%s__A%d_As-%d__DISTRIB_G2-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxG2; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_G2[t]/(double)G2[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_Stoch-%s__A%d_As-%d__DISTRIB_Metaphase-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxM; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_Mitosis[t]/(double)M[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_Stoch-%s__A%d_As-%d__DISTRIB_Telophase-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxT; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_Telophase[t]/(double)T[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            
            sprintf(fn,"%s/%s/_EXP/%s/%s_Stoch-%s__A%d_As-%d__FRACTIONS-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;
                fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        p_hits,
                        100*Statsnow[box].N_cycles/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_deaths/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_endoredupl_G2/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_endoredupl_Aneup/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_endoredupl_T/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_G1/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_R/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_G2/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_M/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_Telo/(double)Statsnow[box].T_live);
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s__A%d_As-%d__Time_Per_cycle-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;
                if (Statsnow[box].N_cycles>0)
                    fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                            p_hits,
                            Statsnow[box].Time_G1/(double)G1[box],
                            Statsnow[box].Time_R/(double)R[box],
                            Statsnow[box].Time_G2/(double)G2[box],
                            Statsnow[box].Time_M/(double)M[box],
                            Statsnow[box].Time_Telo/(double)T[box]);
            }
            fclose(f);
            delete[] Statsnow;Statsnow=NULL;
        }
     }
}

void run_cell_cycle_stochastic_hits(const char InputNode[], double p_hits, longint_ARRAY_pair *Hits,
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                               longint_ARRAY *gl_keep_asis_A, int T_max_cc){
    longint_ARRAY *l;
    Cell_cycle_stats_in_run STnow;
    char fn[400],fn2[400];
    FILE *f;
    unsigned long int G1,G2,R,T,M;
    
    sprintf(fn2, "Stochastic_Hits_%2.lg_",p_hits);
    if((Hits != NULL) && (Hits->N>0)){
        for (unsigned long int l=1; l<=Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD_A->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
        }
    }
    else sprintf(fn2,"%s_WILDTYPE",fn2);
    
  //  unsigned long int id_input=PD_A->BRN->MDAT->get_node_ID(InputNode);
    l=get_Healthy_CellCycle_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            Cell_cycle_stats_in_run Statsnow;
            Statsnow.N_cycles=0;
            Statsnow.N_endoredupl_T=0;
            Statsnow.N_endoredupl_G2=0;
            Statsnow.N_endoredupl_Aneup=0;
            Statsnow.T_live=0;
            Statsnow.Time_G1=0;
            Statsnow.Time_R=0;
            Statsnow.Time_G2=0;
            Statsnow.Time_M=0;
            Statsnow.Time_Telo=0;
            G1=0;R=0;G2=0;M=0;T=0;
            for (int t=0;t<=MAX_ARREST; t++){
                Statsnow.Dist_G1[t]=0;
                Statsnow.Dist_G2[t]=0;
                Statsnow.Dist_Mitosis[t]=0;
                Statsnow.Dist_Replication[t]=0;
                Statsnow.Dist_Telophase[t]=0;
            }
            
            //run_stochastic_HITS_experiment_MOD(id_input,p_hits,Hits,(int)l->A[j],1,PD_A,0);
            
            unsigned long long int t_live_total =0; int deaths =0;
            while (t_live_total<T_max_cc) {
                STnow=count_cell_cycles_stochastic_HITS(InputNode,p_hits,Hits,(int)l->A[j],1,PD_A,0,T_max_cc);
                t_live_total+=STnow.T_live;
                if (STnow.T_live<T_max_cc) {
                    deaths++;
                }
                
                Statsnow.N_cycles+=STnow.N_cycles;
                 Statsnow.N_endoredupl_T+=STnow.N_endoredupl_T;
                Statsnow.N_endoredupl_G2+=STnow.N_endoredupl_G2;
                Statsnow.N_endoredupl_Aneup+=STnow.N_endoredupl_Aneup;
                for (int t=0;t<=MAX_ARREST; t++){
                    Statsnow.Dist_G1[t]+=STnow.Dist_G1[t];
                    if(STnow.Dist_G1[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_G2[t]+=STnow.Dist_G2[t];
                    if(STnow.Dist_G2[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_Mitosis[t]+=STnow.Dist_Mitosis[t];
                    if(STnow.Dist_Mitosis[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_Replication[t]+=STnow.Dist_Replication[t]; if(STnow.Dist_Replication[t]>0) Statsnow.MAX_INTERVAL=t;
                    Statsnow.Dist_Telophase[t]+=STnow.Dist_Telophase[t];
                    if(STnow.Dist_Telophase[t]>0) Statsnow.MAX_INTERVAL=t;
                }
                Statsnow.T_live+=STnow.T_live;
                Statsnow.Time_G1+=STnow.Time_G1;
                Statsnow.Time_R+=STnow.Time_R;
                Statsnow.Time_G2+=STnow.Time_G2;
                Statsnow.Time_M+=STnow.Time_M;
                Statsnow.Time_Telo+=STnow.Time_Telo;
            }
            for (int t=0;t<=MAX_ARREST; t++){
                G1+=Statsnow.Dist_G1[t];
                R+=Statsnow.Dist_Replication[t];
                G2+=Statsnow.Dist_G2[t];
                M+=Statsnow.Dist_Mitosis[t];
                T+=Statsnow.Dist_Telophase[t];
            }
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s__A%d_As-%d__STATS_%d-CC_%d-Dead.txt",MODEL_Directory,PD_A->BRN->path,fn2,fn2,InputNode,(int)l->A[j],1,Statsnow.N_cycles,deaths);
            // printf("%s\n",fn);
            printf("%lg\t%d\t%d\t%lld\n",p_hits,Statsnow.N_cycles,deaths,Statsnow.T_live);
            
            f=fopen(fn,"w");
            for (int t=0;t<=Statsnow.MAX_INTERVAL; t++) {
                fprintf(f, "%d\t%lg\t%lg\t%lg\t%lg\t%lg\n",t,Statsnow.Dist_G1[t]/(double)G1,Statsnow.Dist_Replication[t]/(double)R,Statsnow.Dist_G2[t]/(double)G2,Statsnow.Dist_Mitosis[t]/(double)M,Statsnow.Dist_Telophase[t]/(double)T);

            }
            fclose(f);
        }
    }
}



void run_cell_cycle_stochastic_hits_SWEEP(const char InputNode[],
                                          double p_input,
                                          longint_ARRAY_pair *Hits,
                                          MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                          longint_ARRAY *gl_keep_asis_A,
                                          unsigned long long int T_max_cc,unsigned int G0_CC){
    longint_ARRAY *l;
    Cell_cycle_stats_in_run STnow;
    char fn[400],fn2[400];
    FILE *f;

    unsigned long long int maxG1=0,maxR=0,maxG2=0,maxM=0,maxT=0;
    unsigned long int *G1,*R,*G2,*M,*T;
    
    sprintf(fn2, "Stochastic_Hits_SWEEP_");
    if((Hits != NULL) && (Hits->N>0)){
        for (unsigned long int l=1; l<=Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD_A->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
        }
    }
    else sprintf(fn2,"%s_WILDTYPE",fn2);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD_A->BRN->path,fn2);
    system(fn);
    
  //  unsigned long int id_input=PD_A->BRN->MDAT->get_node_ID(InputNode);
   
    if(G0_CC==1) l=get_Healthy_CellCycle_attractors(PD_A);
    else    l=get_Healthy_G0_attractors(PD_A);

    Cell_cycle_stats_in_run *Statsnow;
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            Statsnow=new Cell_cycle_stats_in_run[BOXN+2];
            G1=new unsigned long int[BOXN+2];
            R=new unsigned long int[BOXN+2];
            G2=new unsigned long int[BOXN+2];
            M=new unsigned long int[BOXN+2];
            T=new unsigned long int[BOXN+2];
            unsigned int box=0;
            double *timetodeath;
            timetodeath=new double[BOXN+2];
            int deathhappens=0;
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;timetodeath[box]=0;
                Statsnow[box].N_cycles=0;
                Statsnow[box].N_endoredupl_T=0;
                Statsnow[box].N_endoredupl_G2=0;
                Statsnow[box].N_endoredupl_Aneup=0;
                Statsnow[box].T_live=0;
                Statsnow[box].N_deaths=0;
                Statsnow[box].Time_G1=0;
                Statsnow[box].Time_R=0;
                Statsnow[box].Time_G2=0;
                Statsnow[box].Time_M=0;
                Statsnow[box].Time_Telo=0;
                G1[box]=0;G2[box]=0;R[box]=0;T[box]=0;M[box]=0;
                for (int t=0;t<=MAX_ARREST; t++){
                    Statsnow[box].Dist_G1[t]=0;
                    Statsnow[box].Dist_G2[t]=0;
                    Statsnow[box].Dist_Mitosis[t]=0;
                    Statsnow[box].Dist_Replication[t]=0;
                    Statsnow[box].Dist_Telophase[t]=0;
                }

               // if(box<=4)
              //  run_stochastic_HITS_experiment_MOD(id_input,p_hits,Hits,(int)l->A[j],1,PD_A,0);

                unsigned long long int t_live_total =0;
                while (t_live_total<T_max_cc) {
                    STnow=count_Stochastic_cell_cycles_Stochastic_HITS(InputNode,p_input,p_hits,Hits,(int)l->A[j],(int)rand_int(1,PD_A->Coupled->Attractor_valleys[(int)l->A[j]]->N_a),PD_A,0,T_max_cc);
                    t_live_total+=STnow.T_live;
                    if (STnow.T_live<T_max_cc) {
                        STnow.N_deaths++;
                        timetodeath[box]+=STnow.T_live;
                        deathhappens=1;
                    }
                    Statsnow[box].N_cycles+=STnow.N_cycles;
                    Statsnow[box].N_endoredupl_T+=STnow.N_endoredupl_T;
                    Statsnow[box].N_endoredupl_G2+=STnow.N_endoredupl_G2;
                    Statsnow[box].N_endoredupl_Aneup+=STnow.N_endoredupl_Aneup;
                    Statsnow[box].T_live+=STnow.T_live;
                    Statsnow[box].N_deaths+=STnow.N_deaths;
                    Statsnow[box].Time_G1+=STnow.Time_G1;
                    Statsnow[box].Time_R+=STnow.Time_R;
                    Statsnow[box].Time_G2+=STnow.Time_G2;
                    Statsnow[box].Time_M+=STnow.Time_M;
                    Statsnow[box].Time_Telo+=STnow.Time_Telo;

                    for (int t=0;t<=MAX_ARREST; t++){
                        Statsnow[box].Dist_G1[t]+=STnow.Dist_G1[t];
                        Statsnow[box].Dist_G2[t]+=STnow.Dist_G2[t];
                        Statsnow[box].Dist_Mitosis[t]+=STnow.Dist_Mitosis[t];
                        Statsnow[box].Dist_Replication[t]+=STnow.Dist_Replication[t];
                        Statsnow[box].Dist_Telophase[t]+=STnow.Dist_Telophase[t];

                        if(STnow.Dist_G1[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_G2[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_Mitosis[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_Replication[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_Telophase[t]>0) Statsnow[box].MAX_INTERVAL=t;
                       
                        if((STnow.Dist_G1[t]>0)&&(t>maxG1)) maxG1=t;
                        if((STnow.Dist_Replication[t]>0)&&(t>maxR)) maxR=t;
                        if((STnow.Dist_G2[t]>0)&&(t>maxG2)) maxG2=t;
                        if((STnow.Dist_Mitosis[t]>0)&&(t>maxM)) maxM=t;
                        if((STnow.Dist_Telophase[t]>0)&&(t>maxT)) maxT=t;

                    }
                }
            }
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;
                for (int t=0;t<=MAX_ARREST; t++){
                    G1[box]+=Statsnow[box].Dist_G1[t];
                    R[box]+=Statsnow[box].Dist_Replication[t];
                    G2[box]+=Statsnow[box].Dist_G2[t];
                    M[box]+=Statsnow[box].Dist_Mitosis[t];
                    T[box]+=Statsnow[box].Dist_Telophase[t];
                }
            }
            if(deathhappens){
                sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__TimetoDeath-BOX_%d.txt",
                        MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                        InputNode,p_input,(int)l->A[j],1,BOXN);
                printf("\t\t%s\n",fn);
                f=fopen(fn,"w");
                unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;
                    if (Statsnow[box].N_deaths>0)
                        fprintf(f, "%lg\t%lg\n",p_hits,timetodeath[box]/(double)Statsnow[box].N_deaths);
                }
                fclose(f);
            }

            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_G1-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            printf("\t\t%s\n",fn);
            f=fopen(fn,"w");
            for (int t=0;t<=maxG1; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_G1[t]/(double)G1[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_Replication-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxR; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_Replication[t]/(double)R[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_G2-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxG2; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_G2[t]/(double)G2[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_Metaphase-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxM; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_Mitosis[t]/(double)M[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_Telophase-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxT; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_Telophase[t]/(double)T[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);

            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__FRACTIONS-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;
                fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        p_hits,
                        100*Statsnow[box].N_cycles/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_deaths/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_endoredupl_G2/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_endoredupl_Aneup/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_endoredupl_T/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_G1/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_R/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_G2/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_M/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_Telo/(double)Statsnow[box].T_live);
            }
            fclose(f);
            
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__Time_Per_cycle-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;
                if (Statsnow[box].N_cycles>0)
                    fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        p_hits,
                        Statsnow[box].Time_G1/(double)G1[box],
                        Statsnow[box].Time_R/(double)R[box],
                        Statsnow[box].Time_G2/(double)T[box],
                        Statsnow[box].Time_M/(double)M[box],
                        Statsnow[box].Time_Telo/(double)T[box]);
            }
            fclose(f);
            delete[] Statsnow;Statsnow=NULL;
        }
    }
}

void run_cell_cycle_stochastic_hits_SWEEP_Assync(const char InputNode[],
                                          double p_input,
                                          longint_ARRAY_pair *Hits,
                                          MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                          longint_ARRAY *gl_keep_asis_A,
                                          unsigned long long int T_max_cc,unsigned int G0_CC){
    longint_ARRAY *l;
    Cell_cycle_stats_in_run STnow;
    char fn[400],fn2[400];
    FILE *f;
    
    unsigned long long int maxG1=0,maxR=0,maxG2=0,maxM=0,maxT=0;
    unsigned long int *G1,*R,*G2,*M,*T;
    
    sprintf(fn2, "Assync_%d_Stochastic_Hits_SWEEP_",ASSYNC_BIAS);
    if((Hits != NULL) && (Hits->N>0)){
        for (unsigned long int l=1; l<=Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD_A->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
        }
    }
    else sprintf(fn2,"%s_Assync_%d_WILDTYPE",fn2,ASSYNC_BIAS);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD_A->BRN->path,fn2);
    system(fn);
    
  //  unsigned long int id_input=PD_A->BRN->MDAT->get_node_ID(InputNode);
    
    if(G0_CC==1) l=get_Healthy_CellCycle_attractors(PD_A);
    else    l=get_Healthy_G0_attractors(PD_A);
    
    Cell_cycle_stats_in_run *Statsnow;
    double *p_hitsl; p_hitsl= new double [Hits->N+1];
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            Statsnow=new Cell_cycle_stats_in_run[BOXN+2];
            G1=new unsigned long int[BOXN+2];
            R=new unsigned long int[BOXN+2];
            G2=new unsigned long int[BOXN+2];
            M=new unsigned long int[BOXN+2];
            T=new unsigned long int[BOXN+2];
            unsigned int box=0;
            double *timetodeath;
            timetodeath=new double[BOXN+2];
            int deathhappens=0;
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;timetodeath[box]=0;
                Statsnow[box].N_cycles=0;
                Statsnow[box].N_endoredupl_T=0;
                Statsnow[box].N_endoredupl_G2=0;
                Statsnow[box].N_endoredupl_Aneup=0;
                Statsnow[box].T_live=0;
                Statsnow[box].N_deaths=0;
                Statsnow[box].Time_G1=0;
                Statsnow[box].Time_R=0;
                Statsnow[box].Time_G2=0;
                Statsnow[box].Time_M=0;
                Statsnow[box].Time_Telo=0;
                G1[box]=0;G2[box]=0;R[box]=0;T[box]=0;M[box]=0;
                for (int t=0;t<=MAX_ARREST; t++){
                    Statsnow[box].Dist_G1[t]=0;
                    Statsnow[box].Dist_G2[t]=0;
                    Statsnow[box].Dist_Mitosis[t]=0;
                    Statsnow[box].Dist_Replication[t]=0;
                    Statsnow[box].Dist_Telophase[t]=0;
                }
                
                // if(box<=4)
                //  run_stochastic_HITS_experiment_MOD(id_input,p_hits,Hits,(int)l->A[j],1,PD_A,0);
                
                for(unsigned int iii=0;iii<=Hits->N;iii++) p_hitsl[iii]=p_hits;
                
                unsigned long long int t_live_total =0;
                while (t_live_total<T_max_cc) {
                    STnow=count_Stochastic_cell_cycles_Stochastic_HITS_Assync(InputNode,p_input,p_hitsl,Hits,(int)l->A[j],(int)rand_int(1,PD_A->Coupled->Attractor_valleys[(int)l->A[j]]->N_a),PD_A,0,T_max_cc);
                    t_live_total+=STnow.T_live;
                    if (STnow.T_live<T_max_cc) {
                        STnow.N_deaths++;
                        timetodeath[box]+=STnow.T_live;
                        deathhappens=1;
                    }
                    Statsnow[box].N_cycles+=STnow.N_cycles;
                    Statsnow[box].N_endoredupl_T+=STnow.N_endoredupl_T;
                    Statsnow[box].N_endoredupl_G2+=STnow.N_endoredupl_G2;
                    Statsnow[box].N_endoredupl_Aneup+=STnow.N_endoredupl_Aneup;
                    Statsnow[box].T_live+=STnow.T_live;
                    Statsnow[box].N_deaths+=STnow.N_deaths;
                    Statsnow[box].Time_G1+=STnow.Time_G1;
                    Statsnow[box].Time_R+=STnow.Time_R;
                    Statsnow[box].Time_G2+=STnow.Time_G2;
                    Statsnow[box].Time_M+=STnow.Time_M;
                    Statsnow[box].Time_Telo+=STnow.Time_Telo;
                    
                    for (int t=0;t<=MAX_ARREST; t++){
                        Statsnow[box].Dist_G1[t]+=STnow.Dist_G1[t];
                        Statsnow[box].Dist_G2[t]+=STnow.Dist_G2[t];
                        Statsnow[box].Dist_Mitosis[t]+=STnow.Dist_Mitosis[t];
                        Statsnow[box].Dist_Replication[t]+=STnow.Dist_Replication[t];
                        Statsnow[box].Dist_Telophase[t]+=STnow.Dist_Telophase[t];
                        
                        if(STnow.Dist_G1[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_G2[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_Mitosis[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_Replication[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        if(STnow.Dist_Telophase[t]>0) Statsnow[box].MAX_INTERVAL=t;
                        
                        if((STnow.Dist_G1[t]>0)&&(t>maxG1)) maxG1=t;
                        if((STnow.Dist_Replication[t]>0)&&(t>maxR)) maxR=t;
                        if((STnow.Dist_G2[t]>0)&&(t>maxG2)) maxG2=t;
                        if((STnow.Dist_Mitosis[t]>0)&&(t>maxM)) maxM=t;
                        if((STnow.Dist_Telophase[t]>0)&&(t>maxT)) maxT=t;
                        
                    }
                }
            }
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;
                for (int t=0;t<=MAX_ARREST; t++){
                    G1[box]+=Statsnow[box].Dist_G1[t];
                    R[box]+=Statsnow[box].Dist_Replication[t];
                    G2[box]+=Statsnow[box].Dist_G2[t];
                    M[box]+=Statsnow[box].Dist_Mitosis[t];
                    T[box]+=Statsnow[box].Dist_Telophase[t];
                }
            }
            if(deathhappens){
                sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__TimetoDeath-BOX_%d.txt",
                        MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                        InputNode,p_input,(int)l->A[j],1,BOXN);
                printf("\t\t%s\n",fn);
                f=fopen(fn,"w");
                unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;
                    if (Statsnow[box].N_deaths>0)
                        fprintf(f, "%lg\t%lg\n",p_hits,timetodeath[box]/(double)Statsnow[box].N_deaths);
                }
                fclose(f);
            }
            
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_G1-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            printf("\t\t%s\n",fn);
            f=fopen(fn,"w");
            for (int t=0;t<=maxG1; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_G1[t]/(double)G1[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_Replication-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxR; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_Replication[t]/(double)R[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_G2-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxG2; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_G2[t]/(double)G2[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_Metaphase-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxM; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_Mitosis[t]/(double)M[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__DISTRIB_Telophase-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            for (int t=0;t<=maxT; t++) {
                fprintf(f, "%d",t); unsigned int box=0;
                for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                    box++;  fprintf(f, "\t%lg",Statsnow[box].Dist_Telophase[t]/(double)T[box]);
                }
                fprintf(f,"\n");
            }
            fclose(f);
            
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__FRACTIONS-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;
                fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        p_hits,
                        100*Statsnow[box].N_cycles/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_deaths/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_endoredupl_G2/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_endoredupl_Aneup/(double)Statsnow[box].T_live,
                        100*Statsnow[box].N_endoredupl_T/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_G1/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_R/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_G2/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_M/(double)Statsnow[box].T_live,
                        100*Statsnow[box].Time_Telo/(double)Statsnow[box].T_live);
            }
            fclose(f);
            
            sprintf(fn,"%s/%s/_EXP/%s/%s_%s-%.2lf__A%d_As-%d__Time_Per_cycle-BOX_%d.txt",
                    MODEL_Directory,PD_A->BRN->path,fn2,fn2,
                    InputNode,p_input,(int)l->A[j],1,BOXN);
            f=fopen(fn,"w");
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOXN){
                box++;
                if (Statsnow[box].N_cycles>0)
                    fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                            p_hits,
                            Statsnow[box].Time_G1/(double)G1[box],
                            Statsnow[box].Time_R/(double)R[box],
                            Statsnow[box].Time_G2/(double)T[box],
                            Statsnow[box].Time_M/(double)M[box],
                            Statsnow[box].Time_Telo/(double)T[box]);
            }
            fclose(f);
            delete[] Statsnow;Statsnow=NULL;
        }
    }
}

void run_Apoptosis_Entry_Delayed_HITS(longint_ARRAY_pair *Hits,
                                 MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                 longint_ARRAY *gl_keep_asis_A,int entry_after_hit){
    longint_ARRAY *l;

    l=get_Healthy_G0_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_Delayed_Hits_Pulse_experiment_MOD(Hits,(int)l->A[j],1,PD_A,0,"Trail",50,entry_after_hit);
        }
    }
}

void run_cell_cycle_Entry_Delayed_HITS(longint_ARRAY_pair *Hits,
                                       MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                       longint_ARRAY *gl_keep_asis_A,int entry_after_hit, int GF_length){
    longint_ARRAY *l;
    
    l=get_Healthy_G0_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_Delayed_Hits_Pulse_experiment_MOD(Hits,(int)l->A[j],1,PD_A,0,"GF_High",GF_length,entry_after_hit);
        }
    }
}

void run_cell_cycle_Pause_Delayed_HITS(longint_ARRAY_pair *Hits,
                                       MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                       longint_ARRAY *gl_keep_asis_A,int entry_after_hit){
    longint_ARRAY *l;
    
    l=get_Healthy_CellCycle_attractors(PD_A);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_Delayed_Hits_Pulse_experiment_MOD(Hits,(int)l->A[j],1,PD_A,0,"Gamma",17,entry_after_hit);
            run_Delayed_Hits_Pulse_experiment_MOD(Hits,(int)l->A[j],1,PD_A,0,"UV",17,entry_after_hit);
        }
    }
}

void run_cell_cycle_Pulse(const char Pulsename[],int pulselength,
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                               longint_ARRAY *gl_keep_asis){
    longint_ARRAY *l;
    
    l=get_Healthy_CellCycle_attractors(PD);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_OneNode_pulse_experiment_MOD(Pulsename,(int)l->A[j],1,PD,0.000000001,pulselength);
        }
    }
}

void build_SUBSET_Model_and_sample_state_space(const char modelname_REF[],
                                               const char modules_name[],
                                               const char subset[],
                                               int subset_mode){
    Boolean_RegNet *A,*B;
    MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,*PD_B;
    
    longint_ARRAY *gl_keep_asis_A,*gl_keep_asis_B;
    
    A=new Boolean_RegNet(modelname_REF);
    gl_keep_asis_A=get_input_node_names(A);

    A->load_Module_IDs_and_Drawing_Order();

//    int m1=A->load_Module_Drawing_order(modules_name);
//    A->load_Module_IDs(modules_name,m1);

    PD_A=build_and_sample(A, gl_keep_asis_A, modules_name,0,A->name,0);
    
    
  //  run_OneNode_pulse_experiment_MOD("GF_High",23,1,PD_A,0.000000001,60);// exit(1);
  //  run_OneNode_pulse_experiment_MOD("GF_High",20,1,PD_A,0.000000001,60);// exit(1);
  //  run_OneNode_pulse_experiment_MOD("GF_High",25,1,PD_A,0.000000001,60);// exit(1);
 //   run_OneNode_pulse_experiment_MOD("GF_High",21,1,PD_A,0.000000001,60);// exit(1);
  //  run_OneNode_pulse_experiment_MOD("GF_High",33,1,PD_A,0.000000001,60);// exit(1);
//    run_OneNode_pulse_experiment_MOD("GF_High",16,1,PD_A,0.000000001,60);// exit(1);
//    run_OneNode_pulse_experiment_MOD("GF_High",17,1,PD_A,0.000000001,60);// exit(1);
 //   run_OneNode_pulse_experiment_MOD("Gamma",13,1,PD_A,0.000000001,60);// exit(1);
 //   run_OneNode_pulse_experiment_MOD("GF_High",25,1,PD_A,0.000000001,60);// exit(1);
  //     run_OneNode_pulse_experiment_MOD("GF",14,1,PD_A,0.000000001,60); exit(1);
////
 //     run_transition_from_cycle(16,"Gamma",PD_A,gl_keep_asis_A);
 //     run_transition_from_cycle(16,"GF_High",PD_A,gl_keep_asis_A);
  //    run_transition_from_cycle(13,"Trail",PD_A,gl_keep_asis_A);
    //     run_transition_from_cycle(13,"GF_High",PD_A,gl_keep_asis_A);
////

//    run_Model_is_Awesome_checklist(PD_A,gl_keep_asis_A);
    
    exit(1);
    
    B=generate_SUBSET_RNB(A,subset,subset_mode);
    gl_keep_asis_B=get_input_node_names(B);
    PD_B=build_and_sample(B, gl_keep_asis_B, modules_name,0,A->name,1);
    
    compare_attractor_graphs(PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
    Export_Attractor_Life_and_Death_table(PD_A,gl_keep_asis_A,1);
    Export_Attractor_Life_and_Death_table(PD_B,gl_keep_asis_B,1);
    COMPARE_Attractor_Life_and_Death_table(PD_A,gl_keep_asis_A,PD_B,gl_keep_asis_B);
    PD_A->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis_B->N)),p_error);
    PD_B->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis_B->N)),p_error);
  
    run_cell_cycle_entry_exit(PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B,3);
    
  //  run_OneNode_pulse_experiment_MOD("GF_High",15,1,PD_B,0.000000001,60);// exit(1);
  //  run_OneNode_pulse_experiment_MOD("Gamma",15,1,PD_B,0.000000001,60);// exit(1);
//    run_Model_is_Awesome_checklist(PD_B,gl_keep_asis_B);
    
}

void Generate_and_RUN_Phyton_script_BAR_GRAPHS(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char nodename[],longint_ARRAY *gl_keep_asis,int lockvalue){
    char fn[400];
    FILE *f;
    
    if(lockvalue==0) sprintf(fn,"%s/%s/_EXP/%s_BARGRAPHS_Knockdown.py",MODEL_Directory,PD->BRN->path,nodename);
    else  sprintf(fn,"%s/%s/_EXP/%s_BARGRAPHS_OverExpression.py",MODEL_Directory,PD->BRN->path,nodename);

    printf("%s\n",fn);//getchar();
    f=fopen(fn, "w");
    fprintf(f, "import matplotlib.pyplot as plt\nimport numpy as np\n");
    fprintf(f, "from mpl_toolkits import mplot3d\nfrom matplotlib.ticker import NullFormatter\n");
 
    fprintf(f, "Directory = '%s/%s/_EXP/'\n",MODEL_Directory,PD->BRN->path);
    fprintf(f, "molec = '%s'\n",nodename);
  
    longint_ARRAY *lcc,*lG0;
    lcc=get_Healthy_CellCycle_attractors(PD);
    lG0=get_Healthy_G0_attractors(PD);

    fprintf(f, "A_CC = '%ld'\n",lcc->A[1]);
    fprintf(f, "A_G0 = '%ld'\n",lG0->A[1]);
   
    fprintf(f, "p0='20'\np1='35'\np2='50'\np3='65'\np4='80'\np5='95'\n");
    
    // GF_High
    
    for (unsigned short int i=0;i<6; i++) {
        fprintf(f,"data = np.loadtxt('%s/%s/_EXP/Stochastic_Hits_SWEEP__%s-%d/Stochastic_Hits_SWEEP__%s-%d_GF_High-%.2lf__A%ld_As-1__FRACTIONS-BOX_10.txt')\n",MODEL_Directory,PD->BRN->path,nodename,lockvalue,nodename,lockvalue,(20+i*15)/100.,lcc->A[1]);
        fprintf(f,"P_lock=data[0:,0]\n");
        fprintf(f,"Normal_cc_%d_0=data[0:,1]\n",i);
        fprintf(f,"Apoptosis_%d_0=data[0:,2]\n",i);
        fprintf(f,"G2_4N_slip_%d_0=data[0:,3]\n",i);
        fprintf(f,"Aneupl_4N_%d_0=data[0:,4]\n",i);
        fprintf(f,"Telophase_4N_%d_0=data[0:,5]\n",i);
    }
    
    fprintf(f,"Normal_cc_row0 = (Normal_cc_0_0,Normal_cc_1_0,Normal_cc_2_0,Normal_cc_3_0,Normal_cc_4_0,Normal_cc_5_0)\n");
    fprintf(f,"Apoptosis_row0 = (Apoptosis_0_0,Apoptosis_1_0,Apoptosis_2_0,Apoptosis_3_0,Apoptosis_4_0,Apoptosis_5_0)\n");
    fprintf(f,"G2_4N_slip_row0 = (G2_4N_slip_0_0,G2_4N_slip_1_0,G2_4N_slip_2_0,G2_4N_slip_3_0,G2_4N_slip_4_0,G2_4N_slip_5_0)\n");
    fprintf(f,"Aneupl_4N_row0 = (Aneupl_4N_0_0,Aneupl_4N_1_0,Aneupl_4N_2_0,Aneupl_4N_3_0,Aneupl_4N_4_0,Aneupl_4N_5_0)\n");
    fprintf(f,"Telophase_4N_row0 = (Telophase_4N_0_0,Telophase_4N_1_0,Telophase_4N_2_0,Telophase_4N_3_0,Telophase_4N_4_0,Telophase_4N_5_0)\n");
    
    // Trail
    for (unsigned short int i=0;i<6; i++) {
        fprintf(f,"data = np.loadtxt('%s/%s/_EXP/Stochastic_Hits_SWEEP__%s-%d/Stochastic_Hits_SWEEP__%s-%d_Trail-%.2lf__A%ld_As-1__FRACTIONS-BOX_10.txt')\n",MODEL_Directory,PD->BRN->path,nodename,lockvalue,nodename,lockvalue,(20+i*15)/100.,lcc->A[1]);
        fprintf(f,"P_lock=data[0:,0]\n");
        fprintf(f,"Normal_cc_%d_1=data[0:,1]\n",i);
        fprintf(f,"Apoptosis_%d_1=data[0:,2]\n",i);
        fprintf(f,"G2_4N_slip_%d_1=data[0:,3]\n",i);
        fprintf(f,"Aneupl_4N_%d_1=data[0:,4]\n",i);
        fprintf(f,"Telophase_4N_%d_1=data[0:,5]\n",i);
    }
    
    fprintf(f,"Normal_cc_row1 = (Normal_cc_0_1,Normal_cc_1_1,Normal_cc_2_1,Normal_cc_3_1,Normal_cc_4_1,Normal_cc_5_1)\n");
    fprintf(f,"Apoptosis_row1 = (Apoptosis_0_1,Apoptosis_1_1,Apoptosis_2_1,Apoptosis_3_1,Apoptosis_4_1,Apoptosis_5_1)\n");
    fprintf(f,"G2_4N_slip_row1 = (G2_4N_slip_0_1,G2_4N_slip_1_1,G2_4N_slip_2_1,G2_4N_slip_3_1,G2_4N_slip_4_1,G2_4N_slip_5_1)\n");
    fprintf(f,"Aneupl_4N_row1 = (Aneupl_4N_0_1,Aneupl_4N_1_1,Aneupl_4N_2_1,Aneupl_4N_3_1,Aneupl_4N_4_1,Aneupl_4N_5_1)\n");
    fprintf(f,"Telophase_4N_row1 = (Telophase_4N_0_1,Telophase_4N_1_1,Telophase_4N_2_1,Telophase_4N_3_1,Telophase_4N_4_1,Telophase_4N_5_1)\n");
    
    
    // Trail, no GF_High
    for (unsigned short int i=0;i<6; i++) {
        fprintf(f,"data = np.loadtxt('%s/%s/_EXP/Stochastic_Hits_SWEEP__%s-%d/Stochastic_Hits_SWEEP__%s-%d_Trail-%.2lf__A%ld_As-1__FRACTIONS-BOX_10.txt')\n",MODEL_Directory,PD->BRN->path,nodename,lockvalue,nodename,lockvalue,(20+i*15)/100.,lG0->A[1]);
        fprintf(f,"P_lock=data[0:,0]\n");
        fprintf(f,"Normal_cc_%d_2=data[0:,1]\n",i);
        fprintf(f,"Apoptosis_%d_2=data[0:,2]\n",i);
        fprintf(f,"G2_4N_slip_%d_2=data[0:,3]\n",i);
        fprintf(f,"Aneupl_4N_%d_2=data[0:,4]\n",i);
        fprintf(f,"Telophase_4N_%d_2=data[0:,5]\n",i);
    }
    
    fprintf(f,"Normal_cc_row2 = (Normal_cc_0_2,Normal_cc_1_2,Normal_cc_2_2,Normal_cc_3_2,Normal_cc_4_2,Normal_cc_5_2)\n");
    fprintf(f,"Apoptosis_row2 = (Apoptosis_0_2,Apoptosis_1_2,Apoptosis_2_2,Apoptosis_3_2,Apoptosis_4_2,Apoptosis_5_2)\n");
    fprintf(f,"G2_4N_slip_row2 = (G2_4N_slip_0_2,G2_4N_slip_1_2,G2_4N_slip_2_2,G2_4N_slip_3_2,G2_4N_slip_4_2,G2_4N_slip_5_2)\n");
    fprintf(f,"Aneupl_4N_row2 = (Aneupl_4N_0_2,Aneupl_4N_1_2,Aneupl_4N_2_2,Aneupl_4N_3_2,Aneupl_4N_4_2,Aneupl_4N_5_2)\n");
    fprintf(f,"Telophase_4N_row2 = (Telophase_4N_0_2,Telophase_4N_1_2,Telophase_4N_2_2,Telophase_4N_3_2,Telophase_4N_4_2,Telophase_4N_5_2)\n");
    
    // GF, no GF_High
    for (unsigned short int i=0;i<6; i++) {
        fprintf(f,"data = np.loadtxt('%s/%s/_EXP/Stochastic_Hits_SWEEP__%s-%d/Stochastic_Hits_SWEEP__%s-%d_GF-%.2lf__A%ld_As-1__FRACTIONS-BOX_10.txt')\n",MODEL_Directory,PD->BRN->path,nodename,lockvalue,nodename,lockvalue,(20+i*15)/100.,lG0->A[1]);
        fprintf(f,"P_lock=data[0:,0]\n");
        fprintf(f,"Normal_cc_%d_3=data[0:,1]\n",i);
        fprintf(f,"Apoptosis_%d_3=data[0:,2]\n",i);
        fprintf(f,"G2_4N_slip_%d_3=data[0:,3]\n",i);
        fprintf(f,"Aneupl_4N_%d_3=data[0:,4]\n",i);
        fprintf(f,"Telophase_4N_%d_3=data[0:,5]\n",i);
    }
    
    fprintf(f,"Normal_cc_row3 = (Normal_cc_0_3,Normal_cc_1_3,Normal_cc_2_3,Normal_cc_3_3,Normal_cc_4_3,Normal_cc_5_3)\n");
    fprintf(f,"Apoptosis_row3 = (Apoptosis_0_3,Apoptosis_1_3,Apoptosis_2_3,Apoptosis_3_3,Apoptosis_4_3,Apoptosis_5_3)\n");
    fprintf(f,"G2_4N_slip_row3 = (G2_4N_slip_0_3,G2_4N_slip_1_3,G2_4N_slip_2_3,G2_4N_slip_3_3,G2_4N_slip_4_3,G2_4N_slip_5_3)\n");
    fprintf(f,"Aneupl_4N_row3 = (Aneupl_4N_0_3,Aneupl_4N_1_3,Aneupl_4N_2_3,Aneupl_4N_3_3,Aneupl_4N_4_3,Aneupl_4N_5_3)\n");
    fprintf(f,"Telophase_4N_row3 = (Telophase_4N_0_3,Telophase_4N_1_3,Telophase_4N_2_3,Telophase_4N_3_3,Telophase_4N_4_3,Telophase_4N_5_3)\n");
    
    
    fprintf(f,"width = 0.05\n");
    fprintf(f,"fig, axes = plt.subplots(4, 6, sharex='all')\n");
    fprintf(f,"for x in range(6):\n");
    for (unsigned short int i=0; i<=3;i++) {
        fprintf(f,"\taxes[%d, x].bar(P_lock, Normal_cc_row%d[x], width)\n",i,i);
        fprintf(f,"\taxes[%d, x].bar(P_lock, G2_4N_slip_row%d[x], width, bottom=Normal_cc_row%d[x])\n",i,i,i);
        fprintf(f,"\taxes[%d, x].bar(P_lock, Aneupl_4N_row%d[x], width, bottom=G2_4N_slip_row%d[x] + Normal_cc_row%d[x])\n",i,i,i,i);
        fprintf(f,"\taxes[%d, x].bar(P_lock, Telophase_4N_row%d[x], width, bottom=Aneupl_4N_row%d[x] + G2_4N_slip_row%d[x] + Normal_cc_row%d[x])\n",i,i,i,i,i);
        fprintf(f,"\taxes[%d, x].bar(P_lock, Apoptosis_row%d[x], width, bottom=Telophase_4N_row%d[x] + Aneupl_4N_row%d[x] + G2_4N_slip_row%d[x] + Normal_cc_row%d[x])\n",i,i,i,i,i,i);
    }
    
    for (unsigned short int i=0; i<6;i++)
        fprintf(f,"axes[0, %d].set_title('0.' + p%d + ' %%')\n",i,i);
    
    fprintf(f,"fig.gca().yaxis.set_minor_formatter(NullFormatter())\n");
    fprintf(f,"fig.subplots_adjust(top=0.92, bottom=0.1, left=0.13, right=0.95, hspace=0.25,wspace=0.45)\n");
    
    if(lockvalue==0) fprintf(f,"axes[3, 2].set_xlabel('%%  %s knockdown')\n",nodename);
    else  fprintf(f,"axes[3, 2].set_xlabel('%%  %s overactivation')\n",nodename);
   
    fprintf(f,"axes[0, 0].set_ylabel('GFH (Trail=0)',fontsize=9)\n");
    fprintf(f,"axes[1, 0].set_ylabel('Trail (GFH=1)',fontsize=9)\n");
    fprintf(f,"axes[2, 0].set_ylabel('Trail (GFH=0)',fontsize=9)\n");
    fprintf(f,"axes[3, 0].set_ylabel('GF (GFH=0)',fontsize=9)\n");
    if(lockvalue==0) fprintf(f,"fig.savefig('%s/%s/_EXP/%s_Knockdown_RESULTS.pdf')\n",MODEL_Directory,PD->BRN->path,nodename);
    else fprintf(f,"fig.savefig('%s/%s/_EXP/%s_Forced_activation_RESULTS.pdf')\n",MODEL_Directory,PD->BRN->path,nodename);
    
    //// Wild type
    
    for (unsigned short int i=0; i<=3;i++){
        fprintf(f,"WT_Normal_cc_row%d = (Normal_cc_0_%d[0],Normal_cc_1_%d[0],Normal_cc_2_%d[0],Normal_cc_3_%d[0],Normal_cc_4_%d[0],Normal_cc_5_%d[0])\n",i,i,i,i,i,i,i);
        fprintf(f,"WT_Apoptosis_row%d = (Apoptosis_0_%d[0],Apoptosis_1_%d[0],Apoptosis_2_%d[0],Apoptosis_3_%d[0],Apoptosis_4_%d[0],Apoptosis_5_%d[0])\n",i,i,i,i,i,i,i);
        fprintf(f,"WT_G2_4N_slip_row%d = (G2_4N_slip_0_%d[0],G2_4N_slip_1_%d[0],G2_4N_slip_2_%d[0],G2_4N_slip_3_%d[0],G2_4N_slip_4_%d[0],G2_4N_slip_5_%d[0])\n",i,i,i,i,i,i,i);
        fprintf(f,"WT_Aneupl_4N_row%d = (Aneupl_4N_0_%d[0],Aneupl_4N_1_%d[0],Aneupl_4N_2_%d[0],Aneupl_4N_3_%d[0],Aneupl_4N_4_%d[0],Aneupl_4N_5_%d[0])\n",i,i,i,i,i,i,i);
        fprintf(f,"WT_Telophase_4N_row%d = (Telophase_4N_0_%d[0],Telophase_4N_1_%d[0],Telophase_4N_2_%d[0],Telophase_4N_3_%d[0],Telophase_4N_4_%d[0],Telophase_4N_5_%d[0])\n",i,i,i,i,i,i,i);
    }
    
    fprintf(f,"Percent_signal = (");
    for (unsigned short int i=0; i<5;i++) fprintf(f,"%.2lf, ",(5+i*15)/100.);
    fprintf(f,"%.2lf)\n",(5+5*15)/100.);
    
    
    fprintf(f,"fig, axes = plt.subplots(1, 4, sharey='all')\n");
    for (unsigned short int i=0; i<=3;i++) {
        fprintf(f,"bottom1 = WT_Normal_cc_row%d\n",i);
        
        fprintf(f,"bottom2 = (Normal_cc_0_%d[0] + G2_4N_slip_0_%d[0], Normal_cc_1_%d[0] + G2_4N_slip_1_%d[0], Normal_cc_2_%d[0]+ G2_4N_slip_2_%d[0], Normal_cc_3_%d[0] + G2_4N_slip_3_%d[0], Normal_cc_4_%d[0]+ G2_4N_slip_4_%d[0], Normal_cc_5_%d[0]+ G2_4N_slip_5_%d[0])\n",i,i,i,i,i,i,i,i,i,i,i,i);
       
        fprintf(f,"bottom3 = (Normal_cc_0_%d[0] + G2_4N_slip_0_%d[0] + Aneupl_4N_0_%d[0], Normal_cc_1_%d[0] + G2_4N_slip_1_%d[0]+ Aneupl_4N_1_%d[0], Normal_cc_2_%d[0]+ G2_4N_slip_2_%d[0]+ Aneupl_4N_2_%d[0], Normal_cc_3_%d[0] + G2_4N_slip_3_%d[0] + Aneupl_4N_3_%d[0], Normal_cc_4_%d[0]+ G2_4N_slip_4_%d[0] + Aneupl_4N_4_%d[0], Normal_cc_5_%d[0] + G2_4N_slip_5_%d[0] + Aneupl_4N_5_%d[0])\n",i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i);
        
        fprintf(f,"bottom4 = (Normal_cc_0_%d[0] + G2_4N_slip_0_%d[0] + Aneupl_4N_0_%d[0] + Telophase_4N_0_%d[0], Normal_cc_1_%d[0] + G2_4N_slip_1_%d[0]+ Aneupl_4N_1_%d[0] + Telophase_4N_1_%d[0], Normal_cc_2_%d[0]+ G2_4N_slip_2_%d[0] + Aneupl_4N_2_%d[0] + Telophase_4N_2_%d[0], Normal_cc_3_%d[0] + G2_4N_slip_3_%d[0] + Aneupl_4N_3_%d[0] + Telophase_4N_3_%d[0], Normal_cc_4_%d[0]+ G2_4N_slip_4_%d[0] + Aneupl_4N_4_%d[0] + Telophase_4N_4_%d[0], Normal_cc_5_%d[0] + G2_4N_slip_5_%d[0] + Aneupl_4N_5_%d[0] + Telophase_4N_5_%d[0])\n",i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i);
        
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Normal_cc_row%d, width)\n",i,i);
        fprintf(f,"axes[%d].bar(Percent_signal, WT_G2_4N_slip_row%d, width, bottom=bottom1)\n",i,i);
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Aneupl_4N_row%d, width, bottom=bottom2)\n",i,i);
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Telophase_4N_row%d, width, bottom=bottom3)\n",i,i);
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Apoptosis_row%d, width, bottom=bottom4)\n",i,i);
    }
    
    fprintf(f,"fig.gca().yaxis.set_minor_formatter(NullFormatter())\n");
    fprintf(f,"fig.subplots_adjust(top=0.92, bottom=0.1, left=0.13, right=0.95, hspace=0.25,wspace=0.45)\n");
    
    fprintf(f,"axes[0].set_xlabel('%% GF_H (Trail=0)',fontsize=9)\n");
    fprintf(f,"axes[1].set_xlabel('%% Trail (GFH=1)',fontsize=9)\n");
    fprintf(f,"axes[2].set_xlabel('%% Trail (GFH=0)',fontsize=9)\n");
    fprintf(f,"axes[3].set_xlabel('%% GF (GFH=0)',fontsize=9)\n");
    
    fprintf(f,"fig.savefig('%s/%s/_EXP/WT_RESULTS.pdf')\n",MODEL_Directory,PD->BRN->path);

    fclose(f);
    
    if(lockvalue==0) sprintf(fn,"python %s/%s/_EXP/%s_BARGRAPHS_Knockdown.py\n",MODEL_Directory_system,PD->BRN->path,nodename);
    else  sprintf(fn,"python %s/%s/_EXP/%s_BARGRAPHS_OverExpression.py\n",MODEL_Directory_system,PD->BRN->path,nodename);
    system(fn);
}

void Generate_and_RUN_Phyton_script_BAR_GRAPHS_Assync(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char nodename[],longint_ARRAY *gl_keep_asis,int lockvalue){
    char fn[400];
    FILE *f;
    sprintf(fn,"%s/%s/_EXP/%s_Assync_%d_BARGRAPHS.py",MODEL_Directory,PD->BRN->path,nodename,ASSYNC_BIAS);
    printf("%s\n",fn);//getchar();
    f=fopen(fn, "w");
    fprintf(f, "import matplotlib.pyplot as plt\nimport numpy as np\n");
    fprintf(f, "from mpl_toolkits import mplot3d\nfrom matplotlib.ticker import NullFormatter\n");
    
    fprintf(f, "Directory = '%s/%s/_EXP/'\n",MODEL_Directory,PD->BRN->path);
    fprintf(f, "molec = '%s'\n",nodename);
    
    longint_ARRAY *lcc,*lG0;
    lcc=get_Healthy_CellCycle_attractors(PD);
    lG0=get_Healthy_G0_attractors(PD);
    
    fprintf(f, "A_CC = '%ld'\n",lcc->A[1]);
    fprintf(f, "A_G0 = '%ld'\n",lG0->A[1]);
    
    fprintf(f, "p0='20'\np1='35'\np2='50'\np3='65'\np4='80'\np5='95'\n");
    for (unsigned short int i=0;i<6; i++) {
        fprintf(f,"data = np.loadtxt('%s/%s/_EXP/Assync_%d_Stochastic_Hits_SWEEP__%s-%d/Assync_%d_Stochastic_Hits_SWEEP__%s-%d_GF_High-%.2lf__A%ld_As-1__FRACTIONS-BOX_10.txt')\n",MODEL_Directory,PD->BRN->path,ASSYNC_BIAS,nodename,lockvalue,ASSYNC_BIAS,nodename,lockvalue,(20+i*15)/100.,lcc->A[1]);
        fprintf(f,"P_lock=data[0:,0]\n");
        fprintf(f,"Normal_cc_%d_0=data[0:,1]\n",i);
        fprintf(f,"Apoptosis_%d_0=data[0:,2]\n",i);
        fprintf(f,"G2_4N_slip_%d_0=data[0:,3]\n",i);
        fprintf(f,"Aneupl_4N_%d_0=data[0:,4]\n",i);
        fprintf(f,"Telophase_4N_%d_0=data[0:,5]\n",i);
    }
    
    fprintf(f,"Normal_cc_row0 = (Normal_cc_0_0,Normal_cc_1_0,Normal_cc_2_0,Normal_cc_3_0,Normal_cc_4_0,Normal_cc_5_0)\n");
    fprintf(f,"Apoptosis_row0 = (Apoptosis_0_0,Apoptosis_1_0,Apoptosis_2_0,Apoptosis_3_0,Apoptosis_4_0,Apoptosis_5_0)\n");
    fprintf(f,"G2_4N_slip_row0 = (G2_4N_slip_0_0,G2_4N_slip_1_0,G2_4N_slip_2_0,G2_4N_slip_3_0,G2_4N_slip_4_0,G2_4N_slip_5_0)\n");
    fprintf(f,"Aneupl_4N_row0 = (Aneupl_4N_0_0,Aneupl_4N_1_0,Aneupl_4N_2_0,Aneupl_4N_3_0,Aneupl_4N_4_0,Aneupl_4N_5_0)\n");
    fprintf(f,"Telophase_4N_row0 = (Telophase_4N_0_0,Telophase_4N_1_0,Telophase_4N_2_0,Telophase_4N_3_0,Telophase_4N_4_0,Telophase_4N_5_0)\n");

    // Trail
    for (unsigned short int i=0;i<6; i++) {
        fprintf(f,"data = np.loadtxt('%s/%s/_EXP/Assync_%d_Stochastic_Hits_SWEEP__%s-%d/Assync_%d_Stochastic_Hits_SWEEP__%s-%d_Trail-%.2lf__A%ld_As-1__FRACTIONS-BOX_10.txt')\n",MODEL_Directory,PD->BRN->path,ASSYNC_BIAS,nodename,lockvalue,ASSYNC_BIAS,nodename,lockvalue,(20+i*15)/100.,lcc->A[1]);
        fprintf(f,"P_lock=data[0:,0]\n");
        fprintf(f,"Normal_cc_%d_1=data[0:,1]\n",i);
        fprintf(f,"Apoptosis_%d_1=data[0:,2]\n",i);
        fprintf(f,"G2_4N_slip_%d_1=data[0:,3]\n",i);
        fprintf(f,"Aneupl_4N_%d_1=data[0:,4]\n",i);
        fprintf(f,"Telophase_4N_%d_1=data[0:,5]\n",i);
    }
    
    fprintf(f,"Normal_cc_row1 = (Normal_cc_0_1,Normal_cc_1_1,Normal_cc_2_1,Normal_cc_3_1,Normal_cc_4_1,Normal_cc_5_1)\n");
    fprintf(f,"Apoptosis_row1 = (Apoptosis_0_1,Apoptosis_1_1,Apoptosis_2_1,Apoptosis_3_1,Apoptosis_4_1,Apoptosis_5_1)\n");
    fprintf(f,"G2_4N_slip_row1 = (G2_4N_slip_0_1,G2_4N_slip_1_1,G2_4N_slip_2_1,G2_4N_slip_3_1,G2_4N_slip_4_1,G2_4N_slip_5_1)\n");
    fprintf(f,"Aneupl_4N_row1 = (Aneupl_4N_0_1,Aneupl_4N_1_1,Aneupl_4N_2_1,Aneupl_4N_3_1,Aneupl_4N_4_1,Aneupl_4N_5_1)\n");
    fprintf(f,"Telophase_4N_row1 = (Telophase_4N_0_1,Telophase_4N_1_1,Telophase_4N_2_1,Telophase_4N_3_1,Telophase_4N_4_1,Telophase_4N_5_1)\n");

    
    // Trail, no GF_High
    for (unsigned short int i=0;i<6; i++) {
        fprintf(f,"data = np.loadtxt('%s/%s/_EXP/Assync_%d_Stochastic_Hits_SWEEP__%s-%d/Assync_%d_Stochastic_Hits_SWEEP__%s-%d_Trail-%.2lf__A%ld_As-1__FRACTIONS-BOX_10.txt')\n",MODEL_Directory,PD->BRN->path,ASSYNC_BIAS,nodename,lockvalue,ASSYNC_BIAS,nodename,lockvalue,(20+i*15)/100.,lG0->A[1]);
        fprintf(f,"P_lock=data[0:,0]\n");
        fprintf(f,"Normal_cc_%d_2=data[0:,1]\n",i);
        fprintf(f,"Apoptosis_%d_2=data[0:,2]\n",i);
        fprintf(f,"G2_4N_slip_%d_2=data[0:,3]\n",i);
        fprintf(f,"Aneupl_4N_%d_2=data[0:,4]\n",i);
        fprintf(f,"Telophase_4N_%d_2=data[0:,5]\n",i);
    }
    
    fprintf(f,"Normal_cc_row2 = (Normal_cc_0_2,Normal_cc_1_2,Normal_cc_2_2,Normal_cc_3_2,Normal_cc_4_2,Normal_cc_5_2)\n");
    fprintf(f,"Apoptosis_row2 = (Apoptosis_0_2,Apoptosis_1_2,Apoptosis_2_2,Apoptosis_3_2,Apoptosis_4_2,Apoptosis_5_2)\n");
    fprintf(f,"G2_4N_slip_row2 = (G2_4N_slip_0_2,G2_4N_slip_1_2,G2_4N_slip_2_2,G2_4N_slip_3_2,G2_4N_slip_4_2,G2_4N_slip_5_2)\n");
    fprintf(f,"Aneupl_4N_row2 = (Aneupl_4N_0_2,Aneupl_4N_1_2,Aneupl_4N_2_2,Aneupl_4N_3_2,Aneupl_4N_4_2,Aneupl_4N_5_2)\n");
    fprintf(f,"Telophase_4N_row2 = (Telophase_4N_0_2,Telophase_4N_1_2,Telophase_4N_2_2,Telophase_4N_3_2,Telophase_4N_4_2,Telophase_4N_5_2)\n");

    // GF, no GF_High
    for (unsigned short int i=0;i<6; i++) {
        fprintf(f,"data = np.loadtxt('%s/%s/_EXP/Assync_%d_Stochastic_Hits_SWEEP__%s-%d/Assync_%d_Stochastic_Hits_SWEEP__%s-%d_GF-%.2lf__A%ld_As-1__FRACTIONS-BOX_10.txt')\n",MODEL_Directory,PD->BRN->path,ASSYNC_BIAS,nodename,lockvalue,ASSYNC_BIAS,nodename,lockvalue,(20+i*15)/100.,lG0->A[1]);
        fprintf(f,"P_lock=data[0:,0]\n");
        fprintf(f,"Normal_cc_%d_3=data[0:,1]\n",i);
        fprintf(f,"Apoptosis_%d_3=data[0:,2]\n",i);
        fprintf(f,"G2_4N_slip_%d_3=data[0:,3]\n",i);
        fprintf(f,"Aneupl_4N_%d_3=data[0:,4]\n",i);
        fprintf(f,"Telophase_4N_%d_3=data[0:,5]\n",i);
    }
    
    fprintf(f,"Normal_cc_row3 = (Normal_cc_0_3,Normal_cc_1_3,Normal_cc_2_3,Normal_cc_3_3,Normal_cc_4_3,Normal_cc_5_3)\n");
    fprintf(f,"Apoptosis_row3 = (Apoptosis_0_3,Apoptosis_1_3,Apoptosis_2_3,Apoptosis_3_3,Apoptosis_4_3,Apoptosis_5_3)\n");
    fprintf(f,"G2_4N_slip_row3 = (G2_4N_slip_0_3,G2_4N_slip_1_3,G2_4N_slip_2_3,G2_4N_slip_3_3,G2_4N_slip_4_3,G2_4N_slip_5_3)\n");
    fprintf(f,"Aneupl_4N_row3 = (Aneupl_4N_0_3,Aneupl_4N_1_3,Aneupl_4N_2_3,Aneupl_4N_3_3,Aneupl_4N_4_3,Aneupl_4N_5_3)\n");
    fprintf(f,"Telophase_4N_row3 = (Telophase_4N_0_3,Telophase_4N_1_3,Telophase_4N_2_3,Telophase_4N_3_3,Telophase_4N_4_3,Telophase_4N_5_3)\n");

    
    fprintf(f,"width = 0.05\n");
    fprintf(f,"fig, axes = plt.subplots(4, 6, sharex='all')\n");
    fprintf(f,"for x in range(6):\n");
    for (unsigned short int i=0; i<=3;i++) {
        fprintf(f,"\taxes[%d, x].bar(P_lock, Normal_cc_row%d[x], width)\n",i,i);
        fprintf(f,"\taxes[%d, x].bar(P_lock, G2_4N_slip_row%d[x], width, bottom=Normal_cc_row%d[x])\n",i,i,i);
        fprintf(f,"\taxes[%d, x].bar(P_lock, Aneupl_4N_row%d[x], width, bottom=G2_4N_slip_row%d[x] + Normal_cc_row%d[x])\n",i,i,i,i);
        fprintf(f,"\taxes[%d, x].bar(P_lock, Telophase_4N_row%d[x], width, bottom=Aneupl_4N_row%d[x] + G2_4N_slip_row%d[x] + Normal_cc_row%d[x])\n",i,i,i,i,i);
        fprintf(f,"\taxes[%d, x].bar(P_lock, Apoptosis_row%d[x], width, bottom=Telophase_4N_row%d[x] + Aneupl_4N_row%d[x] + G2_4N_slip_row%d[x] + Normal_cc_row%d[x])\n",i,i,i,i,i,i);
    }
    for (unsigned short int i=0; i<6;i++)
        fprintf(f,"axes[0, %d].set_title('0.' + p%d + ' %%')\n",i,i);
    
    fprintf(f,"fig.gca().yaxis.set_minor_formatter(NullFormatter())\n");
    fprintf(f,"fig.subplots_adjust(top=0.92, bottom=0.1, left=0.13, right=0.95, hspace=0.25,wspace=0.45)\n");
    
    if(lockvalue==0) fprintf(f,"axes[3, 2].set_xlabel('%%  %s knockdown')\n",nodename);
    else  fprintf(f,"axes[3, 2].set_xlabel('%%  %s overactivation')\n",nodename);
    
    fprintf(f,"axes[0, 0].set_ylabel('GFH (Trail=0)',fontsize=9)\n");
    fprintf(f,"axes[1, 0].set_ylabel('Trail (GFH=1)',fontsize=9)\n");
    fprintf(f,"axes[2, 0].set_ylabel('Trail (GFH=0)',fontsize=9)\n");
    fprintf(f,"axes[3, 0].set_ylabel('GF (GFH=0)',fontsize=9)\n");
    if(lockvalue==0) fprintf(f,"fig.savefig('%s/%s/_EXP/%s_Knockdown_RESULTS_Assync_%d.pdf')\n",MODEL_Directory,PD->BRN->path,nodename,ASSYNC_BIAS);
    else fprintf(f,"fig.savefig('%s/%s/_EXP/%s_Forced_activation_RESULTS_Assync_%d.pdf')\n",MODEL_Directory,PD->BRN->path,nodename,ASSYNC_BIAS);
    
    //// Wild type
    
    for (unsigned short int i=0; i<=3;i++){
        fprintf(f,"WT_Normal_cc_row%d = (Normal_cc_0_%d[0],Normal_cc_1_%d[0],Normal_cc_2_%d[0],Normal_cc_3_%d[0],Normal_cc_4_%d[0],Normal_cc_5_%d[0])\n",i,i,i,i,i,i,i);
        fprintf(f,"WT_Apoptosis_row%d = (Apoptosis_0_%d[0],Apoptosis_1_%d[0],Apoptosis_2_%d[0],Apoptosis_3_%d[0],Apoptosis_4_%d[0],Apoptosis_5_%d[0])\n",i,i,i,i,i,i,i);
        fprintf(f,"WT_G2_4N_slip_row%d = (G2_4N_slip_0_%d[0],G2_4N_slip_1_%d[0],G2_4N_slip_2_%d[0],G2_4N_slip_3_%d[0],G2_4N_slip_4_%d[0],G2_4N_slip_5_%d[0])\n",i,i,i,i,i,i,i);
        fprintf(f,"WT_Aneupl_4N_row%d = (Aneupl_4N_0_%d[0],Aneupl_4N_1_%d[0],Aneupl_4N_2_%d[0],Aneupl_4N_3_%d[0],Aneupl_4N_4_%d[0],Aneupl_4N_5_%d[0])\n",i,i,i,i,i,i,i);
        fprintf(f,"WT_Telophase_4N_row%d = (Telophase_4N_0_%d[0],Telophase_4N_1_%d[0],Telophase_4N_2_%d[0],Telophase_4N_3_%d[0],Telophase_4N_4_%d[0],Telophase_4N_5_%d[0])\n",i,i,i,i,i,i,i);
    }
    
    fprintf(f,"Percent_signal = (");
    for (unsigned short int i=0; i<5;i++) fprintf(f,"%.2lf, ",(5+i*15)/100.);
    fprintf(f,"%.2lf)\n",(5+5*15)/100.);
    
    
    fprintf(f,"fig, axes = plt.subplots(1, 4, sharey='all')\n");
    for (unsigned short int i=0; i<=3;i++) {
        fprintf(f,"bottom1 = WT_Normal_cc_row%d\n",i);
        
        fprintf(f,"bottom2 = (Normal_cc_0_%d[0] + G2_4N_slip_0_%d[0],                                        Normal_cc_1_%d[0] + G2_4N_slip_1_%d[0], Normal_cc_2_%d[0]+ G2_4N_slip_2_%d[0],                 Normal_cc_3_%d[0] + G2_4N_slip_3_%d[0], Normal_cc_4_%d[0]+ G2_4N_slip_4_%d[0], Normal_cc_5_%d[0]+ G2_4N_slip_5_%d[0])\n",i,i,i,i,i,i,i,i,i,i,i,i);
        
        fprintf(f,"bottom3 = (Normal_cc_0_%d[0] + G2_4N_slip_0_%d[0] + Aneupl_4N_0_%d[0],                  Normal_cc_1_%d[0] + G2_4N_slip_1_%d[0]+ Aneupl_4N_1_%d[0],                                      Normal_cc_2_%d[0] + G2_4N_slip_2_%d[0]+ Aneupl_4N_2_%d[0],                                       Normal_cc_3_%d[0] + G2_4N_slip_3_%d[0]+ Aneupl_4N_3_%d[0],                                       Normal_cc_4_%d[0] + G2_4N_slip_4_%d[0]+ Aneupl_4N_4_%d[0],                                       Normal_cc_5_%d[0] + G2_4N_slip_5_%d[0]+ Aneupl_4N_5_%d[0])\n",i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i);
        
        fprintf(f,"bottom4 = (Normal_cc_0_%d[0] + G2_4N_slip_0_%d[0] + Aneupl_4N_0_%d[0] + Telophase_4N_0_%d[0], Normal_cc_1_%d[0] + G2_4N_slip_1_%d[0] + Aneupl_4N_1_%d[0] + Telophase_4N_1_%d[0],              Normal_cc_2_%d[0] + G2_4N_slip_2_%d[0] + Aneupl_4N_2_%d[0] + Telophase_4N_2_%d[0],              Normal_cc_3_%d[0] + G2_4N_slip_3_%d[0] + Aneupl_4N_3_%d[0] + Telophase_4N_3_%d[0],             Normal_cc_4_%d[0] + G2_4N_slip_4_%d[0] + Aneupl_4N_4_%d[0] + Telophase_4N_4_%d[0],             Normal_cc_5_%d[0] + G2_4N_slip_5_%d[0] + Aneupl_4N_5_%d[0] + Telophase_4N_5_%d[0])\n",i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i);
        
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Normal_cc_row%d, width)\n",i,i);
        fprintf(f,"axes[%d].bar(Percent_signal, WT_G2_4N_slip_row%d, width, bottom=bottom1)\n",i,i);
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Aneupl_4N_row%d, width, bottom=bottom2)\n",i,i);
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Telophase_4N_row%d, width, bottom=bottom3)\n",i,i);
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Apoptosis_row%d, width, bottom=bottom4)\n",i,i);
    }
    
    fprintf(f,"fig.gca().yaxis.set_minor_formatter(NullFormatter())\n");
    fprintf(f,"fig.subplots_adjust(top=0.92, bottom=0.1, left=0.13, right=0.95, hspace=0.25,wspace=0.45)\n");
    
    fprintf(f,"axes[0].set_xlabel('%% GF_H (Trail=0)',fontsize=9)\n");
    fprintf(f,"axes[1].set_xlabel('%% Trail (GFH=1)',fontsize=9)\n");
    fprintf(f,"axes[2].set_xlabel('%% Trail (GFH=0)',fontsize=9)\n");
    fprintf(f,"axes[3].set_xlabel('%% GF (GFH=0)',fontsize=9)\n");
    
    fprintf(f,"fig.savefig('%s/%s/_EXP/WT_RESULTS_Assync_%d.pdf')\n",MODEL_Directory,PD->BRN->path,ASSYNC_BIAS);

    fclose(f);
    
    sprintf(fn,"python %s/%s/_EXP/%s_Assync_%d_BARGRAPHS.py\n",MODEL_Directory_system,PD->BRN->path,nodename,ASSYNC_BIAS);
    system(fn);
}

void run_full_partial_lockdown_analysis(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,const char nodename[], int lock_value, longint_ARRAY *gl_keep_asis,int phyton_only){
    
    if (phyton_only==1) {
        if (lock_value==0) Generate_and_RUN_Phyton_script_BAR_GRAPHS(PD, nodename,gl_keep_asis,0);
        else Generate_and_RUN_Phyton_script_BAR_GRAPHS(PD, nodename,gl_keep_asis,1);
        return;
    }
    
    unsigned long int T_sample = 50000;
    longint_ARRAY_pair *Hits;
    
    
    Hits=new longint_ARRAY_pair();
    if (lock_value==0) Hits->add_element(PD->BRN->MDAT->get_node_ID(nodename), 0);
                 else  Hits->add_element(PD->BRN->MDAT->get_node_ID(nodename), 1);

    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.2,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.35,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.50,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.65,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.80,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.95,Hits,PD,gl_keep_asis,T_sample,1);

    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.2,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.35,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.50,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.65,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.80,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.95,Hits,PD,gl_keep_asis,T_sample,1);


    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.2,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.35,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.50,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.65,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.80,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.95,Hits,PD,gl_keep_asis,T_sample,0);

    run_cell_cycle_stochastic_hits_SWEEP("GF",0.2,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.35,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.50,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.65,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.80,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.95,Hits,PD,gl_keep_asis,T_sample,0);

    delete Hits; Hits=NULL;
    if (lock_value==0) Generate_and_RUN_Phyton_script_BAR_GRAPHS(PD, nodename,gl_keep_asis,0);
    else Generate_and_RUN_Phyton_script_BAR_GRAPHS(PD, nodename,gl_keep_asis,1);
}

void run_lockdown_control_async(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,longint_ARRAY *gl_keep_asis,int phyton_only){
    
    unsigned long int T_sample = 50000;
   
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.40,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.60,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.80,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",1.00,NULL,PD,gl_keep_asis,T_sample,1);

    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.40,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.60,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.80,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",1.00,NULL,PD,gl_keep_asis,T_sample,1);

    
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.2,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.35,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.50,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.65,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.80,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.95,NULL,PD,gl_keep_asis,T_sample,1);
    
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.2,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.35,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.50,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.65,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.80,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.95,NULL,PD,gl_keep_asis,T_sample,1);
    
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.2,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.35,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.50,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.65,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.80,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.95,NULL,PD,gl_keep_asis,T_sample,0);
    
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.2,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.35,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.50,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.65,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.80,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.95,NULL,PD,gl_keep_asis,T_sample,0);
}

void run_lockdown_control_Sync(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,longint_ARRAY *gl_keep_asis,int phyton_only){
    
    unsigned long int T_sample = 50000;
    
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.40,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.60,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.80,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",1.00,NULL,PD,gl_keep_asis,T_sample,1);
    
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.40,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.60,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.80,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",1.00,NULL,PD,gl_keep_asis,T_sample,1);
    
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.2,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.35,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.50,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.65,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.80,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.95,NULL,PD,gl_keep_asis,T_sample,1);
    
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.2,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.35,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.50,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.65,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.80,NULL,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.95,NULL,PD,gl_keep_asis,T_sample,1);
    
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.2,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.35,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.50,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.65,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.80,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("Trail",0.95,NULL,PD,gl_keep_asis,T_sample,0);
    
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.2,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.35,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.50,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.65,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.80,NULL,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP("GF",0.95,NULL,PD,gl_keep_asis,T_sample,0);
}


void run_full_partial_lockdown_analysis_asychronous_update(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,const char nodename[], int lock_value, longint_ARRAY *gl_keep_asis,int phyton_only){
    
    if (phyton_only==1) {
        if (lock_value==0) Generate_and_RUN_Phyton_script_BAR_GRAPHS_Assync(PD, nodename,gl_keep_asis,0);
        else Generate_and_RUN_Phyton_script_BAR_GRAPHS_Assync(PD, nodename,gl_keep_asis,1);
        return;
    }
    
    unsigned long int T_sample = 50000;
    longint_ARRAY_pair *Hits;
    
    
    Hits=new longint_ARRAY_pair();
    if (!strcmp(nodename, "FoxO")){
        Hits->add_element(PD->BRN->MDAT->get_node_ID("FoxO1"), lock_value);
        Hits->add_element(PD->BRN->MDAT->get_node_ID("FoxO3"), lock_value);
    }
    else{ Hits->add_element(PD->BRN->MDAT->get_node_ID(nodename), lock_value);
    }

    if (!strcmp(nodename, "PI3K_H")) {
        Hits->add_element(PD->BRN->MDAT->get_node_ID("p110_H"), lock_value);
    }
    
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.40,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.60,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.80,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",1.00,Hits,PD,gl_keep_asis,T_sample,1);

    
    
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.2,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.35,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.50,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.65,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.80,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF_High",0.95,Hits,PD,gl_keep_asis,T_sample,1);

    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.2,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.35,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.50,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.65,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.80,Hits,PD,gl_keep_asis,T_sample,1);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.95,Hits,PD,gl_keep_asis,T_sample,1);


    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.2,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.35,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.50,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.65,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.80,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("Trail",0.95,Hits,PD,gl_keep_asis,T_sample,0);

    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.2,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.35,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.50,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.65,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.80,Hits,PD,gl_keep_asis,T_sample,0);
    run_cell_cycle_stochastic_hits_SWEEP_Assync("GF",0.95,Hits,PD,gl_keep_asis,T_sample,0);

    delete Hits; Hits=NULL;
    if (lock_value==0) Generate_and_RUN_Phyton_script_BAR_GRAPHS_Assync(PD, nodename,gl_keep_asis,0);
    else Generate_and_RUN_Phyton_script_BAR_GRAPHS_Assync(PD, nodename,gl_keep_asis,1);
}

void build_SUBSET_and_SinglegenePerturbation_Model_and_sample_state_space(const char modelname_REF[],
                                            const char modules_name[],
                                            const char subset[],
                                            int subset_mode,
                                            const char konode[],
                                            int setNodeTO){
    Boolean_RegNet *A,*B,*A_mut,*B_mut;
    MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,*PD_B,*PD_A_mut,*PD_B_mut;
    
    longint_ARRAY *gl_keep_asis_A,*gl_keep_asis_B;
    longint_ARRAY_pair *Hits;
    
    A=new Boolean_RegNet(modelname_REF);
    gl_keep_asis_A=get_input_node_names(A);
    A->load_Module_IDs_and_Drawing_Order();

//    int m1=A->load_Module_Drawing_order(modules_name);
//    A->load_Module_IDs(modules_name,m1);

    PD_A=build_and_sample(A, gl_keep_asis_A, modules_name,0,A->name,0); // first 0 or 1 is modules
//            PD_A=new MODULAR_Boolean_Dynamics_PAIRED_TRACKER(A,N_MAX);
//            strcpy(PD_A->Coupled->DIRname,A->name);
//            PD_A->read_Node_Cytoskape_positions();

    run_OneNode_pulse_experiment_MOD("GF_High",52,1,PD_A,0.000000001,30);
    run_OneNode_pulse_experiment_MOD("GF_High",55,1,PD_A,0.000000001,30);

//    run_OneNode_pulse_experiment_MOD("GF_High",18,1,PD_A,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("ROS_ext",27,1,PD_A,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("Gamma",25,1,PD_A,0.000000001,2);
//    run_OneNode_pulse_experiment_MOD("ROS_ext",25,1,PD_A,0.000000001,2);
//
//    run_OneNode_pulse_experiment_MOD("Gamma",50,1,PD_A,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("Gamma",62,1,PD_A,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("Gamma",65,1,PD_A,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("Gamma",56,1,PD_A,0.000000001,150);
 //  run_OneNode_pulse_experiment_MOD("GF_High",39,1,PD_A,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("GF_High",35,1,PD_A,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("GF_High",57,1,PD_A,0.000000001,150);
 //   run_OneNode_pulse_experiment_MOD("GF_High",27,1,PD_A,0.000000001,150);

//    run_OneNode_pulse_experiment_MOD("GF",6,1,PD_A,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("GF",16,1,PD_A,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("GF",37,1,PD_A,0.000000001,150);
  //  run_OneNode_pulse_experiment_MOD("ROS_ext",25,1,PD_A,0.000000001,4);
//
    
//    run_transition_from_cycle(28,"Gamma",PD_A,gl_keep_asis_A);   // exit(1);
//    run_transition_from_cycle(28,"ROS_ext",PD_A,gl_keep_asis_A);
//    exit(1);
////
   //   run_OneNode_pulse_experiment_MOD("GF_High",33,1,PD_A,0.000000001,60); exit(1);
    
//    longint_ARRAY_pair *Hits;
//    Hits=new longint_ARRAY_pair();
//    Hits->add_element(PD_A->BRN->MDAT->get_node_ID(konode), setNodeTO);
//    Hits->add_element(PD_A->BRN->MDAT->get_node_ID("p110_H"), setNodeTO);
//    run_cell_cycle_Delayed_HITS(Hits,PD_A,gl_keep_asis_A);
//    delete Hits; Hits=NULL;
//
    //   run_OneNode_pulse_experiment_MOD("GF_High",19,1,PD_A,0.000000001,60);// exit(1);
    //   run_OneNode_pulse_experiment_MOD("GF_High",40,1,PD_A,0.000000001,60);// exit(1);
    //   run_OneNode_pulse_experiment_MOD("GF",53,1,PD_A,0.000000001,7); exit(1);
    //
   // run_transition_from_cycle(14,"Gamma",PD_A,gl_keep_asis_A);
   // run_transition_from_cycle(14,"GF_High",PD_A,gl_keep_asis_A);
   // run_transition_from_cycle(14,"Trail",PD_A,gl_keep_asis_A);
    
    B=generate_SUBSET_RNB(A,subset,subset_mode);
    gl_keep_asis_B=get_input_node_names(B);
    PD_B=build_and_sample(B, gl_keep_asis_B, modules_name,1,A->name,1); // first 0 or 1 is modules
    char fn[400];
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD_B->BRN->path);
    system(fn);
   
    
  //  Hits=new longint_ARRAY_pair(); run_stochastic_cell_cycle("GF_High", 0.7, Hits, PD_B, gl_keep_asis_B, 100);
    
    run_OneNode_pulse_experiment_MOD("GF_High",49,1,PD_B,0.000000001,30);
    run_OneNode_pulse_experiment_MOD("GF_High",62,1,PD_B,0.000000001,30);
    run_OneNode_pulse_experiment_MOD("GF_High",63,1,PD_B,0.000000001,30);

  //   run_full_partial_lockdown_analysis(PD_B,"Plk1",0,gl_keep_asis_B,0);
  //    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"Plk1",gl_keep_asis_B,1);
   //   run_full_partial_lockdown_analysis_asychronous_update(PD_B,"Plk1",gl_keep_asis_B,0);
    exit(1);
    
//    setup_experiment_with_asynchronous_update(PD_B,7,250,1);
//    setup_experiment_with_asynchronous_update(PD_B,7,251,1);
//    setup_experiment_with_asynchronous_update(PD_B,7,253,1);
 //   setup_experiment_with_asynchronous_update(PD_B,7,254,1);
  
//    setup_experiment_with_asynchronous_update(PD_B,7,3,1000);
//    setup_experiment_with_asynchronous_update(PD_B,7,4,1000);
//    setup_experiment_with_asynchronous_update(PD_B,7,5,1000);
//    //setup_experiment_with_asynchronous_update(PD_B,7,250,1000);

 //   run_full_partial_lockdown_analysis_asychronous_update(PD_B,"Plk1",0,gl_keep_asis_B,0);
 //   run_full_partial_lockdown_analysis_asychronous_update(PD_B,"p110_H",1,gl_keep_asis_B,0);
 //   run_full_partial_lockdown_analysis_asychronous_update(PD_B,"PI3K_H",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"AKT_H",1,gl_keep_asis_B,0);
 //   run_full_partial_lockdown_analysis_asychronous_update(PD_B,"p110_H",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"FoxO1",0,gl_keep_asis_B,0);

    
    //knockout /OE runs
    
  
//    run_full_partial_lockdown_analysis(PD_B,"Ect2",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"Ect2",0,gl_keep_asis_B,0);

//    run_full_partial_lockdown_analysis(PD_B,"MCL-1",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"MCL-1",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"Casp2",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"Casp2",1,gl_keep_asis_B,0);
//
//    run_full_partial_lockdown_analysis(PD_B,"Rheb",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"Rheb",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"FoxO1",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"FoxO1",1,gl_keep_asis_B,0);

    
    run_full_partial_lockdown_analysis(PD_B,"BCLXL",1,gl_keep_asis_B,0);
    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"BCLXL",1,gl_keep_asis_B,0);
    
//    //
//    run_full_partial_lockdown_analysis(PD_B,"p27Kip1",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"p27Kip1",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"Myc",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"Myc",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"CyclinD1",1,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis_asychronous_update(PD_B,"CyclinD1",1,gl_keep_asis_B,0);
//
    
//    run_full_partial_lockdown_analysis(PD_B,"PLCgamma",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"NeddL4",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"Ca2+",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"PDK1",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"ERK",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"TSC2",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"mTORC1",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"mTORC2",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"S6K",0,gl_keep_asis_B,0);
//
//        run_full_partial_lockdown_analysis(PD_B,"Myc",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"CyclinD1",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"p27Kip1",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"pRB",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"p21",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"E2F1",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"CyclinE",0,gl_keep_asis_B,0);
//
//        run_full_partial_lockdown_analysis(PD_B,"FoxM1",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"CyclinA",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"Cdc25A",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"Cdh1",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"Cdc25C",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"Cdc25B",0,gl_keep_asis_B,0);
//
//        run_full_partial_lockdown_analysis(PD_B,"Emi1",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"CyclinB",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"Cdc20",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"Mad2",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"CHK1",0,gl_keep_asis_B,0);
//
//        run_full_partial_lockdown_analysis(PD_B,"BAD",0,gl_keep_asis_B,0);
//        run_full_partial_lockdown_analysis(PD_B,"DR4_5",0,gl_keep_asis_B,0);
//
    run_full_partial_lockdown_analysis(PD_B,"Casp8",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"MCL-1",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"BCL2",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"BCLXL",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"BID",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"BAK",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"BAX",0,gl_keep_asis_B,0);
//    run_full_partial_lockdown_analysis(PD_B,"IAPs",0,gl_keep_asis_B,0);

    
 //   run_lockdown_control_async(PD_B,gl_keep_asis_B,0);
  //  run_lockdown_control_Sync(PD_B,gl_keep_asis_B,0);

   //
   // run_full_partial_lockdown_analysis_asychronous_update(PD_B,"Cdc25A",0,gl_keep_asis_B,0);
   // run_full_partial_lockdown_analysis_asychronous_update(PD_B,"Emi1",0,gl_keep_asis_B,0);
    exit(1);
    
    //   run_full_partial_lockdown_analysis_asychronous_update(PD_B,"p110_H",gl_keep_asis_B,0);
    
//    run_full_partial_lockdown_analysis(PD_B,"FoxM1",gl_keep_asis_B,1);
 //    run_full_partial_lockdown_analysis(PD_B,"FoxO3",gl_keep_asis_B,0);
    
 //    run_full_partial_lockdown_analysis(PD_B,"CyclinA",gl_keep_asis_B,0);
    
 //   PD_B->Coupled->Export_ALL_Attractors_onto_Boolean_Network(1);

   //  run_OneNode_pulse_experiment_MOD("GF_High",7,1,PD_B,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("GF_High",9,1,PD_B,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("GF_High",10,1,PD_B,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("GF_High",15,1,PD_B,0.000000001,150);
//    run_OneNode_pulse_experiment_MOD("GF_High",14,1,PD_B,0.000000001,150);
  
     //   setup_experiment_with_asynchronous_update(PD_B,6,250,1);
    
//     setup_experiment_with_asynchronous_update(PD_B,7,254,1);
//
    
//    run_transition_from_cycle(6,"GF",PD_B,gl_keep_asis_B);
//    run_transition_from_cycle(7,"GF",PD_B,gl_keep_asis_B);
//
    exit(1);
 //    check_gates(PD_B->Coupled); exit(1);
    
    compare_attractor_graphs(PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
    Export_Attractor_Life_and_Death_table(PD_A,gl_keep_asis_A,1);
    Export_Attractor_Life_and_Death_table(PD_B,gl_keep_asis_B,1);
  //  COMPARE_Attractor_Life_and_Death_table(PD_A,gl_keep_asis_A,PD_B,gl_keep_asis_B);
    PD_A->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis_B->N)),p_error);
    PD_B->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis_B->N)),p_error);
    

    run_OneNode_pulse_experiment_MOD("GF_High",19,1,PD_B,0.000000001,50);
 //   run_OneNode_pulse_experiment_MOD("GF",15,1,PD_B,0.000000001,50);
   // run_OneNode_pulse_experiment_MOD("Trail",15,1,PD_B,0.000000001,50);
  //    exit(1);

    run_cell_cycle_entry_exit(PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B,0);
  //  run_cell_cycle_entry_exit(PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B,3);
    
 ///   run_OneNode_pulse_experiment_MOD("Trail",14,1,PD_B,0.000000001,21);
 //   run_OneNode_pulse_experiment_MOD("Trail",14,1,PD_B,0.000000001,22);
//    run_OneNode_pulse_experiment_MOD("Trail",8,1,PD_B,0.000000001,24);
    //  run_cell_cycle_Pulse("Gamma",21,PD_B,gl_keep_asis_B);

 //   run_transition_from_cycle(13,"Trail",PD_B,gl_keep_asis_B);
  //  run_transition_from_cycle(18,"GF_High",PD_B,gl_keep_asis_B);exit(1);
//    run_transition_from_cycle(13,"Trail",PD_B,gl_keep_asis_B);
  //  run_transition_from_cycle(13,"GF",PD_B,gl_keep_asis_B);
//
  //
    
   //
 // run_OneNode_pulse_experiment_MOD("Trail",13,1,PD_B,0.000000001,50);
    
//    Hits=new longint_ARRAY_pair();
//    Hits->add_element(PD_B->BRN->MDAT->get_node_ID("Plk1"), 0);
//    for(unsigned int ll=24;ll<=32;ll++)
//        run_cell_cycle_Entry_Delayed_HITS(Hits,PD_B,gl_keep_asis_B,ll,66);
//    delete Hits; Hits=NULL;
//
//    Hits=new longint_ARRAY_pair();
//    Hits->add_element(PD_B->BRN->MDAT->get_node_ID("p110_H"), 0);
//    for(unsigned int ll=4;ll<=7;ll++){
//        run_cell_cycle_Entry_Delayed_HITS(Hits,PD_B,gl_keep_asis_B,ll,66);
//        run_cell_cycle_Entry_Delayed_HITS(Hits,PD_B,gl_keep_asis_B,ll,71);
//        run_cell_cycle_Entry_Delayed_HITS(Hits,PD_B,gl_keep_asis_B,ll,72);
//    }
//    delete Hits; Hits=NULL;
//    
    exit(1);
//    
    
    
    
    Hits=new longint_ARRAY_pair();
    
//    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.2,Hits,PD_B,gl_keep_asis_B,500000);
//    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.4,Hits,PD_B,gl_keep_asis_B,500000);
//    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.6,Hits,PD_B,gl_keep_asis_B,500000);
//    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.8,Hits,PD_B,gl_keep_asis_B,500000);
//    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.9,Hits,PD_B,gl_keep_asis_B,500000);
//    run_cell_cycle_stochastic_hits_SWEEP("GF_High",1,Hits,PD_B,gl_keep_asis_B,500000);
//
//    run_stochastic_cell_cycle_SWEEP("GF_High",Hits,PD_B,gl_keep_asis_B,500000);
  //  exit(1);

    Hits->add_element(PD_B->BRN->MDAT->get_node_ID(konode), setNodeTO);
    
//    for(unsigned int ll=24;ll<=34;ll++){
//        run_cell_cycle_Entry_Delayed_HITS(Hits,PD_B,gl_keep_asis_B,ll,44);
// //       run_cell_cycle_Entry_Delayed_HITS(Hits,PD_B,gl_keep_asis_B,ll,45);
//    }
//    exit(1);
//
 //   run_stochastic_cell_cycle_SWEEP("GF_High",Hits,PD_B,gl_keep_asis_B,500000);

    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.2,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.4,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.6,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.8,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.9,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",1,Hits,PD_B,gl_keep_asis_B,500000,1);

   exit(1);
    

//    Hits->add_element(PD_B->BRN->MDAT->get_node_ID("PI3K_H"), 1);
  Hits->add_element(PD_B->BRN->MDAT->get_node_ID("FoxO1"), 0);
 //     Hits->add_element(PD_B->BRN->MDAT->get_node_ID("Plk1_H"), 1);

    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.2,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.4,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.6,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.8,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",0.9,Hits,PD_B,gl_keep_asis_B,500000,1);
    run_cell_cycle_stochastic_hits_SWEEP("GF_High",1,Hits,PD_B,gl_keep_asis_B,500000,1);
    //
    run_stochastic_cell_cycle_SWEEP("GF_High",Hits,PD_B,gl_keep_asis_B,500000);

    exit(1);
   
  //  run_Models_are_Awesome_checklist(PD_A,PD_B,gl_keep_asis_A,gl_keep_asis_B);
 
   
    
    A_mut= generate_MUTANT_RNB(A, konode,setNodeTO);
    B_mut= generate_MUTANT_RNB(B, konode,setNodeTO);
    
    PD_A_mut=build_and_sample(A_mut, gl_keep_asis_A, modules_name,0,A->name,1);
    PD_B_mut=build_and_sample(B_mut, gl_keep_asis_B, modules_name,0,A->name,1);
    
    compare_attractor_graphs(PD_A_mut,PD_B_mut,gl_keep_asis_A,gl_keep_asis_B);
    Export_Attractor_Life_and_Death_table(PD_A_mut,gl_keep_asis_A,1);
    Export_Attractor_Life_and_Death_table(PD_B_mut,gl_keep_asis_B,1);
    COMPARE_Attractor_Life_and_Death_table(PD_A_mut,gl_keep_asis_A,PD_B_mut,gl_keep_asis_B);
    
    compare_attractor_graphs(PD_A,PD_A_mut,gl_keep_asis_A,gl_keep_asis_A);
    compare_attractor_graphs(PD_B,PD_B_mut,gl_keep_asis_B,gl_keep_asis_B);
    COMPARE_Attractor_Life_and_Death_table(PD_A,gl_keep_asis_A,PD_A_mut,gl_keep_asis_A);
    COMPARE_Attractor_Life_and_Death_table(PD_B,gl_keep_asis_B,PD_B_mut,gl_keep_asis_B);
    
    PD_A_mut->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis_B->N)),p_error);
    PD_B_mut->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis_B->N)),p_error);
    
}

/* Compare_Boolean_Models_h */
