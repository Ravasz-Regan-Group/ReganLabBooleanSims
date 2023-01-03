//
//  Run_Virtual_Exp_Options.h.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 6/7/19.
//  Copyright © 2019 Regan_Group. All rights reserved.
//

#ifndef Run_Virtual_Exp_Options_h_h
#define Run_Virtual_Exp_Options_h_h

#define T_max_CellCycleSample   50000  //
#define T_max_TransitionsSample    50000 //
#define BOX_ENV    20
#define BOX_KO_OE   10
#define T_MAX_DRAW 650
#define TW  3
#define NameBuff    55
#define SHOW_Init_state 30

#include "Virtual_Exp_classes.h"
#include "Python_Figures_module.h"
#include "Timecourses_and_stats_module.h"



void run_Modules(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    char attr_list[200][5000];
    
    std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(linenow, '(',')');
    char bla2[300]; sprintf(bla2,"%s",word[0].c_str());
    unsigned long int Nf=separate_line_by_string(bla2,attr_list,", ",200);
    if (Nf > 0){
        for(unsigned int i=1;i<=Nf;i++)
            for(unsigned int mid=1;mid<=PD->BRN->MDAT->Module_NR;mid++){
                printf("mid = %d\t %s\n",mid, PD->BRN->Module_Names[mid].c_str());
                if(!strcmp(attr_list[i], PD->BRN->Module_Names[mid].c_str()))
                    work_with_individual_Modules(mid,PD,0.0000001);
            }
    }
  //  printf ("Nice try 1\n"); getchar();
}

void run_Async_Module_Cycles_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    unsigned int STATS=1;
    istringstream iss(linenow);
    string sr; iss >> sr;
    iss >> STATS;
    if (iss.fail() == true) STATS = 1;
    
    for(unsigned int mid=1;mid<=PD->BRN->MDAT->Module_NR;mid++){
        Async_Module_cycle(mid, PD,0.0000001,STATS);
      //  getchar();
    }
   // printf ("Nice try 2\n");
}

longint_ARRAY *get_requested_cell_group(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l; l=NULL;
    
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr;
    iss >> sr;
    switch (sr[0]) {
        case 'P': l=get_Healthy_CellCycle_attractors(PD);
            break;
        case 'Q': l=get_Healthy_G0_attractors(PD);
            break;
        case 'A': l=get_Arrested_G0_attractors(PD);
            break;
        case 'M': l=get_Mesenhymal_attractors(PD);
            break;
        case 'H': l=get_Hybrid_attractors(PD);
            break;
        case 'E': l=get_Epithelial_attractors(PD);
            break;
        case 'S':
            switch (sr[1]) {
                case 'a': l=get_ALL_Senescent_attractors(PD); break;
                case 'A': l=get_ALL_Senescent_attractors(PD); break;
                case 'c': l=get_Chromosomal_Senescent_attractors(PD); break;
                case 'C': l=get_Chromosomal_Senescent_attractors(PD); break;
                case 'm': l=get_Mitochondrial_Senescent_attractors(PD); break;
                case 'M': l=get_Mitochondrial_Senescent_attractors(PD); break;
                case 'd': l=get_Double_Senescent_attractors(PD); break;
                case 'D': l=get_Double_Senescent_attractors(PD); break;
            }
            break;
        default: {
            istringstream iss2(sr);
            unsigned long int Attr_ID;
            iss2 >> Attr_ID;
            if (iss2.fail() == true){
                printf("No options implemented for phenotype group %s\n",sr.c_str()); getchar(); l = NULL;
            }
            else {
                if(Attr_ID <=PD->Coupled->Attractor_NR) {
                    l= new longint_ARRAY(); l->add_element(Attr_ID);
                }
                else {printf("Model only has %lu attractors (Attractor %ld requested)\n",PD->Coupled->Attractor_NR,Attr_ID); getchar(); l = NULL;}
            }
        }
            break;
    }
    return (l);
}

long int draw_timetrace_Multi_INPUT_pulse_MOD(int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                        double p_error,int delay,unsigned long int pulse_length,longint_ARRAY *inp_index,FILE *f,double y0,int W){
    unsigned long int new_id;
    int *s_start;
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    
    s_start= new int[inp_index->N+1];
    for(unsigned int l=1;l<=inp_index->N;l++){
        s_start[l]=PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]-1];
    }
    
    int s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];
    
    for (int t=1;t<=T_MAX_DRAW; t++) {
        if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length))
            for(unsigned int l=1;l<=inp_index->N;l++)
                PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]-1]=48+(49-s_start[l]);
        
        if (t>=SHOW_Init_state+delay+pulse_length) {
            for(unsigned int l=1;l<=inp_index->N;l++)
                PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]-1]=s_start[l];
            
            if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
            if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
        }
        
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length)){
            for(unsigned int l=1;l<=inp_index->N;l++)
                PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]-1]=48+(49-s_start[l]);
        }
        if (t>=SHOW_Init_state+delay+pulse_length) {
            for(unsigned int l=1;l<=inp_index->N;l++)
                PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]-1]=s_start[l];
            if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
            if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
        }
        draw_timestep_MOD(PD,t,f,y0,W, T_MAX_DRAW, TW, NameBuff);
    }
    new_id= PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
    return (new_id);
}


void run_MultiNode_SAMEpulse_experiment_MOD(longint_ARRAY *inp_index,int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                      double p_error,unsigned long int pulse_l){
    int W=8,YTOP,XTOP;
    long int id_new;
    FILE *f;
    char fn[400],fnd[400],fnd2[400];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a+1;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*650+5;
    
    sprintf(fnd,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fnd);
    sprintf(fnd,"mkdir %s/%s/_EXP/Multi-pulse",MODEL_Directory_system,PD->BRN->path);
    for(unsigned int l=1;l<=inp_index->N;l++)
        sprintf(fnd,"%s_%s",fnd,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]]);
    sprintf(fnd,"%s_pulse-%ld\n",fnd,pulse_l);
    system(fnd);
    
    sprintf(fnd,"%s/%s/_EXP/Multi-pulse",MODEL_Directory,PD->BRN->path);
    for(unsigned int l=1;l<=inp_index->N;l++)
        sprintf(fnd,"%s_%s",fnd,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]]);
    sprintf(fnd,"%s_pulse-%ld",fnd,pulse_l);
  
    sprintf(fnd2,"%s/%s/_EXP/Multi-pulse",MODEL_Directory_system,PD->BRN->path);
    for(unsigned int l=1;l<=inp_index->N;l++)
        sprintf(fnd2,"%s_%s",fnd2,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]]);
    sprintf(fnd2,"%s_pulse-%ld",fnd2,pulse_l);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/Multi-pulse",fnd);
        for(unsigned int l=1;l<=inp_index->N;l++)
            sprintf(fn,"%s_%s",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]]);
        sprintf(fn,"%s_pulse-%ld_Mpulse_A%d_As-%d_delay-%d.eps",fn,pulse_l,a_st,a_st_k,delay);
        
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        id_new=draw_timetrace_Multi_INPUT_pulse_MOD(a_st,a_st_k,PD,p_error,delay,pulse_l,inp_index,f,YTOP,W);
        //getchar();
        sprintf(fn,"%s/Multi-pulse",fnd2);
        for(unsigned int l=1;l<=inp_index->N;l++)
            sprintf(fn,"%s_%s",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[inp_index->A[l]]]);
        sprintf(fn,"%s_pulse-%ld_Mpulse_A%d_As-%d_delay-%d",fn,pulse_l,a_st,a_st_k,delay);
        close_EPSDRAW_file(f,fn,2);
    }
    
}


void run_Pulse1_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l;
    unsigned int inputL;
    
    l=get_requested_cell_group(PD,linenow);
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    if (PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str())<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    iss >> inputL;
    if (iss.fail() == true){
        printf("Invalid input length in line %s\n",linenow); getchar(); l = NULL;
    }
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_OneNode_pulse_experiment_MOD(inputname.c_str(),(int)(l->A[j]),1,PD,0.000000001,inputL);
        }
    }
    delete l;l=NULL;
}

void run_Pulse1_Cycle(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l;
    unsigned int inputL_ON,inputL_OFF;
    
    l=get_requested_cell_group(PD,linenow);
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    if (PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str())<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    iss >> inputL_ON;
    if (iss.fail() == true){
        printf("Invalid input ON-length in line %s\n",linenow); getchar(); l = NULL;
    }
    iss >> inputL_OFF;
    if (iss.fail() == true){
        printf("Invalid input OFF-length in line %s\n",linenow); getchar(); l = NULL;
    }
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
          //  run_OneNode_pulse_experiment_MOD(inputname.c_str(),(int)(l->A[j]),1,PD,0.000000001,inputL);
            run_OneNode_pulse_Cycle_experiment_MOD(inputname.c_str(),(int)(l->A[j]),1,PD,0,inputL_ON,inputL_OFF);
        }
    }
    delete l;l=NULL;
}

void run_Pulse2_Cycles(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l;
    unsigned int inputL_ON_1,inputL_OFF_1,inputL_ON_2,inputL_OFF_2;
    
    l=get_requested_cell_group(PD,linenow);
    istringstream iss(linenow);
    string sr,inputname_1,inputname_2;
    iss >> sr; iss >> sr;
    iss >> inputname_1;
    if (PD->Coupled->BRN->MDAT->get_node_ID(inputname_1.c_str())<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname_1.c_str(),linenow);
        return;
    }
    iss >> inputL_ON_1;
    if (iss.fail() == true){
        printf("Invalid input ON_1-length in line %s\n",linenow); getchar(); l = NULL;
    }
    iss >> inputL_OFF_1;
    if (iss.fail() == true){
        printf("Invalid input OFF_1-length in line %s\n",linenow); getchar(); l = NULL;
    }
    iss >> inputname_2;
    if (PD->Coupled->BRN->MDAT->get_node_ID(inputname_2.c_str())<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname_2.c_str(),linenow);
        return;
    }
    iss >> inputL_ON_2;
    if (iss.fail() == true){
        printf("Invalid input ON_2-length in line %s\n",linenow); getchar(); l = NULL;
    }
    iss >> inputL_OFF_2;
    if (iss.fail() == true){
        printf("Invalid input OFF_2-length in line %s\n",linenow); getchar(); l = NULL;
    }
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
          //  run_OneNode_pulse_experiment_MOD(inputname.c_str(),(int)(l->A[j]),1,PD,0.000000001,inputL);
            run_Pulse2_Cycles_experiment_MOD(inputname_1.c_str(),inputname_2.c_str(),(int)(l->A[j]),1,PD,0,inputL_ON_1,inputL_OFF_1,inputL_ON_2,inputL_OFF_2);
        }
    }
    delete l;l=NULL;
}


longint_ARRAY *get_requested_pulsenode_group(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *InputList;
    char attr_list[200][5000];
    
    InputList=new longint_ARRAY();
    
    std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(linenow, '(',')');
    char bla2[300]; sprintf(bla2,"%s",word[0].c_str());
    unsigned long int Nf=separate_line_by_string(bla2,attr_list,", ",200);
    if (Nf > 0){
        for(unsigned int i=1;i<=Nf;i++){
            istringstream iss_in(attr_list[i]);
            string inputname;
            iss_in >> inputname;
            unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
            if (in_id>1){
                unsigned long int i_index=0;
                i_index = PD->Coupled->BRN->active_inputs->check_for_element(in_id);
                if (i_index==0){
                    printf("Node  %s designated for multi-pulse experiment in line %s exists in the model, but it is not an INPUT node! \n\t(ignoring line; press enter to move to next instruction\n",inputname.c_str(),linenow); getchar(); return(NULL);
                }
                else InputList->add_element(i_index);
            }
            else {
                printf("Node %s designated for multi-pulse experiment in line %s is not a valid node in this model! \n\t(ignoring line; press enter to move to next instruction\n",inputname.c_str(),linenow); getchar(); return(NULL);
            }
        }
    }
    else {printf("Multi-pulse experiment input node list is correcly specified in line %s\n\t Please list ALL  input variables you wish to toggle at the same time as illustrated here: \n\t (ECM, Stiff_ECM)\n\t(ignoring line; press enter to move to next instruction\n",linenow); getchar();return(NULL);}
    
   return(InputList);
}

void run_PulseM_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l;
    unsigned long int inputL;
    longint_ARRAY *InputList;
    
    l=get_requested_cell_group(PD,linenow);
    InputList=get_requested_pulsenode_group(PD,linenow);
    
    if ((l!=NULL)&&(InputList!=NULL)){
        string inputLength = getLastWord(linenow);
        istringstream iss(inputLength);
        iss >> inputL;
        if (iss.fail() == true){
            printf("Invalid input length in line %s\n",linenow); getchar(); return;
        }
        
        for(int j=1;j<=l->N;j++){
            run_MultiNode_SAMEpulse_experiment_MOD(InputList,(int)(l->A[j]),1,PD,0.000000001,inputL);
        }
        delete l;l=NULL;
        delete InputList;InputList=NULL;
    }
}


void run_PulseScan_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l;
    
    l=get_requested_cell_group(PD,linenow);
    
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    if (PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str())<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }

    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_transition_from_cycle((int)(l->A[j]),inputname.c_str(),PD);
        }
    }
    delete l;l=NULL;
}

void run_NonSaturating_Draw_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l;
    Environment_and_hits *EnvHit;
    
    l=get_requested_cell_group(PD,linenow);
    EnvHit=extract_environment_and_hits(PD,linenow);
    if((EnvHit==NULL)||(EnvHit->p_inputs[0]<0)) return;
        // this is set when the input is incorrent and virtual_exp line gets ignored
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            for(int a_st_k=1;a_st_k<=PD->Coupled->Attractor_valleys[j]->N_a;a_st_k++){
                run_stochastic_inputs_and_KO_OE_hits_experiment_on_single_attractor_MOD(EnvHit,(int)(l->A[j]), a_st_k,PD, 0);
            }
        }
    }
    delete l;l=NULL;
    delete EnvHit; EnvHit=NULL;
}

void run_Async_NonSaturating_Draw_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l;
    Environment_and_hits *EnvHits;
    
    l=get_requested_cell_group(PD,linenow);
    EnvHits=extract_environment_and_hits(PD,linenow);
    if(EnvHits->p_inputs[0]<0) return; // this is set when the input is incorrent and virtual_exp line gets ignored
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            for(int a_st_k=1;a_st_k<=PD->Coupled->Attractor_valleys[j]->N_a;a_st_k++){
                Asynchronous_run_stochastic_inputs_and_KO_OE_hits_experiment_on_single_attractor_MOD(EnvHits,(int)(l->A[j]), a_st_k,PD, 0);
            }
        }
    }
    delete l;l=NULL;
    delete EnvHits; EnvHits=NULL;
}

/*

void run_NonSaturating_Env_Stats_1D_fixed_mutants_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[],int run_assync){
    // NonSaturating_Stats_1    P    ECM    (GF=1, GF_High=0.9, CellDensity_High=0, CellDensity_Low=0, ECM=1, Stiff_ECM=1) (SNAI1 KD 0.5, ZEB1 OE 0.2)

    
    longint_ARRAY *l;
    Environment_and_hits *EnvHit;
    
    l=get_requested_cell_group(PD,linenow);
    
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    unsigned long int input_id=PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
    if (input_id<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    unsigned long int input_index =PD->Coupled->BRN->active_inputs->check_for_element(input_id);
    if (input_index<1) {
        printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    
    EnvHit = extract_environment_and_hits(PD,linenow);
    if(EnvHit->p_inputs[0]<0) return; // this is set when the input is incorrent and virtual_exp line gets ignored

  
    load_Phenotypes_for_model(PD);
    printf("phenotpyes loaded!\n");//getchar();
    
    bool cc_stats,irrev_stats,rev_stats;
    Python_Grapsh_NonSaturating_Stats_1_env(input_index,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
    
    string datadir,fmain;
    longint_ARRAY *lin;
    lin=new longint_ARRAY();
    lin->add_element(input_index);
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    char fn[400];
    sprintf(fn,"mkdir %s\n",fmain.c_str());
    system(fn);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s\n",datadir.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    delete lin; lin=NULL;

    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        EnvHit->p_inputs[input_index]=p_env;
        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
    }
    
    if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
        Python_Grapsh_NonSaturating_Stats_1_env(input_index,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
   
        if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
            printf("Pyhton script cannot find relevant files immeditately AFTER completion of statistical sampling runs!\n");
            getchar();
        }
    }
}

void run_NonSaturating_Stats_1D_EnvGrid_VaryKDOE_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[],int run_assync){
    // NonSaturating_Stats_1DEnv_KOOE_Scan    P   Gamma  (GF=1, GF_High=0.9, CellDensity_High=0, CellDensity_Low=0, Stiff_ECM=1) (Merlin KD)
  
    
    longint_ARRAY *l; // category of states to run it on
    Environment_and_hits *EnvHit;
    
    l=get_requested_cell_group(PD,linenow);
    
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    unsigned long int input_x=PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
    if (input_x<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    unsigned long int input_index_x =PD->Coupled->BRN->active_inputs->check_for_element(input_x);
    if (input_index_x<1) {
        printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! \n",inputname.c_str(),linenow);
        return;
    }
    
    EnvHit = extract_environment_and_hits(PD,linenow);
    if(EnvHit->p_inputs[0]<0) return; // this is set when the input is incorrent and virtual_exp line gets ignored

    load_Phenotypes_for_model(PD);
    printf("phenotpyes loaded!\n");//getchar();
    
    
     bool cc_stats=false,irrev_stats=false,rev_stats=false;
  //  Python_Grapsh_NonSaturating_Stats_1_env_KOOE_scan(input_index_x,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
  
    
    string datadir,fmain;
    longint_ARRAY *lin;
    lin=new longint_ARRAY();
    lin->add_element(input_index_x);
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    char fn[400];
    sprintf(fn,"mkdir %s\n",fmain.c_str());
    system(fn);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s\n",datadir.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    delete lin; lin=NULL;
    
    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        EnvHit->p_inputs[input_index_x]=p_env_1;
        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
    }
  
  //  if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
    //    Python_Grapsh_NonSaturating_Stats_1_env_KOOE_scan(input_index_x,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);

    //    if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
   //         printf("Pyhton script cannot find relevant files immeditately AFTER completion of statistical sampling runs!\n");
   //         getchar();
   //     }
  //  }
}

void run_NonSaturating_Stats_2D_EnvGrid_VaryKDOE_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[],int run_assync){
    // NonSaturating_Stats_2D    P   ECM Trail  (GF=1, GF_High=0.9, CellDensity_High=0, CellDensity_Low=0, Stiff_ECM=1) (Merlin KD)
    // NonSaturating_Stats_2D    P   ECM Trail  (GF=1, GF_High=0.9, CellDensity_High=0, CellDensity_Low=0, Stiff_ECM=1) (Merlin KD, YAP OE)
  
    
    longint_ARRAY *l; // category of states to run it on
    Environment_and_hits *EnvHit;
    
    l=get_requested_cell_group(PD,linenow);
    
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    unsigned long int input_x=PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
    if (input_x<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    unsigned long int input_index_x =PD->Coupled->BRN->active_inputs->check_for_element(input_x);
    if (input_index_x<1) {
        printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! \n",inputname.c_str(),linenow);
        return;
    }
    iss >> sr; iss >> sr;
    iss >> inputname;
    unsigned long int input_y=PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
    if (input_y<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    unsigned long int input_index_y =PD->Coupled->BRN->active_inputs->check_for_element(input_x);
    if (input_index_y<1) {
        printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! \n",inputname.c_str(),linenow);
        return;
    }
    
    EnvHit = extract_environment_and_hits(PD,linenow);
    if(EnvHit->p_inputs[0]<0) return; // this is set when the input is incorrent and virtual_exp line gets ignored

    load_Phenotypes_for_model(PD);
    printf("phenotpyes loaded!\n");//getchar();
    
    bool cc_stats,irrev_stats,rev_stats;
    Python_Grapsh_NonSaturating_Stats_2_env_KOOE_scan(input_index_x,input_index_y,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
    
    string datadir,fmain;
    longint_ARRAY *lin;
    lin=new longint_ARRAY();
    lin->add_element(input_index_x);
    lin->add_element(input_index_y);
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    char fn[400];
    sprintf(fn,"mkdir %s\n",fmain.c_str());
    system(fn);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s\n",datadir.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    delete lin; lin=NULL;

    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        EnvHit->p_inputs[input_index_x]=p_env_1;
        for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
            EnvHit->p_inputs[input_index_y]=p_env_2;
            if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
            if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
            if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        }
    }
    
    if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
        Python_Grapsh_NonSaturating_Stats_2_env_KOOE_scan(input_index_x,input_index_y,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
   
        if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
            printf("Pyhton script cannot find relevant files immeditately AFTER completion of statistical sampling runs!\n");
            getchar();
        }
    }
}

*/

void NonSaturating_Stats_FIX_KDOE_fn1Env(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[],int run_assync){
   // 2 options!!!
    
    // NonSaturating_Stats_FIX_KDOE_fn1Env    P    ECM    (GF=1, GF_High=0.9, CellDensity_High=0, CellDensity_Low=0, ECM=1, Stiff_ECM=1) (SNAI1 KD 0.5, ZEB1 OE 0.2)
            // 1 env on x axis is ECM, KO at the end fixed to these values, level of KO/OE must be specified!
  
    // NonSaturating_Stats_Scan_1Env_fnKDOE 10   GF_High  (GF=1, GF_High=1, Gamma=0) (p21_mRNA KD)
        // 1 env scan, each x axis is the KO/OE of the group listed at the end, KO at the end fixed to these values, level of KO/OE must be specified!

    longint_ARRAY *l;
    Environment_and_hits *EnvHit;
    
    l=get_requested_cell_group(PD,linenow);
    
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    unsigned long int input_id=PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
    if (input_id<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    unsigned long int input_index =PD->Coupled->BRN->active_inputs->check_for_element(input_id);
    if (input_index<1) {
        printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    
    EnvHit = extract_environment_and_hits(PD,linenow);
    if(EnvHit->p_inputs[0]<0) return; // this is set when the input is incorrent and virtual_exp line gets ignored

  
    load_Phenotypes_for_model(PD);
    printf("phenotpyes loaded!\n");//getchar();
    
    bool cc_stats,irrev_stats,rev_stats;
    Python_Grapsh_NonSaturating_Stats_FIX_KDOE_fn1Env(input_index,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
    
    string datadir,fmain;
    longint_ARRAY *lin;
    lin=new longint_ARRAY();
    lin->add_element(input_index);
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    char fn[400];
    sprintf(fn,"mkdir %s\n",fmain.c_str());
    system(fn);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s\n",datadir.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    delete lin; lin=NULL;

    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        EnvHit->p_inputs[input_index]=p_env;
        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
    }
    
    if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
        Python_Grapsh_NonSaturating_Stats_FIX_KDOE_fn1Env(input_index,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
   
        if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
            printf("Pyhton script cannot find relevant files immeditately AFTER completion of statistical sampling runs!\n");
            getchar();
        }
    }
}

void ModelErrors_NonSaturating_Stats_FIX_KDOE_fn1Env(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[],int run_assync){
   // 2 options!!!
    
    // NonSaturating_Stats_FIX_KDOE_fn1Env    P    ECM    (GF=1, GF_High=0.9, CellDensity_High=0, CellDensity_Low=0, ECM=1, Stiff_ECM=1) (SNAI1 KD 0.5, ZEB1 OE 0.2)
            // 1 env on x axis is ECM, KO at the end fixed to these values, level of KO/OE must be specified!
  
    // NonSaturating_Stats_Scan_1Env_fnKDOE 10   GF_High  (GF=1, GF_High=1, Gamma=0) (p21_mRNA KD)
        // 1 env scan, each x axis is the KO/OE of the group listed at the end, KO at the end fixed to these values, level of KO/OE must be specified!

    longint_ARRAY *l;
    Environment_and_hits *EnvHit;
    
    l=get_requested_cell_group(PD,linenow);
    
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    unsigned long int input_id=PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
    if (input_id<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    unsigned long int input_index =PD->Coupled->BRN->active_inputs->check_for_element(input_id);
    if (input_index<1) {
        printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    
    EnvHit = extract_environment_and_hits(PD,linenow);
    if(EnvHit->p_inputs[0]<0) return; // this is set when the input is incorrent and virtual_exp line gets ignored

  
    load_Phenotypes_for_model(PD);
    printf("phenotpyes loaded!\n");//getchar();
    
    bool cc_stats,irrev_stats,rev_stats;
    Python_Grapsh_NonSaturating_Stats_FIX_KDOE_fn1Env(input_index,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
    
    string datadir,fmain;
    longint_ARRAY *lin;
    lin=new longint_ARRAY();
    lin->add_element(input_index);
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    char fn[400];
    sprintf(fn,"mkdir %s\n",fmain.c_str());
    system(fn);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s\n",datadir.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    delete lin; lin=NULL;

    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        EnvHit->p_inputs[input_index]=p_env;
        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
    }
    
    if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
        Python_Grapsh_NonSaturating_Stats_FIX_KDOE_fn1Env(input_index,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
   
        if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
            printf("Pyhton script cannot find relevant files immeditately AFTER completion of statistical sampling runs!\n");
            getchar();
        }
    }
}

void NonSaturating_Stats_Scan_1KDOE_fn1Env(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[],int run_assync){
    // 2 options!!!
     
     // NonSaturating_Stats_Scan_1KDOE_fn1Env_Code    P    GF_High  (GF=1, GF_High=1, Gamma=0) (p21_mRNA KD)
             // 1 env on x axis is ECM, scan of KO from none to full across (they go in tandem)
    
    longint_ARRAY *l;
    Environment_and_hits *EnvHit;
    
    l=get_requested_cell_group(PD,linenow);
    
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    unsigned long int input_id=PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
    if (input_id<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    unsigned long int input_index =PD->Coupled->BRN->active_inputs->check_for_element(input_id);
    if (input_index<1) {
        printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    
    EnvHit = extract_environment_and_hits(PD,linenow);
    if(EnvHit->p_inputs[0]<0) return; // this is set when the input is incorrent and virtual_exp line gets ignored

  
    load_Phenotypes_for_model(PD);
    printf("phenotpyes loaded!\n");//getchar();
    
    bool cc_stats,irrev_stats,rev_stats;
    Python_Grapsh_NonSaturating_Stats_Scan_1KDOE_fn1Env(input_index,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
    
    string datadir,fmain;
    longint_ARRAY *lin;
    lin=new longint_ARRAY();
    lin->add_element(input_index);
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    char fn[400];
    sprintf(fn,"mkdir %s\n",fmain.c_str());
    system(fn);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s\n",datadir.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    delete lin; lin=NULL;

    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        EnvHit->p_inputs[input_index]=p_env;
        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
    }
    
    if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
        Python_Grapsh_NonSaturating_Stats_Scan_1KDOE_fn1Env(input_index,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
   
        if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
            printf("Pyhton script cannot find relevant files immeditately AFTER completion of statistical sampling runs!\n");
            getchar();
        }
    }
}

void NonSaturating_Stats_Scan_1KDOE_2Env(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[],int run_assync){
     
     //   NonSaturating_Stats_Scan_1KDOE_2Env    P    GF_High  Stiff_ECM (GF=1, GF_High=1, Gamma=0)
     //   NonSaturating_Stats_Scan_1KDOE_2Env    P    GF_High  Stiff_ECM (GF=1, GF_High=1, Gamma=0) (p21_mRNA KD)
             // envs across x and y axis of grid, scan of KO from none to full across x axis of each graph
    
    longint_ARRAY *l;
    Environment_and_hits *EnvHit;
    
    l=get_requested_cell_group(PD,linenow);
    
    istringstream iss(linenow);
    string sr,inputname_1,inputname_2;
    iss >> sr; iss >> sr;
    iss >> inputname_1;
    unsigned long int input_id_1=PD->Coupled->BRN->MDAT->get_node_ID(inputname_1.c_str());
    if (input_id_1<1) {
        printf ("First input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname_1.c_str(),linenow);
        return;
    }
    iss >> inputname_2;
    unsigned long int input_id_2=PD->Coupled->BRN->MDAT->get_node_ID(inputname_2.c_str());
    if (input_id_2<1) {
        printf ("Second input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname_2.c_str(),linenow);
        return;
    }
    
    unsigned long int input_index_1 =PD->Coupled->BRN->active_inputs->check_for_element(input_id_1);
    if (input_index_1<1) {
        printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",inputname_1.c_str(),linenow);
        return;
    }

    
    unsigned long int input_index_2 =PD->Coupled->BRN->active_inputs->check_for_element(input_id_2);
    if (input_index_2<1) {
        printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",inputname_2.c_str(),linenow);
        return;
    }

    EnvHit = extract_environment_and_hits(PD,linenow);
    if(EnvHit->p_inputs[0]<0) return; // this is set when the input is incorrent and virtual_exp line gets ignored


    load_Phenotypes_for_model(PD);
    printf("phenotpyes loaded!\n");//getchar();
    
    bool cc_stats,irrev_stats,rev_stats;
    Python_Grapsh_NonSaturating_Stats_Scan_1KDOE_2Env(input_index_1,input_index_2,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
    
    string datadir,fmain;
    longint_ARRAY *lin;
    lin=new longint_ARRAY();
    lin->add_element(input_index_1);
    lin->add_element(input_index_2);
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    char fn[400];
    sprintf(fn,"mkdir %s\n",fmain.c_str());
    system(fn);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s\n",datadir.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    datadir = get_Figures_directory_now(fmain,lin,PD);
    delete lin; lin=NULL;

    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        EnvHit->p_inputs[input_index_1]=p_env_1;
        for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
            EnvHit->p_inputs[input_index_2]=p_env_2;
            if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
            if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
            if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
        }
    }
    
    if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
        Python_Grapsh_NonSaturating_Stats_Scan_1KDOE_2Env(input_index_1,input_index_2,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);

        if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
            printf("Pyhton script cannot find relevant files immeditately AFTER completion of statistical sampling runs!\n");
            getchar();
        }
    }
}

void MutantSummary_5D_333_22(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[],int run_assync){
    printf ("Nice try mutantsum\n");
    longint_ARRAY *l;
    Environment_and_hits *EnvHit;
    
    l=get_requested_cell_group(PD,linenow);
    
    istringstream iss(linenow);
    string sr,ix1,ix2,iy1,iy2,iz1,iz2,iu,iv,mutant;
    char skipc;
    iss >> sr; iss >> sr;
    iss >> mutant;
    unsigned long int mutant_ID=PD->Coupled->BRN->MDAT->get_node_ID(mutant.c_str());
    if (mutant_ID<1) {
        printf ("Mutant node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",mutant.c_str(),linenow);
        return;
    }
    
    iss >> skipc; iss>>ix1; if (!ix1.empty()) ix1.pop_back();
    unsigned long int ix1_ID=PD->Coupled->BRN->MDAT->get_node_ID(ix1.c_str());
    if (ix1_ID<1) {printf ("X coodinate's first input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",ix1.c_str(),linenow);return;}
    unsigned long int index_x1 =PD->Coupled->BRN->active_inputs->check_for_element(ix1_ID);
    if (index_x1<1) { printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",ix1.c_str(),linenow); return;}

    
    iss>>ix2; if (!ix2.empty()) ix2.pop_back();
    unsigned long int ix2_ID=PD->Coupled->BRN->MDAT->get_node_ID(ix2.c_str());
    if (ix2_ID<1) {printf ("X coodinate's second input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",ix2.c_str(),linenow);return;}
    unsigned long int index_x2 =PD->Coupled->BRN->active_inputs->check_for_element(ix2_ID);
    if (index_x2<1) { printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",ix2.c_str(),linenow); return;}

    
    iss >> skipc; iss>>iy1; if (!iy1.empty()) iy1.pop_back();
    unsigned long int iy1_ID=PD->Coupled->BRN->MDAT->get_node_ID(iy1.c_str());
    if (iy1_ID<1) {
        printf ("Y coodinate's first input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",iy1.c_str(),linenow);
        return;
    }
    unsigned long int index_y1 =PD->Coupled->BRN->active_inputs->check_for_element(iy1_ID);
    if (index_y1<1) { printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",iy1.c_str(),linenow); return;}
   
    iss>>iy2; if (!iy2.empty()) iy2.pop_back();
    unsigned long int iy2_ID=PD->Coupled->BRN->MDAT->get_node_ID(iy2.c_str());
    if (iy2_ID<1) {
        printf ("Y coodinate's second input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",iy2.c_str(),linenow);
        return;
    }
    unsigned long int index_y2 =PD->Coupled->BRN->active_inputs->check_for_element(iy2_ID);
    if (index_y2<1) { printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",iy2.c_str(),linenow); return;}

    iss >> skipc; iss>>iz1; if (!iz1.empty()) iz1.pop_back();
    unsigned long int iz1_ID=PD->Coupled->BRN->MDAT->get_node_ID(iz1.c_str());
    if (iz1_ID<1) {
        printf ("Z coodinate's first input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",iz1.c_str(),linenow);
        return;
    }
    unsigned long int index_z1 =PD->Coupled->BRN->active_inputs->check_for_element(iz1_ID);
    if (index_z1<1) { printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",iz1.c_str(),linenow); return;}

    iss>>iz2; if (!iz2.empty()) iz2.pop_back();
    unsigned long int iz2_ID=PD->Coupled->BRN->MDAT->get_node_ID(iz2.c_str());
    if (iz2_ID<1) {
        printf ("Z coodinate's second input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",iz2.c_str(),linenow);
        return;
    }
    unsigned long int index_z2 =PD->Coupled->BRN->active_inputs->check_for_element(iz2_ID);
    if (index_z2<1) { printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",iz2.c_str(),linenow); return;}

    iss >> iu;
    unsigned long int iu_ID=PD->Coupled->BRN->MDAT->get_node_ID(iu.c_str());
    if (iu_ID<1) {
        printf ("4th coodinate's second input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",iu.c_str(),linenow);
        return;
    }
    unsigned long int index_u =PD->Coupled->BRN->active_inputs->check_for_element(iu_ID);
    if (index_u<1) { printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",iu.c_str(),linenow); return;}

    iss >> iv;
    unsigned long int iv_ID=PD->Coupled->BRN->MDAT->get_node_ID(mutant.c_str());
    if (iv_ID<1) {
        printf ("5th coodinate's second input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",iv.c_str(),linenow);
        return;
    }
    unsigned long int index_v =PD->Coupled->BRN->active_inputs->check_for_element(iv_ID);
    if (index_v<1) { printf ("Input node %s from virtual experiment line %s is a node, but not an INPUT of this Boolean model! (ignoring line)\n",iv.c_str(),linenow); return;}

    string rem(iss.str().substr(iss.tellg()));
    EnvHit = extract_environment_and_hits(PD,rem.c_str());
    if(EnvHit->p_inputs[0]<0) return; // this is set when the input is incorrent and virtual_exp line gets ignored


    load_Phenotypes_for_model(PD);
    printf("phenotpyes loaded!\n");//getchar();
    
    bool cc_stats=0,irrev_stats=0,rev_stats=0;
    unsigned long inputs_5D[8];
    inputs_5D[0]=ix1_ID;inputs_5D[1]=ix2_ID;
    inputs_5D[2]=iy1_ID;inputs_5D[2]=iy2_ID;
    inputs_5D[4]=iz1_ID;inputs_5D[5]=iz2_ID;
    inputs_5D[6]=iu_ID;inputs_5D[7]=iv_ID;

    //  Python_Graphs_MutantSummary_5D_333_22(inputs_5D,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);
    
    string fmain;
    char datadir [400];
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    char fn[400];
    sprintf(fn,"mkdir %s\n",fmain.c_str());
    system(fn);
    sprintf(datadir,"%s/MutSum_%s__%s-%s__%s-%s__%s-%s__%s__%s",fmain.c_str(),mutant.c_str(),ix1.c_str(),ix2.c_str(),iy1.c_str(),iy2.c_str(),iz1.c_str(),iz2.c_str(),iu.c_str(),iv.c_str());
    sprintf(fn,"mkdir %s\n",datadir);
    system(fn);
    
    double eint[5];
    eint[0]=0.05; eint[1]=0.25; eint[2]=0.5; eint[3]=0.75; eint[4]=0.95;
    for(unsigned int x=0;x<5; x++){
        EnvHit->p_inputs[index_x1]=eint[x];
        EnvHit->p_inputs[index_x2]=0;
        
        for(unsigned int y=0;y<5; y++){
            EnvHit->p_inputs[index_y1]=eint[y];
            EnvHit->p_inputs[index_y2]=0;

            for(unsigned int z=0;z<5; z++){
                EnvHit->p_inputs[index_z1]=eint[z];
                EnvHit->p_inputs[index_z2]=0;

                for(unsigned int u=0;u<5; u++){
                    EnvHit->p_inputs[index_u]=eint[u];
                    
                    for(unsigned int v=0;v<5; v++){
                        EnvHit->p_inputs[index_v]=eint[v];
                        
                        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                    }
                }
                printf("Finished x: %.2lf y: %.2lf z: %.2lf\n",eint[x],eint[y],eint[z]);
            }
            for(unsigned int z=0;z<5; z++){
                EnvHit->p_inputs[index_z1]=1;
                EnvHit->p_inputs[index_z2]=eint[z];

                for(unsigned int u=0;u<5; u++){
                    EnvHit->p_inputs[index_u]=eint[u];
                    
                    for(unsigned int v=0;v<5; v++){
                        EnvHit->p_inputs[index_v]=eint[v];
                        
                        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                    }
                }
                printf("Finished x: %.2lf y: %.2lf z: 1-> %.2lf\n",eint[x],eint[y],eint[z]);
            }
        }
        for(unsigned int y=0;y<5; y++){
            EnvHit->p_inputs[index_y1]=1;
            EnvHit->p_inputs[index_y2]=eint[y];

            for(unsigned int z=0;z<5; z++){
                EnvHit->p_inputs[index_z1]=eint[z];
                EnvHit->p_inputs[index_z2]=0;

                for(unsigned int u=0;u<5; u++){
                    EnvHit->p_inputs[index_u]=eint[u];
                    
                    for(unsigned int v=0;v<5; v++){
                        EnvHit->p_inputs[index_v]=eint[v];
                        
                        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                    }
                }
                printf("Finished x: %.2lf y: 1-> %.2lf z: %.2lf\n",eint[x],eint[y],eint[z]);
            }
            for(unsigned int z=0;z<5; z++){
                EnvHit->p_inputs[index_z1]=1;
                EnvHit->p_inputs[index_z2]=eint[z];

                for(unsigned int u=0;u<5; u++){
                    EnvHit->p_inputs[index_u]=eint[u];
                    
                    for(unsigned int v=0;v<5; v++){
                        EnvHit->p_inputs[index_v]=eint[v];
                        
                        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                    }
                }
                printf("Finished x: %.2lf y: 1-> %.2lf z: 1-> %.2lf\n",eint[x],eint[y],eint[z]);
            }
        }
    }
    for(unsigned int x=0;x<5; x++){
        EnvHit->p_inputs[index_x1]=1;
        EnvHit->p_inputs[index_x2]=eint[x];
        
        for(unsigned int y=0;y<5; y++){
            EnvHit->p_inputs[index_y1]=eint[y];
            EnvHit->p_inputs[index_y2]=0;

            for(unsigned int z=0;z<5; z++){
                EnvHit->p_inputs[index_z1]=eint[z];
                EnvHit->p_inputs[index_z2]=0;

                for(unsigned int u=0;u<5; u++){
                    EnvHit->p_inputs[index_u]=eint[u];
                    
                    for(unsigned int v=0;v<5; v++){
                        EnvHit->p_inputs[index_v]=eint[v];
                        
                        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                    }
                }
                printf("Finished x: 1-> %.2lf y: %.2lf z: %.2lf\n",eint[x],eint[y],eint[z]);
            }
            for(unsigned int z=0;z<5; z++){
                EnvHit->p_inputs[index_z1]=1;
                EnvHit->p_inputs[index_z2]=eint[z];

                for(unsigned int u=0;u<5; u++){
                    EnvHit->p_inputs[index_u]=eint[u];
                    
                    for(unsigned int v=0;v<5; v++){
                        EnvHit->p_inputs[index_v]=eint[v];
                        
                        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                    }
                }
                printf("Finished x: 1-> %.2lf y: %.2lf z: 1-> %.2lf\n",eint[x],eint[y],eint[z]);
            }
        }
        for(unsigned int y=0;y<5; y++){
            EnvHit->p_inputs[index_y1]=1;
            EnvHit->p_inputs[index_y2]=eint[y];

            for(unsigned int z=0;z<5; z++){
                EnvHit->p_inputs[index_z1]=eint[z];
                EnvHit->p_inputs[index_z2]=0;

                for(unsigned int u=0;u<5; u++){
                    EnvHit->p_inputs[index_u]=eint[u];
                    
                    for(unsigned int v=0;v<5; v++){
                        EnvHit->p_inputs[index_v]=eint[v];
                        
                        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                    }
                }
                printf("Finished x: 1-> %.2lf y: 1-> %.2lf z: %.2lf\n",eint[x],eint[y],eint[z]);
            }
            for(unsigned int z=0;z<5; z++){
                EnvHit->p_inputs[index_z1]=1;
                EnvHit->p_inputs[index_z2]=eint[z];

                for(unsigned int u=0;u<5; u++){
                    EnvHit->p_inputs[index_u]=eint[u];
                    
                    for(unsigned int v=0;v<5; v++){
                        EnvHit->p_inputs[index_v]=eint[v];
                        
                        if(cc_stats!=Phen_STATs.CellCycleErrors) CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(irrev_stats!=Phen_STATs.Irrev) Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                        if(rev_stats!=Phen_STATs.Rev) RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(datadir,EnvHit,l,PD, run_assync);
                    }
                }
                printf("Finished x: 1-> %.2lf y: 1-> %.2lf z: 1-> %.2lf\n",eint[x],eint[y],eint[z]);
            }
        }
    }
    if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
       // Python_Grapsh_NonSaturating_Stats_Scan_1KDOE_2Env(input_index_1,input_index_2,EnvHit,l,PD, run_assync,&cc_stats,&irrev_stats,&rev_stats);

        if(((cc_stats!=Phen_STATs.CellCycleErrors))||((irrev_stats!=Phen_STATs.Irrev))||((rev_stats!=Phen_STATs.Rev))) {
            printf("Pyhton script cannot find relevant files immeditately AFTER completion of statistical sampling runs!\n");
            getchar();
        }
    }
}
 
void run_Timing_Cycles_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    printf ("Nice try 7\n");
}
void run_Timing_Transition_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    printf ("Nice try 9\n");
}
void run_Async_Draw_SyncCycles_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    printf ("Nice try 13\n");
}
void run_Async_Pulse1_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l;
    unsigned int inputL;
    
    l=get_requested_cell_group(PD,linenow);
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    if (PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str())<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    iss >> inputL;
    if (iss.fail() == true){
        printf("Invalid input length in line %s\n",linenow); getchar(); l = NULL;
    }
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_OneNode_pulse_experiment_ASSYNC(inputname.c_str(),(int)(l->A[j]),1,PD,0.000000001,inputL);
        }
    }
    delete l;l=NULL;
}

void run_ModelErrors_Pulse1_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l;
    unsigned int inputL,nr_errors;
    
    l=get_requested_cell_group(PD,linenow);
    istringstream iss(linenow);
    string sr,inputname;
    iss >> sr; iss >> sr;
    iss >> inputname;
    if (PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str())<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    iss >> inputL;
    if (iss.fail() == true){
        printf("Invalid input length in line %s\n",linenow); getchar(); l = NULL;
    }
    
    iss >> nr_errors;
    if (iss.fail() == true){
        printf("Invalid number of random link removals / gate flips in line %s\n",linenow); getchar(); l = NULL;
    }
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            run_OneNode_pulse_experiment_ModelErrors(inputname.c_str(),(int)(l->A[j]),1,PD,inputL,nr_errors);
        }
    }
    delete l;l=NULL;
}


void run_Async_PulseM_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    printf ("Nice try 14\n");
}

void run_Async_Timing_Cycles_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    printf ("Nice try 18\n");
}
void run_Async_Timing_Transition_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    printf ("Nice try 20\n");
}

void run_ModulePulse1_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    unsigned int inputL;
    string sr,modulename,inputname;
    
    istringstream iss(linenow);
    iss >> sr;
    iss >> modulename;
    iss >> inputname;
    if (PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str())<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    if (PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str())<1) {
        printf ("Input node %s from virtual experiment line %s is not in Boolean model! (ignoring line)\n",inputname.c_str(),linenow);
        return;
    }
    iss >> inputL;
    if (iss.fail() == true){
        printf("Invalid input length in line %s\n",linenow); getchar();
    }
    
    for(unsigned int mID=1;mID<=PD->BRN->MDAT->Module_NR;mID++){
        if(modulename == PD->BRN->Module_Names[mID]){
            if( PD->BRN_Module[mID]->N<=20){
                Boolean_Dynamics *D;
                D=new Boolean_Dynamics(PD->BRN_Module[mID]);
                D->generate_full_state_graph_synchronously(0);
                for (unsigned int k=1; k<=D->Attractor_NR; k++) {
                    run_OneNode_pulse_experiment_MOD(inputname.c_str(),k,1,D,0.000000001,inputL);
                }
                delete D;D=NULL;
            }
            else {
                Boolean_Dynamics_TRACKER *D2;
                D2=new Boolean_Dynamics_TRACKER(PD->BRN_Module[mID],500);
                D2->record_timetraces_from_collection_of_random_inputs(100,100,0.01);
                for (unsigned int k=1; k<=D2->Attractor_NR; k++) {
                    run_OneNode_pulse_experiment_MOD(inputname.c_str(),k,1,D2,0.000000001,inputL);
                }
                delete D2;D2=NULL;
            }
        
        }
        
    }
    
    printf ("Nice try 24\n");
}

void run_ModulePulseM_Code(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    printf ("Nice try 25\n");
}

void run_CycleDrivingNetwork(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    longint_ARRAY *l, cg;
    char fnd[600];
    l=get_requested_cell_group(PD,linenow);
    
    sprintf(fnd,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fnd);
    sprintf(fnd,"mkdir %s/%s/_EXP/CyclingAttractor_SCCs",MODEL_Directory_system,PD->BRN->path);
    system(fnd);
    sprintf(fnd,"mkdir %s/%s/_EXP/CyclingAttractor_SCCs",MODEL_Directory_system,PD->BRN->path);
    system(fnd);
 

    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            if (PD->Coupled->Attractor_valleys[l->A[j]]->N_a>1){
                // find nodes that cycle
                short int *stable;
                stable = new short int[PD->BRN->N+1];
                for(int i=1;i<=PD->BRN->N;i++)
                    stable[i] = PD->Coupled->Attractor_valleys[l->A[j]]->Basin_Stored[PD->Coupled->Attractor_valleys[l->A[j]]->Attractor[1]]->s[i-1];
                for(int a_st_k=1;a_st_k<=PD->Coupled->Attractor_valleys[l->A[j]]->N_a;a_st_k++)
                    for(int i=1;i<=PD->BRN->N;i++)
                        if (stable[i] != PD->Coupled->Attractor_valleys[l->A[j]]->Basin_Stored[PD->Coupled->Attractor_valleys[l->A[j]]->Attractor[a_st_k]]->s[i-1])
                            stable[i]  = -1;
               
                // make graph of nodes that cycle
                longint_ARRAY *cg; cg=new longint_ARRAY();
                for(int i=1;i<=PD->BRN->N;i++)
                    if (stable[i]==-1) {
                        cg->add_element(i);
                     //   printf("cycling nr %ld, id %ld, name %s\n",cg->N,i,PD->BRN->MDAT->Node_name[i]);
                    }
                
                sgraph *net,*gcc,*gcc_broken; char nev[400];
                sprintf(nev,"%s_A-%ld_cycling_node_SSC",PD->BRN->path,l->A[j]);
                net=new sgraph(0,nev);
             //   printf("before nodes, node nr %ld\n",net->N);getchar();
                unsigned long int id1;
                for(int i=1;i<=cg->N;i++){
                    id1 = net->add_node(PD->BRN->MDAT->Node_name[cg->A[i]]);
                    net->node[i].ID = cg->A[i];
                //    printf("cg node %ld, ID %ld, name %s\n",i,net->node[id1].ID,net->node[id1].node_draw_props->node_label);
                }
              //  getchar();
             //   printf("before links, node nr %ld\n",net->N);getchar();
                for(int i=1;i<=cg->N;i++)
                    for(int k=1;k<=cg->N;k++)
                        if((PD->BRN->BN->node[cg->A[i]].Li.get_dirlink(cg->A[k])) && (PD->Coupled->cascading_change(l->A[j],cg->A[k],cg->A[i])))
                            net->add_link(i,k,PD->BRN->BN->get_dir_link(cg->A[i],cg->A[k]));
                
                //printf("links added, node nr %ld\n",net->N);getchar();
                
                gcc=get_Strongly_Connected_Components_linked(net);
                PD->Coupled->export_subgraph_to_GML(gcc, "CyclingAttractor_SCCs");
                gcc_broken=get_Strongly_Connected_Components_Disconneceted(net);
                PD->Coupled->export_subgraph_to_GML(gcc_broken, "CyclingAttractor_SCCs");
              //  write_network_draw_attr(gcc,nev);
            }
        }
    }
    delete l;l=NULL;
}



#endif /* Run_Virtual_Exp_Options_h_h */
