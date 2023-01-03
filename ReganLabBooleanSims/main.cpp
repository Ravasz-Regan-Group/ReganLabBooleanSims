//
//  main.cpp
//  ReganLabBooleanSims
//
//  Created by Erzsó Ravasz Regan on 9/6/16.
//  Copyright © 2016 Regan_Group. All rights reserved.
//
#define whoamI  "eregan"
#define CODE_Directory "/Users/eregan/Dropbox/_CODE/Boolean code_for_IS"
#define MODEL_Directory "/Users/eregan/Dropbox/__RESEARCH/MODEL_Folders/BooleanModels_dmms"
#define MODEL_Directory_system "/Users/eregan/Dropbox/__RESEARCH/MODEL_Folders/BooleanModels_dmms"

#define DEFAULT_DELIMITER 9


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <time.h>
#include <limits.h>

char Model_type[150] = "Life_and_Death";

char GF_low_name[150],GF_high_name[150];
char CD_low_name[150] = "CellDensity_Low",CD_high_name[150] = "CellDensity_High";

using std::cout;
using std::endl;

struct equal_string_for_google_hash
{   bool operator()(const char* s1, const char* s2) const
    { return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
    }
};


#include "rand_fg.h"
#include "defs.h"
#include "base_convert.h"
#include "Epslib_withoutX.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include <gsl/gsl_cblas.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_bspline.h>

#include "logar.h"
#include "NetBasics.h"
#include "Net_IO.h"


//#include "B_spline_enthropy.h"
//#include "Hiearachical_tree_on_matrix.h"
//#include "studentttests.h"
//#include "K_means_datastructure.h"
//#include "nrutil.c"

//extern "C" {
//#include"fisher2.h"
//}

//#define ASSYNC_BIAS 0
//#define ASSYNC_BIAS 1 // fast vs. slow for all TFs and 4N DNA and cytokinesis
#define ASSYNC_BIAS 30
/// PI3K paper is going in with update rule 30


#include "Boolean_modeling_central.h"
    //#include "Dynamical_Modularity.h"
    //#include "Boolean_modeling_visualize_landscape.h"

#include "Boolean_Model_sampling_2.h"
#include "Dynamical_Mod_Paired_Sampling.h"
#include "Asynchronous_version.h"

#include "Model_specific_info.h"
#include "Modules.h"
#include "Analyze_Boolean_Model.h"
#include "Verify_gates.h"
#include "Compare_Boolean_Models.h"
#include "Run_Virtual_Exp_Options.h"

void run_small_model(const char mname[]){
        Boolean_Dynamics *D;
        Boolean_RegNet *A;
        //  snode::slinklist::slink *sc2;
        double p_error=0.02;
    
        double beta=0.5 *log(1./p_error-1.);
    
        A=new Boolean_RegNet(mname);
    
        D=new Boolean_Dynamics(A);
        A->export_Boolean_RegNet();
        A->export_Boolean_RegNet_with_names_to_Cytoskape();
        
        longint_ARRAY *SCC;
        SCC=A->Find_Strongly_Connected_Component();
        for(unsigned int i=1;i<=A->N;i++){
            if (SCC->check_for_element(i)==0) {
                printf("Node %s is not in the module's SCC\n",A->MDAT->Node_name[i]);
            }
        }
        //   getchar();
        
        D->generate_full_state_graph_synchronously(0);
        
        D->calculate_noisy_landscape(beta);
        
        D->export_state_Graph_to_Cytoskape_no_noise();
        D->export_state_Graph_to_Cytoskape();
        //export_Module_Attractors_onto_Boolean_network__PDF(mID, PD,D);
        delete D;D=NULL;
        exit(1);
}

#define Comment_Code 0
#define Modules_Code 1
#define Async_Module_Cycles_Code 2
#define Pulse1_Code 3
#define PulseM_Code 12
#define PulseScan_Code 4
#define NonSaturating_Draw_Code 5
#define NonSaturating_Stats_FIX_KDOE_fn1Env_Code 6
#define NonSaturating_Stats_Scan_1KDOE_fn1Env_Code 7
#define NonSaturating_Stats_Scan_2KDOE_fn1Env_Code 8
#define NonSaturating_Stats_FIX_Env_fnKDOE_Code 9
#define NonSaturating_Stats_Scan_1Env_fnKDOE_Code 10
#define NonSaturating_Stats_Scan_2Env_fnKDOE_Code 11
#define Async_PulseM_Code    13
#define Async_Draw_SyncCycles_Code   14
#define Async_PulseScan_Code 15
#define Async_Pulse1_Code    16


#define Async_NonSaturating_Draw_Code    17

#define Timing_Cycles_Code  19
#define Timing_Transition_Code 20
#define Async_Timing_Cycles_Code  21
#define Async_Timing_Transition_Code 22

#define ModulePulse1_Code    24
#define ModulePulseM_Code    25

#define Async_NonSaturating_Stats_FIX_KDOE_fn1Env_Code 36
#define Async_NonSaturating_Stats_Scan_1KDOE_fn1Env_Code 37
#define Async_NonSaturating_Stats_Scan_2KDOE_fn1Env_Code 38
#define Async_NonSaturating_Stats_FIX_Env_fnKDOE_Code 39
#define Async_NonSaturating_Stats_Scan_1Env_fnKDOE_Code 40
#define Async_NonSaturating_Stats_Scan_2Env_fnKDOE_Code 41
#define AsyncBias_NonSaturating_Stats_FIX_KDOE_fn1Env_Code 42
#define AsyncBias_NonSaturating_Stats_Scan_1KDOE_fn1Env_Code 43
#define AsyncBias_NonSaturating_Stats_Scan_2KDOE_fn1Env_Code 44
#define AsyncBias_NonSaturating_Stats_FIX_Env_fnKDOE_Code 45
#define AsyncBias_NonSaturating_Stats_Scan_1Env_fnKDOE_Code 46
#define AsyncBias_NonSaturating_Stats_Scan_2Env_fnKDOE_Code 47

#define Pulse1_Cycle 50
#define Pulse2_Cycles 51

#define NonSaturating_Stats_Scan_1KDOE_2Env_Code 52
#define MutantSummary_5D_333_22_Code 53

#define CycleDrivingNetwork 54

#define ModelErrors_Pulse1 103

unsigned int get_Exp_Code (string w){
    if (w == "--")  return 0;
    if (w == "Modules")  return Modules_Code;
    if (w == "Async_Module_Cycles")  return Async_Module_Cycles_Code;
    if (w == "Pulse1")  return Pulse1_Code;
    if (w == "PulseM")  return PulseM_Code;
    if (w == "PulseScan")  return PulseScan_Code;
    if (w == "NonSaturating_Draw")  return NonSaturating_Draw_Code;
    
    if (w == "NonSaturating_Stats_FIX_KDOE_fn1Env")  return NonSaturating_Stats_FIX_KDOE_fn1Env_Code;
    if (w == "Async_NonSaturating_Stats_FIX_KDOE_fn1Env")  return Async_NonSaturating_Stats_FIX_KDOE_fn1Env_Code;
    if (w == "AsyncBias_NonSaturating_Stats_FIX_KDOE_fn1Env")  return AsyncBias_NonSaturating_Stats_FIX_KDOE_fn1Env_Code;
   
    if (w == "NonSaturating_Stats_Scan_1Env_fnKDOE")  return NonSaturating_Stats_Scan_1Env_fnKDOE_Code;
    if (w == "NonSaturating_Stats_Scan_1KDOE_2Env")  return NonSaturating_Stats_Scan_1KDOE_2Env_Code;
    if (w == "MutantSummary_5D_333_22")  return MutantSummary_5D_333_22_Code;
    
    if (w == "Async_NonSaturating_Stats_Scan_1Env_fnKDOE")  return Async_NonSaturating_Stats_Scan_1Env_fnKDOE_Code;
    if (w == "AsyncBias_NonSaturating_Stats_Scan_1Env_fnKDOE")  return AsyncBias_NonSaturating_Stats_Scan_1Env_fnKDOE_Code;

    if (w == "NonSaturating_Stats_Scan_1KDOE_fn1Env")  return NonSaturating_Stats_Scan_1KDOE_fn1Env_Code;
    if (w == "Async_NonSaturating_Stats_Scan_1KDOE_fn1Env")  return Async_NonSaturating_Stats_Scan_1KDOE_fn1Env_Code;
    if (w == "AsyncBias_NonSaturating_Stats_Scan_1KDOE_fn1Env")  return AsyncBias_NonSaturating_Stats_Scan_1KDOE_fn1Env_Code;
    
    if (w == "Timing_Cycles")  return Timing_Cycles_Code;
    if (w == "Timing_Transition")  return Timing_Transition_Code;
    if (w == "Async_Timing_Cycles")  return Async_Timing_Cycles_Code;
    if (w == "Async_Timing_Transition")  return Async_Timing_Transition_Code;

    if (w == "Pulse1_Cycle")  return Pulse1_Cycle;
    if (w == "Pulse2_Cycles")  return Pulse2_Cycles;

    if (w == "Async_Draw_SyncCycles")  return Async_Draw_SyncCycles_Code;
    if (w == "Async_PulseScan")  return Async_PulseScan_Code;
    if (w == "Async_Pulse1")  return Async_Pulse1_Code;
    if (w == "Async_PulseM")  return Async_PulseM_Code;
    if (w == "Async_NonSaturating_Draw")  return Async_NonSaturating_Draw_Code;
    if (w == "Async_NonSaturating_Stats_FIX_KDOE_fn1Env")  return Async_NonSaturating_Stats_FIX_KDOE_fn1Env_Code;
    if (w == "ModulePulse1")  return ModulePulse1_Code;
    if (w == "ModulePulseM")  return ModulePulseM_Code;
    
    if (w == "CycleDrivingNetwork")  return CycleDrivingNetwork;
    
    if (w == "ModelErrors_Pulse1")  return ModelErrors_Pulse1;
   
    return 0;
}

void run_InSilico_Experimetnal_Package(const char  ModelName[], const char VexpName[]){
    char fn[300];
    FILE *f;
    Boolean_RegNet *A;
    unsigned int N_trace,N_RND,N_MAX;
    double p_error;
    MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD = NULL;
    char attr_list[200][5000],attr_list2[200][5000];
    
    N_trace=5;N_RND=100;
    p_error=0.02;N_MAX=1000;
    
    sprintf(fn,"%s/%s/%s",MODEL_Directory, ModelName,VexpName);
    f=fopen(fn,"r");
    if(f==NULL) {printf("Virtual Experiment File %s not found in Model directory %s \n \tC urrent model folder: %s\n%s\n",VexpName,ModelName,MODEL_Directory_system,fn); exit(1);}
    
    // first line: set up model
    if(fgets(GSE_line,MAX_LINE,f)!=NULL){
        istringstream iss(GSE_line);
        string mname,gf,gfh; iss >> mname;
        if (mname != ModelName) {
            printf("Virtual Experiment File must start with the correct model name %s, followed by an optional list of the two nodes corresponding to low and high growth factor conditions\n",ModelName);
            exit(1);
        }
        iss >> gf;
        if (iss.fail() != true) {
            iss >> gfh;
            if (iss.fail() != true) {
                strcpy(GF_low_name,gf.c_str());strcpy(GF_high_name,gfh.c_str());
            }
            else { //  3-level GF inputs need to come in pairs
                printf("Please specify two node names for LOW and HIGH growth factor conditions, if applicable (if GF levels are binary, the first line should only contain the model name).\n");
                exit(1);
            }
        }
        else { strcpy(GF_low_name,"");strcpy(GF_high_name,"");
            // there are no 3-level GF inputs; read model and determine all inputs from file, as well as all silenced inputs
        }
        
        A=new Boolean_RegNet(ModelName);
        A->load_Module_IDs_and_Drawing_Order();
        A->silenced_inputs = new longint_ARRAY();
        
        if(fgets(GSE_line,MAX_LINE,f)!=NULL){
            istringstream iss(GSE_line);
            string sr; iss >> sr;
            if (sr == "Sampling") {
                iss >> N_RND;
                if (iss.fail() == true) {
                    printf("Please the number of random initial conditions to be sampled for each unique environment as a valid integer on line 2.\n");
                    exit(1);
                }
                iss >> sr;
                if (sr[0] == '('){
                    std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(GSE_line, '(',')');
                    if (word.size() > 1) {
                        printf("Line %s not correctly formatted; code needs a single list of input node names that should be kept locked into a single state throughout the sampling. They need to be enclosed in ( ... ) brackets with no space on either side before / after the first /last nodename. Nodes sgould be separatd by comma and space: e.g.: (Trail=0, UV=0, ROS_ext=0) \n",GSE_line);
                        getchar(); exit (1);
                    }
                    char bla2[600]; sprintf(bla2,"%s",word[0].c_str());
                    unsigned long int Nf=separate_line_by_string(bla2,attr_list,", ",200);
                    string sr2;
                    if (Nf > 0){
                        for(unsigned int i=1;i<=Nf;i++){
                            unsigned long int Nf2=separate_line_by_string(attr_list[i],attr_list2,"=",10);
                            if (Nf2 > 0){
                                istringstream iss_in(attr_list2[1]);
                                string inputname;
                                iss_in >> inputname;
                                unsigned int in_id = A->MDAT->get_node_ID(inputname.c_str());
                                if (in_id>1){
                                    istringstream iss_p(attr_list2[2]);
                                    unsigned int v = -1;
                                    iss_p >> v;
                                    if ((iss_p.fail() == true) || (v<0) || (v>1)){
                                        printf("Environmental setting %s for node %s in is not a valid input of  0 and 1!\n",attr_list2[2], inputname.c_str()); getchar(); exit(1);
                                    }
                                    else {
                                        A->silenced_inputs->add_element(in_id);
                                        for(unsigned long long kk=0;kk<A->Gate[in_id]->N;kk++)
                                            A->Gate[in_id]->O[kk]=v;
                                    }
                                }
                                else {
                                    printf("Input node %s not in this model, please check list of inputs you wish to lock in the Sampling options of the %s file \n",inputname.c_str(),VexpName);
                                    getchar(); exit(1);
                                }
                            }
                            else {
                                printf("Line %s not correctly formatted; code needs a single list of input node names that should be kept locked into a single state throughout the sampling. They need to be enclosed in ( ... ) brackets with no space on either side before / after the first /last nodename. Nodes sgould be separatd by comma and space: e.g.: (Trail=0, UV=0, ROS_ext=0) \n",GSE_line);
                                getchar(); exit(1);
                            }
                        }
                    }
                }
            }
            else {fclose(f); f=fopen(fn,"r"); fgets(GSE_line,MAX_LINE,f); // re-read first line to reset buffer
            }
        }
        
        A->identify_input_nodes();
        printf("GF_low = %s\t GF_High = %s\n",GF_low_name,GF_high_name);
        A->input_node_report();
        //getchar();
       
        PD=new MODULAR_Boolean_Dynamics_PAIRED_TRACKER(A,N_MAX);
        PD->BRN->export_Boolean_RegNet_with_names_to_Cytoskape();
       // PD->read_Node_Cytoskape_positions();
        PD->read_Node_positions_from_DMMS_gDetails();
        //PD->Coupled->export_GML();
        
       
        PD->Coupled->record_timetraces_from_collection_of_random_inputs_per_inputBOX(N_trace,N_RND,p_error,A->active_inputs);
        
        for(int i=1;i<=PD->Coupled->Attractor_NR;i++){
            write_link_for_All_environmental_changes_from_ONE_attractor(A->active_inputs,i,PD);
        }
        PD->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,A->active_inputs->N)),p_error);
        Export_Attractor_Life_and_Death_table(PD,A->active_inputs,1);
      //  PD->Coupled->export_Attractor_onto_Boolean_network__PDF(1);
        
        longint_ARRAY *active_inputs_minusCD;
        active_inputs_minusCD=new longint_ARRAY();
        active_inputs_minusCD->add_longint_ARRAY(A->active_inputs);
        
        longint_ARRAY *active_inputs_noGF_minusCD;
        active_inputs_noGF_minusCD=new longint_ARRAY();
        active_inputs_noGF_minusCD->add_longint_ARRAY(A->active_inputs_noGF);
        
        unsigned int CDL_id = A->MDAT->get_node_ID(CD_low_name);
        if(CDL_id>0) {
            active_inputs_minusCD->delete_value(CDL_id);
            active_inputs_noGF_minusCD->delete_value(CDL_id);
        }
        printf("running 3D map with:\n");
        printf("active_inputs_minusCD input names:\t");
        for (int i=1; i<=active_inputs_minusCD->N; i++)
            printf("%s\t",A->MDAT->Node_name[active_inputs_minusCD->A[i]]);
        printf("\n");
        printf("active_inputs_noGF_minusCD input names:\t");
        for (int i=1; i<=active_inputs_noGF_minusCD->N; i++)
            printf("%s\t",A->MDAT->Node_name[active_inputs_noGF_minusCD->A[i]]);
        printf("\n");
       // getchar();
        
        if(active_inputs_minusCD->N > active_inputs_noGF_minusCD ->N) {
            for (int i=1; i<=active_inputs_minusCD->N; i++) {
                make_3D_map_of_steady_states_new(PD,active_inputs_noGF_minusCD,active_inputs_minusCD,i,0);
            }
//            run_OnePulse_analysis(active_inputs_minusCD,PD);
//            longint_ARRAY *l_sim;
//            l_sim= Find_Similar_Attractors__Life_and_Death(PD,A->active_inputs,1);
//            delete l_sim;
        }
        PD->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,A->active_inputs->N)),p_error);
        //PD->Coupled->Export_ALL_Attractors_onto_Boolean_Network(1);
        PD->Coupled->BRN->print_gates_to_File();
    }
    

    if (PD == NULL) exit(1);

    set_up_module_MetaData(PD);

    while (fgets(GSE_line,MAX_LINE,f)!=NULL) {
        char tagnow [300]; strcpy(tagnow,"");
        sscanf(GSE_line, "%s",tagnow);
      //  printf("GSE_line = %s\n",GSE_line); getchar();
        if (strlen(tagnow) > 0) {
            switch (get_Exp_Code(tagnow)) {
                case Modules_Code: run_Modules(PD,GSE_line); break;
                case Async_Module_Cycles_Code: run_Async_Module_Cycles_Code(PD,GSE_line); break;
                case Pulse1_Code: run_Pulse1_Code(PD,GSE_line); break;
                case PulseM_Code: run_PulseM_Code(PD,GSE_line); break;
                case PulseScan_Code: run_PulseScan_Code(PD,GSE_line); break;
                case NonSaturating_Draw_Code: run_NonSaturating_Draw_Code(PD,GSE_line); break;
                case NonSaturating_Stats_FIX_KDOE_fn1Env_Code: NonSaturating_Stats_FIX_KDOE_fn1Env(PD,GSE_line,-1); break;
                case Async_NonSaturating_Stats_FIX_KDOE_fn1Env_Code: NonSaturating_Stats_FIX_KDOE_fn1Env(PD,GSE_line,0); break;
                case AsyncBias_NonSaturating_Stats_FIX_KDOE_fn1Env_Code: NonSaturating_Stats_FIX_KDOE_fn1Env(PD,GSE_line,ASSYNC_BIAS); break;
               
                case NonSaturating_Stats_Scan_1Env_fnKDOE_Code: NonSaturating_Stats_FIX_KDOE_fn1Env(PD,GSE_line,-1); break;
                case Async_NonSaturating_Stats_Scan_1Env_fnKDOE_Code: NonSaturating_Stats_FIX_KDOE_fn1Env(PD,GSE_line,0); break;
                case AsyncBias_NonSaturating_Stats_Scan_1Env_fnKDOE_Code: NonSaturating_Stats_FIX_KDOE_fn1Env(PD,GSE_line,ASSYNC_BIAS); break;
              
                case NonSaturating_Stats_Scan_1KDOE_fn1Env_Code: NonSaturating_Stats_Scan_1KDOE_fn1Env(PD,GSE_line,-1); break;
                case NonSaturating_Stats_Scan_1KDOE_2Env_Code: NonSaturating_Stats_Scan_1KDOE_2Env(PD,GSE_line,-1); break;
               
                case MutantSummary_5D_333_22_Code: MutantSummary_5D_333_22(PD,GSE_line,-1); break;
               
                case Pulse1_Cycle: run_Pulse1_Cycle(PD,GSE_line); break;
                case Pulse2_Cycles: run_Pulse2_Cycles(PD,GSE_line); break;
                case CycleDrivingNetwork: run_CycleDrivingNetwork(PD,GSE_line); break;
               
                case Async_Draw_SyncCycles_Code: run_Async_Draw_SyncCycles_Code(PD,GSE_line); break;
                case Async_Pulse1_Code: run_Async_Pulse1_Code(PD,GSE_line); break;
                case Async_PulseM_Code: run_Async_PulseM_Code(PD,GSE_line); break;
                case Async_NonSaturating_Draw_Code: run_Async_NonSaturating_Draw_Code(PD,GSE_line); break;
                case Async_Timing_Cycles_Code: run_Async_Timing_Cycles_Code(PD,GSE_line); break;
                case Async_Timing_Transition_Code: run_Async_Timing_Transition_Code(PD,GSE_line); break;
                case ModulePulse1_Code: run_ModulePulse1_Code(PD,GSE_line); break;
                case ModulePulseM_Code: run_ModulePulseM_Code(PD,GSE_line); break;
                    
                case ModelErrors_Pulse1: run_ModelErrors_Pulse1_Code(PD,GSE_line); break;
               
                default: break;
            }
         //   printf("Line %s DONE\n",GSE_line); getchar();
        }
    }
    fclose(f);
    PD->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,A->active_inputs->N)),p_error);
}

int main (int argc, const char * argv[])
{
    init_random(0); setcolors_PAJEK();
   
    //  run_InSilico_Experimetnal_Package("EMT_Mechanosensing", "VirtualExp_1.txt");
     run_InSilico_Experimetnal_Package("EMT_Mechanosensing_TGFbeta", "VirtualExp_1.txt");
  
    return 0;
}
