//
//  Python_Figures_module.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 10/23/19.
//  Copyright © 2019 Regan_Group. All rights reserved.
//

#ifndef Python_Figures_module_h
#define Python_Figures_module_h

#define PLOT_WIDTH 0.025
//#define PLOT_top 0.92
//#define PLOT_bottom 0.1
//#define PLOT_left 0.13
//#define PLOT_right 0.95
//#define PLOT_hspace 0.1
//#define PLOT_wspace 0.45
#define PLOT_label_fontsize 30

using namespace std;

void initialize_Python_script(FILE *f,const char pt[]){
    fprintf(f, "import matplotlib.pyplot as plt\nimport numpy as np\n");
   // fprintf(f, "from mpl_toolkits import mplot3d\nfrom matplotlib.ticker import NullFormatter\n");
    fprintf(f, "from mpl_toolkits import mplot3d\n");
    fprintf(f, "from matplotlib.backends.backend_pdf import PdfPages\n");
    fprintf(f, "import seaborn as sns\n");
    fprintf(f, "Directory = '%s/%s/_EXP/'\n",MODEL_Directory,pt);

}

string get_Figures_directory_main(Environment_and_hits *EnvHit,
                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                            int run_assync){
    char fn2[400];

    sprintf(fn2, "%s/%s/_EXP/Stochastic_SWEEP_",MODEL_Directory,PD->BRN->path);
    if(run_assync != -1) sprintf(fn2, "%s_ASYNC_",fn2);
    if(((EnvHit->BKR_Hits != NULL)&&(EnvHit->BKR_Hits->N > 0)) || (EnvHit->Scan_Hits->N>0)){
        for (unsigned long int l=1; l<=EnvHit->BKR_Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld-%.2lf",fn2,PD->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[l]],EnvHit->BKR_Hits->B[l],EnvHit->p_hits[l]);
        }
        for (unsigned long int l=1; l<=EnvHit->Scan_Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[l]],EnvHit->Scan_Hits->B[l]);
        }
    }
    else {
        sprintf(fn2,"%s_WILDTYPE",fn2);
        if(run_assync != -1)  sprintf(fn2, "%s_ASYNC_",fn2);
    }
    return(fn2);
}
string get_Figures_directory_main_system(Environment_and_hits *EnvHit,
                                  MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                  int run_assync){
    char fn2[400];
    
    sprintf(fn2, "%s/%s/_EXP/Stochastic_SWEEP_",MODEL_Directory_system,PD->BRN->path);
    if(run_assync != -1)  sprintf(fn2, "%s_ASYNC_",fn2);
    if(((EnvHit->BKR_Hits != NULL)&&(EnvHit->BKR_Hits->N > 0)) || (EnvHit->Scan_Hits->N>0)){
        for (unsigned long int l=1; l<=EnvHit->BKR_Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld-%.2lf",fn2,PD->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[l]],EnvHit->BKR_Hits->B[l],EnvHit->p_hits[l]);
        }
        for (unsigned long int l=1; l<=EnvHit->Scan_Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[l]],EnvHit->Scan_Hits->B[l]);
        }
    }
    else {
        sprintf(fn2,"%s_WILDTYPE",fn2);
        if(run_assync != -1)  sprintf(fn2, "%s_ASYNC_",fn2);
    }
    return(fn2);
}

string get_Figures_directory_now(string maindir,longint_ARRAY *input_index_list,
                                 MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){

    string datadir;
    datadir= maindir + "/";
    for(unsigned int l=1;l<=input_index_list->N;l++){
        if (l>1)  datadir = datadir + "_";
        datadir = datadir + PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_list->A[l]]];
    }
    return (datadir);
}

string get_CC_filename(string fn2,Environment_and_hits *EnvHit,
                       longint_ARRAY *l,
                       unsigned long int lindex,
                       MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    char fn[400];

    sprintf(fn,"%s/CC-FRACTIONS_E-",fn2.c_str());
    for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
        sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
    if(EnvHit->BKR_Hits->N>0) {
        sprintf(fn,"%s___BKG-",fn);
        for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
            sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
    }
    if(EnvHit->Scan_Hits->N>0){
        sprintf(fn,"%s___JointScan-",fn);
        for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
            if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
            else sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
    }
    sprintf(fn,"%s__A%d__BOX_%d.txt",fn,(int)l->A[lindex],BOX_KO_OE);
    return(fn);
}

string get_Irrev_filename(string fn2,Environment_and_hits *EnvHit,
                          longint_ARRAY *l,
                          unsigned long int lindex,
                          MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                          unsigned long int irrev_index){
    char fn[400];
    sprintf(fn,"%s/Irrev_Stats_%s-",fn2.c_str(),Phen_STATs.Irreversible_STATS[irrev_index].c_str());
    for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
        sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
    if(EnvHit->BKR_Hits->N>0) {
        sprintf(fn,"%s___BKG-",fn);
        for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
            sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
    }
    if(EnvHit->Scan_Hits->N>0){
        sprintf(fn,"%s___JointScan-",fn);
        for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
            if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
            else sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
    }
    sprintf(fn,"%s__A%d__BOX_%d.txt",fn,(int)l->A[lindex],BOX_KO_OE);
    return(fn);
}

string get_Rev_filename(string fn2,Environment_and_hits *EnvHit,
                          longint_ARRAY *l,
                          unsigned long int lindex,
                          MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                          unsigned long int rev_index){
    char fn[400];
    sprintf(fn,"%s/RevTransition_Stats_%s_vs_%s-",fn2.c_str(),
            Phen_STATs.Reversible_Transitions_1[rev_index].c_str(),
            Phen_STATs.Reversible_Transitions_2[rev_index].c_str());
    for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
        sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
    if(EnvHit->BKR_Hits->N>0) {
        sprintf(fn,"%s___BKG-",fn);
        for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
            sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
    }
    if(EnvHit->Scan_Hits->N>0){
        sprintf(fn,"%s___JointScan-",fn);
        for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
            if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
            else sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
    }
    sprintf(fn,"%s__A%d__BOX_%d.txt",fn,(int)l->A[lindex],BOX_KO_OE);
    return(fn);
}

bool check_all_files(unsigned long int input_index,
                     string fn2,Environment_and_hits *EnvHit,
                     longint_ARRAY *l,
                     MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                     int run_assync,
                     bool *cc_stats,bool *irrev_stats,bool *rev_stats){
    string fn;
    FILE *f;
    
    *cc_stats=Phen_STATs.CellCycleErrors;
    *irrev_stats=Phen_STATs.Irrev;
    *rev_stats=Phen_STATs.Rev;
    
    for(int j=1;j<=l->N;j++){
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            EnvHit->p_inputs[input_index]=p_env;
            
            fn=get_CC_filename(fn2,EnvHit,l,j,PD);
            f=fopen(fn.c_str(),"r");
            if(f==NULL) *cc_stats = false;
            fclose(f);
        
            for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
                fn=get_Irrev_filename(fn2,EnvHit,l,j,PD,vi);
                f=fopen(fn.c_str(),"r");
                if(f==NULL) *irrev_stats = false;
                fclose(f);
            }
        
            for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
                fn=get_Rev_filename(fn2,EnvHit,l,j,PD,vi);
                f=fopen(fn.c_str(),"r");
                if(f==NULL) *rev_stats = false;
                fclose(f);
            }
        }
    }
    if((*cc_stats!=Phen_STATs.CellCycleErrors)||(*irrev_stats!=Phen_STATs.Irrev)||(*rev_stats!=Phen_STATs.Rev)) return false;
    else return true;
}


bool check_all_files_2D(unsigned long int input_index_1,
                        unsigned long int input_index_2,
                        string fn2,Environment_and_hits *EnvHit,
                        longint_ARRAY *l,
                        MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                        int run_assync,
                        bool *cc_stats,bool *irrev_stats,bool *rev_stats){
    string fn;
    FILE *f;
    
    *cc_stats=Phen_STATs.CellCycleErrors;
    *irrev_stats=Phen_STATs.Irrev;
    *rev_stats=Phen_STATs.Rev;
    
    for(int j=1;j<=l->N;j++){
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            EnvHit->p_inputs[input_index_1]=p_env_1;
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                EnvHit->p_inputs[input_index_2]=p_env_2;
          
                fn=get_CC_filename(fn2,EnvHit,l,j,PD);
                f=fopen(fn.c_str(),"r");
                if(f==NULL) *cc_stats = false;
                fclose(f);
            
                for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
                    fn=get_Irrev_filename(fn2,EnvHit,l,j,PD,vi);
                    f=fopen(fn.c_str(),"r");
                    if(f==NULL) *irrev_stats = false;
                    fclose(f);
                }
            
                for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
                    fn=get_Rev_filename(fn2,EnvHit,l,j,PD,vi);
                    f=fopen(fn.c_str(),"r");
                    if(f==NULL) *rev_stats = false;
                    fclose(f);
                }
            }
        }
    }
    if((*cc_stats!=Phen_STATs.CellCycleErrors)||(*irrev_stats!=Phen_STATs.Irrev)||(*rev_stats!=Phen_STATs.Rev)) return false;
    else return true;
}


void print_CC_file_to_PythonScript(FILE *f,const char fn[],unsigned int ebox_s,unsigned int row){
    fprintf(f,"data_%d = np.loadtxt('%s',skiprows=1)\n",ebox_s, fn);
    fprintf(f,"P_lock=data_%d[:,0]\n",ebox_s);
    fprintf(f,"Normal_CC_%d_%d=data_%d[:,1]\n",ebox_s,row,ebox_s);
    fprintf(f,"G2_4N_slip_%d_%d=data_%d[:,2]\n",ebox_s,row,ebox_s);
    fprintf(f,"Aneupl_4N_%d_%d=data_%d[:,3]\n",ebox_s,row,ebox_s);
    fprintf(f,"Telophase_4N_%d_%d=data_%d[:,4]\n",ebox_s,row,ebox_s);
}
void print_CC_file_to_PythonScript_2D(FILE *f,const char fn[],unsigned int ebox1,unsigned int ebox2,unsigned int row){
    fprintf(f,"data_%d_%d = np.loadtxt('%s',skiprows=1)\n",ebox1,ebox2, fn);
    fprintf(f,"P_lock=data_%d_%d[:,0]\n",ebox1,ebox2);
    fprintf(f,"Normal_CC_%d_%d_%d=data_%d_%d[:,1]\n",ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"G2_4N_slip_%d_%d_%d=data_%d_%d[:,2]\n",ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"Aneupl_4N_%d_%d_%d=data_%d_%d[:,3]\n",ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"Telophase_4N_%d_%d_%d=data_%d_%d[:,4]\n",ebox1,ebox2,row,ebox1,ebox2);
}

void print_Irrev_file_to_PythonScript(FILE *f,const char fn[],unsigned int vi, unsigned int ebox_s,unsigned int row){
    fprintf(f,"Irrev_data_%d = np.loadtxt('%s',skiprows=1)\n",ebox_s, fn);
    fprintf(f,"P_lock=Irrev_data_%d[:,0]\n",ebox_s);
    fprintf(f,"rate_rel_CC_IR%d_%d_%d=Irrev_data_%d[:,1]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"AV_time_to_IR%d_%d_%d=Irrev_data_%d[:,2]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"SD_time_to_IR%d_%d_%d=Irrev_data_%d[:,3]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"HalfLife_IR%d_%d_%d=Irrev_data_%d[:,4]\n", vi, ebox_s,row,ebox_s);
}

void print_Irrev_file_to_PythonScript_2D(FILE *f,const char fn[],unsigned int vi, unsigned int ebox1,unsigned int ebox2,unsigned int row){
    fprintf(f,"Irrev_data_%d_%d = np.loadtxt('%s',skiprows=1)\n",ebox1,ebox2, fn);
    fprintf(f,"P_lock=Irrev_data_%d_%d[:,0]\n",ebox1,ebox2);
    fprintf(f,"rate_rel_CC_IR%d_%d_%d_%d=Irrev_data_%d_%d[:,1]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"AV_time_to_IR%d_%d_%d_%d=Irrev_data_%d_%d[:,2]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"SD_time_to_IR%d_%d_%d_%d=Irrev_data_%d_%d[:,3]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"HalfLife_IR%d_%d_%d_%d=Irrev_data_%d_%d[:,4]\n", vi, ebox1,ebox2,row,ebox1,ebox2);
}

void print_Rev_file_to_PythonScript(FILE *f,const char fn[],unsigned int vi, unsigned int ebox_s,unsigned int row){
    fprintf(f,"Rev_data_%d = np.loadtxt('%s',skiprows=1)\n",ebox_s, fn);
    fprintf(f,"P_lock=Rev_data_%d[:,0]\n",ebox_s);
    fprintf(f,"Time_in_1_Rpair%d_%d_%d=Rev_data_%d[:,1]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"rate_to_1_rel_CC_Rpair%d_%d_%d=Rev_data_%d[:,2]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"AV_window_in_1_Rpair%d_%d_%d=Rev_data_%d[:,3]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"SD_window_in_1_Rpair%d_%d_%d=Rev_data_%d[:,4]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"HalfLife_1_Rpair%d_%d_%d=Rev_data_%d[:,5]\n", vi, ebox_s,row,ebox_s);
    fprintf(f,"Time_in_2_Rpair%d_%d_%d=Rev_data_%d[:,6]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"rate_to_2_rel_CC_Rpair%d_%d_%d=Rev_data_%d[:,7]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"AV_window_in_2_Rpair%d_%d_%d=Rev_data_%d[:,8]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"SD_window_in_2_Rpair%d_%d_%d=Rev_data_%d[:,9]\n",vi, ebox_s,row,ebox_s);
    fprintf(f,"HalfLife_2_Rpair%d_%d_%d=Rev_data_%d[:,10]\n", vi, ebox_s,row,ebox_s);
}
void print_Rev_file_to_PythonScript_2D(FILE *f,const char fn[],unsigned int vi, unsigned int ebox1,unsigned int ebox2,unsigned int row){
    fprintf(f,"Rev_data_%d_%d = np.loadtxt('%s',skiprows=1)\n",ebox1,ebox2, fn);
    fprintf(f,"P_lock=Rev_data_%d_%d[:,0]\n",ebox1,ebox2);
    fprintf(f,"Time_in_1_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,1]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"rate_to_1_rel_CC_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,2]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"AV_window_in_1_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,3]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"SD_window_in_1_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,4]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"HalfLife_1_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,5]\n", vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"Time_in_2_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,6]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"rate_to_2_rel_CC_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,7]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"AV_window_in_2_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,8]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"SD_window_in_2_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,9]\n",vi, ebox1,ebox2,row,ebox1,ebox2);
    fprintf(f,"HalfLife_2_Rpair%d_%d_%d_%d=Rev_data_%d_%d[:,10]\n", vi, ebox1,ebox2,row,ebox1,ebox2);
}


void collate_CC_KO_data_for_increasing_input_in_PythonScript(FILE *f,const char fn[],unsigned int row){
    
    unsigned int ebox=0;
    fprintf(f,"Normal_CC_row%d = (",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"Normal_CC_%d_%d",ebox,row);
    }
    fprintf(f,")\n");

    ebox=0;
    fprintf(f,"G2_4N_slip_row%d = (",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"G2_4N_slip_%d_%d",ebox,row);
    }
    fprintf(f,")\n");
 
    ebox=0;
    fprintf(f,"Aneupl_4N_row%d = (",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"Aneupl_4N_%d_%d",ebox,row);
    }
    fprintf(f,")\n");

    ebox=0;
    fprintf(f,"Telophase_4N_row%d = (",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"Telophase_4N_%d_%d",ebox,row);
    }
    fprintf(f,")\n");
}

void collate_CC_KO_data_for_increasing_input_in_PythonScript_2D(FILE *f,const char fn[],unsigned int row){
    
    unsigned int ebox1=0,ebox2=0;
    fprintf(f,"Normal_CC_row%d = (",row);
    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        ebox1++;
        if(ebox1>1)  fprintf(f,",");
        fprintf(f,"(");
        ebox2=0;
        for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
            ebox2++;
            if(ebox2>1)  fprintf(f,",");
            fprintf(f,"Normal_CC_%d_%d_%d",ebox1,ebox2,row);
        }
        fprintf(f,")");
    }
    fprintf(f,")\n");

    ebox1=0;ebox2=0;
    fprintf(f,"G2_4N_slip_row%d = (",row);
    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        ebox1++;
        if(ebox1>1)  fprintf(f,",");
        fprintf(f,"(");
        ebox2=0;
        for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
            ebox2++;
            if(ebox2>1)  fprintf(f,",");
            fprintf(f,"G2_4N_slip_%d_%d_%d",ebox1,ebox2,row);
        }
        fprintf(f,")");
    }
    fprintf(f,")\n");
    
    ebox1=0;ebox2=0;
    fprintf(f,"Aneupl_4N_row%d = (",row);
    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        ebox1++;
        if(ebox1>1)  fprintf(f,",");
        fprintf(f,"(");
        ebox2=0;
        for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
            ebox2++;
            if(ebox2>1)  fprintf(f,",");
            fprintf(f,"Aneupl_4N_%d_%d_%d",ebox1,ebox2,row);
        }
        fprintf(f,")");
    }
    fprintf(f,")\n");
    
    ebox1=0;ebox2=0;
    fprintf(f,"Telophase_4N_row%d = (",row);
    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        ebox1++;
        if(ebox1>1)  fprintf(f,",");
        fprintf(f,"(");
        ebox2=0;
        for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
            ebox2++;
            if(ebox2>1)  fprintf(f,",");
            fprintf(f,"Telophase_4N_%d_%d_%d",ebox1,ebox2,row);
        }
        fprintf(f,")");
    }
    fprintf(f,")\n");

}

void collate_CC_KO_data_for_increasing_input_in_PythonScript_1D(FILE *f,const char fn[],unsigned int row){
    
    unsigned int ebox=0;
    fprintf(f,"Normal_CC_row%d = (",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"Normal_CC_%d_%d",ebox,row);
    }
    fprintf(f,")\n");

    ebox=0;
    fprintf(f,"G2_4N_slip_row%d = (",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"G2_4N_slip_%d_%d",ebox,row);
    }
    fprintf(f,")\n");
 
    ebox=0;
    fprintf(f,"Aneupl_4N_row%d = (",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"Aneupl_4N_%d_%d",ebox,row);
    }
    fprintf(f,")\n");

    ebox=0;
    fprintf(f,"Telophase_4N_row%d = (",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"Telophase_4N_%d_%d",ebox,row);
    }
    fprintf(f,")\n");
}

void collate_Irrev_KO_data_for_increasing_input_in_PythonScript(FILE *f,const char fn[],unsigned int row){
    unsigned int ebox=0;
    for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
        ebox=0;
        fprintf(f,"rate_rel_CC_IR%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"rate_rel_CC_IR%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");

        ebox=0;
        fprintf(f,"AV_time_to_IR%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"AV_time_to_IR%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox=0;
        fprintf(f,"SD_time_to_IR%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"SD_time_to_IR%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox=0;
        fprintf(f,"HalfLife_IR%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"HalfLife_IR%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
    }
}

void collate_Irrev_KO_data_for_increasing_input_in_PythonScript_2D(FILE *f,const char fn[],unsigned int row){
    unsigned int ebox1=0,ebox2;
    for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
        ebox1=0;
        fprintf(f,"rate_rel_CC_IR%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"rate_rel_CC_IR%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")\n");
        }
        fprintf(f,")\n");

        ebox1=0;
        fprintf(f,"AV_time_to_IR%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"AV_time_to_IR%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")\n");
        }
        fprintf(f,")\n");
        
        ebox1=0;
        fprintf(f,"SD_time_to_IR%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"SD_time_to_IR%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")\n");
        }
        fprintf(f,")\n");
        
        ebox1=0;
        fprintf(f,"HalfLife_IR%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"HalfLife_IR%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");
    }
}

void collate_Rev_KO_data_for_increasing_input_in_PythonScript(FILE *f,const char fn[],unsigned int row){
    unsigned int ebox=0;
    for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
        ebox=0;
        fprintf(f,"Time_in_1_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"Time_in_1_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");

        ebox=0;
        fprintf(f,"rate_to_1_rel_CC_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"rate_to_1_rel_CC_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox=0;
        fprintf(f,"AV_window_in_1_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"AV_window_in_1_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox=0;
        fprintf(f,"SD_window_in_1_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"SD_window_in_1_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox=0;
        fprintf(f,"HalfLife_1_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"HalfLife_1_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox=0;
        fprintf(f,"Time_in_2_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"Time_in_2_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox=0;
        fprintf(f,"rate_to_2_rel_CC_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"rate_to_2_rel_CC_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");

        ebox=0;
        fprintf(f,"AV_window_in_2_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"AV_window_in_2_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox=0;
        fprintf(f,"SD_window_in_2_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"SD_window_in_2_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox=0;
        fprintf(f,"HalfLife_2_Rpair%d_row%d = (",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"HalfLife_2_Rpair%d_%d_%d",vi,ebox,row);
        }
        fprintf(f,")\n");
    }
}

void collate_Rev_KO_data_for_increasing_input_in_PythonScript_2D(FILE *f,const char fn[],unsigned int row){
    unsigned int ebox1=0,ebox2=0;
   
    for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
        ebox1=0;
        fprintf(f,"Time_in_1_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"Time_in_1_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");

        ebox1=0;
        fprintf(f,"rate_to_1_rel_CC_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"rate_to_1_rel_CC_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");
        
        ebox1=0;
        fprintf(f,"AV_window_in_1_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"AV_window_in_1_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");
        
        ebox1=0;
        fprintf(f,"SD_window_in_1_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"SD_window_in_1_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");
        
        ebox1=0;
        fprintf(f,"HalfLife_1_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"HalfLife_1_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");
        
        ebox1=0;
        fprintf(f,"Time_in_2_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"Time_in_2_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");
        
        ebox1=0;
        fprintf(f,"rate_to_2_rel_CC_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"rate_to_2_rel_CC_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");
        
        ebox1=0;
        fprintf(f,"AV_window_in_2_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"AV_window_in_2_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");

        ebox1=0;
        fprintf(f,"SD_window_in_2_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"SD_window_in_2_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");
        
        ebox1=0;
        fprintf(f,"HalfLife_2_Rpair%d_row%d = (",vi,row);
        for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
            ebox1++;
            if(ebox1>1)  fprintf(f,",");
            ebox2=0;
            fprintf(f,"(");
            for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
                ebox2++;
                if(ebox2>1)  fprintf(f,",");
                fprintf(f,"HalfLife_2_Rpair%d_%d_%d_%d",vi,ebox1,ebox2,row);
            }
            fprintf(f,")");
        }
        fprintf(f,")\n");
    }
}


void collapse_CC_data_when_KO_is_absent_in_PythonScript(FILE *f,const char fn[],unsigned int row){
    
    unsigned int box=0, ebox=0;
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        fprintf(f,"WT_Normal_CC_%d_%d = ",ebox,row);
        box=0;
        for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
            box++;
            if(box>1)  fprintf(f," + ");
            fprintf(f,"Normal_CC_%d_%d[%d]",ebox,row,box-1);
        }
        fprintf(f,"\n");
        fprintf(f,"WT_Normal_CC_%d_%d = WT_Normal_CC_%d_%d / %d\n",ebox,row,ebox,row,box);
    }
    
    ebox=0;
    fprintf(f,"WT_Normal_CC_row%d =(",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"WT_Normal_CC_%d_%d",ebox,row);
    }
    fprintf(f,")\n");
    
    box=0, ebox=0;
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        fprintf(f,"WT_G2_4N_slip_%d_%d = ",ebox,row);
        box=0;
        for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
            box++;
            if(box>1)  fprintf(f," + ");
            fprintf(f,"G2_4N_slip_%d_%d[%d]",ebox,row,box-1);
        }
        fprintf(f,"\n");
        fprintf(f,"WT_G2_4N_slip_%d_%d = WT_G2_4N_slip_%d_%d / %d\n",ebox,row,ebox,row,box);
    }
    ebox=0;
    fprintf(f,"WT_G2_4N_slip_row%d =(",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"WT_G2_4N_slip_%d_%d",ebox,row);
    }
    fprintf(f,")\n");
    
    box=0, ebox=0;
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        fprintf(f,"WT_Aneupl_4N_%d_%d = ",ebox,row);
        box=0;
        for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
            box++;
            if(box>1)  fprintf(f," + ");
            fprintf(f,"Aneupl_4N_%d_%d[%d]",ebox,row,box-1);
        }
        fprintf(f,"\n");
        fprintf(f,"WT_Aneupl_4N_%d_%d = WT_Aneupl_4N_%d_%d / %d\n",ebox,row,ebox,row,box);
    }
    ebox=0;
    fprintf(f,"WT_Aneupl_4N_row%d = (",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"WT_Aneupl_4N_%d_%d",ebox,row);
    }
    fprintf(f,")\n");
    
    box=0, ebox=0;
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        fprintf(f,"WT_Telophase_4N_%d_%d = ",ebox,row);
        box=0;
        for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
            box++;
            if(box>1)  fprintf(f," + ");
            fprintf(f,"Telophase_4N_%d_%d[%d]",ebox,row,box-1);
        }
        fprintf(f,"\n");
        fprintf(f,"WT_Telophase_4N_%d_%d = WT_Telophase_4N_%d_%d / %d\n",ebox,row,ebox,row,box);
    }
    ebox=0;
    fprintf(f,"WT_Telophase_4N_row%d =(",row);
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        if(ebox>1)  fprintf(f,",");
        fprintf(f,"WT_Telophase_4N_%d_%d",ebox,row);
    }
    fprintf(f,")\n");
}

void collapse_Irrev_data_when_KO_is_absent_in_PythonScript(FILE *f,const char fn[],unsigned int row){
    
    
    unsigned int box=0, ebox=0;
    for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_rate_rel_CC_IR%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"rate_rel_CC_IR%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_rate_rel_CC_IR%d_%d_%d = WT_rate_rel_CC_IR%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_rate_rel_CC_IR%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_rate_rel_CC_IR%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_AV_time_to_IR%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"AV_time_to_IR%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_AV_time_to_IR%d_%d_%d = WT_AV_time_to_IR%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_AV_time_to_IR%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_AV_time_to_IR%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_SD_time_to_IR%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"SD_time_to_IR%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_SD_time_to_IR%d_%d_%d = WT_SD_time_to_IR%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_SD_time_to_IR%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_SD_time_to_IR%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_HalfLife_IR%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"HalfLife_IR%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_HalfLife_IR%d_%d_%d = WT_HalfLife_IR%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_HalfLife_IR%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_HalfLife_IR%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
    }
}

void collapse_Rev_data_when_KO_is_absent_in_PythonScript(FILE *f,const char fn[],unsigned int row){
    
    
    unsigned int box=0, ebox=0;
    for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_Time_in_1_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"Time_in_1_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_Time_in_1_Rpair%d_%d_%d = WT_Time_in_1_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_Time_in_1_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_Time_in_1_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_rate_to_1_rel_CC_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"rate_to_1_rel_CC_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_rate_to_1_rel_CC_Rpair%d_%d_%d = WT_rate_to_1_rel_CC_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_rate_to_1_rel_CC_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_rate_to_1_rel_CC_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_AV_window_in_1_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"AV_window_in_1_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_AV_window_in_1_Rpair%d_%d_%d = WT_AV_window_in_1_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_AV_window_in_1_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_AV_window_in_1_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_SD_window_in_1_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"SD_window_in_1_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_SD_window_in_1_Rpair%d_%d_%d = WT_SD_window_in_1_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_SD_window_in_1_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_SD_window_in_1_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_HalfLife_1_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"HalfLife_1_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_HalfLife_1_Rpair%d_%d_%d = WT_HalfLife_1_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_HalfLife_1_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_HalfLife_1_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_Time_in_2_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"Time_in_2_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_Time_in_2_Rpair%d_%d_%d = WT_Time_in_2_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_Time_in_2_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_Time_in_2_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_rate_to_2_rel_CC_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"rate_to_2_rel_CC_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_rate_to_2_rel_CC_Rpair%d_%d_%d = WT_rate_to_2_rel_CC_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_rate_to_2_rel_CC_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_rate_to_2_rel_CC_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_AV_window_in_2_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"AV_window_in_2_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_AV_window_in_2_Rpair%d_%d_%d = WT_AV_window_in_2_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_AV_window_in_2_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_AV_window_in_2_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_SD_window_in_2_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"SD_window_in_2_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_SD_window_in_2_Rpair%d_%d_%d = WT_SD_window_in_2_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_SD_window_in_2_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_SD_window_in_2_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
        
        box=0, ebox=0;
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            fprintf(f,"WT_HalfLife_2_Rpair%d_%d_%d = ",vi,ebox,row);
            box=0;
            for(double p=0;p<=1.; p+=1./(double)BOX_KO_OE){
                box++;
                if(box>1)  fprintf(f," + ");
                fprintf(f,"HalfLife_2_Rpair%d_%d_%d[%d]",vi,ebox,row,box-1);
            }
            fprintf(f,"\n");
            fprintf(f,"WT_HalfLife_2_Rpair%d_%d_%d = WT_HalfLife_2_Rpair%d_%d_%d  / %d\n",vi,ebox,row,vi,ebox,row,box);
        }
        ebox=0;
        fprintf(f,"WT_HalfLife_2_Rpair%d_row%d =(",vi,row);
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,",");
            fprintf(f,"WT_HalfLife_2_Rpair%d_%d_%d ",vi,ebox,row);
        }
        fprintf(f,")\n");
    }
}

void read_CC_datafiles_in_Python_script(FILE *f,
                                     unsigned long int input_index,
                                     string fn2,Environment_and_hits *EnvHit,
                                     longint_ARRAY *l,
                                     unsigned long int lindex,
                                     MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    string fn;
    unsigned int ebox=0;
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        EnvHit->p_inputs[input_index]=p_env;
        
        fn=get_CC_filename(fn2,EnvHit,l,lindex,PD);
        print_CC_file_to_PythonScript(f,fn.c_str(),ebox,0);
    }
    collate_CC_KO_data_for_increasing_input_in_PythonScript(f,fn.c_str(),0);
}

void read_CC_datafiles_in_Python_script_2D(FILE *f,
                                     unsigned long int input_index_1,
                                     unsigned long int input_index_2,
                                     string fn2,Environment_and_hits *EnvHit,
                                     longint_ARRAY *l,
                                     unsigned long int lindex,
                                     MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    string fn;
    unsigned int ebox1=0, ebox2=0;
    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        ebox1++;
        EnvHit->p_inputs[input_index_1]=p_env_1;
        ebox2=0;
        for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
            ebox2++;
            EnvHit->p_inputs[input_index_2]=p_env_2;
      
            fn=get_CC_filename(fn2,EnvHit,l,lindex,PD);
            print_CC_file_to_PythonScript_2D(f,fn.c_str(),ebox1,ebox2,0);
        }
    }
    collate_CC_KO_data_for_increasing_input_in_PythonScript_2D(f,fn.c_str(),0);
}

void read_CC_datafiles_in_Python_script_1D(FILE *f,
                                     unsigned long int input_index_x,
                                     string fn2,Environment_and_hits *EnvHit,
                                     longint_ARRAY *l,
                                     unsigned long int lindex,
                                     MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    string fn;
    unsigned int ebox=0;
    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        EnvHit->p_inputs[input_index_x]=p_env_1;
        ebox++;
        fn=get_CC_filename(fn2,EnvHit,l,lindex,PD);
        print_CC_file_to_PythonScript(f,fn.c_str(),ebox,0);
    }
    collate_CC_KO_data_for_increasing_input_in_PythonScript_1D(f,fn.c_str(),0);
}


void read_Irrev_datafiles_in_Python_script(FILE *f,
                                        unsigned long int input_index,
                                        string fn2,Environment_and_hits *EnvHit,
                                        longint_ARRAY *l,
                                        unsigned long int lindex,
                                        MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    string fn;
    unsigned int ebox=0;
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        EnvHit->p_inputs[input_index]=p_env;
        for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
            fn=get_Irrev_filename(fn2,EnvHit,l,lindex,PD,vi);
            print_Irrev_file_to_PythonScript(f,fn.c_str(),vi,ebox,0);
        }
    }
    collate_Irrev_KO_data_for_increasing_input_in_PythonScript(f,fn.c_str(),0);

}

void read_Rev_datafiles_in_Python_script(FILE *f,
                                           unsigned long int input_index,
                                           string fn2,Environment_and_hits *EnvHit,
                                           longint_ARRAY *l,
                                           unsigned long int lindex,
                                           MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    string fn;
    unsigned int ebox=0;
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        ebox++;
        EnvHit->p_inputs[input_index]=p_env;
        for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
            fn=get_Rev_filename(fn2,EnvHit,l,lindex,PD,vi);
            print_Rev_file_to_PythonScript(f,fn.c_str(),vi,ebox,0);
        }
    }
    collate_Rev_KO_data_for_increasing_input_in_PythonScript(f,fn.c_str(),0);
}

void read_Irrev_datafiles_in_Python_script_2D(FILE *f,
                                        unsigned long int input_index_1,
                                        unsigned long int input_index_2,
                                        string fn2,Environment_and_hits *EnvHit,
                                        longint_ARRAY *l,
                                        unsigned long int lindex,
                                        MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    string fn;
    unsigned int ebox1=0,ebox2=0;
    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        ebox1++;
        EnvHit->p_inputs[input_index_1]=p_env_1;
        ebox2=0;
        for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
            ebox2++;
            EnvHit->p_inputs[input_index_2]=p_env_2;
            for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
                fn=get_Irrev_filename(fn2,EnvHit,l,lindex,PD,vi);
                print_Irrev_file_to_PythonScript_2D(f,fn.c_str(),vi,ebox1,ebox2,0);
            }
        }
    }
    collate_Irrev_KO_data_for_increasing_input_in_PythonScript_2D(f,fn.c_str(),0);

}

void read_Rev_datafiles_in_Python_script_2D(FILE *f,
                                           unsigned long int input_index_1,
                                           unsigned long int input_index_2,
                                           string fn2,Environment_and_hits *EnvHit,
                                           longint_ARRAY *l,
                                           unsigned long int lindex,
                                           MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    string fn;
    unsigned int ebox1=0,ebox2=0;
    for(double p_env_1=0;p_env_1<=1.; p_env_1+=1./(double)BOX_ENV){
        ebox1++;
        EnvHit->p_inputs[input_index_1]=p_env_1;
        ebox2=0;
        for(double p_env_2=0;p_env_2<=1.; p_env_2+=1./(double)BOX_ENV){
            ebox2++;
            EnvHit->p_inputs[input_index_2]=p_env_2;
            for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
                fn=get_Rev_filename(fn2,EnvHit,l,lindex,PD,vi);
                print_Rev_file_to_PythonScript_2D(f,fn.c_str(),vi,ebox1,ebox2,0);
            }
        }
    }
    collate_Rev_KO_data_for_increasing_input_in_PythonScript_2D(f,fn.c_str(),0);
}


void draw_NOHITs_stat_1input_panel(FILE *f,
                                   unsigned long int input_index,
                                   string fn2,Environment_and_hits *EnvHit,
                                   longint_ARRAY *l,
                                   unsigned long int lindex,
                                   MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                   unsigned int row){

    fprintf(f,"Percent_signal = (");
    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
        if(p_env>0)  fprintf(f,", ");
        fprintf(f,"%.3lf",p_env);
    }
    fprintf(f,")\n");
        
    fprintf(f,"width = %lf\n",PLOT_WIDTH);
    fprintf(f,"plt.tick_params(axis='both', which='major', labelsize=%d)\n",(int)(PLOT_label_fontsize*0.75));
    fprintf(f,"fig, axes = plt.subplots(1, 2+%lu, sharey='none')\n",Phen_STATs.Reversible_Transitions_1.size());
    
    int box_drawn=0;
    if(Phen_STATs.CellCycleErrors){
       fprintf(f,"axes[%d].xaxis.set_tick_params(which='both', labelsize=%d, labelbottom=True)\n",box_drawn,(int)(PLOT_label_fontsize*0.75));
        fprintf(f,"axes[%d].yaxis.set_tick_params(which='both', labelsize=%d, labelbottom=True)\n",box_drawn,(int)(PLOT_label_fontsize*0.75));
        //fprintf(f,"axes[%d].set_ylim([0,1])\n",box_drawn);
        box_drawn++;
    }
    if(Phen_STATs.Irrev){
        fprintf(f,"axes[%d].xaxis.set_tick_params(which='both', labelsize=%d, labelbottom=True)\n",box_drawn,(int)(PLOT_label_fontsize*0.75));
        fprintf(f,"axes[%d].yaxis.set_tick_params(which='both', labelsize=%d, labelbottom=True)\n",box_drawn,(int)(PLOT_label_fontsize*0.75));
        //if(Phen_STATs.CellCycleErrors) fprintf(f,"axes[%d].set_ylim([0,1])\n",box_drawn);
        box_drawn++;
    }
    if(Phen_STATs.Irrev){
        for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
            fprintf(f,"axes[%d+%d].xaxis.set_tick_params(which='both', labelsize=%d, labelbottom=True)\n",box_drawn,vi,(int)(PLOT_label_fontsize*0.75));
            fprintf(f,"axes[%d+%d].yaxis.set_tick_params(which='both', labelsize=%d, labelbottom=True)\n",box_drawn,vi,(int)(PLOT_label_fontsize*0.75));
            //if(Phen_STATs.CellCycleErrors) fprintf(f,"axes[%d+%d].set_ylim([0,1])\n",box_drawn,vi);
        }
    }
    fprintf(f,"fig.set_size_inches(10*(%d+%lu),10)\n",box_drawn, Phen_STATs.Reversible_Transitions_1.size());
    
    box_drawn=0;
    if(Phen_STATs.CellCycleErrors){
        fprintf(f,"bottom1 = WT_Normal_CC_row%d\n",row);
        
        unsigned int ebox = 0;
        fprintf(f,"bottom2 = (");
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,", ");
            fprintf(f,"WT_Normal_CC_%d_%d + WT_G2_4N_slip_%d_%d",ebox,row,ebox,row);
        }
        fprintf(f,")\n");
        
        ebox = 0;
        fprintf(f,"bottom3 = (");
        for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
            ebox++;
            if(ebox>1)  fprintf(f,", ");
            fprintf(f,"WT_Normal_CC_%d_%d + WT_G2_4N_slip_%d_%d + WT_Aneupl_4N_%d_%d",ebox,row,ebox,row,ebox,row);
        }
        fprintf(f,")\n");
        
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Normal_CC_row%d, 2*width, align='center', label=\"%s\")\n",box_drawn,row,"Normal cycle");
        fprintf(f,"axes[%d].bar(Percent_signal, WT_G2_4N_slip_row%d, 2*width, align='center', bottom=bottom1, label=\"%s\")\n",box_drawn,row,"G2->4N");
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Aneupl_4N_row%d, 2*width, align='center', bottom=bottom2, label=\"%s\")\n",box_drawn,row,"Aneuploidy->4N");
        fprintf(f,"axes[%d].bar(Percent_signal, WT_Telophase_4N_row%d, 2*width, align='center', bottom=bottom3, label=\"%s\")\n",box_drawn,row,"Telophase->4N");
        
        box_drawn++;
    }
    
    if(Phen_STATs.Irrev){
        fprintf(f,"barGap = %lf\n",(0.25/(double)BOX_ENV) /(double)Phen_STATs.Irreversible_STATS.size() );
        
        //fprintf(f,"Percent_signal_sh = Percent_signal\n");
        fprintf(f,"Percent_signal_sh = [x - barGap * %lf + barGap/2. for x in Percent_signal]\n",Phen_STATs.Irreversible_STATS.size()/2.);
         //fprintf(f,"bt_irr=0\n");
         for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
           if(Phen_STATs.CellCycleErrors){
               if(vi==0) {fprintf(f,"axes[%d].bar(Percent_signal_sh, WT_rate_rel_CC_IR%d_row%d, 2*width, align='center', label=\"%s\")\n",box_drawn,vi,row,Phen_STATs.Irreversible_STATS[vi].c_str());
                   fprintf(f,"bt_irr = WT_rate_rel_CC_IR%d_row%d\n", vi,row);
               }
               else {
                   fprintf(f,"axes[%d].bar(Percent_signal_sh, WT_rate_rel_CC_IR%d_row%d, 2*width, align='center', bottom=bt_irr, label=\"%s\")\n",box_drawn,vi,row,Phen_STATs.Irreversible_STATS[vi].c_str());
                   unsigned int ebox = 0;
                   fprintf(f,"bt_irr = (");
                   for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
                       ebox++;
                       if(ebox>1)  fprintf(f,", ");
                       fprintf(f,"bt_irr[%d] + WT_rate_rel_CC_IR%d_row%d[%d]",ebox-1,vi,row,ebox-1);
                   }
                   fprintf(f,")\n");
               }
            }
            else {
                if(vi==0) {
                    fprintf(f,"axes[%d].bar(Percent_signal_sh, WT_AV_time_to_IR%d_row%d, 2*width, align='center', label=\"%s\")\n",box_drawn,vi,row,Phen_STATs.Irreversible_STATS[vi].c_str());
                    fprintf(f,"bt_irr = WT_AV_time_to_IR%d_row%d\n", vi,row);
                }
                else {fprintf(f,"axes[%d].bar(Percent_signal_sh, WT_AV_time_to_IR%d_row%d, 2*width, align='center', bottom=bt_irr, label=\"%s\")\n",box_drawn,vi,row,Phen_STATs.Irreversible_STATS[vi].c_str());
                    unsigned int ebox = 0;
                    fprintf(f,"bt_irr = (");
                    for(double p_env=0;p_env<=1.; p_env+=1./(double)BOX_ENV){
                        ebox++;
                        if(ebox>1)  fprintf(f,", ");
                        fprintf(f,"bt_irr[%d] + WT_rate_rel_CC_IR%d_row%d[%d]",ebox-1,vi,row,ebox-1);
                    }
                    fprintf(f,")\n");
                }
            }
            fprintf(f,"Percent_signal_sh = [x  + barGap for x in Percent_signal_sh]\n");
        }
        box_drawn++;
    }
    //fprintf(f,"Percent_signal_sh = [x - barGap * 0.25 for x in Percent_signal]\n");
    //fprintf(f,"Percent_signal_sh = Percent_signal\n");
//    fprintf(f,"Percent_signal_sh = [x - barGap * %lf + barGap/2. for x in Percent_signal]\n",Phen_STATs.Reversible_Transitions_1.size()/2.);
//        fprintf(f,"axes[2].bar(Percent_signal_sh, WT_rate_to_1_rel_CC_Rpair%d_row%d, width, align='center', label=\"%s\")\n",vi,row,Phen_STATs.Reversible_Transitions_1[vi].c_str());
//        fprintf(f,"axes[2].bar(Percent_signal_sh, WT_rate_to_2_rel_CC_Rpair%d_row%d, width, align='center',  bottom=WT_rate_to_1_rel_CC_Rpair%d_row%d, label=\"%s\")\n",vi,row,vi,row,Phen_STATs.Reversible_Transitions_2[vi].c_str());
//        fprintf(f,"Percent_signal_sh = [x  + barGap for x in Percent_signal_sh]\n");
//    }

    if(Phen_STATs.Rev){
        for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
            fprintf(f,"Percent_signal_sh = [x - barGap * %lf + barGap/2. for x in Percent_signal]\n",Phen_STATs.Reversible_Transitions_1.size()/2.);
            if(Phen_STATs.CellCycleErrors){
                fprintf(f,"axes[%d+%d].bar(Percent_signal_sh, WT_Time_in_1_Rpair%d_row%d, width, align='center', label=\"%s\")\n",box_drawn,vi,vi,row,Phen_STATs.Reversible_Transitions_1[vi].c_str());
                fprintf(f,"axes[%d+%d].bar(Percent_signal_sh, WT_Time_in_2_Rpair%d_row%d, width, align='center',  bottom=WT_Time_in_1_Rpair%d_row%d, label=\"%s\")\n",box_drawn,vi,vi,row,vi,row,Phen_STATs.Reversible_Transitions_2[vi].c_str());
            }
            else{
                fprintf(f,"axes[%d+%d].bar(Percent_signal_sh, WT_AV_window_in_1_Rpair%d_row%d, width, align='center', label=\"%s\")\n",box_drawn,vi,vi,row,Phen_STATs.Reversible_Transitions_1[vi].c_str());
                fprintf(f,"axes[%d+%d].bar(Percent_signal_sh, WT_AV_window_in_1_Rpair%d_row%d, width, align='center',  bottom=WT_Time_in_1_Rpair%d_row%d, label=\"%s\")\n",box_drawn,vi,vi,row,vi,row,Phen_STATs.Reversible_Transitions_2[vi].c_str());
            }
            fprintf(f,"Percent_signal_sh = [x  + barGap for x in Percent_signal_sh]\n");
        }
        box_drawn++;
    }
 //   fprintf(f,"fig.gca().yaxis.set_minor_formatter(NullFormatter())\n");
 //   fprintf(f,"fig.subplots_adjust(top=%lf, bottom=%lf, left=%lf, right=%lf, hspace=%lf,wspace=%lf)\n",PLOT_top,PLOT_bottom,PLOT_left,PLOT_right,PLOT_hspace,PLOT_wspace);

    box_drawn=0;
    if(Phen_STATs.CellCycleErrors){
        fprintf(f,"axes[%d].set_xlabel('%% %s',fontsize=%d)\n",box_drawn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]],PLOT_label_fontsize);
        fprintf(f,"axes[%d].legend(prop={'size': %d})\n",box_drawn,(int)(PLOT_label_fontsize*0.75));
        fprintf(f,"axes[%d].set_title('Cell cycle', fontsize=%d)\n",box_drawn,PLOT_label_fontsize);
        fprintf(f,"axes[%d].set_ylabel('Cell cycle events / WT cell cycle',fontsize=%d)\n",box_drawn,PLOT_label_fontsize);
        box_drawn++;
    }
    if(Phen_STATs.Irrev){
        fprintf(f,"axes[%d].set_xlabel('%% %s',fontsize=%d)\n",box_drawn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]],PLOT_label_fontsize);
        if(Phen_STATs.CellCycleErrors)
             fprintf(f,"axes[%d].set_ylabel('%% time in states',fontsize=%d)\n",box_drawn,PLOT_label_fontsize);
        else fprintf(f,"axes[%d].set_ylabel('Av. time to transition',fontsize=%d)\n",box_drawn,PLOT_label_fontsize);
        fprintf(f,"axes[%d].legend(prop={'size': %d})\n",box_drawn,(int)(PLOT_label_fontsize*0.75));
        fprintf(f,"axes[%d].set_title('Irreversible transitions', fontsize=%d)\n",box_drawn,PLOT_label_fontsize);
        box_drawn++;
    }
    if(Phen_STATs.Rev){
        for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
            fprintf(f,"axes[%d+%d].set_xlabel('%% %s',fontsize=%d)\n",box_drawn,vi,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]],PLOT_label_fontsize);
            if(Phen_STATs.CellCycleErrors)
                  fprintf(f,"axes[%d+%d].set_ylabel('%% time in states',fontsize=%d)\n",box_drawn,vi,PLOT_label_fontsize);
            else  fprintf(f,"axes[%d+%d].set_ylabel('Av. time-interval in states',fontsize=%d)\n",box_drawn,vi,PLOT_label_fontsize);
            fprintf(f,"axes[%d+%d].legend(prop={'size': %d})\n",box_drawn,vi,(int)(PLOT_label_fontsize*0.75));
            fprintf(f,"axes[%d+%d].set_title('%% time in %s vs. %s', fontsize=%d)\n",box_drawn,vi,Phen_STATs.Reversible_Transitions_1[vi].c_str(),Phen_STATs.Reversible_Transitions_2[vi].c_str(),PLOT_label_fontsize);
        }
    }
    fprintf(f,"plt.tight_layout()\n");
    //fprintf(f,"plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)\n");
}

void draw_HITs_stat_1input_panel(FILE *f,
                                   unsigned long int input_index,
                                   string fn2,Environment_and_hits *EnvHit,
                                   MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                   unsigned int row){
    
    fprintf(f,"P_lock = (");
    for(double ph=0;ph<=1.; ph+=1./(double)BOX_KO_OE){
        if(ph>0)  fprintf(f,", ");
        fprintf(f,"%.3lf",ph);
    }
    fprintf(f,")\n");
    
    fprintf(f,"width = %lf\n",PLOT_WIDTH);
    fprintf(f,"fig, axes = plt.subplots(2+%lu, %d, sharex=True, sharey='row')\n",Phen_STATs.Reversible_Transitions_1.size(),BOX_ENV/2);
    fprintf(f,"fig.set_size_inches(10*%d,10*(2+%lu))\n",BOX_ENV/2,Phen_STATs.Reversible_Transitions_1.size());
    fprintf(f,"plt.tick_params(axis='both', which='major', labelsize=%d)\n",(int)(PLOT_label_fontsize*0.75));
   // fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
   // if(Phen_STATs.WT_cell_cycle>1){
   //     fprintf(f,"\taxes[0, x].set_ylim([0,1])\n");
   //     fprintf(f,"\taxes[1, x].set_ylim([0,1])\n");
   // }
    
  /*
    if(Phen_STATs.Rev){
        fprintf(f,"for i in range(2+%lu):\n",Phen_STATs.Reversible_Transitions_1.size());
        fprintf(f,"\t\taxes[i, x].tick_params(axis=\"x\", labelsize=%d,labelbottom=True)\n",(int)(PLOT_label_fontsize*0.75));
        fprintf(f,"\t\taxes[i, x].tick_params(axis=\"y\", labelsize=%d,labelleft=True)\n",(int)(PLOT_label_fontsize*0.75));
    }
*/
    ////
    if(Phen_STATs.CellCycleErrors){
       fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
        fprintf(f,"\taxes[0, x].bar(P_lock, Normal_CC_row%d[2*x+1], 2*width, align='center', label=\"%s\")\n",row,"Normal cycle");
        fprintf(f,"\taxes[0, x].bar(P_lock, G2_4N_slip_row%d[2*x+1], 2*width, align='center', bottom=Normal_CC_row%d[2*x+1], label=\"%s\")\n",row,row,"G2->4N");
        fprintf(f,"\taxes[0, x].bar(P_lock, Aneupl_4N_row%d[2*x+1], 2*width, align='center', bottom=G2_4N_slip_row%d[2*x+1] + Normal_CC_row%d[2*x+1], label=\"%s\")\n",row,row,row,"Aneuploidy->4N");
        fprintf(f,"\taxes[0, x].bar(P_lock, Telophase_4N_row%d[2*x+1], 2*width, align='center', bottom=Aneupl_4N_row%d[2*x+1] + G2_4N_slip_row%d[2*x+1] + Normal_CC_row%d[2*x+1], label=\"%s\")\n",row,row,row,row,"Telophase->4N");
        fprintf(f,"\taxes[0, x].tick_params(axis='both', which='major', labelsize=22)\n");
    }
    
    if(Phen_STATs.Irrev){
        fprintf(f,"barGap = %lf\n",(1./(double)BOX_ENV) /(double)Phen_STATs.Irreversible_STATS.size() );
        
        //fprintf(f,"Percent_signal_sh = Percent_signal\n");
        fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
            for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
                fprintf(f,"\tPercent_signal_sh = [a - barGap * %lf + barGap/2. + %d * barGap for a in P_lock]\n",Phen_STATs.Irreversible_STATS.size()/2.,vi);
                fprintf(f,"\taxes[1, x].bar(Percent_signal_sh, rate_rel_CC_IR%d_row%d[2*x+1], 2*width, align='center', label=\"%s\")\n",vi,row,Phen_STATs.Irreversible_STATS[vi].c_str());
                fprintf(f,"\taxes[1, x].tick_params(axis='both', which='major', labelsize=22)\n");
                //fprintf(f,"\tPercent_signal_sh = [y  + barGap for y in Percent_signal_sh]\n");
        }
    }
        
    //fprintf(f,"Percent_signal_sh = [x - barGap * 0.25 for x in Percent_signal]\n");
    //fprintf(f,"Percent_signal_sh = Percent_signal\n");
//    fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
//    fprintf(f,"\tPercent_signal_sh = [a - barGap * %lf + barGap/2. for a in P_lock]\n",Phen_STATs.Reversible_Transitions_1.size()/2.);
//    for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
//        fprintf(f,"\taxes[2, x].bar(Percent_signal_sh, rate_to_1_rel_CC_Rpair%d_row%d[2*x+1], width, align='center', label=\"%s\")\n",vi,row,Phen_STATs.Reversible_Transitions_1[vi].c_str());
//        fprintf(f,"\taxes[2, x].bar(Percent_signal_sh, rate_to_2_rel_CC_Rpair%d_row%d[2*x+1], width, align='center',  bottom=rate_to_1_rel_CC_Rpair%d_row%d[2*x+1], label=\"%s\")\n",vi,row,vi,row,Phen_STATs.Reversible_Transitions_2[vi].c_str());
//        fprintf(f,"\tPercent_signal_sh = [y  + barGap for y in Percent_signal_sh]\n");
//    }
 
    if(Phen_STATs.Rev){
        fprintf(f,"barGap = %lf\n",(1./(double)BOX_ENV) /(double)Phen_STATs.Reversible_Transitions_1.size() );
       
        for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
            fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
            fprintf(f,"\tPercent_signal_sh = [a - barGap * %lf + barGap/2. for a in P_lock]\n",Phen_STATs.Reversible_Transitions_1.size()/2.);
            fprintf(f,"\taxes[%d, x].bar(Percent_signal_sh, Time_in_1_Rpair%d_row%d[2*x+1], width, align='center', label=\"%s\")\n",vi+2,vi,row,Phen_STATs.Reversible_Transitions_1[vi].c_str());
            fprintf(f,"\taxes[%d, x].bar(Percent_signal_sh, Time_in_2_Rpair%d_row%d[2*x+1], width, align='center',  bottom=Time_in_1_Rpair%d_row%d[2*x+1], label=\"%s\")\n",vi+2,vi,row,vi,row,Phen_STATs.Reversible_Transitions_2[vi].c_str());
            fprintf(f,"\tPercent_signal_sh = [y  + barGap for y in Percent_signal_sh]\n");
            fprintf(f,"axes[%d, x].tick_params(axis='both', which='major', labelsize=22)\n",vi+2);
        }
    }
    
    unsigned int ebox_s=0;
    char xlabel[400] = "";
    for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
        if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
        sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
        if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
        else sprintf (xlabel, "%s%s",xlabel,"OE");
    }
    ebox_s=0;
    for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
     //   fprintf(f,"axes[0, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
     //   fprintf(f,"axes[1, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
        fprintf(f,"axes[%lu+1, %d].set_xlabel('%s',fontsize=%d)\n",Phen_STATs.Reversible_Transitions_1.size(),ebox_s,xlabel,PLOT_label_fontsize);
        fprintf(f,"axes[0, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,
                p_env, PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]],PLOT_label_fontsize);
        ebox_s+=1;
    }
    
    fprintf(f,"axes[0, 0].set_ylabel('Cell cycle events / WT cell cycle',fontsize=%d)\n",PLOT_label_fontsize);
    fprintf(f,"axes[1, 0].set_ylabel('Transitions / WT cell cycle',fontsize=%d)\n",PLOT_label_fontsize);
    
    for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++)
        fprintf(f,"axes[%d+2, 0].set_ylabel('%% time in %s vs. %s',fontsize=%d)\n",
                vi,
                Phen_STATs.Reversible_Transitions_1[vi].c_str(),
                Phen_STATs.Reversible_Transitions_2[vi].c_str(),
                PLOT_label_fontsize);

    ebox_s=0;
    for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
        fprintf(f,"axes[0, %d].legend(prop={'size': %d})\n",ebox_s,(int)(PLOT_label_fontsize*0.75));
        fprintf(f,"axes[1, %d].legend(prop={'size': %d})\n",ebox_s,(int)(PLOT_label_fontsize*0.75));
        for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++)
            fprintf(f,"axes[2+%d, %d].legend(prop={'size': %d})\n",vi,ebox_s,(int)(PLOT_label_fontsize*0.75));
        ebox_s++;
    }
    fprintf(f,"plt.tight_layout()\n");
    //fprintf(f,"plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)\n");
}

void draw_HITs_stat_1input_panel_inverted(FILE *f,
                                   unsigned long int input_index,
                                   string fn2,Environment_and_hits *EnvHit,
                                   MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                   unsigned int row){
    
    fprintf(f,"P_lock = (");
    for(double ph=0;ph<=1.; ph+=1./(double)BOX_KO_OE){
        if(ph>0)  fprintf(f,", ");
        fprintf(f,"%.3lf",ph);
    }
    fprintf(f,")\n");
    
    fprintf(f,"width = %lf\n",PLOT_WIDTH);
    fprintf(f,"fig, axes = plt.subplots(2+%lu, %d, sharex=True, sharey='row')\n",Phen_STATs.Reversible_Transitions_1.size(),BOX_ENV/2);
    fprintf(f,"fig.set_size_inches(10*%d,10*(2+%lu))\n",BOX_ENV/2,Phen_STATs.Reversible_Transitions_1.size());
    fprintf(f,"plt.tick_params(axis='both', which='major', labelsize=%d)\n",(int)(PLOT_label_fontsize*0.75));
    fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
  //  if(Phen_STATs.WT_cell_cycle>1){
  //      fprintf(f,"\taxes[0, x].set_ylim([0,1])\n");
   //     fprintf(f,"\taxes[1, x].set_ylim([0,1])\n");
   // }
    fprintf(f,"\tfor i in range(2+%lu):\n",Phen_STATs.Reversible_Transitions_1.size());
    fprintf(f,"\t\taxes[i, x].tick_params(axis=\"x\", labelsize=%d,labelbottom=True)\n",(int)(PLOT_label_fontsize*0.75));
    fprintf(f,"\t\taxes[i, x].tick_params(axis=\"y\", labelsize=%d,labelleft=True)\n",(int)(PLOT_label_fontsize*0.75));


    ////
    fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
    fprintf(f,"\taxes[0, x].bar(P_lock, Normal_CC_row%d[2*x+1], 2*width, align='center', label=\"%s\")\n",row,"Normal cycle");
    fprintf(f,"\taxes[0, x].bar(P_lock, G2_4N_slip_row%d[2*x+1], 2*width, align='center', bottom=Normal_CC_row%d[2*x+1], label=\"%s\")\n",row,row,"G2->4N");
    fprintf(f,"\taxes[0, x].bar(P_lock, Aneupl_4N_row%d[2*x+1], 2*width, align='center', bottom=G2_4N_slip_row%d[2*x+1] + Normal_CC_row%d[2*x+1], label=\"%s\")\n",row,row,row,"Aneuploidy->4N");
    fprintf(f,"\taxes[0, x].bar(P_lock, Telophase_4N_row%d[2*x+1], 2*width, align='center', bottom=Aneupl_4N_row%d[2*x+1] + G2_4N_slip_row%d[2*x+1] + Normal_CC_row%d[2*x+1], label=\"%s\")\n",row,row,row,row,"Telophase->4N");
    
    fprintf(f,"barGap = %lf\n",(0.25/(double)BOX_ENV) /(double)Phen_STATs.Irreversible_STATS.size() );
    
    //fprintf(f,"Percent_signal_sh = Percent_signal\n");
    fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
    fprintf(f,"\tPercent_signal_sh = [a - barGap * %lf + barGap/2. for a in P_lock]\n",Phen_STATs.Irreversible_STATS.size()/2.);
    for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
        fprintf(f,"\taxes[1, x].bar(Percent_signal_sh, rate_rel_CC_IR%d_row%d[2*x+1], 2*width, align='center', label=\"%s\")\n",vi,row,Phen_STATs.Irreversible_STATS[vi].c_str());
        fprintf(f,"\tPercent_signal_sh = [y  + barGap for y in Percent_signal_sh]\n");
    }
    
    //fprintf(f,"Percent_signal_sh = [x - barGap * 0.25 for x in Percent_signal]\n");
    //fprintf(f,"Percent_signal_sh = Percent_signal\n");
//    fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
//    fprintf(f,"\tPercent_signal_sh = [a - barGap * %lf + barGap/2. for a in P_lock]\n",Phen_STATs.Reversible_Transitions_1.size()/2.);
//    for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
//        fprintf(f,"\taxes[2, x].bar(Percent_signal_sh, rate_to_1_rel_CC_Rpair%d_row%d[2*x+1], width, align='center', label=\"%s\")\n",vi,row,Phen_STATs.Reversible_Transitions_1[vi].c_str());
//        fprintf(f,"\taxes[2, x].bar(Percent_signal_sh, rate_to_2_rel_CC_Rpair%d_row%d[2*x+1], width, align='center',  bottom=rate_to_1_rel_CC_Rpair%d_row%d[2*x+1], label=\"%s\")\n",vi,row,vi,row,Phen_STATs.Reversible_Transitions_2[vi].c_str());
//        fprintf(f,"\tPercent_signal_sh = [y  + barGap for y in Percent_signal_sh]\n");
//    }
 
    for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
        fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
        fprintf(f,"\tPercent_signal_sh = [a - barGap * %lf + barGap/2. for a in P_lock]\n",Phen_STATs.Reversible_Transitions_1.size()/2.);
        fprintf(f,"\taxes[%d, x].bar(Percent_signal_sh, Time_in_1_Rpair%d_row%d[2*x+1], width, align='center', label=\"%s\")\n",vi+2,vi,row,Phen_STATs.Reversible_Transitions_1[vi].c_str());
        fprintf(f,"\taxes[%d, x].bar(Percent_signal_sh, Time_in_2_Rpair%d_row%d[2*x+1], width, align='center',  bottom=Time_in_1_Rpair%d_row%d[2*x+1], label=\"%s\")\n",vi+2,vi,row,vi,row,Phen_STATs.Reversible_Transitions_2[vi].c_str());
        fprintf(f,"\tPercent_signal_sh = [y  + barGap for y in Percent_signal_sh]\n");
    }

    unsigned int ebox_s=0;
    char xlabel[400] = "";
    for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
        if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
        sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
        if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
        else sprintf (xlabel, "%s%s",xlabel,"OE");
    }
    ebox_s=0;
    for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
     //   fprintf(f,"axes[0, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
     //   fprintf(f,"axes[1, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
        fprintf(f,"axes[%lu+1, %d].set_xlabel('%s',fontsize=%d)\n",Phen_STATs.Reversible_Transitions_1.size(),ebox_s,xlabel,PLOT_label_fontsize);
        fprintf(f,"axes[0, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,
                p_env, PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]],PLOT_label_fontsize);
        ebox_s+=1;
    }
    
    fprintf(f,"axes[0, 0].set_ylabel('Cell cycle events / WT cell cycle',fontsize=%d)\n",PLOT_label_fontsize);
    fprintf(f,"axes[1, 0].set_ylabel('Transitions / WT cell cycle',fontsize=%d)\n",PLOT_label_fontsize);
    
    for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++)
        fprintf(f,"axes[%d+2, 0].set_ylabel('%% time in %s vs. %s',fontsize=%d)\n",
                vi,
                Phen_STATs.Reversible_Transitions_1[vi].c_str(),
                Phen_STATs.Reversible_Transitions_2[vi].c_str(),
                PLOT_label_fontsize);

    ebox_s=0;
    for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
        fprintf(f,"axes[0, %d].legend(prop={'size': %d})\n",ebox_s,(int)(PLOT_label_fontsize*0.75));
        fprintf(f,"axes[1, %d].legend(prop={'size': %d})\n",ebox_s,(int)(PLOT_label_fontsize*0.75));
        for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++)
            fprintf(f,"axes[2+%d, %d].legend(prop={'size': %d})\n",vi,ebox_s,(int)(PLOT_label_fontsize*0.75));
        ebox_s++;
    }
    fprintf(f,"plt.tight_layout()\n");
    //fprintf(f,"plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)\n");
}


void draw_HITs_stat_2input_panel_KOx_CellCycle(FILE *f,
                                   unsigned long int input_index_1,
                                   unsigned long int input_index_2,
                                   string fn2,Environment_and_hits *EnvHit,
                                   MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                   unsigned int row){
    
    fprintf(f,"P_lock = (");
    for(double ph=0;ph<=1.; ph+=1./(double)BOX_KO_OE){
        if(ph>0)  fprintf(f,", ");
        fprintf(f,"%.3lf",ph);
    }
    fprintf(f,")\n");
    fprintf(f,"width = %lf\n",PLOT_WIDTH);
    
    if(Phen_STATs.CellCycleErrors){
        fprintf(f,"plt.tick_params(axis='both', which='major', labelsize=%d)\n",(int)(PLOT_label_fontsize*0.75));
        
        fprintf(f,"fig_cc, axes = plt.subplots(%d, %d, sharex=True, sharey=True)\n", BOX_ENV/2,BOX_ENV/2);
        fprintf(f,"fig_cc.set_size_inches(5*%d,5*%d)\n",BOX_ENV/2,BOX_ENV/2);
   
        fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
            fprintf(f,"\tfor y in range(%d):\n",BOX_ENV/2);
                fprintf(f,"\t\taxes[y, x].tick_params(axis=\"x\", labelsize=%d,labelbottom=True)\n",(int)(PLOT_label_fontsize*0.75));
                fprintf(f,"\t\taxes[y, x].tick_params(axis=\"y\", labelsize=%d,labelleft=True)\n",(int)(PLOT_label_fontsize*0.75));

        
        fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
            fprintf(f,"\tfor y in range(%d):\n",BOX_ENV/2);
                fprintf(f,"\t\taxes[y, x].bar(P_lock, Normal_CC_row%d[2*x+1][2*y+1], 2*width, align='center', label=\"%s\")\n",row,"Normal cycle");
                fprintf(f,"\t\taxes[y, x].bar(P_lock, G2_4N_slip_row%d[2*x+1][2*y+1], 2*width, align='center', bottom=Normal_CC_row%d[2*x+1][2*y+1], label=\"%s\")\n",row,row,"G2->4N");
                fprintf(f,"\t\taxes[y, x].bar(P_lock, Aneupl_4N_row%d[2*x+1][2*y+1], 2*width, align='center', bottom=G2_4N_slip_row%d[2*x+1][2*y+1] + Normal_CC_row%d[2*x+1][2*y+1], label=\"%s\")\n",row,row,row,"Aneuploidy->4N");
                fprintf(f,"\t\taxes[y, x].bar(P_lock, Telophase_4N_row%d[2*x+1][2*y+1], 2*width, align='center', bottom=Aneupl_4N_row%d[2*x+1][2*y+1] + G2_4N_slip_row%d[2*x+1][2*y+1] + Normal_CC_row%d[2*x+1][2*y+1], label=\"%s\")\n",row,row,row,row,"Telophase->4N");

        
        unsigned int ebox_s=0;
        char xlabel[400] = "";
        for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
            if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
            sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
            if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
            else sprintf (xlabel, "%s%s",xlabel,"OE");
        }
        ebox_s=0;
        for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
         //   fprintf(f,"axes[0, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
         //   fprintf(f,"axes[1, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
            fprintf(f,"axes[0, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,
                    p_env, PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize);
            ebox_s+=1;
        }
        ebox_s=0;
        for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
         //   fprintf(f,"axes[0, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
         //   fprintf(f,"axes[1, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
            fprintf(f,"axes[%d, 0].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,
                    p_env, PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize);
            ebox_s+=1;
        }

        fprintf(f,"axes[0, 0].set_ylabel('Cell cycle events / WT cell cycle',fontsize=%d)\n",PLOT_label_fontsize);
    
        fprintf(f,"axes[0, 0].legend(prop={'size': %d})\n",(int)(PLOT_label_fontsize*0.75));
    
        fprintf(f,"plt.tight_layout()\n");
    //fprintf(f,"plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)\n");
        
        
        
        
        
        
        fprintf(f,"P_lock = (");
        for(double ph=0;ph<=1.; ph+=1./(double)BOX_ENV){
            if(ph>0)  fprintf(f,", ");
            fprintf(f,"%.2lf",ph);
        }
        fprintf(f,")\n");
        
        fprintf(f,"sns.set(font_scale = 0.2)\n");
        
        fprintf(f,"fig_cc_hm, axes = plt.subplots(4, %d+2, sharex=True, sharey=True)\n", BOX_KO_OE);
        fprintf(f,"cbar_ax = fig_cc_hm.add_axes([.91, .3, .03, .4])\n");

        
        fprintf(f,"maxnow=max(np.max(np.array(Normal_CC_row%d)),np.max(np.array(G2_4N_slip_row%d)),np.max(np.array(Aneupl_4N_row%d)),np.max(np.array(Telophase_4N_row%d)))\n",row,row,row,row);
        ebox_s=0;
        fprintf(f,"data_now=np.array(Normal_CC_row%d)\n",row);
        
        for(double p_hit=0;p_hit<=1.; p_hit+=1./(double)BOX_KO_OE){
            fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[0, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
            
            char xlabel[400] = "";
            for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                else sprintf (xlabel, "%s%s",xlabel,"OE");
            }
            
            fprintf(f,"axes[0, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
           
            fprintf(f,"axes[0, %d].set_xlabel('%% %s', fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
            if(ebox_s==0)  fprintf(f,"axes[0, %d].set_ylabel(\"Normal cycle \\n %% %s\", fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
            else fprintf(f,"axes[0, %d].set_ylabel(\"\")\n",ebox_s);
           
            ebox_s+=1;
        }
    
        
        ebox_s=0;
        fprintf(f,"data_now=np.array(G2_4N_slip_row%d)\n",row);
        for(double p_hit=0;p_hit<=1.; p_hit+=1./(double)BOX_KO_OE){
            fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[1, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
            
            char xlabel[400] = "";
            for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                else sprintf (xlabel, "%s%s",xlabel,"OE");
            }
            
            fprintf(f,"axes[1, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
            fprintf(f,"axes[1, %d].set_xlabel('%% %s', fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
            if(ebox_s==0) fprintf(f,"axes[1, %d].set_ylabel(\"G2 -> 4N slip \\n %% %s\", fontsize=%d)\n", ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
            else fprintf(f,"axes[1, %d].set_ylabel(\"\")\n",ebox_s);
            ebox_s+=1;
        }
        
        ebox_s=0;
        fprintf(f,"data_now=np.array(Aneupl_4N_row%d)\n",row);
        for(double p_hit=0;p_hit<=1.; p_hit+=1./(double)BOX_KO_OE){
            fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[2, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
            
            char xlabel[400] = "";
            for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                else sprintf (xlabel, "%s%s",xlabel,"OE");
            }
            
            fprintf(f,"axes[2, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
            fprintf(f,"axes[2, %d].set_xlabel(\"Aneuploidy -> 4N slip \\n %% %s\", fontsize=%d)\n", ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
            if(ebox_s==0) fprintf(f,"axes[2, %d].set_ylabel('%% %s', fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
            else fprintf(f,"axes[2, %d].set_ylabel(\"\")\n",ebox_s);
            ebox_s+=1;
        }

        ebox_s=0;
        fprintf(f,"data_now=np.array(Telophase_4N_row%d)\n",row);
        for(double p_hit=0;p_hit<=1.; p_hit+=1./(double)BOX_KO_OE){
            fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[3, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
            fprintf(f,"axes[3, %d].figure.axes[-1].yaxis.label.set_size(%d)\n",ebox_s,PLOT_label_fontsize/10);

            char xlabel[400] = "";
            for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                else sprintf (xlabel, "%s%s",xlabel,"OE");
            }
            
            fprintf(f,"axes[3, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
            fprintf(f,"axes[3, %d].set_xlabel(\"Telophase -> 4N slip \\n %% %s\", fontsize=%d)\n", ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
            if(ebox_s==0) fprintf(f,"axes[3, %d].set_ylabel('%% %s', fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
            else fprintf(f,"axes[3, %d].set_ylabel(\"\")\n",ebox_s);
           
            ebox_s+=1;
        }
        ebox_s-=1;
        //fprintf(f,"axes[3, %d].figure.axes[-1].yaxis.label.set_size(%d)\n",ebox_s-1,PLOT_label_fontsize/10);

        fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[3, %d], square=True, cbar_ax=cbar_ax, cbar_kws={'label': \"Events / WT cell cycle\"})\n", ebox_s,ebox_s);
        
        fprintf(f,"plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6)\n");
        
        
        
        
        fprintf(f,"P_lock = (");
        for(double ph=0;ph<=1.; ph+=1./(double)BOX_ENV){
            if(ph>0)  fprintf(f,", ");
            fprintf(f,"%.2lf",ph);
        }
        fprintf(f,")\n");
        
        fprintf(f,"sns.set(font_scale = 0.2)\n");
        
        fprintf(f,"fig_cc_hm_error, axes = plt.subplots(4, %d+2, sharex=True, sharey=True)\n", BOX_KO_OE);
        fprintf(f,"cbar_ax = fig_cc_hm_error.add_axes([.91, .3, .03, .4])\n");

       // fprintf(f,"fig_cc_hm.set_size_inches(20*4,20*%d)\n",BOX_ENV/2);
        
        fprintf(f,"maxnow=max(np.max(np.array(G2_4N_slip_row%d)),np.max(np.array(Aneupl_4N_row%d)),np.max(np.array(Telophase_4N_row%d)))\n",row,row,row);
       // fprintf(f,"maxnow=max(np.max(np.array(Normal_CC_row%d)),np.max(np.array(G2_4N_slip_row%d)),np.max(np.array(Aneupl_4N_row%d)),np.max(np.array(Telophase_4N_row%d)))\n",row,row,row,row);
        ebox_s=0;
        fprintf(f,"data_now=np.array(Normal_CC_row%d)\n",row);
        
        for(double p_hit=0;p_hit<=1.; p_hit+=1./(double)BOX_KO_OE){
             fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[0, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
            
            char xlabel[400] = "";
            for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                else sprintf (xlabel, "%s%s",xlabel,"OE");
            }
            
            fprintf(f,"axes[0, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
           
            fprintf(f,"axes[0, %d].set_xlabel('%% %s', fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
            if(ebox_s==0)  fprintf(f,"axes[0, %d].set_ylabel(\"Normal cycle \\n %% %s\", fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
            else fprintf(f,"axes[0, %d].set_ylabel(\"\")\n",ebox_s);
           
            ebox_s+=1;
        }
    
        
        ebox_s=0;
        fprintf(f,"data_now=np.array(G2_4N_slip_row%d)\n",row);
        for(double p_hit=0;p_hit<=1.; p_hit+=1./(double)BOX_KO_OE){
           fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[1, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
            
            char xlabel[400] = "";
            for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                else sprintf (xlabel, "%s%s",xlabel,"OE");
            }
            
            fprintf(f,"axes[1, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
            fprintf(f,"axes[1, %d].set_xlabel('%% %s', fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
            if(ebox_s==0) fprintf(f,"axes[1, %d].set_ylabel(\"G2 -> 4N slip \\n %% %s\", fontsize=%d)\n", ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
            else fprintf(f,"axes[1, %d].set_ylabel(\"\")\n",ebox_s);
            ebox_s+=1;
        }
        
        ebox_s=0;
        fprintf(f,"data_now=np.array(Aneupl_4N_row%d)\n",row);
        for(double p_hit=0;p_hit<=1.; p_hit+=1./(double)BOX_KO_OE){
            fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[2, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
            
            char xlabel[400] = "";
            for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                else sprintf (xlabel, "%s%s",xlabel,"OE");
            }
            
            fprintf(f,"axes[2, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
            fprintf(f,"axes[2, %d].set_xlabel(\"Aneuploidy -> 4N slip \\n %% %s\", fontsize=%d)\n", ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
            if(ebox_s==0) fprintf(f,"axes[2, %d].set_ylabel('%% %s', fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
            else fprintf(f,"axes[2, %d].set_ylabel(\"\")\n",ebox_s);
            ebox_s+=1;
        }

        ebox_s=0;
        fprintf(f,"data_now=np.array(Telophase_4N_row%d)\n",row);
        for(double p_hit=0;p_hit<=1.; p_hit+=1./(double)BOX_KO_OE){
            fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[3, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
            fprintf(f,"axes[3, %d].figure.axes[-1].yaxis.label.set_size(%d)\n",ebox_s,PLOT_label_fontsize/10);

            char xlabel[400] = "";
            for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                else sprintf (xlabel, "%s%s",xlabel,"OE");
            }
            
            fprintf(f,"axes[3, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
            fprintf(f,"axes[3, %d].set_xlabel(\"Telophase -> 4N slip \\n %% %s\", fontsize=%d)\n", ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
            if(ebox_s==0) fprintf(f,"axes[3, %d].set_ylabel('%% %s', fontsize=%d)\n",ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
            else fprintf(f,"axes[3, %d].set_ylabel(\"\")\n",ebox_s);
            ebox_s+=1;
        }
        ebox_s-=1;
        //fprintf(f,"axes[3, %d].figure.axes[-1].yaxis.label.set_size(%d)\n",ebox_s-1,PLOT_label_fontsize/10);

        fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[3, %d], square=True, cbar_ax=cbar_ax, cbar_kws={'label': \"Events / WT cell cycle\"})\n", ebox_s,ebox_s);
        
        fprintf(f,"plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6)\n");
    }
    
    
    if(Phen_STATs.Irrev){
        fprintf(f,"P_lock = (");
        for(double ph=0;ph<=1.; ph+=1./(double)BOX_KO_OE){
            if(ph>0)  fprintf(f,", ");
            fprintf(f,"%.3lf",ph);
        }
        fprintf(f,")\n");
        
        fprintf(f,"fig_irrev, axes = plt.subplots(%d, %d, sharex=True, sharey=True)\n", BOX_ENV/2,BOX_ENV/2);
        fprintf(f,"fig_irrev.set_size_inches(5*%d,5*%d)\n",BOX_ENV/2,BOX_ENV/2);
        fprintf(f,"plt.tick_params(axis='both', which='major', labelsize=%d)\n",(int)(PLOT_label_fontsize*0.75));

        fprintf(f,"barGap = %lf\n",(1./(double)BOX_ENV) /(double)Phen_STATs.Irreversible_STATS.size() );
       
        fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
            fprintf(f,"\tfor y in range(%d):\n",BOX_ENV/2);
                fprintf(f,"\t\taxes[y, x].tick_params(axis=\"x\", labelsize=%d,labelbottom=True)\n",(int)(PLOT_label_fontsize*0.75));
                fprintf(f,"\t\taxes[y, x].tick_params(axis=\"y\", labelsize=%d,labelleft=True)\n",(int)(PLOT_label_fontsize*0.75));

        
        fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
            fprintf(f,"\tfor y in range(%d):\n",BOX_ENV/2);
                for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
                    fprintf(f,"\t\tPercent_signal_sh = [a - barGap * %lf + barGap/2. + %d * barGap for a in P_lock]\n",Phen_STATs.Irreversible_STATS.size()/2.,vi);
                    fprintf(f,"\t\taxes[y, x].bar(Percent_signal_sh, rate_rel_CC_IR%d_row%d[2*x+1][2*y+1], 2*width, align='center', label=\"%s\")\n",vi,row,Phen_STATs.Irreversible_STATS[vi].c_str());
                //                //fprintf(f,"\tPercent_signal_sh = [y  + barGap for y in Percent_signal_sh]\n");
                    }
 
        unsigned int ebox_s=0;
        char xlabel[400] = "";
        for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
            if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
            sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
            if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
            else sprintf (xlabel, "%s%s",xlabel,"OE");
        }
        ebox_s=0;
        for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
         //   fprintf(f,"axes[0, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
         //   fprintf(f,"axes[1, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
            fprintf(f,"axes[0, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,
                    p_env, PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize);
            ebox_s+=1;
        }
        ebox_s=0;
        for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
         //   fprintf(f,"axes[0, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
         //   fprintf(f,"axes[1, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
            fprintf(f,"axes[%d, 0].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,
                    p_env, PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize);
            ebox_s+=1;
        }

        fprintf(f,"axes[0, 0].set_ylabel('Transitions / WT cell cycle',fontsize=%d)\n",PLOT_label_fontsize);
    
        fprintf(f,"axes[0, 0].legend(prop={'size': %d})\n",(int)(PLOT_label_fontsize*0.75));
    
        fprintf(f,"plt.tight_layout()\n");
        
        
        
        
        
        fprintf(f,"P_lock = (");
        for(double ph=0;ph<=1.; ph+=1./(double)BOX_ENV){
            if(ph>0)  fprintf(f,", ");
            fprintf(f,"%.2lf",ph);
        }
        fprintf(f,")\n");
        
        fprintf(f,"sns.set(font_scale = 0.2)\n");
        
        fprintf(f,"fig_irrev_hm, axes = plt.subplots(%lu, %d+2)\n", Phen_STATs.Irreversible_STATS.size(),BOX_KO_OE);
      
        //fprintf(f,"maxnow=max(np.max(np.array(G2_4N_slip_row%d)),np.max(np.array(Aneupl_4N_row%d)),np.max(np.array(Telophase_4N_row%d)))\n",row,row,row);
        
        
        if(Phen_STATs.Irreversible_STATS.size()>1){
            fprintf(f,"shax = axes[0, 0].get_shared_x_axes()\n");
            fprintf(f,"shay = axes[0, 0].get_shared_y_axes()\n");
            
            fprintf(f,"for ax in axes[:, :-1].ravel():\n");
                fprintf(f,"\tshax.join(axes[0, 0], ax)\n");
                fprintf(f,"\tshay.join(axes[0, 0], ax)\n");
        }
        else{
            fprintf(f,"shax = axes[0].get_shared_x_axes()\n");
            fprintf(f,"shay = axes[0].get_shared_y_axes()\n");
            fprintf(f,"for ax in axes[:-2].ravel():\n");
                fprintf(f,"\tshax.join(axes[0], ax)\n");
            fprintf(f,"for ax in axes[:-1].ravel():\n");
                fprintf(f,"\tshay.join(axes[0], ax)\n");
        }
        for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
            ebox_s=0;
            fprintf(f,"nmaxnow=np.max(np.array(rate_rel_CC_IR%d_row%d))\n",vi,row);
            fprintf(f,"data_now=np.array(rate_rel_CC_IR%d_row%d)\n",vi,row);
            
//            fprintf(f,"cbar_ax_%d = fig_irrev_hm.add_axes([.91, .3, .03, .4])\n",vi);

            for(double p_hit=0;p_hit<=1.; p_hit+=1./(double)BOX_KO_OE){
               // if(p_hit+2./(double)BOX_KO_OE >=1.) fprintf(f,"cbar_flag=True\n"); else fprintf(f,"cbar_flag=False\n");
             
                if(Phen_STATs.Irreversible_STATS.size()>1)
                     fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[%u, %d], square=True, cbar=False)\n", ebox_s,vi,ebox_s);
                else fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[%d],       square=True, cbar=False)\n", ebox_s,ebox_s);
              
                char xlabel[400] = "";
                for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                    if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                    sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                    if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                    else sprintf (xlabel, "%s%s",xlabel,"OE");
                }
                
                if(Phen_STATs.Irreversible_STATS.size()>1){
                    fprintf(f,"axes[%u, %d].set_title('%.2lf %% %s', fontsize=%d)\n",vi,ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
                   
                    fprintf(f,"axes[%u, %d].set_xlabel('%% %s', fontsize=%d)\n", vi,ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
                    if(ebox_s==0) fprintf(f,"axes[%u, %d].set_ylabel(\"Transitions / WT cell cycle \\n %% %s\", fontsize=%d)\n", vi,ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
                    else fprintf(f,"axes[%u, %d].set_ylabel(\"\")\n",vi,ebox_s);
                }
                else{
                    fprintf(f,"axes[%d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
                   
                    fprintf(f,"axes[%d].set_xlabel('%% %s', fontsize=%d)\n", ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
                    if(ebox_s==0) fprintf(f,"axes[%d].set_ylabel(\"%s / WT cell cycle \\n %% %s\", fontsize=%d)\n", ebox_s,
                            Phen_STATs.Irreversible_STATS[vi].c_str(),
                            PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
                    else fprintf(f,"axes[%d].set_ylabel(\"\")\n",ebox_s);
                   
                }
                ebox_s+=1;
            }
            
            fprintf(f,"cbar = np.array([np.flip(np.arange(0, maxnow, maxnow/%lf))]).transpose()\n",(float)BOX_ENV);
           // fprintf(f,"cbar = np.array([np.arange(0, maxnow, maxnow/%lf)]).transpose()\n",(float)BOX_KO_OE);
            //fprintf(f,"cbar = np.array([np.arange(0, maxnow, maxnow/20.)]).transpose()\n");
            if(Phen_STATs.Irreversible_STATS.size()>1)
                 fprintf(f,"sns.heatmap(cbar, ax=axes[%lu, %d], square=True, cbar=False)\n",  Phen_STATs.Irreversible_STATS.size()-1,ebox_s);
            else fprintf(f,"sns.heatmap(cbar, ax=axes[%d], square=True, cbar=True)\n", ebox_s);
        }
        fprintf(f,"plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6)\n");
        
    }
     
    if(Phen_STATs.Rev){
        
        for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
            fprintf(f,"P_lock = (");
            for(double ph=0;ph<=1.; ph+=1./(double)BOX_KO_OE){
                if(ph>0)  fprintf(f,", ");
                fprintf(f,"%.3lf",ph);
            }
            fprintf(f,")\n");
            
            fprintf(f,"fig_rev_%d, axes = plt.subplots(%d, %d, sharex=True, sharey=True)\n",vi, BOX_ENV/2,BOX_ENV/2);
            fprintf(f,"fig_rev_%d.set_size_inches(5*%d,5*%d)\n",vi, BOX_ENV/2,BOX_ENV/2);
            fprintf(f,"plt.tick_params(axis='both', which='major', labelsize=%d)\n",(int)(PLOT_label_fontsize*0.75));

            fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
                fprintf(f,"\tfor y in range(%d):\n",BOX_ENV/2);
                    fprintf(f,"\t\taxes[y, x].tick_params(axis=\"x\", labelsize=%d,labelbottom=True)\n",(int)(PLOT_label_fontsize*0.75));
                    fprintf(f,"\t\taxes[y, x].tick_params(axis=\"y\", labelsize=%d,labelleft=True)\n",(int)(PLOT_label_fontsize*0.75));

            
            fprintf(f,"for x in range(%d):\n",BOX_ENV/2);
                fprintf(f,"\tfor y in range(%d):\n",BOX_ENV/2);
                    fprintf(f,"\t\taxes[y, x].bar(Percent_signal_sh, Time_in_1_Rpair%d_row%d[2*x+1][2*y+1], width, align='center', label=\"%s\")\n",vi,row,Phen_STATs.Reversible_Transitions_1[vi].c_str());
                    fprintf(f,"\t\taxes[y, x].bar(Percent_signal_sh, Time_in_2_Rpair%d_row%d[2*x+1][2*y+1], width, align='center',  bottom=Time_in_1_Rpair%d_row%d[2*x+1][2*y+1], label=\"%s\")\n",vi,row,vi,row,Phen_STATs.Reversible_Transitions_2[vi].c_str());
            //fprintf(f,"\t\tPercent_signal_sh = [y  + barGap for y in Percent_signal_sh]\n");
            
            
            unsigned int ebox_s=0;
            char xlabel[400] = "";
            for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                else sprintf (xlabel, "%s%s",xlabel,"OE");
            }
            ebox_s=0;
            for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
             //   fprintf(f,"axes[0, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
             //   fprintf(f,"axes[1, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
                fprintf(f,"axes[0, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,
                        p_env, PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize);
                ebox_s+=1;
            }
            ebox_s=0;
            for(double p_env=1./(double)BOX_ENV;p_env<1.; p_env+=2./(double)BOX_ENV){
             //   fprintf(f,"axes[0, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
             //   fprintf(f,"axes[1, %d].set_xlabel('%s',fontsize=%d)\n",ebox_s,xlabel,PLOT_label_fontsize);
                fprintf(f,"axes[%d, 0].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,
                        p_env, PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize);
                ebox_s+=1;
            }

            fprintf(f,"axes[0, 0].set_ylabel('%% time in %s vs. %s',fontsize=%d)\n",
                    Phen_STATs.Reversible_Transitions_1[vi].c_str(),
                    Phen_STATs.Reversible_Transitions_2[vi].c_str(),
                    PLOT_label_fontsize);
            
            fprintf(f,"axes[0, 0].legend(prop={'size': %d})\n",(int)(PLOT_label_fontsize*0.75));
        
            fprintf(f,"plt.tight_layout()\n");
        
        
            
            fprintf(f,"P_lock = (");
            for(double ph=0;ph<=1.; ph+=1./(double)BOX_ENV){
                if(ph>0)  fprintf(f,", ");
                fprintf(f,"%.2lf",ph);
            }
            fprintf(f,")\n");
            
            fprintf(f,"sns.set(font_scale = 0.2)\n");
            
            fprintf(f,"fig_rev_hm_%d, axes = plt.subplots(2, %d+2, sharey=True)\n", vi,BOX_KO_OE);
            fprintf(f,"cbar_ax = fig_rev_hm_%d.add_axes([.91, .3, .03, .4])\n",vi);
            
            
            fprintf(f,"maxnow = 0\n");
            fprintf(f,"maxnow=max(maxnow,np.max(np.array(Time_in_1_Rpair%d_row%d)))\n",vi,row);
            fprintf(f,"maxnow=max(maxnow,np.max(np.array(Time_in_2_Rpair%d_row%d)))\n",vi,row);
                
            ebox_s=0;
            fprintf(f,"data_now=np.array(Time_in_1_Rpair%d_row%d)\n",vi,row);
            
            for(double p_hit=0;p_hit<1.; p_hit+=1./(double)BOX_KO_OE){
                fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[0, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
                
                char xlabel[400] = "";
                for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                    if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                    sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                    if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                    else sprintf (xlabel, "%s%s",xlabel,"OE");
                }
                
                fprintf(f,"axes[0, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
               
                fprintf(f,"axes[0, %d].set_xlabel('%% %s', fontsize=%d)\n", ebox_s,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
                if(ebox_s==0)  fprintf(f,"axes[0, %d].set_ylabel(\"Time in %s \\n %% %s\", fontsize=%d)\n", ebox_s,Phen_STATs.Reversible_Transitions_1[vi].c_str(), PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
                else fprintf(f,"axes[0, %d].set_ylabel(\"\")\n",ebox_s);
                
                ebox_s+=1;
            }
                 
            fprintf(f,"cbar = np.array([np.flip(np.arange(0, maxnow, maxnow/%lf))]).transpose()\n",(float)BOX_ENV);
            fprintf(f,"if maxnow > 0:\n\tsns.heatmap(cbar, ax=axes[0, %d], square=True, cbar=True, cbar_kws={'label': \"%% Time in Cell state\"})\n",ebox_s);
   
            ebox_s=0;
            fprintf(f,"data_now=np.array(Time_in_2_Rpair%d_row%d)\n",vi,row);
            
            for(double p_hit=0;p_hit<1.; p_hit+=1./(double)BOX_KO_OE){
                fprintf(f,"sns.heatmap(data_now[:,:,%d], vmin=0, vmax=maxnow, yticklabels=P_lock, xticklabels=P_lock, ax=axes[1, %d], square=True, cbar=False)\n", ebox_s,ebox_s);
                
                char xlabel[400] = "";
                for(unsigned int i=1;i<=EnvHit->Scan_Hits->N;i++){
                    if(i>1)  sprintf (xlabel, "%s%s ",xlabel,", ");
                    sprintf (xlabel, "%s%s ",xlabel, PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[i]]);
                    if (EnvHit->Scan_Hits->B[i]==0)  sprintf (xlabel, "%s%s",xlabel,"KD");
                    else sprintf (xlabel, "%s%s",xlabel,"OE");
                }
                
                fprintf(f,"axes[1, %d].set_title('%.2lf %% %s', fontsize=%d)\n",ebox_s,p_hit,xlabel,PLOT_label_fontsize/10);
               
                fprintf(f,"axes[1, %d].set_xlabel('%% %s', fontsize=%d)\n", ebox_s, PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]],PLOT_label_fontsize/10);
                if(ebox_s==0)  fprintf(f,"axes[1, %d].set_ylabel(\"Time in %s \\n %% %s\", fontsize=%d)\n", ebox_s,Phen_STATs.Reversible_Transitions_2[vi].c_str(), PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PLOT_label_fontsize/10);
                else fprintf(f,"axes[1, %d].set_ylabel(\"\")\n",ebox_s);
                ebox_s+=1;
            }
            
            fprintf(f,"cbar = np.array([np.flip(np.arange(0, maxnow, maxnow/%lf))]).transpose()\n",(float)BOX_ENV);
            fprintf(f,"if maxnow > 0:\n\tsns.heatmap(cbar, ax=axes[1, %d], square=True, cbar=True, cbar_kws={'label': \"%% Time in Cell state\"})\n",ebox_s);
         
        }
        fprintf(f,"plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6)\n");
     
    }
    
      
  
    
}

void Python_Grapsh_NonSaturating_Stats_FIX_KDOE_fn1Env(unsigned long int input_index,
                                                Environment_and_hits *EnvHit,
                                                longint_ARRAY *l,
                                                MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                int run_assync,
                                                bool *cc_stats,bool *irrev_stats,bool *rev_stats){
    string fn2,fneps,fmain;
    char fn[400],fn_comm[400],fnpdf[400];
    longint_ARRAY *lin;
    FILE *f;
    
    lin=new longint_ARRAY();
    lin->add_element(input_index);
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
   
    if(!check_all_files(input_index,fn2,EnvHit,l,PD,run_assync,cc_stats,irrev_stats,rev_stats)) return;
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s/_FIGURES\n",fn2.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    
    for(int j=1;j<=l->N;j++){
        sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
        for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
            if(e!=input_index)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
        sprintf(fn,"%s___SCAN-%s___",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]]);
        if(EnvHit->BKR_Hits->N>0) {
            sprintf(fn,"%s___BKG-",fn);
            for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
        }
        if(EnvHit->Scan_Hits->N>0){
            sprintf(fn,"%s___JointScan-",fn);
            for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
                if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
                else sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
        }
        sprintf(fnpdf,"%s__A%d__BOX_%d",fn,(int)l->A[j],BOX_KO_OE);
        sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);
        
        f=fopen(fn,"w");
        initialize_Python_script(f,PD->BRN->path);
        fprintf(f, "molec = '%s'\n",PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]]);
     
        fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
        fn2 = get_Figures_directory_now(fmain,lin,PD);
        
        if(Phen_STATs.CellCycleErrors) read_CC_datafiles_in_Python_script(f,input_index,fn2,EnvHit,l,j,PD);
        if(Phen_STATs.Irrev) read_Irrev_datafiles_in_Python_script(f,input_index,fn2,EnvHit,l,j,PD);
        if(Phen_STATs.Rev) read_Rev_datafiles_in_Python_script(f,input_index,fn2,EnvHit,l,j,PD);

        if(EnvHit->Scan_Hits->N==0){
            // 1 row, 3 graphs for the 3 kids of results; x axis in each is input_index probs.
            // data averaged over all intake!
            if(Phen_STATs.CellCycleErrors) collapse_CC_data_when_KO_is_absent_in_PythonScript(f, fn2.c_str(), 0);
            if(Phen_STATs.Irrev) collapse_Irrev_data_when_KO_is_absent_in_PythonScript(f, fn2.c_str(), 0);
            if(Phen_STATs.Rev) collapse_Rev_data_when_KO_is_absent_in_PythonScript(f, fn2.c_str(), 0);

            draw_NOHITs_stat_1input_panel(f,input_index,fn2,EnvHit,l,j,PD,0);
            fneps = (string)fnpdf + ".pdf";
            fprintf(f,"fig.savefig('%s')\n",fneps.c_str());
        }
        else{
             // 3 rows for the 3 kids of results
            // each row: input_index goes up for each graph, x axis inside graphs is KO/OE p_now.
            // all data counts
            draw_HITs_stat_1input_panel(f,input_index,fn2,EnvHit,PD,0);
            fneps = (string)fnpdf + ".pdf";
            fprintf(f,"fig.savefig('%s')\n",fneps.c_str());
       }
        fclose(f);
    }
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    
    for(int j=1;j<=l->N;j++){
        sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
        for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
            if(e!=input_index)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
        sprintf(fn,"%s___SCAN-%s___",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]]);
        if(EnvHit->BKR_Hits->N>0) {
            sprintf(fn,"%s___BKG-",fn);
            for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
        }
        if(EnvHit->Scan_Hits->N>0){
            sprintf(fn,"%s___JointScan-",fn);
            for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
                if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
                else sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
        }
        sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);
    
        //sprintf(fn_comm,"/opt/homebrew/bin//python3 %s\n",fn);
        sprintf(fn_comm,"python3 %s\n",fn);
        // printf(fn_comm); getchar();
        system(fn_comm);
    }
    delete lin; lin=NULL;
    
  //  getchar();
}

void Python_Grapsh_NonSaturating_Stats_Scan_1KDOE_fn1Env(unsigned long int input_index,
                                                Environment_and_hits *EnvHit,
                                                 longint_ARRAY *l,
                                                MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                int run_assync,
                                                bool *cc_stats,bool *irrev_stats,bool *rev_stats){
    string fn2,fneps,fmain;
    char fn[400],fn_comm[400],fnpdf[400];
    longint_ARRAY *lin;
    FILE *f;
    
    lin=new longint_ARRAY();
    lin->add_element(input_index);
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
   
    if(!check_all_files(input_index,fn2,EnvHit,l,PD,run_assync,cc_stats,irrev_stats,rev_stats)) return;
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s/_FIGURES\n",fn2.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    
    for(int j=1;j<=l->N;j++){
        sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
        for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
            if(e!=input_index)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
        sprintf(fn,"%s___SCAN_Xaxis-%s___",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]]);
        if(EnvHit->BKR_Hits->N>0) {
            sprintf(fn,"%s___BKG-",fn);
            for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
        }
        if(EnvHit->Scan_Hits->N>0){
            sprintf(fn,"%s___JointScan-",fn);
            for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
                if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
                else sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
        }
        sprintf(fnpdf,"%s__A%d__BOX_%d",fn,(int)l->A[j],BOX_KO_OE);
        sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);
        
        f=fopen(fn,"w");
        initialize_Python_script(f,PD->BRN->path);
        fprintf(f, "molec = '%s'\n",PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]]);
     
        fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
        fn2 = get_Figures_directory_now(fmain,lin,PD);
        
        if(Phen_STATs.CellCycleErrors) read_CC_datafiles_in_Python_script(f,input_index,fn2,EnvHit,l,j,PD);
        if(Phen_STATs.Irrev) read_Irrev_datafiles_in_Python_script(f,input_index,fn2,EnvHit,l,j,PD);
        if(Phen_STATs.Rev) read_Rev_datafiles_in_Python_script(f,input_index,fn2,EnvHit,l,j,PD);

        if(EnvHit->Scan_Hits->N==0){
            printf("This function should not be called; this situation is handled by the old one; check inut and mandate a KO/OE molecule!\n"); getchar(); return;
        }
        else{
            //draw_HITs_stat_1input_panel(f,input_index,fn2,EnvHit,PD,0);
            draw_HITs_stat_1input_panel_inverted(f,input_index,fn2,EnvHit,PD,0);
            fneps = (string)fnpdf + ".pdf";
            fprintf(f,"fig.savefig('%s')\n",fneps.c_str());
       }
        fclose(f);
    }
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    
    for(int j=1;j<=l->N;j++){
        sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
        for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
            if(e!=input_index)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
        sprintf(fn,"%s___SCAN_Xaxis-%s___",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]]);
        if(EnvHit->BKR_Hits->N>0) {
            sprintf(fn,"%s___BKG-",fn);
            for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
        }
        if(EnvHit->Scan_Hits->N>0){
            sprintf(fn,"%s___JointScan-",fn);
            for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
                if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
                else sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
        }
        sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);
    
        //sprintf(fn_comm,"/opt/homebrew/bin//python3 %s\n",fn);
        sprintf(fn_comm,"python3 %s\n",fn);
        system(fn_comm);
    }
    delete lin; lin=NULL;
    
  //  getchar();
}



void Python_Grapsh_NonSaturating_Stats_Scan_1KDOE_2Env(unsigned long int input_index_1,
                                                       unsigned long int input_index_2,
                                                       Environment_and_hits *EnvHit,
                                                       longint_ARRAY *l,
                                                       MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                       int run_assync,
                                                       bool *cc_stats,bool *irrev_stats,bool *rev_stats){
    string fn2,fneps,fmain;
    char fn[400],fn_comm[400],fnpdf[400];
    longint_ARRAY *lin;
    FILE *f;
    
    lin=new longint_ARRAY();
    lin->add_element(input_index_1);
    lin->add_element(input_index_2);
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
   
    if(!check_all_files_2D(input_index_1,input_index_2,fn2,EnvHit,l,PD,run_assync,cc_stats,irrev_stats,rev_stats)) return;
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s/_FIGURES\n",fn2.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    
    for(int j=1;j<=l->N;j++){
        sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
        for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
            if((e!=input_index_1)&&(e!=input_index_2))
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
        sprintf(fn,"%s___SCAN_Xaxis-%s___SCAN_Yaxis-%s___",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]]);
        if(EnvHit->BKR_Hits->N>0) {
            sprintf(fn,"%s___BKG-",fn);
            for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
        }
        if(EnvHit->Scan_Hits->N>0){
            sprintf(fn,"%s___JointScan-",fn);
            for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
                if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
                else sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
        }
        sprintf(fnpdf,"%s__A%d__BOX_%d",fn,(int)l->A[j],BOX_KO_OE);
        sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);
        
        f=fopen(fn,"w");
        initialize_Python_script(f,PD->BRN->path);
        fprintf(f, "molec1 = '%s'\n",PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]]);
        fprintf(f, "molec2 = '%s'\n",PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]]);
     
        fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
        fn2 = get_Figures_directory_now(fmain,lin,PD);
        
        if(Phen_STATs.CellCycleErrors) read_CC_datafiles_in_Python_script_2D(f,input_index_1,input_index_2,fn2,EnvHit,l,j,PD);
        if(Phen_STATs.Irrev) read_Irrev_datafiles_in_Python_script_2D(f,input_index_1,input_index_2,fn2,EnvHit,l,j,PD);
        if(Phen_STATs.Rev) read_Rev_datafiles_in_Python_script_2D(f,input_index_1,input_index_2,fn2,EnvHit,l,j,PD);

        if(EnvHit->Scan_Hits->N==0){
            printf("This function should not be called; this situation is handled by the old one; check inut and mandate a KO/OE molecule!\n"); getchar(); return;
        }
        else{
            //draw_HITs_stat_1input_panel(f,input_index,fn2,EnvHit,PD,0);
           // draw_HITs_stat_1input_panel_inverted(f,input_index,fn2,EnvHit,PD,0);
            draw_HITs_stat_2input_panel_KOx_CellCycle(f,input_index_1,input_index_2,fn2,EnvHit,PD,0);
            fneps = (string)fnpdf + ".pdf";
         //   fprintf(f,"fig.savefig('%s')\n",fneps.c_str());
            
            fprintf(f,"def save_multi_image(filename):\n\tpp = PdfPages(filename)\n\tfig_nums = plt.get_fignums()\n");
            fprintf(f,"\tfigs = [plt.figure(n) for n in fig_nums]\n\tfor fig in figs:\n\t\tfig.savefig(pp, format='pdf')\n\tpp.close()\n");

            fprintf(f,"save_multi_image(\"%s\")",fneps.c_str());
       }
        fclose(f);
    }
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    
    for(int j=1;j<=l->N;j++){
        sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
        for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
            if((e!=input_index_1)&&(e!=input_index_2))
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
        sprintf(fn,"%s___SCAN_Xaxis-%s___SCAN_Yaxis-%s___",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_1]],PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index_2]]);
        if(EnvHit->BKR_Hits->N>0) {
            sprintf(fn,"%s___BKG-",fn);
            for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
        }
        if(EnvHit->Scan_Hits->N>0){
            sprintf(fn,"%s___JointScan-",fn);
            for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
                if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
                else sprintf(fn,"%s_%s-KD",fn,PD->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
        }
        sprintf(fnpdf,"%s__A%d__BOX_%d",fn,(int)l->A[j],BOX_KO_OE);
        sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);
    
        //sprintf(fn_comm,"/opt/homebrew/bin//python3 %s\n",fn);
        sprintf(fn_comm,"python3 %s\n",fn);
        system(fn_comm);
    }
    delete lin; lin=NULL;
    
  //  getchar();
}

void Python_Grapsh_NonSaturating_Stats_1_env_KOOE_scan(unsigned long int input_index,
                                                       Environment_and_hits *EnvHit,
                                                       longint_ARRAY *l,
                                                       MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                       int run_assync,
                                                       bool *cc_stats,bool *irrev_stats,bool *rev_stats){
           string fn2,fneps,fmain;
           char fn[400];
           longint_ARRAY *lin;

           
           lin=new longint_ARRAY();
           lin->add_element(input_index);
           fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
           fn2 = get_Figures_directory_now(fmain,lin,PD);
          
           //if(!check_all_files_1_env_KOOE_scan(input_index,fn2,EnvHit,l,PD,run_assync,cc_stats,irrev_stats,rev_stats)) return;
           if(!check_all_files(input_index,fn2,EnvHit,l,PD,run_assync,cc_stats,irrev_stats,rev_stats)) return;
    
           fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
           fn2 = get_Figures_directory_now(fmain,lin,PD);
           sprintf(fn,"mkdir %s/_FIGURES\n",fn2.c_str());
           system(fn);
           
           fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
           fn2 = get_Figures_directory_now(fmain,lin,PD);
           /*
           for(int j=1;j<=l->N;j++){
               sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
               for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
                   if(e!=input_index)
                       sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
               sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);

               f=fopen(fn,"w");
               initialize_Python_script(f,PD->BRN->path);
               fprintf(f, "molec = '%s'\n",PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_index]]);
            
               fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
               fn2 = get_Figures_directory_now(fmain,lin,PD);
               
               if(Phen_STATs.CellCycleErrors) read_CC_datafiles_in_Python_script(f,input_index,fn2,EnvHit,l,j,PD);
               if(Phen_STATs.Irrev) read_Irrev_datafiles_in_Python_script(f,input_index,fn2,EnvHit,l,j,PD);
               if(Phen_STATs.Rev) read_Rev_datafiles_in_Python_script(f,input_index,fn2,EnvHit,l,j,PD);

               if(EnvHit->Scan_Hits->N==0){
                   // 1 row, 3 graphs for the 3 kids of results; x axis in each is input_index probs.
                   // data averaged over all intake!
                   if(Phen_STATs.CellCycleErrors) collapse_CC_data_when_KO_is_absent_in_PythonScript(f, fn2.c_str(), 0);
                   if(Phen_STATs.Irrev) collapse_Irrev_data_when_KO_is_absent_in_PythonScript(f, fn2.c_str(), 0);
                   if(Phen_STATs.Rev) collapse_Rev_data_when_KO_is_absent_in_PythonScript(f, fn2.c_str(), 0);

                   draw_NOHITs_stat_1input_panel(f,input_index,fn2,EnvHit,l,j,PD,0);
                   fneps = (string)fn + ".pdf";
                   fprintf(f,"fig.savefig('%s')\n",fneps.c_str());
               }
               else{
                    // 3 rows for the 3 kids of results
                   // each row: input_index goes up for each graph, x axis inside graphs is KO/OE p_now.
                   // all data counts
                   draw_HITs_stat_1input_panel(f,input_index,fn2,EnvHit,PD,0);
                   fneps = (string)fn + ".pdf";
                   fprintf(f,"fig.savefig('%s')\n",fneps.c_str());
              }
               fclose(f);
           }
           
           fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
           fn2 = get_Figures_directory_now(fmain,lin,PD);
           
           for(int j=1;j<=l->N;j++){
               sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
               for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
                   if(e!=input_index)
                       sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
               sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);
           
               sprintf(fn_comm,"/opt/homebrew/bin//python3 %s\n",fn);
               system(fn_comm);
           }
           delete lin; lin=NULL;
           */
         //  getchar();
       }

void Python_Grapsh_NonSaturating_Stats_2_env_KOOE_scan (unsigned long int input_x,unsigned long int input_y,
                                                Environment_and_hits *EnvHit,
                                                longint_ARRAY *l,
                                                MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                int run_assync,
                                                bool *cc_stats,bool *irrev_stats,bool *rev_stats){
    string fn2,fneps,fmain;
    char fn[400],fn_comm[400];
    longint_ARRAY *lin;
   // FILE *f;
    
    lin=new longint_ARRAY();
    lin->add_element(input_x);
    lin->add_element(input_y);
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
   
    if(!check_all_files_2D(input_x,input_y,fn2,EnvHit,l,PD,run_assync,cc_stats,irrev_stats,rev_stats)) return;
    
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    sprintf(fn,"mkdir %s/_FIGURES\n",fn2.c_str());
    system(fn);
    
    fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    /*
    for(int j=1;j<=l->N;j++){
        sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
        for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
            if((e!=input_x)&&(e!=input_y))
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
        sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);

        f=fopen(fn,"w");
        initialize_Python_script(f,PD->BRN->path);
        fprintf(f, "molec_x = '%s'\n",PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_x]]);
        fprintf(f, "molec_y = '%s'\n",PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[input_y]]);
     
        fmain=get_Figures_directory_main(EnvHit,PD,run_assync);
        fn2 = get_Figures_directory_now(fmain,lin,PD);
        
        if(Phen_STATs.CellCycleErrors) read_CC_datafiles_in_Python_script_2D(f,input_x,input_y,fn2,EnvHit,l,j,PD);
        if(Phen_STATs.Irrev) read_Irrev_datafiles_in_Python_script_2D(f,input_x,input_y,fn2,EnvHit,l,j,PD);
        if(Phen_STATs.Rev) read_Rev_datafiles_in_Python_script_2D(f,input_x,input_y,fn2,EnvHit,l,j,PD);

        if(EnvHit->Scan_Hits->N==0){
            // 1 row, 3 graphs for the 3 kids of results; x axis in each is input_index probs.
            // data averaged over all intake!
            if(Phen_STATs.CellCycleErrors) collapse_CC_data_when_KO_is_absent_in_PythonScript_2D(f, fn2.c_str(), 0);
            if(Phen_STATs.Irrev) collapse_Irrev_data_when_KO_is_absent_in_PythonScript_2D(f, fn2.c_str(), 0);
            if(Phen_STATs.Rev) collapse_Rev_data_when_KO_is_absent_in_PythonScript_2D(f, fn2.c_str(), 0);

            draw_NOHITs_stat_2_inputs_panel(f,input_x,input_y,fn2,EnvHit,l,j,PD,0);
            fneps = (string)fn + ".pdf";
            fprintf(f,"fig.savefig('%s')\n",fneps.c_str());
        }
        else{
             // 3 rows for the 3 kids of results
            // each row: input_index goes up for each graph, x axis inside graphs is KO/OE p_now.
            // all data counts
            draw_HITs_stat_2_inputs_panel(f,input_x,input_y,fn2,EnvHit,PD,0);
            fneps = (string)fn + ".pdf";
            fprintf(f,"fig.savefig('%s')\n",fneps.c_str());
       }
        fclose(f);
    }
    */
    fmain=get_Figures_directory_main_system(EnvHit,PD,run_assync);
    fn2 = get_Figures_directory_now(fmain,lin,PD);
    
    for(int j=1;j<=l->N;j++){
        sprintf(fn,"%s/_FIGURES/BAR-GRAPHS-",fn2.c_str());
        for(unsigned int e=1;e<=PD->Coupled->BRN->active_inputs->N;e++)
            if((e!=input_x)&&(e!=input_y))
                sprintf(fn,"%s_%s-%.2lf",fn,PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
       /*
        if(EnvHit->BKR_Hits->N>0) {
            sprintf(fn,"%s___BKG-",fn);
            for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++)
                sprintf(fn,"%s_%s-%.2lf",fn,PD_A->Coupled->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[t]],EnvHit->p_hits[t]);
        }
        if(EnvHit->Scan_Hits->N>0){
            sprintf(fn,"%s___JointScan-",fn);
            for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++)
                if(EnvHit->Scan_Hits->B[t]==0) sprintf(fn,"%s_%s-KD",fn,PD_A->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
                else sprintf(fn,"%s_%s-KD",fn,PD_A->Coupled->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[t]]);
        }
        else sprintf(fn,"%s___WType",fn);
        */
        
        sprintf(fn,"%s__A%d__BOX_%d.py",fn,(int)l->A[j],BOX_KO_OE);
    
        //sprintf(fn_comm,"/opt/homebrew/bin//python3 %s\n",fn);
        sprintf(fn_comm,"python3 %s\n",fn);
        system(fn_comm);
    }
    delete lin; lin=NULL;
    
  //  getchar();
}

#endif /* Python_Figures_module_h */
