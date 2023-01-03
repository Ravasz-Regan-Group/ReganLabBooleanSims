//
//  Analyze_Boolean_Model.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 11/21/17.
//  Copyright © 2017 Regan_Group. All rights reserved.
//

#ifndef Analyze_Boolean_Model_h
#define Analyze_Boolean_Model_h
#define MAX_ARREST 100000
#define TIMESTEPS_DEFAULT 400

int DRAW_PDFS =0;

struct Cell_cycle_stats_in_run {
    unsigned int    N_cycles, N_deaths,
    N_endoredupl_T,N_endoredupl_G2,N_endoredupl_Aneup,Cell_cycle_stats_in_run;
    
    unsigned long int Time_G1,Time_R,Time_G2,Time_M,Time_Telo;
    unsigned long int N_G1,N_G2,N_R,N_M,N_T;
    unsigned long long int   T_live,MAX_INTERVAL=0;
    unsigned int    Dist_G1[MAX_ARREST+1],
    Dist_G2[MAX_ARREST+1],
    Dist_Mitosis[MAX_ARREST+1],
    Dist_Replication[MAX_ARREST+1],
    Dist_Telophase[MAX_ARREST+1];
};



void draw_timestep_MOD(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int t,
                       FILE *f,double y0,int W,int T_MAX,int TW,int NameBuff){
    
    long int i;
    int poz=0;
    
    for(int m=1;m<=PD->BRN->Module_Order->N;m++){
        if((PD->BRN->Module_Order->A[m]>0)&&(PD->BRN->Module_Order->A[m]<=PD->BRN->MDAT->Module_NR)){
            // printf("Module %d\t",m);
            // printf("Going for it %ld \n",PD->BRN->Module_Order->A[m]);
            // printf("\tNodes: %ld \n",PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N);
            for(int nn=1;nn<=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N;nn++){
                i=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->A[nn];
                //   printf("\t nn=%d node = %ld %s\n",nn,i,PD->BRN->MDAT->Node_name[i]);
                if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
                EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                EpsSetRgb(f,0,0,0);
                EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
                poz++;
            }
        }
        else{
            i=PD->BRN->Module_Order->B[m];
            if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
            EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
            EpsSetRgb(f,0,0,0);
            EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
            poz++;
        }
        poz++;
    }
}

void draw_timestep_ASSYNC_Average(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, double *state_average, int t,
                       FILE *f,double y0,int W,int T_MAX,int TW,int NameBuff){
    
    long int i;
    int poz=0;
    
    for(int m=1;m<=PD->BRN->Module_Order->N;m++){
        if((PD->BRN->Module_Order->A[m]>0)&&(PD->BRN->Module_Order->A[m]<=PD->BRN->MDAT->Module_NR)){
            // printf("Module %d\t",m);
            // printf("Going for it %ld \n",PD->BRN->Module_Order->A[m]);
            // printf("\tNodes: %ld \n",PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N);
            for(int nn=1;nn<=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N;nn++){
                i=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->A[nn];
                //   printf("\t nn=%d node = %ld %s\n",nn,i,PD->BRN->MDAT->Node_name[i]);
               // if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
                if (state_average[i-1]<0.5) EpsSetRgb(f,0,0,1-state_average[i-1]);
                    else EpsSetRgb(f,state_average[i-1],state_average[i-1]/2.,0);
                EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                EpsSetRgb(f,0,0,0);
                EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
                poz++;
            }
        }
        else{
            i=PD->BRN->Module_Order->B[m];
            if (state_average[i-1]<0.5) EpsSetRgb(f,0,0,1-state_average[i-1]);
                else EpsSetRgb(f,state_average[i-1],state_average[i-1]/2.,0);
            EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
            EpsSetRgb(f,0,0,0);
            EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
            poz++;
        }
        poz++;
    }
}

void draw_timestep_MOD(longint_ARRAY_pair *Hits,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int t,
                       FILE *f,double y0,int W,int T_MAX,int TW,int NameBuff){
    
    long int i;
    int poz=0;
    
    for(int m=1;m<=PD->BRN->Module_Order->N;m++){
        if((PD->BRN->Module_Order->A[m]>0)&&(PD->BRN->Module_Order->A[m]<=PD->BRN->MDAT->Module_NR)){
            // printf("Module %d\t",m);
            // printf("Going for it %ld \n",PD->BRN->Module_Order->A[m]);
            // printf("\tNodes: %ld \n",PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N);
            for(int nn=1;nn<=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N;nn++){
                i=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->A[nn];
                //   printf("\t nn=%d node = %ld %s\n",nn,i,PD->BRN->MDAT->Node_name[i]);
                if(Hits->check_for_element_A(i)){
                    if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,0); else EpsSetRgb(f,1,1,0);
                }
                else { if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);}
                EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                EpsSetRgb(f,0,0,0);
                EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
                poz++;
            }
        }
        else{
            i=PD->BRN->Module_Order->B[m];
            if(Hits->check_for_element_A(i)){
                if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,0); else EpsSetRgb(f,1,1,1);
            }
            else { if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);}
            EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
            EpsSetRgb(f,0,0,0);
            EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
            poz++;
        }
        poz++;
    }
}

void draw_timestep_MOD(Boolean_Dynamics_TRACKER *D,int t,
                       FILE *f,double y0,int W,int T_MAX,int TW,int NameBuff){
    
    long int i;
    int poz=0;
    
    if(D->BRN->Module_Order!=NULL){
        for(int m=1;m<=D->BRN->Module_Order->N;m++){
            if(D->BRN->Module_Order->A[m]>0){
                for(int nn=1;nn<=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[m]]->N;nn++){
                    i=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[m]]->A[nn];
                    if (D->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
                    EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                    EpsSetRgb(f,0,0,0);
                    EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
                    poz++;
                }
            }
            else{
                i=D->BRN->Module_Order->B[m];
                if (D->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
                EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                EpsSetRgb(f,0,0,0);
                EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
                poz++;
            }
            poz++;
        }
    }
    else{
        for(int nn=1;nn<=D->BRN->N;nn++){
            i=nn;
            if (D->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
            EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
            EpsSetRgb(f,0,0,0);
            EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
            poz++;
        }
    }
}

void draw_timestep_MOD(Boolean_Dynamics *D,int t,
                       FILE *f,double y0,int W,int T_MAX,int TW,int NameBuff){
    
    long int i;
    int poz=0;
    
    if(D->BRN->Module_Order!=NULL){
        for(int m=1;m<=D->BRN->Module_Order->N;m++){
            if(D->BRN->Module_Order->A[m]>0){
                for(int nn=1;nn<=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[m]]->N;nn++){
                    i=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[m]]->A[nn];
                    if (D->s[i]==0) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
                    EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                    EpsSetRgb(f,0,0,0);
                    EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
                    poz++;
                }
            }
            else{
                i=D->BRN->Module_Order->B[m];
                if (D->s[i]==0) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
                EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                EpsSetRgb(f,0,0,0);
                EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
                poz++;
            }
            poz++;
        }
    }
    else{
        for(int nn=1;nn<=D->BRN->N;nn++){
            i=nn;
            if (D->s[i]==0) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
            EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
            EpsSetRgb(f,0,0,0);
            EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
            poz++;
        }
    }
}


long int draw_timetrace_INPUT_pulse_MOD(int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                        double p_error,int delay,int pulse_length,int id_pulse,FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    long int  new_id;
    int SHOW_Init_state = 50;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);

    int s_start=PD->Coupled->s[id_pulse-1];
    int s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];

    for (int t=1;t<=T_MAX; t++) {
        if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length))
            PD->Coupled->s[id_pulse-1]=48+(49-s_start);
        if (t>=SHOW_Init_state+delay+pulse_length) {
            if(PD->Coupled->BRN->active_inputs->check_for_element(id_pulse)>0)
                PD->Coupled->s[id_pulse-1]=s_start;
            if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
            if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
        }
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length))
            PD->Coupled->s[id_pulse-1]=48+(49-s_start);
        if (t>=SHOW_Init_state+delay+pulse_length) {
                            PD->Coupled->s[id_pulse-1]=s_start;
            if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
            if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
        }
        draw_timestep_MOD(PD,t,f,y0,W, T_MAX, TW, NameBuff);
    }
    new_id= PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
    return (new_id);
}


long int draw_timetrace_INPUT_pulse_ASSYNC(int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int delay,int pulse_length,int id_pulse,FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    long int  new_id;
    int SHOW_Init_state = 50;
    int Average_over=1000;
    doublepointer *av_state;
    
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
  
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);

    int s_start=PD->Coupled->s[id_pulse-1];
    int s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];
    
    av_state = new doublepointer[T_MAX+1];
    for (int t=1;t<=T_MAX; t++) {
        av_state[t] = new double[PD->Coupled->BRN->N];
        for (unsigned int snow=0; snow<PD->Coupled->BRN->N; snow++) av_state[t][snow]=0;
    }
    
    for (unsigned long int stat=1; stat<=Average_over; stat++) {
        PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
        
        for (int t=1;t<=T_MAX; t++) {
            if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length))
                PD->Coupled->s[id_pulse-1]=48+(49-s_start);
            if (t>=SHOW_Init_state+delay+pulse_length) {
                if(PD->Coupled->BRN->active_inputs->check_for_element(id_pulse)>0)
                    PD->Coupled->s[id_pulse-1]=s_start;
                if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
            }
            
            if(ASSYNC_BIAS == 0 ) PD->Coupled->Asynchronously_update_all_Gates();
            if(ASSYNC_BIAS == 30 ) PD->Coupled->Asynchronously_update_all_Gates(ASSYNC_BIAS,NULL,NULL);
            
            if(ASSYNC_BIAS == 1 ) PD->Coupled->Asynchronously_update_all_Gates_Twoscale(NULL,NULL);
           
            if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length))
                PD->Coupled->s[id_pulse-1]=48+(49-s_start);
            if (t>=SHOW_Init_state+delay+pulse_length) {
                                PD->Coupled->s[id_pulse-1]=s_start;
                if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
            }
            
            for (unsigned int snow=0; snow<PD->Coupled->BRN->N; snow++){
                av_state[t][snow]+=PD->Coupled->s[snow]-48;
            }
        }
    
    }
    for (int t=1;t<=T_MAX; t++) {
        for (unsigned int snow=0; snow<PD->Coupled->BRN->N; snow++) av_state[t][snow]=av_state[t][snow]/(double) Average_over;
      //  printf("t=%d\t %s = %lg\t\t %s = %lg\n",t, PD->Coupled->BRN->MDAT->Node_name[id_pulse],av_state[t][id_pulse-1],PD->Coupled->BRN->MDAT->Node_name[id_pulse+4],av_state[t][id_pulse+3]);
        draw_timestep_ASSYNC_Average(PD,av_state[t],t,f,y0,W, T_MAX, TW, NameBuff);
    }
    new_id= PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
    return (new_id);
}

void draw_timetrace_INPUT_pulse_ERROR(const char errortype[],int howmanyerrors,int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int delay,int pulse_length,int id_pulse,FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    long int  new_id;
    int SHOW_Init_state = 50;
    int Average_over=1000;
    doublepointer *av_state;
    Boolean_Dynamics_TRACKER *Mutant;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
  
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);

    int s_start=PD->Coupled->s[id_pulse-1];
    int s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];
    
    av_state = new doublepointer[T_MAX+1];
    for (int t=1;t<=T_MAX; t++) {
        av_state[t] = new double[PD->Coupled->BRN->N];
        for (unsigned int snow=0; snow<PD->Coupled->BRN->N; snow++) av_state[t][snow]=0;
    }
    
    for (unsigned long int stat=1; stat<=Average_over; stat++) {
        Mutant = PD->Coupled->mutant_dynamics(errortype,howmanyerrors);
        Mutant->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
        
        for (int t=1;t<=T_MAX; t++) {
            if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length))
                Mutant->s[id_pulse-1]=48+(49-s_start);
            if (t>=SHOW_Init_state+delay+pulse_length) {
                if(PD->Coupled->BRN->active_inputs->check_for_element(id_pulse)>0)
                    Mutant->s[id_pulse-1]=s_start;
                if(id_GF>0) Mutant->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) Mutant->s[id_CDL-1]=s_start_CDL;
            }
            
            Mutant->synchronously_update_all_Gates();

            if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length))
                Mutant->s[id_pulse-1]=48+(49-s_start);
            if (t>=SHOW_Init_state+delay+pulse_length) {
                Mutant->s[id_pulse-1]=s_start;
                if(id_GF>0) Mutant->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) Mutant->s[id_CDL-1]=s_start_CDL;
            }
            
            for (unsigned int snow=0; snow<Mutant->BRN->N; snow++){
                av_state[t][snow]+=Mutant->s[snow]-48;
            }
        }

        delete Mutant->BRN;  Mutant->BRN = NULL;
        delete Mutant; Mutant = NULL;
    }
    
    for (int t=1;t<=T_MAX; t++) {
        for (unsigned int snow=0; snow<PD->Coupled->BRN->N; snow++) av_state[t][snow]=av_state[t][snow]/(double) Average_over;
        draw_timestep_ASSYNC_Average(PD,av_state[t],t,f,y0,W, T_MAX, TW, NameBuff);
    }
}



void draw_timetrace_INPUT_pulse_Cycle_MOD(int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                        double p_error,int delay,int pulse_l_ON,int pulse_l_OFF,int id_pulse,FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    int SHOW_Init_state = 50;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);

    int s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];
    int per = (pulse_l_ON + pulse_l_OFF);
    
    for (int t=1;t<=T_MAX; t++) {
       // n_cycles = (unsigned long int) floor(t-SHOW_Init_state /  (pulse_l_ON + pulse_l_OFF));
        if (t>=SHOW_Init_state){
            if ((t-SHOW_Init_state+delay) % per < pulse_l_ON){
              //  printf("ON: t=%d\t t-SHOW_Init_state+delay = %d   per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
             //   getchar();
                PD->Coupled->s[id_pulse-1]=49;
            }
            else {
               // printf("OFF: t=%d\t t-SHOW_Init_state+delay = %d    per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
               // getchar();
                if(PD->Coupled->BRN->active_inputs->check_for_element(id_pulse)>0)
                    PD->Coupled->s[id_pulse-1]=48;
                if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
            }
        }
        
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        if (t>=SHOW_Init_state){
            if ((t-SHOW_Init_state+delay) % per < pulse_l_ON){
               PD->Coupled->s[id_pulse-1]=49;
            }
            else {
                if(PD->Coupled->BRN->active_inputs->check_for_element(id_pulse)>0)
                    PD->Coupled->s[id_pulse-1]=48;
                if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
            }
        }
        draw_timestep_MOD(PD,t,f,y0,W, T_MAX, TW, NameBuff);
    }
}

void draw_timetrace_INPUT_pulse2_Cycles_MOD(int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error,
                                            int delay,int offset_p,
                                            int pulse_l_ON_1,int pulse_l_OFF_1,int id_pulse_1,
                                            int pulse_l_ON_2,int pulse_l_OFF_2,int id_pulse_2,
                                            FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    int SHOW_Init_state = 50;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);

    int s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];
    int per1 = (pulse_l_ON_1 + pulse_l_OFF_1);
    int per2 = (pulse_l_ON_2 + pulse_l_OFF_2);
    
    for (int t=1;t<=T_MAX; t++) {
       // n_cycles = (unsigned long int) floor(t-SHOW_Init_state /  (pulse_l_ON + pulse_l_OFF));
        if (t>=SHOW_Init_state){
            if ((t-SHOW_Init_state+delay) % per1 < pulse_l_ON_1){
              //  printf("ON: t=%d\t t-SHOW_Init_state+delay = %d   per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
             //   getchar();
                PD->Coupled->s[id_pulse_1-1]=49;
            }
            else {
               // printf("OFF: t=%d\t t-SHOW_Init_state+delay = %d    per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
               // getchar();
                if(PD->Coupled->BRN->active_inputs->check_for_element(id_pulse_1)>0)
                    PD->Coupled->s[id_pulse_1-1]=48;
                if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
            }
            if ((t-SHOW_Init_state+delay+offset_p) % per2 < pulse_l_ON_2){
              //  printf("ON: t=%d\t t-SHOW_Init_state+delay = %d   per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
             //   getchar();
                PD->Coupled->s[id_pulse_2-1]=49;
            }
            else {
               // printf("OFF: t=%d\t t-SHOW_Init_state+delay = %d    per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
               // getchar();
                if(PD->Coupled->BRN->active_inputs->check_for_element(id_pulse_2)>0)
                    PD->Coupled->s[id_pulse_2-1]=48;
                if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
            }
        }
        
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        if (t>=SHOW_Init_state){
            if ((t-SHOW_Init_state+delay) % per1 < pulse_l_ON_1){
              //  printf("ON: t=%d\t t-SHOW_Init_state+delay = %d   per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
             //   getchar();
                PD->Coupled->s[id_pulse_1-1]=49;
            }
            else {
               // printf("OFF: t=%d\t t-SHOW_Init_state+delay = %d    per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
               // getchar();
                if(PD->Coupled->BRN->active_inputs->check_for_element(id_pulse_1)>0)
                    PD->Coupled->s[id_pulse_1-1]=48;
                if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
            }
            if ((t-SHOW_Init_state+delay+offset_p) % per2 < pulse_l_ON_2){
              //  printf("ON: t=%d\t t-SHOW_Init_state+delay = %d   per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
             //   getchar();
                PD->Coupled->s[id_pulse_2-1]=49;
            }
            else {
               // printf("OFF: t=%d\t t-SHOW_Init_state+delay = %d    per=%d    t-SHOW_Init_state+delay mod (pulse_l_ON + pulse_l_OFF) = %d\n",t,t-SHOW_Init_state+delay,per,(t-SHOW_Init_state+delay) % per);
               // getchar();
                if(PD->Coupled->BRN->active_inputs->check_for_element(id_pulse_2)>0)
                    PD->Coupled->s[id_pulse_2-1]=48;
                if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
                if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
            }
        }
        draw_timestep_MOD(PD,t,f,y0,W, T_MAX, TW, NameBuff);
    }
}

long int draw_timetrace_INPUT_pulse_MOD(int a_st,int a_st_k,Boolean_Dynamics_TRACKER *D,
                                        double p_error,int delay,int pulse_length,int id_pulse,FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    long int  new_id;
    int SHOW_Init_state = 20;
    
    D->set_state(D->Attractor_valleys[a_st]->Basin_Stored[D->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int id_GF=D->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=D->BRN->MDAT->get_node_ID(CD_low_name);

    int s_start=D->s[id_pulse-1];
    int s_start_GF=D->s[id_GF-1];
    int s_start_CDL=D->s[id_CDL-1];
    
    for (int t=1;t<=T_MAX; t++) {
        if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length))
            D->s[id_pulse-1]=48+(49-s_start);
        if (t>=SHOW_Init_state+delay+pulse_length) {
            if(D->BRN->active_inputs->check_for_element(id_pulse)>0)
                D->s[id_pulse-1]=s_start;
            if(id_GF>0) D->s[id_GF-1]=s_start_GF;
            if(id_CDL>0) D->s[id_CDL-1]=s_start_CDL;
        }
        D->synchronously_update_all_Gates_with_noise(p_error);
        if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length))
            D->s[id_pulse-1]=48+(49-s_start);
        if (t>=SHOW_Init_state+delay+pulse_length) {
            if(D->BRN->active_inputs->check_for_element(id_pulse)>0)
                D->s[id_pulse-1]=s_start;
            if(id_GF>0) D->s[id_GF-1]=s_start_GF;
            if(id_CDL>0) D->s[id_CDL-1]=s_start_CDL;
        }
        draw_timestep_MOD(D,t,f,y0,W, T_MAX, TW, NameBuff);
    }
    new_id= D->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(D->s);
    return (new_id);
}

void draw_timetrace_INPUT_pulse_MOD(int a_st,int a_st_k,Boolean_Dynamics *D,
                                        double p_error,int delay,int pulse_length,int id_pulse,FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    int SHOW_Init_state = 20;
    
    D->set_state(D->get_attractor_state_string(a_st,a_st_k));

    int id_GF=D->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=D->BRN->MDAT->get_node_ID(CD_low_name);
    int s_start=D->s[id_pulse];
    int s_start_GF=D->s[id_GF];
    int s_start_CDL=D->s[id_CDL];
    
    for (int t=1;t<=T_MAX; t++) {
        if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length)) D->s[id_pulse]=(1-s_start);
        if (t>=SHOW_Init_state+delay+pulse_length) {
            if(D->BRN->active_inputs->check_for_element(id_pulse)>0) D->s[id_pulse]=s_start;
            if(id_GF>0) D->s[id_GF]=s_start_GF;
            if(id_CDL>0) D->s[id_CDL]=s_start_CDL;
        }
        D->synchronously_update_all_Gates_with_noise(p_error);
        if ((t>=SHOW_Init_state+delay)&&(t<SHOW_Init_state+delay+pulse_length)) D->s[id_pulse]=(1-s_start);
        if (t>=SHOW_Init_state+delay+pulse_length) {
            if(D->BRN->active_inputs->check_for_element(id_pulse)>0) D->s[id_pulse]=s_start;
            if(id_GF>0) D->s[id_GF]=s_start_GF;
            if(id_CDL>0) D->s[id_CDL]=s_start_CDL;
        }
        draw_timestep_MOD(D,t,f,y0,W, T_MAX, TW, NameBuff);
    }
}

long int draw_timetrace_Delayed_Knockdown_MOD(int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                              double p_error,int delay,int id_KO,FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    long int  new_id;
    int SHOW_Init_state = 20;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    for (int t=1;t<=T_MAX; t++) {
        if (t>=SHOW_Init_state+delay) PD->Coupled->s[id_KO-1]=48;
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        if (t>=SHOW_Init_state+delay) PD->Coupled->s[id_KO-1]=48;
        draw_timestep_MOD(PD,t,f,y0,W, T_MAX, TW, NameBuff);
    }
    new_id= PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
    return (new_id);
}

long int draw_timetrace_Delayed_OverExpression_MOD(int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                   double p_error,int delay,int id_KO,FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    long int  new_id;
    int SHOW_Init_state = 20;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    for (int t=1;t<=T_MAX; t++) {
        if (t>=SHOW_Init_state+delay) PD->Coupled->s[id_KO-1]=49;
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        if (t>=SHOW_Init_state+delay) PD->Coupled->s[id_KO-1]=49;
        draw_timestep_MOD(PD,t,f,y0,W, T_MAX, TW, NameBuff);
    }
    new_id= PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
    return (new_id);
}

void draw_timetrace_Delayed_MultiHIT_MOD(longint_ARRAY_pair *Hits,
                                         int a_st,int a_st_k,
                                         MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                         double p_error,int HIT_delay,
                                         int pulse_length,int pulseID,int pulse_delay,
                                         FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    int SHOW_Init_state = 20;
    int id_GF=0,s_start1=0,s_start_GF=0;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    if(pulseID>0){
        id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
       // int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
        s_start1=PD->Coupled->s[pulseID-1];
        s_start_GF=PD->Coupled->s[id_GF-1];
    }
    
    for (int t=1;t<=T_MAX; t++) {
        
        if(pulseID>0){
            if (t>=SHOW_Init_state + pulse_delay) {
                PD->Coupled->s[pulseID-1]=48+(49-s_start1);
            }
            if (t>=SHOW_Init_state + pulse_delay + pulse_length) {
                PD->Coupled->s[pulseID-1]=s_start1;
            }
        }
        
        if (t>=SHOW_Init_state + pulse_delay + HIT_delay) {
            for(unsigned long int l=1;l<=Hits->N;l++){
                PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
            }
        }
        
        if (t>=SHOW_Init_state + pulse_delay + HIT_delay)
            draw_timestep_MOD(Hits,PD,t,f,y0,W, T_MAX, TW, NameBuff);
        else draw_timestep_MOD(PD,t,f,y0,W, T_MAX, TW, NameBuff);
        
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        
        for(unsigned long int l=1;l<=Hits->N;l++)
            if (t>=SHOW_Init_state + HIT_delay) {
                for(unsigned long int l=1;l<=Hits->N;l++){
                    PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
                }
            }
    }
}

void draw_timetrace_Stochastic_MultiHIT_MOD(unsigned long int id_input,double p_input,longint_ARRAY_pair *Hits,
                                            int a_st,int a_st_k,
                                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                            double p_error,
                                            FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    int id_GF=0,s_start1=0,s_start_GF=0;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    s_start1=PD->Coupled->s[id_input-1];
    s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL];

    for (int t=1;t<=T_MAX; t++) {
        for(unsigned long int l=1;l<=Hits->N;l++){
            PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
        }
        double r=rand_real(0,1);
        if(r<p_input) PD->Coupled->s[id_input-1]=49;
        else PD->Coupled->s[id_input-1]=48;
        
        draw_timestep_MOD(Hits,PD,t,f,y0,W, T_MAX, TW, NameBuff);
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        PD->Coupled->s[id_GF-1]=s_start_GF;
        PD->Coupled->s[id_CDL-1]=s_start_CDL;
    }
}

void draw_timetrace_Stochastic_HITS_MultiHIT_MOD(unsigned long int id_input,
                                                 double p_hits,longint_ARRAY_pair *Hits,
                                                 int a_st,int a_st_k,
                                                 MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                 double p_error,
                                                 FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    int id_GF=0,s_start1=0,s_start_GF=0,s_start_CDL=0;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    s_start1   =PD->Coupled->s[id_input-1];
    s_start_GF =PD->Coupled->s[id_GF-1];
    s_start_CDL=PD->Coupled->s[id_CDL-1];
    
    for (int t=1;t<=T_MAX; t++) {
        PD->Coupled->s[id_input-1]=49;
        double r=rand_real(0,1);
        if(r<p_hits){
            for(unsigned long int l=1;l<=Hits->N;l++){
                PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
                // else PD->Coupled->s[Hits->A[l]-1]=48+(1-Hits->B[l]);
            }
            draw_timestep_MOD(Hits,PD,t,f,y0,W, T_MAX, TW, NameBuff);
        }
        else draw_timestep_MOD(PD,t,f,y0,W, T_MAX, TW, NameBuff);
        
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        PD->Coupled->s[id_GF-1]=s_start_GF;
        PD->Coupled->s[id_CDL-1]=s_start_CDL;
    }
}

void draw_timetrace_Stochastic_MultiHIT_MOD(unsigned long int id_input,double p_input,
                                            unsigned long int id_input2,double p_input2,
                                            longint_ARRAY_pair *Hits,
                                            int a_st,int a_st_k,
                                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                            double p_error,
                                            FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    int id_GF=0,s_start1=0,s_start2=0,s_start_GF=0;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    s_start1=PD->Coupled->s[id_input-1];
    s_start2=PD->Coupled->s[id_input2-1];
    s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];
    
    for (int t=1;t<=T_MAX; t++) {
        for(unsigned long int l=1;l<=Hits->N;l++){
            PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
        }
        double r=rand_real(0,1);
        if(r<p_input) PD->Coupled->s[id_input-1]=49;
        else PD->Coupled->s[id_input-1]=48;
        
        r=rand_real(0,1);
        if(r<p_input2) PD->Coupled->s[id_input2-1]=49;
        else PD->Coupled->s[id_input2-1]=48;
        
        draw_timestep_MOD(Hits,PD,t,f,y0,W, T_MAX, TW, NameBuff);
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        PD->Coupled->s[id_GF-1]=s_start_GF;
        PD->Coupled->s[id_CDL-1]=s_start_CDL;
    }
}


long int draw_timetrace_INPUT_twopulse_MOD(int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                           double p_error,int delay,int pulse_length,int id_pulse1,int id_pulse2,
                                           FILE *f,double y0,int W){
    int T_MAX=TIMESTEPS_DEFAULT,TW=3,NameBuff=55;
    long int  new_id;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    int s_start1=PD->Coupled->s[id_pulse1-1];
    int s_start2=PD->Coupled->s[id_pulse2-1];
    int s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];
    for (int t=1;t<=T_MAX; t++) {
        if (t>=30+delay) {
            PD->Coupled->s[id_pulse1-1]=48+(49-s_start1);
            PD->Coupled->s[id_pulse2-1]=48+(49-s_start2);
        }
        if (t>=30+delay+pulse_length) {
            PD->Coupled->s[id_pulse1-1]=s_start1;
            PD->Coupled->s[id_pulse2-1]=s_start2;
            PD->Coupled->s[id_GF-1]=s_start_GF;
            PD->Coupled->s[id_CDL-1]=s_start_CDL;
        }
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        if (t>=30+delay) {
            PD->Coupled->s[id_pulse1-1]=48+(49-s_start1);
            PD->Coupled->s[id_pulse2-1]=48+(49-s_start2);
        }
        if (t>=30+delay+pulse_length) {
            PD->Coupled->s[id_pulse1-1]=s_start1;
            PD->Coupled->s[id_pulse2-1]=s_start2;
            PD->Coupled->s[id_GF-1]=s_start_GF;
            PD->Coupled->s[id_CDL-1]=s_start_CDL;
        }
        draw_timestep_MOD(PD,t,f,y0,W, T_MAX, TW, NameBuff);
    }
    new_id= PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
    return (new_id);
}

void run_OneNode_pulse_experiment_MOD(const char nodename[],int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                      double p_error,int pulse_l){
    int W=8,YTOP,XTOP;
    long int id_new;
    FILE *f;
    char fn[400];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    int id_Input=PD->BRN->MDAT->get_node_ID(nodename);
    if (id_Input == 0) {
        printf("Node altered during pulse experiment %s does not exist in network %s\n",nodename,PD->Coupled->BRN->name);
        getchar();
        return;
    }
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s_pulse_%d\n",MODEL_Directory_system,PD->BRN->path,nodename,pulse_l);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d.eps",MODEL_Directory,PD->BRN->path,nodename,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        id_new=draw_timetrace_INPUT_pulse_MOD(a_st,a_st_k,PD,p_error,delay,pulse_l,id_Input,f,YTOP,W);
        printf("When flipping %s for %d steps - delay %d, attractor %d went to %ld\n",
               nodename,pulse_l,delay,a_st,id_new);
        //getchar();
        sprintf(fn,"%s/%s/_EXP/%s_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d",MODEL_Directory_system,PD->BRN->path,nodename,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
        close_EPSDRAW_file(f,fn,2);
    }
    
}

void run_OneNode_pulse_experiment_ASSYNC(const char nodename[],int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                      double p_error,int pulse_l){
    int W=8,YTOP,XTOP;
    long int id_new;
    FILE *f;
    char fn[400];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    int id_Input=PD->BRN->MDAT->get_node_ID(nodename);
    if (id_Input == 0) {
        printf("Node altered during pulse experiment %s does not exist in network %s\n",nodename,PD->Coupled->BRN->name);
        getchar();
        return;
    }
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s_ASSYNC-Bias%d_pulse_%d\n",MODEL_Directory_system,PD->BRN->path,nodename,ASSYNC_BIAS,pulse_l);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s_ASSYNC-Bias%d_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d.eps",MODEL_Directory,PD->BRN->path,nodename,ASSYNC_BIAS,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        id_new=draw_timetrace_INPUT_pulse_ASSYNC(a_st,a_st_k,PD,delay,pulse_l,id_Input,f,YTOP,W);
        printf("When flipping %s for %d steps - delay %d, attractor %d went to %ld\n",
               nodename,pulse_l,delay,a_st,id_new);
        //getchar();
        sprintf(fn,"%s/%s/_EXP/%s_ASSYNC-Bias%d_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d",MODEL_Directory,PD->BRN->path,nodename,ASSYNC_BIAS,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
        close_EPSDRAW_file(f,fn,2);
    }
}

void run_OneNode_pulse_experiment_ModelErrors(const char nodename[],int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int pulse_l,int error_strength){
        int W=8,YTOP,XTOP;
        FILE *f;
        char fn[400];
        
        PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
        // PD->Coupled->print_state();//getchar();
        
        int id_Input=PD->BRN->MDAT->get_node_ID(nodename);
        if (id_Input == 0) {
            printf("Node altered during pulse experiment %s does not exist in network %s\n",nodename,PD->Coupled->BRN->name);
            getchar();
            return;
        }
        long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
        int N_lines=PD->BRN->N+10;
        YTOP=(W*N_lines+10*W)+5;
        XTOP=2*TIMESTEPS_DEFAULT+5;
        
        sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
        system(fn);
    
        sprintf(fn,"mkdir %s/%s/_EXP/%s_LINK-%d-ERRORS_pulse_%d\n",MODEL_Directory_system,PD->BRN->path,nodename,error_strength,pulse_l);
        system(fn);
        for (int delay=1; delay<=DMAX; delay++) {
            sprintf(fn,"%s/%s/_EXP/%s_LINK-%d-ERRORS_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d.eps",MODEL_Directory,PD->BRN->path,nodename,error_strength,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
            printf("%s\n",fn);
            f=fopen(fn,"w");
            EpsInit(f,-5,-5,XTOP,YTOP);
            EpsSetFont(f,"Times-Roman",W);
            draw_timetrace_INPUT_pulse_ERROR("remove_links",error_strength,a_st,a_st_k,PD,delay,pulse_l,id_Input,f,YTOP,W);
            sprintf(fn,"%s/%s/_EXP/%s_LINK-%d-ERRORS_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d",MODEL_Directory,PD->BRN->path,nodename,error_strength,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
            close_EPSDRAW_file(f,fn,2);
        }
      
        sprintf(fn,"mkdir %s/%s/_EXP/%s_GATE-%d-ERRORS_pulse_%d\n",MODEL_Directory_system,PD->BRN->path,nodename,error_strength,pulse_l);
        system(fn);
        for (int delay=1; delay<=DMAX; delay++) {
            sprintf(fn,"%s/%s/_EXP/%s_GATE-%d-ERRORS_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d.eps",MODEL_Directory,PD->BRN->path,nodename,error_strength,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
            printf("%s\n",fn);
            f=fopen(fn,"w");
            EpsInit(f,-5,-5,XTOP,YTOP);
            EpsSetFont(f,"Times-Roman",W);
            draw_timetrace_INPUT_pulse_ERROR("gate_oneoutput_flips",error_strength, a_st,a_st_k,PD,delay,pulse_l,id_Input,f,YTOP,W);
            sprintf(fn,"%s/%s/_EXP/%s_GATE-%d-ERRORS_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d",MODEL_Directory,PD->BRN->path,nodename,error_strength,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
            close_EPSDRAW_file(f,fn,2);
        }
       
        sprintf(fn,"mkdir %s/%s/_EXP/%s_NODELOCK-%d-ERRORS_pulse_%d\n",MODEL_Directory_system,PD->BRN->path,nodename,error_strength,pulse_l);
        system(fn);
        for (int delay=1; delay<=DMAX; delay++) {
            sprintf(fn,"%s/%s/_EXP/%s_NODELOCK-%d-ERRORS_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d.eps",MODEL_Directory,PD->BRN->path,nodename,error_strength,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
            printf("%s\n",fn);
            f=fopen(fn,"w");
            EpsInit(f,-5,-5,XTOP,YTOP);
            EpsSetFont(f,"Times-Roman",W);
            draw_timetrace_INPUT_pulse_ERROR("lock_nodes_ONOFF",error_strength, a_st,a_st_k,PD,delay,pulse_l,id_Input,f,YTOP,W);
            sprintf(fn,"%s/%s/_EXP/%s_NODELOCK-%d-ERRORS_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d",MODEL_Directory,PD->BRN->path,nodename,error_strength,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
            close_EPSDRAW_file(f,fn,2);
        }
}

void run_OneNode_pulse_Cycle_experiment_MOD(const char nodename[],int a_st,int a_st_k,
                                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                            double p_error,int pulse_l_ON,int pulse_l_OFF){
    int W=8,YTOP,XTOP;
    //long int id_new;
    FILE *f;
    char fn[400];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    int id_Input=PD->BRN->MDAT->get_node_ID(nodename);
    if (id_Input == 0) {
        printf("Node altered during pulse experiment %s does not exist in network %s\n",nodename,PD->Coupled->BRN->name);
        getchar();
        return;
    }
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s_pulse_Cycle_%d-ON_%d-OFF\n",MODEL_Directory_system,PD->BRN->path,nodename,pulse_l_ON,pulse_l_OFF);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s_pulse_Cycle_%d-ON_%d-OFF/%s_pulse_experiment_A%d_As-%d__pulse_Cycle_%d-ON_%d-OFF_delay-%d.eps",
                MODEL_Directory,PD->BRN->path,nodename,pulse_l_ON,pulse_l_OFF,nodename,a_st,a_st_k,pulse_l_ON,pulse_l_OFF,delay);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
     //   draw_timetrace_INPUT_pulse_MOD(a_st,a_st_k,PD,p_error,delay,pulse_l,id_Input,f,YTOP,W);
        draw_timetrace_INPUT_pulse_Cycle_MOD(a_st,a_st_k,PD,p_error,delay,pulse_l_ON,pulse_l_OFF,id_Input,f,YTOP,W);
        //getchar();
        sprintf(fn,"%s/%s/_EXP/%s_pulse_Cycle_%d-ON_%d-OFF/%s_pulse_experiment_A%d_As-%d__pulse_Cycle_%d-ON_%d-OFF_delay-%d",
                MODEL_Directory,PD->BRN->path,nodename,pulse_l_ON,pulse_l_OFF,nodename,a_st,a_st_k,pulse_l_ON,pulse_l_OFF,delay);
        close_EPSDRAW_file(f,fn,2);
    }
    
}

void run_Pulse2_Cycles_experiment_MOD(const char nodename_1[],const char nodename_2[],int a_st,int a_st_k,
                                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                            double p_error,int pulse_l_ON_1,int pulse_l_OFF_1,int pulse_l_ON_2,int pulse_l_OFF_2){
    int W=8,YTOP,XTOP;
    //long int id_new;
    FILE *f;
    char fn[400];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    int id_Input_1=PD->BRN->MDAT->get_node_ID(nodename_1);
    if (id_Input_1 == 0) {
        printf("Node altered during pulse experiment %s does not exist in network %s\n",nodename_1,PD->Coupled->BRN->name);
        getchar();
        return;
    }
    int id_Input_2=PD->BRN->MDAT->get_node_ID(nodename_2);
    if (id_Input_2 == 0) {
        printf("Node altered during pulse experiment %s does not exist in network %s\n",nodename_2,PD->Coupled->BRN->name);
        getchar();
        return;
    }
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/Pulse2_Cycles_%s_%d-ON_%d-OFF__%s_%d-ON_%d-OFF\n",MODEL_Directory_system,PD->BRN->path,nodename_1,pulse_l_ON_1,pulse_l_OFF_1,nodename_2,pulse_l_ON_2,pulse_l_OFF_2);
    system(fn);

    for (int delay=1; delay<=DMAX; delay++) {
        for (int offset_p=0; offset_p<pulse_l_ON_1+pulse_l_OFF_1+pulse_l_ON_2+pulse_l_OFF_2; offset_p++) {
            sprintf(fn,"%s/%s/_EXP/Pulse2_Cycles_%s_%d-ON_%d-OFF__%s_%d-ON_%d-OFF/Pulse2_experiment_A%d_As-%d__%s_%d-ON_%d-OFF__%s_%d-ON_%d-OFF_offset-%d__delay-%d.eps",
                MODEL_Directory,PD->BRN->path,nodename_1,pulse_l_ON_1,pulse_l_OFF_1,nodename_2,pulse_l_ON_2,pulse_l_OFF_2,a_st,a_st_k,nodename_1,pulse_l_ON_1,pulse_l_OFF_1,nodename_2,pulse_l_ON_2,pulse_l_OFF_2,offset_p,delay);
            printf("%s\n",fn);
            f=fopen(fn,"w");
            
            EpsInit(f,-5,-5,XTOP,YTOP);
            EpsSetFont(f,"Times-Roman",W);
         //   draw_timetrace_INPUT_pulse_MOD(a_st,a_st_k,PD,p_error,delay,pulse_l,id_Input,f,YTOP,W);
            draw_timetrace_INPUT_pulse2_Cycles_MOD(a_st,a_st_k,PD,p_error,delay,offset_p,pulse_l_ON_1,pulse_l_OFF_1,id_Input_1,pulse_l_ON_2,pulse_l_OFF_2,id_Input_2,f,YTOP,W);
            //getchar();
            sprintf(fn,"%s/%s/_EXP/Pulse2_Cycles_%s_%d-ON_%d-OFF__%s_%d-ON_%d-OFF/Pulse2_experiment_A%d_As-%d__%s_%d-ON_%d-OFF__%s_%d-ON_%d-OFF_offset-%d__delay-%d",
                MODEL_Directory,PD->BRN->path,nodename_1,pulse_l_ON_1,pulse_l_OFF_1,nodename_2,pulse_l_ON_2,pulse_l_OFF_2,a_st,a_st_k,nodename_1,pulse_l_ON_1,pulse_l_OFF_1,nodename_2,pulse_l_ON_2,pulse_l_OFF_2,offset_p,delay);
           close_EPSDRAW_file(f,fn,2);
        }
    }
    
}
void run_OneNode_Delayed_KO_experiment_MOD(const char nodename[],int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    int W=8,YTOP,XTOP;
    long int id_new;
    FILE *f;
    char fn[400];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    int id_KO=PD->BRN->MDAT->get_node_ID(nodename);
    
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s_KnockDown\n",MODEL_Directory_system,PD->BRN->path,nodename);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s_KnockDown/%s_KnockDown_experiment_A%d_As-%d__delay-%d.eps",MODEL_Directory,PD->BRN->path,nodename,nodename,a_st,a_st_k,delay);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        id_new=draw_timetrace_Delayed_Knockdown_MOD(a_st,a_st_k,PD,p_error,delay,id_KO,f,YTOP,W);
        printf("When knocking down %s at delay %d, attractor %d went to %ld\n",
               nodename,delay,a_st,id_new);
        //getchar();
        sprintf(fn,"%s/%s/_EXP/%s_KnockDown/%s_KnockDown_experiment_A%d_As-%d__delay-%d",MODEL_Directory_system,PD->BRN->path,nodename,nodename,a_st,a_st_k,delay);
        close_EPSDRAW_file(f,fn,2);
    }
}

void run_OneNode_Delayed_OE_experiment_MOD(const char nodename[],int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    int W=8,YTOP,XTOP;
    long int id_new;
    FILE *f;
    char fn[400];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    int id_KO=PD->BRN->MDAT->get_node_ID(nodename);
    
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s_OverExpression\n",MODEL_Directory_system,PD->BRN->path,nodename);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s_OverExpression/%s_OverExpression_experiment_A%d_As-%d__delay-%d.eps",MODEL_Directory,PD->BRN->path,nodename,nodename,a_st,a_st_k,delay);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        id_new=draw_timetrace_Delayed_OverExpression_MOD(a_st,a_st_k,PD,p_error,delay,id_KO,f,YTOP,W);
        printf("When knocking down %s at delay %d, attractor %d went to %ld\n",
               nodename,delay,a_st,id_new);
        //getchar();
        sprintf(fn,"%s/%s/_EXP/%s_OverExpression/%s_OverExpression_experiment_A%d_As-%d__delay-%d",MODEL_Directory_system,PD->BRN->path,nodename,nodename,a_st,a_st_k,delay);
        close_EPSDRAW_file(f,fn,2);
    }
}

void run_Delayed_Hits_experiment_MOD(longint_ARRAY_pair *Hits,int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    int W=8,YTOP,XTOP;
    FILE *f;
    char fn[400],fn2[300];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn2, "Hits");
    for (unsigned long int l=1; l<=Hits->N; l++) {
        sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
    }
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD->BRN->path,fn2);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s/%s_experiment_A%d_As-%d__delay-%d.eps",MODEL_Directory,PD->BRN->path,fn2,fn2,a_st,a_st_k,delay);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        draw_timetrace_Delayed_MultiHIT_MOD(Hits,a_st,a_st_k,PD,p_error,delay,0,0,0,f,YTOP,W);
        //getchar();
        sprintf(fn,"%s/%s/_EXP/%s/%s_experiment_A%d_As-%d__delay-%d",MODEL_Directory_system,PD->BRN->path,fn2,fn2,a_st,a_st_k,delay);
        close_EPSDRAW_file(f,fn,2);
    }
}

void run_stochastic_Input_experiment_MOD(unsigned long int id_input,double p_input,longint_ARRAY_pair *Hits,int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    int W=8,YTOP,XTOP;
    FILE *f;
    char fn[400],fn2[300];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    //  long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn2, "Hits");
    if(Hits->N>0){
        for (unsigned long int l=1; l<=Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
        }
    }
    else sprintf(fn2,"%s_WILDTYPE",fn2);
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD->BRN->path,fn2);
    system(fn);
    
    sprintf(fn,"%s/%s/_EXP/%s/%s_Stochastic_%s-%.2lf__A%d_As-%d.eps",MODEL_Directory,PD->BRN->path,fn2,fn2,PD->BRN->MDAT->Node_name[id_input],p_input,a_st,a_st_k);
    printf("%s\n",fn);
    f=fopen(fn,"w");
    
    EpsInit(f,-5,-5,XTOP,YTOP);
    EpsSetFont(f,"Times-Roman",W);
    draw_timetrace_Stochastic_MultiHIT_MOD(id_input,p_input,Hits,a_st,a_st_k,PD,p_error,f,YTOP,W);
    //getchar();
    sprintf(fn,"%s/%s/_EXP/%s/%s_Stochastic_%s-%.2lf__A%d_As-%d",MODEL_Directory_system,PD->BRN->path,fn2,fn2,PD->BRN->MDAT->Node_name[id_input],p_input,a_st,a_st_k);
    close_EPSDRAW_file(f,fn,2);
}

void run_stochastic_HITS_experiment_MOD(unsigned long int id_input,double p_hits,longint_ARRAY_pair *Hits,int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    int W=8,YTOP,XTOP;
    FILE *f;
    char fn[400],fn2[300];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    //  long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn2, "Stochastic_Hits_%.2lg_",p_hits);
    if(Hits->N>0){
        for (unsigned long int l=1; l<=Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
        }
    }
    else sprintf(fn2,"%s_WILDTYPE",fn2);
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD->BRN->path,fn2);
    system(fn);
    
    sprintf(fn,"%s/%s/_EXP/%s/%s_%s__A%d_As-%d.eps",MODEL_Directory,PD->BRN->path,fn2,fn2,PD->BRN->MDAT->Node_name[id_input],a_st,a_st_k);
    printf("%s\n",fn);
    f=fopen(fn,"w");
    
    EpsInit(f,-5,-5,XTOP,YTOP);
    EpsSetFont(f,"Times-Roman",W);
    draw_timetrace_Stochastic_HITS_MultiHIT_MOD(id_input,p_hits,Hits,a_st,a_st_k,PD,p_error,f,YTOP,W);
    //getchar();
    sprintf(fn,"%s/%s/_EXP/%s/%s_%s__A%d_As-%d",MODEL_Directory_system,PD->BRN->path,fn2,fn2,PD->BRN->MDAT->Node_name[id_input],a_st,a_st_k);
    close_EPSDRAW_file(f,fn,2);
}

void run_stochastic_Input_experiment_MOD(unsigned long int id_input,double p_input,
                                         unsigned long int id_input2,double p_input2,
                                         longint_ARRAY_pair *Hits,int a_st,int a_st_k,
                                         MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    int W=8,YTOP,XTOP;
    FILE *f;
    char fn[400],fn2[300];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    //  long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn2, "Hits");
    if(Hits->N>0){
        for (unsigned long int l=1; l<=Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
        }
    }
    else sprintf(fn2,"%s_WILDTYPE",fn2);
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD->BRN->path,fn2);
    system(fn);
    
    sprintf(fn,"%s/%s/_EXP/%s/%s_Stochastic_%s-%.2lf_%s-%.2lf__A%d_As-%d.eps",MODEL_Directory,PD->BRN->path,fn2,fn2,PD->BRN->MDAT->Node_name[id_input],p_input,PD->BRN->MDAT->Node_name[id_input2],p_input2,a_st,a_st_k);
    printf("%s\n",fn);
    f=fopen(fn,"w");
    
    EpsInit(f,-5,-5,XTOP,YTOP);
    EpsSetFont(f,"Times-Roman",W);
    draw_timetrace_Stochastic_MultiHIT_MOD(id_input,p_input,id_input2,p_input2,Hits,a_st,a_st_k,PD,p_error,f,YTOP,W);
    //getchar();
    sprintf(fn,"%s/%s/_EXP/%s/%s_Stochastic_%s-%.2lf_%s-%.2lf__A%d_As-%d",MODEL_Directory_system,PD->BRN->path,fn2,fn2,PD->BRN->MDAT->Node_name[id_input],p_input,PD->BRN->MDAT->Node_name[id_input2],p_input2,a_st,a_st_k);
    close_EPSDRAW_file(f,fn,2);
}

void run_Delayed_Hits_Pulse_experiment_MOD(longint_ARRAY_pair *Hits,int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error,const char pulsename[], int lpulse,int shift_pulse){
    int W=8,YTOP,XTOP;
    FILE *f;
    char fn[400],fn2[300];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    int id_pulse=PD->Coupled->BRN->MDAT->get_node_ID(pulsename);
    
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn2, "Hits");
    for (unsigned long int l=1; l<=Hits->N; l++) {
        sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[Hits->A[l]],Hits->B[l]);
    }
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD->BRN->path,fn2);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s/%s_experiment_A%d_As-%d__%s_for_%d__delay-%d_Hitat-%d.eps",MODEL_Directory,PD->BRN->path,fn2,fn2,a_st,a_st_k,PD->BRN->MDAT->Node_name[id_pulse],lpulse,delay,shift_pulse);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        draw_timetrace_Delayed_MultiHIT_MOD(Hits,a_st,a_st_k,PD,p_error,delay+shift_pulse,
                                            lpulse,id_pulse,delay,f,YTOP,W);
        //getchar();
        sprintf(fn,"%s/%s/_EXP/%s/%s_experiment_A%d_As-%d__%s_for_%d__delay-%d_Hitat-%d",MODEL_Directory_system,PD->BRN->path,fn2,fn2,a_st,a_st_k,PD->BRN->MDAT->Node_name[id_pulse],lpulse,delay,shift_pulse);
        close_EPSDRAW_file(f,fn,2);
    }
}



void run_OneNode_pulse_experiment_MOD(const char nodename[],int a_st,int a_st_k,Boolean_Dynamics_TRACKER *D,
                                      double p_error,int pulse_l){
    int W=8,YTOP,XTOP;
    long int id_new;
    FILE *f;
    char fn[400];
    
    D->set_state(D->Attractor_valleys[a_st]->Basin_Stored[D->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    int id_Input=D->BRN->MDAT->get_node_ID(nodename);
    if (id_Input == 0) {
        printf("Node altered during pulse experiment %s does not exist in network %s\n",nodename,D->BRN->name);
        getchar();
        return;
    }
    long int DMAX=D->Attractor_valleys[a_st]->N_a+1;
    int N_lines=D->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,D->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s_pulse_%d\n",MODEL_Directory_system,D->BRN->path,nodename,pulse_l);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d.eps",MODEL_Directory,D->BRN->path,nodename,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        id_new=draw_timetrace_INPUT_pulse_MOD(a_st,a_st_k,D,p_error,delay,pulse_l,id_Input,f,YTOP,W);
        printf("When flipping %s for %d steps - delay %d, attractor %d went to %ld\n",
               nodename,pulse_l,delay,a_st,id_new);
        //getchar();
        sprintf(fn,"%s/%s/_EXP/%s_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d",MODEL_Directory_system,D->BRN->path,nodename,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
        close_EPSDRAW_file(f,fn,2);
    }
    
}

void run_OneNode_pulse_experiment_MOD(const char nodename[],int a_st,int a_st_k,Boolean_Dynamics *D,
                                      double p_error,int pulse_l){
    int W=8,YTOP,XTOP;
    FILE *f;
    char fn[400];
    
    D->set_state(D->get_attractor_state_string(a_st,a_st_k));

    unsigned int id_Input=D->BRN->MDAT->get_node_ID(nodename);
    if (id_Input == 0) {
        printf("Node altered during pulse experiment %s does not exist in network %s\n",nodename,D->BRN->name);
        getchar();
        return;
    }
    long int DMAX=D->Attractor_Cycle[a_st]->N+1;
    
    int N_lines=D->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,D->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s_pulse_%d\n",MODEL_Directory_system,D->BRN->path,nodename,pulse_l);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d.eps",MODEL_Directory,D->BRN->path,nodename,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        draw_timetrace_INPUT_pulse_MOD(a_st,a_st_k,D,p_error,delay,pulse_l,id_Input,f,YTOP,W);
        sprintf(fn,"%s/%s/_EXP/%s_pulse_%d/%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d",MODEL_Directory_system,D->BRN->path,nodename,pulse_l,nodename,a_st,a_st_k,pulse_l,delay);
        close_EPSDRAW_file(f,fn,2);
    }
    
}

int chech_time_of_death_following_OneNode_pulse_experiment_MOD(const char nodename[],int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    
    int T_MAX=100;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int id_pulse=PD->BRN->MDAT->get_node_ID(nodename);
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    int id_dead=PD->Coupled->BRN->MDAT->get_node_ID("CAD");
    
    int s_start=PD->Coupled->s[id_pulse-1];
    int s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];
    int p_min=-1;
    for (int pulse_l=1;pulse_l<=T_MAX; pulse_l++){
        PD->Coupled->s[id_pulse-1]=48+(49-s_start);
        for (int t=1;t<=T_MAX; t++) {
            PD->Coupled->synchronously_update_all_Gates();
            if (t>=pulse_l) {
                PD->Coupled->s[id_pulse-1]=s_start;
                if(id_pulse != id_GF) PD->Coupled->s[id_GF-1]=s_start_GF;
                if(id_pulse != id_CDL) PD->Coupled->s[id_CDL-1]=s_start_CDL;
            }
            if(p_min<0) {
                if (PD->Coupled->s[id_dead-1]==49) {
                    p_min=pulse_l;
                    t=T_MAX+10;
                    pulse_l=T_MAX+10;
                }
            }
        }
    }
    return (p_min);
}


int count_cell_cycles_following_OneNode_pulse_experiment_MOD(const char nodename[],int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int pulse_l){
    int nc=0;
    int T_MAX=TIMESTEPS_DEFAULT;
    
    if(PD->Coupled->Attractor_valleys[a_st]->N_a>1) return (-1);
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    //   PD->Coupled->print_state();//getchar();
    
    int id_pulse=PD->BRN->MDAT->get_node_ID(nodename);
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    int s_start=PD->Coupled->s[id_pulse-1];
    int s_start_GF=PD->Coupled->s[id_GF-1];
    int s_start_CDL=PD->Coupled->s[id_CDL-1];
    int idm=PD->Coupled->BRN->MDAT->get_node_ID("f4N_DNA");
    if (idm==0) idm=PD->BRN->MDAT->get_node_ID("4N_DNA");
    
    int r0=0;
    PD->Coupled->s[id_pulse-1]=48+(49-s_start);
    for (int t=1;t<=T_MAX; t++) {
        PD->Coupled->synchronously_update_all_Gates();
        if (t>=pulse_l) {
            PD->Coupled->s[id_pulse-1]=s_start;
            PD->Coupled->s[id_GF-1]=s_start_GF;
            PD->Coupled->s[id_CDL-1]=s_start_CDL;
        }
        if ((r0==0)&&(PD->Coupled->s[idm-1]==49)) {r0=1;}
        if ((r0==1)&&(PD->Coupled->s[idm-1]==48)) {r0=0;nc++; printf("t=%d\t ncycles =%d\n",t,nc);}
    }
    return(nc);
}

Cell_cycle_stats_in_run count_cell_cycles_stochastic_HITS(const char inputnode[],double p_hits,
                                                          longint_ARRAY_pair *Hits,
                                                          int a_st,int a_st_k,
                                                          MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                          double p_error,unsigned long long int T_MAX){
    Cell_cycle_stats_in_run Statsnow;
    
    
    Statsnow.N_cycles=0;
    Statsnow.N_endoredupl_T=0;
    Statsnow.N_endoredupl_G2=0;
    Statsnow.N_endoredupl_Aneup=0;
    Statsnow.MAX_INTERVAL=0;
    Statsnow.T_live=T_MAX;
    Statsnow.N_deaths=0;
    Statsnow.Time_G1=0;
    Statsnow.Time_R=0;
    Statsnow.Time_G2=0;
    Statsnow.Time_M=0;
    Statsnow.Time_Telo=0;
    Statsnow.N_G1=0;
    Statsnow.N_R=0;
    Statsnow.N_G2=0;
    Statsnow.N_M=0;
    Statsnow.N_T=0;
    for (int t=0;t<=MAX_ARREST; t++){
        Statsnow.Dist_G1[t]=0;
        Statsnow.Dist_G2[t]=0;
        Statsnow.Dist_Mitosis[t]=0;
        Statsnow.Dist_Replication[t]=0;
        Statsnow.Dist_Telophase[t]=0;
    }
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    //   PD->Coupled->print_state();//getchar();
    
    //  int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    // int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    //   int s_start_GF=PD->Coupled->s[id_GF-1];
    // int s_start_CDL=PD->Coupled->s[id_CDL-1];
    //int s_start=PD->Coupled->s[id_input-1];
    
    int id_input=PD->Coupled->BRN->MDAT->get_node_ID(inputnode);
    int id_r=PD->Coupled->BRN->MDAT->get_node_ID("Replication");
    int id_4N=PD->Coupled->BRN->MDAT->get_node_ID("4N_DNA");
    if (id_4N == 0) id_4N=PD->Coupled->BRN->MDAT->get_node_ID("f4N_DNA");
    int id_M=PD->Coupled->BRN->MDAT->get_node_ID("U_Kinetochores");
    int id_SAC=PD->Coupled->BRN->MDAT->get_node_ID("A_Kinetochores");
    int id_cyto=PD->Coupled->BRN->MDAT->get_node_ID("Ect2");
    int id_Dead=PD->Coupled->BRN->MDAT->get_node_ID("CAD");
    
    int r0=0;int t_G1=1;int t_G2=0;int t_M=0;int t_R=0;int t_T=0;
    
    for (unsigned long long int t=1;t<=T_MAX; t++) {
        for(unsigned long int l=1;l<=Hits->N;l++){
            double r=rand_real(0,1);
            if(r<p_hits) PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
            // else PD->Coupled->s[Hits->A[l]-1]=48+(1-Hits->B[l]);
        }
        
        PD->Coupled->s[id_input-1]=49;
        //        PD->Coupled->s[id_GF-1]=s_start_GF;
        //        PD->Coupled->s[id_CDL-1]=s_start_CDL;
        
        
        if(PD->Coupled->s[id_Dead-1]==49) {
            if(t_M>0) {
                Statsnow.Dist_Mitosis[t_M]++;
                Statsnow.Time_M+=t_M;
                Statsnow.N_M++;
            }
            Statsnow.T_live=t;
            return(Statsnow);
        }
        
        if (r0==0){
            if(PD->Coupled->s[id_r-1]==49) {
                r0=1;
                Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
            else t_G1++;
        }
        if (r0==1){ // in Replication
            if(PD->Coupled->s[id_r-1]==49) {t_R++; }
            if(PD->Coupled->s[id_4N-1]==49) {
                r0=2;
                Statsnow.Dist_Replication[t_R]++; Statsnow.Time_R+=t_R; Statsnow.N_R++;
                //  printf("t=%d \t t_R=%d recorded\n",t,t_R);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
        }
        
        if(r0==2){ // 4N DNA accomplished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {t_G2++; }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                //  Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                //  Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
            }
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        
        if(r0==3){ // mitosis, not finished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_T++;
            }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                        r0=3;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        
        if(r0==4){
            if(PD->Coupled->s[id_SAC-1]==48){
                if(t_M>0) {
                    Statsnow.Dist_Mitosis[t_M]++; Statsnow.Time_M+=t_M; Statsnow.N_M++;
                    //printf("t=%d \t t_M=%d recorded\n",t,t_M);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
                if(PD->Coupled->s[id_r-1]==49) {
                    r0=1;
                    if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    t_R++;
                    //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                    //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                    Statsnow.N_endoredupl_T++;
                }
                else t_T++;
            }
            else t_M++;
            
            if((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==49)) {
                r0=4;
                if(t_T>0) {
                    Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;
                    // printf("t=%d \t t_T=%d recorded\n",t,t_T);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
            }
        }
        if (r0==4){
            if ((PD->Coupled->s[id_4N-1]==48)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==48)){
                r0=0;
                Statsnow.N_cycles++;
                // getchar();
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_G1=1;
            }
        }
        PD->Coupled->synchronously_update_all_Gates();
    }
    
    return(Statsnow);
}


Cell_cycle_stats_in_run count_Stochastic_cell_cycles_Stochastic_HITS(const char inputnode[],
                                                                     double p_input,
                                                                     double p_hits,
                                                                     longint_ARRAY_pair *Hits,
                                                                     int a_st,int a_st_k,
                                                                     MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                                     double p_error,unsigned long long int T_MAX){
    Cell_cycle_stats_in_run Statsnow;
    
    Statsnow.N_cycles=0;
    Statsnow.N_endoredupl_T=0;
    Statsnow.N_endoredupl_G2=0;
    Statsnow.N_endoredupl_Aneup=0;
    
    Statsnow.MAX_INTERVAL=0;
    Statsnow.T_live=T_MAX;
    Statsnow.N_deaths=0;
    Statsnow.Time_G1=0;
    Statsnow.Time_R=0;
    Statsnow.Time_G2=0;
    Statsnow.Time_M=0;
    Statsnow.Time_Telo=0;
    Statsnow.N_G1=0;
    Statsnow.N_R=0;
    Statsnow.N_G2=0;
    Statsnow.N_M=0;
    Statsnow.N_T=0;
    for (int t=0;t<=MAX_ARREST; t++){
        Statsnow.Dist_G1[t]=0;
        Statsnow.Dist_G2[t]=0;
        Statsnow.Dist_Mitosis[t]=0;
        Statsnow.Dist_Replication[t]=0;
        Statsnow.Dist_Telophase[t]=0;
    }
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    //   PD->Coupled->print_state();//getchar();
    
    //  int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    // int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    //   int s_start_GF=PD->Coupled->s[id_GF-1];
    // int s_start_CDL=PD->Coupled->s[id_CDL-1];
    //int s_start=PD->Coupled->s[id_input-1];
    
    int id_input=PD->Coupled->BRN->MDAT->get_node_ID(inputnode);
    int id_r=PD->Coupled->BRN->MDAT->get_node_ID("Replication");
    int id_4N=PD->Coupled->BRN->MDAT->get_node_ID("4N_DNA");
    if (id_4N == 0) id_4N=PD->Coupled->BRN->MDAT->get_node_ID("f4N_DNA");
    int id_M=PD->Coupled->BRN->MDAT->get_node_ID("U_Kinetochores");
    int id_SAC=PD->Coupled->BRN->MDAT->get_node_ID("A_Kinetochores");
    int id_cyto=PD->Coupled->BRN->MDAT->get_node_ID("Ect2");
    if (id_cyto==0) id_cyto=PD->Coupled->BRN->MDAT->get_node_ID("Cytokinesis");
    int id_Dead=PD->Coupled->BRN->MDAT->get_node_ID("CAD");
    
    int r0=0;int t_G1=1;int t_G2=0;int t_M=0;int t_R=0;int t_T=0;
    double r;
    
    for (unsigned long long int t=1;t<=T_MAX; t++) {
        if(Hits!= NULL){
            for(unsigned long int l=1;l<=Hits->N;l++){
                r=rand_real(0,1);
                if(r<p_hits) PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
                // else PD->Coupled->s[Hits->A[l]-1]=48+(1-Hits->B[l]);
            }
        }
        r=rand_real(0,1);
        if(r<p_input)
            PD->Coupled->s[id_input-1]=49;
        else PD->Coupled->s[id_input-1]=48;
        
        
        if(PD->Coupled->s[id_Dead-1]==49) {
            Statsnow.T_live=t;
            return(Statsnow);
        }
        
        if (r0==0){
            if(PD->Coupled->s[id_r-1]==49) {
                r0=1;
                Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
            else t_G1++;
        }
        if (r0==1){ // in Replication
            if(PD->Coupled->s[id_r-1]==49) {t_R++; }
            if(PD->Coupled->s[id_4N-1]==49) {
                r0=2;
                Statsnow.Dist_Replication[t_R]++; Statsnow.Time_R+=t_R; Statsnow.N_R++;
                //  printf("t=%d \t t_R=%d recorded\n",t,t_R);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
        }
        
        if(r0==2){ // 4N DNA accomplished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {t_G2++; }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                //  Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                //  Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        if(r0==3){ // mitosis, not finished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_T++;
            }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                        r0=3;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        if(r0==4){
            if(PD->Coupled->s[id_SAC-1]==48){
                if(t_M>0) {
                    Statsnow.Dist_Mitosis[t_M]++; Statsnow.Time_M+=t_M; Statsnow.N_M++;
                    //printf("t=%d \t t_M=%d recorded\n",t,t_M);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
                if(PD->Coupled->s[id_r-1]==49) {
                    r0=1;
                    if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    t_R++;
                    //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                    //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                    Statsnow.N_endoredupl_T++;
                }
                else t_T++;
            }
            else t_M++;
            
            if((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==49)) {
                r0=4;
                if(t_T>0) {
                    Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;
                    // printf("t=%d \t t_T=%d recorded\n",t,t_T);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
            }
        }
        if (r0==4){
            if ((PD->Coupled->s[id_4N-1]==48)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==48)){
                r0=0;
                Statsnow.N_cycles++;
                // getchar();
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_G1=1;
            }
        }
        PD->Coupled->synchronously_update_all_Gates();
    }
    return(Statsnow);
}

Cell_cycle_stats_in_run count_Stochastic_cell_cycles_Stochastic_HITS_Assync(const char inputnode[],
                                                                            double p_input,
                                                                            double p_hits[],
                                                                            longint_ARRAY_pair *Hits,
                                                                            int a_st,int a_st_k,
                                                                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                                            double p_error,unsigned long long int T_MAX){
    Cell_cycle_stats_in_run Statsnow;
    
    Statsnow.N_cycles=0;
    Statsnow.N_endoredupl_T=0;
    Statsnow.N_endoredupl_G2=0;
    Statsnow.N_endoredupl_Aneup=0;
    
    Statsnow.MAX_INTERVAL=0;
    Statsnow.T_live=T_MAX;
    Statsnow.N_deaths=0;
    Statsnow.Time_G1=0;
    Statsnow.Time_R=0;
    Statsnow.Time_G2=0;
    Statsnow.Time_M=0;
    Statsnow.Time_Telo=0;
    Statsnow.N_G1=0;
    Statsnow.N_R=0;
    Statsnow.N_G2=0;
    Statsnow.N_M=0;
    Statsnow.N_T=0;
    for (int t=0;t<=MAX_ARREST; t++){
        Statsnow.Dist_G1[t]=0;
        Statsnow.Dist_G2[t]=0;
        Statsnow.Dist_Mitosis[t]=0;
        Statsnow.Dist_Replication[t]=0;
        Statsnow.Dist_Telophase[t]=0;
    }
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    //   PD->Coupled->print_state();//getchar();
    
    //  int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    //  int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    //   int s_start_GF=PD->Coupled->s[id_GF-1];
    // int s_start_CDL=PD->Coupled->s[id_CDL-1];
    //int s_start=PD->Coupled->s[id_input-1];
    
    int id_input=PD->Coupled->BRN->MDAT->get_node_ID(inputnode);
    int id_r=PD->Coupled->BRN->MDAT->get_node_ID("Replication");
    int id_4N=PD->Coupled->BRN->MDAT->get_node_ID("4N_DNA");
    if (id_4N == 0) id_4N=PD->Coupled->BRN->MDAT->get_node_ID("f4N_DNA");

    int id_M=PD->Coupled->BRN->MDAT->get_node_ID("U_Kinetochores");
    int id_SAC=PD->Coupled->BRN->MDAT->get_node_ID("A_Kinetochores");
    int id_cyto=PD->Coupled->BRN->MDAT->get_node_ID("Ect2");
    if (id_cyto==0) id_cyto=PD->Coupled->BRN->MDAT->get_node_ID("Cytokinesis");
    int id_Dead=PD->Coupled->BRN->MDAT->get_node_ID("CAD");
    
    int r0=0;int t_G1=1;int t_G2=0;int t_M=0;int t_R=0;int t_T=0;
    double r;
    
    for (unsigned long long int t=1;t<=T_MAX; t++) {
        if((Hits != NULL) && (Hits->N>0)){
            for(unsigned long int l=1;l<=Hits->N;l++){
                r=rand_real(0,1);
                if(r<p_hits[l]) PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
                // else PD->Coupled->s[Hits->A[l]-1]=48+(1-Hits->B[l]);
            }
        }
        r=rand_real(0,1);
        if(r<p_input)
            PD->Coupled->s[id_input-1]=49;
        else PD->Coupled->s[id_input-1]=48;
        
        
        if(PD->Coupled->s[id_Dead-1]==49) {
            Statsnow.T_live=t;
            return(Statsnow);
        }
        
        if (r0==0){
            if(PD->Coupled->s[id_r-1]==49) {
                r0=1;
                Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
            else t_G1++;
        }
        if (r0==1){ // in Replication
            if(PD->Coupled->s[id_r-1]==49) {t_R++; }
            if(PD->Coupled->s[id_4N-1]==49) {
                r0=2;
                Statsnow.Dist_Replication[t_R]++; Statsnow.Time_R+=t_R; Statsnow.N_R++;
                //  printf("t=%d \t t_R=%d recorded\n",t,t_R);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
        }
        
        if(r0==2){ // 4N DNA accomplished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {t_G2++; }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                //  Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                //  Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        if(r0==3){ // mitosis, not finished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_T++;
            }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                        r0=3;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        if(r0==4){
            if(PD->Coupled->s[id_SAC-1]==48){
                if(t_M>0) {
                    Statsnow.Dist_Mitosis[t_M]++; Statsnow.Time_M+=t_M; Statsnow.N_M++;
                    //printf("t=%d \t t_M=%d recorded\n",t,t_M);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
                if(PD->Coupled->s[id_r-1]==49) {
                    r0=1;
                    if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    t_R++;
                    //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                    //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                    Statsnow.N_endoredupl_T++;
                }
                else t_T++;
            }
            else t_M++;
            
            if((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==49)) {
                r0=4;
                if(t_T>0) {
                    Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;
                    // printf("t=%d \t t_T=%d recorded\n",t,t_T);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
            }
        }
        if (r0==4){
            if ((PD->Coupled->s[id_4N-1]==48)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==48)){
                r0=0;
                Statsnow.N_cycles++;
                // getchar();
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_G1=1;
            }
        }
        if(ASSYNC_BIAS == 0 ) PD->Coupled->Asynchronously_update_all_Gates();
        else PD->Coupled->Asynchronously_update_all_Gates(ASSYNC_BIAS,Hits,p_hits);
    }
    return(Statsnow);
}




Cell_cycle_stats_in_run count_cell_cycles_stochastic_input(const char inputnode[],double p_input,
                                                           longint_ARRAY_pair *Hits,
                                                           int a_st,int a_st_k,
                                                           MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                           double p_error,unsigned long long int T_MAX){
    Cell_cycle_stats_in_run Statsnow;
    
    Statsnow.N_cycles=0;
    Statsnow.N_endoredupl_T=0;
    Statsnow.N_endoredupl_G2=0;
    Statsnow.N_endoredupl_Aneup=0;
    
    Statsnow.MAX_INTERVAL=0;
    Statsnow.T_live=T_MAX;
    Statsnow.N_deaths=0;
    Statsnow.Time_G1=0;
    Statsnow.Time_R=0;
    Statsnow.Time_G2=0;
    Statsnow.Time_M=0;
    Statsnow.Time_Telo=0;
    Statsnow.N_G1=0;
    Statsnow.N_R=0;
    Statsnow.N_G2=0;
    Statsnow.N_M=0;
    Statsnow.N_T=0;
    for (int t=0;t<=MAX_ARREST; t++){
        Statsnow.Dist_G1[t]=0;
        Statsnow.Dist_G2[t]=0;
        Statsnow.Dist_Mitosis[t]=0;
        Statsnow.Dist_Replication[t]=0;
        Statsnow.Dist_Telophase[t]=0;
    }
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    //   PD->Coupled->print_state();//getchar();
    
    //  int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    //  int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    //   int s_start_GF=PD->Coupled->s[id_GF-1];
    // int s_start_CDL=PD->Coupled->s[id_CDL-1];
    //int s_start=PD->Coupled->s[id_input-1];
    
    int id_input=PD->Coupled->BRN->MDAT->get_node_ID(inputnode);
    int id_r=PD->Coupled->BRN->MDAT->get_node_ID("Replication");
    int id_4N=PD->Coupled->BRN->MDAT->get_node_ID("4N_DNA");
    if (id_4N == 0) id_4N=PD->Coupled->BRN->MDAT->get_node_ID("f4N_DNA");

    int id_M=PD->Coupled->BRN->MDAT->get_node_ID("U_Kinetochores");
    int id_SAC=PD->Coupled->BRN->MDAT->get_node_ID("A_Kinetochores");
    int id_cyto=PD->Coupled->BRN->MDAT->get_node_ID("Ect2");
    int id_Dead=PD->Coupled->BRN->MDAT->get_node_ID("CAD");
    
    int r0=0;int t_G1=1;int t_G2=0;int t_M=0;int t_R=0;int t_T=0;
    
    for (unsigned long long int t=1;t<=T_MAX; t++) {
        for(unsigned long int l=1;l<=Hits->N;l++){
            PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
        }
        double r=rand_real(0,1);
        if(r<p_input) PD->Coupled->s[id_input-1]=49;
        else PD->Coupled->s[id_input-1]=48;
        //        PD->Coupled->s[id_GF-1]=s_start_GF;
        //        PD->Coupled->s[id_CDL-1]=s_start_CDL;
        
        
        if(PD->Coupled->s[id_Dead-1]==49) {
            Statsnow.T_live=t;
            return(Statsnow);
        }
        
        if (r0==0){
            if(PD->Coupled->s[id_r-1]==49) {
                r0=1;
                Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
            else t_G1++;
        }
        if (r0==1){ // in Replication
            if(PD->Coupled->s[id_r-1]==49) {t_R++; }
            if(PD->Coupled->s[id_4N-1]==49) {
                r0=2;
                Statsnow.Dist_Replication[t_R]++; Statsnow.Time_R+=t_R; Statsnow.N_R++;
                //  printf("t=%d \t t_R=%d recorded\n",t,t_R);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
        }
        
        if(r0==2){ // 4N DNA accomplished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {t_G2++; }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                //  Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                //Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
                t_R++;
            }
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        if(r0==3){ // mitosis, not finished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_T++;
            }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                        r0=3;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        if(r0==4){
            if(PD->Coupled->s[id_SAC-1]==48){
                if(t_M>0) {
                    Statsnow.Dist_Mitosis[t_M]++; Statsnow.Time_M+=t_M; Statsnow.N_M++;
                    //printf("t=%d \t t_M=%d recorded\n",t,t_M);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
                if(PD->Coupled->s[id_r-1]==49) {
                    r0=1;
                    if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    t_R++;
                    //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                    //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                    Statsnow.N_endoredupl_T++;
                }
                else t_T++;
            }
            else t_M++;
            
            if((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==49)) {
                r0=4;
                if(t_T>0) {
                    Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;
                    // printf("t=%d \t t_T=%d recorded\n",t,t_T);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
            }
        }
        if (r0==4){
            if ((PD->Coupled->s[id_4N-1]==48)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==48)){
                r0=0;
                Statsnow.N_cycles++;
                // getchar();
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_G1++;
            }
        }
        PD->Coupled->synchronously_update_all_Gates();
    }
    return(Statsnow);
}

Cell_cycle_stats_in_run count_cell_cycles_stochastic_input(const char inputnode[],double p_input,
                                                           const char inputnode2[],double p_input2,
                                                           longint_ARRAY_pair *Hits,
                                                           int a_st,int a_st_k,
                                                           MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                           double p_error,unsigned long long int T_MAX){
    Cell_cycle_stats_in_run Statsnow;
    
    Statsnow.N_cycles=0;
    Statsnow.N_endoredupl_T=0;
    Statsnow.N_endoredupl_G2=0;
    Statsnow.N_endoredupl_Aneup=0;
    
    Statsnow.MAX_INTERVAL=0;
    Statsnow.T_live=T_MAX;
    Statsnow.N_deaths=0;
    Statsnow.Time_G1=0;
    Statsnow.Time_R=0;
    Statsnow.Time_G2=0;
    Statsnow.Time_M=0;
    Statsnow.Time_Telo=0;
    Statsnow.N_G1=0;
    Statsnow.N_R=0;
    Statsnow.N_G2=0;
    Statsnow.N_M=0;
    Statsnow.N_T=0;
    for (int t=0;t<=MAX_ARREST; t++){
        Statsnow.Dist_G1[t]=0;
        Statsnow.Dist_G2[t]=0;
        Statsnow.Dist_Mitosis[t]=0;
        Statsnow.Dist_Replication[t]=0;
        Statsnow.Dist_Telophase[t]=0;
    }
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    //   PD->Coupled->print_state();//getchar();
    
    //  int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    // int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    //   int s_start_GF=PD->Coupled->s[id_GF-1];
    // int s_start_CDL=PD->Coupled->s[id_CDL-1];
    //int s_start=PD->Coupled->s[id_input-1];
    
    int id_input=PD->Coupled->BRN->MDAT->get_node_ID(inputnode);
    int id_input2=PD->Coupled->BRN->MDAT->get_node_ID(inputnode2);
    int id_r=PD->Coupled->BRN->MDAT->get_node_ID("Replication");
    int id_4N=PD->Coupled->BRN->MDAT->get_node_ID("4N_DNA");
    if (id_4N == 0) id_4N=PD->Coupled->BRN->MDAT->get_node_ID("f4N_DNA");

    int id_M=PD->Coupled->BRN->MDAT->get_node_ID("U_Kinetochores");
    int id_SAC=PD->Coupled->BRN->MDAT->get_node_ID("A_Kinetochores");
    int id_cyto=PD->Coupled->BRN->MDAT->get_node_ID("Ect2");
    int id_Dead=PD->Coupled->BRN->MDAT->get_node_ID("CAD");
    
    int r0=0;int t_G1=1;int t_G2=0;int t_M=0;int t_R=0;int t_T=0;
    
    for (unsigned long long int t=1;t<=T_MAX; t++) {
        for(unsigned long int l=1;l<=Hits->N;l++){
            PD->Coupled->s[Hits->A[l]-1]=48+Hits->B[l];
        }
        double r=rand_real(0,1);
        if(r<p_input) PD->Coupled->s[id_input-1]=49;
        else PD->Coupled->s[id_input-1]=48;
        
        r=rand_real(0,1);
        if(r<p_input2) PD->Coupled->s[id_input2-1]=49;
        else PD->Coupled->s[id_input2-1]=48;
        //        PD->Coupled->s[id_GF-1]=s_start_GF;
        //        PD->Coupled->s[id_CDL-1]=s_start_CDL;
        
        
        if(PD->Coupled->s[id_Dead-1]==49) {
            Statsnow.T_live=t;
            return(Statsnow);
        }
        
        if (r0==0){
            if(PD->Coupled->s[id_r-1]==49) {
                r0=1;
                Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
            else t_G1++;
        }
        if (r0==1){ // in Replication
            if(PD->Coupled->s[id_r-1]==49) {t_R++; }
            if(PD->Coupled->s[id_4N-1]==49) {
                r0=2;
                Statsnow.Dist_Replication[t_R]++; Statsnow.Time_R+=t_R; Statsnow.N_R++;
                //  printf("t=%d \t t_R=%d recorded\n",t,t_R);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
        }
        
        if(r0==2){ // 4N DNA accomplished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {t_G2++; }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                //  Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                        r0=3;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        if(r0==3){ // mitosis, not finished
            if((PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==48)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48))  {
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_T++;
            }
            
            if((t_M>0) &&(PD->Coupled->s[id_4N-1]==49)
               &&(PD->Coupled->s[id_r-1]==49)
               &&(PD->Coupled->s[id_M-1]==48)
               &&(PD->Coupled->s[id_SAC-1]==48)) {
                r0=1;
                if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)
                &&(PD->Coupled->s[id_r-1]==48)
                &&((PD->Coupled->s[id_M-1]==49)
                   ||(PD->Coupled->s[id_SAC-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                        r0=3;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[id_r-1]==49)) {
                r0=1;
                Statsnow.Dist_G2[t_G2]++; Statsnow.Time_G2+=t_G2; Statsnow.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==49)) {r0=4;}
        }
        
        
        if(r0==4){
            if(PD->Coupled->s[id_SAC-1]==48){
                if(t_M>0) {
                    Statsnow.Dist_Mitosis[t_M]++; Statsnow.Time_M+=t_M; Statsnow.N_M++;
                    //printf("t=%d \t t_M=%d recorded\n",t,t_M);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
                if(PD->Coupled->s[id_r-1]==49) {
                    r0=1;
                    if(t_T>0)  {Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;}
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    t_R++;
                    //Statsnow.Dist_G1[t_G1]++; Statsnow.Time_G1+=t_G1; Statsnow.N_G1++;
                    //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                    Statsnow.N_endoredupl_T++;
                }
                else t_T++;
            }
            else t_M++;
            
            if((PD->Coupled->s[id_4N-1]==49)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==49)) {
                r0=4;
                if(t_T>0) {
                    Statsnow.Dist_Telophase[t_T]++; Statsnow.Time_Telo+=t_T; Statsnow.N_T++;
                    // printf("t=%d \t t_T=%d recorded\n",t,t_T);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
            }
        }
        if (r0==4){
            if ((PD->Coupled->s[id_4N-1]==48)&&(PD->Coupled->s[id_SAC-1]==48)&&(PD->Coupled->s[id_cyto-1]==48)){
                r0=0;
                Statsnow.N_cycles++;
                // getchar();
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_G1++;
            }
        }
        PD->Coupled->synchronously_update_all_Gates();
    }
    return(Statsnow);
}

void run_TwoNode_pulse_experiment_MOD(const char nodename1[],const char nodename2[],int a_st,int a_st_k,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                      double p_error,int pulse_l){
    int W=8,YTOP,XTOP;
    long int id_new;
    FILE *f;
    char fn[400];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    //   PD->Coupled->print_state();//getchar();
    
    int id_node_1=PD->BRN->MDAT->get_node_ID(nodename1);
    int id_node_2=PD->BRN->MDAT->get_node_ID(nodename2);
    
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+3*W)+5;
    XTOP=2*TIMESTEPS_DEFAULT+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->name);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s_and_%s_pulse_%d\n",MODEL_Directory_system,PD->BRN->name,nodename1,nodename2,pulse_l);
    system(fn);
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/%s_and_%s_pulse_%d/%s_and_%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d.eps",MODEL_Directory,PD->BRN->name,nodename1,nodename2,
                pulse_l,nodename1,nodename2,a_st,a_st_k,pulse_l,delay);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
        id_new=draw_timetrace_INPUT_twopulse_MOD(a_st,a_st_k,PD,p_error,delay,pulse_l,id_node_1,id_node_2,f,YTOP,W);
        printf("When flipping %s and %s for %d steps - delay %d, attractor %d went to %ld\n",
               nodename1,nodename2,pulse_l,delay,a_st,id_new);
        //getchar();
        sprintf(fn,"%s/%s/_EXP/%s_and_%s_pulse_%d/%s_and_%s_pulse_experiment_A%d_As-%d__pulse-%d_delay-%d",MODEL_Directory_system,PD->BRN->name,nodename1,nodename2,pulse_l,nodename1,nodename2,a_st,a_st_k,pulse_l,delay);
        close_EPSDRAW_file(f,fn,2);
    }
    
}

void run_OnePulse_analysis_on_ONE_attractor(const char nodename[],int attr_id,
                                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Boolean_Dynamics_TRACKER *D;
    Phenotype_now P1,P2;
    
    D=PD->Coupled;
    
    int MAX_PULSE=60;
    
    int id_pulse=PD->BRN->MDAT->get_node_ID(nodename);
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    int *b_new;
    b_new=new int[MAX_PULSE+1];
    
    //for(unsigned int b=1;b<=D->Attractor_NR;b++) {
    unsigned int b=attr_id;
    for (int a_st_k=1; a_st_k<=D->Attractor_valleys[b]->N_a; a_st_k++) {
        PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[a_st_k]]->s);
        //    PD->Coupled->print_state();//getchar();
        P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        
        int s_start=PD->Coupled->s[id_pulse-1];
        int s_start_GF=PD->Coupled->s[id_GF-1];
        //int s_start_CDL=PD->Coupled->s[id_CDL-1];
        
        int olda=b;
        for (int pulse_length=1; pulse_length<=MAX_PULSE; pulse_length++) {
            PD->Coupled->s[id_pulse-1]=48+(49-s_start);
            PD->Coupled->synchronously_update_all_Gates();
            
            PD->Coupled->s[id_pulse-1]=s_start;
            PD->Coupled->s[id_GF-1]=s_start_GF;
            
            b_new[pulse_length]=(int)PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
            P2=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b_new[pulse_length],1);
            
            if (( b_new[pulse_length]!=b) && ( b_new[pulse_length]!=olda)) {
                if (s_start==48) printf("%s ON:",nodename); else printf("%s OFF:",nodename);
                printf(" %d (%d) --> %d (pulse=%d)\t\t",b,a_st_k, b_new[pulse_length],pulse_length);
                if (P1.Dead>P2.Dead) printf("Resurrected!\t");
                if (P1.Dead<P2.Dead) printf("Died\t");
                if (P1.Sen>P2.Sen) printf("Recovered from Senescence!\t");
                if (P1.Sen<P2.Sen) printf("Senesced\t");
                
                if (P1.CC!=P2.CC) {
                    switch (P1.CC) {
                        case 0:printf("Cell_Cycle ---> ");break;
                        case 1:printf("G0/G1 ---> ");break;
                        case 2:printf("G2 ---> ");break;
                        case 3:printf("SAC ---> ");break;
                        case 4:printf("Broken_Cycle ---> ");break;
                        case 5:printf("PHS_on_barrier ---> ");break;
                        case 6:printf("?? ---> ");break;
                        default: break;
                    }
                    
                    switch (P2.CC) {
                        case 0:printf("Cell_Cycle\t");break;
                        case 1:printf("G0/G1\t");break;
                        case 2:printf("G2\t");break;
                        case 3:printf("SAC\t");break;
                        case 4:printf("Broken_Cycle\t");break;
                        case 5:printf("PHS_on_barrier\t");break;
                        case 6:printf("??\t");break;
                        default: break;
                    }
                }
                
                if (P1.OtherC>P2.OtherC) printf("Stopped OtherCycle!\t");
                if (P1.OtherC<P2.OtherC) printf("Started OtherCycle\t");
                if (P1.ESC!=P2.ESC) {
                    switch (P1.ESC) {
                        case 0:printf("Naive ESC ---> ");break;
                        case 1:printf("Primed ESC ---> ");break;
                        case 2:printf("Diff_1 ---> ");break;
                        case 3:printf("Diff_2 ---> ");break;
                        case 5:printf("Neuron ---> ");break;
                        default: break;
                    }
                    switch (P2.ESC) {
                        case 0:printf("Naive ESC\n");break;
                        case 1:printf("Primed ESC\n");break;
                        case 2:printf("Diff_1\n");break;
                        case 3:printf("Diff_2\n");break;
                        case 5:printf("Neuron\n");break;
                        default: break;
                    }
                }
                else printf("\n");
            }
            olda=b_new[pulse_length];
        }
    }
}


void run_All_environmental_changes_from_ONE_attractor(longint_ARRAY *gl_keep_asis, int attr_id,
                                                      MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Boolean_Dynamics_TRACKER *D;
    Phenotype_now P1,P2;
    
    D=PD->Coupled;
    
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
   // int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    int *b_new;
    b_new=new int[gl_keep_asis->N+1];
    
    unsigned int b=attr_id;
    for (int a_st_k=1; a_st_k<=D->Attractor_valleys[b]->N_a; a_st_k++) {
        PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[a_st_k]]->s);
        //    PD->Coupled->print_state();//getchar();
        P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        int olda=b;
        
        for (int idp=1; idp<=gl_keep_asis->N;idp++) {
            PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[a_st_k]]->s);
            //    PD->Coupled->print_state();//getchar();
            int s_start=PD->Coupled->s[gl_keep_asis->A[idp]-1];
            int s_start_GF=PD->Coupled->s[id_GF-1];
           // int s_start_CDL=PD->Coupled->s[id_CDL-1];
            
            PD->Coupled->s[gl_keep_asis->A[idp]-1]=48+(49-s_start);
            PD->Coupled->synchronously_update_all_Gates();
            if(gl_keep_asis->A[idp]!=id_GF) PD->Coupled->s[id_GF-1]=s_start_GF;
            
            b_new[idp]=(int)PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
            P2=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b_new[idp],1);
            
            if (( b_new[idp]!=b) && ( b_new[idp]!=olda)) {
                if (s_start==48) printf("%s ON:",PD->BRN->MDAT->Node_name[gl_keep_asis->A[idp]]);
                else printf("%s OFF:",PD->BRN->MDAT->Node_name[gl_keep_asis->A[idp]]);
                printf(" %d (%d) --> %d \t\t",b,a_st_k,b_new[idp]);
                if (P1.Dead>P2.Dead) printf("Resurrected!\t");
                if (P1.Dead<P2.Dead) printf("Died\t");
                if (P1.Sen>P2.Sen) printf("Recovered from Senescence!\t");
                if (P1.Sen<P2.Sen) printf("Senesced\t");
                
                if (P1.CC!=P2.CC) {
                    switch (P1.CC) {
                        case 0:printf("Cell_Cycle ---> ");break;
                        case 1:printf("G0/G1 ---> ");break;
                        case 2:printf("G2 ---> ");break;
                        case 3:printf("SAC ---> ");break;
                        case 4:printf("Broken_Cycle ---> ");break;
                        case 5:printf("PHS_on_barrier ---> ");break;
                        case 6:printf("?? ---> ");break;
                        default: break;
                    }
                    
                    switch (P2.CC) {
                        case 0:printf("Cell_Cycle\t");break;
                        case 1:printf("G0/G1\t");break;
                        case 2:printf("G2\t");break;
                        case 3:printf("SAC\t");break;
                        case 4:printf("Broken_Cycle\t");break;
                        case 5:printf("PHS_on_barrier\t");break;
                        case 6:printf("??\t");break;
                        default: break;
                    }
                }
                
                if (P1.OtherC>P2.OtherC) printf("Stopped OtherCycle!\t");
                if (P1.OtherC<P2.OtherC) printf("Started OtherCycle\t");
                if (P1.ESC!=P2.ESC) {
                    switch (P1.ESC) {
                        case 0:printf("Naive ESC ---> ");break;
                        case 1:printf("Primed ESC ---> ");break;
                        case 2:printf("Endoderm ---> ");break;
                        case 3:printf("Trophoectoderm ---> ");break;
                        case 4:printf("Ectoderm ---> ");break;
                        case 5:printf("Neuron ---> ");break;
                        case 6:printf("Weird, undiff. ---> ");break;
                        default: break;
                    }
                    switch (P2.ESC) {
                        case 0:printf("Naive ESC\n");break;
                        case 1:printf("Primed ESC\n");break;
                        case 2:printf("Endoderm\n");break;
                        case 3:printf("Trophoectoderm\n");break;
                        case 4:printf("Ectoderm\n");break;
                        case 5:printf("Neuron\n");break;
                        case 6:printf("Weird, undiff.\n");break;
                        default: break;
                    }
                }
                else printf("\n");
            }
        }
    }
}


class Attractors_on_MAP {
public:
    unsigned short int *shape,*col,*border,*bar,*dot,*celltype;
    longint_ARRAY *al;
    unsigned int *a_order,*at_poz;
    
    Attractors_on_MAP(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,unsigned long int b){
        al=new longint_ARRAY();
        al->add_element(b);
        
        a_order=NULL;at_poz=NULL;
        shape = new unsigned short int [al->N+1];
        col = new unsigned short int [al->N+1];
        border = new unsigned short int [al->N+1];
        bar = new unsigned short int [al->N+1];
        dot = new unsigned short int [al->N+1];
        celltype =new unsigned short int [al->N+1];
        for(int bb=1;bb<=al->N;bb++){
            shape[bb]=0;
            col[bb]=0;
            border[bb]=0;
            bar[bb]=0;
            dot[bb]=0;
            celltype[bb]=0;
        }
        
        Phenotype_now P1;
        for(int bb=1;bb<=al->N;bb++){
            P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,al->A[bb],1);
            
            if (P1.Dead==1) border[bb]=1; else border[bb]=0;
            if (P1.G2_struck==1) bar[bb]=1; else bar[bb]=0;
            dot[bb]=P1.DNA_damage;
            celltype[bb]=P1.ESC;
          //  if (P1.Migr==1) col[bb]=10; break; // migration
          //  if (P1.Migr==2) col[bb]=11; break; // migration
            if (P1.Sen>0) {
                switch (P1.Sen) {
                    case 3:col[bb]=1; break; // chrom
                    case 2:col[bb]=2; break;// mito
                    case 1:col[bb]=3; break; // both
                    default: col[bb]=6; break; // weird..., should not happen
                }
            }
            switch (P1.CC) {
                case 0: {shape[bb]=1;  if(P1.Sen==0) col[bb]=7;} break; //CC
                case 1: {shape[bb]=2;  if(P1.Sen==0) col[bb]=4;} break; //G0_1
                case 2: {shape[bb]=3;  if(P1.Sen==0) col[bb]=5;} break;    // G2
                case 3: {shape[bb]=4;  if(P1.Sen==0) col[bb]=6;} break; // SAC
                case 4: {shape[bb]=5;  if(P1.Sen==0) col[bb]=8;} break; // Broken_Cycle
                case 5: {shape[bb]=6;  if(P1.Sen==0) col[bb]=4;} break; //PHS_on_barrier
                case 6: {shape[bb]=5;  if(P1.Sen==0) col[bb]=8;} break; // ??
                default:{shape[bb]=5;  if(P1.Sen==0) col[bb]=8;} break;
            }
        }
    }
    
    Attractors_on_MAP(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,longint_ARRAY *alin){
        al=new longint_ARRAY();
        al->add_longint_ARRAY(alin);
        
        a_order=NULL;at_poz=NULL;
        shape = new unsigned short int [al->N+1];
        col = new unsigned short int [al->N+1];
        border = new unsigned short int [al->N+1];
        bar = new unsigned short int [al->N+1];
        dot = new unsigned short int [al->N+1];
        celltype =new unsigned short int [al->N+1];
        for(int bb=1;bb<=al->N;bb++){
            shape[bb]=0;
            col[bb]=0;
            border[bb]=0;
            bar[bb]=0;
            dot[bb]=0;
            celltype[bb]=0;
        }
        
        Phenotype_now P1;
        for(int bb=1;bb<=al->N;bb++){
            P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,al->A[bb],1);
            
            if (P1.Dead==1) border[bb]=1; else border[bb]=0;
            if (P1.G2_struck==1) bar[bb]=1; else bar[bb]=0;
            dot[bb]=P1.DNA_damage;
            celltype[bb]=P1.ESC;
           // if (P1.Migr==1) col[bb]=10; break; // migration
           // if (P1.Migr==2) col[bb]=11; break; // migration
            if (P1.Sen>0) {
                switch (P1.Sen) {
                    case 3:col[bb]=1; break; // chrom
                    case 2:col[bb]=2; break;// mito
                    case 1:col[bb]=3; break; // both
                    default: col[bb]=6; break; // weird..., should not happen
                }
            }
            switch (P1.CC) {
                case 0: {shape[bb]=1;  if(P1.Sen==0) col[bb]=7;} break; //CC
                case 1: {shape[bb]=2;  if(P1.Sen==0) col[bb]=4;} break; //G0_1
                case 2: {shape[bb]=3;  if(P1.Sen==0) col[bb]=5;} break;    // G2
                case 3: {shape[bb]=4;  if(P1.Sen==0) col[bb]=6;} break; // SAC
                case 4: {shape[bb]=5;  if(P1.Sen==0) col[bb]=8;} break; // Broken_Cycle
                case 5: {shape[bb]=6;  if(P1.Sen==0) col[bb]=4;} break; //PHS_on_barrier
                case 6: {shape[bb]=5;  if(P1.Sen==0) col[bb]=8;} break; // ??
                default:{shape[bb]=5;  if(P1.Sen==0) col[bb]=8;} break;
            }
        }
    }
    
    Attractors_on_MAP(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,longint_ARRAY *yztu_inputs,int id_GF,int id_GFH,int i_now,int j_now,int k_now,int l_now,int m_now){
        Phenotype_now P1;
        
        al=new longint_ARRAY();
        for(int b=1;b<=PD->Coupled->Attractor_NR;b++){
            switch (i_now) {
                case 0:
                    if ((PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GF-1]==48)&&(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GFH-1]==48))
                        al->add_element(b);
                    break;
                case 1:
                    if ((PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GF-1]==49)&&(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GFH-1]==48))
                        al->add_element(b);
                    break;
                case 2:{
                    if ((PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GF-1]==48)&&(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GFH-1]==49))  al->add_element(b);
                    if ((PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GF-1]==49)&&(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GFH-1]==49))  al->add_element(b);
                }
                    break;
                    
                default:
                    break;
            }
            // printf("b=%d\t So far: %ld\n",b,al->N);//getchar();
        }
        
        for(int b=1;b<=PD->Coupled->Attractor_NR;b++){
            if (j_now==0)  if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[1]-1]==49) al->delete_value(b);
            if (j_now==1)  if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[1]-1]==48) al->delete_value(b);
            
            // printf("b=%d\t So far: %ld\n",b,al->N);//getchar();
            
            if(yztu_inputs->N>=2){
                if (k_now==0)  if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[2]-1]==49) al->delete_value(b);
                if (k_now==1)  if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[2]-1]==48) al->delete_value(b);
            }
            // printf("b=%d\t So far: %ld\n",b,al->N);//getchar();
            
            if(yztu_inputs->N>=3){
                if (l_now==0)  if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[3]-1]==49) al->delete_value(b);
                if (l_now==1)  if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[3]-1]==48) al->delete_value(b);
            }
            // printf("b=%d\t So far: %ld\n",b,al->N);//getchar();
            
            if(yztu_inputs->N==4){
                if (m_now==0)  if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[4]-1]==49) al->delete_value(b);
                if (m_now==1)  if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[4]-1]==48) al->delete_value(b);
            }
            // printf("b=%d\t So far: %ld\n",b,al->N);//getchar();
            
        }
       
        if(al->N ==0) {
            printf("%lu attractors found in i_now = %d, j_now = %d, k_now = %d, l_now = %d, m_now = %d corner!\n",al->N,i_now, j_now, k_now, l_now, m_now);
            getchar();
        }
        a_order=NULL;at_poz=NULL;
        shape = new unsigned short int [al->N+1];
        col = new unsigned short int [al->N+1];
        border = new unsigned short int [al->N+1];
        bar = new unsigned short int [al->N+1];
        dot = new unsigned short int [al->N+1];
        celltype =new unsigned short int [al->N+1];
        for(int bb=1;bb<=al->N;bb++){
            shape[bb]=0;
            col[bb]=0;
            border[bb]=0;
            bar[bb]=0;
            dot[bb]=0;
            celltype[bb]=0;
        }
        
        for(int bb=1;bb<=al->N;bb++){
            P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,al->A[bb],1);
            
            if (P1.Dead==1) border[bb]=1; else border[bb]=0;
            if (P1.G2_struck==1) bar[bb]=1; else bar[bb]=0;
            dot[bb]=P1.DNA_damage;
            celltype[bb]=P1.ESC;
           // if (P1.Migr==1) col[bb]=10; break; // migration
           // if (P1.Migr==2) col[bb]=11; break; // migration
            if (P1.Sen>0) {
                switch (P1.Sen) {
                    case 3:col[bb]=1; break; // chrom
                    case 2:col[bb]=2; break;// mito
                    case 1:col[bb]=3; break; // both
                    default: col[bb]=6; break; // weird..., should not happen
                }
            }
            
            switch (P1.CC) {
                case 0: {shape[bb]=1;  if(P1.Sen==0) col[bb]=7;} break; //CC
                case 1: {shape[bb]=2;  if(P1.Sen==0) col[bb]=4;} break; //G0_1
                case 2: {shape[bb]=3;  if(P1.Sen==0) col[bb]=5;} break;    // G2
                case 3: {shape[bb]=4;  if(P1.Sen==0) col[bb]=6;} break; // SAC
                case 4: {shape[bb]=5;  if(P1.Sen==0) col[bb]=8;} break; // Broken_Cycle
                case 5: {shape[bb]=6;  if(P1.Sen==0) col[bb]=4;} break; //PHS_on_barrier
                case 6: {shape[bb]=5;  if(P1.Sen==0) col[bb]=8;} break; // ??
                default:{shape[bb]=5;  if(P1.Sen==0) col[bb]=8;} break;
            }
        }
        sort_attractors();
    }
    
    int swap_two(int b1,int b2){
        
        if(border[b1]<border[b2]) return (1);
        if(border[b1]==border[b2]){
            if(col[b1]>col[b2]) return (1);
            if(col[b1]==col[b2]){
                if(celltype[b1]<celltype[b2]) return (1);
                if(celltype[b1]==celltype[b2]){
                    if(shape[b1]>shape[b2]) return (1);
                    if(shape[b1]==shape[b2]){
                        if(bar[b1]>bar[b2]) return(1);
                    }
                }
            }
        }
        return (0);
    }
    
    void sort_attractors(){
        // adaptation of bubble sort for attractors; ordering operation is "swap_two"
        a_order=new unsigned int[al->N+1];
        at_poz=new unsigned int[al->N+1];
        
        for(int x=1; x<=al->N; x++) {
            a_order[x]=x;
            at_poz[x]=x;
        }
        int swapnow=1;
        while (swapnow) {
            swapnow=0;
            for(int y=1; y<=al->N-1; y++){
                if(swap_two(at_poz[y], at_poz[y+1])){
                    int temp = at_poz[y+1];
                    at_poz[y+1] = at_poz[y];
                    at_poz[y] = temp;
                    swapnow=1;
                    //  printf("\nswapped poz %d (%d -> a %ld) with %d (%d -> a %ld)\n",
                    //         y,at_poz[y],al->A[at_poz[y]],
                    //         y+1,at_poz[y+1],al->A[at_poz[y]+1]);//getchar();
                }
            }
        }
        //  for(int x=1; x<=al->N; x++)
        //     printf("poz:%d -> index =%u (Attr=%ld)\n",x, at_poz[x],al->A[at_poz[x]]);
        //getchar();
        
    }
    
    ~Attractors_on_MAP(){
        delete[] shape;     shape=NULL;     delete[] col;   col =NULL;
        delete[] border;    border=NULL;    delete[] bar;   bar=NULL;
        if(al!=NULL) {delete  al;al=NULL;}
        if(a_order!=NULL) {delete a_order;a_order=NULL;}
        if(at_poz!=NULL) {delete at_poz;at_poz=NULL;}
    }
};



void make_3D_map_of_steady_states(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,longint_ARRAY *yztu_inputs,longint_ARRAY *gl_keep_asis, int inputnow){
    char fn[300], bla[100];
    FILE *f;
    double a ,alpha, beta, TW,xd,yd,aa,bb;
    
    Boolean_Dynamics_TRACKER *D;
    Phenotype_now P1,P2;
    
    D=PD->Coupled;
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    //int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    int id_GFH=PD->BRN->MDAT->get_node_ID(GF_high_name);
   // int id_CDH=PD->BRN->MDAT->get_node_ID(CD_high_name);
    int id_dead=PD->BRN->MDAT->get_node_ID("CAD");
    int *b_new;
    b_new=new int[gl_keep_asis->N+1];
    
    
    a=100;
    alpha=a*0.9;beta=a*0.65;
    xd=0.4;yd=0.55;
    TW=30;
    aa=a/10.; bb=0.5*a/10.;
    double xnow, ynow;
    int shape=0;
    double sc_node=0.4;
    
    switch (yztu_inputs->N) {
        case 1:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s__Change_%s.eps", MODEL_Directory,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 2:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s__Change_%s.eps", MODEL_Directory,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 3:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s_%s__Change_%s.eps", MODEL_Directory,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[3]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 4:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s_%s_%s__Change_%s.eps", MODEL_Directory,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[3]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[4]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        default: exit(1);
            break;
    }
    
    
    printf("%s\n",fn);
    
    f=fopen(fn,"w");
    
    COLOR_SCHEME=GYR;
    setcolors_PAJEK();
    
    double XMAX=5*a+2*alpha+2*TW;
    double YMAX=3*a+2*beta+2*TW;
    
    //printf("XMIN=%lg\tXMAX=%lg\nYMIN=%lg\tYMAX=%lg\n",XMIN,XMAX,YMIN,YMAX);
    //getchar();
    EpsInit(f,-25,-25,XMAX+25,YMAX+25);
    
    EpsSetRgb(f,0.2,0.2,0.2);
    EpsSetLinewidth(f, 1);
    EpsSetFont(f,"Times-Roman",10);
    
    fivepointer *x,*y;
    x=new fivepointer[3];y=new fivepointer[3];
    for (int i=0; i<3; i++) {
        x[i]=new fourpointer[2];y[i]=new fourpointer[2];
        for (int j=0; j<2; j++) {
            x[i][j]=new triplepointer[2];y[i][j]=new triplepointer[2];
            for (int k=0; k<2; k++) {
                x[i][j][k]=new doublepointer[2];y[i][j][k]=new doublepointer[2];
                for (int l=0; l<2; l++) {
                    x[i][j][k][l]=new double[2];y[i][j][k][l]=new double[2];
                }
            }
        }
    }
    // shift by x coordinate
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            for (int k=0; k<2; k++) {
                for (int l=0; l<2; l++) {
                    x[i][j][k][l][0]=TW+i*a + j*xd*a       + l*(2*a+alpha);
                    x[i][j][k][l][1]=TW+i*a + j*xd*a       + l*(2*a+alpha) + alpha/2.;
                    y[i][j][k][l][0]=TW     + j*yd*a + k*a + l* beta;
                    y[i][j][k][l][1]=TW     + j*yd*a + k*a + l* beta + 1.6*a+ beta;
                }
            }
        }
    }
    
    // horizontals
    for (int j=0; j<2; j++) {
        for (int k=0; k<2; k++) {
            for (int l=0; l<2; l++) {
                for (int m=0; m<2; m++){
                    EpsDrawLine(f, x[0][j][k][l][m],y[0][j][k][l][m],x[2][j][k][l][m],y[2][j][k][l][m]);
                    if ((j==0)&&(k==0)) {
                        EpsDrawString_Centered(f, x[0][j][k][l][m],y[0][j][k][l][m]-a/5., "No GF");
                        EpsDrawString_Centered(f, x[1][j][k][l][m],y[0][j][k][l][m]-a/5., "Low GF");
                        EpsDrawString_Centered(f, x[2][j][k][l][m],y[0][j][k][l][m]-a/5., "High GF");
                    }
                }
            }
        }
    }
    // verticals
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            for (int l=0; l<2; l++) {
                for (int m=0; m<2; m++){
                    EpsDrawLine(f, x[i][j][0][l][m],y[i][j][0][l][m],x[i][j][1][l][m],y[i][j][1][l][m]);
                    if ((i==0)&&(j==0)) {
                        EpsDrawString_Centered(f,  (x[0][0][0][l][m]+x[0][0][1][l][m])/2.-a/6.,
                                               (y[0][0][0][l][m]+y[0][0][1][l][m])/2.,PD->BRN->MDAT->Node_name[yztu_inputs->A[1]]);
                    }
                }
            }
        }
    }
    
    if (yztu_inputs->N>1){
        // diagonals
        for (int i=0; i<3; i++) {
            for (int k=0; k<2; k++) {
                for (int l=0; l<2; l++) {
                    for (int m=0; m<2; m++){
                        EpsDrawLine(f, x[i][0][k][l][m],y[i][0][k][l][m],x[i][1][k][l][m],y[i][1][k][l][m]);
                        if ((i==0)&&(k==0)) {
                            EpsDrawString_Centered(f,  (x[0][0][0][l][m]+x[0][1][0][l][m])/2.-a/6.,
                                                   (y[0][0][0][l][m]+y[0][1][0][l][m])/2., PD->BRN->MDAT->Node_name[yztu_inputs->A[2]]);
                        }
                    }
                }
            }
        }
    }
    
    if (yztu_inputs->N>2){
        EpsSetLinewidth(f, 1.2);
        EpsDrawLine(f,  x[0][1][1][0][1]-6*aa,     y[0][1][1][0][1]+2*aa,
                    x[0][1][1][0][1]+2*aa+2*a, y[0][1][1][0][1]+2*aa);
        sprintf(bla, "No %s",PD->BRN->MDAT->Node_name[yztu_inputs->A[3]]);
        EpsDrawString_Centered(f, (x[0][1][1][0][1]-6*aa+x[0][1][1][0][1]+2*aa+2*a)/2.,
                               (y[0][1][1][0][1]+2*aa+y[0][1][1][0][1]+2*aa)/2.+aa, bla);
        
        EpsDrawLine(f,  x[0][1][1][1][1]-6*aa,     y[0][1][1][1][1]+2*aa,
                    x[0][1][1][1][1]+2*aa+2*a, y[0][1][1][1][1]+2*aa);
        EpsDrawString_Centered(f, (x[0][1][1][1][1]-6*aa+x[0][1][1][1][1]+2*aa+2*a)/2.,
                               (y[0][1][1][1][1]+2*aa+y[0][1][1][1][1]+2*aa)/2.+aa, PD->BRN->MDAT->Node_name[yztu_inputs->A[3]]);
    }
    
    if (yztu_inputs->N==4) {
        EpsDrawLine(f, x[0][0][0][0][0]-3*aa, y[0][0][0][0][0]-2*aa, x[0][0][1][0][0]-3*aa, y[0][0][1][0][0]+4*aa);
        sprintf(bla, "No %s",PD->BRN->MDAT->Node_name[yztu_inputs->A[4]]);
        EpsDrawString_Centered(f, (x[0][0][0][0][0]+x[0][0][1][0][0]-6*aa)/2.-2*aa,
                               (y[0][0][0][0][0]+y[0][0][1][0][0]+2*aa)/2., bla);
        
        // EpsDrawCircle(f, x[0][0][0][0][1], y[0][0][0][0][1], 20);
        EpsDrawLine(f, x[0][0][0][0][1]-3*aa, y[0][0][0][0][1]-2*aa, x[0][0][1][0][1]-3*aa, y[0][0][1][0][1]+4*aa);
        EpsDrawString_Centered(f,   (x[0][0][0][0][1]+x[0][0][1][0][1]-6*aa)/2.-2*aa,
                               (y[0][0][0][0][1]+y[0][0][1][0][1]+2*aa)/2., PD->BRN->MDAT->Node_name[yztu_inputs->A[4]]);
    }
    
    int i=0,j=0,k=0,l=0,m=0;
    
    for(int b=1;b<=PD->Coupled->Attractor_NR;b++){
        P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        
        if ((PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GF-1]==48)&&(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GFH-1]==48)) i=0;
        if ((PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GF-1]==48)&&(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GFH-1]==49)) i=2;
        if ((PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GF-1]==49)&&(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GFH-1]==49)) i=2;
        if ((PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GF-1]==49)&&(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[id_GFH-1]==48)) i=1;
        
        if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[1]-1]==48) k=0; else k=1;
        
        if(yztu_inputs->N>1){
            if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[2]-1]==48) j=0;
            else j=1;
        }
        if(yztu_inputs->N>2){
            if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[3]-1]==48) l=0;
            else l=1;
        }
        
        if(yztu_inputs->N==4){
            if (PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[yztu_inputs->A[4]-1]==48) m=0;
            else m=1;
        }
        else m=0;
        
        xnow=x[i][j][k][l][m];ynow=y[i][j][k][l][m];
        
        EpsSetRgb(f,0,0,0);
        EpsSetLinewidth(f, 0.25);
        EpsDrawRectangle(f,xnow-1.5*aa,ynow-1.5*aa,xnow+1.5*aa,ynow+1.5*aa);
        
        if (P1.Dead==1)     {xnow-=bb/4.;   ynow+=5*aa/6.;  shape=0;}
        else if (P1.Sen==1) {xnow+=bb/4.;   ynow+=5*aa/6.; shape=1;}
        switch (P1.CC) {
            case 0: {shape=2;} break; //CC
            case 1: if(P1.Dead==0) {xnow-=bb/4.;shape=3;}  break; //G0_1
            case 2: if(P1.Dead==0) {ynow-=5*aa/6.; shape=4;} break;    // G2
            case 3: if(P1.Dead==0) {xnow+=bb/4.; ynow-=5*aa/6.;shape=5;} break;
            case 4: if(P1.Dead==0) {xnow+=bb/4.; ynow-=5*aa/6.;  shape=5;} break;
            case 5: if(P1.Dead==0) {xnow+=bb/4.; ynow-=5*aa/6.;  shape=5;} break;
            case 6: if(P1.Dead==0) {xnow+=bb/4.; ynow-=5*aa/6.; shape=5;} break;
            default: break;
        }
        
        //   if (AL[b]->Phenotype.OtherC==1) printf("OtherCycle; ");
        
        if (P1.ESC!=-1){
            switch (P1.ESC) {
                case 0:{ ynow+=3*bb/3.;xnow+=3*bb/2.;
                    EpsSetColor_PAJEK(f, szinnev[17]);} break; //Naive_ESC // "LightOrange")
                case 1:{ ynow+=bb/3.; xnow+=bb/2.;
                    EpsSetColor_PAJEK(f, szinnev[10]);} break; //Primed_ESC // Plum
                case 2:{ ynow-=bb/3.; xnow-=bb/2.;
                    EpsSetColor_PAJEK(f, szinnev[2]);} break; //blue //Endoderm
                case 3:{ ynow-=5*aa/3.;xnow-=5*aa/2.;
                    EpsSetColor_PAJEK(f, szinnev[3]); } break; // Green //Trophectoderm
                case 4:{ ynow-=3*bb/3.;xnow-=3*bb/2.;
                    EpsSetColor_PAJEK(f, szinnev[1]);} break; //red //Ectoderm
                case 5:{ ynow-=7*bb/3.;xnow-=7*bb/2.;
                    EpsSetRgb(f, 0.8, 0.8, 0.8);} break; //weird diff
                default: EpsSetRgb(f, 0.9, 0.9, 0.9); break;
            }
        }
        
        PD->Coupled->Attractor_valleys[b]->x=xnow;
        PD->Coupled->Attractor_valleys[b]->y=ynow;
        
        //printf("Dead = %d \t shape =%d \t x= %lg  y = %lg\n",P1.Dead,shape,xnow,ynow);//getchar();
        switch (shape) {
            case 0: {EpsFillTriangle(f, xnow-aa*sc_node,ynow-bb*sc_node,xnow+aa*sc_node,ynow+bb*sc_node,xnow,ynow+bb*sc_node);} //dead
                break;
                
            case 1: {EpsFillRectangle(f, xnow-aa*sc_node,ynow-bb*sc_node/2.,xnow+aa*sc_node,ynow+bb*sc_node/2.);} //senescent
                break;
                
            case 4: {EpsFillRectangle(f, xnow-bb*sc_node,ynow-bb*sc_node,xnow+bb*sc_node,ynow+bb*sc_node);}
                //G2
                break;
                
            case 3: {EpsFillCircle(f, xnow, ynow, aa*sc_node);} //G0_G1
                break;
                
            case 2: {EpsFillEllipse(f, xnow, ynow, aa*sc_node, bb*sc_node);} //Cellcycle
                break;
                
            case 5: {EpsFillRectangle(f, xnow-bb*sc_node,ynow-aa*sc_node,xnow+bb*sc_node,ynow+aa*sc_node);} //sen
                break;
                
            default:
                break;
        }
        EpsSetFont(f,"Times-Roman",3);
        EpsSetRgb(f, 0,0,0);
        sprintf(bla, "%d",b);
        EpsDrawString_Centered(f, PD->Coupled->Attractor_valleys[b]->x, PD->Coupled->Attractor_valleys[b]->y, bla);
        
    }
    
    for(int b=1;b<=PD->Coupled->Attractor_NR;b++){
        P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        for (int a_st_k=1; a_st_k<=D->Attractor_valleys[b]->N_a; a_st_k++) {
            PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[a_st_k]]->s);
            //    PD->Coupled->print_state();//getchar();
            int olda=b;
            
            int idp=inputnow;
            //           for (int idp=1; idp<=gl_keep_asis->N;idp++) {
            PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[a_st_k]]->s);
            //    PD->Coupled->print_state();//getchar();
            int s_start=PD->Coupled->s[gl_keep_asis->A[idp]-1];
            int s_start_GF=PD->Coupled->s[id_GF-1];
           // int s_start_CDL=PD->Coupled->s[id_CDL-1];
            
            PD->Coupled->s[gl_keep_asis->A[idp]-1]=48+(49-s_start);
            PD->Coupled->synchronously_update_all_Gates();
            if(gl_keep_asis->A[idp]!=id_GF) PD->Coupled->s[id_GF-1]=s_start_GF;
            
            b_new[idp]=(int)PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
            P2=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b_new[idp],1);
            
            if (( b_new[idp]!=b) && ( b_new[idp]!=olda)) {
                // if (s_start==48) printf("%s ON:",PD->BRN->MDAT->Node_name[gl_keep_asis->A[idp]]);
                // else printf("%s OFF:",PD->BRN->MDAT->Node_name[gl_keep_asis->A[idp]]);
                // printf(" %d \t-->\t %d \t %d\n",b,b_new[idp],a_st_k);
                // getchar();
                if ((PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[i]->Attractor[1]]->s[id_dead-1]==48)||(PD->Coupled->Attractor_valleys[i]->Basin_Stored[PD->Coupled->Attractor_valleys[b_new[idp]]->Attractor[1]]->s[id_dead-1]==48)) {
                    EpsSetColor_PAJEK(f, szinnev[idp]);
                    EpsDrawLine(f, PD->Coupled->Attractor_valleys[b]->x, PD->Coupled->Attractor_valleys[b]->y, PD->Coupled->Attractor_valleys[b_new[idp]]->x, PD->Coupled->Attractor_valleys[b_new[idp]]->y);
                    EpsFillCircle(f, PD->Coupled->Attractor_valleys[b]->x - 3*(PD->Coupled->Attractor_valleys[b]->x - PD->Coupled->Attractor_valleys[b_new[idp]]->x)/4., PD->Coupled->Attractor_valleys[b]->y - 3*(PD->Coupled->Attractor_valleys[b]->y - PD->Coupled->Attractor_valleys[b_new[idp]]->y)/4., 1.2);
                }
            }
            //}
        }
    }
    
    switch (yztu_inputs->N) {
        case 1:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s__Change_%s", MODEL_Directory_system,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 2:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s__Change_%s", MODEL_Directory_system,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 3:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s_%s__Change_%s", MODEL_Directory_system,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[3]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 4:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s_%s_%s__Change_%s", MODEL_Directory_system,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[3]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[4]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        default: exit(1);
            break;
    }
    
    
    close_EPSDRAW_file(f,fn,2);
}

void set_attractor_color(FILE *f,int col){
    switch (col) {
        case 1: EpsSetRgb_255(f, 0, 145, 145); break;// both kinds of senesence -- dark green
        case 2: EpsSetRgb_255(f, 148, 208, 255); break;//  mitochondrial senesence -- light blue
        case 3: EpsSetRgb_255(f, 0, 0, 255); break;//  chromosomal of senesence -- blue
            
        case 4: EpsSetRgb_255(f, 255, 128, 0); break;// G0/G1 -- orange
        case 5: EpsSetRgb_255(f, 128, 128, 0); break;// G2 --  asparagus
        case 6: EpsSetRgb_255(f, 128, 64, 0); break;// SAC -- brown
            
        case 7: EpsSetRgb_255(f, 255, 0, 0); break;// cell cycle -- red
        case 8: EpsSetRgb_255(f, 230, 230, 230); break;// other (on barrier, weird, broken) - light grey

        case 10: EpsSetRgb_255(f, 148, 82, 0); break;// sometimes migratory
        case 11: EpsSetRgb_255(f, 148, 17, 0); break;// other (on barrier, weird, broken) - light grey
        default: EpsSetRgb(f, 0, 0, 0); break;
    }
    
}

void set_celltype_color(FILE *f,int col){
    switch (col) {
        case -1: {EpsSetRgb_255(f, 0, 0, 0);return;} break;// nothing
        case 0: {EpsSetRgb_255(f, 255*0.7, 0, 0);return;} break;// naive
        case 1: {EpsSetRgb_255(f, 255*0.7, 255*0.3, 0);return;} break;//  primed
        case 2: {EpsSetRgb_255(f, 0, 255*0.5, 255*0.5);return;} break;//  Endoderm - needs renaming
        case 3: {EpsSetRgb_255(f, 0, 255*0.75, 0);return;} break;// Trophectoderm
        case 4: {EpsSetRgb_255(f, 0, 255*1, 255*1);return;} break;// Ectoderm ? doesn't really work in model
        case 5: {EpsSetRgb_255(f, 79., 143., 0);return;} break;// Neuron
        case 6: {EpsSetRgb_255(f, 128., 64., 0);return;} break;// Differentiated, weird
        default: {EpsSetRgb(f, 0, 0, 0);return;} break; // ?? - model does not include cell type
    }
}

unsigned long int get_env_ID_from_env_Index(int i,int j,int k,int l,int m){
    return (m+l*10+k*100+j*1000+i*10000);
}

void  draw_attractor_symbol(FILE *f, Attractors_on_MAP *AMAP,unsigned int index_a,double xn,double yn,
                            double aa, double bb,double sc_node){
    
    // aa space for sybmols in x? is x symbol size
    // bb is ??? space for sybmols in y? is x symbol size
    // sc_node is  ?? scaling factor
    char bla[300];
    
    EpsSetFont(f,"Times-Roman",5);
    set_attractor_color(f,AMAP->col[index_a]);
    EpsSetLinewidth(f, 0.2);
    switch (AMAP->shape[index_a]) {
        case 1:{// cell cycle
            EpsSetLinewidth(f, 0);
            EpsFillEllipse(f, xn, yn, aa*sc_node, bb*sc_node);
            EpsSetRgb(f, 0.3, 0.3,0.3);
            if(AMAP->border[index_a]==1){
                EpsSetLinewidth(f, 1.2);
                EpsDrawEllipse(f, xn, yn, aa*sc_node, bb*sc_node);
            }
            EpsSetRgb(f, 0, 0,0);
            if(AMAP->bar[index_a]) EpsDrawLine(f, xn, yn+bb*sc_node, xn, yn-bb*sc_node);
            EpsDrawCircle(f, xn, yn,AMAP->dot[index_a]/3.);
            sprintf(bla, "%ld",AMAP->al->A[index_a]);
            EpsDrawString(f, 0,xn+aa*sc_node*0.2,yn-aa*sc_node*0.2, bla);
        }
            break;
        case 2: {//G0_G1
            EpsSetLinewidth(f, 0);
            EpsFillCircle(f, xn, yn, aa*sc_node/2.);
            EpsSetRgb(f, 0.3, 0.3,0.3);
            if(AMAP->border[index_a]==1){
                EpsSetLinewidth(f, 1.2);
                EpsDrawCircle(f, xn, yn, aa*sc_node/2.);
            }
            EpsSetRgb(f, 0, 0,0);
            if(AMAP->bar[index_a])  EpsDrawLine(f, xn, yn+aa*sc_node, xn, yn-aa*sc_node/2.);
            EpsDrawCircle(f, xn, yn,AMAP->dot[index_a]/3.);
            sprintf(bla, "%ld",AMAP->al->A[index_a]);
            EpsDrawString(f, 0,xn+aa*sc_node*0.2,yn-aa*sc_node*0.2, bla);
        }
            break;
        case 3: {//G2
            EpsSetLinewidth(f, 0);
            EpsFillRectangle(f, xn-bb*sc_node,yn-bb*sc_node,xn+bb*sc_node,yn+bb*sc_node);
            EpsSetRgb(f, 0.3, 0.3,0.3);
            if(AMAP->border[index_a]==1){
                EpsSetLinewidth(f, 1.2);
                EpsDrawRectangle(f, xn-bb*sc_node,yn-bb*sc_node,xn+bb*sc_node,yn+bb*sc_node);
            }
            EpsSetRgb(f, 0, 0,0);
            if(AMAP->bar[index_a])  EpsDrawLine(f, xn, yn+bb*sc_node, xn, yn-bb*sc_node);
            EpsDrawCircle(f, xn, yn,AMAP->dot[index_a]/3.);
            sprintf(bla, "%ld",AMAP->al->A[index_a]);
            EpsDrawString(f, 0,xn+aa*sc_node*0.2,yn-aa*sc_node*0.2, bla);
        }
            break;
        case 4: {//SAC
            EpsSetLinewidth(f, 0);
            EpsFillRectangle(f, xn-aa*sc_node,yn-bb*sc_node/2.,xn+aa*sc_node,yn+bb*sc_node/2.);
            EpsSetRgb(f, 0.3, 0.3,0.3);
            if(AMAP->border[index_a]==1){
                EpsSetLinewidth(f, 1.2);
                EpsDrawRectangle(f, xn-aa*sc_node,yn-bb*sc_node/2.,xn+aa*sc_node,yn+bb*sc_node/2.);
            }
            EpsSetRgb(f, 0, 0,0);
            if(AMAP->bar[index_a]) EpsDrawLine(f, xn, yn+bb*sc_node/2., xn, yn-bb*sc_node/2.);
            EpsDrawCircle(f, xn, yn,AMAP->dot[index_a]/3.);
            sprintf(bla, "%ld",AMAP->al->A[index_a]);
            EpsDrawString(f, 0,xn+aa*sc_node*0.2,yn-aa*sc_node*0.2, bla);
        }
            break;
        case 5: {//broken cycle or ??
            EpsSetLinewidth(f, 0);
            EpsFillTriangle(f, xn-aa*sc_node,yn-bb*sc_node,xn+aa*sc_node,yn+bb*sc_node,xn,yn+bb*sc_node);
            EpsSetRgb(f, 0.3, 0.3,0.3);
            if(AMAP->border[index_a]==1){
                EpsSetLinewidth(f, 1.2);
                EpsDrawTriangle(f, xn-aa*sc_node,yn-bb*sc_node,xn+aa*sc_node,yn+bb*sc_node,xn,yn+bb*sc_node);
            }
            EpsSetRgb(f, 0, 0,0);
            if(AMAP->bar[index_a]) EpsDrawLine(f, xn, yn+bb*sc_node, xn, yn);
            EpsDrawCircle(f, xn, yn,AMAP->dot[index_a]/3.);
            sprintf(bla, "%ld",AMAP->al->A[index_a]);
            EpsDrawString(f, 0,xn+aa*sc_node*0.2,yn-aa*sc_node*0.2, bla);
        }
            break;
        case 6: {//Phase switch poised
            EpsSetLinewidth(f, 0);
            EpsFillEllipse(f, xn, yn, bb*sc_node, aa*sc_node);
            EpsSetRgb(f, 0.3, 0.3,0.3);
            if(AMAP->border[index_a]==1){
                EpsSetLinewidth(f, 1.2);
                EpsDrawEllipse(f, xn, yn, bb*sc_node, aa*sc_node);
            }
            EpsSetRgb(f, 0, 0,0);
            if(AMAP->bar[index_a])   EpsDrawLine(f, xn, yn+aa*sc_node, xn, yn-aa*sc_node);
            EpsDrawCircle(f, xn, yn,AMAP->dot[index_a]/3.);
            sprintf(bla, "%ld",AMAP->al->A[index_a]);
            EpsDrawString(f, 0,xn+aa*sc_node*0.2,yn-aa*sc_node*0.2, bla);
        }
            break;
        default:{//other, ??
            EpsSetLinewidth(f, 0);
            EpsFillEllipse(f, xn, yn, bb*sc_node/10, aa*sc_node/10);
            EpsSetRgb(f, 0.3, 0.3,0.3);
            if(AMAP->border[index_a]==1){
                EpsSetLinewidth(f, 1.2);
                EpsDrawEllipse(f, xn, yn, bb*sc_node/10, aa*sc_node/10);
            }
            EpsSetRgb(f, 0, 0,0);
            if(AMAP->bar[index_a]) EpsDrawLine(f, xn, yn+aa*sc_node/10, xn, yn-aa*sc_node/10);
            EpsDrawCircle(f, xn, yn,AMAP->dot[index_a]/3.);
            sprintf(bla, "%ld",AMAP->al->A[index_a]);
            EpsDrawString(f, 0,xn+aa*sc_node*0.2,yn-aa*sc_node*0.2, bla);
        }
            break;
    }
}

void draw_all_Attractors_in_environment(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                        longint_ARRAY *gl_keep_asis,
                                        longint_ARRAY *yztu_inputs,int id_GF,int id_GFH,
                                        int i,int j,int k,int l,int m, FILE *f,
                                        double xnow,double ynow,
                                        double aa, double bb,double sc_node){
    
    Attractors_on_MAP *AMAP;
    
    AMAP=new Attractors_on_MAP(PD,yztu_inputs,id_GF,id_GFH,i,j,k,l,m);
    
    char fn[500];
//    sprintf(fn,"mkdir %s/%s/Cytoskape/Inputs",MODEL_Directory_system,PD->Coupled->DIRname);
//    for (int kkk=2;kkk<= gl_keep_asis->N; kkk++) {
//        sprintf(fn,"%s_%s",fn, PD->Coupled->BRN->MDAT->Node_name[gl_keep_asis->A[kkk]]);
//        switch (kkk) {
//            case 1:sprintf(fn,"%s-%d",fn,i);break;
//            case 2:sprintf(fn,"%s-%d",fn,j);break;
//            case 3:sprintf(fn,"%s-%d",fn,k);break;
//            case 4:sprintf(fn,"%s-%d",fn,l);break;
//            case 5:sprintf(fn,"%s-%d",fn,m);break;
//            default:
//                break;
//        }
//    }
//    sprintf(fn,"%s\n",fn);
//    system(fn);
    
    sprintf(fn,"Inputs");
    for (int kkk=2;kkk<= gl_keep_asis->N; kkk++) {
        sprintf(fn,"%s_%s",fn, PD->Coupled->BRN->MDAT->Node_name[gl_keep_asis->A[kkk]]);
        switch (kkk) {
            case 1:sprintf(fn,"%s-%d",fn,i);break;
            case 2:sprintf(fn,"%s-%d",fn,j);break;
            case 3:sprintf(fn,"%s-%d",fn,k);break;
            case 4:sprintf(fn,"%s-%d",fn,l);break;
            case 5:sprintf(fn,"%s-%d",fn,m);break;
            default:
                break;
        }
    }
    
    // printf("i=%d j=%d k=%d l=%d m = %d    %d attractors\n",i,j,k,l,m,AMAP->al->N);
    // getchar();
    if (AMAP->al->N>0) {
        EpsSetRgb(f,0,0,0);
        EpsSetLinewidth(f, 0.25);
        EpsDrawRectangle(f,xnow-2*aa,ynow-2*aa,xnow+2*aa,ynow+2*aa);
        
        for(unsigned long int am=1;am<=AMAP->al->N;am++){
            for(unsigned int ll=1;ll<=gl_keep_asis->N;ll++){
                PD->Coupled->Attr_Transitions[ll]->node[AMAP->al->A[AMAP->at_poz[am]]].node_props=new node_attr(0,0,0);
                PD->Coupled->Attr_Transitions[ll]->node[AMAP->al->A[AMAP->at_poz[am]]].node_props->basin = get_env_ID_from_env_Index(i,j,k,l,m);
            }
            
            PD->Coupled->color[AMAP->al->A[AMAP->at_poz[am]]]=AMAP->col[AMAP->at_poz[am]];
            
            // array over positions!
            int xn= xnow-2*aa+4*aa*am/(double)AMAP->al->N;
            int yn= ynow-aa+2*aa*am/(double)AMAP->al->N;
            
            PD->Coupled->Attractor_valleys[AMAP->al->A[AMAP->at_poz[am]]]->x=xn;
            PD->Coupled->Attractor_valleys[AMAP->al->A[AMAP->at_poz[am]]]->y=yn;
            //  printf("Attractor %ld: index in list=%u, position=%lu\ncolor %d; shape %d; border %d; bar %d\n",
            //               AMAP->al->A[AMAP->at_poz[am]],AMAP->at_poz[am],am,
            //             AMAP->col[AMAP->at_poz[am]],AMAP->shape[AMAP->at_poz[am]],
            //           AMAP->border[AMAP->at_poz[am]],AMAP->bar[AMAP->at_poz[am]]);
            //  getchar();
            draw_attractor_symbol(f,AMAP,AMAP->at_poz[am],xn,yn,aa, bb, sc_node);
            if(DRAW_PDFS) PD->Coupled->export_Attractor_onto_Boolean_network__PDF(AMAP->al->A[AMAP->at_poz[am]],fn);
            //  if(DRAW_PDFS) PD->Coupled->export_Attractor_onto_Boolean_network__PDF(AMAP->al->A[AMAP->at_poz[am]]);
        }
        
        for(unsigned long int am=1;am<=AMAP->al->N;am++){
            // array over positions!
            int xn= xnow-2*aa+4*aa*am/(double)AMAP->al->N;
            int yn= ynow-aa+2*aa*am/(double)AMAP->al->N;
            EpsSetFont(f,"Times-Roman",7);
            set_celltype_color(f,AMAP->celltype[AMAP->at_poz[am]]);
            EpsDrawString(f, 0, xn-aa*sc_node*0.8, yn+aa*sc_node*0.2, celltype_name(AMAP->celltype[AMAP->at_poz[am]]));
            //  EpsDrawLine(f, xn-bb*sc_node-1, yn, xn-bb*sc_node-1, yn+bb*sc_node);
        }
    }
    delete AMAP;AMAP=NULL;
    EpsSetLinewidth(f, 0.2);
}

void make_3D_map_of_steady_states_new(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,longint_ARRAY *yztu_inputs,longint_ARRAY *gl_keep_asis, int inputnow,int pairlines){
    char fn[300], bla[100];
    FILE *f;
    double a ,alpha, beta, TW,xd,yd,aa,bb;
    
    Boolean_Dynamics_TRACKER *D;
    Phenotype_now P1,P2;
    
    D=PD->Coupled;
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    int id_GFH=PD->BRN->MDAT->get_node_ID(GF_high_name);
    
    int *b_new;
    b_new=new int[gl_keep_asis->N+1];
    
    D->Attr_Transitions = new sgraph_array[gl_keep_asis->N+1];
    D->Attr_Transitions[0]=NULL;
    for (unsigned int i=1; i<=gl_keep_asis->N; i++) {
        D->Attr_Transitions[i] = new sgraph(D->Attractor_NR);
    }
    
    a=300;
    alpha=a*0.9;beta=a*0.65;
    xd=0.4;yd=0.55;
    TW=150;
    aa=a/7.; bb=0.5*a/7.;
    double xnow = 0.0, ynow = 0.0;
    double sc_node=0.2;
    
    switch (yztu_inputs->N) {
        case 1:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s__Change_%s.eps", MODEL_Directory,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 2:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s__Change_%s.eps", MODEL_Directory,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 3:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s_%s__Change_%s.eps", MODEL_Directory,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[3]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 4:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s_%s_%s__Change_%s.eps", MODEL_Directory,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[3]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[4]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        default: exit(1);
            break;
    }
    
    
    printf("%s\n",fn);
    
    f=fopen(fn,"w");
    
    COLOR_SCHEME=GYR;
    setcolors_PAJEK();
    
    double XMAX=5*a+2*alpha+2*TW;
    double YMAX=3*a+2*beta+2*TW;
    
    //printf("XMIN=%lg\tXMAX=%lg\nYMIN=%lg\tYMAX=%lg\n",XMIN,XMAX,YMIN,YMAX);
    //getchar();
    EpsInit(f,-25,-25,XMAX+25,YMAX+25);
    
    EpsSetRgb(f,0.2,0.2,0.2);
    EpsSetLinewidth(f, 1);
    EpsSetFont(f,"Times-Roman",10);
    
    fivepointer *x,*y;
    x=new fivepointer[3];y=new fivepointer[3];
    for (int i=0; i<3; i++) {
        x[i]=new fourpointer[2];y[i]=new fourpointer[2];
        for (int j=0; j<2; j++) {
            x[i][j]=new triplepointer[2];y[i][j]=new triplepointer[2];
            for (int k=0; k<2; k++) {
                x[i][j][k]=new doublepointer[2];y[i][j][k]=new doublepointer[2];
                for (int l=0; l<2; l++) {
                    x[i][j][k][l]=new double[2];y[i][j][k][l]=new double[2];
                }
            }
        }
    }
    // shift by x coordinate
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            for (int k=0; k<2; k++) {
                for (int l=0; l<2; l++) {
                    x[i][j][k][l][0]=TW+i*a + k*xd*a       + l*(2*a+alpha);
                    x[i][j][k][l][1]=TW+i*a + k*xd*a       + l*(2*a+alpha)   + alpha/2.;
                    y[i][j][k][l][0]=TW     + k*yd*a + j*a + l* beta;
                    y[i][j][k][l][1]=TW     + k*yd*a + j*a + l* beta + 1.6*a + beta;
                }
            }
        }
    }
    
    // horizontals
    for (int j=0; j<2; j++) {
        for (int k=0; k<2; k++) {
            for (int l=0; l<2; l++) {
                for (int m=0; m<2; m++){
                    EpsDrawLine(f, x[0][j][k][l][m],y[0][j][k][l][m],x[2][j][k][l][m],y[2][j][k][l][m]);
                    if ((j==0)&&(k==0)) {
                        EpsDrawString_Centered(f, x[0][j][k][l][m],y[0][j][k][l][m]-a/5., "No GF");
                        EpsDrawString_Centered(f, x[1][j][k][l][m],y[0][j][k][l][m]-a/5., "Low GF");
                        EpsDrawString_Centered(f, x[2][j][k][l][m],y[0][j][k][l][m]-a/5., "High GF");
                    }
                }
            }
        }
    }
    // verticals
    
    if (yztu_inputs->N>1){
        for (int i=0; i<3; i++) {
            for (int j=0; j<2; j++) {
                for (int l=0; l<2; l++) {
                    for (int m=0; m<2; m++){
                        EpsDrawLine(f, x[i][j][0][l][m],y[i][j][0][l][m],x[i][j][1][l][m],y[i][j][1][l][m]);
                        if ((i==0)&&(j==0)) {
                            EpsDrawString_Centered(f,  (x[0][0][0][l][m]+x[0][0][1][l][m])/2.-a/6.,
                                                   (y[0][0][0][l][m]+y[0][0][1][l][m])/2.,PD->BRN->MDAT->Node_name[yztu_inputs->A[2]]);
                        }
                    }
                }
            }
        }
    }
    // diagonals
    for (int i=0; i<3; i++) {
        for (int k=0; k<2; k++) {
            for (int l=0; l<2; l++) {
                for (int m=0; m<2; m++){
                    EpsDrawLine(f, x[i][0][k][l][m],y[i][0][k][l][m],x[i][1][k][l][m],y[i][1][k][l][m]);
                    if ((i==0)&&(k==0)) {
                        EpsDrawString_Centered(f,  (x[0][0][0][l][m]+x[0][1][0][l][m])/2.-a/6.,
                                               (y[0][0][0][l][m]+y[0][1][0][l][m])/2., PD->BRN->MDAT->Node_name[yztu_inputs->A[1]]);
                    }
                }
            }
        }
    }
    
    
    
    if (yztu_inputs->N==3){
        EpsSetLinewidth(f, 1.2);
        EpsDrawLine(f,  x[0][1][1][0][1]-6*aa,     y[0][1][1][0][1]+2*aa,
                    x[0][1][1][0][1]+2*aa+2*a, y[0][1][1][0][1]+2*aa);
        sprintf(bla, "No %s",PD->BRN->MDAT->Node_name[yztu_inputs->A[3]]);
        EpsDrawString_Centered(f, (x[0][1][1][0][1]-6*aa+x[0][1][1][0][1]+2*aa+2*a)/2.,
                               (y[0][1][1][0][1]+2*aa+y[0][1][1][0][1]+2*aa)/2.+aa, bla);
        
        EpsDrawLine(f,  x[0][1][1][1][1]-6*aa,     y[0][1][1][1][1]+2*aa,
                    x[0][1][1][1][1]+2*aa+2*a, y[0][1][1][1][1]+2*aa);
        EpsDrawString_Centered(f, (x[0][1][1][1][1]-6*aa+x[0][1][1][1][1]+2*aa+2*a)/2.,
                               (y[0][1][1][1][1]+2*aa+y[0][1][1][1][1]+2*aa)/2.+aa, PD->BRN->MDAT->Node_name[yztu_inputs->A[3]]);
    }
    
    if (yztu_inputs->N==4) {
        EpsDrawLine(f, x[0][0][0][0][0]-3*aa, y[0][0][0][0][0]-2*aa,
                    x[0][0][0][0][0]-3*aa, y[0][0][0][0][0]+4*aa);
        sprintf(bla, "No %s",PD->BRN->MDAT->Node_name[yztu_inputs->A[4]]);
        EpsDrawString_Centered(f, (x[0][0][0][0][0]+x[0][0][1][0][0]-6*aa)/2.-2*aa,
                               (y[0][0][0][0][0]+y[0][0][1][0][0]+2*aa)/2., bla);
        
        // EpsDrawCircle(f, x[0][0][0][0][1], y[0][0][0][0][1], 20);
        EpsDrawLine(f, x[0][0][0][0][1]-3*aa, y[0][0][0][0][1]-2*aa,
                    x[0][0][0][0][1]-3*aa, y[0][0][0][0][1]+4*aa);
        EpsDrawString_Centered(f,   (x[0][0][0][0][1]+x[0][0][1][0][1]-6*aa)/2.-2*aa,
                               (y[0][0][0][0][1]+y[0][0][1][0][1]+2*aa)/2., PD->BRN->MDAT->Node_name[yztu_inputs->A[4]]);
    }
    
    
    EpsSetFont(f,"Times-Roman",5);
    EpsSetLinewidth(f, 0.2);
    
    if(PD->Coupled->color==NULL) {
        PD->Coupled->color=new unsigned int[PD->Coupled->Attractor_NR+1];
        for (int i=0; i<=PD->Coupled->Attractor_NR; i++) PD->Coupled->color[i]=0;
    }
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            if(yztu_inputs->N>1){
                for (int k=0; k<2; k++) {
                    if(yztu_inputs->N>2){
                        for (int l=0; l<2; l++) {
                            if(yztu_inputs->N>3) for (int m=0; m<2; m++){
                                xnow=x[i][j][k][l][m];
                                ynow=y[i][j][k][l][m];
                                draw_all_Attractors_in_environment(PD,gl_keep_asis, yztu_inputs, id_GF, id_GFH,i,j,k,l,m,f,xnow,ynow,aa,bb,sc_node);
                            }
                            else {
                                xnow=x[i][j][k][l][0];
                                ynow=y[i][j][k][l][0];
                                draw_all_Attractors_in_environment(PD, gl_keep_asis,yztu_inputs, id_GF, id_GFH,i,j,k,l,0,f,xnow,ynow,aa,bb,sc_node);
                            }
                        }
                    }
                    else {   xnow=x[i][j][k][0][0];
                        ynow=y[i][j][k][0][0];
                        draw_all_Attractors_in_environment(PD, gl_keep_asis,yztu_inputs, id_GF, id_GFH,i,j,k,0,0,f,xnow,ynow,aa,bb,sc_node);
                    }
                }
            }
            else {   xnow=x[i][j][0][0][0];
                ynow=y[i][j][0][0][0];
                draw_all_Attractors_in_environment(PD, gl_keep_asis,yztu_inputs, id_GF, id_GFH,i,j,0,0,0,f,xnow,ynow,aa,bb,sc_node);
            }
        }
    }
    
    for(int b=1;b<=PD->Coupled->Attractor_NR;b++){
        P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        for (int a_st_k=1; a_st_k<=D->Attractor_valleys[b]->N_a; a_st_k++) {
            int olda=b;
            int idp=inputnow;
            PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[a_st_k]]->s);
            //    PD->Coupled->print_state();//getchar();
            int s_start=PD->Coupled->s[gl_keep_asis->A[idp]-1];
            int s_start_GF=PD->Coupled->s[id_GF-1];
            
            PD->Coupled->s[gl_keep_asis->A[idp]-1]=48+(49-s_start);
            PD->Coupled->synchronously_update_all_Gates();
            if(gl_keep_asis->A[idp]!=id_GF) PD->Coupled->s[id_GF-1]=s_start_GF;
            
            b_new[idp]=(int)PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
            
            P2=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b_new[idp],1);
            
            //    printf("Transition when %s is flipped from %d to %d: Arrt %ld -> %ld\n",
            //          PD->BRN->MDAT->Node_name[gl_keep_asis->A[idp]],
            //          s_start,PD->Coupled->s[gl_keep_asis->A[idp]-1],
            //          b,b_new[idp]);
            
            if (( b_new[idp]!=b) && ( b_new[idp]!=olda)) {
                if(D->Attractor_valleys[b_new[idp]]->N_a==1){
                    olda=b_new[idp];
                    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[olda]->Basin_Stored[PD->Coupled->Attractor_valleys[olda]->Attractor[1]]->s);
                    int s_start=PD->Coupled->s[gl_keep_asis->A[idp]-1];
                    int s_start_GF=PD->Coupled->s[id_GF-1];
                    
                    PD->Coupled->s[gl_keep_asis->A[idp]-1]=48+(49-s_start);
                    PD->Coupled->synchronously_update_all_Gates();
                    if(gl_keep_asis->A[idp]!=id_GF) PD->Coupled->s[id_GF-1]=s_start_GF;
                    int newback=(int)PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
                    if ((pairlines==1)||(newback!=b)) {
                        set_attractor_color(f,PD->Coupled->color[b]);
                        EpsDrawLine(f,  PD->Coupled->Attractor_valleys[b]->x,
                                    PD->Coupled->Attractor_valleys[b]->y,
                                    PD->Coupled->Attractor_valleys[b_new[idp]]->x,
                                    PD->Coupled->Attractor_valleys[b_new[idp]]->y);
                        EpsFillCircle(f,
                                      PD->Coupled->Attractor_valleys[b]->x - 3*(PD->Coupled->Attractor_valleys[b]->x - PD->Coupled->Attractor_valleys[b_new[idp]]->x)/4.,
                                      PD->Coupled->Attractor_valleys[b]->y - 3*(PD->Coupled->Attractor_valleys[b]->y - PD->Coupled->Attractor_valleys[b_new[idp]]->y)/4., 1.2);
                    }
                }
                else { set_attractor_color(f,PD->Coupled->color[b]);
                    //EpsSetColor_PAJEK(f, szinnev[0]);
                    EpsDrawLine(f,  PD->Coupled->Attractor_valleys[b]->x,
                                PD->Coupled->Attractor_valleys[b]->y,
                                PD->Coupled->Attractor_valleys[b_new[idp]]->x,
                                PD->Coupled->Attractor_valleys[b_new[idp]]->y);
                    EpsFillCircle(f,
                                  PD->Coupled->Attractor_valleys[b]->x - 3.5*(PD->Coupled->Attractor_valleys[b]->x - PD->Coupled->Attractor_valleys[b_new[idp]]->x)/4.,
                                  PD->Coupled->Attractor_valleys[b]->y - 3.5*(PD->Coupled->Attractor_valleys[b]->y - PD->Coupled->Attractor_valleys[b_new[idp]]->y)/4., 1.2);
                }
            }
        }
    }
    
    switch (yztu_inputs->N) {
        case 1:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s__Change_%s", MODEL_Directory_system,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 2:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s__Change_%s", MODEL_Directory_system,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 3:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s_%s__Change_%s", MODEL_Directory_system,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[3]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        case 4:
            sprintf(fn,"%s/%s/%s_Steady_state_MAP_GF_%s_%s_%s_%s__Change_%s", MODEL_Directory_system,
                    PD->Coupled->DIRname,PD->BRN->name,
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[1]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[2]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[3]],
                    PD->BRN->MDAT->Node_name[yztu_inputs->A[4]],PD->BRN->MDAT->Node_name[gl_keep_asis->A[inputnow]]);
            break;
        default: exit(1);
            break;
    }
    
    close_EPSDRAW_file(f,fn,2);
}

void write_link_for_All_environmental_changes_from_ONE_attractor(longint_ARRAY *gl_keep_asis, int attr_id,
                                                                 MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Boolean_Dynamics_TRACKER *D;
    Phenotype_now P1,P2;
    
    D=PD->Coupled;
    
    int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    //int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    int *b_new;
    b_new=new int[gl_keep_asis->N+1];
    
    unsigned int b=attr_id;
    for (int a_st_k=1; a_st_k<=D->Attractor_valleys[b]->N_a; a_st_k++) {
        PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[a_st_k]]->s);
        //    PD->Coupled->print_state();//getchar();
        P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        
        for (int idp=1; idp<=gl_keep_asis->N;idp++) {
            PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[a_st_k]]->s);
            //    PD->Coupled->print_state();//getchar();
            int s_start=PD->Coupled->s[gl_keep_asis->A[idp]-1];
            int s_start_GF=PD->Coupled->s[id_GF-1];
            //int s_start_CDL=PD->Coupled->s[id_CDL-1];
            
            PD->Coupled->s[gl_keep_asis->A[idp]-1]=48+(49-s_start);
            PD->Coupled->synchronously_update_all_Gates();
            if(gl_keep_asis->A[idp]!=id_GF) PD->Coupled->s[id_GF-1]=s_start_GF;
            
            b_new[idp]=(int)PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
            P2=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b_new[idp],1);
            
            //if (( b_new[idp]!=b) && ( b_new[idp]!=olda)) {
            //    printf(" %d \t-->\t %d \t %d\n",b,b_new[idp],a_st_k);
            //}
        }
    }
}

void run_OnePulse_analysis(longint_ARRAY *gl_keep_asis, MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Boolean_Dynamics_TRACKER *D;
    Phenotype_now P1,P2;
    
    D=PD->Coupled;
    
    int MAX_PULSE=80;
    unsigned long int oldan=D->Attractor_NR;
    
    for(unsigned int ll = 1; ll<=gl_keep_asis->N ; ll++){
        unsigned int id_pulse=(unsigned int)gl_keep_asis->A[ll];
        unsigned int id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
       // unsigned int id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
        unsigned int *b_new,bfinal;
        b_new=new unsigned int[MAX_PULSE+1];
        
        for(unsigned int b=1;b<=D->Attractor_NR;b++) {
            for (int a_st_k=1; a_st_k<=D->Attractor_valleys[b]->N_a; a_st_k++) {
                PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[a_st_k]]->s);
                //    PD->Coupled->print_state();//getchar();
                P1=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
                
                unsigned int s_start=PD->Coupled->s[id_pulse-1];
                unsigned int s_start_GF=PD->Coupled->s[id_GF-1];
               // if(id_CDL>0) int s_start_CDL=PD->Coupled->s[id_CDL-1];
              //  else int s_start_CDL=0;
                
                unsigned int olda=b;
                unsigned int irrev=0;
                bfinal=0;
                for (unsigned int pulse_length=1; pulse_length<=MAX_PULSE; pulse_length++) {
                    PD->Coupled->s[id_pulse-1]=48+(49-s_start);
                    PD->Coupled->synchronously_update_all_Gates();
                    
                    if(id_pulse != id_GF) PD->Coupled->s[id_GF-1]=s_start_GF;
                    if(pulse_length==MAX_PULSE){
                        bfinal=(int)PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
                    }
                    PD->Coupled->s[id_pulse-1]=s_start;
                    
                    // b_new is AFTER the pulse is over and env. returned to original!
                    b_new[pulse_length]=(int)PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
                    P2=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b_new[pulse_length],1);
                    if (( b_new[pulse_length]!=b) && ( b_new[pulse_length]!=olda)) {
                        if (irrev==0) irrev=pulse_length;
                        if ((b_new[pulse_length]<=PD->Coupled->Attr_Transitions[ll]->N)&&(b <= PD->Coupled->Attr_Transitions[ll]->N))
                            if((b<=oldan)&&(b_new[pulse_length]<=oldan)) PD->Coupled->Attr_Transitions[ll]->add_link(b, b_new[pulse_length], pulse_length,a_st_k,id_pulse);
                    }
                    olda=b_new[pulse_length];
                }
                if (irrev==0){
                    if((b<=oldan)&&(bfinal<=oldan))
                        PD->Coupled->Attr_Transitions[ll]->add_link(b, bfinal, 0.,a_st_k,id_pulse);
                }
                else   {
                    if((b<=oldan)&&(bfinal<=oldan))
                        PD->Coupled->Attr_Transitions[ll]->add_link(b, bfinal, -1.,a_st_k,id_pulse);
                }
            }
        }
    }
}



void draw_mitogen_dose_Reprogramming(int a_st, int a_st_k, MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, double p_error, double p_GF, int idgf, FILE *f, double y0, int W) {
    int T_MAX = 450, TW = 5, NameBuff = 30;
    double r;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    int YAM = PD->BRN->MDAT->get_node_ID("Yamanaka");
    // int Twoi = PD->BRN->MDAT->get_node_ID("2i");
    
    for (int t = 1; t <= T_MAX; t++) {
        if (t <= 30) { PD->Coupled->s[YAM - 1] = 48;}
        else { PD->Coupled->s[YAM - 1] = 49;}
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        r = rand_real(0, 1);
        if (r <= p_GF) PD->Coupled->s[idgf - 1] = 49;
        else PD->Coupled->s[idgf - 1] = 48;
        draw_timestep_MOD(PD, t, f, y0, W, T_MAX, TW, NameBuff);
    }
}

void draw_mitogen_dose_series_Reprogramming(int a_st, int a_st_k, MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, double p_error) {
    int W = 8, YTOP, XTOP;
    FILE *f;
    char fn[400];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    PD->Coupled->print_state();//getchar();
    
    int idgf = PD->BRN->MDAT->get_node_ID(GF_high_name);
    
    //int idtrack=PD->BRN->MDAT->get_node_ID(markernode);
    double pGF = 0.5;
    sprintf(fn,"%s/%s/_EXP/Mitogen_dose_experiment.eps_%lg.eps",MODEL_Directory,PD->BRN->path,pGF);
    f = fopen(fn, "w");
    
    //int DMAX = 25;
    int N_lines = PD->BRN->N + 10;
    YTOP = (W*N_lines + 3 * W) + 5;
    XTOP = 2 * 800 + 5;
    
    EpsInit(f, -5, -5, XTOP, YTOP);
    EpsSetFont(f, "Times-Roman", W);
    
    draw_mitogen_dose_Reprogramming(a_st, a_st_k, PD, p_error, pGF, idgf, f, YTOP, W);
    sprintf(fn,"mkdir %s/%s/_EXP/Mitogen_dose_experiment.eps_%lg",MODEL_Directory_system,PD->BRN->path,pGF);
    close_EPSDRAW_file(f,fn,2);
    
}

unsigned long int attractor_after_pulse(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                        unsigned long int b,unsigned int st_a,
                                        unsigned long int id_pulse,unsigned long int id_GF, unsigned long int id_CDL, unsigned int pulse){
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[st_a]]->s);
    
    unsigned int s_start=PD->Coupled->s[id_pulse-1];
    unsigned int s_start_GF=PD->Coupled->s[id_GF-1];
    unsigned int s_start_CDL=PD->Coupled->s[id_CDL-1];
    
    unsigned int b_back=0;
    for (unsigned int pulse_length=1; pulse_length<=pulse; pulse_length++) {
        PD->Coupled->s[id_pulse-1]=48+(49-s_start);
        PD->Coupled->synchronously_update_all_Gates();
        if(id_pulse != id_GF) PD->Coupled->s[id_GF-1]=s_start_GF;
        if(id_pulse != id_CDL) PD->Coupled->s[id_CDL-1]=s_start_CDL;
    }
    PD->Coupled->s[id_pulse-1]=s_start;
    
    b_back=(int)PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
    return(b_back);
}


void run_transition_from_cycle(unsigned long int a_ID,const char nodename[],MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A){
    FILE *f;
    Attractors_on_MAP *AMAP;
    double X_BOX=20.,Y_BOX=20.;
    Boolean_Dynamics_TRACKER *D;
    int ll_now;
    char bla[300];
    
    unsigned long int pulse_ID = PD_A->BRN->MDAT->get_node_ID(nodename);
    ll_now=0;
    for(unsigned int ll=1;ll<=PD_A->BRN->active_inputs->N;ll++){
        if (pulse_ID == PD_A->BRN->active_inputs->A[ll]) {
            ll_now=ll;
        }
    }
    if (ll_now==0) {
        printf("Pulse of %s does now work, %s is not one of the pre-set imput nodes!\n",nodename,nodename);
        return;
    }
    unsigned long int GF_id = PD_A->BRN->MDAT->get_node_ID(GF_low_name);
    unsigned long int id_CDL= PD_A->BRN->MDAT->get_node_ID(CD_low_name);
    D=PD_A->Coupled;
    
    double pmax=2*(D->Attractor_valleys[a_ID]->N_a+1);
   /* if(D->Attr_Transitions[ll_now]!=NULL){
        snode::slinklist::slink *sc2;
        sc2=D->Attr_Transitions[ll_now]->node[a_ID].Li.first_link_pt;
        while(sc2!=NULL){
            printf("\tlink from attractor %ld to %llu :ID = %d  label = %d (pulseid= %ld)  w= %lg\n",
                   a_ID,sc2->sneighbor,
                   sc2->link_props->ID,
                   sc2->link_props->distinc_link_labels,pulse_ID,
                   sc2->link_props->lw);//getchar();
            if((sc2->link_props->distinc_link_labels == pulse_ID)&&
               (sc2->link_props->lw > pmax)) pmax=(int)sc2->link_props->lw;
            sc2=sc2->next(sc2);
        }
    }
    else {printf("PD_A->Coupled->Attr_Transitions[%d] does not exist!!!\n",ll_now);getchar();}
    */
    char fn[400];
    sprintf(fn,"%s/%s/_EXP/%s_Response_of_ATTR-%ld.eps",MODEL_Directory,PD_A->BRN->path,nodename,a_ID);
    printf("%s\npmax = %lf",fn,pmax);
    
    int draw=1;
    if(pmax==0) {draw=0;}
    f=fopen(fn,"w");
    EpsInit(f,-5,-5,X_BOX*(pmax+2),Y_BOX*(D->Attractor_valleys[a_ID]->N_a+2));
    EpsSetFont(f,"Times-Roman",5);
    
    double y;
    for (unsigned int pp=1; pp<=pmax+5; pp++) {
        sprintf(bla, "%d",pp);
        EpsDrawString(f,0, pp*X_BOX,(D->Attractor_valleys[a_ID]->N_a+1)*Y_BOX, bla);
        //   if(draw) run_OneNode_pulse_experiment_MOD(nodename,(int)(a_ID),1,PD_A,0.000000001,pp);
        for(int i=1;i<=D->Attractor_valleys[a_ID]->N_a;i++){
            unsigned long int bback=attractor_after_pulse(PD_A,a_ID,i,pulse_ID,GF_id,id_CDL,pp);
            AMAP = new Attractors_on_MAP(PD_A,bback);
            y=D->Attractor_valleys[a_ID]->N_a - (i-1);
            draw_attractor_symbol(f, AMAP,1,pp*X_BOX,y*Y_BOX,30,30,0.2);
            delete AMAP; AMAP=NULL;
            sprintf(bla, "%d",i);
            EpsDrawString(f,0, 0 ,y*Y_BOX, bla);
        }
    }
    if(draw) run_OneNode_pulse_experiment_MOD(nodename,(int)(a_ID),1,PD_A,0.000000001,pmax);
    
    sprintf(fn,"%s/%s/_EXP/%s_Response_of_ATTR-%ld",MODEL_Directory_system,PD_A->BRN->path,nodename,a_ID);
    close_EPSDRAW_file(f,fn,2);
}


void build_Model_and_sample_state_space(const char modelname[],const char modules_name[]){
    Boolean_RegNet *A;
    //Boolean_Dynamics_TRACKER *D;
    unsigned int N_trace,N_RND,N_MAX;
    double p_error;
    MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD;
    longint_ARRAY *gl_keep_asis;
    
    N_trace=5;N_RND=100;
    
    p_error=0.02;N_MAX=1000;
    
    double x_int[2];x_int[0]=8;x_int[1]=10000;
    
    A=new Boolean_RegNet(modelname);
    
    A->load_Module_IDs_and_Drawing_Order();
    //    int m=A->load_Module_Drawing_order(modules_name);
    //    A->load_Module_IDs(modules_name,m);
    
    gl_keep_asis=get_input_node_names(A);
    
    PD=new MODULAR_Boolean_Dynamics_PAIRED_TRACKER(A,N_MAX);
    PD->BRN->export_Boolean_RegNet_with_names_to_Cytoskape();
    
    PD->read_Node_Cytoskape_positions();
    
    // PD->Unconnected->record_timetraces_from_collection_of_random_inputs_per_inputBOX(N_trace,N_RND,p_error,gl_keep_asis);
    // PD->Unconnected->Export_ALL_Attractors_onto_Boolean_Network_ALIVE();
    // Export_Attractor_Life_and_Death_table(PD,gl_keep_asis,0);
    
    PD->Coupled->record_timetraces_from_collection_of_random_inputs_per_inputBOX(N_trace,N_RND,p_error,gl_keep_asis);
    
    
    for(int i=1;i<=PD->Coupled->Attractor_NR;i++){
        write_link_for_All_environmental_changes_from_ONE_attractor(gl_keep_asis,i,PD);
    }
    PD->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis->N)),p_error);
    Export_Attractor_Life_and_Death_table(PD,gl_keep_asis,1);
    
    set_up_module_MetaData(PD);
    //  work_with_individual_Modules(1,PD,0.0000001); // R-SW
    //  work_with_individual_Modules(2,PD,0.0000001); // PH-SW
    //   work_with_individual_Modules(9,PD,0.0000001); // ORC
    //   work_with_individual_Modules(4,PD_A,0.0000001); // Apoptosis
    work_with_individual_Modules(3,PD,0.0000001); // AKT
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
    
    
    run_OneNode_pulse_experiment_MOD("CellDensity_high",27,1,PD,0.000000001,100);
    run_OneNode_pulse_experiment_MOD("CellDensity_high",28,1,PD,0.000000001,100);
    run_OneNode_pulse_experiment_MOD("CellDensity_high",25,1,PD,0.000000001,100);
    
    run_OneNode_pulse_experiment_MOD("ECM",27,1,PD,0.000000001,100);
    run_OneNode_pulse_experiment_MOD("ECM",28,1,PD,0.000000001,100);
    run_OneNode_pulse_experiment_MOD("ECM",25,1,PD,0.000000001,100);
    
    //  run_OneNode_pulse_experiment_MOD("ECM",26,1,PD,0.000000001,100);
    
    //   run_OneNode_pulse_experiment_MOD("NDD",83,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("EGF_High",65,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("DLL1",65,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("NDD",22,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("NDD",48,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("NDD",67,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("NDD",66,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("NDD",48,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("EGF_High",57,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("NGF",40,1,PD,0.000000001,100);
    //    run_OneNode_pulse_experiment_MOD("DLL1",57,1,PD,0.000000001,100);
    
    //   exit(1);
    //
    // PD->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis->N)),p_error);
    
    //  PD->Coupled->Export_ALL_Attractors_onto_Boolean_Network(0);
    
    //PD->Coupled->Export_ALL_Attractors_onto_Boolean_Network_ALIVE(0);
    
    //  PD->Coupled->Export_ALL_Attractors_onto_Boolean_Network_ALIVE_or_cycle(0);
    
    
    
    longint_ARRAY *yztu_inputs;
    yztu_inputs=get_input_nodes_for_5input_drawing(A);
    for (int i=1; i<=gl_keep_asis->N; i++) {
        make_3D_map_of_steady_states_new(PD,yztu_inputs,gl_keep_asis,i,0);
    }
    run_OnePulse_analysis(gl_keep_asis,PD);
    
    longint_ARRAY *l_sim;
    l_sim= Find_Similar_Attractors__Life_and_Death(PD,gl_keep_asis,1);
    
    //   run_transition_from_cycle(51,"NGF",PD,gl_keep_asis);
    
    exit(1);
    
    //PD->Coupled->export_Attractor_onto_Boolean_network__PDF(31);
    
    
    run_TwoNode_pulse_experiment_MOD("Gamma","GF_High",1,1,PD,0.000000001,80);
    
    //    run_OneNode_pulse_experiment_MOD("GF_High",36,1,PD,0.000000001,40);
    //    run_OneNode_pulse_experiment_MOD("GF_High",38,1,PD,0.000000001,40);
    //   run_OneNode_pulse_experiment_MOD("GF_High",39,1,PD,0.000000001,40);
    
    //run_OneNode_pulse_experiment_MOD("Gamma",47,1,PD,0.000000001,50);
    //run_OneNode_pulse_experiment_MOD("Gamma",56,1,PD,0.000000001,50);
    //run_OneNode_pulse_experiment_MOD("Gamma",58,1,PD,0.000000001,50);
    //run_OneNode_pulse_experiment_MOD("Gamma",64,1,PD,0.000000001,50);
    //run_OneNode_pulse_experiment_MOD("Gamma",66,1,PD,0.000000001,50);
    //run_OneNode_pulse_experiment_MOD("Gamma",54,1,PD,0.000000001,50);
    
    // run_OneNode_pulse_experiment_MOD("UV",59,1,PD,0.000000001,3);
    //  run_OneNode_pulse_experiment_MOD("GF_High",56,1,PD,0.000000001,40);
    
    //  exit(1);
    
    longint_ARRAY *l;
    l=get_Healthy_CellCycle_attractors(PD);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){ // for every healthy cell cycle it finds
            for (int i=1;i<=2; i++) { // for pulse lengths from 1 to 2
                run_OneNode_pulse_experiment_MOD("Gamma",(int)(l->A[j]),1,PD,0.000000001,i);
                run_OneNode_pulse_experiment_MOD("ROS_ext",(int)(l->A[j]),1,PD,0.000000001,i);
            }
            run_OneNode_pulse_experiment_MOD("Gamma",(int)(l->A[j]),1,PD,0.000000001,50);
            run_OneNode_pulse_experiment_MOD("ROS_ext",(int)(l->A[j]),1,PD,0.000000001,50);
        }
    }
    delete l;l=NULL;
    //  exit(1);
    
    l=get_Healthy_G0_attractors(PD);
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            for (int i=24;i<=24; i++) {
                //           run_OneNode_pulse_experiment_MOD("GF_High",(int)(l->A[j]),1,PD,0.000000001,i);
            }
            for (int i=1;i<=2; i++) {
                run_OneNode_pulse_experiment_MOD("Gamma",(int)(l->A[j]),1,PD,0.000000001,i);
            }
            run_OneNode_pulse_experiment_MOD("GF_High",(int)(l->A[j]),1,PD,0.000000001,50);
            run_OneNode_pulse_experiment_MOD("Gamma",(int)(l->A[j]),1,PD,0.000000001,50);
            run_OneNode_pulse_experiment_MOD("GF",(int)(l->A[j]),1,PD,0.000000001,50);
            run_OneNode_pulse_experiment_MOD("ROS_ext",(int)(l->A[j]),1,PD,0.000000001,50);
        }
    }
    
    
    PD->Coupled->export_state_of_tracking(N_trace,(int)(N_RND*pow(2,gl_keep_asis->N)),p_error);
    
    Export_Attractor_Life_and_Death_table(PD,gl_keep_asis,1);
    
    printf("DONE FOR REAL!\n");
    exit(1);
    
}


#endif /* Analyze_Boolean_Model_h */
