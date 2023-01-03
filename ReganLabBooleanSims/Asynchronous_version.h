//
//  Asynchronous_version.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 10/2/18.
//  Copyright © 2018 Regan_Group. All rights reserved.
//

#ifndef Asynchronous_version_h
#define Asynchronous_version_h

#define forced_onoff 0  // input nodes, random prob p of ON  and 1 - p for OFF
#define forced_on_reg 1 // knockdown , random prob p of OFF  and 1 - p for whatever the Boolean rules dictate
#define forced_off_reg 2 // overexpression , random prob p of ON  and 1 - p for whatever the Boolean rules dictate
#define forced_on 3  // input nodes, ON
#define forced_off 4  // input nodes, OFF

struct pulse {
    char nodename[100];
    unsigned int t_st,t_end;
    unsigned short int pulsetype;
    double p_force;
};

typedef pulse * pulselist;

void draw_timestep_MOD_Assync(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int t,
                              doublepointer *A,
                              FILE *f,double y0,int W,int T_MAX){
    
    long int i;
    int poz=0;
    int TW=3;
    int NameBuff=55;
    
    for(int m=1;m<=PD->BRN->Module_Order->N;m++){
        if((PD->BRN->Module_Order->A[m]>0)&&(PD->BRN->Module_Order->A[m]<=PD->BRN->MDAT->Module_NR)){
            // printf("Module %d\t",m);
            // printf("Going for it %ld \n",PD->BRN->Module_Order->A[m]);
            // printf("\tNodes: %ld \n",PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N);
            for(int nn=1;nn<=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N;nn++){
                i=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->A[nn];
                //   printf("\t nn=%d node = %ld %s\n",nn,i,PD->BRN->MDAT->Node_name[i]);
                EpsSetRgb(f,A[i][t],0.5*A[i][t],1-A[i][t]);
//                if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
                EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                EpsSetRgb(f,0,0,0);
                EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
                poz++;
            }
        }
        else{
            i=PD->BRN->Module_Order->B[m];
            EpsSetRgb(f,A[i][t],0.5*A[i][t],1.-A[i][t]);
         //   if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
            EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
            EpsSetRgb(f,0,0,0);
            EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
            poz++;
        }
        poz++;
    }
}

void draw_timestep_MOD_Assync(Boolean_Dynamics_TRACKER *D,int t,
                              doublepointer *A,
                              FILE *f,double y0,int W,int T_MAX){
    
    long int i;
    int poz=0;
    int TW=3;
    int NameBuff=55;
    
    for(int m=1;m<=D->BRN->Module_Order->N;m++){
        if((D->BRN->Module_Order->A[m]>0)&&(D->BRN->Module_Order->A[m]<=D->BRN->MDAT->Module_NR)){
            // printf("Module %d\t",m);
            // printf("Going for it %ld \n",PD->BRN->Module_Order->A[m]);
            // printf("\tNodes: %ld \n",PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N);
            for(int nn=1;nn<=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[m]]->N;nn++){
                i=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[m]]->A[nn];
                //   printf("\t nn=%d node = %ld %s\n",nn,i,PD->BRN->MDAT->Node_name[i]);
                EpsSetRgb(f,A[i][t],0.5*A[i][t],1-A[i][t]);
                //                if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
                EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                EpsSetRgb(f,0,0,0);
                EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
                poz++;
            }
        }
        else{
            i=D->BRN->Module_Order->B[m];
            EpsSetRgb(f,A[i][t],0.5*A[i][t],1.-A[i][t]);
            //   if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
            EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
            EpsSetRgb(f,0,0,0);
            EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
            poz++;
        }
        poz++;
    }
}

void draw_timestep_MOD_Assync(Boolean_Dynamics *D,int t,
                              doublepointer *A,
                              FILE *f,double y0,int W,int T_MAX){
    
    long int i;
    int poz=0;
    int TW=3;
    int NameBuff=55;
    
    for(int m=1;m<=D->BRN->Module_Order->N;m++){
        if((D->BRN->Module_Order->A[m]>0)&&(D->BRN->Module_Order->A[m]<=D->BRN->MDAT->Module_NR)){
            // printf("Module %d\t",m);
            // printf("Going for it %ld \n",PD->BRN->Module_Order->A[m]);
            // printf("\tNodes: %ld \n",PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N);
            for(int nn=1;nn<=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[m]]->N;nn++){
                i=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[m]]->A[nn];
                //   printf("\t nn=%d node = %ld %s\n",nn,i,PD->BRN->MDAT->Node_name[i]);
                EpsSetRgb(f,A[i][t],0.5*A[i][t],1-A[i][t]);
                //                if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
                EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
                EpsSetRgb(f,0,0,0);
                EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
                poz++;
            }
        }
        else{
            i=D->BRN->Module_Order->B[m];
            EpsSetRgb(f,A[i][t],0.5*A[i][t],1.-A[i][t]);
            //   if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1); else EpsSetRgb(f,1,0.5,0);
            EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
            EpsSetRgb(f,0,0,0);
            EpsDrawString(f,0,0,y0-poz*W-W,D->BRN->MDAT->Node_name[i]);
            poz++;
        }
        poz++;
    }
}

void timecourse_with_asynchronous_update_Overlay(Boolean_Dynamics_TRACKER *D,
                                                 unsigned int a_st,unsigned int a_st_k,
                                                 unsigned int STATS){
    
    int W=8,YTOP,XTOP,T_MAX=300;
    FILE *f;
    char fn[600];
  //  unsigned int *pl_ID,*pl_st;
    doublepointer *Act_av;
    
    D->set_state(D->Attractor_valleys[a_st]->Basin_Stored[D->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int N_lines=D->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*650+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,D->BRN->path);
    system(fn);
    
    sprintf(fn,"mkdir %s/%s/_EXP/Assync_%d_TimeCourses\n",MODEL_Directory_system,D->BRN->path,ASSYNC_BIAS);
    system(fn);
    
    Act_av = new doublepointer [D->BRN->N+1];
    for (int i=0; i<=D->BRN->N; i++){
        Act_av[i]=new double[T_MAX+1];
        for (int j=0; j<=T_MAX; j++){
            Act_av[i][j]=0;
        }
    }
    
    sprintf(fn,"%s/%s/_EXP/Assync_%d_TimeCourses/Attr_%u-%u__S-%d.eps",
            MODEL_Directory,D->BRN->path,ASSYNC_BIAS,a_st,a_st_k,STATS);
   // printf("%s\n",fn);
    f=fopen(fn,"w");
    
    EpsInit(f,-5,-5,XTOP,YTOP);
    EpsSetFont(f,"Times-Roman",W);
    
    for (unsigned int s=1;s<=STATS; s++) {
        D->set_state(D->Attractor_valleys[a_st]->Basin_Stored[D->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
        
        for (int t=1;t<=T_MAX; t++) {
            if(ASSYNC_BIAS == 0)    D->Asynchronously_update_all_Gates();
            else D->Asynchronously_update_all_Gates(ASSYNC_BIAS,NULL,0);
            
            for (int i=1; i<=D->BRN->N; i++){
                Act_av[i][t]+=D->s[i-1]-48;
            }
        }
    }
    
    for (int i=1; i<=D->BRN->N; i++){
        for (int t=1;t<=T_MAX; t++){ Act_av[i][t]/=(double)STATS;}
    }
    
    for (int t=1;t<=T_MAX; t++){
        draw_timestep_MOD_Assync(D,t,Act_av,f,YTOP,W,T_MAX);
    }
    sprintf(fn,"%s/%s/_EXP/Assync_%d_TimeCourses/Attr_%u-%u__S-%d",
            MODEL_Directory,D->BRN->path,ASSYNC_BIAS,a_st,a_st_k,STATS);
    close_EPSDRAW_file(f,fn,2);
    
}

void timecourse_with_asynchronous_update_Overlay(Boolean_Dynamics *D,
                                                 unsigned int a_st,unsigned int a_st_k,
                                                 unsigned int STATS){
    
    int W=8,YTOP,XTOP,T_MAX=300;
    FILE *f;
    char fn[600];
    //  unsigned int *pl_ID,*pl_st;
    doublepointer *Act_av;
    
    D->set_state(D->get_attractor_state_string(a_st,a_st_k));
    
   // int id_GF=D->BRN->MDAT->get_node_ID(GF_low_name);
    
    int N_lines=D->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*650+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,D->BRN->path);
    system(fn);
    
    sprintf(fn,"mkdir %s/%s/_EXP/Assync_%d_TimeCourses\n",MODEL_Directory_system,D->BRN->path,ASSYNC_BIAS);
    system(fn);
    
    Act_av = new doublepointer [D->BRN->N+1];
    for (int i=0; i<=D->BRN->N; i++){
        Act_av[i]=new double[T_MAX+1];
        for (int j=0; j<=T_MAX; j++){
            Act_av[i][j]=0;
        }
    }
    
    sprintf(fn,"%s/%s/_EXP/Assync_%d_TimeCourses/%s_Attr_%u-%u__S-%d.eps",
            MODEL_Directory,D->BRN->path,ASSYNC_BIAS,D->BRN->name,a_st,a_st_k,STATS);
    // printf("%s\n",fn);
    f=fopen(fn,"w");
    
    EpsInit(f,-5,-5,XTOP,YTOP);
    EpsSetFont(f,"Times-Roman",W);
    
    for (unsigned int s=1;s<=STATS; s++) {
        D->set_state(D->get_attractor_state_string(a_st,a_st_k));
        
        for (int t=1;t<=T_MAX; t++) {
            if(ASSYNC_BIAS == 0)    D->Asynchronously_update_all_Gates();
            else D->Asynchronously_update_all_Gates(ASSYNC_BIAS,NULL,0);
            
            for (int i=1; i<=D->BRN->N; i++){
                //Act_av[i][t]+=D->s[i-1]-48;
                Act_av[i][t]+=D->s[i];
            }
        }
    }
    
    for (int i=1; i<=D->BRN->N; i++){
        for (int t=1;t<=T_MAX; t++){ Act_av[i][t]/=(double)STATS;}
    }
    
    for (int t=1;t<=T_MAX; t++){
        draw_timestep_MOD_Assync(D,t,Act_av,f,YTOP,W,T_MAX);
    }
    sprintf(fn,"%s/%s/_EXP/Assync_%d_TimeCourses/%s_Attr_%u-%u__S-%d",
            MODEL_Directory_system,D->BRN->path,ASSYNC_BIAS,D->BRN->name,a_st,a_st_k,STATS);
    close_EPSDRAW_file(f,fn,2);
}


void perform_experiment_with_asynchronous_update(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                 unsigned int a_st,unsigned int a_st_k,
                                                 unsigned int PL_nr,pulselist *pl,
                                                 unsigned int STATS){
    
    int W=8,YTOP,XTOP,T_MAX=300;
    FILE *f;
    char fn[600];
    unsigned int *pl_ID,*pl_st;
    doublepointer *Act_av;
    double p;
    // PD->Coupled->print_state();//getchar();
  
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
   // int id_GF=PD->BRN->MDAT->get_node_ID("GF");
    //int s_start_GF=PD->Coupled->s[id_GF-1];
    pl_ID = new unsigned int[PL_nr+1];
    pl_st= new unsigned int[PL_nr+1];
    for (unsigned int i=1;i<=PL_nr; i++) {
        pl_ID[i]=PD->BRN->MDAT->get_node_ID(pl[i]->nodename);
        if(pl_ID[i]<=0) getchar();
        pl_st[i]=PD->Coupled->s[pl_ID[i]-1];
    }
    
    long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a+1;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*650+5;
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    
    char pulsecollection[400];
    sprintf(pulsecollection,"%s_%d_t-%d",pl[1]->nodename,pl[1]->pulsetype,pl[1]->t_end-pl[1]->t_st);
    for (unsigned int i=2;i<=PL_nr; i++) {
        sprintf(pulsecollection,"%s__%s_%u_t-%d",pulsecollection,pl[i]->nodename,pl[i]->pulsetype,pl[i]->t_end-pl[i]->t_st);
    }
    sprintf(fn,"mkdir %s/%s/_EXP/Assync_%d_Pulses_%s\n",MODEL_Directory_system,PD->BRN->path,ASSYNC_BIAS,pulsecollection);
    system(fn);
   
    char namecollection[400];
    sprintf(namecollection,"%s-%d-%d_%.3lf",pl[1]->nodename,pl[1]->t_st,pl[1]->t_end,pl[1]->p_force);
    for (unsigned int i=2;i<=PL_nr; i++) {
        sprintf(namecollection,"%s__%s-%d-%d_%.3lf",
                namecollection,pl[i]->nodename,pl[i]->t_st,pl[i]->t_end,pl[i]->p_force);
    }
    
    Act_av = new doublepointer [PD->Coupled->BRN->N+1];
    for (int i=0; i<=PD->Coupled->BRN->N; i++){
        Act_av[i]=new double[T_MAX+1];
        for (int j=0; j<=T_MAX; j++){
            Act_av[i][j]=0;
        }
    }
    
    for (int delay=1; delay<=DMAX; delay++) {
        sprintf(fn,"%s/%s/_EXP/Assync_%d_Pulses_%s/%s_Delay-%d__S-%d.eps",
                   MODEL_Directory,PD->BRN->path,ASSYNC_BIAS,pulsecollection,namecollection,delay,STATS);
        printf("%s\n",fn);
        f=fopen(fn,"w");
        
        EpsInit(f,-5,-5,XTOP,YTOP);
        EpsSetFont(f,"Times-Roman",W);
  
        for (unsigned int s=1;s<=STATS; s++) {
            PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
            
            for (int t=1;t<=T_MAX; t++) {
               // int resetGF=1;
              //  for (unsigned int l=1; l<=PL_nr; l++)
              //      if(pl_ID[l] == id_GF) resetGF=0;
              //  if(resetGF==1) PD->Coupled->s[id_GF-1]=s_start_GF;
               
                for (unsigned int l=1; l<=PL_nr; l++) {
                    if (t>=pl[l]->t_st+delay) {
                        p=rand_real(0,1);
                        if(p<=pl[l]->p_force){
                            if(pl[l]->pulsetype==forced_on) PD->Coupled->s[pl_ID[l]-1]=49;
                            if(pl[l]->pulsetype==forced_off) PD->Coupled->s[pl_ID[l]-1]=48;
                            if(pl[l]->pulsetype==forced_on_reg) PD->Coupled->s[pl_ID[l]-1]=49;
                            if(pl[l]->pulsetype==forced_off_reg) PD->Coupled->s[pl_ID[l]-1]=48;
                        }
                        else {if(pl[l]->pulsetype==forced_on) PD->Coupled->s[pl_ID[l]-1]=48;}
                    }
                    if (t>=pl[l]->t_end+delay){
                        PD->Coupled->s[pl_ID[l]-1]=pl_st[l];
                    }
                }

                if(ASSYNC_BIAS == 0)    PD->Coupled->Asynchronously_update_all_Gates();
                else PD->Coupled->Asynchronously_update_all_Gates(ASSYNC_BIAS,NULL,0);

                for (int i=1; i<=PD->Coupled->BRN->N; i++){
                    Act_av[i][t]+=PD->Coupled->s[i-1]-48;
                }
            }
        }
        
        for (int i=1; i<=PD->Coupled->BRN->N; i++){
            for (int t=1;t<=T_MAX; t++){ Act_av[i][t]/=(double)STATS;}
        }
        
        for (int t=1;t<=T_MAX; t++){
             draw_timestep_MOD_Assync(PD,t,Act_av,f,YTOP,W,T_MAX);
        }
        sprintf(fn,"%s/%s/_EXP/Assync_%d_Pulses_%s/%s_Delay-%d__S-%d",
                MODEL_Directory_system,PD->BRN->path,ASSYNC_BIAS,pulsecollection,namecollection,delay,STATS);
        close_EPSDRAW_file(f,fn,2);
    }
    
}


void setup_experiment_with_asynchronous_update(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int a_st,int pulseL,int stat){
    
    pulselist *P;
    P=new pulselist[1];
    P[1]=new pulse;
    sprintf(P[1]->nodename,"GF_High");
    P[1]->p_force = 1;
    P[1]->t_st =20;
    P[1]->t_end =P[1]->t_st + pulseL;
    P[1]->pulsetype =forced_on;

    perform_experiment_with_asynchronous_update(PD,a_st,1,1,P,stat);
    delete[] P; P=NULL;
   
    P=new pulselist[1];
    P[1]=new pulse;
    sprintf(P[1]->nodename,"Trail");
    P[1]->p_force = 1;
    P[1]->t_st =20;
    P[1]->t_end =P[1]->t_st + pulseL;
    P[1]->pulsetype =forced_on;
    
    perform_experiment_with_asynchronous_update(PD,a_st,1,1,P,stat);
    delete[] P; P=NULL;
    
    P=new pulselist[1];
    P[1]=new pulse;
    sprintf(P[1]->nodename,"GF");
    P[1]->p_force = 1;
    P[1]->t_st =20;
    P[1]->t_end =P[1]->t_st + pulseL;
    P[1]->pulsetype =forced_off;
    
    perform_experiment_with_asynchronous_update(PD,a_st,1,1,P,stat);
    delete[] P; P=NULL;
}

#endif /* Asynchronous_version_h */
