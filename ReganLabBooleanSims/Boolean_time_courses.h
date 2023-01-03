//
//  Boolean_time_courses.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 9/21/16.
//  Copyright © 2016 Regan_Group. All rights reserved.
//

#ifndef Boolean_time_courses_h
#define Boolean_time_courses_h


void draw_timesptep(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                        longint_ARRAY_list *node_order,
                        int t, FILE *f,double y0,double y02,int W,int T_MAX,int TW,int NameBuff){
    
    int i;
    int poz=0;
    for (int m=0; m<=PD->Coupled->BRN->MDAT->Module_NR; m++) {
        for(k=1;k<=node_order[m]->N;k++){
            i=node_order[m]->A[k];
            if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,0.8); else EpsSetRgb(f,1,0,0);
            EpsFillRectangle(f,NameBuff+t*TW,y0-poz*W-W,NameBuff+(t+1)*TW,y0-poz*W);
      //      EpsSetRgb(f,0,0,0);
      //     EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
            poz++;
        }
        poz++;
    }
}

void draw_nodenames(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                        longint_ARRAY_list *node_order,
                        FILE *f,double y0,double y02,int W,int T_MAX,int TW,int NameBuff){
    
    int i;
    int poz=0;
    for (int m=0; m<=PD->Coupled->BRN->MDAT->Module_NR; m++) {
        for(k=1;k<=node_order[m]->N;k++){
            i=node_order[m]->A[k];
            EpsSetRgb(f,0,0,0);
            EpsDrawString(f,0,0,y0-poz*W-W,PD->BRN->MDAT->Node_name[i]);
            poz++;
        }
        poz++;
    }
}

void draw_timetrace_ON_OFF_single_INPUT_pulse(int a_st,MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error,int delay,int id_input,FILE *f,double y0,FILE *f2,double y02,int W){
    // delay is the length of the
    int T_MAX=200,TW=3,NameBuff=30;
    longint_ARRAY_list *node_order;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[1]]->s);
    
    node_order=Generate_drawing_node_order(PD->Coupled->BRN);

    draw_nodenames(PD,node_order,idgf,f,y0,y02,W, T_MAX, TW, NameBuff);

    for (int t=1;t<=T_MAX; t++) {
        if (t<=10) PD->Coupled->s[id_input-1]=48;
        if ((t>10)&&(t<=10+delay)) PD->Coupled->s[id_input-1]=49;
        if (t>10+delay) PD->Coupled->s[id_input-1]=48;
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        draw_timesptep(PD,node_order,t,idgf,f,y0,y02,W, T_MAX, TW, NameBuff);
    }
    
    for (int m=0; m<=PD->Coupled->BRN->MDAT->Module_NR; m++) {
        delete node_order[m]; node_order[m]=NULL;
    }
    delete[] node_order; node_order=NULL;
}



#endif /* Boolean_time_courses_h */
