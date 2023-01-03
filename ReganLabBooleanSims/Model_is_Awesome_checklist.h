//
//  Model_is_Awesome_checklist.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 12/23/17.
//  Copyright © 2017 Regan_Group. All rights reserved.
//

#ifndef Model_is_Awesome_checklist_h
#define Model_is_Awesome_checklist_h

#define SEN_CUT 1


unsigned short int similar_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                      longint_ARRAY *gl_keep_asis,
                                      Attractor_char_list *AL){
    unsigned long int n=0;
    longint_ARRAY *l_sim;
    l_sim= Find_Similar_Attractors__Life_and_Death(PD,gl_keep_asis,1);
    n=l_sim->N;
    delete l_sim; l_sim=NULL;
    return(n);
}

unsigned short int dead_when_needed(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                        longint_ARRAY *gl_keep_asis,
                                        Attractor_char_list *AL){
    unsigned long int n=0;
    unsigned long int GFi=PD->BRN->MDAT->get_node_ID("GF");
    unsigned long int UV=PD->BRN->MDAT->get_node_ID("UV");
    unsigned long int Trail=PD->BRN->MDAT->get_node_ID("Trail");
    
    unsigned long int ll_GF=0,ll_UV=0,ll_Trail=0;
    for(unsigned int ll=1;ll<=gl_keep_asis->N;ll++){
        if (GFi == gl_keep_asis->A[ll]) ll_GF=ll;
        if (UV == gl_keep_asis->A[ll]) ll_UV=ll;
        if (Trail == gl_keep_asis->A[ll]) ll_Trail=ll;
    }
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        if ((ll_GF>0)&&(AL[b]->env_now[ll_GF] == 0) && (AL[b]->Phenotype.Dead==0)){
            printf("\t\t %u is Alive in no GF environment!\n",b);
            n++;
        }
        if ((ll_UV>0)&&(AL[b]->env_now[ll_UV] == 1) && (AL[b]->Phenotype.Dead==0)){
            printf("\t\t %u is Alive in UV!\n",b);
            n++;
        }
        if ((ll_Trail>0)&&(AL[b]->env_now[ll_Trail] == 1) && (AL[b]->Phenotype.Dead==0)
            && (AL[b]->Phenotype.Sen==0)){
            printf("\t\t %u is Alive in Trail, not Senescent!\n",b);
            n++;
        }
    }
    return(n);
}

unsigned short int alive_when_needed(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                    longint_ARRAY *gl_keep_asis,
                                    Attractor_char_list *AL){
    unsigned long int n=0;
    unsigned long int GFi=PD->BRN->MDAT->get_node_ID("GF");
    unsigned long int Gamma=PD->BRN->MDAT->get_node_ID("Gamma");
    unsigned long int UV=PD->BRN->MDAT->get_node_ID("UV");
    unsigned long int Trail=PD->BRN->MDAT->get_node_ID("Trail");
    unsigned long int ROS=PD->BRN->MDAT->get_node_ID("ROS_ext");
    unsigned long int GFH=PD->BRN->MDAT->get_node_ID("GF_High");
    
    unsigned long int ll_GF=0,ll_Gamma=0,ll_UV=0,ll_Trail=0,ll_ROS=0,ll_GFH=0;
    for(unsigned int ll=1;ll<=gl_keep_asis->N;ll++){
        if (GFi == gl_keep_asis->A[ll]) ll_GF=ll;
        if (Gamma == gl_keep_asis->A[ll]) ll_Gamma=ll;
        if (UV == gl_keep_asis->A[ll]) ll_UV=ll;
        if (Trail == gl_keep_asis->A[ll]) ll_Trail=ll;
        if (ROS == gl_keep_asis->A[ll]) ll_ROS=ll;
        if (GFH == gl_keep_asis->A[ll]) ll_GFH=ll;
    }
    
    int ok_lowGF=0,ok_lowGF_Gamma=0,ok_highGF_Gamma=0;
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        if ((ll_GF>0)&&(AL[b]->env_now[ll_GF] == 1) && (AL[b]->Phenotype.Dead==0)) ok_lowGF++;
        if (((ll_Gamma==0)||(AL[b]->env_now[ll_Gamma] == 0))&&
            ((ll_UV==0)||(AL[b]->env_now[ll_UV] == 0))&&
            ((ll_Trail==0)||(AL[b]->env_now[ll_Trail] == 0))&&
            ((ll_ROS==0)||(AL[b]->env_now[ll_ROS] == 0))&&
            ((ll_Gamma==0)||(AL[b]->env_now[ll_Gamma] == 0))){// none of the stressors are on
           
            if(PD->Coupled->Attr_Transitions[ll_GFH]!=NULL){
                snode::slinklist::slink *sc2;
                sc2=PD->Coupled->Attr_Transitions[ll_GFH]->node[b].Li.first_link_pt;
                while(sc2!=NULL){
                    if(AL[b]->Phenotype.Dead + AL[sc2->sneighbor]->Phenotype.Dead ==1 ){
                        printf("\t\t Dies going back and forth from healthy GF High / low: %d <-> %llu !\n",
                               b,sc2->sneighbor);
                        n++;
                    }
                    sc2=sc2->next(sc2);
                }
            }
            else {printf("PD_A->Coupled->Attr_Transitions[%lu] does not exist!!!\n",ll_GF);getchar();}
        }
        
        if (((ll_Gamma>0)&&(AL[b]->env_now[ll_Gamma] == 1))&&
            ((ll_GF>0)&&(AL[b]->env_now[ll_GF] != 0))){
                if ((ll_GFH>0)&&(AL[b]->env_now[ll_GFH] == 1)) ok_highGF_Gamma++;
                if ((ll_GFH>0)&&(AL[b]->env_now[ll_GFH] == 0)) ok_lowGF_Gamma++;
        }
    }
    
    if (n==0) printf("\t\t\t Alive when going back and forth between low-high GF -- OK!\n");
    
    int go_gamma_ok=1;
    longint_ARRAY *l;
    l=get_Healthy_G0_attractors(PD);
    for(int j=1;j<=l->N;j++){
        if(PD->Coupled->Attr_Transitions[ll_Gamma]!=NULL){
            snode::slinklist::slink *sc2;
            sc2=PD->Coupled->Attr_Transitions[ll_Gamma]->node[l->A[j]].Li.first_link_pt;
            while(sc2!=NULL){
                if((AL[l->A[j]]->Phenotype.Dead ==0)&&(AL[sc2->sneighbor]->Phenotype.Dead == 1)){
                    printf("\t\t Dies when exposed to Gamma from healthy GO! %ld -> %llu !\n",
                           l->A[j],sc2->sneighbor);
                    n++;go_gamma_ok=0;
                }
                sc2=sc2->next(sc2);
            }
        }
        else {printf("PD_A->Coupled->Attr_Transitions[%lu] does not exist!!!\n",ll_Gamma);getchar();}
    }
    
    if (ok_lowGF==0) {printf("\t\t No alive states in low GF!\n");n++;}
        else printf("\t\t\t Alive in low GF -- OK!\n");
    if (ok_lowGF_Gamma==0) {printf("\t\t No alive states in low GF and Gamma!\n");n++;}
        else printf("\t\t\t Alive in low GF & Gamma -- OK!\n");
    if (ok_highGF_Gamma==0) {printf("\t\t No alive states in High GF and Gamma!\n");n++;}
        else printf("\t\t\t Alive in high GF & Gamma -- OK!\n");
    if (go_gamma_ok==1) printf("\t\t\t Alive when hit with gamma in G0 -- OK!\n");
    return(n);
}

unsigned short int irreversible_senesence_3way(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                     longint_ARRAY *gl_keep_asis,
                                     Attractor_char_list *AL){
    unsigned long int n=0;
 /*   unsigned long int GFi=PD->BRN->MDAT->get_node_ID("GF");
    unsigned long int Gamma=PD->BRN->MDAT->get_node_ID("Gamma");
    unsigned long int UV=PD->BRN->MDAT->get_node_ID("UV");
    unsigned long int Trail=PD->BRN->MDAT->get_node_ID("Trail");
    unsigned long int ROS=PD->BRN->MDAT->get_node_ID("ROS_ext");
    unsigned long int GFH=PD->BRN->MDAT->get_node_ID("GF_High");
    
    unsigned long int ll_GF=0,ll_Gamma=0,ll_UV=0,ll_Trail=0,ll_ROS=0,ll_GFH=0;
    for(unsigned int ll=1;ll<=gl_keep_asis->N;ll++){
        if (GF == gl_keep_asis->A[ll]) ll_GF=ll;
        if (Gamma == gl_keep_asis->A[ll]) ll_Gamma=ll;
        if (UV == gl_keep_asis->A[ll]) ll_UV=ll;
        if (Trail == gl_keep_asis->A[ll]) ll_Trail=ll;
        if (ROS == gl_keep_asis->A[ll]) ll_ROS=ll;
        if (GFH == gl_keep_asis->A[ll]) ll_GFH=ll;
    }
    */
    
    int Mok=0,Chok=0,both=0;
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        if(AL[b]->Phenotype.Sen==1) Mok++;
        if(AL[b]->Phenotype.Sen==2) Chok++;
        if(AL[b]->Phenotype.Sen==3) both++;

        if(AL[b]->Phenotype.Sen>0){
            for(unsigned int ll=1;b<=gl_keep_asis->N;ll++) {
                if(PD->Coupled->Attr_Transitions[ll]!=NULL){
                    snode::slinklist::slink *sc2;
                    sc2=PD->Coupled->Attr_Transitions[ll]->node[b].Li.first_link_pt;
                    while(sc2!=NULL){
                        if((AL[sc2->sneighbor]->Phenotype.Dead == 0) &&(AL[b]->Phenotype.Sen==0)){
                            switch (AL[b]->Phenotype.Sen) {
                                case 1:{printf("\t\t Mitochondrial senesence was revesible!: %d -> %llu with a %s pulse length %lg!\n",
                                               b,sc2->sneighbor,
                                               PD->Coupled->BRN->MDAT->Node_name[sc2->link_props->distinc_link_labels], sc2->link_props->lw);
                                        n++;
                                }  break;
                                case 2:{printf("\t\t Chromosomal senesence was revesible!: %d -> %llu with a %s pulse length %lg!\n",
                                               b,sc2->sneighbor,
                                               PD->Coupled->BRN->MDAT->Node_name[sc2->link_props->distinc_link_labels], sc2->link_props->lw);
                                    n++;
                                } break;
                                case 3:{printf("\t\t Double senesence was revesible!: %d -> %llu with a %s pulse length %lg!\n",
                                               b,sc2->sneighbor,
                                               PD->Coupled->BRN->MDAT->Node_name[sc2->link_props->distinc_link_labels], sc2->link_props->lw);
                                    n++;
                                } break;
                                default:
                                    break;
                            }
                            
                            n++;
                        }
                        sc2=sc2->next(sc2);
                    }
                }
                else {printf("PD_A->Coupled->Attr_Transitions[%u] does not exist!!!\n",ll);getchar();}
            }
        }
    }
    if (Mok==0) {printf("\t\t NO Mitochondrial Senesence found!\n");n++;}
    if (Chok==0) {printf("\t\t NO Chromosomal Senesence found!\n");n++;}
    if (both==0) {printf("\t\t NO Double Senesence found!\n");n++;}
    
    if (n==0) printf("\t\t\t Senesence exists in 3 ways and it is irreversible -- OK!\n");
    return(n);
}


unsigned short int normal_cell_Cycle_only(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                     longint_ARRAY *gl_keep_asis,
                                     Attractor_char_list *AL){
    unsigned long int n=0;
    unsigned long int GFi=PD->BRN->MDAT->get_node_ID("GF");
    unsigned long int Gamma=PD->BRN->MDAT->get_node_ID("Gamma");
    unsigned long int UV=PD->BRN->MDAT->get_node_ID("UV");
    unsigned long int Trail=PD->BRN->MDAT->get_node_ID("Trail");
    unsigned long int ROS=PD->BRN->MDAT->get_node_ID("ROS_ext");
    unsigned long int GFH=PD->BRN->MDAT->get_node_ID("GF_High");
    
    unsigned long int ll_GF=0,ll_Gamma=0,ll_UV=0,ll_Trail=0,ll_ROS=0,ll_GFH=0;
    for(unsigned int ll=1;ll<=gl_keep_asis->N;ll++){
        if (GFi == gl_keep_asis->A[ll]) ll_GF=ll;
        if (Gamma == gl_keep_asis->A[ll]) ll_Gamma=ll;
        if (UV == gl_keep_asis->A[ll]) ll_UV=ll;
        if (Trail == gl_keep_asis->A[ll]) ll_Trail=ll;
        if (ROS == gl_keep_asis->A[ll]) ll_ROS=ll;
        if (GFH == gl_keep_asis->A[ll]) ll_GFH=ll;
    }
    
   // int ok_lowGF=0,ok_lowGF_Gamma=0,ok_highGF_Gamma=0;
    int ncc=0;
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        if (((ll_GFH>0)&&(AL[b]->env_now[ll_GFH] == 1))&&
            ((ll_Gamma==0)||(AL[b]->env_now[ll_Gamma] == 0))&&
            ((ll_UV==0)||(AL[b]->env_now[ll_UV] == 0))&&
            ((ll_Trail==0)||(AL[b]->env_now[ll_Trail] == 0))&&
            ((ll_ROS==0)||(AL[b]->env_now[ll_ROS] == 0))&&
            ((ll_Gamma==0)||(AL[b]->env_now[ll_Gamma] == 0))){// none of the stressors are on, GF is high
            
            switch (AL[b]->Phenotype.CC) {
                case 0:{ncc++;} break;
                case 1: break;
                case 2: break;
                case 3: break;
                case 4:{printf("\t\t Broken cycle found! -- %u\n",b);n++;} break;
                case 5:{printf("\t\t One CC swithc stuck on a barrier (not known arrest)! -- %u\n",b);n++;} break;
                case 6:{printf("\t\t CC state could not be determined (not known arrest)! -- %u\n",b);n++;} break;
                default: break;
            }
        }
    }
    if (ncc==0) {printf("\t\t NO cell cycle found!\n");n++;}
    if (ncc>1) {printf("\t\t More than one cell cycle found: %d!\n",ncc);n++;}
    if (n==0) printf("\t\t\t Only one cell cycle -- OK!\n");
    return(n);
}
unsigned short int normal_cell_Cycle_entry_from_1_G0(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                          longint_ARRAY *gl_keep_asis,
                                          Attractor_char_list *AL){
    unsigned long int n=0;
    longint_ARRAY *l;
    l=get_Healthy_G0_attractors(PD);
    if (l->N==0) {printf("\t\t NO healthy G0 cycle found in no-stress low GF!\n");n++;}
    if (l->N>1) {printf("\t\t More than one healthy G0 cycle found in no-stress low GF!: %lu!\n",l->N);n++;}
    
    unsigned long int HGF=PD->BRN->MDAT->get_node_ID("GF_High");
    unsigned long int ll_GFH=0;
    for(unsigned int ll=1;ll<=gl_keep_asis->N;ll++){
        if (HGF == gl_keep_asis->A[ll]) ll_GFH=ll;
    }
    
    for(int j=1;j<=l->N;j++){
        int ok=0;
        if(PD->Coupled->Attr_Transitions[ll_GFH]!=NULL){
            snode::slinklist::slink *sc2;
            sc2=PD->Coupled->Attr_Transitions[ll_GFH]->node[l->A[j]].Li.first_link_pt;
            while(sc2!=NULL){
                if(AL[sc2->sneighbor]->Phenotype.CC == 0){
                    ok++;
                }
                sc2=sc2->next(sc2);
            }
        }
        else {printf("PD_A->Coupled->Attr_Transitions[%lu] does not exist!!!\n",ll_GFH);getchar();}
        if (ok==0) {printf("\t\t Healthy G0 %ld cannot enter the cell cycle!\n",l->A[j]);n++;}
    }
    if (n==0) printf("\t\t\t Only one healthy G0, and it enters the CC -- OK!\n");
    return (n);
}

unsigned short int modules_test(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                     longint_ARRAY *gl_keep_asis,
                                                     Attractor_char_list *AL){
    unsigned long int n=0;
//    longint_ARRAY *l;
//    l=get_Healthy_G0_attractors(PD,gl_keep_asis);
//    if (l->N==0) {printf("\t\t NO healthy G0 cycle found in no-stress low GF!\n");n++;}
//    if (l->N>1) {printf("\t\t More than one healthy G0 cycle found in no-stress low GF!: %lu!\n",l->N);n++;}
//
//    unsigned long int HGF=PD->BRN->MDAT->get_node_ID("GF_High");
//    unsigned long int ll_GFH=0;
//    for(unsigned int ll=1;ll<=gl_keep_asis->N;ll++){
//        if (HGF == gl_keep_asis->A[ll]) ll_GFH=ll;
//    }
//
//    for(int j=1;j<=l->N;j++){
//        int ok=0;
//        if(PD->Coupled->Attr_Transitions[ll_GFH]!=NULL){
//            snode::slinklist::slink *sc2;
//            sc2=PD->Coupled->Attr_Transitions[ll_GFH]->node[l->A[j]].Li.first_link_pt;
//            while(sc2!=NULL){
//                if(AL[sc2->sneighbor]->Phenotype.CC == 0){
//                    ok++;
//                }
//                sc2=sc2->next(sc2);
//            }
//        }
//        else {printf("PD_A->Coupled->Attr_Transitions[%lu] does not exist!!!\n",ll_GFH);getchar();}
//        if (ok==0) {printf("\t\t Healthy G0 %ld cannot enter the cell cycle!\n",l->A[j]);n++;}
//    }
//    if (n==0) printf("\t\t\t Only one healthy G0, and it enters the CC -- OK!\n");
    return (n);
}



unsigned short int run_Model_is_Awesome_checklist(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                    longint_ARRAY *gl_keep_asis,int mutant){
    Attractor_char_list *AL;
    unsigned short int not_ok, failed;
    failed=0;
    
    AL=grab_attractor_charact_Coupled(PD,gl_keep_asis);
    
    not_ok=similar_attractors(PD,gl_keep_asis,AL); if (not_ok>0) {printf("\t Fail -- similar attractors: %d\n",not_ok); failed++;}
    
    not_ok=dead_when_needed(PD,gl_keep_asis,AL);  if (not_ok>0) {printf("\t Failed dead test: %d\n",not_ok); failed++;}

    not_ok=alive_when_needed(PD,gl_keep_asis,AL);  if (not_ok>0) {printf("\t Failed alive test: %d\n",not_ok); failed++;}
  
    not_ok=normal_cell_Cycle_only(PD,gl_keep_asis,AL);  if (not_ok>0) {printf("\t Failed Cell cycle test: %d\n",not_ok); failed++;}
  
    not_ok=normal_cell_Cycle_entry_from_1_G0(PD,gl_keep_asis,AL);  if (not_ok>0) {printf("\t Failed G0 test: %d\n",not_ok); failed++;}
    
    if (mutant!=SEN_CUT) {
        not_ok=irreversible_senesence_3way(PD,gl_keep_asis,AL);  if (not_ok>0) {printf("\t Failed Senesence test: %d\n",not_ok); failed++;}
    }
    
    not_ok=modules_test(PD,gl_keep_asis,AL);  if (not_ok>0) {printf("\t Failed modules test: %d\n",not_ok); failed++;}
    
    return (failed);
}

void run_Models_are_Awesome_checklist(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                    MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_B,
                                    longint_ARRAY *gl_keep_asis_A,
                                    longint_ARRAY *gl_keep_asis_B){
    unsigned short int failed_L, failed_S;
    
    printf("\n\n\n Models Are Awesome Checklists\n");
    printf("\nLarge Model\n");
    failed_L = run_Model_is_Awesome_checklist(PD_A,gl_keep_asis_A,0);
    printf("\n\tFailed %d checks!\n",failed_L);
    getchar();
    
    printf("\nSubset Model\n");
    failed_S = run_Model_is_Awesome_checklist(PD_B,gl_keep_asis_B,SEN_CUT);
    printf("\n\tFailed %d checks!\n",failed_S);
     
    getchar();
}

void run_Model_is_Awesome_checklist(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                      longint_ARRAY *gl_keep_asis_A){
    unsigned short int failed_L;
    
    printf("\n\n\n Models Are Awesome Checklists\n");
    printf("\nLarge Model\n");
    failed_L = run_Model_is_Awesome_checklist(PD_A,gl_keep_asis_A,0);
    printf("\n\tFailed %d checks!\n",failed_L);
   // getchar();
    
}

#endif /* Model_is_Awesome_checklist_h */
