//
//  Modules.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 1/3/18.
//  Copyright © 2018 Regan_Group. All rights reserved.
//

#ifndef Modules_h
#define Modules_h


void export_Module_Attractors_onto_Boolean_network__PDF(unsigned int mID, MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,Boolean_Dynamics *D){
    unsigned long int i;
    FILE *f;
    char fn[300];
    double XMAX,YMAX,XMIN, YMIN,e_min,e_max,ch_min,ch_max,RMAX=15;
    unsigned long int s;
    snode::slinklist::slink *sc2;
    
    for (unsigned long int b=1; b<=D->Attractor_NR; b++) {
         for(i=1;i<=D->Attractor_Cycle[b]->N;i++) {
            D->set_state(D->Attractor_Cycle[b]->A[i]);
            
             sprintf(fn,"%s/%s/MODULE_%d_%s/%s_Attr_%ld_C-%lu.eps", MODEL_Directory,PD->Coupled->DIRname,mID,PD->BRN->Module_Names[mID].c_str(),PD->Coupled->BRN->name,b,i);
            printf("%s\n",fn);
            
            f=fopen(fn,"w");
            
            COLOR_SCHEME=GYR;
            setcolors_PAJEK();
            
            XMAX=maximum(D->BRN->N, D->BRN->MDAT->x);
            XMIN=minimum(D->BRN->N, D->BRN->MDAT->x);
            YMAX=maximum(D->BRN->N, D->BRN->MDAT->y);
            YMIN=minimum(D->BRN->N, D->BRN->MDAT->y);
            //printf("XMIN=%lg\tXMAX=%lg\nYMIN=%lg\tYMAX=%lg\n",XMIN,XMAX,YMIN,YMAX);
            //getchar();
            EpsInit(f,-25,-25,XMAX-XMIN+25,YMAX-YMIN+25);
            e_min=2000000;e_max=0;
            ch_min=2000000;ch_max=-2000000;
            
            EpsSetRgb(f,0.2,0.2,0.2);
            EpsSetLinewidth(f, 1);
            EpsSetFont(f,"Times-Roman",RMAX);
            int linksign;
            for(unsigned int s=1;s<=D->BRN->N;s++){
                int k=0;
                sc2=D->BRN->BN->node[s].Li.first_link_pt;
                while(sc2!=NULL){
                    k++;
                    // printf("Link %s (%d) -> %s (%llu) (input # %d)\n",BRN->MDAT->Node_name[s],s, BRN->MDAT->Node_name[sc2->sneighbor],sc2->sneighbor,k);
                    linksign=D->BRN->Gate[s]->check_overall_sign_of_input(k);
                    switch (linksign) {
                        case 1:  EpsSetRgb(f,0,0,0.8);break;
                        case -1: EpsSetRgb(f,0.8,0,0);break;
                        default: EpsSetRgb(f,0.2,0.2,0.2); break;
                    }
                    EpsDrawLine(f,D->BRN->MDAT->x[s]-XMIN,
                                   YMAX-D->BRN->MDAT->y[s],
                                   D->BRN->MDAT->x[sc2->sneighbor]-XMIN,
                                    YMAX-D->BRN->MDAT->y[sc2->sneighbor]);
                    EpsFillCircle(f,D->BRN->MDAT->x[s]-XMIN - (D->BRN->MDAT->x[s] - D->BRN->MDAT->x[sc2->sneighbor])/4.,
                                  YMAX-D->BRN->MDAT->y[s] + (D->BRN->MDAT->y[s] - D->BRN->MDAT->y[sc2->sneighbor])/4.,3);
                    sc2=sc2->next(sc2);
                }
            }
            
            //all nodes with energy color, size and attracror border
            
            EpsSetLinewidth(f,3);
            
            for(s=1;s<=D->BRN->N;s++){
                switch(D->s[s]){
                        //      case 48:  EpsSetRgb(f,0,1,0);break; //0
                        //      case 49:  EpsSetRgb(f,1,0,0);break; //1
                    case 48:  EpsSetRgb(f,0.5,0.5,1);break; //0
                    case 49:  EpsSetRgb(f,1,0.5,0);break; //1
                    case 0:  EpsSetRgb(f,0.5,0.5,1);break; //0
                    case 1:  EpsSetRgb(f,1,0.5,0);break; //1
                    default:  EpsSetRgb(f,1,1,1);break;
                }
                //printf("node %s\t color %d\n",BRN->MDAT->Node_name[s],Attractor_valleys[b]->Basin_Stored[Attractor_valleys[b]->Attractor[i]]->s[s]);
                
                EpsFillCircle(f,D->BRN->MDAT->x[s]-XMIN,YMAX-D->BRN->MDAT->y[s],RMAX);
                
                if (D->BRN->MDAT->Module_ID!=NULL) {
                    switch(D->BRN->MDAT->Module_ID[s]){
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
                EpsDrawCircle(f,D->BRN->MDAT->x[s]-XMIN,YMAX-D->BRN->MDAT->y[s],RMAX);
                EpsSetRgb(f,0,0,0);
                EpsDrawString_Centered(f, D->BRN->MDAT->x[s]-XMIN,YMAX-D->BRN->MDAT->y[s], D->BRN->MDAT->Node_name[s]);
            }
            
            sprintf(fn,"%s/%s/MODULE_%d_%s/%s_Attr_%ld_C-%lu", MODEL_Directory_system,PD->Coupled->DIRname,mID,PD->BRN->Module_Names[mID].c_str(),PD->Coupled->BRN->name,b,i);
            close_EPSDRAW_file(f,fn,2);
        }
    }
    //  printf("\ndone\n");
}

void set_up_module_MetaData(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    
    for(int j=1;j<=PD->BRN->N;j++){
        printf("j=%d  Mod ID=%d ",j,PD->Module_ID[j]);
        printf("ID_in_Mod=%d\n",PD->ID_in_Module[j]);
        
        if(PD->Module_ID[j]>0) {
            PD->BRN_Module[PD->Module_ID[j]]->MDAT->x[PD->ID_in_Module[j]]=PD->BRN->MDAT->x[j];
            PD->BRN_Module[PD->Module_ID[j]]->MDAT->y[PD->ID_in_Module[j]]=PD->BRN->MDAT->y[j];
            strcpy(PD->BRN_Module[PD->Module_ID[j]]->MDAT->Node_name[PD->ID_in_Module[j]],PD->BRN->MDAT->Node_name[j]);
        }
    }
    for(int m=1;m<=PD->Coupled->BRN->Module_Order->N;m++){
        if((PD->Coupled->BRN->Module_Order->A[m]>0)&&(PD->Coupled->BRN->Module_Order->A[m]<=PD->Coupled->BRN->MDAT->Module_NR)){
            printf("Dealing with module order %d  - id %ld\n",m,PD->Coupled->BRN->Module_Order->A[m]);
            long int module_now = PD->Coupled->BRN->Module_Order->A[m];
            PD->BRN_Module[module_now]->Module_Order =  new longint_ARRAY_pair();
            PD->BRN_Module[module_now]->Module_Order->add_element(1,0);
            PD->BRN_Module[module_now]->MDAT->Module_NR=1;
            PD->BRN_Module[module_now]->MDAT->Module_ID = new unsigned int[PD->BRN_Module[module_now]->N+1];
            for(int nn=1;nn<=PD->BRN_Module[module_now]->N;nn++) PD->BRN_Module[module_now]->MDAT->Module_ID[nn]=1;
            //printf("\t adding (1 0) to module order: better match (%ld 0)\n",PD->Coupled->BRN->Module_Order->A[m]);
            PD->BRN_Module[module_now]->Draw_in_Modules=new longint_ARRAY_list[2];
            PD->BRN_Module[module_now]->Draw_in_Modules[1]=new longint_ARRAY();
            PD->BRN_Module[module_now]->active_inputs=new longint_ARRAY();
            PD->BRN_Module[module_now]->silenced_inputs=new longint_ARRAY();
            PD->BRN_Module[module_now]->active_inputs_noGF=new longint_ARRAY();

            for(int nn=1;nn<=PD->Coupled->BRN->Draw_in_Modules[module_now]->N;nn++){
                long int i=PD->Coupled->BRN->Draw_in_Modules[module_now]->A[nn];
                long int i_mod = PD->BRN_Module[module_now]->MDAT->get_node_ID(PD->Coupled->BRN->MDAT->Node_name[i]);
                PD->BRN_Module[module_now]->Draw_in_Modules[1]->add_element(i_mod);
                printf("\t %ld %ld: adding node %s (%s) to module drawing order, poz %d\n",i, i_mod,PD->Coupled->BRN->MDAT->Node_name[i],
                       PD->BRN_Module[module_now]->MDAT->Node_name[i_mod],nn);
            }
        }
    }
   // getchar();
}

void work_with_individual_Modules(int mID, MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    Boolean_Dynamics *D;
  //  snode::slinklist::slink *sc2;
    
    double beta=0.5 *log(1./p_error-1.);
    printf("p_err=%lg\t beta=%lg\n",p_error,beta);//getchar();
    
    char fn[300];
    sprintf(fn,"mkdir %s/%s/MODULE_%d_%s\n",MODEL_Directory_system,PD->Coupled->DIRname,mID,PD->BRN->Module_Names[mID].c_str());
    system(fn);
    
    D=new Boolean_Dynamics(PD->BRN_Module[mID]);
    PD->BRN_Module[mID]->export_Boolean_RegNet();
   // sprintf(fn, "%s/%s/MODULE_%d_%s/Module_%d_%s_GATES.txt",MODEL_Directory,PD->Coupled->DIRname,mID,PD->BRN->Module_Names[mID].c_str(),mID,PD->BRN->Module_Names[mID].c_str());
    PD->BRN_Module[mID]->print_gates_to_File();
    PD->BRN_Module[mID]->export_Boolean_RegNet_with_names_to_Cytoskape();
    
    longint_ARRAY *SCC;
    SCC=PD->BRN_Module[mID]->Find_Strongly_Connected_Component();
    for(unsigned int i=1;i<=PD->BRN_Module[mID]->N;i++){
        if (SCC->check_for_element(i)==0) {
            printf("Node %s is not in the module's SCC\n",PD->BRN_Module[mID]->MDAT->Node_name[i]);
        }
    }
 //   getchar();
    
    if(D->BRN->N<=20){
        D->generate_full_state_graph_synchronously(0);
//        D->calculate_noisy_landscape(beta);
        D->export_state_Graph_to_Cytoskape_no_noise();
        export_Module_Attractors_onto_Boolean_network__PDF(mID, PD,D);
        delete D;D=NULL;
    }
    else {
        Boolean_Dynamics_TRACKER *D2;
        char fn[400];
        
        D2=new Boolean_Dynamics_TRACKER(PD->BRN_Module[mID],500);
        sprintf(fn,"mkdir %s/%s\n",MODEL_Directory_system,D2->DIRname);
        system(fn);
        D2->record_timetraces_from_collection_of_random_inputs(100,100,0.01);
        D2->BRN->export_Boolean_RegNet_with_names_to_Cytoskape();
        D2->Export_ALL_Attractors_onto_Boolean_Network(0);
        sprintf(fn,"rm -r %s/%s/%s\n",MODEL_Directory_system,PD->Coupled->DIRname,D2->DIRname);
        system(fn);
        sprintf(fn,"mv %s/%s %s/%s\n",MODEL_Directory_system,D2->DIRname, MODEL_Directory_system,PD->Coupled->DIRname);
        system(fn);
        //printf("%s",fn);getchar();
    }

}


void Async_Module_cycle(int mID, MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error,unsigned int STATS){
    
    double beta=0.5 *log(1./p_error-1.);
    printf("p_err=%lg\t beta=%lg\n",p_error,beta);//getchar();
    
    char fn[300];
    sprintf(fn,"mkdir %s/%s/MODULE_%d_%s\n",MODEL_Directory_system,PD->Coupled->DIRname,mID,PD->BRN->Module_Names[mID].c_str());
    system(fn);
    
    PD->BRN_Module[mID]->export_Boolean_RegNet();
    //sprintf(fn, "%s/%s/MODULE_%d_%s/Module_%d_%s_GATES.txt",MODEL_Directory,PD->Coupled->DIRname,mID,PD->BRN->Module_Names[mID].c_str(),mID,PD->BRN->Module_Names[mID].c_str());
   
    // PD->BRN_Module[mID]->print_gates_to_File(fn);
    PD->BRN_Module[mID]->print_gates_to_File();
    PD->BRN_Module[mID]->export_Boolean_RegNet_with_names_to_Cytoskape();
    
    longint_ARRAY *SCC;
    SCC=PD->BRN_Module[mID]->Find_Strongly_Connected_Component();
    for(unsigned int i=1;i<=PD->BRN_Module[mID]->N;i++){
        if (SCC->check_for_element(i)==0) {
            printf("Node %s is not in the module's SCC\n",PD->BRN_Module[mID]->MDAT->Node_name[i]);
        }
    }
    //   getchar();
    
    if( PD->BRN_Module[mID]->N<=20){
        Boolean_Dynamics *D;
        D=new Boolean_Dynamics(PD->BRN_Module[mID]);
        D->generate_full_state_graph_synchronously(0);
        //        D->calculate_noisy_landscape(beta);
        D->export_state_Graph_to_Cytoskape_no_noise();
        export_Module_Attractors_onto_Boolean_network__PDF(mID, PD,D);
        for (unsigned int k=1; k<=D->Attractor_NR; k++) {
            printf("k=%d\t cycle: %llu\n",k,D->Attractor_Cycle[k]->N);
            if(D->Attractor_Cycle[k]->N>1){
                printf("D->BRN->Module_Order->A[1] = %lu\n",D->BRN->Module_Order->A[1]);
                if((D->BRN->Module_Order->A[1]>0)&&(D->BRN->Module_Order->A[1]<=D->BRN->MDAT->Module_NR)){
                    printf("\tNodes: %ld \n",D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[1]]->N);
                    for(int nn=1;nn<=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[1]]->N;nn++){
                        long int i=D->BRN->Draw_in_Modules[D->BRN->Module_Order->A[1]]->A[nn];
                        printf("\t\t: %lu %s \n",i,D->BRN->MDAT->Node_name[i]);
                    }
                }
              //  getchar();
                timecourse_with_asynchronous_update_Overlay(D,k,1,STATS);
            }
        }
        delete D;D=NULL;
    }
    else {
        Boolean_Dynamics_TRACKER *D2;
        char fn[400];
        
        D2=new Boolean_Dynamics_TRACKER(PD->BRN_Module[mID],500);
        sprintf(fn,"mkdir %s/%s\n",MODEL_Directory_system,D2->DIRname);
        system(fn);
        D2->record_timetraces_from_collection_of_random_inputs(100,100,0.01);
        D2->BRN->export_Boolean_RegNet_with_names_to_Cytoskape();
        D2->Export_ALL_Attractors_onto_Boolean_Network(0);
        for (unsigned int k=1; k<=D2->Attractor_NR; k++) {
            if(D2->Attractor_valleys[k]->N_a>1)
                timecourse_with_asynchronous_update_Overlay(D2,k,1,STATS);
        }
        delete D2;D2=NULL;
    }
    
}

#endif /* Modules_h */
