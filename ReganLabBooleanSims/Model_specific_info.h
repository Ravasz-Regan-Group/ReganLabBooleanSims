//
//  Model_specific_info.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 9/7/16.
//  Copyright © 2016 Regan_Group. All rights reserved.
//

#ifndef Model_specific_info_h
#define Model_specific_info_h

struct Phenotype_now {
    int Dead,Sen,CC,OtherC,G2_struck,ESC,DNA_damage,Migr;
};

struct Attractor_char_ALL {
    Phenotype_now Phenotype;
    short int env_now[50];
};

int identical_phenotype(Phenotype_now a, Phenotype_now b){
    if(a.Dead!=b.Dead) return(0);
    if(a.Sen!=b.Sen) return(0);
    if(a.CC!=b.CC) return(0);
    if(a.OtherC!=b.OtherC) return(0);
    if(a.G2_struck!=b.G2_struck) return(0);
    if(a.ESC!=b.ESC) return(0);
    if(a.DNA_damage!=b.DNA_damage) return(0);
  //  if(a.Migr!=b.Migr) return(0);
    return(1);
}
typedef  Attractor_char_ALL * Attractor_char_list;


// setting the nodes that represent the environment
longint_ARRAY *get_input_node_names(Boolean_RegNet *A){
    longint_ARRAY *gl_keep_asis;
    
    gl_keep_asis=new longint_ARRAY();
    int idnow;
    
    if (!strcmp(Model_type,"Life_and_Death")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("GF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("GF_High");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    if (!strcmp(Model_type,"Life_and_Death_CIP")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("GF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("GF_High");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("ECM");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("CellDensity_high");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    
    
    if (!strcmp(Model_type,"Life_and_Death_Sen_UV_Gamma")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("GF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("GF_High");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Gamma");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("UV");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("ROS_ext");
        if(idnow>0) gl_keep_asis->add_element(idnow);

        //idnow=A->MDAT->get_node_ID("UV");
        //if(idnow>0) gl_keep_asis->add_element(idnow);
        //idnow=A->MDAT->get_node_ID("Onc_Ras");
        //if(idnow>0) gl_keep_asis->add_element(idnow);
        
    }
    
    if (!strcmp(Model_type,"Neuron")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("EGF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("EGF_High");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("NGF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("DLL1");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("NDD");
        if(idnow>0) gl_keep_asis->add_element(idnow);
     //   idnow=A->MDAT->get_node_ID("Trail");
     //   if(idnow>0) gl_keep_asis->add_element(idnow);
//        idnow=A->MDAT->get_node_ID("Neuron");
//        if(idnow>0) gl_keep_asis->add_element(idnow);
//        idnow=A->MDAT->get_node_ID("ROS_ext");
//        if(idnow>0) gl_keep_asis->add_element(idnow);

        //idnow=A->MDAT->get_node_ID("UV");
        //if(idnow>0) gl_keep_asis->add_element(idnow);
        //idnow=A->MDAT->get_node_ID("Onc_Ras");
        //if(idnow>0) gl_keep_asis->add_element(idnow);
        
    }
    
    
    if (!strcmp(Model_type,"Life_and_Death_UV_Gamma")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("GF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("GF_High");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("UV");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Gamma");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    if (!strcmp(Model_type,"Life_and_Death")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("GF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("GF_High");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
   
    
    if (!strcmp(Model_type,"Life_and_Death_Sen_UV_Gamma_ESC")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("GF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("GF_High");
        if(idnow>0) gl_keep_asis->add_element(idnow);
       // idnow=A->MDAT->get_node_ID("Gamma");
      //  if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("LIF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("2i");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Yamanaka");
        if(idnow>0) gl_keep_asis->add_element(idnow);

        /*idnow=A->MDAT->get_node_ID("2i");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Gata6_On");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Cdx2_On");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("BMP");
        if(idnow>0) gl_keep_asis->add_element(idnow);
         */
    }

    if (!strcmp(Model_type,"Life_and_Death_ESC_Reprogramming")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("GF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("GF_High");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Two_i");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("LIF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Gata6_On");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Cdx2_On");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("BMP");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Yamanaka");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    

    

    
    if (!strcmp(Model_type,"Senescence")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("GF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("GF_High");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("UV");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Gamma");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    if (!strcmp(Model_type,"Mito_only")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("Stress");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    return (gl_keep_asis);
}


longint_ARRAY *get_input_nodes_for_5input_drawing(Boolean_RegNet *A){
    longint_ARRAY *gl_keep_asis;
    
    gl_keep_asis=new longint_ARRAY();
    int idnow;
    
    if (!strcmp(Model_type,"Life_and_Death")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    if (!strcmp(Model_type,"Life_and_Death_CIP")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("ECM");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("CellDensity_high");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    if (!strcmp(Model_type,"Life_and_Death_Sen_UV_Gamma")) { // if this is the model you run now
      //  idnow=A->MDAT->get_node_ID("Onc_Ras");
      //  if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Gamma");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("UV");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("ROS_ext");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    if (!strcmp(Model_type,"Neuron")) { // if this is the model you run now
    //    idnow=A->MDAT->get_node_ID("Gamma");
    //    if(idnow>0) gl_keep_asis->add_element(idnow);
    //    idnow=A->MDAT->get_node_ID("Trail");
     //   if(idnow>0) gl_keep_asis->add_element(idnow);
//         idnow=A->MDAT->get_node_ID("Neuron");
//         if(idnow>0) gl_keep_asis->add_element(idnow);
         idnow=A->MDAT->get_node_ID("NGF");
         if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("DLL1");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("NDD");
        if(idnow>0) gl_keep_asis->add_element(idnow);
//         idnow=A->MDAT->get_node_ID("ROS_ext");
//         if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    if (!strcmp(Model_type,"Life_and_Death_UV_Gamma")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("UV");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Gamma");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    
    if (!strcmp(Model_type,"Life_and_Death_Sen_UV_Gamma_ESC")) { // if this is the model you run now
        //  idnow=A->MDAT->get_node_ID("Onc_Ras");
        //  if(idnow>0) gl_keep_asis->add_element(idnow);
      //  idnow=A->MDAT->get_node_ID("Gamma");
      //  if(idnow>0) gl_keep_asis->add_element(idnow);
        //  idnow=A->MDAT->get_node_ID("UV");
        //  if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("LIF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("2i");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Yamanaka");
        if(idnow>0) gl_keep_asis->add_element(idnow);
  }
    
    if (!strcmp(Model_type,"Life_and_Death_ESC")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("Two_i");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("LIF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Gata6_On");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Cdx2_On");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("BMP");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    
    if (!strcmp(Model_type,"Life_and_Death_ESC_Reprogramming")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("Two_i");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("LIF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
      //  idnow=A->MDAT->get_node_ID("Gata6_On");
       // if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Cdx2_On");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Yamanaka");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    if (!strcmp(Model_type,"Life_and_Death_ESC_Neuro")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("Two_i");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("LIF");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("BMP");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }

    if (!strcmp(Model_type,"Senescence")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("Trail");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("UV");
        if(idnow>0) gl_keep_asis->add_element(idnow);
        idnow=A->MDAT->get_node_ID("Gamma");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    if (!strcmp(Model_type,"Mito_only")) { // if this is the model you run now
        idnow=A->MDAT->get_node_ID("Stress");
        if(idnow>0) gl_keep_asis->add_element(idnow);
    }
    
    return (gl_keep_asis);
}

unsigned short int cell_is_Dead(Boolean_Dynamics_TRACKER *D,long int b){
    int idnow;
    idnow=D->BRN->MDAT->get_node_ID("CAD");
    if(idnow==0) return (0);
    for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow-1]==49) return(1);
    }
    return (0);
}

unsigned short int cell_is_Senescent(Boolean_Dynamics_TRACKER *D,long int b){
    int idnow1,idnow2,id3,id4,id5;
    int chs=0,ms=0;
    idnow1=D->BRN->MDAT->get_node_ID("SAHF");
    idnow2=D->BRN->MDAT->get_node_ID("MFN1_2");
    id3=D->BRN->MDAT->get_node_ID("Hyperfused");
    id4=D->BRN->MDAT->get_node_ID("Mito_Pot_Low");
    if(id4==0) id4=D->BRN->MDAT->get_node_ID("MP_Low");
    id5=D->BRN->MDAT->get_node_ID("AMPK_a1");
    if(id5==0) id5=D->BRN->MDAT->get_node_ID("AMPK");
    if((idnow1==0)&&(idnow2==0)) return (0);
    for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow1-1]==49) chs=1;
        if (
           // (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow2-1]==49) &&
            (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id3-1]==49) &&
            (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id4-1]==49) &&
            (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id5-1]==49)
            )            ms++;
    }
    if (ms!=D->Attractor_valleys[b]->N_a) ms=0; else ms=1;
    if((ms==0)&&(chs==0)) return (0);
    if((ms==1)&&(chs==0)) return (1);  // mitochondrial only
    if((ms==0)&&(chs==1)) return (2);   // chromosomal only
    return (3);                         // both
}

unsigned short int Double_DNA(Boolean_Dynamics_TRACKER *D,long int b){
    int idnow;
    
    idnow=D->BRN->MDAT->get_node_ID("Reduplicated_DNA");
    if(idnow!=0) {
        for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
            if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow-1]==49) return(2);
        }
    }
    
    idnow=D->BRN->MDAT->get_node_ID("4N_DNA");
    if(idnow==0) {
        idnow=D->BRN->MDAT->get_node_ID("f4N_DNA");
        if(idnow==0) return (0);
    }
    for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow-1]==49) return(1);
    }
    
    return (0);
}

unsigned short int DNA_damage(Boolean_Dynamics_TRACKER *D,long int b){
    // if it's always on : DSB = 1
    // if it's always on : Mild = 2
    // if it's always on : Severe = 3
    // if it's always on : both_m = 4
    // if it's always on : both_s = 5
    // intermittent 6
    int idnow_1,idnow_2,idnow_3;
    idnow_1=D->BRN->MDAT->get_node_ID("DSBF");
    idnow_2=D->BRN->MDAT->get_node_ID("Mild_DNA_Lesions");
    idnow_3=D->BRN->MDAT->get_node_ID("DNA_Lesions.csv");
    if((idnow_1==0)&&(idnow_2==0)&&(idnow_3==0)) return (0);
    
    int all_ds=1; int all_ss=1; int all_ss_bad=1;int alloff=1;
    for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow_1-1]==49) alloff=0;
        else all_ds = 0;
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow_2-1]==49) alloff=0;
        else all_ss = 0;
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow_3-1]==49) alloff=0;
        else all_ss_bad = 0;
    }
    if(all_ds==1){
        if ((all_ss==1)&&(all_ss_bad==0)) return(4);
        if ((all_ss==1)&&(all_ss_bad==1)) return(5);
        return (1);
    }
    if ((all_ss==1)&&(all_ss_bad==0)) return(2);
    if ((all_ss==1)&&(all_ss_bad==1)) return(3);
    if (alloff!=1) return (6);
    return (0);
}


unsigned short int non_CC_rhythm(Boolean_Dynamics_TRACKER *D,long int b){
    int idnow,M_RS,M_PHS,M_DNADam;
    if (D->Attractor_valleys[b]->N_a==1) return (0);
    
    idnow=D->BRN->MDAT->get_node_ID("CyclinB");
    M_PHS=D->BRN->MDAT->Module_ID[idnow];
    idnow=D->BRN->MDAT->get_node_ID("CyclinE");
    M_RS=D->BRN->MDAT->Module_ID[idnow];
    
    idnow=D->BRN->MDAT->get_node_ID("Ph_p53");
    M_DNADam=D->BRN->MDAT->Module_ID[idnow];

    int dnadam_only=1;
    int c=0;
    for(int s=1;s<=D->BRN->N;s++)
        if ((D->BRN->MDAT->Module_ID[s]!=M_RS)&&(D->BRN->MDAT->Module_ID[s]!=M_PHS)) {
            for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
                if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[s-1]!=
                    D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[s-1]
                    ) {
                    if (D->BRN->MDAT->Module_ID[s]!=M_DNADam) dnadam_only=0;
                    c=1;
                }
            }
        }
    if (c==1) {
        if(dnadam_only) return (2);
        else return (1);
    }
    return (0);
}

unsigned short int ESC_status(Boolean_Dynamics_TRACKER *D,long int b){
    int idnow;
    
    idnow=D->BRN->MDAT->get_node_ID("Nanog");
    int idnow1=D->BRN->MDAT->get_node_ID("Oct4");
    int idnow2=D->BRN->MDAT->get_node_ID("Oct4_Sox2");
    int idnow_NR=D->BRN->MDAT->get_node_ID("Neuron");
  //  if (idnow*idnow1*idnow2 + idnow_NR==0)  return (-1);
  
    int idnow3=D->BRN->MDAT->get_node_ID("Gata6");
    int idnow4=D->BRN->MDAT->get_node_ID("Cdx2");
    int idnow5=D->BRN->MDAT->get_node_ID("Zic1_3");
    
    for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
        if((idnow>0)&& (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow-1]==49)) return(0); // naive
    }
    int diff=1;
    
    int endod=0;int trop=0;int ect=0;int neuron = 0;
    for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
        if(idnow3>0) {
            if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow3-1]==49)
            endod=1;
        }
        if(idnow4>0) {
            if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow4-1]==49)
            trop=1;
        }
        if(idnow5>0) {
            if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow5-1]==49)
            ect=1;
        }
        if(idnow_NR>0) {
            if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow_NR-1]==49)
                neuron=1;
        }
    }
    
    for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow1-1]==49)
            diff=0;
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow2-1]==49)
            diff=0;
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow3-1]==48)
            endod=0;
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow4-1]==48)
            trop=0;
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow5-1]==48)
            ect=0;
    }
    
    if (diff)   {
        if (endod) return (2); // endoderm
        if (trop) return (3); // trophectoderm
        if (ect) return (4); // ectoderm
        if (neuron) return (5); // neuron
       
        int idnow1=D->BRN->MDAT->get_node_ID("Ecadherin_mRNA");
        int idnow2=D->BRN->MDAT->get_node_ID("miR_200");
        int idnow3=D->BRN->MDAT->get_node_ID("miR_34");
        int idnow4=D->BRN->MDAT->get_node_ID("ZEB1_H");
        int idnow5=D->BRN->MDAT->get_node_ID("Fast_Migration");

        if(idnow1*idnow2*idnow3*idnow4*idnow5){
            int mesench = 1; int mig = 1; int partial = 1;
            
            for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
                if ((D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow1-1]==49)
                    && (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow2-1]==49)
                    && (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow3-1]==49)){
                        mesench=0; partial = 0;
                }
                if ((D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow1-1]==49)
                    && (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow2-1]==49)
                    && (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow3-1]==48))
                    mesench=0;
                
                if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow5-1]==48)
                    mig=0;
            }
            if (mesench) return (7); // mesenchymal
            if (partial){
                if (mig) return (8);// partial, mig
                else return (8); // partial, nonmig
            }
            return (10); // epithelial
        }
        return (6); //weird undidentified differentiated state
    }
    else {
        return (1);// primed
    }
}

unsigned short int Migration_status(Boolean_Dynamics_TRACKER *D,long int b){
    int idnow;
    
    idnow=D->BRN->MDAT->get_node_ID("Fast_Migration");
    if(idnow>0){
        int allon=1, all_off=1;
        for(unsigned int i=1;i<=D->Attractor_valleys[b]->N_a;i++){
            if(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idnow-1]==48) allon = 0;
            else all_off =0;
        }
        if(allon) return (2);// always migrating
        if(all_off) return (0);// never migrating
        return (1);// sometimes migrating
    }
    return(0); // never migrating, no migr module
}

const char * celltype_name(unsigned int i){
    switch (i) {
        case 0:
            return ("Na");
            break;
        case 1:
            return ("Pr");
            break;
        case 2:
            return ("En");
            break;
        case 3:
            return ("Tr");
            break;
        case 4:
            return ("Ec");
            break;
        case 5:
            return ("Neuron");
            break;
        case 6:
            return ("Di");
            break;
        case 7:
            return ("M");
            break;
        case 8:
            return ("H_Mig");
            break;
        case 9:
            return ("H_NoMig");
            break;
        case 10:
            return ("E");
            break;
        default:
            return ("?");
            break;
    }
    return ("?");
}
unsigned short int cell_cycle_or_phase(Boolean_Dynamics_TRACKER *D,long int b){
    int id1,id2,id3,id4,id5,id6,g0,g2,sac;
    
    id1=D->BRN->MDAT->get_node_ID("Cdh1");
    id2=D->BRN->MDAT->get_node_ID("CyclinA");
    id3=D->BRN->MDAT->get_node_ID("UbcH10");
    id4=D->BRN->MDAT->get_node_ID("Cdk1");
    id5=D->BRN->MDAT->get_node_ID("Mad2");
    id6=D->BRN->MDAT->get_node_ID("pAPC");
    
    //printf("id6=%d\n",id6);getchar();
    
    //0 - fprintf(f,"Cell_Cycle\t");
    //1 - fprintf(f,"G0/G1\t");
    //2 - fprintf(f,"G2\t");
    //3 - fprintf(f,"SAC\t");
    //4 - fprintf(f,"Broken_Cycle\t");
    //5 - fprintf(f,"PHS_on_barrier\t");
    
    
    if (D->Attractor_valleys[b]->N_a==1) {
        g0= (1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id1-1]-49))*
            (1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id2-1]-48))*
            (1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id4-1]-48))*
            (1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id5-1]-48))*
            (1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id6-1]-48));
        
        g2=(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id2-1]-48)*
        (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id3-1]-48);
        
        sac=(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id4-1]-48)*
        (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id5-1]-48)*
        (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[id6-1]-48);
        
        if ((g0==1)&&(sac==0)&&(g2==0)) return (1); // G0/G1
        if ((g0==0)&&(sac==0)&&(g2==1)) return (2); //G2
        if ((g0==0)&&(sac==1)&&(g2==0)) return (3); //SAC
        return(5); // PHS_on_barrier
    }
    else{
        int rflip1=0;
        int idm=D->BRN->MDAT->get_node_ID("Replication");
        for (int i=2; i<=D->Attractor_valleys[b]->N_a; i++) {
            if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idm-1]!=
                D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i-1]]->s[idm-1])
                rflip1++;
        }
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[idm-1]!=
            D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[D->Attractor_valleys[b]->N_a]]->s[idm-1])
            rflip1++;
        if ((rflip1!=2)&&(rflip1!=0)) {
          //  printf("Broken_Cycle: ");
           // printf("b=%d\trflip1=%d\n",b,rflip1);//getchar();
            return (4);
        } // Broken_Cycle
        
        int rflip2=0;
        idm=D->BRN->MDAT->get_node_ID("4N_DNA");
        if(idm==0) idm=D->BRN->MDAT->get_node_ID("f4N_DNA");
        
        for (int i=2; i<=D->Attractor_valleys[b]->N_a; i++) {
            if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[idm-1]!=
                D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i-1]]->s[idm-1])
                rflip2++;
        }
        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[idm-1]!=
            D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[D->Attractor_valleys[b]->N_a]]->s[idm-1])
            rflip2++;
        if ((rflip2!=2)&&(rflip2!=0)) {
            printf("Broken_Cycle: ");
            printf("b=%ld\trflip1=%d\trflip1=%d\n",b,rflip1,rflip2);//getchar();
            return (4);
        } // Broken_Cycle
        if ((rflip1==2)&&(rflip2==2)) return (0); // Cell_Cycle
        else {
            g0=0;g2=0;sac=0;
            for (int i=1; i<=D->Attractor_valleys[b]->N_a; i++) {
                g0+=(1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id1-1]-49))*
                (1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id2-1]-48))*
                (1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id4-1]-48))*
                (1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id5-1]-48))*
                (1-(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id6-1]-48));
                
                g2+=(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id2-1]-48)*
                (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id3-1]-48);
                sac+=(D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id4-1]-48)*
                (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id5-1]-48)*
                (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[i]]->s[id6-1]-48);
            }
            // printf("b=%d\tg0=%d\tg2=%d\tsac=%d\n",b,g0,g2,sac); getchar();
            if ((g0==D->Attractor_valleys[b]->N_a)&&(sac==0)&&(g2==0)) return (1); // G0/G1
            if ((g0==0)&&(sac==0)&&(g2==D->Attractor_valleys[b]->N_a)) return (2); // G2
            if ((g0==0)&&(sac==D->Attractor_valleys[b]->N_a)&&(g2==0)) return (3); // SAC
            return(4); // Broken_Cycle
        }
    }
    return (-1);
}


Phenotype_now Grab_Phenotype_of_Attractor_Life_Death_ESC(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,long int b, int Coup){
    Phenotype_now PN;
    
    Boolean_Dynamics_TRACKER *D;
    if(Coup==1) { D=PD->Coupled; }
    else {D=PD->Unconnected;}

    if (cell_is_Dead(D,b)) PN.Dead=1;  else  PN.Dead=0;
    PN.Sen=cell_is_Senescent(D,b);
    PN.CC=cell_cycle_or_phase(D,b);
    PN.OtherC=non_CC_rhythm(D,b);
    PN.DNA_damage=DNA_damage(D,b);
    if (Double_DNA(D,b)==1) PN.G2_struck=1;
        else  {if (Double_DNA(D,b)==2) PN.G2_struck=2; else PN.G2_struck=0;}
    PN.ESC=ESC_status(D,b);
    PN.Migr=Migration_status(D,b);
    return (PN);
}

void write_attractor_line(FILE *f, Boolean_Dynamics_TRACKER *D, unsigned int b,longint_ARRAY *gl_inp){

    for(unsigned int k=1;k<=gl_inp->N;k++) {
        fprintf(f,"%d\t",D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[gl_inp->A[k]-1]-48);
    }
    if (cell_is_Dead(D,b)) fprintf(f,"YES\t"); else fprintf(f,"NO\t");
    
    int a= cell_is_Senescent(D,b);
    switch (a) {
        case 0:
            fprintf(f,"NO\t");
            break;
        case 1:
            fprintf(f,"Mitochondrial\t");
            break;
        case 2:
            fprintf(f,"SAHF\t");
            break;
        case 3:
            fprintf(f,"Mit_SAHF\t");
            break;
        default: fprintf(f,"NO\t");
            break;
    }
    
    a=cell_cycle_or_phase(D,b);
    switch (a) {
        case 0:
            fprintf(f,"Cell_Cycle\t");
            break;
        case 1:
            fprintf(f,"G0/G1\t");
            break;
        case 2:
            fprintf(f,"G2\t");
            break;
        case 3:
            fprintf(f,"SAC\t");
            break;
        case 4:
            fprintf(f,"Broken_Cycle\t");
            break;
        case 5:
            fprintf(f,"PHS_on_barrier\t");
            break;
        default: fprintf(f,"??\t");
            break;
    }
    
    if (non_CC_rhythm(D,b)) fprintf(f,"YES\t"); else fprintf(f,"NO\t");
    
    if (Double_DNA(D,b)==1) fprintf(f,"4N\t");
    else { if (Double_DNA(D,b)==2) fprintf(f,"Reduplicated\t");
            else  fprintf(f,"NO\t");fprintf(f,"NO\t");
    }
    
    a=DNA_damage(D,b);
    switch (a) {
        case 0:{fprintf(f,"No_DNA_damage\t");} break;
        case 1:{fprintf(f,"DS_breaks\t");} break;
        case 2:{fprintf(f,"Mild_ss_breaks\t");} break;
        case 3:{fprintf(f,"Severe_ss_breaks\t");} break;
        case 4:{fprintf(f,"DS_Mild_SS_breaks\t");} break;
        case 5:{fprintf(f,"DS_Severe_SS_breaks\t");} break;
        case 6:{fprintf(f,"Intermittent_damage\t");} break;
        default: break;
    }
    
    a=ESC_status(D,b);
    switch (a) {
        case 0:
            fprintf(f,"Naive"); break;
        case 1:
            fprintf(f,"Primed"); break;
        case 2:
            fprintf(f,"Endoderm"); break;
        case 3:
            fprintf(f,"Trophectoderm"); break;
        case 4:
            fprintf(f,"Ectoderm"); break;
        case 5:
            fprintf(f,"Neuron"); break;
        case 6:
            fprintf(f,"Differentiated, weird"); break;
        default: fprintf(f,"??");
            break;
    }
}

void Export_Attractor_Life_and_Death_table(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,longint_ARRAY *gl_inp,int Coup){
    char fn[300];
    FILE *f;
    Boolean_Dynamics_TRACKER *D;
    
    
    if(Coup==1) {
        D=PD->Coupled;
        sprintf(fn,"%s/%s/%s_Attractor_Profile",MODEL_Directory,D->DIRname,D->SAMPLEname);
        if(PD->BRN->silenced_inputs->N>0)
            for (unsigned k=1; k<=PD->BRN->silenced_inputs->N; k++) {
                sprintf(fn,"%s_%s",fn,PD->BRN->MDAT->Node_name[PD->BRN->silenced_inputs->A[k]]);
            }
        sprintf(fn, "%s.txt",fn);
    }
    else {D=PD->Unconnected;
        sprintf(fn,"%s/%s/%s_UNCONNECTED_Attractor_Profile.txt",MODEL_Directory,D->DIRname,D->SAMPLEname);
    }
    f=fopen(fn,"w");
    fprintf(f,"A_ID\t");
    for(unsigned int b=1;b<=gl_inp->N;b++) {
        fprintf(f,"%s\t",D->BRN->MDAT->Node_name[gl_inp->A[b]]);
    }
    if(D->Attr_Transitions==NULL) fprintf(f,"Dead\tSenescent\tCellCycle\tOtherCycle\t4N_DNA\tDNA_Damage\tESC\n");
    else fprintf(f,"Dead\tSenescent\tCellCycle\tOtherCycle\t4N_DNA\tDNA_Damage\tESC\t\t\tOther_ID\tMatch\n");
    // YES/NO; YES/NO; CC, G1, G2, SAC
    
    for(unsigned int b=1;b<=D->Attractor_NR;b++) {
        fprintf(f,"%d\t",b);
        write_attractor_line(f,D,b,gl_inp);
        if(D->Attr_Transitions==NULL) fprintf(f,"\n");
        else fprintf(f, "\t\t%lld\t%lg\n",D->Attr_Transitions[1]->node[b].ID,D->Attr_Transitions[1]->node[b].node_props->nw);
    }
    
    if(D->Attr_Transitions!=NULL){
        fprintf(f,"\n\n\n");
        snode::slinklist::slink *sc2;
        for(unsigned int ll=1;ll<=gl_inp->N;ll++){
            for(unsigned long int i=1;i<=D->Attr_Transitions[ll]->N;i++){
                sc2=D->Attr_Transitions[ll]->node[i].Li.first_link_pt;
                while(sc2!=NULL){
                    fprintf(f,"%lu\t%lld\t\t%d\t%d\t%lf\n",
                                 i,sc2->sneighbor,
                                    sc2->link_props->ID,
                                    sc2->link_props->distinc_link_labels,
                                    sc2->link_props->lw);
                    sc2=sc2->next(sc2);
                }
            }
        }
    }
    fclose(f);
}


void COMPARE_Attractor_Life_and_Death_table(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                            longint_ARRAY *gl_inp_A,
                                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_B,
                                            longint_ARRAY *gl_inp_B
                                            ){
    char fn[300];
    FILE *f;
    Boolean_Dynamics_TRACKER *D_A, *D_B;
    
    if(gl_inp_B->N!=gl_inp_A->N){
        printf("Number of environment-setting variables need to match!\n"); exit(1);
    }
    for(unsigned long int l=1;l<=gl_inp_B->N;l++)
        if(PD_B->BRN->MDAT->BIG_ID[gl_inp_B->A[l]]!= gl_inp_A->A[l]){
            printf("Order of environment-setting variables need to match!\n\t Mismatch at poz %ld: %s (index %ld) of full network vs.  %s (index %ld, bigID %d) of subset\n",l,PD_A->Coupled->BRN->MDAT->Node_name[gl_inp_A->A[l]],gl_inp_A->A[l],PD_B->Coupled->BRN->MDAT->Node_name[gl_inp_B->A[l]],gl_inp_B->A[l],PD_B->BRN->MDAT->BIG_ID[gl_inp_B->A[l]]); exit(1);
        }

    D_A=PD_A->Coupled;  if(D_A->Attr_Transitions==NULL) return;
    D_B=PD_B->Coupled;   if(D_B->Attr_Transitions==NULL) return;
    
    sprintf(fn,"%s/%s/COMP_%s_vs_%s_Attractor_Profile.txt",MODEL_Directory,D_A->DIRname,
                                                            D_A->SAMPLEname,D_B->SAMPLEname);
    
    f=fopen(fn,"w");
    fprintf(f,"COMPARISON\tA_ID\tSub_ID\t");
    for(unsigned int b=1;b<=gl_inp_A->N;b++) {
        fprintf(f,"%s\t",D_A->BRN->MDAT->Node_name[gl_inp_A->A[b]]);
    }
    fprintf(f,"Dead\tSenescent\tCellCycle\tOtherCycle\t4N_DNA\tDNA_Damage\tESC\tMatch\n");
    
    for(unsigned int b=1;b<=D_A->Attractor_NR;b++){
        if(D_A->Attr_Transitions[1]->node[b].ID!=0){
            if(D_A->Attr_Transitions[1]->node[b].node_props->nw > D_B->BRN->N-0.5)
                fprintf(f,"IDENTICAL\t%d\t%lld\t",b,D_A->Attr_Transitions[1]->node[b].ID);
            else {
                if((D_A->Attr_Transitions[1]->node[b].node_props->nw < D_B->BRN->N-0.5)
                 &&(D_A->Attr_Transitions[1]->node[b].node_props->nw > 0))
                    fprintf(f,"ALTERED\t%d\t%lld\t",b,D_A->Attr_Transitions[1]->node[b].ID);
                else {
                    if(  (D_A->Attr_Transitions[1]->node[b].node_props->nw < 0)
                       &&(D_A->Attr_Transitions[1]->node[b].node_props->nw > -99))
                        fprintf(f,"SAME_LENGHT_CYCLCE\t%d\t%lld\t",b,D_A->Attr_Transitions[1]->node[b].ID);
                    else {
                        if( (D_A->Attr_Transitions[1]->node[b].node_props->nw < -99)
                           &&(D_A->Attr_Transitions[1]->node[b].node_props->nw > -101))
                            fprintf(f,"ONE_DIFF_LENGHT_CYCLCE\t%d\t-100\t",b);
                        else {
                            if((D_A->Attr_Transitions[1]->node[b].node_props->nw < -101))
                                fprintf(f,"SEVERAL_DIFF_LENGHT_CYCLCES\t%d\t%lg\t",b,D_A->Attr_Transitions[1]->node[b].node_props->nw);
                            else fprintf(f,"??\t%d\t%lg\t",b,D_A->Attr_Transitions[1]->node[b].node_props->nw);
                        }
                    }
                }
            }
        }
        else fprintf(f,"LOST\t%d\t0\t",b);
        write_attractor_line(f,D_A,b,gl_inp_A);
        fprintf(f,"\n");
    }
    for(unsigned int b=1;b<=D_B->Attractor_NR;b++){
        if(D_B->Attr_Transitions[1]->node[b].ID==0){
            fprintf(f,"GAINED\t0\t%d\t",b);
            write_attractor_line(f,D_B,b,gl_inp_B);
            fprintf(f,"\n");
         }
    }
    
    /*
    snode::slinklist::slink *sc2;
    snode::slinklist::slink *sc3;
    
    for(unsigned long int l=1;l<=gl_inp_A->N;l++) {
         if(D_A->Attr_Transitions[l]!=NULL){
            fprintf(f,"\n\n\n");
            for(unsigned long int i=1;i<=D_A->Attr_Transitions[l]->N;i++){
                sc2=D_A->Attr_Transitions[l]->node[i].Li.first_link_pt;
                while(sc2!=NULL){
                    if(D_A->Attr_Transitions[l]->node[i].ID!=0){ // start maps to a reduced state
                        if(D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID!=0){ // end  maps to a reduced state
                            sc3=D_B->Attr_Transitions[l]->node[D_A->Attr_Transitions[l]->node[i].ID].Li.first_link_pt;  // scroll through links of the start map !
                            int nothingfound=1;
                            while(sc3!=NULL){
                                if (sc3->sneighbor == D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID){
                                    nothingfound=0;
                                       // there is one where both ends match!
                                    if (fabs(sc2->link_props->lw-sc3->link_props->lw)<0.1){
                                            // same type of transtion and/or same length
                                        fprintf(f,"PRESERVED\t%lu\t%llu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                i,D_A->Attr_Transitions[l]->node[i].ID,
                                                sc2->sneighbor,D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID,
                                                sc2->link_props->ID,sc3->link_props->ID,
                                                sc2->link_props->distinc_link_labels,
                                                sc2->link_props->lw,sc3->link_props->lw);
                                    }
                                    else {
                                       if((sc2->link_props->lw>-0.5)&&(sc2->link_props->lw<0.5)&& // 0
                                          ((sc3->link_props->lw<-0.5)||(sc3->link_props->lw>0.5))) // not 0
                                          fprintf(f,"Turned_Irrev\t%lu\t%llu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                   i,D_A->Attr_Transitions[l]->node[i].ID,
                                                   sc2->sneighbor,D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID,
                                                   sc2->link_props->ID,sc3->link_props->ID,
                                                   sc2->link_props->distinc_link_labels,
                                                   sc2->link_props->lw,sc3->link_props->lw);
                                        else if((sc3->link_props->lw>-0.5)&&(sc3->link_props->lw<0.5)&& //  0
                                           ((sc2->link_props->lw<-0.5)||(sc2->link_props->lw>0.5))) // not 0
                                           fprintf(f,"Turned_Rev\t%lu\t%llu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                    i,D_A->Attr_Transitions[l]->node[i].ID,
                                                    sc2->sneighbor,D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID,
                                                    sc2->link_props->ID,sc3->link_props->ID,
                                                    sc2->link_props->distinc_link_labels,
                                                    sc2->link_props->lw,sc3->link_props->lw);
                                        else if((sc2->link_props->lw>0.5)&&(sc3->link_props->lw>0.5))
                                            fprintf(f,"Altered_duration\t%lu\t%llu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                    i,D_A->Attr_Transitions[l]->node[i].ID,
                                                    sc2->sneighbor,D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID,
                                                    sc2->link_props->ID,sc3->link_props->ID,
                                                    sc2->link_props->distinc_link_labels,
                                                    sc2->link_props->lw,sc3->link_props->lw);
                                        else fprintf(f,"???_not_cool\t%lu\t%llu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                    i,D_A->Attr_Transitions[l]->node[i].ID,
                                                    sc2->sneighbor,D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID,
                                                    sc2->link_props->ID,sc3->link_props->ID,
                                                    sc2->link_props->distinc_link_labels,
                                                    sc2->link_props->lw,sc3->link_props->lw);
                                    }
                                }
                                sc2=sc2->next(sc2);
                            }
                            if(nothingfound==1){
                                // no match on other end!
                                fprintf(f,"NO_Matching_Transition_in_SM\t%lu\t%llu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t-\n",
                                        i,D_A->Attr_Transitions[l]->node[i].ID,
                                        sc2->sneighbor,D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID,
                                        sc2->link_props->ID,0,
                                        sc2->link_props->distinc_link_labels,
                                        sc2->link_props->lw);
                            }
                        }
                        else{
                            fprintf(f,"NO_Matching_END_state_in_SM\t%lu\t%llu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t-\n",
                                    i,D_A->Attr_Transitions[l]->node[i].ID,
                                    sc2->sneighbor,D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID,
                                    sc2->link_props->ID,0,
                                    sc2->link_props->distinc_link_labels,
                                    sc2->link_props->lw);
                        }
                    }
                    else {
                        fprintf(f,"NO_Matching_START_state_in_SM\t%lu\t%llu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t-\n",
                                    i,D_A->Attr_Transitions[l]->node[i].ID,
                                    sc2->sneighbor,D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID,
                                    sc2->link_props->ID,0,
                                    sc2->link_props->distinc_link_labels,
                                    sc2->link_props->lw);
                    }
                    sc2=sc2->next(sc2);
                }
            }
        }
     }
    
    for(unsigned long int l=1;l<=gl_inp_A->N;l++) {
        if(D_B->Attr_Transitions[l]!=NULL){
            for(unsigned long int i=1;i<=D_B->Attr_Transitions[l]->N;i++){
                sc2=D_B->Attr_Transitions[l]->node[i].Li.first_link_pt;
                while(sc2!=NULL){
                    if(D_B->Attr_Transitions[l]->node[i].ID!=0){ // start maps to a full state
                        
                        if(D_A->Attr_Transitions[l]->node[sc2->sneighbor].ID!=0){ // end  maps to a full state
                            sc3=D_A->Attr_Transitions[l]->node[D_B->Attr_Transitions[l]->node[i].ID].Li.first_link_pt;  // scroll through links of the start map !
                            int nothingfound=1;
                            while(sc3!=NULL){
                                if (sc3->sneighbor == D_B->Attr_Transitions[l]->node[sc2->sneighbor].ID){
                                    nothingfound=0;
                                    // there is one where both ends match!
                                    if (fabs(sc2->link_props->lw-sc3->link_props->lw)<0.1){
                                        // same type of transtion and/or same length
                                        fprintf(f,"PRESERVED\t%llu\t%lu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                D_B->Attr_Transitions[l]->node[i].ID,i,
                                                D_B->Attr_Transitions[l]->node[sc2->sneighbor].ID,sc2->sneighbor,
                                                sc3->link_props->ID,sc2->link_props->ID,
                                                sc2->link_props->distinc_link_labels,
                                                sc3->link_props->lw,sc2->link_props->lw);
                                    }
                                    else {
                                        if((sc2->link_props->lw>-0.5)&&(sc2->link_props->lw<0.5)&& // 0
                                           ((sc3->link_props->lw<-0.5)||(sc3->link_props->lw>0.5))) // not 0
                                            fprintf(f,"Turned_Rev\t%llu\t%lu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                    D_B->Attr_Transitions[l]->node[i].ID,i,
                                                    D_B->Attr_Transitions[l]->node[sc2->sneighbor].ID,sc2->sneighbor,
                                                    sc3->link_props->ID,sc2->link_props->ID,
                                                    sc2->link_props->distinc_link_labels,
                                                    sc3->link_props->lw,sc2->link_props->lw);
                                        else if((sc3->link_props->lw>-0.5)&&(sc3->link_props->lw<0.5)&& //  0
                                                ((sc2->link_props->lw<-0.5)||(sc2->link_props->lw>0.5))) // not 0
                                            fprintf(f,"Turned_Irrev\t%llu\t%lu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                    D_B->Attr_Transitions[l]->node[i].ID,i,
                                                    D_B->Attr_Transitions[l]->node[sc2->sneighbor].ID,sc2->sneighbor,
                                                    sc3->link_props->ID,sc2->link_props->ID,
                                                    sc2->link_props->distinc_link_labels,
                                                    sc3->link_props->lw,sc2->link_props->lw);
                                        else if((sc2->link_props->lw>0.5)&&(sc3->link_props->lw>0.5))
                                            fprintf(f,"Altered_duration\t%llu\t%lu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                    D_B->Attr_Transitions[l]->node[i].ID,i,
                                                    D_B->Attr_Transitions[l]->node[sc2->sneighbor].ID,sc2->sneighbor,
                                                    sc3->link_props->ID,sc2->link_props->ID,
                                                    sc2->link_props->distinc_link_labels,
                                                    sc3->link_props->lw,sc2->link_props->lw);
                                        else fprintf(f,"???_not_cool\t%llu\t%lu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t%lf\n",
                                                     D_B->Attr_Transitions[l]->node[i].ID,i,
                                                     D_B->Attr_Transitions[l]->node[sc2->sneighbor].ID,sc2->sneighbor,
                                                     sc3->link_props->ID,sc2->link_props->ID,
                                                     sc2->link_props->distinc_link_labels,
                                                     sc3->link_props->lw,sc2->link_props->lw);
                                    }
                                }
                                sc2=sc2->next(sc2);
                            }
                            if(nothingfound==1){
                                // no match on other end!
                                fprintf(f,"NO_Matching_Transition_in_FULL\t%llu\t%lu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t-\n",
                                        D_B->Attr_Transitions[l]->node[i].ID,i,
                                        D_B->Attr_Transitions[l]->node[sc2->sneighbor].ID,sc2->sneighbor,
                                        sc3->link_props->ID,sc2->link_props->ID,
                                        sc2->link_props->distinc_link_labels,
                                        sc2->link_props->lw);
                            }
                        }
                        else{
                            fprintf(f,"NO_Matching_END_state_in_FULL\t%llu\t%lu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t-\n",
                                    D_B->Attr_Transitions[l]->node[i].ID,i,
                                    D_B->Attr_Transitions[l]->node[sc2->sneighbor].ID,sc2->sneighbor,
                                    0,sc2->link_props->ID,
                                    sc2->link_props->distinc_link_labels,
                                    sc2->link_props->lw);
                        }
                    }
                    else {
                        fprintf(f,"NO_Matching_START_state_in_FULL\t%llu\t%lu\t%lld\t%lld\t\t%d\t%d\t%d\t%lf\t-\n",
                                D_B->Attr_Transitions[l]->node[i].ID,i,
                                D_B->Attr_Transitions[l]->node[sc2->sneighbor].ID,sc2->sneighbor,
                                0,sc2->link_props->ID,
                                sc2->link_props->distinc_link_labels,
                                sc2->link_props->lw);
                    }
                    sc2=sc2->next(sc2);
                }
            }
        }
    }
     */
    fclose(f);
}

longint_ARRAY *get_Healthy_G0_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        int ok=1;
        int idnow=PD->BRN->MDAT->get_node_ID(GF_low_name);
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==48) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID(GF_high_name);
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID("Onc_Ras");
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID("UV");
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID("Gamma");
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID("4N_DNA");
        if(idnow == 0)   idnow=PD->BRN->MDAT->get_node_ID("f4N_DNA");
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID("ROS_ext");
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID("Trail");
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if (at.Phenotype.Dead==1) ok=0;
        if (at.Phenotype.Sen>0) ok=0;
        if (at.Phenotype.DNA_damage>0) ok=0;
        if(ok==1) l->add_element(b);
    }
    return (l);
}
longint_ARRAY *get_Arrested_G0_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        int ok=1;
        int idnow=PD->BRN->MDAT->get_node_ID(GF_high_name);
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==48) ok=0; // first condition: has to occur in high GF env!
        
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if (at.Phenotype.Dead==1) ok=0;
        if (at.Phenotype.Sen>0) ok=0;
       // if (at.Phenotype.DNA_damage>0) ok=0;
        if(ok==1) l->add_element(b);
    }
    return (l);
}

longint_ARRAY *get_ALL_Senescent_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        int ok=1;
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if (at.Phenotype.Dead==1) ok=0;
        if (at.Phenotype.Sen==0) ok=0;
        if(ok==1) l->add_element(b);
    }
    return (l);
}

longint_ARRAY *get_Chromosomal_Senescent_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        int ok=1;
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if (at.Phenotype.Dead==1) ok=0;
        if (at.Phenotype.Sen!=2) ok=0;
        if(ok==1) l->add_element(b);
    }
    return (l);
}

longint_ARRAY *get_Mesenhymal_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if ((at.Phenotype.ESC==7)&&(at.Phenotype.Dead==0)) l->add_element(b);
    }
    return (l);
}

longint_ARRAY *get_Hybrid_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if ((at.Phenotype.ESC==8)&&(at.Phenotype.Dead==0)) l->add_element(b);
    }
    return (l);
}

longint_ARRAY *get_Epithelial_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if ((at.Phenotype.ESC==10)&&(at.Phenotype.Dead==0)) l->add_element(b);
    }
    return (l);
}


longint_ARRAY *get_Mitochondrial_Senescent_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        int ok=1;
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if (at.Phenotype.Dead==1) ok=0;
        if (at.Phenotype.Sen!=1) ok=0;
        if(ok==1) l->add_element(b);
    }
    return (l);
}

longint_ARRAY *get_Double_Senescent_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        int ok=1;
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if (at.Phenotype.Dead==1) ok=0;
        if (at.Phenotype.Sen!=3) ok=0;
        if(ok==1) l->add_element(b);
    }
    return (l);
}

longint_ARRAY *get_Healthy_CellCycle_attractors(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    Attractor_char_ALL at;
    longint_ARRAY *l;
    l=new longint_ARRAY();
    
    for(unsigned int b=1;b<=PD->Coupled->Attractor_NR;b++) {
        int ok=1;
        int idnow=PD->BRN->MDAT->get_node_ID(GF_low_name);
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==48) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID(GF_high_name);
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==48) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID("Onc_Ras");
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID("UV");
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        idnow=PD->BRN->MDAT->get_node_ID("Gamma");
        if(idnow>0) if(PD->Coupled->Attractor_valleys[b]->Basin_Stored[PD->Coupled->Attractor_valleys[b]->Attractor[1]]->s[idnow-1]==49) ok=0;
        
        at.Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        if (at.Phenotype.Dead==1) ok=0;
        if (at.Phenotype.Sen>0) ok=0;
        if (at.Phenotype.CC>0) ok=0;
        if(ok==1) l->add_element(b);
    }
    return (l);
}

longint_ARRAY *Find_Similar_Attractors__Life_and_Death(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,longint_ARRAY *gl_inp,int Coup){
    Boolean_Dynamics_TRACKER *D;
    Attractor_char_list * AL;
    FILE *f;
    char fn[300];
    longint_ARRAY *l_sim;
    l_sim=new longint_ARRAY();
    
    if(Coup==1) {
        D=PD->Coupled;
        sprintf(fn,"%s/%s/%s_Similar_attractors_Report",MODEL_Directory,D->DIRname,D->SAMPLEname);
        if(PD->BRN->silenced_inputs->N>0)
            for (unsigned k=1; k<=PD->BRN->silenced_inputs->N; k++) {
                sprintf(fn,"%s_%s",fn,PD->BRN->MDAT->Node_name[PD->BRN->silenced_inputs->A[k]]);
            }
        sprintf(fn,"%s.txt",fn);
    }
    else {D=PD->Unconnected;
        sprintf(fn,"%s/%s/%s_UNCONNECTED_Similar_attractors_Report",MODEL_Directory,D->DIRname,D->SAMPLEname);
        if(PD->BRN->silenced_inputs->N>0)
            for (unsigned k=1; k<=PD->BRN->silenced_inputs->N; k++) {
                sprintf(fn,"%s_%s",fn,PD->BRN->MDAT->Node_name[PD->BRN->silenced_inputs->A[k]]);
            }
        sprintf(fn,"%s.txt",fn);
    }
    f=fopen(fn,"w");
    
    AL=new Attractor_char_list[D->Attractor_NR+1];
    for(unsigned int b=1;b<=D->Attractor_NR;b++) {
        AL[b]=new Attractor_char_ALL;
        for(unsigned int k=1;k<=gl_inp->N;k++) {
            AL[b]->env_now[k]=D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[gl_inp->A[k]-1]-48;
            AL[b]->Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        }
    }
    
    for(unsigned int b=2;b<=D->Attractor_NR;b++)
        for(unsigned int b2=1;b2<b;b2++){
            int same=1;
            for(unsigned int k=1;k<=gl_inp->N;k++)
                if (AL[b]->env_now[k]!=AL[b2]->env_now[k])  same=0;
            if (same==1) {
                if (AL[b]->Phenotype.Dead!=AL[b2]->Phenotype.Dead) same=0;
                if (AL[b]->Phenotype.Sen!=AL[b2]->Phenotype.Sen) same=0;
                if (AL[b]->Phenotype.CC!=AL[b2]->Phenotype.CC) same=0;
                if (AL[b]->Phenotype.OtherC!=AL[b2]->Phenotype.OtherC) same=0;
                if (AL[b]->Phenotype.G2_struck!=AL[b2]->Phenotype.G2_struck) same=0;
                if (AL[b]->Phenotype.ESC!=AL[b2]->Phenotype.ESC) same=0;
                if (AL[b]->Phenotype.DNA_damage!=AL[b2]->Phenotype.DNA_damage) same=0;
            }
            if (same==1){
                printf("%u\t%u\t",b,b2);
              //  D->export_Attractor_onto_Boolean_network__PDF(b);
              //  D->export_Attractor_onto_Boolean_network__PDF(b2);
                l_sim->add_element(b);
                l_sim->add_element(b2);
                fprintf(f,"%u\t%u\t",b,b2);
                if (AL[b]->Phenotype.Dead==1) {printf("Dead; ");fprintf(f,"Dead; ");}
                if (AL[b]->Phenotype.Sen==1) {printf("Senescent; ");fprintf(f,"Senescent; ");}
                switch (AL[b]->Phenotype.CC) {
                    case 0:{printf("in_Cell_Cycle; ");fprintf(f,"in_Cell_Cycle; ");} break;
                    case 1:{printf("in_G0/G1; ");fprintf(f,"in_G0/G1; ");} break;
                    case 2:{printf("in_G2; ");fprintf(f,"in_G2; ");} break;
                    case 3:{printf("at_SAC; ");fprintf(f,"at_SAC; ");} break;
                    case 4:{printf("Broken_Cycle; ");fprintf(f,"Broken_Cycle; ");} break;
                    case 5:{printf("On_barrier; ");fprintf(f,"On_barrier; ");} break;
                    case 6:{printf("CC=??;");fprintf(f,"CC=??;");} break;
                    default: break;
                }
                switch (AL[b]->Phenotype.DNA_damage) {
                    case 0:{printf("No DNA damage; ");fprintf(f,"No_DNA_damage; ");} break;
                    case 1:{printf("DS_breaks; ");fprintf(f,"DS_breaks; ");} break;
                    case 2:{printf("Mild_ss_breaks; ");fprintf(f,"Mild_ss_breaks; ");} break;
                    case 3:{printf("Severe_ss_breaks; ");fprintf(f,"Severe_ss_breaks; ");} break;
                    case 4:{printf("DS_Mild_SS_breaks; ");fprintf(f,"DS_Mild_SS_breaks; ");} break;
                    case 5:{printf("DS_Severe_SS_breaks; ");fprintf(f,"DS_Severe_SS_breaks; ");} break;
                    case 6:{printf("Intermittent_damage;");fprintf(f,"Intermittent_damage;");} break;
                    default: break;
                }
                
                if (AL[b]->Phenotype.OtherC==1) printf("OtherCycle; ");
                if (AL[b]->Phenotype.OtherC==2) printf("DNADamage_cycle; ");

                if (AL[b]->Phenotype.ESC!=-1){
                    switch (AL[b]->Phenotype.ESC) {
                        case 0:{printf("Naive_ESC;");fprintf(f,"Naive_ESC;");} break;
                        case 1:{printf("Primed_ESC;");fprintf(f,"Primed_ESC;");} break;
                        case 2:{printf("Endoderm;");fprintf(f,"Endoderm;");} break;
                        case 3:{printf("Trophectoderm;");fprintf(f,"Trophectoderm;");} break;
                        case 4:{printf("Ectoderm;");fprintf(f,"Ectoderm;");} break;
                        case 5:{printf("Differentiated, weird;");fprintf(f,"Differentiated, weird;");} break;
                        default: break;
                    }
                }
                
                if ((D->Attractor_valleys[b]->N_a==1)&&(D->Attractor_valleys[b2]->N_a==1)) {
                    for (int i=1; i<=D->BRN->N; i++) {
                        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[i-1]!=D->Attractor_valleys[b2]->Basin_Stored[D->Attractor_valleys[b2]->Attractor[1]]->s[i-1]) {printf("\t%s",D->BRN->MDAT->Node_name[i]);
                                    fprintf(f,"\t%s",D->BRN->MDAT->Node_name[i]);
                        }
                    }
                }
                else {
                        int k1_best=0;int k2_best=0;int n_min=D->BRN->N;
                        for (int k1=1; k1<=D->Attractor_valleys[b]->N_a; k1++) {
                            for (int k2=1; k2<=D->Attractor_valleys[b2]->N_a; k2++){
                                int nn=0;
                                for (int i=1; i<=D->BRN->N; i++)
                                    if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[k1]]->s[i-1]!=D->Attractor_valleys[b2]->Basin_Stored[D->Attractor_valleys[b2]->Attractor[k2]]->s[i-1]) nn++;
                                if (nn<n_min) {
                                    n_min=nn;k1_best=k1;k2_best=k2;
                                }
                            }
                        }
                        for (int i=1; i<=D->BRN->N; i++) {
                            if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[k1_best]]->s[i-1]!=D->Attractor_valleys[b2]->Basin_Stored[D->Attractor_valleys[b2]->Attractor[k2_best]]->s[i-1]) {
                                    printf("\t%s",D->BRN->MDAT->Node_name[i]);
                                    fprintf(f,"\t%s",D->BRN->MDAT->Node_name[i]);
                            }
                        }
                }
                printf("\n");
                fprintf(f,"\n");//getchar();
            }
        }
    fclose(f);
    return (l_sim);
}

Attractor_char_list *grab_attractor_charact_Coupled(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                    longint_ARRAY *gl_inp){
    Attractor_char_list *AL;
    Boolean_Dynamics_TRACKER *D;
    
    D=PD->Coupled;
    AL=new Attractor_char_list[D->Attractor_NR+1];
    for(unsigned int b=1;b<=D->Attractor_NR;b++) {
        AL[b]=new Attractor_char_ALL;
        for(unsigned int k=1;k<=gl_inp->N;k++) {
            AL[b]->env_now[k]=D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[gl_inp->A[k]-1]-48;
            AL[b]->Phenotype=Grab_Phenotype_of_Attractor_Life_Death_ESC(PD,b,1);
        }
    }
    return (AL);
}

longint_ARRAY *Find_Similar_Attractors__Life_and_Death(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,longint_ARRAY *gl_inp,Attractor_char_list * AL){
   
    Boolean_Dynamics_TRACKER *D;
    FILE *f;
    char fn[300];
    longint_ARRAY *l_sim;
    l_sim=new longint_ARRAY();
    
    D=PD->Coupled;
    sprintf(fn,"%s/%s/%s_Similar_attractors_Report",MODEL_Directory,D->DIRname,D->SAMPLEname);
    if(PD->BRN->silenced_inputs->N>0)
        for (unsigned k=1; k<=PD->BRN->silenced_inputs->N; k++) {
            sprintf(fn,"%s_%s",fn,PD->BRN->MDAT->Node_name[PD->BRN->silenced_inputs->A[k]]);
        }
    sprintf(fn,"%s.txt",fn);

    f=fopen(fn,"w");
    
    for(unsigned int b=2;b<=D->Attractor_NR;b++)
        for(unsigned int b2=1;b2<b;b2++){
            int same=1;
            for(unsigned int k=1;k<=gl_inp->N;k++)
                if (AL[b]->env_now[k]!=AL[b2]->env_now[k])  same=0;
            if (same==1) {
                if (AL[b]->Phenotype.Dead!=AL[b2]->Phenotype.Dead) same=0;
                if (AL[b]->Phenotype.Sen!=AL[b2]->Phenotype.Sen) same=0;
                if (AL[b]->Phenotype.CC!=AL[b2]->Phenotype.CC) same=0;
                if (AL[b]->Phenotype.OtherC!=AL[b2]->Phenotype.OtherC) same=0;
                if (AL[b]->Phenotype.G2_struck!=AL[b2]->Phenotype.G2_struck) same=0;
                if (AL[b]->Phenotype.ESC!=AL[b2]->Phenotype.ESC) same=0;
                if (AL[b]->Phenotype.DNA_damage!=AL[b2]->Phenotype.DNA_damage) same=0;
            }
            if (same==1){
                printf("\t\t\t%u\t%u\t",b,b2);
                //  D->export_Attractor_onto_Boolean_network__PDF(b);
                //  D->export_Attractor_onto_Boolean_network__PDF(b2);
                l_sim->add_element(b);
                l_sim->add_element(b2);
                fprintf(f,"%u\t%u\t",b,b2);
                if (AL[b]->Phenotype.Dead==1) {printf("Dead; ");fprintf(f,"Dead; ");}
                if (AL[b]->Phenotype.Sen>0) {printf("Senescent; ");fprintf(f,"Senescent; ");}
                switch (AL[b]->Phenotype.CC) {
                    case 0:{printf("in_Cell_Cycle; ");fprintf(f,"in_Cell_Cycle; ");} break;
                    case 1:{printf("in_G0/G1; ");fprintf(f,"in_G0/G1; ");} break;
                    case 2:{printf("in_G2; ");fprintf(f,"in_G2; ");} break;
                    case 3:{printf("at_SAC; ");fprintf(f,"at_SAC; ");} break;
                    case 4:{printf("Broken_Cycle; ");fprintf(f,"Broken_Cycle; ");} break;
                    case 5:{printf("On_barrier; ");fprintf(f,"On_barrier; ");} break;
                    case 6:{printf("CC=??;");fprintf(f,"CC=??;");} break;
                    default: break;
                }
                switch (AL[b]->Phenotype.DNA_damage) {
                    case 0:{printf("No DNA damage; ");fprintf(f,"No_DNA_damage; ");} break;
                    case 1:{printf("DS_breaks; ");fprintf(f,"DS_breaks; ");} break;
                    case 2:{printf("Mild_ss_breaks; ");fprintf(f,"Mild_ss_breaks; ");} break;
                    case 3:{printf("Severe_ss_breaks; ");fprintf(f,"Severe_ss_breaks; ");} break;
                    case 4:{printf("DS_Mild_SS_breaks; ");fprintf(f,"DS_Mild_SS_breaks; ");} break;
                    case 5:{printf("DS_Severe_SS_breaks; ");fprintf(f,"DS_Severe_SS_breaks; ");} break;
                    case 6:{printf("Intermittent_damage;");fprintf(f,"Intermittent_damage;");} break;
                    default: break;
                }
                
                if (AL[b]->Phenotype.OtherC==1) printf("OtherCycle; ");
                if (AL[b]->Phenotype.OtherC==2) printf("DNADamage_cycle; ");
                
                if (AL[b]->Phenotype.ESC!=-1){
                    switch (AL[b]->Phenotype.ESC) {
                        case 0:{printf("Naive_ESC;");fprintf(f,"Naive_ESC;");} break;
                        case 1:{printf("Primed_ESC;");fprintf(f,"Primed_ESC;");} break;
                        case 2:{printf("Endoderm;");fprintf(f,"Endoderm;");} break;
                        case 3:{printf("Trophectoderm;");fprintf(f,"Trophectoderm;");} break;
                        case 4:{printf("Ectoderm;");fprintf(f,"Ectoderm;");} break;
                        case 5:{printf("Differentiated, weird;");fprintf(f,"Differentiated, weird;");} break;
                        default: break;
                    }
                }
                
                if ((D->Attractor_valleys[b]->N_a==1)&&(D->Attractor_valleys[b2]->N_a==1)) {
                    for (int i=1; i<=D->BRN->N; i++) {
                        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[1]]->s[i-1]!=D->Attractor_valleys[b2]->Basin_Stored[D->Attractor_valleys[b2]->Attractor[1]]->s[i-1]) {printf("\t%s",D->BRN->MDAT->Node_name[i]);
                            fprintf(f,"\t%s",D->BRN->MDAT->Node_name[i]);
                        }
                    }
                }
                else {
                    int k1_best=0;int k2_best=0;int n_min=D->BRN->N;
                    for (int k1=1; k1<=D->Attractor_valleys[b]->N_a; k1++) {
                        for (int k2=1; k2<=D->Attractor_valleys[b2]->N_a; k2++){
                            int nn=0;
                            for (int i=1; i<=D->BRN->N; i++)
                                if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[k1]]->s[i-1]!=D->Attractor_valleys[b2]->Basin_Stored[D->Attractor_valleys[b2]->Attractor[k2]]->s[i-1]) nn++;
                            if (nn<n_min) {
                                n_min=nn;k1_best=k1;k2_best=k2;
                            }
                        }
                    }
                    for (int i=1; i<=D->BRN->N; i++) {
                        if (D->Attractor_valleys[b]->Basin_Stored[D->Attractor_valleys[b]->Attractor[k1_best]]->s[i-1]!=D->Attractor_valleys[b2]->Basin_Stored[D->Attractor_valleys[b2]->Attractor[k2_best]]->s[i-1]) {
                            printf("\t%s",D->BRN->MDAT->Node_name[i]);
                            fprintf(f,"\t%s",D->BRN->MDAT->Node_name[i]);
                        }
                    }
                }
                printf("\n");
                fprintf(f,"\n");//getchar();
            }
        }
    fclose(f);
    return (l_sim);
}



#endif /* Model_specific_info_h */
