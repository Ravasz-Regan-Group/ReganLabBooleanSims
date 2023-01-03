//
//  Timecourses_and_stats_module.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 10/17/19.
//  Copyright © 2019 Regan_Group. All rights reserved.
//

#ifndef Timecourses_and_stats_module_h
#define Timecourses_and_stats_module_h

class One_Env_Mut_profile_of_Attr{
public:
    double *CC_KD,*Irrev_KD,*Rev_KD;
    double *CC_OE,*Irrev_OE,*Rev_OE;
    // CC: HealthyCC, Endoredupl_G2, Endoredupl_Aneup, Endoredupl_T;
    double Mark_Threshold=0.1; // marked as 0 below 10% difference
    
    One_Env_Mut_profile_of_Attr(int irr,int rpairs){
        CC_KD=new double[4];
        Irrev_KD = new double[irr];
        Rev_KD = new double[2*rpairs];
        CC_OE=new double[4];
        Irrev_OE = new double[irr];
        Rev_OE = new double[2*rpairs];
        for(int j=0;j<4;j++) {CC_KD[j]=0;CC_OE[j]=0;}
        for(int j=0;j<irr;j++) {Irrev_KD[j]=0;Irrev_OE[j]=0;}
        for(int j=0;j<2*rpairs;j++) {Rev_KD[j]=0;Rev_OE[j]=0;}
    }
    ~One_Env_Mut_profile_of_Attr(){
        delete[] CC_KD;CC_KD=NULL;delete[] CC_OE;CC_OE=NULL;
        delete[] Irrev_KD;Irrev_KD=NULL; delete[] Irrev_OE;Irrev_OE=NULL;
        delete Rev_KD;Rev_KD=NULL;delete Rev_OE;Rev_OE=NULL;
    }
};

typedef  One_Env_Mut_profile_of_Attr * One_Env_Mut_profile_of_Attr_1D;
typedef  One_Env_Mut_profile_of_Attr_1D * One_Env_Mut_profile_of_Attr_2D;
typedef  One_Env_Mut_profile_of_Attr_2D * One_Env_Mut_profile_of_Attr_3D;
typedef  One_Env_Mut_profile_of_Attr_3D * One_Env_Mut_profile_of_Attr_4D;
typedef  One_Env_Mut_profile_of_Attr_4D * One_Env_Mut_profile_of_Attr_5D;

class Mut_profile_of_Attr_5D{
public:
    One_Env_Mut_profile_of_Attr_5D *Mut_Profile;
    // CC: HealthyCC, Endoredupl_G2, Endoredupl_Aneup, Endoredupl_T;
    double Mark_Threshold=0.1; // marked as 0 below 10% difference
    unsigned int env_breakdown;
    
    Mut_profile_of_Attr_5D(int e,int irr,int rpairs){
        env_breakdown=e;
        Mut_Profile = new One_Env_Mut_profile_of_Attr_5D[2*env_breakdown];
        for(unsigned int i=0;i<env_breakdown*2;i++){    //GF fixed
            Mut_Profile[i]=new One_Env_Mut_profile_of_Attr_4D[2*env_breakdown];
            for(unsigned int j=0;j<env_breakdown*2;j++){ // GF neighbors fixed
                Mut_Profile[i][j]=new One_Env_Mut_profile_of_Attr_3D[2*env_breakdown];
                for(unsigned int k=0;k<env_breakdown*2;k++){ // GF neighbors ECM fixed
                    Mut_Profile[i][j][k]=new One_Env_Mut_profile_of_Attr_2D[env_breakdown];
                    for(unsigned int l=0;l<env_breakdown;l++){ // GF neighbors ECM fixed
                        Mut_Profile[i][j][k][l]=new One_Env_Mut_profile_of_Attr_1D[env_breakdown];
                        for(unsigned int m=0;m<env_breakdown;m++){ // GF neighbors ECM TGF fixed
                            Mut_Profile[i][j][k][l][m]=new One_Env_Mut_profile_of_Attr(irr,rpairs);
                        }
                    }
                }
            }
        }
    }
    ~Mut_profile_of_Attr_5D(){
        for(unsigned int i=0;i<env_breakdown*2;i++){    //GF fixed
            for(unsigned int j=0;j<env_breakdown*2;j++){ // GF neighbors fixed
                for(unsigned int k=0;k<env_breakdown*2;k++){ // GF neighbors ECM fixed
                    for(unsigned int l=0;l<env_breakdown;l++){ // GF neighbors ECM fixed
                        for(unsigned int m=0;m<env_breakdown;m++){ // GF neighbors ECM TGF fixed
                            delete Mut_Profile[i][j][k][l][m]; Mut_Profile[i][j][k][l][m]=NULL;
                        }
                        delete[] Mut_Profile[i][j][k][l]; Mut_Profile[i][j][k][l]=NULL;
                    }
                    delete[] Mut_Profile[i][j][k]; Mut_Profile[i][j][k]=NULL;
                }
                delete[] Mut_Profile[i][j];Mut_Profile[i][j]=NULL;
            }
            delete[] Mut_Profile[i];Mut_Profile[i]=NULL;
        }
       delete Mut_Profile;Mut_Profile=NULL;
    }
};

Environment_and_hits *extract_environment_and_hits(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, const char linenow[]){
    char attr_list[200][5000];
    char attr_list_condition[200][5000];
    Environment_and_hits *EnvHit;
    EnvHit=new Environment_and_hits((unsigned int)PD->Coupled->BRN->active_inputs->N);
    
    std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(linenow, '(',')');
    
    char bla2[300];unsigned long int Nf=0;
    if(word.size()>0) {
        sprintf(bla2,"%s",word[0].c_str());
        Nf=separate_line_by_string(bla2,attr_list,", ",200);
    }
    if (Nf > 0){
        for(unsigned int i=1;i<=Nf;i++){
            unsigned long int Nf2=separate_line_by_string(attr_list[i],attr_list_condition,"=",10);
            if (Nf2 > 0){
                istringstream iss_in(attr_list_condition[1]);
                string inputname;
                iss_in >> inputname;
                unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
                if (in_id>=1){
                    unsigned long int i_index=0;
                    i_index = PD->Coupled->BRN->active_inputs->check_for_element(in_id);
                    if (i_index==0){
                        printf("Environmental condition %s in line %s exists in the model, but it is either not an INPUT node, or it is currenly locked ON/OFF! \n\t(ignoring line; press enter to move to next instruction\n",inputname.c_str(),linenow); getchar(); return(NULL);
                    }
                    else {
                        istringstream iss_p(attr_list_condition[2]);
                        double input_p = -1;
                        iss_p >> input_p;
                        if ((iss_p.fail() == true) || (input_p<0) || (input_p>1)){
                            printf("Environmental condition %s for node %s in line %s cannot be read as a real number between 0 and 1! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[2], inputname.c_str(),linenow); getchar(); return(NULL);
                        }
                        else EnvHit->p_inputs[i_index] =input_p;
                    }
                }
                else {
                    printf("Environmental condition %s in line %s is not a valid node in this model! \n\t(ignoring line; press enter to move to next instruction\n",inputname.c_str(),linenow); getchar(); return(NULL);
                }
            }
            else {
                printf("Environmental conditions are not correcly specified in line %s\n\t Please list ALL active input variables and their desired levels between 0 and 1 as illustrated here for cell exposed to medium-high growth factor stimualtion in the Sizek etal PI3K paper: \n\t (GF=1, GF_High=0.5, Trail=0)\n\t(ignoring line; press enter to move to next instruction\n",linenow); getchar(); return(NULL);
            }
        }
    }
 //   else {printf("Environmental conditions are not correcly specified in line %s\n\t Please list ALL active input variables and their desired levels between 0 and 1 as illustrated here for cell exposed to medium-high growth factor stimualtion in the Sizek etal PI3K paper: \n\t (GF=1, GF_High=0.5, Trail=0)\n\t(ignoring line; press enter to move to next instruction\n",linenow); getchar();return(NULL);}
    
    for(int j=1;j<=PD->Coupled->BRN->active_inputs->N;j++){
        if(EnvHit->p_inputs[j]<0) {
            printf("Environmental conditions are not set for model input %s in line %s\n\t settting this input to 0. \n",PD->Coupled->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[j]],linenow);
            EnvHit->p_inputs[j]=0;
        }
    }
    EnvHit->p_inputs[0]=1;
    
    EnvHit->BKR_Hits=new longint_ARRAY_pair();
    EnvHit->Scan_Hits=new longint_ARRAY_pair();
    
    if(word.size()>=2){
        sprintf(bla2,"%s",word[1].c_str());
        Nf=separate_line_by_string(bla2,attr_list,", ",200);
        if (Nf > 0){
            // test which one:
            unsigned long int Nf2=separate_line_by_string(attr_list[1],attr_list_condition," ",10);
            if (Nf2 == 3) EnvHit->p_hits=new double[Nf+1];
            
            for(unsigned int i=1;i<=Nf;i++){
                unsigned long int Nf2=separate_line_by_string(attr_list[i],attr_list_condition," ",10);
                if (Nf2 == 3){
                    istringstream iss_in(attr_list_condition[1]);
                    string hitname;
                    iss_in >> hitname;
                    unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(hitname.c_str());
                    if (in_id>=1){
                            istringstream iss_p(attr_list_condition[3]);
                            double hit_p = -1;
                            iss_p >> hit_p;
                            if ((iss_p.fail() == true) || (hit_p<0) || (hit_p>1)){
                                printf("Knockout / overexpression condition %s for node %s in line %s cannot be read as a real number between 0 and 1! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[3], hitname.c_str(),linenow); getchar(); return(NULL);
                            }
                            else {
                                if(!strcmp(attr_list_condition[2],"KD")) {EnvHit->BKR_Hits->add_element(in_id, 0);EnvHit->p_hits[i]=hit_p;}
                                else if(!strcmp(attr_list_condition[2],"OE")) {EnvHit->BKR_Hits->add_element(in_id, 1);EnvHit->p_hits[i]=hit_p;}
                                else{
                                    printf("Knockout / overexpression condition %s for node %s in line %s is neigher KD for knockdown or OE for overexpression! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[2], hitname.c_str(),linenow); getchar(); return(NULL);
                                }
                            }
                        }
                    else {
                        printf("Knockout / overexpression condition %s in line %s is not a valid node in this model! \n\t(ignoring line; press enter to move to next instruction\n",hitname.c_str(),linenow); getchar(); return(NULL);
                    }
                }
                else {
                    if (Nf2 == 2){
                        istringstream iss_in(attr_list_condition[1]);
                        string hitname;
                        iss_in >> hitname;
                        unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(hitname.c_str());
                        if (in_id>=1){
                            if(!strcmp(attr_list_condition[2],"KD")) EnvHit->Scan_Hits->add_element(in_id, 0);
                                else if(!strcmp(attr_list_condition[2],"OE")) EnvHit->Scan_Hits->add_element(in_id, 1);
                                else{
                                    printf("Knockout / overexpression scan condition %s for node %s in line %s is neigher KD for knockdown or OE for overexpression! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[2], hitname.c_str(),linenow); getchar(); return(NULL);
                                }
                            }
                        else {
                            printf("Knockout / overexpression condition %s in line %s is not a valid node in this model! \n\t(ignoring line; press enter to move to next instruction\n",hitname.c_str(),linenow); getchar(); return(NULL);
                        }
                    }
                    else{
                        printf("Knockout / overexpression conditions are not correcly specified in line %s\n\t Please specify one or both of the following:\n\t nodes that are partially knocked down or overexpressed at specific desired levels between 0 and 1:\t (SNAI1 KD 0.5, ZEB1 OE 0.2)\n\t list of nodes that are knocked down/overexpressed in tandem at levels scanned from 0 to 1:\t(Merlin KD, YAP OE) \n\t (ignoring line; press enter to move to next instruction\n",linenow); getchar(); return(NULL);
                    }
                }
            }
        }
    }
    else{
        EnvHit->p_hits=new double[1];
        EnvHit->p_hits[0]=0;
    }
    
    if(word.size()==3){
        sprintf(bla2,"%s",word[2].c_str());
        Nf=separate_line_by_string(bla2,attr_list,", ",200);
        if (Nf > 0){
            unsigned long int Nf2=separate_line_by_string(attr_list[1],attr_list_condition," ",10);
            if (Nf2 == 3) EnvHit->p_hits=new double[Nf+1];
            for(unsigned int i=1;i<=Nf;i++){
                unsigned long int Nf2=separate_line_by_string(attr_list[i],attr_list_condition," ",10);
                if ((Nf2 == 3)&&(EnvHit->BKR_Hits->N==0)){
                    istringstream iss_in(attr_list_condition[1]);
                    string hitname;
                    iss_in >> hitname;
                    unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(hitname.c_str());
                    if (in_id>=1){
                        istringstream iss_p(attr_list_condition[3]);
                        double hit_p = -1;
                        iss_p >> hit_p;
                        if ((iss_p.fail() == true) || (hit_p<0) || (hit_p>1)){
                            printf("Knockout / overexpression condition %s for node %s in line %s cannot be read as a real number between 0 and 1! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[3], hitname.c_str(),linenow); getchar(); return(NULL);
                        }
                        else {
                            if(!strcmp(attr_list_condition[2],"KD")) {EnvHit->BKR_Hits->add_element(in_id, 0);EnvHit->p_hits[i]=hit_p;}
                            else if(!strcmp(attr_list_condition[2],"OE")) {EnvHit->BKR_Hits->add_element(in_id, 1);EnvHit->p_hits[i]=hit_p;}
                            else{
                                printf("Knockout / overexpression condition %s for node %s in line %s is neigher KD for knockdown or OE for overexpression! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[2], hitname.c_str(),linenow); getchar(); return(NULL);
                            }
                        }
                    }
                    else {
                        printf("Knockout / overexpression condition %s in line %s is not a valid node in this model! \n\t(ignoring line; press enter to move to next instruction\n",hitname.c_str(),linenow); getchar(); return(NULL);
                    }
                }
                else {
                    if ((Nf2 == 2)&&(EnvHit->Scan_Hits->N==0)){
                        istringstream iss_in(attr_list_condition[1]);
                        string hitname;
                        iss_in >> hitname;
                        unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(hitname.c_str());
                        if (in_id>=1){
                            if(!strcmp(attr_list_condition[2],"KD")) EnvHit->Scan_Hits->add_element(in_id, 0);
                            else if(!strcmp(attr_list_condition[2],"OE")) EnvHit->Scan_Hits->add_element(in_id, 1);
                            else{
                                printf("Knockout / overexpression scan condition %s for node %s in line %s is neigher KD for knockdown or OE for overexpression! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[2], hitname.c_str(),linenow); getchar(); return(NULL);
                            }
                        }
                        else {
                            printf("Knockout / overexpression condition %s in line %s is not a valid node in this model! \n\t(ignoring line; press enter to move to next instruction\n",hitname.c_str(),linenow); getchar(); return(NULL);
                        }
                    }
                    else{
                        printf("Knockout / overexpression conditions are not correcly specified in line %s\n\t Please specify one or both of the following:\n\t nodes that are partially knocked down or overexpressed at specific desired levels between 0 and 1:\t (SNAI1 KD 0.5, ZEB1 OE 0.2)\n\t list of nodes that are knocked down/overexpressed in tandem at levels scanned from 0 to 1:\t(Merlin KD, YAP OE) \n\t (ignoring line; press enter to move to next instruction\n",linenow); getchar(); return(NULL);
                    }
                }
            }
        }
    }
    else{
        if(word.size()>3){
            printf("Too many knockout / overexpression conditions in line %s\n\t Please specify one or both of the following:\n\t nodes that are partially knocked down or overexpressed at specific desired levels between 0 and 1:\t (SNAI1 KD 0.5, ZEB1 OE 0.2)\n\t list of nodes that are knocked down/overexpressed in tandem at levels scanned from 0 to 1:\t(Merlin KD, YAP OE) \n\t (ignoring line; press enter to move to next instruction\n",linenow); getchar(); return(NULL);
        }
    }
    return (EnvHit);
}

Environment_and_hits *extract_Hits(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD, unsigned long int a_st, const char linenow[]){
    char attr_list[200][5000];
    char attr_list_condition[200][5000];
    Environment_and_hits *EnvHit;
    EnvHit=new Environment_and_hits((unsigned int)PD->Coupled->BRN->active_inputs->N);
    
    std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(linenow, '(',')');
    
    char bla2[300];unsigned long int Nf=0;
    
    for(int j=1;j<=PD->Coupled->BRN->active_inputs->N;j++){
        EnvHit->p_inputs[j]=PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[1]]->s[PD->Coupled->BRN->active_inputs->A[j]-1];
    }
    EnvHit->p_inputs[0]=1;
    
    EnvHit->BKR_Hits=new longint_ARRAY_pair();
    EnvHit->Scan_Hits=new longint_ARRAY_pair();
    
    if(word.size()>=1){
        sprintf(bla2,"%s",word[0].c_str());
        Nf=separate_line_by_string(bla2,attr_list,", ",200);
        if (Nf > 0){
            // test which one:
            unsigned long int Nf2=separate_line_by_string(attr_list[1],attr_list_condition," ",10);
            for(unsigned int i=1;i<=Nf;i++){
                unsigned long int Nf2=separate_line_by_string(attr_list[i],attr_list_condition," ",10);
                if (Nf2 == 2){
                        istringstream iss_in(attr_list_condition[1]);
                        string hitname;
                        iss_in >> hitname;
                        unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(hitname.c_str());
                        if (in_id>=1){
                            if(!strcmp(attr_list_condition[2],"KD")) EnvHit->Scan_Hits->add_element(in_id, 0);
                                else if(!strcmp(attr_list_condition[2],"OE")) EnvHit->Scan_Hits->add_element(in_id, 1);
                                else{
                                    printf("Knockout / overexpression scan condition %s for node %s in line %s is neigher KD for knockdown or OE for overexpression! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[2], hitname.c_str(),linenow); getchar(); return(NULL);
                                }
                            }
                        else {
                            printf("Knockout / overexpression condition %s in line %s is not a valid node in this model! \n\t(ignoring line; press enter to move to next instruction\n",hitname.c_str(),linenow); getchar(); return(NULL);
                        }
                }
                else{
                    printf("Knockout / overexpression conditions are not correcly specified in line %s\n\t Please specify the list of nodes that are knocked down/overexpressed in tandem at different points in time along a network state or trajectory:\t(Merlin KD, YAP OE) \n\t (ignoring line; press enter to move to next instruction\n",linenow); getchar(); return(NULL);
                    }
            }
        }
    }
    else{
        EnvHit->p_hits=new double[1];
        EnvHit->p_hits[0]=0;
    }
    
    if(word.size()==3){
        sprintf(bla2,"%s",word[2].c_str());
        Nf=separate_line_by_string(bla2,attr_list,", ",200);
        if (Nf > 0){
            unsigned long int Nf2=separate_line_by_string(attr_list[1],attr_list_condition," ",10);
            if (Nf2 == 3) EnvHit->p_hits=new double[Nf+1];
            for(unsigned int i=1;i<=Nf;i++){
                unsigned long int Nf2=separate_line_by_string(attr_list[i],attr_list_condition," ",10);
                if ((Nf2 == 3)&&(EnvHit->BKR_Hits->N==0)){
                    istringstream iss_in(attr_list_condition[1]);
                    string hitname;
                    iss_in >> hitname;
                    unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(hitname.c_str());
                    if (in_id>=1){
                        istringstream iss_p(attr_list_condition[3]);
                        double hit_p = -1;
                        iss_p >> hit_p;
                        if ((iss_p.fail() == true) || (hit_p<0) || (hit_p>1)){
                            printf("Knockout / overexpression condition %s for node %s in line %s cannot be read as a real number between 0 and 1! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[3], hitname.c_str(),linenow); getchar(); return(NULL);
                        }
                        else {
                            if(!strcmp(attr_list_condition[2],"KD")) {EnvHit->BKR_Hits->add_element(in_id, 0);EnvHit->p_hits[i]=hit_p;}
                            else if(!strcmp(attr_list_condition[2],"OE")) {EnvHit->BKR_Hits->add_element(in_id, 1);EnvHit->p_hits[i]=hit_p;}
                            else{
                                printf("Knockout / overexpression condition %s for node %s in line %s is neigher KD for knockdown or OE for overexpression! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[2], hitname.c_str(),linenow); getchar(); return(NULL);
                            }
                        }
                    }
                    else {
                        printf("Knockout / overexpression condition %s in line %s is not a valid node in this model! \n\t(ignoring line; press enter to move to next instruction\n",hitname.c_str(),linenow); getchar(); return(NULL);
                    }
                }
                else {
                    if ((Nf2 == 2)&&(EnvHit->Scan_Hits->N==0)){
                        istringstream iss_in(attr_list_condition[1]);
                        string hitname;
                        iss_in >> hitname;
                        unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(hitname.c_str());
                        if (in_id>=1){
                            if(!strcmp(attr_list_condition[2],"KD")) EnvHit->Scan_Hits->add_element(in_id, 0);
                            else if(!strcmp(attr_list_condition[2],"OE")) EnvHit->Scan_Hits->add_element(in_id, 1);
                            else{
                                printf("Knockout / overexpression scan condition %s for node %s in line %s is neigher KD for knockdown or OE for overexpression! \n\t(ignoring line; press enter to move to next instruction\n",attr_list_condition[2], hitname.c_str(),linenow); getchar(); return(NULL);
                            }
                        }
                        else {
                            printf("Knockout / overexpression condition %s in line %s is not a valid node in this model! \n\t(ignoring line; press enter to move to next instruction\n",hitname.c_str(),linenow); getchar(); return(NULL);
                        }
                    }
                    else{
                        printf("Knockout / overexpression conditions are not correcly specified in line %s\n\t Please specify one or both of the following:\n\t nodes that are partially knocked down or overexpressed at specific desired levels between 0 and 1:\t (SNAI1 KD 0.5, ZEB1 OE 0.2)\n\t list of nodes that are knocked down/overexpressed in tandem at levels scanned from 0 to 1:\t(Merlin KD, YAP OE) \n\t (ignoring line; press enter to move to next instruction\n",linenow); getchar(); return(NULL);
                    }
                }
            }
        }
    }
    else{
        if(word.size()>3){
            printf("Too many knockout / overexpression conditions in line %s\n\t Please specify one or both of the following:\n\t nodes that are partially knocked down or overexpressed at specific desired levels between 0 and 1:\t (SNAI1 KD 0.5, ZEB1 OE 0.2)\n\t list of nodes that are knocked down/overexpressed in tandem at levels scanned from 0 to 1:\t(Merlin KD, YAP OE) \n\t (ignoring line; press enter to move to next instruction\n",linenow); getchar(); return(NULL);
        }
    }
    return (EnvHit);
}

longint_ARRAY_pair *extract_sequence_of_states(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,string inparantheses,bool *cyc){
    char attr_list[200][5000];
    char attr_list_condition[200][5000];
    longint_ARRAY_pair *statepairs;
    *cyc =true;

    statepairs=new longint_ARRAY_pair();
    char bla2[300]; sprintf(bla2,"%s",inparantheses.c_str());
    
    unsigned long int Nf=separate_line_by_string(bla2,attr_list,", ",200);
    if (Nf == 0)  Nf=separate_line_by_string(bla2,attr_list," -> ",200); else *cyc =false;
    if (Nf > 0){
        for(unsigned int i=1;i<=Nf;i++){
            unsigned long int Nf2=separate_line_by_string(attr_list[i],attr_list_condition,"=",10);
            if (Nf2 > 0){
                istringstream iss_in(attr_list_condition[1]);
                string inputname;
                iss_in >> inputname;
                unsigned int in_id = PD->Coupled->BRN->MDAT->get_node_ID(inputname.c_str());
                if (in_id>=1){
                    istringstream iss_p(attr_list_condition[2]);
                    unsigned int nstate = -1;
                    iss_p >> nstate;
                    if ((iss_p.fail() == true) || (nstate<0) || (nstate>1)){
                        printf("Node state %s for node %s in group of nodestates %s cannot be read as an integer of 0 and 1! \n\t(ignoring group of conditions; press enter to move to next instruction\n",attr_list_condition[2], inputname.c_str(),inparantheses.c_str()); getchar();
                        delete statepairs; statepairs=new longint_ARRAY_pair();
                        return(statepairs);
                    }
                    else statepairs->add_element(in_id, nstate);
                }
                else {
                    printf("Node %s in line %s is not a valid node in this model! \n\t(ignoring group of conditions; press enter to move to next instruction\n",inputname.c_str(),inparantheses.c_str()); getchar();
                    delete statepairs; statepairs=new longint_ARRAY_pair();
                    return(statepairs);
                }
            }
            else {
                printf("Node states are not correcly specified in: %s\n\t Please list ALL nodes required to identify a phenotype as:\n\t a stable list for fixed points, e.g.: (SNAI1=0, Twist=0, miR_34=1, Ecadherin_mRNA=1)\n\t a list of node-states that must follow each other in an ordered cylce, e.g.: (p110_H=1 -> AKT_H=1 -> Nedd4L=1 -> p110_H=0 -> FoxO3=1)\n\t\t(ignoring line; press enter to move to next instruction\n",inparantheses.c_str()); getchar();return(statepairs);
            }
        }
    }
    else {
        printf("Node states are not correcly specified in: %s\n\t Please list ALL nodes required to identify a phenotype as:\n\t a stable list for fixed points, e.g.: (SNAI1=0, Twist=0, miR_34=1, Ecadherin_mRNA=1)\n\t a list of node-states that must follow each other in an ordered cylce, e.g.: (p110_H=1 -> AKT_H=1 -> Nedd4L=1 -> p110_H=0 -> FoxO3=1)\n\t\t(ignoring line; press enter to move to next instruction\n",inparantheses.c_str());getchar();
        return(statepairs);
    }
    return (statepairs);
}

void load_Phenotypes_for_model(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD){
    char fn[300];
    FILE *f;
    longint_ARRAY_pair *a,aa;
    longint_ARRAY_pair *b,bb;
    
    if(Phen_STATs.WT_cell_cycle>0) return;
    
    string CC_block = "CellCycleError_Points:";
    string I_CC_block = "Irreversible_Do_Not_Count_CellCycle_Errors:";
    string I_STAT_block = "Irreversible_STATS:";
    string R_STAT_block = "Reversible_Transitions:";
    
    sprintf(fn,"%s/%s/VirtualExp_Phenotypes_for_STATS.txt",MODEL_Directory, PD->Coupled->SAMPLEname);
    f=fopen(fn,"r");
    if(f==NULL) {printf("VirtualExp_Phenotypes_for_STATS.txt File is not found in Model directory %s \n \tCurrent model folder: %s\n%s\n",PD->Coupled->SAMPLEname,MODEL_Directory_system,PD->Coupled->SAMPLEname); exit(1);}
    
    while(fgets(GSE_line,MAX_LINE,f)!=NULL){
        istringstream iss(GSE_line);
        string categoryname;
        iss >> categoryname;
        if (categoryname == CC_block) {
            unsigned int ccl = 0;
            iss >> ccl;
            if ((iss.fail() == true) || (ccl<0)){
                printf("In line %s please specify length of wild-type cell cycle under saturating growth factor conditions\n\t e.g.:  CellCycleError_Points:    21\n\t\t -- if more than one cycle exist, please choose the fastes one you wish to use for reference when calculating the relative rate of other phenotype transitions\n",GSE_line); getchar();
            }
            else {Phen_STATs.WT_cell_cycle = ccl; Phen_STATs.CellCycleErrors=true;}
            for(unsigned int i=1;i<=5;i++)
                if(fgets(GSE_line,MAX_LINE,f)!=NULL){
                    istringstream iss2(GSE_line);
                    string cc,nn,bla;
                    iss2 >> cc; iss2 >> bla; iss2 >> nn;
                    if(cc == "DNA_Replication_START"){
                        Phen_STATs.DNA_Replication_START = PD->Coupled->BRN->MDAT->get_node_ID(nn.c_str());
                        if(Phen_STATs.DNA_Replication_START == 0) {
                            printf("Node %s, listed as a marker of %s is not in the model!\n",nn.c_str(),cc.c_str());
                            getchar();
                        }
                    }
                    if(cc == "4N_DNA"){
                        Phen_STATs.DNA_Replication_DONE = PD->Coupled->BRN->MDAT->get_node_ID(nn.c_str());
                        if(Phen_STATs.DNA_Replication_DONE == 0) {
                            printf("Node %s, listed as a marker of %s is not in the model!\n",nn.c_str(),cc.c_str());
                            getchar();
                        }
                    }
                    if(cc == "Metaphase_START"){
                        Phen_STATs.Metaphase_START = PD->Coupled->BRN->MDAT->get_node_ID(nn.c_str());
                        if(Phen_STATs.Metaphase_START == 0) {
                            printf("Node %s, listed as a marker of %s is not in the model!\n",nn.c_str(),cc.c_str());
                            getchar();
                        }
                    }
                    if(cc == "SAC_Cleared"){
                        Phen_STATs.SAC_Cleared = PD->Coupled->BRN->MDAT->get_node_ID(nn.c_str());
                        if(Phen_STATs.SAC_Cleared == 0) {
                            printf("Node %s, listed as a marker of %s is not in the model!\n",nn.c_str(),cc.c_str());
                            getchar();
                        }
                    }
                    if(cc == "Cytokinesis"){
                        Phen_STATs.Cytokinesis = PD->Coupled->BRN->MDAT->get_node_ID(nn.c_str());
                        if(Phen_STATs.Cytokinesis == 0) {
                            printf("Node %s, listed as a marker of %s is not in the model!\n",nn.c_str(),cc.c_str());
                            getchar();
                        }
                    }
                }
        }
        
        if (categoryname == I_CC_block) {
            while ((fgets(GSE_line,MAX_LINE,f)!=NULL)&&((GSE_line[0]=='\t')||(GSE_line[0]=='-'))) {
                istringstream iss2(GSE_line);
                string cc;
                iss2 >> cc;
                if(cc != "--"){
                    Phen_STATs.Irreversible_Do_Not_Count_CellCycle_Errors.push_back(cc); // marks each irreversible phen that should stop cell cycle counts
                    longint_ARRAY_pair *a;
                    bool c = false;
                    std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(GSE_line, '(',')');
                    a=extract_sequence_of_states(PD,word[0], &c); // this is a list of node-state pairs belonging to a single phneotype
                    if (a->N>0){
                         Phen_STATs.Irreversible_Do_Not_Count_CellCycle_Errors_NodeStates.push_back(*a); // the whole package is added as a vector element, there is one for each phenotye
                         Phen_STATs.Irreversible_Do_Not_Count_CellCycle_Errors_LimitCycle.push_back(c);
                    }
                    else {printf("Incorrectly formatted  phenotype string in the Irreversible_Do_Not_Count_CellCycle_Errors: block %s\n\t code requires one phenotype name, then a single = character bracketed by two spaces, then a list of stable states or cycles; e.g.:\n\tEpithelial = (SNAI1=0, Twist=0, miR_34=1, Ecadherin_mRNA=1)\n\tPI3KOscillation = (p110_H=1 -> AKT_H=1 -> Nedd4L=1 -> p110_H=0 -> FoxO3=1)\t skipping line, press ENTER to continue.\n",GSE_line);
                        getchar ();}
                    
                }
            }
        }
        
        if (categoryname == I_STAT_block) {
            while ((fgets(GSE_line,MAX_LINE,f)!=NULL)&&((GSE_line[0]=='\t')||(GSE_line[0]=='-'))) {
                istringstream iss2(GSE_line);
                string cc;
                iss2 >> cc;
                if(cc != "--"){
                    Phen_STATs.Irreversible_STATS.push_back(cc);
                    bool c = false;
                    std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(GSE_line, '(',')');
                    a=extract_sequence_of_states(PD,word[0], &c);
                    if(a->N>0){
                        Phen_STATs.Irreversible_STATS_NodeStates.push_back(*a);
                        Phen_STATs.Irreversible_STATS_LimitCycle.push_back(c);
                        Phen_STATs.Irrev =true;
                    }
                    else {printf("Incorrectly formatted  phenotype string in the Irreversible_STATS: block %s\n\t code requires one phenotype name, then a single = character bracketed by two spaces, then a list of stable states or cycles; e.g.:\n\tEpithelial = (SNAI1=0, Twist=0, miR_34=1, Ecadherin_mRNA=1)\n\tPI3KOscillation = (p110_H=1 -> AKT_H=1 -> Nedd4L=1 -> p110_H=0 -> FoxO3=1)\t skipping line, press ENTER to continue.\n",GSE_line);
                        getchar ();
                    }
                    cc = getLastWord(GSE_line);
                    if(cc == "STOP_ALL") Phen_STATs.Irreversible_STATS_GO_ON.push_back(0);
                    else {if(cc == "STOP_REVERSIBLE_COUNTS") Phen_STATs.Irreversible_STATS_GO_ON.push_back(1);
                        else {if(cc == "CONTINUE_REVERSIBLE_COUNTS") Phen_STATs.Irreversible_STATS_GO_ON.push_back(2);
                            else Phen_STATs.Irreversible_STATS_GO_ON.push_back(0);
                        }
                    }
                  }
            }
        }
        if (categoryname == R_STAT_block) {
            while ((fgets(GSE_line,MAX_LINE,f)!=NULL)&&((GSE_line[0]=='\t')||(GSE_line[0]=='-'))) {
                istringstream iss2(GSE_line);
                string p1,p2, bla;
                iss2 >> p1;
                iss2 >> p2;
                iss2 >> bla;
                if((p1.c_str()[0] != '-') && (bla == "=")){
                    bool c1 = false;
                    std::vector<std::string> word = Extract_String_Between_Delimiters_On_Same_Level(GSE_line, '(',')');
                    a=extract_sequence_of_states(PD,word[0], &c1);
                    if (a->N>0){
                        bool c2 = false;
                        b=extract_sequence_of_states(PD,word[1], &c2);
                        if (b->N>0){
                            Phen_STATs.Reversible_Transitions_1.push_back(p1);
                            Phen_STATs.Reversible_Transitions_2.push_back(p2);
                            Phen_STATs.Reversible_Transitions_1_NodeStates.push_back(new longint_ARRAY_pair(a));
                            Phen_STATs.Reversible_Transitions_1_LimitCycle.push_back(c1);
                            Phen_STATs.Reversible_Transitions_2_NodeStates.push_back(new longint_ARRAY_pair(b));
                            Phen_STATs.Reversible_Transitions_2_LimitCycle.push_back(c2);
                            Phen_STATs.Rev =true;
                        }
                        else {printf("Incorrectly formatted second phenotype string in the Reversible_Transitions: block %s\n\t code requires two phenotype names, then a single = character bracketed by two spaces, then two lists of stable states or cycles; e.g.:\n\tEpithelial Mesenchymal = (SNAI1=0, Twist=0, miR_34=1, Ecadherin_mRNA=1) (SNAI1=1, Twist=1, N_bcatenin_H=0, ZEB1_H=0, Ecadherin_mRNA=1)\n\t\t skipping line, press ENTER to continue.\n",GSE_line);
                            getchar ();}
                    }
                    else {
                        printf("Incorrectly formatted first phenotype string in the Reversible_Transitions: block %s\n\t code requires two phenotype names, then a single = character bracketed by two spaces, then two lists of stable states or cycles; e.g.:\n\tEpithelial Mesenchymal = (SNAI1=0, Twist=0, miR_34=1, Ecadherin_mRNA=1) (SNAI1=1, Twist=1, N_bcatenin_H=0, ZEB1_H=0, Ecadherin_mRNA=1)\n\t\t skipping line, press ENTER to continue.\n",GSE_line);
                        getchar ();
                    }
                 }
                else {
                    if(p1.c_str()[0] != '-') {
                        printf("Incorrectly formatted line in the Reversible_Transitions: block: %s\n\t code requires two phenotype names, then a single = character bracketed by two spaces, then two lists of stable states or cycles; e.g.:\n\tEpithelial Mesenchymal = (SNAI1=0, Twist=0, miR_34=1, Ecadherin_mRNA=1) (SNAI1=1, Twist=1, N_bcatenin_H=0, ZEB1_H=0, Ecadherin_mRNA=1)\n\t\t skipping line, press ENTER to continue.\n",GSE_line);
                        getchar ();
                    }
                }
            }
        }
    }
    
    if(Phen_STATs.CellCycleErrors==false) Phen_STATs.WT_cell_cycle=1;
}

void draw_timestep_MOD_RND_Hit(longint_ARRAY_pair *Hits,
                               unsigned int *gothit,
                               MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int t,
                               FILE *f,double y0,int W){
    
    long int i;
    int poz=0;
    unsigned long int hitid;
    
    for(int m=1;m<=PD->BRN->Module_Order->N;m++){
        if((PD->BRN->Module_Order->A[m]>0)&&(PD->BRN->Module_Order->A[m]<=PD->BRN->MDAT->Module_NR)){
            // printf("Module %d\t",m);
            // printf("Going for it %ld \n",PD->BRN->Module_Order->A[m]);
            // printf("\tNodes: %ld \n",PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N);
            for(int nn=1;nn<=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N;nn++){
                i=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->A[nn];
                //   printf("\t nn=%d node = %ld %s\n",nn,i,PD->BRN->MDAT->Node_name[i]);
                hitid=Hits->check_for_element_A(i);
                if(hitid>0){
                    if(gothit[hitid]==0){
                        if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,1);
                        else EpsSetRgb(f,1,0.5,0);
                    }
                    else{
                        if (PD->Coupled->s[i-1]==48) EpsSetRgb(f,0,0,0);
                        else EpsSetRgb(f,1,1,0);
                    }
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

void draw_timetrace_stochastic_inputs_and_KO_OE_hits_experiment_on_single_attractor_MOD(
                                            Environment_and_hits *EnvHit,
                                            int a_st,int a_st_k,
                                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                            double p_error,
                                            FILE *f,double y0,int W){
    int id_GF=0,id_CDL=0, s_start_GF=0,s_start_CDL=0;
    unsigned int *gothit;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    
    if(id_GF>0) s_start_GF=PD->Coupled->s[id_GF-1];
    if(id_CDL>0) s_start_CDL=PD->Coupled->s[id_CDL-1];

    gothit = new unsigned int[EnvHit->BKR_Hits->N+1];
    for (int t=1;t<=T_MAX_DRAW; t++) {
        for (unsigned int inp=1; inp<=PD->Coupled->BRN->active_inputs->N; inp++) {
            double r=rand_real(0,1);
            if(r<EnvHit->p_inputs[inp]) PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp]-1]=49;
            else PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp]-1]=48;
        }
        
        for(unsigned long int l=1;l<=EnvHit->BKR_Hits->N;l++){
            double r=rand_real(0,1);
            if(r<EnvHit->p_hits[l]) {
                PD->Coupled->s[EnvHit->BKR_Hits->A[l]-1]=48+EnvHit->BKR_Hits->B[l];
                gothit[l]=1;
            }
            else gothit[l]=0;
        }
        draw_timestep_MOD_RND_Hit(EnvHit->BKR_Hits,gothit,PD,t,f,y0,W);
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
     
        if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
        if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
    }
}

void Asynchronous_draw_timetrace_stochastic_inputs_and_KO_OE_hits_experiment_on_single_attractor_MOD(
                                            Environment_and_hits *EnvHit,
                                            int a_st,int a_st_k,
                                            MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                            double p_error,
                                            FILE *f,double y0,int W){
    int id_GF=0,id_CDL=0,s_start_GF=0,s_start_CDL=0;
    unsigned int *gothit;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
    id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
    
    if(id_GF>0) s_start_GF=PD->Coupled->s[id_GF-1];
    if(id_CDL>0) s_start_CDL=PD->Coupled->s[id_CDL-1];
    
    gothit = new unsigned int[EnvHit->BKR_Hits->N+1];
    for (int t=1;t<=T_MAX_DRAW; t++) {
        // hits random forcing with p_hits is embedded in asynch update, since these nodes are updated at the top of the order of a timestep and then removed from the randomized list
       for(unsigned long int l=1;l<=EnvHit->BKR_Hits->N;l++){
            double r=rand_real(0,1);
            if(r<EnvHit->p_hits[l]) {
                PD->Coupled->s[EnvHit->BKR_Hits->A[l]-1]=48+EnvHit->BKR_Hits->B[l];
                gothit[l]=1;
            }
            else gothit[l]=0;
        }
        
        draw_timestep_MOD_RND_Hit(EnvHit->BKR_Hits,gothit,PD,t,f,y0,W);
        if(ASSYNC_BIAS == 0 ) PD->Coupled->Asynchronously_update_all_Gates(EnvHit->BKR_Hits,EnvHit->p_hits);
        else PD->Coupled->Asynchronously_update_all_Gates(ASSYNC_BIAS,EnvHit->BKR_Hits,EnvHit->p_hits);
        
        if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
        if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
    }
}


void run_stochastic_inputs_and_KO_OE_hits_experiment_on_single_attractor_MOD(
                                         Environment_and_hits *EnvHit,
                                         int a_st,int a_st_k,
                                         MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    int W=8,YTOP,XTOP;
    FILE *f;
    char fn[400],fn1[400],fn2[300],fn22[300];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    //  long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a+1;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*650+5;
    
    if(EnvHit->BKR_Hits->N>0){
        sprintf(fn2, "DRAW_Hits");
        for (unsigned long int l=1; l<=EnvHit->BKR_Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[l]],EnvHit->BKR_Hits->B[l]);
        }
    }
    else sprintf(fn2,"WILDTYPE");
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD->BRN->path,fn2);
    system(fn);
    
    if(EnvHit->BKR_Hits->N>0){
        sprintf(fn22, "DRAW_Hits");
        for (unsigned long int l=1; l<=EnvHit->BKR_Hits->N; l++) {
            sprintf(fn22,"%s_%s-%ld-%.2lf",fn22,PD->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[l]],EnvHit->BKR_Hits->B[l],EnvHit->p_hits[l]);
        }
    }
    else sprintf(fn22,"WILDTYPE");

    sprintf(fn,"%s/%s/_EXP/%s/%s_Stochastic_A%d_As-%d__",MODEL_Directory,PD->BRN->path,fn2,fn22,a_st,a_st_k);
    for (unsigned int i=1; i<=PD->Coupled->BRN->active_inputs->N; i++) {
        sprintf(fn,"%s_%s-%.2lf",fn,PD->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[i]],EnvHit->p_inputs[i]);
    }
    sprintf(fn1,"%s.eps",fn);

    printf("%s\n",fn1);
    f=fopen(fn1,"w");
    
    EpsInit(f,-5,-5,XTOP,YTOP);
    EpsSetFont(f,"Times-Roman",W);
    draw_timetrace_stochastic_inputs_and_KO_OE_hits_experiment_on_single_attractor_MOD(EnvHit,a_st,a_st_k,PD,p_error,f,YTOP,W);
    
    //getchar();
    sprintf(fn,"%s/%s/_EXP/%s/%s_Stochastic_A%d_As-%d__",MODEL_Directory_system,PD->BRN->path,fn2,fn22,a_st,a_st_k);
    for (unsigned int i=1; i<=PD->Coupled->BRN->active_inputs->N; i++) {
        sprintf(fn,"%s_%s-%.2lf",fn,PD->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[i]],EnvHit->p_inputs[i]);
    }
    close_EPSDRAW_file(f,fn,2);
}

void Asynchronous_run_stochastic_inputs_and_KO_OE_hits_experiment_on_single_attractor_MOD(
                                        Environment_and_hits *EnvHit,
                                        int a_st,int a_st_k,
                                        MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    int W=8,YTOP,XTOP;
    FILE *f;
    char fn[400],fn1[400],fn2[300],fn22[300];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    //  long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a+1;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*650+5;
    
    if(EnvHit->BKR_Hits->N>0){
        sprintf(fn2, "DRAW_Hits");
        for (unsigned long int l=1; l<=EnvHit->BKR_Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[l]],EnvHit->BKR_Hits->B[l]);
        }
    }
    else sprintf(fn2,"WILDTYPE");

    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD->BRN->path,fn2);
    system(fn);

    
    if(EnvHit->BKR_Hits->N>0){
        sprintf(fn22, "DRAW_Hits");
        for (unsigned long int l=1; l<=EnvHit->BKR_Hits->N; l++) {
            sprintf(fn22,"%s_%s-%ld-%.2lf",fn22,PD->BRN->MDAT->Node_name[EnvHit->BKR_Hits->A[l]],EnvHit->BKR_Hits->B[l],EnvHit->p_hits[l]);
        }
    }
    else sprintf(fn2,"WILDTYPE");

    
    sprintf(fn,"%s/%s/_EXP/%s/%s_ASYNCR_Stochastic_A%d_As-%d__",MODEL_Directory,PD->BRN->path,fn2,fn22,a_st,a_st_k);
    for (unsigned int i=1; i<=PD->Coupled->BRN->active_inputs->N; i++) {
        sprintf(fn,"%s_%s-%.2lf",fn,PD->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[i]],EnvHit->p_inputs[i]);
    }
    sprintf(fn1,"%s.eps",fn);
    
    printf("%s\n",fn1);
    f=fopen(fn1,"w");
    
    EpsInit(f,-5,-5,XTOP,YTOP);
    EpsSetFont(f,"Times-Roman",W);
 Asynchronous_draw_timetrace_stochastic_inputs_and_KO_OE_hits_experiment_on_single_attractor_MOD(EnvHit,a_st,a_st_k,PD,p_error,f,YTOP,W);
    
    //getchar();
    sprintf(fn,"%s/%s/_EXP/%s/%s_ASYNCR_Stochastic_A%d_As-%d__",MODEL_Directory_system,PD->BRN->path,fn2,fn22,a_st,a_st_k);
    for (unsigned int i=1; i<=PD->Coupled->BRN->active_inputs->N; i++) {
        sprintf(fn,"%s_%s-%.2lf",fn,PD->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[i]],EnvHit->p_inputs[i]);
    }
    close_EPSDRAW_file(f,fn,2);
}

unsigned short int check_network_state_against_stable_phenotype(char *s,longint_ARRAY_pair a){
    
    unsigned int ok=1;
    for(unsigned int i=1;i<=a.N;i++) // each node on the list for this phenotype
        if((s[a.A[i] - 1] - 48) != a.B[i]) ok=0;
    return(ok);
}

void init_stats(Cell_cycle_Irrev_Rev_stats_in_run *Statsnow){
   // Cell_cycle_Irrev_Rev_stats_in_run Statsnow;
    Statsnow->CC.N_cycles=0;
    Statsnow->CC.N_endoredupl_T=0;
    Statsnow->CC.N_endoredupl_G2=0;
    Statsnow->CC.N_endoredupl_Aneup=0;
    
    Statsnow->CC.MAX_INTERVAL=0;
    Statsnow->CC.T_live=0;
    Statsnow->CC.Time_G1=0;
    Statsnow->CC.Time_R=0;
    Statsnow->CC.Time_G2=0;
    Statsnow->CC.Time_M=0;
    Statsnow->CC.Time_Telo=0;
    Statsnow->CC.N_G1=0;
    Statsnow->CC.N_R=0;
    Statsnow->CC.N_G2=0;
    Statsnow->CC.N_M=0;
    Statsnow->CC.N_T=0;
    for (int t=0;t<=MAX_ARREST; t++){
        Statsnow->CC.Dist_G1[t]=0;
        Statsnow->CC.Dist_G2[t]=0;
        Statsnow->CC.Dist_Mitosis[t]=0;
        Statsnow->CC.Dist_Replication[t]=0;
        Statsnow->CC.Dist_Telophase[t]=0;
    }
    for (int t=0;t<Phen_STATs.Irreversible_STATS.size(); t++){
        Irrev_Phenotyope_stats a;
        a.N_events=0;a.T_live=0;
        a.av_time_to_event=0;a.SD_time_to_event=0; a.Half_life=0;
        Statsnow->IRev.push_back(a);
        for (int i=0;i<=MAX_ARREST; i++){
            Statsnow->IRev[t].Dist_Time_to_event[i]=0;
        }
    }
    for (int t=0;t<Phen_STATs.Reversible_Transitions_1.size(); t++){
        Rev_Transition_stats a;
        a.N_event_to_1=0;a.N_event_to_2=0;a.T_live=0;
        a.av_window_in_1=0;a.av_window_in_2=0;
        a.SD_window_in_1=0;a.SD_window_in_2=0;
        a.Half_life_1=0;a.Half_life_2=0;
        a.T_in_1=0;a.T_in_2=0;
        Statsnow->Rev.push_back(a);
        for (int i=0;i<=MAX_ARREST; i++){
            Statsnow->Rev[t].Dist_Time_in_1[i]=0;
            Statsnow->Rev[t].Dist_Time_in_2[i]=0;
        }
    }
}

Cell_cycle_Irrev_Rev_stats_in_run count_for_stochastic_cell_cycles_Stochastic_HITS(Environment_and_hits *EnvHit,
                                                                                   double p_now,
                                                                        int a_st,int a_st_k,
                                                                        MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                                        double p_error,int run_async){
    // 1 single time-course with a fixed set of forced KO/OE locks going in sync, and a fixed set of intermediate env inputs
    
    Cell_cycle_Irrev_Rev_stats_in_run Statsnow;
    
    init_stats(&Statsnow);
    
    longint_ARRAY_pair *allhits;
    allhits=new longint_ARRAY_pair();
    allhits->add_longint_ARRAY_pair(EnvHit->BKR_Hits);
    allhits->add_longint_ARRAY_pair(EnvHit->Scan_Hits);
    double *allp; allp=new double[EnvHit->BKR_Hits->N+EnvHit->Scan_Hits->N+1];
    for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++) allp[t]=EnvHit->p_hits[t]; // BKR_Hits: these we do not scan for in tandem, they are fixed
    for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++) allp[EnvHit->BKR_Hits->N+t]=p_now; // Scan_Hits: these we planned to scan for (outside this function)

    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);

    int r0=0;int t_G1=1;int t_G2=0;int t_M=0;int t_R=0;int t_T=0;
    double r;

    
    for (unsigned long long int t=1;t<=T_max_CellCycleSample; t++) { // these are the stochastic hits being scanned and/or set, collated so all are set randomly here
  
        for(unsigned long int inp=1;inp<=PD->Coupled->BRN->active_inputs->N;inp++){ // these are all the stochastic inputs (across all inputs)
            r=rand_real(0,1);
            if(r<EnvHit->p_inputs[inp])
                    PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp]-1]=49;
            else PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp]-1]=48;
        }
        for(unsigned long int l=1;l<=allhits->N;l++){
            r=rand_real(0,1);
            if(r<allp[l]) PD->Coupled->s[allhits->A[l]-1]=48+allhits->B[l];
            // else PD->Coupled->s[Hits->A[l]-1]=48+(1-Hits->B[l]);
        }

        for(int ii=0; ii < Phen_STATs.Irreversible_Do_Not_Count_CellCycle_Errors_NodeStates.size(); ii++){ //for each phenptype that stiops the cc count
            if(check_network_state_against_stable_phenotype(PD->Coupled->s,Phen_STATs.Irreversible_Do_Not_Count_CellCycle_Errors_NodeStates[ii])){
                Statsnow.CC.T_live=t;
                return(Statsnow);
            }
        }
        
        if (r0==0){
            if(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==49) {
                r0=1;
                Statsnow.CC.Dist_G1[t_G1]++; Statsnow.CC.Time_G1+=t_G1; Statsnow.CC.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
            else t_G1++;
        }
        if (r0==1){ // in Replication
            if(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==49) {t_R++; }
            if(PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49) {
                r0=2;
                Statsnow.CC.Dist_Replication[t_R]++; Statsnow.CC.Time_R+=t_R; Statsnow.CC.N_R++;
                //  printf("t=%d \t t_R=%d recorded\n",t,t_R);
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
            }
        }

        if(r0==2){ // 4N DNA accomplished
            if((PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49)
                &&(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==48)
                &&(PD->Coupled->s[Phen_STATs.Metaphase_START-1]==48)
                &&(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==48))  {t_G2++; }
            
            if((t_M>0) &&(PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49)
                &&(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==49)
                &&(PD->Coupled->s[Phen_STATs.Metaphase_START-1]==48)
                &&(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==48)) {
                r0=1;
                //  Statsnow.CC.Dist_G2[t_G2]++; Statsnow.CC.Time_G2+=t_G2; Statsnow.CC.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.CC.Dist_G1[t_G1]++; Statsnow.CC.Time_G1+=t_G1; Statsnow.CC.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.CC.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49)
                &&(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==48)
                &&((PD->Coupled->s[Phen_STATs.Metaphase_START-1]==49)
                    ||(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.CC.Dist_G2[t_G2]++; Statsnow.CC.Time_G2+=t_G2; Statsnow.CC.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==49)) {
                r0=1;
                //  Statsnow.CC.Dist_G2[t_G2]++; Statsnow.CC.Time_G2+=t_G2; Statsnow.CC.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.CC.Dist_G1[t_G1]++; Statsnow.CC.Time_G1+=t_G1; Statsnow.CC.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.CC.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49)&&(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==49)) {r0=4;}
        }
        if(r0==3){ // mitosis, not finished
            if((PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49)
                &&(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==48)
                &&(PD->Coupled->s[Phen_STATs.Metaphase_START-1]==48)
                &&(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==48))  {
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_T++;
            }
            
            if((t_M>0) &&(PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49)
                &&(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==49)
                &&(PD->Coupled->s[Phen_STATs.Metaphase_START-1]==48)
                &&(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==48)) {
                r0=1;
                if(t_T>0)  {Statsnow.CC.Dist_Telophase[t_T]++; Statsnow.CC.Time_Telo+=t_T; Statsnow.CC.N_T++;}
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.CC.Dist_G1[t_G1]++; Statsnow.CC.Time_G1+=t_G1; Statsnow.CC.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.CC.N_endoredupl_Aneup++;
            }
            
            if ((PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49)
                &&(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==48)
                &&((PD->Coupled->s[Phen_STATs.Metaphase_START-1]==49)
                    ||(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==49))) {
                    if(t_G2>0) {
                        Statsnow.CC.Dist_G2[t_G2]++; Statsnow.CC.Time_G2+=t_G2; Statsnow.CC.N_G2++;
                        // printf("t=%d \t t_G2=%d recorded\n",t,t_G2);
                        t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                        r0=3;
                    }
                    t_M++;
                }
            
            if((t_G2>0)&&(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==49)) {
                r0=1;
                Statsnow.CC.Dist_G2[t_G2]++; Statsnow.CC.Time_G2+=t_G2; Statsnow.CC.N_G2++;
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_R++;
                //Statsnow.CC.Dist_G1[t_G1]++; Statsnow.CC.Time_G1+=t_G1; Statsnow.CC.N_G1++;
                //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                Statsnow.CC.N_endoredupl_G2++;
            }
            
            if ((PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49)&&(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==49)) {r0=4;}
        }
        if(r0==4){
            if(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==48){
                if(t_M>0) {
                    Statsnow.CC.Dist_Mitosis[t_M]++; Statsnow.CC.Time_M+=t_M; Statsnow.CC.N_M++;
                    //printf("t=%d \t t_M=%d recorded\n",t,t_M);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
                if(PD->Coupled->s[Phen_STATs.DNA_Replication_START-1]==49) {
                    r0=1;
                    if(t_T>0)  {Statsnow.CC.Dist_Telophase[t_T]++; Statsnow.CC.Time_Telo+=t_T; Statsnow.CC.N_T++;}
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                    t_R++;
                    //Statsnow.CC.Dist_G1[t_G1]++; Statsnow.CC.Time_G1+=t_G1; Statsnow.CC.N_G1++;
                    //  printf("t=%d \t t_G1=%d recorded\n",t,t_G1);
                    Statsnow.CC.N_endoredupl_T++;
                }
                else t_T++;
            }
            else t_M++;
            
            if((PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==49)&&(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==48)&&(PD->Coupled->s[Phen_STATs.Cytokinesis-1]==49)) {
                r0=4;
                if(t_T>0) {
                    Statsnow.CC.Dist_Telophase[t_T]++; Statsnow.CC.Time_Telo+=t_T; Statsnow.CC.N_T++;
                    // printf("t=%d \t t_T=%d recorded\n",t,t_T);
                    t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                }
            }
        }
        if (r0==4){
            if ((PD->Coupled->s[Phen_STATs.DNA_Replication_DONE-1]==48)&&(PD->Coupled->s[Phen_STATs.SAC_Cleared-1]==48)&&(PD->Coupled->s[Phen_STATs.Cytokinesis-1]==48)){
                r0=0;
                Statsnow.CC.N_cycles++;
                // getchar();
                t_G1=0;t_R=0;t_M=0;t_G2=0;t_T=0;
                t_G1=1;
            }
        }
        if (run_async == -1 ) PD->Coupled->synchronously_update_all_Gates();
        else {
            if(run_async == 0 ) PD->Coupled->Asynchronously_update_all_Gates(allhits,allp);
            else PD->Coupled->Asynchronously_update_all_Gates(run_async,allhits,allp);
        }
    }
    Statsnow.CC.T_live=T_max_CellCycleSample;
    return(Statsnow);
}


void CC_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(
                                          string datadir,
                                          Environment_and_hits *EnvHit,
                                          longint_ARRAY *l,
                                          MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                          int run_assync){
    Cell_cycle_Irrev_Rev_stats_in_run STnow;
    char fn[400],fn2[400];
    FILE *f;
    
    unsigned long long int maxG1=0,maxR=0,maxG2=0,maxM=0,maxT=0;
    unsigned long int *G1,*R,*G2,*M,*T;
   
    strcpy(fn2,datadir.c_str());
    Cell_cycle_Irrev_Rev_stats_in_run *Statsnow;
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){ // for every attractor on the list of states that match the call
            Statsnow=new Cell_cycle_Irrev_Rev_stats_in_run[BOX_KO_OE+2];
            G1=new unsigned long int[BOX_KO_OE+2];
            R=new unsigned long int[BOX_KO_OE+2];
            G2=new unsigned long int[BOX_KO_OE+2];
            M=new unsigned long int[BOX_KO_OE+2];
            T=new unsigned long int[BOX_KO_OE+2];
            unsigned int box=0;
            
            box=0;
            for(double p_now=0;p_now<=1.; p_now+=1./(double)BOX_KO_OE){ // KO_OE sweep
                box++;
                //timetodeath[box]=0;
                init_stats(&(Statsnow[box]));
                
                G1[box]=0;G2[box]=0;R[box]=0;T[box]=0;M[box]=0;
            
                unsigned long long int t_live_total =0;
                while (t_live_total<T_max_CellCycleSample) {
                    STnow=count_for_stochastic_cell_cycles_Stochastic_HITS(EnvHit,p_now,(int)l->A[j],(int)rand_int(1,PD_A->Coupled->Attractor_valleys[(int)l->A[j]]->N_a),PD_A,0,run_assync);
                    t_live_total+=STnow.CC.T_live;
                    
                    Statsnow[box].CC.T_live+=STnow.CC.T_live;
                    Statsnow[box].CC.N_cycles+=STnow.CC.N_cycles;
                    Statsnow[box].CC.N_endoredupl_T+=STnow.CC.N_endoredupl_T;
                    Statsnow[box].CC.N_endoredupl_G2+=STnow.CC.N_endoredupl_G2;
                    Statsnow[box].CC.N_endoredupl_Aneup+=STnow.CC.N_endoredupl_Aneup;
                    Statsnow[box].CC.Time_G1+=STnow.CC.Time_G1;
                    Statsnow[box].CC.Time_R+=STnow.CC.Time_R;
                    Statsnow[box].CC.Time_G2+=STnow.CC.Time_G2;
                    Statsnow[box].CC.Time_M+=STnow.CC.Time_M;
                    Statsnow[box].CC.Time_Telo+=STnow.CC.Time_Telo;
                    
                    for (int t=0;t<=MAX_ARREST; t++){
                        Statsnow[box].CC.Dist_G1[t]+=STnow.CC.Dist_G1[t];
                        Statsnow[box].CC.Dist_G2[t]+=STnow.CC.Dist_G2[t];
                        Statsnow[box].CC.Dist_Mitosis[t]+=STnow.CC.Dist_Mitosis[t];
                        Statsnow[box].CC.Dist_Replication[t]+=STnow.CC.Dist_Replication[t];
                        Statsnow[box].CC.Dist_Telophase[t]+=STnow.CC.Dist_Telophase[t];
                        
                        if(STnow.CC.Dist_G1[t]>0) Statsnow[box].CC.MAX_INTERVAL=t;
                        if(STnow.CC.Dist_G2[t]>0) Statsnow[box].CC.MAX_INTERVAL=t;
                        if(STnow.CC.Dist_Mitosis[t]>0) Statsnow[box].CC.MAX_INTERVAL=t;
                        if(STnow.CC.Dist_Replication[t]>0) Statsnow[box].CC.MAX_INTERVAL=t;
                        if(STnow.CC.Dist_Telophase[t]>0) Statsnow[box].CC.MAX_INTERVAL=t;
                        
                        if((STnow.CC.Dist_G1[t]>0)&&(t>maxG1)) maxG1=t;
                        if((STnow.CC.Dist_Replication[t]>0)&&(t>maxR)) maxR=t;
                        if((STnow.CC.Dist_G2[t]>0)&&(t>maxG2)) maxG2=t;
                        if((STnow.CC.Dist_Mitosis[t]>0)&&(t>maxM)) maxM=t;
                        if((STnow.CC.Dist_Telophase[t]>0)&&(t>maxT)) maxT=t;
                        
                    }
                }
            }
            box=0;
            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOX_KO_OE){ // KO_OE sweep
                box++;
                for (int t=0;t<=MAX_ARREST; t++){
                    G1[box]+=Statsnow[box].CC.Dist_G1[t];
                    R[box]+=Statsnow[box].CC.Dist_Replication[t];
                    G2[box]+=Statsnow[box].CC.Dist_G2[t];
                    M[box]+=Statsnow[box].CC.Dist_Mitosis[t];
                    T[box]+=Statsnow[box].CC.Dist_Telophase[t];
                }
            }
            
            sprintf(fn,"%s/CC-FRACTIONS_E-",fn2);
            for(unsigned int e=1;e<=PD_A->Coupled->BRN->active_inputs->N;e++)
                sprintf(fn,"%s_%s-%.2lf",fn,PD_A->Coupled->BRN->MDAT->Node_name[PD_A->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
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
            sprintf(fn,"%s__A%d__BOX_%d.txt",fn,(int)l->A[j],BOX_KO_OE);
            f=fopen(fn,"w");
            box=0;
            fprintf(f, "p_hits\tNormal\tN_endoredupl_G2\tN_endoredupl_Aneup\tN_endoredupl_T\tTime_G1\tTime_R\tTime_G2\tTime_M\tTime_Telo\n");

            for(double p_hits=0;p_hits<=1.; p_hits+=1./(double)BOX_KO_OE){
                box++;
                fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                        p_hits,
                        Phen_STATs.WT_cell_cycle * Statsnow[box].CC.N_cycles/((double)Statsnow[box].CC.T_live),
                        //100*Statsnow[box].CC.N_deaths/(double)Statsnow[box].CC.T_live,
                        Phen_STATs.WT_cell_cycle * Statsnow[box].CC.N_endoredupl_G2/((double)Statsnow[box].CC.T_live),
                        Phen_STATs.WT_cell_cycle * Statsnow[box].CC.N_endoredupl_Aneup/((double)Statsnow[box].CC.T_live),
                        Phen_STATs.WT_cell_cycle * Statsnow[box].CC.N_endoredupl_T/((double)Statsnow[box].CC.T_live),
                        Phen_STATs.WT_cell_cycle * Statsnow[box].CC.Time_G1/((double)Statsnow[box].CC.T_live),
                        Phen_STATs.WT_cell_cycle * Statsnow[box].CC.Time_R/((double)Statsnow[box].CC.T_live),
                        Phen_STATs.WT_cell_cycle * Statsnow[box].CC.Time_G2/((double)Statsnow[box].CC.T_live),
                        Phen_STATs.WT_cell_cycle * Statsnow[box].CC.Time_M/((double)Statsnow[box].CC.T_live),
                        Phen_STATs.WT_cell_cycle * Statsnow[box].CC.Time_Telo/((double)Statsnow[box].CC.T_live));
            }
            fclose(f);
            delete[] Statsnow;Statsnow=NULL;
        }
    }
}

Cell_cycle_Irrev_Rev_stats_in_run count_irreversible_transitions_Stochastic_HITS(unsigned int vi, Environment_and_hits *EnvHit,
                                                                                   double p_now,
                                                                                   int a_st,int a_st_k,
                                                                                   MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                                                   double p_error,int run_async){
    Cell_cycle_Irrev_Rev_stats_in_run Statsnow;
    init_stats(&Statsnow);
    
    longint_ARRAY_pair *allhits;
    allhits=new longint_ARRAY_pair();
    allhits->add_longint_ARRAY_pair(EnvHit->BKR_Hits);
    allhits->add_longint_ARRAY_pair(EnvHit->Scan_Hits);
    double *allp; allp=new double[EnvHit->BKR_Hits->N+EnvHit->Scan_Hits->N+1];
    for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++) allp[t]=EnvHit->p_hits[t];
    for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++) allp[EnvHit->BKR_Hits->N+t]=p_now;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    double r;
    
    for (unsigned long long int t=1;t<=T_max_CellCycleSample; t++) {
        for(unsigned long int inp=1;inp<=PD->Coupled->BRN->active_inputs->N;inp++){
            r=rand_real(0,1);
            if(r<EnvHit->p_inputs[inp])
                PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp]-1]=49;
            else PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp]-1]=48;
        }
       
        for(unsigned long int l=1;l<=allhits->N;l++){
            r=rand_real(0,1);
            if(r<allp[l]) PD->Coupled->s[allhits->A[l]-1]=48+allhits->B[l];
            // else PD->Coupled->s[Hits->A[l]-1]=48+(1-Hits->B[l]);
        }
         
        for(int ii=0; ii < Phen_STATs.Irreversible_STATS.size(); ii++){
            if((Phen_STATs.Irreversible_STATS_GO_ON[ii] == 0) || ii == vi){
                if(Phen_STATs.Irreversible_STATS_LimitCycle[ii]==0){
                    if(check_network_state_against_stable_phenotype(PD->Coupled->s,Phen_STATs.Irreversible_STATS_NodeStates[ii])){
                        if(ii==vi){
                            Statsnow.IRev[ii].N_events++;
                            Statsnow.IRev[ii].av_time_to_event+=t;
                            Statsnow.IRev[ii].SD_time_to_event+=t*t;
                            if(t<MAX_ARREST) Statsnow.IRev[ii].Dist_Time_to_event[t]++;
                         //   printf("cgecking  %s\t",Phen_STATs.Irreversible_STATS[ii].c_str());
                        //    printf("SAHF = %d\t",PD->Coupled->s[PD->Coupled->BRN->MDAT->get_node_ID("SAHF")-1]);
                        //    printf(" Match!\n");getchar();
                            Statsnow.IRev[ii].T_live=t;
                            return(Statsnow);
                        }
                    }
                }
                else {//// catch cycle
                }
            }
        }
        
        if (run_async == -1 ) PD->Coupled->synchronously_update_all_Gates();
        else {
            if(run_async == 0 ) PD->Coupled->Asynchronously_update_all_Gates(allhits,allp);
            else PD->Coupled->Asynchronously_update_all_Gates(run_async,allhits,allp);
        }
    }
    Statsnow.IRev[vi].T_live=T_max_CellCycleSample;
    return(Statsnow);
}


void Irrev_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(
                                                                              string datadir,
                                                                           Environment_and_hits *EnvHit,
                                                                           longint_ARRAY *l,
                                                                           MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                                                           int run_assync){
    Cell_cycle_Irrev_Rev_stats_in_run STnow;
    char fn[400],fn2[400];
    FILE *f;
    unsigned int boxp=0;
    
    strcpy(fn2,datadir.c_str());
    Cell_cycle_Irrev_Rev_stats_in_run *Statsnow;
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            Statsnow=new Cell_cycle_Irrev_Rev_stats_in_run[BOX_KO_OE+2];
            //printf("boxp=%u\t vect = %lu\n",boxp,Phen_STATs.Irreversible_STATS.size());
          //  getchar();
            for(unsigned int vi=0;vi<Phen_STATs.Irreversible_STATS.size();vi++){
                boxp=0;
                for(double p_now=0;p_now<=1.; p_now+=1./(double)BOX_KO_OE){
                    boxp++;
                    init_stats(&(Statsnow[boxp]));
                    
                    unsigned long long int t_live_total =0;
                    unsigned int n_runs=0;
                    while (t_live_total<T_max_TransitionsSample) {
                        STnow=count_irreversible_transitions_Stochastic_HITS(vi,EnvHit,p_now,(int)l->A[j],(int)rand_int(1,PD_A->Coupled->Attractor_valleys[(int)l->A[j]]->N_a),PD_A,0,run_assync);
                        t_live_total+=STnow.IRev[vi].T_live;
                        n_runs++;
                        if (STnow.IRev[vi].T_live<T_max_CellCycleSample) {
                            Statsnow[boxp].IRev[vi].T_live+=STnow.IRev[vi].T_live;
                            Statsnow[boxp].IRev[vi].N_events+=STnow.IRev[vi].N_events;
                            Statsnow[boxp].IRev[vi].av_time_to_event+=STnow.IRev[vi].av_time_to_event;
                            Statsnow[boxp].IRev[vi].SD_time_to_event+=STnow.IRev[vi].SD_time_to_event;
                            for (int t=0;t<=MAX_ARREST; t++)
                                Statsnow[boxp].IRev[vi].Dist_Time_to_event[t]+=STnow.IRev[vi].Dist_Time_to_event[t];
                        }
                    }

                    if(Statsnow[boxp].IRev[vi].N_events>0)
                        Statsnow[boxp].IRev[vi].av_time_to_event/=(double) Statsnow[boxp].IRev[vi].N_events;
                    Statsnow[boxp].IRev[vi].SD_time_to_event=sigma_from_average_and_sum_of_squares(Statsnow[boxp].IRev[vi].av_time_to_event,Statsnow[boxp].IRev[vi].SD_time_to_event,Statsnow[boxp].IRev[vi].N_events);
                    for (int t=0;t<=MAX_ARREST; t++)
                        if(Statsnow[boxp].IRev[vi].Dist_Time_to_event[t]>=(double)(n_runs)/2.) Statsnow[boxp].IRev[vi].Half_life = t;
                }
                
                sprintf(fn,"%s/Irrev_Stats_%s-",fn2,Phen_STATs.Irreversible_STATS[vi].c_str());
                for(unsigned int e=1;e<=PD_A->Coupled->BRN->active_inputs->N;e++)
                    sprintf(fn,"%s_%s-%.2lf",fn,PD_A->Coupled->BRN->MDAT->Node_name[PD_A->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
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
                sprintf(fn,"%s__A%d__BOX_%d.txt",fn,(int)l->A[j],BOX_KO_OE);
            
                f=fopen(fn,"w");
                fprintf(f,"p_hit\t%s_rate_rel_CC\tAV_time_to_%s\tSD_time_to_%s\tHalfLife_%s\n",
                        Phen_STATs.Irreversible_STATS[vi].c_str(),
                        Phen_STATs.Irreversible_STATS[vi].c_str(),
                        Phen_STATs.Irreversible_STATS[vi].c_str(),
                        Phen_STATs.Irreversible_STATS[vi].c_str());
                boxp=0;
                for(double p_now=0;p_now<=1.; p_now+=1./(double)BOX_KO_OE){
                    boxp++;
                    //if (Statsnow[box].IRev[vi].N_events>0)
                    if(Statsnow[boxp].IRev[vi].T_live==0) {
                        Statsnow[boxp].IRev[vi].N_events=0;
                        Statsnow[boxp].IRev[vi].T_live=1;
                    }
                    fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\n",
                            p_now,
                            Phen_STATs.WT_cell_cycle * Statsnow[boxp].IRev[vi].N_events/((double)Statsnow[boxp].IRev[vi].T_live),
                            Statsnow[boxp].IRev[vi].av_time_to_event,
                            Statsnow[boxp].IRev[vi].SD_time_to_event,
                            Statsnow[boxp].IRev[vi].Half_life);
                    //else
                }
                fclose(f);
            }
            
            delete[] Statsnow;Statsnow=NULL;
        }
    }
}

Cell_cycle_Irrev_Rev_stats_in_run count_Reversible_transitions_Stochastic_HITS(unsigned int vi,
                                                                               Environment_and_hits *EnvHit,
                                                                                 double p_now,
                                                                                 int a_st,int a_st_k,
                                                                                 MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                                                 double p_error,int run_async){
    Cell_cycle_Irrev_Rev_stats_in_run Statsnow;
    init_stats(&Statsnow);
    
    longint_ARRAY_pair *allhits;
    allhits=new longint_ARRAY_pair();
    allhits->add_longint_ARRAY_pair(EnvHit->BKR_Hits);
    allhits->add_longint_ARRAY_pair(EnvHit->Scan_Hits);
    double *allp; allp=new double[EnvHit->BKR_Hits->N+EnvHit->Scan_Hits->N+1];
    for (unsigned long int t=1;t<=EnvHit->BKR_Hits->N; t++) allp[t]=EnvHit->p_hits[t];
    for (unsigned long int t=1;t<=EnvHit->Scan_Hits->N; t++) allp[EnvHit->BKR_Hits->N+t]=p_now;
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    double r;
    bool *first_transition_ok;
    unsigned int *first_state,*previous_state;
    first_transition_ok=new bool [Phen_STATs.Reversible_Transitions_1.size()];
    first_state=new unsigned int [Phen_STATs.Reversible_Transitions_1.size()];
    previous_state=new unsigned int [Phen_STATs.Reversible_Transitions_1.size()];
    for(int ii=0; ii < Phen_STATs.Reversible_Transitions_1.size(); ii++) {
        first_transition_ok[ii]=0;
        first_state[ii]=0;
        previous_state[ii]=0;
    }
    
    unsigned long long int t_last_transition = 1;
    for (unsigned long long int t=1;t<=T_max_CellCycleSample; t++) {
        for(unsigned long int inp=1;inp<=PD->Coupled->BRN->active_inputs->N;inp++){
            r=rand_real(0,1);
            if(r<EnvHit->p_inputs[inp])
                PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp]-1]=49;
            else PD->Coupled->s[PD->Coupled->BRN->active_inputs->A[inp]-1]=48;
        }
       
        for(unsigned long int l=1;l<=allhits->N;l++){
            r=rand_real(0,1);
            if(r<allp[l]) PD->Coupled->s[allhits->A[l]-1]=48+allhits->B[l];
            // else PD->Coupled->s[Hits->A[l]-1]=48+(1-Hits->B[l]);
        }
        
        for(int ii=0; ii < Phen_STATs.Irreversible_STATS.size(); ii++){
            if(Phen_STATs.Irreversible_STATS_GO_ON[ii] < 2) {
                if(Phen_STATs.Irreversible_STATS_LimitCycle[ii]==0){
                    if(check_network_state_against_stable_phenotype(PD->Coupled->s,Phen_STATs.Irreversible_STATS_NodeStates[ii])){
                            Statsnow.Rev[vi].T_live=t;
                            return(Statsnow);
                        }
                    }
                }
                else {//// catch cycle
            }
        }
        
        if(Phen_STATs.Reversible_Transitions_1_LimitCycle[vi]==0){
            if(check_network_state_against_stable_phenotype(PD->Coupled->s,Phen_STATs.Reversible_Transitions_1_NodeStates[vi])){
                Statsnow.Rev[vi].T_in_1++;
                if(first_transition_ok[vi]==true){
                    if(previous_state[vi]==2){
                        Statsnow.Rev[vi].N_event_to_1++;
                        Statsnow.Rev[vi].av_window_in_2+=(t-t_last_transition);
                        Statsnow.Rev[vi].SD_window_in_2+=(t-t_last_transition)*(t-t_last_transition);
                        if(t<MAX_ARREST) Statsnow.Rev[vi].Dist_Time_in_2[t]++;
                        previous_state[vi]=1;
                        t_last_transition = t;
                    }
                }
                else{
                    if(first_state[vi]==0) {first_state[vi]=1;previous_state[vi]=1;}
                    else if(first_transition_ok[vi]==false) first_transition_ok[vi]=true;
                }
            }
        }
        else {//// catch cycle
        }
        
        if(Phen_STATs.Reversible_Transitions_2_LimitCycle[vi]==0){
            if(check_network_state_against_stable_phenotype(PD->Coupled->s,Phen_STATs.Reversible_Transitions_2_NodeStates[vi])){
                Statsnow.Rev[vi].T_in_2++;
                if(first_transition_ok[vi]==true){
                    if(previous_state[vi]==1){
                        Statsnow.Rev[vi].N_event_to_2++;
                        Statsnow.Rev[vi].av_window_in_1+=(t-t_last_transition);
                        Statsnow.Rev[vi].SD_window_in_1+=(t-t_last_transition)*(t-t_last_transition);
                        if(t<MAX_ARREST) Statsnow.Rev[vi].Dist_Time_in_1[t]++;
                        previous_state[vi]=2;
                        t_last_transition = t;
                    }
                }
                else{
                    if(first_state[vi]==0) {first_state[vi]=2;previous_state[vi]=2;}
                    else if(first_transition_ok[vi]==false) first_transition_ok[vi]=true;
                }
            }
        }
        else {//// catch cycle
        }
        
        
        if (run_async == -1 ) PD->Coupled->synchronously_update_all_Gates();
        else {
            if(run_async == 0 ) PD->Coupled->Asynchronously_update_all_Gates(allhits,allp);
            else PD->Coupled->Asynchronously_update_all_Gates(ASSYNC_BIAS,allhits,allp);
        }
    }
    Statsnow.Rev[vi].T_live=T_max_CellCycleSample;
    return(Statsnow);
}

void RevTrans_run_stochastic_inputs_and_KO_OE_hits_STATS_SWEEP_on_attractor_list(
                                                                                 string datadir,
                                                                              Environment_and_hits *EnvHit,
                                                                              longint_ARRAY *l,
                                                                              MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD_A,
                                                                              int run_assync){
    Cell_cycle_Irrev_Rev_stats_in_run STnow;
    char fn[400],fn2[400];
    FILE *f;
    
    strcpy(fn2,datadir.c_str());
    
    Cell_cycle_Irrev_Rev_stats_in_run *Statsnow;
    
    if (l!=NULL){
        for(int j=1;j<=l->N;j++){
            Statsnow=new Cell_cycle_Irrev_Rev_stats_in_run[BOX_KO_OE+2];
            unsigned int boxp=0;
            
            for(unsigned int vi=0;vi<Phen_STATs.Reversible_Transitions_1.size();vi++){
                boxp=0;
                for(double p_now=0;p_now<=1.; p_now+=1./(double)BOX_KO_OE){
                    boxp++;
                    init_stats(&(Statsnow[boxp]));
                    
                    unsigned long long int t_live_total =0;
                    unsigned int n_runs=0;
                    while (t_live_total<T_max_TransitionsSample) {
                        STnow=count_Reversible_transitions_Stochastic_HITS(vi,EnvHit,p_now,(int)l->A[j],(int)rand_int(1,PD_A->Coupled->Attractor_valleys[(int)l->A[j]]->N_a),PD_A,0,run_assync);
                        t_live_total+=STnow.Rev[vi].T_live;
                        n_runs++;
                       // printf("vi=%d\t time in 1 = %u\t  time in 2 = %u\n",vi,STnow.Rev[vi].T_in_1,STnow.Rev[vi].T_in_2);
                       // getchar();
                        
                        Statsnow[boxp].Rev[vi].T_live+=STnow.Rev[vi].T_live;
                        Statsnow[boxp].Rev[vi].N_event_to_1+=STnow.Rev[vi].N_event_to_1;
                        Statsnow[boxp].Rev[vi].N_event_to_2+=STnow.Rev[vi].N_event_to_2;
                        Statsnow[boxp].Rev[vi].T_in_1+=STnow.Rev[vi].T_in_1;
                        Statsnow[boxp].Rev[vi].T_in_2+=STnow.Rev[vi].T_in_2;
                        Statsnow[boxp].Rev[vi].av_window_in_1+=STnow.Rev[vi].av_window_in_1;
                        Statsnow[boxp].Rev[vi].av_window_in_2+=STnow.Rev[vi].av_window_in_2;
                        Statsnow[boxp].Rev[vi].SD_window_in_1+=STnow.Rev[vi].SD_window_in_1;
                        Statsnow[boxp].Rev[vi].SD_window_in_2+=STnow.Rev[vi].SD_window_in_2;
                        for (int i=0;i<=MAX_ARREST; i++){
                            Statsnow[boxp].Rev[vi].Dist_Time_in_1[i]+=STnow.Rev[vi].Dist_Time_in_1[i];
                            Statsnow[boxp].Rev[vi].Dist_Time_in_2[i]+=STnow.Rev[vi].Dist_Time_in_2[i];
                        }
                    }
                    
                    if(Statsnow[boxp].Rev[vi].N_event_to_2>0)
                        Statsnow[boxp].Rev[vi].av_window_in_1/=(double) Statsnow[boxp].Rev[vi].N_event_to_2;
                    if(Statsnow[boxp].Rev[vi].N_event_to_1>0)
                        Statsnow[boxp].Rev[vi].av_window_in_2/=(double) Statsnow[boxp].Rev[vi].N_event_to_1;
                    Statsnow[boxp].Rev[vi].SD_window_in_1=sigma_from_average_and_sum_of_squares(
                                               Statsnow[boxp].Rev[vi].av_window_in_1,
                                               Statsnow[boxp].Rev[vi].SD_window_in_1,
                                               Statsnow[boxp].Rev[vi].N_event_to_2);
                    Statsnow[boxp].Rev[vi].SD_window_in_2=sigma_from_average_and_sum_of_squares(
                                               Statsnow[boxp].Rev[vi].av_window_in_2,
                                               Statsnow[boxp].Rev[vi].SD_window_in_2,
                                               Statsnow[boxp].Rev[vi].N_event_to_1);
                    for (int t=0;t<=MAX_ARREST; t++){
                        if(Statsnow[boxp].Rev[vi].N_event_to_2>0)
                            if(Statsnow[boxp].Rev[vi].Dist_Time_in_1[t]>=(double)(Statsnow[boxp].Rev[vi].N_event_to_2)/2.)
                                Statsnow[boxp].Rev[vi].Half_life_1 = t;
                        if(Statsnow[boxp].Rev[vi].N_event_to_1>0)
                            if(Statsnow[boxp].Rev[vi].Dist_Time_in_2[t]>=(double)(Statsnow[boxp].Rev[vi].N_event_to_1)/2.)
                                Statsnow[boxp].Rev[vi].Half_life_2 = t;
                    }
                }
                
                sprintf(fn,"%s/RevTransition_Stats_%s_vs_%s-",fn2,Phen_STATs.Reversible_Transitions_1[vi].c_str(),Phen_STATs.Reversible_Transitions_2[vi].c_str());
                for(unsigned int e=1;e<=PD_A->Coupled->BRN->active_inputs->N;e++)
                    sprintf(fn,"%s_%s-%.2lf",fn,PD_A->Coupled->BRN->MDAT->Node_name[PD_A->Coupled->BRN->active_inputs->A[e]],EnvHit->p_inputs[e]);
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
                sprintf(fn,"%s__A%d__BOX_%d.txt",fn,(int)l->A[j],BOX_KO_OE);
                
                f=fopen(fn,"w");
                fprintf(f,"p_hit\tT_in_%s\trate_to_%s_rel_CC\tAV_window_in_%s\tSD_window_in_%s\tHalfLife_%s\tT_in_%s\trate_to_%s_rel_CC\tAV_window_in_%s\tSD_window_in_%s\tHalfLife_%s\n",
                        Phen_STATs.Reversible_Transitions_1[vi].c_str(),
                        Phen_STATs.Reversible_Transitions_1[vi].c_str(),
                        Phen_STATs.Reversible_Transitions_1[vi].c_str(),
                        Phen_STATs.Reversible_Transitions_1[vi].c_str(),
                        Phen_STATs.Reversible_Transitions_1[vi].c_str(),
                        Phen_STATs.Reversible_Transitions_2[vi].c_str(),
                        Phen_STATs.Reversible_Transitions_2[vi].c_str(),
                        Phen_STATs.Reversible_Transitions_2[vi].c_str(),
                        Phen_STATs.Reversible_Transitions_2[vi].c_str(),
                        Phen_STATs.Reversible_Transitions_2[vi].c_str());
                unsigned int boxp=0;
                
               for(double p_now=0;p_now<=1.; p_now+=1./(double)BOX_KO_OE){
                   boxp++;
                   if(Statsnow[boxp].Rev[vi].T_live==0) {
                       Statsnow[boxp].Rev[vi].N_event_to_1=0;
                       Statsnow[boxp].Rev[vi].N_event_to_2=0;
                       Statsnow[boxp].Rev[vi].T_live=1;
                   }
                   
                   if(Statsnow[boxp].Rev[vi].N_event_to_1 * Statsnow[boxp].Rev[vi].N_event_to_2 >0){
                       fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                                p_now,
                               Statsnow[boxp].Rev[vi].T_in_1/(double)(Statsnow[boxp].Rev[vi].T_live),
                               Phen_STATs.WT_cell_cycle * Statsnow[boxp].Rev[vi].N_event_to_1/(Statsnow[boxp].Rev[vi].av_window_in_2 * Statsnow[boxp].Rev[vi].N_event_to_2),
                                Statsnow[boxp].Rev[vi].av_window_in_1,
                                Statsnow[boxp].Rev[vi].SD_window_in_1,
                                Statsnow[boxp].Rev[vi].Half_life_1,
                                Statsnow[boxp].Rev[vi].T_in_2/(double)Statsnow[boxp].Rev[vi].T_live,
                               Phen_STATs.WT_cell_cycle * Statsnow[boxp].Rev[vi].N_event_to_2/(Statsnow[boxp].Rev[vi].av_window_in_1 * Statsnow[boxp].Rev[vi].N_event_to_1),
                                Statsnow[boxp].Rev[vi].av_window_in_2,
                                Statsnow[boxp].Rev[vi].SD_window_in_2,
                                Statsnow[boxp].Rev[vi].Half_life_2);
                   }
                   else  fprintf(f, "%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                                 p_now,
                                 Statsnow[boxp].Rev[vi].T_in_1/(double)Statsnow[boxp].Rev[vi].T_live,
                                 0.,
                                 Statsnow[boxp].Rev[vi].av_window_in_1,
                                 Statsnow[boxp].Rev[vi].SD_window_in_1,
                                 Statsnow[boxp].Rev[vi].Half_life_1,
                                 Statsnow[boxp].Rev[vi].T_in_2/(double)Statsnow[boxp].Rev[vi].T_live,
                                 0.,
                                 Statsnow[boxp].Rev[vi].av_window_in_2,
                                 Statsnow[boxp].Rev[vi].SD_window_in_2,
                                 Statsnow[boxp].Rev[vi].Half_life_2);
                }
                fclose(f);
            }
            delete[] Statsnow;Statsnow=NULL;
        }
    }
}

void draw_timetrace_KO_OE_experiment_on_single_attractor_MOD(
                                                             Environment_and_hits *EnvHit,
                                                             int a_st,int a_st_k,
                                                             MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                                                             double p_error,
                                                             FILE *f,double y0,int W){
     int id_GF=0,id_CDL=0, s_start_GF=0,s_start_CDL=0;
     unsigned int *gothit;
     
     PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[1]]->s);
     
     id_GF=PD->BRN->MDAT->get_node_ID(GF_low_name);
     id_CDL=PD->BRN->MDAT->get_node_ID(CD_low_name);
     
     if(id_GF>0) s_start_GF=PD->Coupled->s[id_GF-1];
     if(id_CDL>0) s_start_CDL=PD->Coupled->s[id_CDL-1];

     gothit = new unsigned int[EnvHit->Scan_Hits->N+1];
    
// Initial state"
    for (int t=1;t<=maximum(40,2*PD->Coupled->Attractor_valleys[a_st]->N_a); t++) {
        draw_timestep_MOD_RND_Hit(EnvHit->Scan_Hits,gothit,PD,t,f,y0,W);
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);

    }
    unsigned int t0=maximum(40,2*PD->Coupled->Attractor_valleys[a_st]->N_a);
// knockdown at state  a_st_k along cycle
     for (int t=t0+1;t<=t0+PD->Coupled->Attractor_valleys[a_st]->N_a; t++) {
         if(t>=t0+a_st_k){
             for(unsigned long int l=1;l<=EnvHit->Scan_Hits->N;l++){
                 PD->Coupled->s[EnvHit->Scan_Hits->A[l]-1]=48+EnvHit->Scan_Hits->B[l];
                 gothit[l]=1;
             }
             draw_timestep_MOD_RND_Hit(EnvHit->Scan_Hits,gothit,PD,t,f,y0,W);
             PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
             
             if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
             if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
         }
         else{
             draw_timestep_MOD_RND_Hit(EnvHit->Scan_Hits,gothit,PD,t,f,y0,W);
             PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
         }
     }
    t0+=PD->Coupled->Attractor_valleys[a_st]->N_a;
// Followup with knockdown in place state"
    for (int t=t0+1;t<=t0+maximum(40,2*PD->Coupled->Attractor_valleys[a_st]->N_a); t++) {
        for(unsigned long int l=1;l<=EnvHit->Scan_Hits->N;l++){
            PD->Coupled->s[EnvHit->Scan_Hits->A[l]-1]=48+EnvHit->Scan_Hits->B[l];
            gothit[l]=1;
        }
        draw_timestep_MOD_RND_Hit(EnvHit->Scan_Hits,gothit,PD,t,f,y0,W);
        PD->Coupled->synchronously_update_all_Gates_with_noise(p_error);
        
        if(id_GF>0) PD->Coupled->s[id_GF-1]=s_start_GF;
        if(id_CDL>0) PD->Coupled->s[id_CDL-1]=s_start_CDL;
    }
 }


void run_timed_full_KO_OE_experiment_on_single_attractor_MOD(
                                                             Environment_and_hits *EnvHit,
                                                             int a_st,int a_st_k,
                                                             MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,double p_error){
    int W=8,YTOP,XTOP;
    FILE *f;
    char fn[400],fn1[400],fn2[300],fn22[300];
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    // PD->Coupled->print_state();//getchar();
    
    //  long int DMAX=PD->Coupled->Attractor_valleys[a_st]->N_a+1;
    int N_lines=PD->BRN->N+10;
    YTOP=(W*N_lines+10*W)+5;
    XTOP=2*650+5;
    
    if(EnvHit->Scan_Hits->N>0){
        sprintf(fn2, "DRAW_Hits");
        for (unsigned long int l=1; l<=EnvHit->Scan_Hits->N; l++) {
            sprintf(fn2,"%s_%s-%ld",fn2,PD->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[l]],EnvHit->Scan_Hits->B[l]);
        }
    }
    else sprintf(fn2,"WILDTYPE");
    
    sprintf(fn,"mkdir %s/%s/_EXP\n",MODEL_Directory_system,PD->BRN->path);
    system(fn);
    sprintf(fn,"mkdir %s/%s/_EXP/%s\n",MODEL_Directory_system,PD->BRN->path,fn2);
    system(fn);
    
    if(EnvHit->Scan_Hits->N>0){
        sprintf(fn22, "DRAW_Hits");
        for (unsigned long int l=1; l<=EnvHit->Scan_Hits->N; l++) {
            sprintf(fn22,"%s_%s-%ld",fn22,PD->BRN->MDAT->Node_name[EnvHit->Scan_Hits->A[l]],EnvHit->Scan_Hits->B[l]);
        }
    }
    else sprintf(fn22,"WILDTYPE");

    sprintf(fn,"%s/%s/_EXP/%s/%s_TIMED_A%d_As-%d__",MODEL_Directory,PD->BRN->path,fn2,fn22,a_st,a_st_k);
    for (unsigned int i=1; i<=PD->Coupled->BRN->active_inputs->N; i++) {
        sprintf(fn,"%s_%s-%.2lf",fn,PD->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[i]],EnvHit->p_inputs[i]);
    }
    sprintf(fn1,"%s.eps",fn);

    printf("%s\n",fn1);
    f=fopen(fn1,"w");
    
    EpsInit(f,-5,-5,XTOP,YTOP);
    EpsSetFont(f,"Times-Roman",W);
    draw_timetrace_KO_OE_experiment_on_single_attractor_MOD(EnvHit,a_st,a_st_k,PD,p_error,f,YTOP,W);
   
    //getchar();
    sprintf(fn,"%s/%s/_EXP/%s/%s_TIMED_A%d_As-%d__",MODEL_Directory_system,PD->BRN->path,fn2,fn22,a_st,a_st_k);
    for (unsigned int i=1; i<=PD->Coupled->BRN->active_inputs->N; i++) {
        sprintf(fn,"%s_%s-%.2lf",fn,PD->BRN->MDAT->Node_name[PD->Coupled->BRN->active_inputs->A[i]],EnvHit->p_inputs[i]);
    }
    close_EPSDRAW_file(f,fn,2);
}



#endif /* Timecourses_and_stats_module_h */
