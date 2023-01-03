//
//  DynamicalModularity.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 3/13/20.
//  Copyright © 2020 Regan_Group. All rights reserved.
//
using namespace std;

class DM_Phenotype {
public:
    string DM_Phen_Name;
    vector<string> DM_Phen_StateCycle; // 00100110, 00100101, ...
    vector<string> DM_Phen_PIN; // __1__11_, __1__10_, ...
    vector<string> Modules_to_stop_transition_Stats_time_for;
        // the names of these modules indicate which module-transitions I stop
        // tracking once this phenotyope is observed.
        // look for inconsistency if those other modules do still flip!
        // also look for this phenotype reversing

    /*
     NEED:
        -- constructors with both types of module state sampling outcome
            -- read new pins & stops file
            -- compare and match incoming attracvtors & pins, fill named phenoytpes and pins file
     
        --
     */
    DM_Phenotype(GSE_line,vector<string> DM_Node_Names){
        istringstream iss(GSE_line);
        if (GSE_line[0]!='\t') {
            printf("Module phenotype line expected;  place tab in front of phenotype name in line %s. Please resolve before moving on with module analysis.\n",GSE_line); exit(1);
        }
        iss >> DM_Phen_Name;
        
        vector<string> word = Extract_String_Between_Delimiters_On_Same_Level(GSE_line, '(',')');
        for(unsigned int k=0;k<word.size();k++)
            DM_Phen_PIN.push_back(extract_sequence_of_PINS(word[k],DM_Node_Names));

    }
    
    string extract_sequence_of_PINS(string inparantheses,vector<string> DM_Node_Names){
        // wrong place; can't access module node name list!!!
        
        /// not edited yet!!!
        
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
                    if (in_id>1){
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
};

class Dynamical_Module {
public:
    string DM_Name,Coupled_Name;
    unsigned int DM_Phenotype_NR,DM_Node_NR;
    vector<string> DM_Node_Names;
    vector<DM_Phenotype> DM_Phen;

    Dynamical_Module(const char modname, const char modelname,Boolean_Dynamics *D){
        DM_Name = modname;
        Coupled_Name = modelname;
        fill_DM_Node_Names(D);
        create_phenotypes_based_on_PINS(D);
        match_PINS_to_attractors(D);
    }
    
    fill_DM_Node_Names(Boolean_Dynamics *D){
        DM_Node_NR=D->BRN->N;
        for(unsigned int s=1;s<=D->BRN->N;s++)
            DM_Node_Names.push_back(D->BRN->MDAT->Node_name[s]);
    }
    
    create_phenotypes_based_on_PINS(Boolean_Dynamics *D){
        char fn[400];
        FILE *f;
        sprint(fn,"%s/%s/Module_Phenotype_Info.txt",MODEL_Directory_system,Coupled_Name.c_str());
        f=fopen(fn,"r");
        if(f==NULL) {printf("Module_Phenotype_Info.txt File is not found in Model directory %s \n \tCurrent model folder: %s\n%s\n",Coupled_Name.c_str(),MODEL_Directory_system,Coupled_Name.c_str()); exit(1);}
        
        while(fgets(GSE_line,MAX_LINE,f)!=NULL){
            istringstream iss(GSE_line);
            if ((GSE_line[0]!='\t')&&(GSE_line[0]!='\n')&&(GSE_line[0]!='\r')&&(GSE_line[0]!=' ') {
                string mname;
                iss >> mname;
                 if (mname == DM_Name) {
                     unsigned int mattrnr = 0;
                     iss >> mattrnr;
                     if (iss.fail() == true){
                         printf("After module name %s please specify number of expected attractors (integer expected)\n",DM_Name.c_str()); exit(1);
                     }
                     if(mattrnr!=D->Attractor_NR){
                         printf("Module %s has %d attractors when isolated from the coupled network. Module_Phenotype_Info.txt file is expecting %d attractors. Please resolve before moving on with module analysis.\n",DM_Name.c_str(),D->Attractor_NR,mattrnr); exit(1);
                     }
                     else {
                         DM_Phenotype_NR = D->Attractor_NR;
                         for(unsigned int ia=1;ia<=DM_Phenotype_NR;ia++)
                             if(fgets(GSE_line,MAX_LINE,f)!=NULL) {
                                 DM_Phen.push_back(DM_Phenotype(GSE_line),DM_Node_Names);
                             }
                             else {printf("Could not read expected phenotype with index %f for module %s. Please resolve before moving on with module analysis.\n",ia, DM_Name.c_str()); exit(1);}
                     }
                 }
            }
        }
        
        fclose(f);
    }
};

/* DynamicalModularity_h */
