//
//  Virtual_Exp_classes.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 10/23/19.
//  Copyright © 2019 Regan_Group. All rights reserved.
//

#ifndef Virtual_Exp_classes_h
#define Virtual_Exp_classes_h

class  Phenotypes_for_Stats {
public:
    unsigned long int WT_cell_cycle = 0;
    bool CellCycleErrors = 0, Irrev = 0 , Rev = 0;
    unsigned long int DNA_Replication_START, DNA_Replication_DONE,Metaphase_START,SAC_Cleared,Cytokinesis;
    std::vector<std::string> Irreversible_Do_Not_Count_CellCycle_Errors;
    std::vector<bool> Irreversible_Do_Not_Count_CellCycle_Errors_LimitCycle;
    
    std::vector<std::string> Irreversible_STATS;
    std::vector<longint_ARRAY_pair> Irreversible_Do_Not_Count_CellCycle_Errors_NodeStates;
    std::vector<longint_ARRAY_pair> Irreversible_STATS_NodeStates;
    std::vector<bool> Irreversible_STATS_LimitCycle;
    std::vector<unsigned short int> Irreversible_STATS_GO_ON;
    
    std::vector<std::string> Reversible_Transitions_1,Reversible_Transitions_2;
    std::vector<longint_ARRAY_pair> Reversible_Transitions_1_NodeStates,Reversible_Transitions_2_NodeStates;
    std::vector<bool> Reversible_Transitions_1_LimitCycle,Reversible_Transitions_2_LimitCycle;

    Phenotypes_for_Stats (){}
    
} Phen_STATs;

struct Cell_cycle_stats {
    unsigned int    N_cycles,
    N_endoredupl_T,N_endoredupl_G2,N_endoredupl_Aneup,Cell_cycle_stats_in_run;
    
    unsigned long int Time_G1,Time_R,Time_G2,Time_M,Time_Telo;
    unsigned long int N_G1,N_G2,N_R,N_M,N_T;
    unsigned long long int   T_live,MAX_INTERVAL=0;
    unsigned int Dist_G1[MAX_ARREST+1],
    Dist_G2[MAX_ARREST+1],
    Dist_Mitosis[MAX_ARREST+1],
    Dist_Replication[MAX_ARREST+1],
    Dist_Telophase[MAX_ARREST+1];
};

struct Irrev_Phenotyope_stats{
    unsigned long int N_events;
    unsigned long long int T_live;
    double av_time_to_event,SD_time_to_event,Half_life;
    unsigned int Dist_Time_to_event[MAX_ARREST+1];
    
};

struct Rev_Transition_stats{
    unsigned long int N_event_to_1,N_event_to_2;
    unsigned long long int T_live;
    double av_window_in_1,SD_window_in_1,Half_life_1;
    double av_window_in_2,SD_window_in_2,Half_life_2;
    unsigned int T_in_1,T_in_2;
    unsigned int Dist_Time_in_1[MAX_ARREST+1],Dist_Time_in_2[MAX_ARREST+1];
};

struct Cell_cycle_Irrev_Rev_stats_in_run {
    Cell_cycle_stats CC;
    std::vector<Irrev_Phenotyope_stats> IRev;
    std::vector<Rev_Transition_stats> Rev;
};


class Environment_and_hits{
public:
    double *p_inputs,*p_hits;
    longint_ARRAY_pair *BKR_Hits;
    longint_ARRAY_pair *Scan_Hits;
    
    Environment_and_hits(unsigned int ni){
        p_inputs=new double[ni+1];
        for(int j=0;j<=ni;j++)
            p_inputs[j]=-1;
    }
    ~Environment_and_hits(){
        delete[] p_inputs;p_inputs=NULL;
        delete[] p_hits;p_hits=NULL;
        delete BKR_Hits;BKR_Hits=NULL;
        delete Scan_Hits;Scan_Hits=NULL;
    }
};


#endif /* Virtual_Exp_classes_h */
