//
//  Interactive_timecourse.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 1/22/18.
//  Copyright © 2018 Regan_Group. All rights reserved.
//

#ifndef Interactive_timecourse_h
#define Interactive_timecourse_h

#include "SDL2/SDL.h"
#include "SDL2_ttf/SDL_ttf.h"
#define T_W 10
#define G_W 12
#define Name_W 30
#define T_MAX 150

using namespace std;
SDL_Window* window;
SDL_Renderer* renderer;

int Setup_Timecourse_window(unsigned int N_row,unsigned int N_t) {
    window = SDL_CreateWindow("Timecourse!", 5, 5, Name_W + T_W * N_t, G_W * N_row, 0);
    
    if (window == NULL){
        printf("Failed to create window : %s\n ",SDL_GetError());
        return -1;
    }
    
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    
    if (renderer == nullptr){
        printf("Failed to create renderer : %s\n", SDL_GetError());
        return -1;
    }
    
    SDL_RenderSetLogicalSize(renderer, Name_W + T_W * N_t, G_W * N_row);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);
    return (1);
}


void drawText(char* string,
              int size,
              int x, int y,
              unsigned char fR, unsigned char fG, unsigned char fB){
    
    char fn[300];
    sprintf(fn, "/Users/%s/Dropbox/_CODE/freefont-ttf/sfd/arial.ttf",whoamI);
    TTF_Font* Sans = TTF_OpenFont(fn, G_W);
   // TTF_Font* Sans = TTF_OpenFont("Arial.ttf", G_W);
    SDL_Color textColor = { fR, fG, fB };
    SDL_Surface * surface = TTF_RenderText_Solid(Sans,string,textColor);
    SDL_Texture * texture = SDL_CreateTextureFromSurface(renderer, surface);
    
    int texW = 0;
    int texH = 0;
    SDL_QueryTexture(texture, NULL, NULL, &texW, &texH);
    SDL_Rect dstrect = { x, y, texW, texH };
    SDL_RenderCopy(renderer, texture, NULL, &dstrect);
    
    SDL_DestroyTexture(texture);
    SDL_FreeSurface(surface);
    TTF_CloseFont(Sans);
}

int node_ID_from_poz(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,int poz_g){
    
    int poz=1;
    for(int m=1;m<=PD->BRN->Module_Order->N;m++){
        if(PD->BRN->Module_Order->A[m]>0){
            for(int nn=1;nn<=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N;nn++){
                int i=(int)PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->A[nn];
                poz++;
                if(poz==poz_g) return(i);
            }
        }
        else{
            int i=(int)PD->BRN->Module_Order->B[m];
            poz++;
            if(poz==poz_g) return(i);
        }
        poz++;
    }
    return(0);
}

char *get_inputchange_info(int id_pulse,
                           int a_st,int a_st_k,
                           MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                           int pulse_l,int delay,
                           int poz_g,int t_now){
    char *whatsup;
    int nodeID=node_ID_from_poz(PD,poz_g);
    int *inputnodes;
    int *inputstates_before, *inputstates_after;
    int SHOW_Init_state = 50;
    
    whatsup=NULL;
    
    if (nodeID>0) {
        PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
        
        int s_start=PD->Coupled->s[id_pulse-1];
        int id_GF=PD->BRN->MDAT->get_node_ID("GF");
        int s_start_GF=PD->Coupled->s[id_GF-1];
        
        for (int t=1;t<t_now; t++) {
            if (t>=SHOW_Init_state+delay) PD->Coupled->s[id_pulse-1]=48+(49-s_start);
            if (t>=SHOW_Init_state+delay + pulse_l) {
                PD->Coupled->s[id_pulse-1]=s_start;
                PD->Coupled->s[id_GF-1]=s_start_GF;
            }
            PD->Coupled->synchronously_update_all_Gates();
        }
        
        inputnodes=new int[PD->Coupled->BRN->BN->node[nodeID].k_in+1];
        inputstates_before=new int[PD->Coupled->BRN->BN->node[nodeID].k_in+1];
        
        snode::slinklist::slink *sc2;
        sc2=PD->Coupled->BRN->BN->node[nodeID].Li.first_link_pt;
        int ii=1;
        while(sc2!=NULL){
            inputnodes[ii]=(int)sc2->sneighbor;
            inputstates_before[ii]= PD->Coupled->s[sc2->sneighbor-1];
            ii++;
            sc2=sc2->next(sc2);
        }
        // timestep that got us to NOW
        
        if (t_now>=SHOW_Init_state+delay) PD->Coupled->s[id_pulse-1]=48+(49-s_start);
        if (t_now>=SHOW_Init_state+delay + pulse_l) {
            PD->Coupled->s[id_pulse-1]=s_start;
            PD->Coupled->s[id_GF-1]=s_start_GF;
        }
        PD->Coupled->synchronously_update_all_Gates();

        inputstates_after=new int[PD->Coupled->BRN->BN->node[nodeID].k_out+1];
        sc2=PD->Coupled->BRN->BN->node[nodeID].Li.first_link_pt;
        ii=1;
        while(sc2!=NULL){
            inputstates_after[ii]= PD->Coupled->s[sc2->sneighbor-1];
            ii++;
            sc2=sc2->next(sc2);
        }
        
        whatsup=new char[2000];
        sprintf(whatsup, "%s from last timestep -   ",PD->Coupled->BRN->MDAT->Node_name[nodeID]);
        sc2=PD->Coupled->BRN->BN->node[nodeID].Li.first_link_pt;
        ii=1;
        while(sc2!=NULL){
            if((inputstates_after[sc2->sneighbor]==0 ) &&( inputstates_before[sc2->sneighbor] ==1))
                sprintf(whatsup,"%s %s turned OFF  ",whatsup,PD->Coupled->BRN->MDAT->Node_name[sc2->sneighbor]);
            if((inputstates_after[sc2->sneighbor]==1 ) &&( inputstates_before[sc2->sneighbor] ==0))
                sprintf(whatsup,"%s %s turned ON ",whatsup,PD->Coupled->BRN->MDAT->Node_name[sc2->sneighbor]);
            sc2=sc2->next(sc2);
        }
    }
    return(whatsup);
}

void Draw_Timecourse(MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                     unsigned int N_row,unsigned int N_t,
                     us_shortintpointer *s,charpointer *nodename){
    SDL_Rect r;
   
    SDL_Init(SDL_INIT_EVERYTHING);
    TTF_Init();

    if (Setup_Timecourse_window(N_row,N_t) == -1) exit(1);

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);
    
    for (int i = 1; i <= N_t; i++) {
        for (int j = 1; j <= N_row; j++) {
            r.x = Name_W + T_W * i;
            r.y = G_W * j;
            r.w = T_W;
            r.h = G_W;
            if(s[i][j]==48) SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255); // blue
            if(s[i][j]==49) SDL_SetRenderDrawColor(renderer, 255, 128, 0, 255); // orange
            SDL_RenderFillRect(renderer, &r);
        }
    }
    

    for (unsigned long int j = 1; j <= N_row; j++) {
        drawText(nodename[j],G_W,1,G_W * j,255,0,0);
    }

    SDL_RenderPresent(renderer);
    
}

void animate_timecourse(const char nodename[],
                        int a_st,int a_st_k,
                        MODULAR_Boolean_Dynamics_PAIRED_TRACKER *PD,
                        int pulse_l,int delay){
 //   long int id_new;
    us_shortintpointer *s;
    charpointer *nodenamelist;
    int SHOW_Init_state = 50;
    unsigned int N_row = 0;
    
    for(int m=1;m<=PD->BRN->Module_Order->N;m++){
        if(PD->BRN->Module_Order->A[m]>0){
            N_row +=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N;
        }
        else{ N_row++;}
        N_row++;
    }
    
    nodenamelist=new charpointer[N_row+1];
    for (int j=1; j<=N_row; j++) {
        nodenamelist[j]=new char[50];
        sprintf(nodenamelist[j],"");
    }
    
    int poz=1;
    for(int m=1;m<=PD->BRN->Module_Order->N;m++){
        if(PD->BRN->Module_Order->A[m]>0){
            for(int nn=1;nn<=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N;nn++){
                unsigned long int i=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->A[nn];
                sprintf(nodenamelist[poz],"%s",PD->BRN->MDAT->Node_name[i]);
                poz++;
            }
        }
        else{
            unsigned long int i=PD->BRN->Module_Order->B[m];
            sprintf(nodenamelist[poz],"%s",PD->BRN->MDAT->Node_name[i]);
            poz++;
        }
        poz++;
    }
    
    s=new us_shortintpointer[T_MAX+1];
    for(int m=1;m<=T_MAX;m++) {
        s[m]=new unsigned short int[N_row+1];
        for (int j=1; j<=N_row; j++) {
            s[m][j]=-1;
        }
    }
    //int T_MAX=650,TW=5,NameBuff=35;
    
    int id_pulse=PD->BRN->MDAT->get_node_ID(nodename);
    
    PD->Coupled->set_state(PD->Coupled->Attractor_valleys[a_st]->Basin_Stored[PD->Coupled->Attractor_valleys[a_st]->Attractor[a_st_k]]->s);
    
    int s_start=PD->Coupled->s[id_pulse-1];
    int id_GF=PD->BRN->MDAT->get_node_ID("GF");
    int s_start_GF=PD->Coupled->s[id_GF-1];
    
    for (int t=1;t<=T_MAX; t++) {
        if (t>=SHOW_Init_state+delay) PD->Coupled->s[id_pulse-1]=48+(49-s_start);
        if (t>=SHOW_Init_state+delay + pulse_l) {
            PD->Coupled->s[id_pulse-1]=s_start;
            PD->Coupled->s[id_GF-1]=s_start_GF;
        }
        PD->Coupled->synchronously_update_all_Gates();
        int poz=1;
        
        for(int m=1;m<=PD->BRN->Module_Order->N;m++){
            if(PD->BRN->Module_Order->A[m]>0){
                for(int nn=1;nn<=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->N;nn++){
                    unsigned long int i=PD->BRN->Draw_in_Modules[PD->BRN->Module_Order->A[m]]->A[nn];
                    s[t][poz]=PD->Coupled->s[i-1];
                    poz++;
                }
            }
            else{
                unsigned long int i=PD->BRN->Module_Order->B[m];
                s[t][poz]=PD->Coupled->s[i-1];
                poz++;
            }
            poz++;
        }
    }
    
    Draw_Timecourse(PD,N_row,(unsigned int)T_MAX, s,nodenamelist);
    
   // new_id= PD->Coupled->find_attractor_basin_for_state_without_altering_tracking_if_in_know_attractor(PD->Coupled->s);
    
    
    
    SDL_Event event;
    char *whatsup;
    
    /* A bool to check if the program has exited */
    int quit = 0;
    
    /* While the program is running */
    while (!quit){
        /* Check for new events */
        while(SDL_PollEvent(&event)){
            /* If a quit event has been sent */
            if (event.type == SDL_QUIT) {
                /* Quit the application */
                quit = 1;
            }
            
            /* If a button on the mouse is pressed. */
            if (event.type == SDL_MOUSEBUTTONDOWN){
                /* If the left button was pressed. */
                if (event.button.button == SDL_BUTTON_LEFT){
                    whatsup=NULL;
                    whatsup=get_inputchange_info(id_pulse, a_st, a_st_k,PD,pulse_l, delay,
                                             event.motion.y / G_W +1,(event.motion.x - Name_W) / T_W );
                    if(whatsup!=NULL){
                        printf("%s\n",whatsup);
                        SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_INFORMATION,"Gene report",whatsup, window);
                        delete[] whatsup;whatsup=NULL;
                    }
                }
            }
        }
    }
    
    getchar();
    TTF_Quit();
    SDL_DestroyWindow(window);
    SDL_Quit();
}

#endif /* Interactive_timecourse_h */
