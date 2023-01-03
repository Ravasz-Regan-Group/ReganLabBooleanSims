//
//  Ingrid_NGS_Li_Lu_Ki_LPS_Tie2.h
//  Erzso_Code_Central
//
//  Created by Erzsó Ravasz Regan on 6/26/17.
//  Copyright © 2017 Beth Israel Deaconess Medical Center. All rights reserved.
//

#ifndef Ingrid_NGS_Li_Lu_Ki_LPS_Tie2_h
#define Ingrid_NGS_Li_Lu_Ki_LPS_Tie2_h

struct Gene {
    unsigned long int ID;
    char ENSMUSG[100];
    char Symbol[100];
    unsigned long int ENTREZ;
    char GeneName[1000];
};

struct Gene_Results {
    longint_ARRAY  *gene_IDs;
    double *FC,*padj;
};

unsigned int get_NR_up(Gene_Results *gr){
    unsigned int n=0;
    for(unsigned int i=1;i<=gr->gene_IDs->N;i++)
        if(gr->FC[i]>0) n++;
    return(n);
}

unsigned int get_NR_down(Gene_Results *gr){
    unsigned int n=0;
    for(unsigned int i=1;i<=gr->gene_IDs->N;i++)
        if(gr->FC[i]<0) n++;
    return(n);
}

unsigned long int N_GENES = 0;
typedef Gene * GeneList;
GeneList *G_names;


typedef Gene_Results * Gene_Results_List;
typedef Gene_Results_List * Gene_Results_Matrix;

struct Exp_Design{
    unsigned int NT,NC,NG;
    charpointer *Tissues, *Conditions,*Genotypes;
};

unsigned long int get_gene_ID(const char *gE){
    for(unsigned long int i=1;i<=N_GENES;i++){
        if(strcmp (G_names[i]->ENSMUSG,gE) == 0 ) return (i);
        
    }
    printf("gene %s not found!\n\n\n",gE);
    //printf("__%s\n__%s__\n",gE,G_names[1]->ENSMUSG);getchar();
    return(0);
}

Gene_Results *get_gene_list(const char *fn ,double FC,double padj){
    FILE *f;
    Gene_Results *GR;
    
    GR=new Gene_Results;
    GR->gene_IDs=new longint_ARRAY();
    
    double fc = 0.0, p = 1.0;
    unsigned long int  k;
    char bla1[300];
    std::string ss;
    char attr_list[100][5000];
    unsigned long int Nf;
    
    sprintf(bla1,"%s.csv",fn);
    f=fopen(bla1,"r");
    printf("%s\n",fn);//getchar();
    fgets(GSE_line,MAX_LINE,f);
    
    for (unsigned long int i=1; i<=N_GENES; i++) {
        if ( fgets(GSE_line,MAX_LINE,f) != NULL ){
           // printf("%s\n",GSE_line);
            Nf=separate_line(GSE_line,attr_list,',',5000);
            
           // fscanf(f,"%s\t%lg\t%lg\t%lg\t%lg\t%s\t%s",bla1,&m,&fc,&sd,&s,bla2,bla3);
            ss=attr_list[1];
            ss.erase(
                    remove( ss.begin(), ss.end(), '\"' ),
                    ss.end()
                    );

            if(attr_list[3][0]!='N') fc=atof(attr_list[3]); else fc=0;
            if(attr_list[7][0]!='N') p=atof(attr_list[7]); else p=1;
            
           // printf("%s\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",ss.c_str(),atof(attr_list[2]),fc,
           //        atof(attr_list[4]),atof(attr_list[5]),atof(attr_list[6]),padj);
           // getchar();
            
            if((p<=padj)&&(fabs(fc)>=FC)){
                unsigned long int ID = get_gene_ID(ss.c_str());
                if(ID>0) GR->gene_IDs->add_element(ID);
            }
        }
    }
    fclose(f);
    printf("%ld genes pass the test\n",GR->gene_IDs->N);//getchar();
    
    GR->FC=new double[GR->gene_IDs->N+1];
    GR->padj=new double[GR->gene_IDs->N+1];
    
    f=fopen(bla1,"r");
    fgets(GSE_line,MAX_LINE,f);
    k=0;
    for (unsigned long int i=1; i<=N_GENES; i++) {
        if ( fgets(GSE_line,MAX_LINE,f) != NULL ){
            Nf=separate_line(GSE_line,attr_list,',',5000);
            if(attr_list[3][0]!='N') fc=atof(attr_list[3]); else fc=0;
            if(attr_list[7][0]!='N') p=atof(attr_list[7]); else p=1;

            if((p<=padj)&&(fabs(fc)>=FC)){
                k++;
                GR->FC[k]=fc;
                GR->padj[k]=p;
            }
        }
    }
    fclose(f);
    
    return(GR);
}

void export_subset_for_heatmap(Gene_Results *gr,const char path[],const char fn[],int PG,
                               int Cin,int Gin,int t,int excl){
                                    // value of -1 means ALL in, 0 or avobe pics the one to include, except if excl=1 for exluding 1 tissue...
    FILE *f,*f_aeag,*f_aetg,*f_ceag,*f_cetg;
    char fname[300];
    std::string ss;
    char attr_list[100][5000];
    unsigned long int Nf,ID,k_gene;
    double *FCord, FC_top_cut=0.0;
    int Exp_include[100];
    
    if(PG<gr->gene_IDs->N){
        FCord=new double[gr->gene_IDs->N+1];
        for(unsigned int i=0;i<=gr->gene_IDs->N;i++) FCord[i] = fabs(gr->FC[i]);
        hpsort((int)gr->gene_IDs->N, FCord);
        FC_top_cut=FCord[gr->gene_IDs->N-PG+1];
    }
    for (unsigned int i=0;i<100;i++) Exp_include[i]=0;
    for (unsigned int i=1;i<=72;i++) Exp_include[i]=1;
    if(Cin==0){ // no LPS wanted!
        for (unsigned int i=1;i<=6;i++) Exp_include[i]=0;
        for (unsigned int i=13;i<=18;i++) Exp_include[i]=0;
        for (unsigned int i=25;i<=30;i++) Exp_include[i]=0;
        for (unsigned int i=37;i<=42;i++) Exp_include[i]=0;
        for (unsigned int i=49;i<=54;i++) Exp_include[i]=0;
        for (unsigned int i=61;i<=66;i++) Exp_include[i]=0;
    }
    if(Cin==1){ // no Control wanted!
        for (unsigned int i=7;i<=12;i++) Exp_include[i]=0;
        for (unsigned int i=19;i<=24;i++) Exp_include[i]=0;
        for (unsigned int i=31;i<=36;i++) Exp_include[i]=0;
        for (unsigned int i=43;i<=48;i++) Exp_include[i]=0;
        for (unsigned int i=55;i<=60;i++) Exp_include[i]=0;
        for (unsigned int i=67;i<=72;i++) Exp_include[i]=0;
    }
    if(Gin==0){ // no WT wanted!
        for (unsigned int i=13;i<=24;i++) Exp_include[i]=0;
        for (unsigned int i=37;i<=54;i++) Exp_include[i]=0;
        for (unsigned int i=61;i<=72;i++) Exp_include[i]=0;
    }
    if(Gin==1){ // no Tie2 wanted!
        for (unsigned int i=1;i<=12;i++) Exp_include[i]=0;
        for (unsigned int i=25;i<=36;i++) Exp_include[i]=0;
        for (unsigned int i=55;i<=60;i++) Exp_include[i]=0;
        }
    if(t==0){
        if(excl==0) {for (unsigned int i=25;i<=72;i++) Exp_include[i]=0;}
        else {for (unsigned int i=1;i<=24;i++) Exp_include[i]=0;}
    }
    if(t==1){
        if(excl==0) {
            for (unsigned int i=1;i<=24;i++) Exp_include[i]=0;
            for (unsigned int i=49;i<=72;i++) Exp_include[i]=0;
        }
        else {for (unsigned int i=25;i<=48;i++) Exp_include[i]=0;}
    }
    if(t==2){
        if(excl==0) {for (unsigned int i=1;i<=48;i++) Exp_include[i]=0;}
        else {for (unsigned int i=49;i<=72;i++) Exp_include[i]=0;}
    }
    
    sprintf(fname,"%s/_ALL.csv",path);
   // printf("%s\n",fname);getchar();
    f=fopen(fname, "r");
    sprintf(fname,"%s_Show_AllExp_AllSigGenes.csv",fn);         f_aeag=fopen(fname, "w");
    if(PG<gr->gene_IDs->N)
        sprintf(fname,"%s_Show_AllExp_TOPGenes.csv",fn);         f_aetg=fopen(fname, "w");
    sprintf(fname,"%s_Show_Contrast_AllSigGenes.csv",fn);  f_ceag=fopen(fname, "w");
    if(PG<gr->gene_IDs->N)
        sprintf(fname,"%s_Show_Contrast_TOPGenes.csv",fn);     f_cetg=fopen(fname, "w");
    
    fgets(GSE_line,MAX_LINE,f);
    Nf=separate_line(GSE_line,attr_list,',',5000);
    fprintf(f_aeag,"%s",GSE_line);
    fprintf(f_aetg,"%s",GSE_line);
    
    fprintf(f_ceag,"%s",attr_list[1]);
    fprintf(f_cetg,"%s",attr_list[1]);
    
    for(unsigned int k=1;k<=72;k++)
        if(Exp_include[k]==1){
            fprintf(f_ceag, "\t%s",attr_list[k+1]);
            fprintf(f_cetg, "\t%s",attr_list[k+1]);
        }
    fprintf(f_ceag,"\n");
    fprintf(f_cetg,"\n");
    
    for (unsigned long int i=1; i<=N_GENES; i++) {
        if ( fgets(GSE_line,MAX_LINE,f) != NULL ){
            Nf=separate_line(GSE_line,attr_list,',',5000);
            ID = get_gene_ID(attr_list[1]);
            k_gene=gr->gene_IDs->check_for_element(ID);
            if(k_gene>0){
                fprintf(f_aeag,"%s",GSE_line);
                fprintf(f_ceag,"%s",attr_list[1]);
                if(PG<gr->gene_IDs->N)
                     if(fabs(gr->FC[k_gene])>=FC_top_cut){
                         fprintf(f_aetg,"%s",GSE_line);
                         fprintf(f_cetg,"%s",attr_list[1]);
                    }
                for(unsigned int k=1;k<=72;k++)
                    if(Exp_include[k]==1){
                        fprintf(f_ceag, "\t%s",attr_list[k+1]);
                        if(PG<gr->gene_IDs->N)
                            if(fabs(gr->FC[k_gene])>=FC_top_cut)
                                 fprintf(f_cetg, "\t%s",attr_list[k+1]);
                    }
                fprintf(f_ceag,"\n");
                fprintf(f_cetg,"\n");
            }
        }
    }
    fclose(f);
    fclose(f_aeag);
    fclose(f_ceag);
    fclose(f_aetg);
    fclose(f_cetg);
}

void process_tissue_cuts(Exp_Design ED,const char path[],double FC, double padj){
    Gene_Results_Matrix *GR;
    GR=new Gene_Results_Matrix[3];
    GR[0] =new Gene_Results_List[7+1]; // Liver
    GR[1] =new Gene_Results_List[7+1]; // Lung
    GR[2] =new Gene_Results_List[7+1]; // Kidney
    
    int PlotGene=200;
    char fn[300];

    for(int i=0;i<ED.NT;i++){
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Genotypes",path,ED.Tissues[i],ED.Tissues[i],
                ED.Conditions[0],ED.Conditions[1]);
        GR[i][0]=get_gene_list(fn,FC,padj);
        //   export_subset_for_heatmap(Gene_Results *gr, char path[],char fn[],int c,int t,int g)
        export_subset_for_heatmap(GR[i][0],path,fn,PlotGene,-1,-1,i,0);
        // int Cin,int Gin,int Lungin,int Liverin,int Kidneyin
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Tissues[i],ED.Tissues[i],
                ED.Conditions[0],ED.Conditions[1],ED.Genotypes[0]);
        GR[i][1]=get_gene_list(fn,FC,padj);
        export_subset_for_heatmap(GR[i][1],path,fn,PlotGene,-1,0,i,0);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Tissues[i],ED.Tissues[i],
                ED.Conditions[0],ED.Conditions[1],ED.Genotypes[1]);
        GR[i][2]=get_gene_list(fn,FC,padj);
        export_subset_for_heatmap(GR[i][2],path,fn,PlotGene,-1,1,i,0);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Conditions",path,ED.Tissues[i],ED.Tissues[i],
                ED.Genotypes[0],ED.Genotypes[1]);
        GR[i][3]=get_gene_list(fn,FC,padj);
        export_subset_for_heatmap(GR[i][3],path,fn,PlotGene,-1,-1,i,0);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Tissues[i],ED.Tissues[i],
                ED.Genotypes[0],ED.Genotypes[1],ED.Conditions[0]);
        GR[i][4]=get_gene_list(fn,FC,padj);
        export_subset_for_heatmap(GR[i][4],path,fn,PlotGene,0,-1,i,0);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Tissues[i],ED.Tissues[i],
                ED.Genotypes[0],ED.Genotypes[1],ED.Conditions[1]);
        GR[i][5]=get_gene_list(fn,FC,padj);
        export_subset_for_heatmap(GR[i][5],path,fn,PlotGene,1,-1,i,0);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s_over_%s",path,ED.Tissues[i],ED.Tissues[i],
                ED.Genotypes[0],ED.Genotypes[1],ED.Conditions[1],ED.Conditions[0]);
        GR[i][6]=get_gene_list(fn,FC,padj);
        export_subset_for_heatmap(GR[i][6],path,fn,PlotGene,-1,-1,i,0);
        
    }
}

Gene_Results_Matrix *process_Genotype_cut(Exp_Design ED,const char path[],double FC, double padj,int exp_OK,int textok){
    Gene_Results_Matrix *GR;
    GR=new Gene_Results_Matrix[3];
    GR[0] =new Gene_Results_List[10];
            // LPS vs. control all for 0,
            // 1-3 : -- 3 tissues for index
            // 4-6 : LPS vs. control : 3 differential responses
            // 6-9: 3 tissue pair comparisons
    GR[1] =new Gene_Results_List[10];
            // same for Tie2
    int PlotGene=200;
    
    char fn[300];
    
    for(int i=0;i<ED.NG;i++){
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Tissues",path,ED.Genotypes[i],ED.Genotypes[i],
                ED.Conditions[0],ED.Conditions[1]); // e.g., WT_Control_to_LPS_all_Tissues.csv
        GR[i][0]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][0],path,fn,PlotGene,-1,0,-1,0);
        
        
        // single-tissue LPS vs. ct.
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Genotypes[i],ED.Genotypes[i],
                ED.Conditions[0],ED.Conditions[1],ED.Tissues[0]);
                        //e.g., WT_Control_to_LPS_in_Kidney.csv
        GR[i][1]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][1],path,fn,PlotGene,-1,0,0,0);
      
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Genotypes[i],ED.Genotypes[i],
                ED.Conditions[0],ED.Conditions[1],ED.Tissues[1]);
        GR[i][2]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][2],path,fn,PlotGene,-1,0,1,0);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Genotypes[i],ED.Genotypes[i],
                ED.Conditions[0],ED.Conditions[1],ED.Tissues[2]);
        GR[i][3]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][3],path,fn,PlotGene,-1,0,2,0);
        
        // LPS vs. ct. differential in tissues
        
         
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s_over_%s",path,ED.Genotypes[i],ED.Genotypes[i],
                 ED.Conditions[0],ED.Conditions[1],ED.Tissues[0],ED.Tissues[1]);
                    // WT_Control_to_LPS_in_Liver_over_Kidney.csv
        GR[i][4]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][4],path,fn,PlotGene,-1,0,2,1);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s_over_%s",path,ED.Genotypes[i],ED.Genotypes[i],
                ED.Conditions[0],ED.Conditions[1],ED.Tissues[0],ED.Tissues[2]);
        GR[i][5]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][5],path,fn,PlotGene,-1,0,1,1);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s_over_%s",path,ED.Genotypes[i],ED.Genotypes[i],
                ED.Conditions[0],ED.Conditions[1],ED.Tissues[1],ED.Tissues[2]);
        GR[i][6]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][6],path,fn,PlotGene,-1,0,0,1);
        
        // difference across tissues (all treatments)
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Conditions",path,ED.Genotypes[i],ED.Genotypes[i],
                ED.Tissues[1],ED.Tissues[0]);
                    // WT_Liver_to_Kidney_all_Conditions.csv
        GR[i][7]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][7],path,fn,PlotGene,-1,0,2,1);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Conditions",path,ED.Genotypes[i],ED.Genotypes[i],
                ED.Tissues[0],ED.Tissues[2]);
        GR[i][8]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][8],path,fn,PlotGene,-1,0,1,1);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Conditions",path,ED.Genotypes[i],ED.Genotypes[i],
                ED.Tissues[1],ED.Tissues[2]);
        GR[i][9]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][9],path,fn,PlotGene,-1,0,0,1);
        
    }
    
    if(textok) sprintf(fn,"%s/Genotype_summary_text.eps", path);
    else sprintf(fn,"%s/Genotype_summary_raw.eps", path);
    FILE *f;
    f=fopen(fn,"w");
    
    COLOR_SCHEME=GYR;
    setcolors_PAJEK();
    int XMAX,XMIN,YMAX, YMIN;
    
    XMAX=600;
    XMIN=0;
    YMAX=600;
    YMIN=0;
    double scalef = XMAX/4000.;
    
    EpsInit(f,-10,-15,XMAX-XMIN+10,YMAX-YMIN+10);
    
    EpsSetFont(f,"Times-Roman",9);
    double dim=1.5;
    EpsSetLinewidth(f, scalef*10);

    //WT
  //  EpsSetRgb(f,0.9,0.9,0.9);
  //  EpsFillCircle(f,XMAX/4.,YMAX/4.,scalef*GR[0][0]->gene_IDs->N);
  //  EpsSetRgb(f,0.7,0.7,0.7);
  //  EpsFillCircle(f,XMAX/4.,YMAX/4.,scalef*get_NR_down(GR[0][0]));
    printf("WT: %u U %u  D \n",get_NR_up(GR[0][0]),get_NR_down(GR[0][0]));
 
    longint_ARRAY *ao,*bo,*co1,*ab,*ac,*bc,*abc;
   
    // Liver
    EpsSetRgb(f,0,219/255.,152/255.);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*GR[0][1]->gene_IDs->N*sin(-M_PI/3.-M_PI/20.),
                                        YMAX/4.+ scalef*GR[0][1]->gene_IDs->N*cos(-M_PI/3.-M_PI/20.) ,
                                        XMAX/4.+ scalef*GR[0][1]->gene_IDs->N*sin(-M_PI/3.+M_PI/20.),
                                        YMAX/4.+ scalef*GR[0][1]->gene_IDs->N*cos(-M_PI/3.+M_PI/20.) );
    EpsSetRgb(f,0,(219/255.)/dim,2*(152/255.)/dim);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*get_NR_down(GR[0][1])*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[0][1])*cos(-M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*get_NR_down(GR[0][1])*sin(-M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[0][1])*cos(-M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[0],
                GR[0][1]->gene_IDs->N,get_NR_up(GR[0][1]),get_NR_down(GR[0][1]));
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*GR[0][1]->gene_IDs->N*sin(-M_PI/3.+M_PI/20.)+5,
                        YMAX/4.+ scalef*GR[0][1]->gene_IDs->N*cos(-M_PI/3.+M_PI/20.)+5,bla1);
    // Lung
    EpsSetRgb(f,1,147/255.,137/255.);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*GR[0][2]->gene_IDs->N*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*GR[0][2]->gene_IDs->N*cos(M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*GR[0][2]->gene_IDs->N*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*GR[0][2]->gene_IDs->N*cos(M_PI/3.+M_PI/20.) );
    EpsSetRgb(f,1./dim,(147/255.)/dim,(137/255.)/dim);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*get_NR_down(GR[0][2])*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[0][2])*cos(M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*get_NR_down(GR[0][2])*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[0][2])*cos(M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[1],
            GR[0][2]->gene_IDs->N,get_NR_up(GR[0][2]),get_NR_down(GR[0][2]));
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*GR[0][2]->gene_IDs->N*sin(M_PI/3.+M_PI/20.)+5,
                        YMAX/4.+ scalef*GR[0][2]->gene_IDs->N*cos(M_PI/3.+M_PI/20.)+5,bla1);
    
    // Kidney
    EpsSetRgb(f,1,138/255.,1);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*GR[0][3]->gene_IDs->N*sin(M_PI-M_PI/20.),
                    YMAX/4.+ scalef*GR[0][3]->gene_IDs->N*cos(M_PI-M_PI/20.) ,
                    XMAX/4.+ scalef*GR[0][3]->gene_IDs->N*sin(M_PI+M_PI/20.),
                    YMAX/4.+ scalef*GR[0][3]->gene_IDs->N*cos(M_PI+M_PI/20.) );
    EpsSetRgb(f,1/dim,(138/255.)/dim,1/dim);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*get_NR_down(GR[0][3])*sin(M_PI-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[0][3])*cos(M_PI-M_PI/20.) ,
                    XMAX/4.+ scalef*get_NR_down(GR[0][3])*sin(M_PI+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[0][3])*cos(M_PI+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[2],
            GR[0][3]->gene_IDs->N,get_NR_up(GR[0][3]),get_NR_down(GR[0][3]));
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*GR[0][3]->gene_IDs->N*sin(M_PI+M_PI/20.),
                        YMAX/4.+ scalef*GR[0][3]->gene_IDs->N*cos(M_PI+M_PI/20.)-10,bla1);
    
    // comparative
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[0][1]->gene_IDs,GR[0][2]->gene_IDs,GR[0][4]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Li_o=%ld   Lung_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Lu=%ld   Li_and_Diff=%ld   Lu_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef*10);
    EpsSetRgb(f,168/255.,1,3/255.);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/20.) ,
                    XMAX/4.+ scalef*(ab->N+abc->N)*sin(+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(+M_PI/20.) );
    EpsDrawEllipseArc(f, XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      30+180/40.,150-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*abc->N*sin(-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-M_PI/20.) ,
                    XMAX/4.+ scalef*abc->N*sin(+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(+M_PI/20.) );
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N)*sin(+M_PI/20.) -25,
                        YMAX/4.+ scalef*(ab->N+abc->N)*cos(+M_PI/20.) + 5,bla1);
    
    EpsSetLinewidth(f, scalef*(double)GR[0][0]->gene_IDs->N*5);
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f, XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/3.+M_PI/20.),
                   YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/3.+M_PI/20.),
                   XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.),
                   YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.));
    EpsDrawLine(f, XMAX/4.+ scalef*(ab->N+abc->N)*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(M_PI/3.-M_PI/20.),
                    XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.));
    EpsSetRgb(f,168/255.,1,3/255.);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[1],ac->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.)+5,
                  YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.) + 5,bla1);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[0],bc->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.)+5,
                  YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.) + 5,bla1);
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f, XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(-M_PI/3.+M_PI/20.));
   
    EpsDrawLine(f, XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(M_PI/3.-M_PI/20.));
    
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc;delete abc;
   
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[0][1]->gene_IDs,GR[0][3]->gene_IDs,GR[0][5]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Li_o=%ld   Ki_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Ki=%ld   Li_and_Diff=%ld   Ki_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
   
    EpsSetLinewidth(f, scalef*10);
    EpsSetRgb(f,1,0,0);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.+M_PI/20.) );
    EpsDrawEllipseArc(f, XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      150+180/40.,270-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*abc->N*sin(-2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-2*M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*abc->N*sin(-2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-2*M_PI/3.+M_PI/20.) );
    EpsSetLinewidth(f, scalef*(double)GR[0][0]->gene_IDs->N*5);
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.+M_PI/20.) - 50,
                  YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.+M_PI/20.)+5,bla1);
    
    
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/3.-M_PI/20.),
                    XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.));
     EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI+M_PI/20.),
                    XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.));
    EpsSetRgb(f,1,0,0);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[2],ac->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.)-5,
                        YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.) - 5,bla1);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[0],bc->N);
    if(textok) EpsDrawString(f, 0,  XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.) - 15,
                          YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.) - 25,bla1);
    
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.),
                    XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(-M_PI/3.-M_PI/20.));
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.),
                    XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(-M_PI+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(-M_PI+M_PI/20.));
    
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc;delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[0][2]->gene_IDs,GR[0][3]->gene_IDs,GR[0][6]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Lu_o=%ld   Ki_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Lu_and_Ki=%ld   Lu_and_Diff=%ld   Ki_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef*10);
    EpsSetRgb(f,0,102/255.,1);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.+M_PI/20.) );
    EpsDrawEllipseArc(f, XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      270+180/40.,30-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*abc->N*sin(2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(2*M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*abc->N*sin(2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(2*M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.+M_PI/20.)+10,
                  YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.+M_PI/20.)+5,bla1);
    
    EpsSetLinewidth(f, scalef*(double)GR[0][0]->gene_IDs->N*5);
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N)*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(M_PI/3.+M_PI/20.),
                    XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.));
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI-M_PI/20.),
                    XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.));
    EpsSetRgb(f,0,102/255.,1);
    sprintf(bla1, "%ld NOT in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[2],ac->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.) - 5,
                        YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.) - 5,bla1);
    sprintf(bla1, "%ld NOT in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[1],bc->N);
    if(textok) EpsDrawString(f, 0,  XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.)+15,
                         YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.) + 5,bla1);
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.),
                    XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(M_PI/3.+M_PI/20.));
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.),
                    XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(-M_PI-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(-M_PI-M_PI/20.));
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc; delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[0][1]->gene_IDs,GR[0][2]->gene_IDs,GR[0][3]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Li_o=%ld   Ku_o=%ld   Ki_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Lu=%ld   Li_and_Ki=%ld   Lu_and_Ki=%ld   \n",ab->N,ac->N,bc->N);
    printf("All_3_Organs=%ld\n",abc->N);//getchar();
    EpsSetRgb(f, 0,0,0);
    EpsDrawCircle(f, XMAX/4., YMAX/4., scalef*abc->N);
    //getchar();
    
  
    
    
    EpsSetLinewidth(f, scalef*2);
    
    // Tie2
    //  EpsSetRgb(f,0.9,0.9,0.9);
    //  EpsFillCircle(f,XMAX/4.,YMAX/4.,scalef*GR[0][0]->gene_IDs->N);
    //  EpsSetRgb(f,0.7,0.7,0.7);
    //  EpsFillCircle(f,XMAX/4.,YMAX/4.,scalef*get_NR_down(GR[0][0]));
    printf("WT: %u U %u  D \n",get_NR_up(GR[1][0]),get_NR_down(GR[1][0]));
    EpsSetFont(f,"Times-Roman",5);
    
    // Liver
    EpsSetRgb(f,0,219/255.,152/255.);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*GR[1][1]->gene_IDs->N*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*GR[1][1]->gene_IDs->N*cos(-M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*GR[1][1]->gene_IDs->N*sin(-M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*GR[1][1]->gene_IDs->N*cos(-M_PI/3.+M_PI/20.) );
    EpsSetRgb(f,0,(219/255.)/dim,2*(152/255.)/dim);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*get_NR_down(GR[1][1])*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[1][1])*cos(-M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*get_NR_down(GR[1][1])*sin(-M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[1][1])*cos(-M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[0],
            GR[1][1]->gene_IDs->N,get_NR_up(GR[1][1]),get_NR_down(GR[1][1]));
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*GR[1][1]->gene_IDs->N*sin(-M_PI/3.+M_PI/20.)+5,
                  YMAX/4.+ scalef*GR[1][1]->gene_IDs->N*cos(-M_PI/3.+M_PI/20.)+5,bla1);
    // Lung
    EpsSetRgb(f,1,147/255.,137/255.);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*GR[1][2]->gene_IDs->N*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*GR[1][2]->gene_IDs->N*cos(M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*GR[1][2]->gene_IDs->N*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*GR[1][2]->gene_IDs->N*cos(M_PI/3.+M_PI/20.) );
    EpsSetRgb(f,1./dim,(147/255.)/dim,(137/255.)/dim);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*get_NR_down(GR[1][2])*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[1][2])*cos(M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*get_NR_down(GR[1][2])*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[1][2])*cos(M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[1],
            GR[1][2]->gene_IDs->N,get_NR_up(GR[1][2]),get_NR_down(GR[1][2]));
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*GR[1][2]->gene_IDs->N*sin(M_PI/3.+M_PI/20.)+5,
                  YMAX/4.+ scalef*GR[1][2]->gene_IDs->N*cos(M_PI/3.+M_PI/20.)+5,bla1);
    
    // Kidney
    EpsSetRgb(f,1,138/255.,1);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*GR[1][3]->gene_IDs->N*sin(M_PI-M_PI/20.),
                    YMAX/4.+ scalef*GR[1][3]->gene_IDs->N*cos(M_PI-M_PI/20.) ,
                    3*XMAX/4.+ scalef*GR[1][3]->gene_IDs->N*sin(M_PI+M_PI/20.),
                    YMAX/4.+ scalef*GR[1][3]->gene_IDs->N*cos(M_PI+M_PI/20.) );
    EpsSetRgb(f,1/dim,(138/255.)/dim,1/dim);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*get_NR_down(GR[1][3])*sin(M_PI-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[1][3])*cos(M_PI-M_PI/20.) ,
                    3*XMAX/4.+ scalef*get_NR_down(GR[1][3])*sin(M_PI+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_down(GR[1][3])*cos(M_PI+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[2],
            GR[1][3]->gene_IDs->N,get_NR_up(GR[1][3]),get_NR_down(GR[1][3]));
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*GR[1][3]->gene_IDs->N*sin(M_PI+M_PI/20.),
                  YMAX/4.+ scalef*GR[1][3]->gene_IDs->N*cos(M_PI+M_PI/20.)-10,bla1);
    
    // comparative
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[1][1]->gene_IDs,GR[1][2]->gene_IDs,GR[1][4]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Li_o=%ld   Lung_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Lu=%ld   Li_and_Diff=%ld   Lu_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef*2);
    EpsSetRgb(f,168/255.,1,3/255.);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/20.) ,
                    3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(+M_PI/20.) );
    EpsDrawEllipseArc(f, 3*XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      30+180/40.,150-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*abc->N*sin(-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-M_PI/20.) ,
                    3*XMAX/4.+ scalef*abc->N*sin(+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(+M_PI/20.) );
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(+M_PI/20.) -25,
                  YMAX/4.+ scalef*(ab->N+abc->N)*cos(+M_PI/20.) + 5,bla1);
    
    EpsSetLinewidth(f, scalef*(double)GR[1][0]->gene_IDs->N*2);
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/3.+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.));
    EpsDrawLine(f, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(M_PI/3.-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.));
    EpsSetRgb(f,168/255.,1,3/255.);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[1],ac->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.)+5,
                  YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.) + 5,bla1);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[0],bc->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.)+5,
                  YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.) + 5,bla1);
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f, 3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(-M_PI/3.+M_PI/20.));
    
    EpsDrawLine(f, 3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(M_PI/3.-M_PI/20.));
    
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc;delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[1][1]->gene_IDs,GR[1][3]->gene_IDs,GR[1][5]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Li_o=%ld   Ki_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Ki=%ld   Li_and_Diff=%ld   Ki_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef*2);
    EpsSetRgb(f,1,0,0);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.+M_PI/20.) );
    EpsDrawEllipseArc(f, 3*XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      150+180/40.,270-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*abc->N*sin(-2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-2*M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*abc->N*sin(-2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-2*M_PI/3.+M_PI/20.) );
    
    EpsSetLinewidth(f, scalef*(double)GR[1][0]->gene_IDs->N*2);
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.+M_PI/20.) - 50,
                  YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.+M_PI/20.)+5,bla1);
    
    
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/3.-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.));
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.));
    EpsSetRgb(f,1,0,0);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[2],ac->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.)-5,
                  YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.) - 5,bla1);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[0],bc->N);
    if(textok) EpsDrawString(f, 0,  3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.) - 15,
                  YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.) - 25,bla1);
    
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(-M_PI/3.-M_PI/20.));
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(-M_PI+M_PI/20.));
    
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc;delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[1][2]->gene_IDs,GR[1][3]->gene_IDs,GR[1][6]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Lu_o=%ld   Ki_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Lu_and_Ki=%ld   Lu_and_Diff=%ld   Ki_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef*2);
    EpsSetRgb(f,0,102/255.,1);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.+M_PI/20.) );
    EpsDrawEllipseArc(f, 3*XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      270+180/40.,30-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*abc->N*sin(2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(2*M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*abc->N*sin(2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(2*M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.+M_PI/20.)+10,
                  YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.+M_PI/20.)+5,bla1);
    
    EpsSetLinewidth(f, scalef*(double)GR[1][0]->gene_IDs->N*2);
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(M_PI/3.+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.));
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.));
    EpsSetRgb(f,0,102/255.,1);
    sprintf(bla1, "%ld NOT in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[2],ac->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.) - 5,
                  YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.) - 5,bla1);
    sprintf(bla1, "%ld NOT in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[1],bc->N);
    if(textok) EpsDrawString(f, 0,  3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.)+15,
                  YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.) + 5,bla1);
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(M_PI/3.+M_PI/20.));
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(-M_PI-M_PI/20.));
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc; delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[1][1]->gene_IDs,GR[1][2]->gene_IDs,GR[1][3]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Li_o=%ld   Ku_o=%ld   Ki_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Lu=%ld   Li_and_Ki=%ld   Lu_and_Ki=%ld   \n",ab->N,ac->N,bc->N);
    printf("All_3_Organs=%ld\n",abc->N);//getchar();
    EpsSetRgb(f, 0,0,0);
    EpsDrawCircle(f, 3*XMAX/4., YMAX/4., scalef*abc->N);
    //getchar();
    
    if(textok) sprintf(fn,"%s/Genotype_summary_text", path);
    else sprintf(fn,"%s/Genotype_summary_raw", path);
    close_EPSDRAW_file(f,fn,2);

    ao=new longint_ARRAY(); bo=new longint_ARRAY(); ab=new longint_ARRAY();
    VENN_DIAGRAM(GR[0][1]->gene_IDs,GR[1][1]->gene_IDs,ao,bo,ab);
    printf("Liver WT_only=%ld   Tie2_only=%ld   Both=%ld   \n",ao->N,bo->N,ab->N);getchar();
    delete ao;delete bo;delete ab;
        
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); ab=new longint_ARRAY();
    VENN_DIAGRAM(GR[0][2]->gene_IDs,GR[1][2]->gene_IDs,ao,bo,ab);
    printf("Lung WT_only=%ld   Tie2_only=%ld   Both=%ld   \n",ao->N,bo->N,ab->N);getchar();
    delete ao;delete bo;delete ab;
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); ab=new longint_ARRAY();
    VENN_DIAGRAM(GR[0][3]->gene_IDs,GR[1][3]->gene_IDs,ao,bo,ab);
    printf("Kidney WT_only=%ld   Tie2_only=%ld   Both=%ld   \n",ao->N,bo->N,ab->N);getchar();
    delete ao;delete bo;delete ab;
    
    return(GR);
}

Gene_Results_Matrix *process_Condition_cut(Exp_Design ED,const char path[],double FC, double padj,int exp_OK,int textok){
    Gene_Results_Matrix *GR;
    GR=new Gene_Results_Matrix[3];
    GR[0] =new Gene_Results_List[10];
    // LPS vs. control all for 0,
    // 1-3 : -- 3 tissues for index
    // 4-6 : LPS vs. control : 3 differential responses
    // 6-9: 3 tissue pair comparisons
    GR[1] =new Gene_Results_List[10];
    // same for Tie2
    int PlotGene=200;
    
    char fn[300];
    
    for(int i=0;i<ED.NC;i++){
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Tissues",path,ED.Conditions[i],ED.Conditions[i],
                ED.Genotypes[0],ED.Genotypes[1]); // e.g., WT_Control_to_LPS_all_Tissues.csv
        GR[i][0]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][0],path,fn,PlotGene,-1,0,-1,0);
        
        
        // single-tissue LPS vs. ct.
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Conditions[i],ED.Conditions[i],
                ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[0]);
        //e.g., WT_Control_to_LPS_in_Kidney.csv
        GR[i][1]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][1],path,fn,PlotGene,-1,0,0,0);
        
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Conditions[i],ED.Conditions[i],
                ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[1]);
        GR[i][2]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][2],path,fn,PlotGene,-1,0,1,0);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s",path,ED.Conditions[i],ED.Conditions[i],
                ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[2]);
        GR[i][3]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][3],path,fn,PlotGene,-1,0,2,0);
        
        // LPS vs. ct. differential in tissues
        
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s_over_%s",path,ED.Conditions[i],ED.Conditions[i],
                ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[0],ED.Tissues[1]);
        // WT_Control_to_LPS_in_Liver_over_Kidney.csv
        GR[i][4]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][4],path,fn,PlotGene,-1,0,2,1);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s_over_%s",path,ED.Conditions[i],ED.Conditions[i],
                ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[0],ED.Tissues[2]);
        GR[i][5]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][5],path,fn,PlotGene,-1,0,1,1);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_in_%s_over_%s",path,ED.Conditions[i],ED.Conditions[i],
                ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[1],ED.Tissues[2]);
        GR[i][6]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][6],path,fn,PlotGene,-1,0,0,1);
        
        // difference across tissues (all treatments)
       /*
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Genotypes",path,ED.Conditions[i],ED.Conditions[i],
                ED.Tissues[1],ED.Tissues[0]);
        // WT_Liver_to_Kidney_all_Conditions.csv
        GR[i][7]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][7],path,fn,PlotGene,-1,0,2,1);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Genotypes",path,ED.Conditions[i],ED.Conditions[i],
                ED.Tissues[0],ED.Tissues[2]);
        GR[i][8]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][8],path,fn,PlotGene,-1,0,1,1);
        
        sprintf(fn,"%s/%s/%s_%s_to_%s_all_Genotypes",path,ED.Conditions[i],ED.Conditions[i],
                ED.Tissues[1],ED.Tissues[2]);
        GR[i][9]=get_gene_list(fn,FC,padj);
        if(exp_OK) export_subset_for_heatmap(GR[i][9],path,fn,PlotGene,-1,0,0,1);
        */
    }
    
    if(textok) sprintf(fn,"%s/Condition_summary_text.eps", path);
    else sprintf(fn,"%s/Condition_summary_raw.eps", path);
    FILE *f;
    f=fopen(fn,"w");
    
    COLOR_SCHEME=GYR;
    setcolors_PAJEK();
    int XMAX,XMIN,YMAX, YMIN;
    
    XMAX=600;
    XMIN=0;
    YMAX=600;
    YMIN=0;
    double scalef = XMAX/400.;
    
    EpsInit(f,-10,-15,XMAX-XMIN+10,YMAX-YMIN+10);
    
    EpsSetFont(f,"Times-Roman",9);
    double dim=1.5;
    EpsSetLinewidth(f, scalef);
    
    //WT
    //  EpsSetRgb(f,0.9,0.9,0.9);
    //  EpsFillCircle(f,XMAX/4.,YMAX/4.,scalef*GR[0][0]->gene_IDs->N);
    //  EpsSetRgb(f,0.7,0.7,0.7);
    //  EpsFillCircle(f,XMAX/4.,YMAX/4.,scalef*get_NR_down(GR[0][0]));
    printf("WT: %u U %u  D \n",get_NR_up(GR[0][0]),get_NR_down(GR[0][0]));
    
    longint_ARRAY *ao,*bo,*co1,*ab,*ac,*bc,*abc;
    
    // Liver
    EpsSetRgb(f,0,219/255.,152/255.);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*GR[0][1]->gene_IDs->N*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*GR[0][1]->gene_IDs->N*cos(-M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*GR[0][1]->gene_IDs->N*sin(-M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*GR[0][1]->gene_IDs->N*cos(-M_PI/3.+M_PI/20.) );
    EpsSetRgb(f,0,(219/255.)/dim,2*(152/255.)/dim);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*get_NR_up(GR[0][1])*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[0][1])*cos(-M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*get_NR_up(GR[0][1])*sin(-M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[0][1])*cos(-M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[0],
            GR[0][1]->gene_IDs->N,get_NR_down(GR[0][1]),get_NR_up(GR[0][1]));
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*GR[0][1]->gene_IDs->N*sin(-M_PI/3.+M_PI/20.)+5,
                             YMAX/4.+ scalef*GR[0][1]->gene_IDs->N*cos(-M_PI/3.+M_PI/20.)+5,bla1);
    // Lung
    EpsSetRgb(f,1,147/255.,137/255.);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*GR[0][2]->gene_IDs->N*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*GR[0][2]->gene_IDs->N*cos(M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*GR[0][2]->gene_IDs->N*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*GR[0][2]->gene_IDs->N*cos(M_PI/3.+M_PI/20.) );
    EpsSetRgb(f,1./dim,(147/255.)/dim,(137/255.)/dim);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*get_NR_up(GR[0][2])*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[0][2])*cos(M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*get_NR_up(GR[0][2])*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[0][2])*cos(M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[1],
            GR[0][2]->gene_IDs->N,get_NR_down(GR[0][2]),get_NR_up(GR[0][2]));
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*GR[0][2]->gene_IDs->N*sin(M_PI/3.+M_PI/20.)+5,
                             YMAX/4.+ scalef*GR[0][2]->gene_IDs->N*cos(M_PI/3.+M_PI/20.)+5,bla1);
    
    // Kidney
    EpsSetRgb(f,1,138/255.,1);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*GR[0][3]->gene_IDs->N*sin(M_PI-M_PI/20.),
                    YMAX/4.+ scalef*GR[0][3]->gene_IDs->N*cos(M_PI-M_PI/20.) ,
                    XMAX/4.+ scalef*GR[0][3]->gene_IDs->N*sin(M_PI+M_PI/20.),
                    YMAX/4.+ scalef*GR[0][3]->gene_IDs->N*cos(M_PI+M_PI/20.) );
    EpsSetRgb(f,1/dim,(138/255.)/dim,1/dim);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*get_NR_up(GR[0][3])*sin(M_PI-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[0][3])*cos(M_PI-M_PI/20.) ,
                    XMAX/4.+ scalef*get_NR_up(GR[0][3])*sin(M_PI+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[0][3])*cos(M_PI+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[2],
            GR[0][3]->gene_IDs->N,get_NR_down(GR[0][3]),get_NR_up(GR[0][3]));
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*GR[0][3]->gene_IDs->N*sin(M_PI+M_PI/20.),
                             YMAX/4.+ scalef*GR[0][3]->gene_IDs->N*cos(M_PI+M_PI/20.)-10,bla1);
    
    // comparative
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[0][1]->gene_IDs,GR[0][2]->gene_IDs,GR[0][4]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Condition --  Control, Tie2 vs WT:   \n");
    printf("Li_o=%ld   Lung_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Lu=%ld   Li_and_Diff=%ld   Lu_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef);
    EpsSetRgb(f,168/255.,1,3/255.);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/20.) ,
                    XMAX/4.+ scalef*(ab->N+abc->N)*sin(+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(+M_PI/20.) );
    EpsDrawEllipseArc(f, XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      30+180/40.,150-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*abc->N*sin(-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-M_PI/20.) ,
                    XMAX/4.+ scalef*abc->N*sin(+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(+M_PI/20.) );
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N)*sin(+M_PI/20.) -25,
                             YMAX/4.+ scalef*(ab->N+abc->N)*cos(+M_PI/20.) + 5,bla1);
    
    EpsSetLinewidth(f, scalef*(double)GR[0][0]->gene_IDs->N*5);
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f, XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/3.+M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.));
    EpsDrawLine(f, XMAX/4.+ scalef*(ab->N+abc->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(M_PI/3.-M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.));
    EpsSetRgb(f,168/255.,1,3/255.);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[1],ac->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.)+5,
                             YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.) + 5,bla1);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[0],bc->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.)+5,
                             YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.) + 5,bla1);
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f, XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(-M_PI/3.+M_PI/20.));
    
    EpsDrawLine(f, XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(M_PI/3.-M_PI/20.));
    
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc;delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[0][1]->gene_IDs,GR[0][3]->gene_IDs,GR[0][5]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Condition --  Control, Tie2 vs WT:   \n");
    printf("Li_o=%ld   Ki_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Ki=%ld   Li_and_Diff=%ld   Ki_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef);
    EpsSetRgb(f,1,0,0);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.+M_PI/20.) );
    EpsDrawEllipseArc(f, XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      150+180/40.,270-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*abc->N*sin(-2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-2*M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*abc->N*sin(-2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-2*M_PI/3.+M_PI/20.) );
    EpsSetLinewidth(f, scalef*(double)GR[0][0]->gene_IDs->N*5);
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.+M_PI/20.) - 50,
                             YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.+M_PI/20.)+5,bla1);
    
    
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/3.-M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.));
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI+M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.));
    EpsSetRgb(f,1,0,0);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[2],ac->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.)-5,
                             YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.) - 5,bla1);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[0],bc->N);
    if(textok) EpsDrawString(f, 0,  XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.) - 15,
                             YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.) - 25,bla1);
    
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(-M_PI/3.-M_PI/20.));
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(-M_PI+M_PI/20.));
    
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc;delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[0][2]->gene_IDs,GR[0][3]->gene_IDs,GR[0][6]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Condition --  Control, Tie2 vs WT:   \n");
    printf("Lu_o=%ld   Ki_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Lu_and_Ki=%ld   Lu_and_Diff=%ld   Ki_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef);
    EpsSetRgb(f,0,102/255.,1);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.+M_PI/20.) );
    EpsDrawEllipseArc(f, XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      270+180/40.,30-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, XMAX/4.,YMAX/4., XMAX/4.+ scalef*abc->N*sin(2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(2*M_PI/3.-M_PI/20.) ,
                    XMAX/4.+ scalef*abc->N*sin(2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(2*M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.+M_PI/20.)+10,
                             YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.+M_PI/20.)+5,bla1);
    
    EpsSetLinewidth(f, scalef*(double)GR[0][0]->gene_IDs->N*5);
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(M_PI/3.+M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.));
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI-M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.));
    EpsSetRgb(f,0,102/255.,1);
    sprintf(bla1, "%ld NOT in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[2],ac->N);
    if(textok) EpsDrawString(f, 0, XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.) - 5,
                             YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.) - 5,bla1);
    sprintf(bla1, "%ld NOT in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[1],bc->N);
    if(textok) EpsDrawString(f, 0,  XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.)+15,
                             YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.) + 5,bla1);
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(M_PI/3.+M_PI/20.));
    EpsDrawLine(f,  XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.),
                XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(-M_PI-M_PI/20.));
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc; delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[0][1]->gene_IDs,GR[0][2]->gene_IDs,GR[0][3]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Condition --  Control, Tie2 vs WT:   \n");
    printf("Li_o=%ld   Lu_o=%ld   Ki_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Lu=%ld   Li_and_Ki=%ld   Lu_and_Ki=%ld   \n",ab->N,ac->N,bc->N);
    printf("All_3_Organs=%ld\n",abc->N);//getchar();
    EpsSetRgb(f, 0,0,0);
    EpsDrawCircle(f, XMAX/4., YMAX/4., scalef*abc->N);
    //getchar();
    
    
    
    
    EpsSetLinewidth(f, scalef*2);
    
    // LPS
    //  EpsSetRgb(f,0.9,0.9,0.9);
    //  EpsFillCircle(f,XMAX/4.,YMAX/4.,scalef*GR[0][0]->gene_IDs->N);
    //  EpsSetRgb(f,0.7,0.7,0.7);
    //  EpsFillCircle(f,XMAX/4.,YMAX/4.,scalef*get_NR_down(GR[0][0]));
    printf("LPS: %u U %u  D \n",get_NR_up(GR[1][0]),get_NR_down(GR[1][0]));
    EpsSetFont(f,"Times-Roman",5);
    
    // Liver
    EpsSetRgb(f,0,219/255.,152/255.);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*GR[1][1]->gene_IDs->N*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*GR[1][1]->gene_IDs->N*cos(-M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*GR[1][1]->gene_IDs->N*sin(-M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*GR[1][1]->gene_IDs->N*cos(-M_PI/3.+M_PI/20.) );
    EpsSetRgb(f,0,(219/255.)/dim,2*(152/255.)/dim);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*get_NR_up(GR[1][1])*sin(-M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[1][1])*cos(-M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*get_NR_up(GR[1][1])*sin(-M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[1][1])*cos(-M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[0],
            GR[1][1]->gene_IDs->N,get_NR_down(GR[1][1]),get_NR_up(GR[1][1]));
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*GR[1][1]->gene_IDs->N*sin(-M_PI/3.+M_PI/20.)+5,
                             YMAX/4.+ scalef*GR[1][1]->gene_IDs->N*cos(-M_PI/3.+M_PI/20.)+5,bla1);
    // Lung
    EpsSetRgb(f,1,147/255.,137/255.);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*GR[1][2]->gene_IDs->N*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*GR[1][2]->gene_IDs->N*cos(M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*GR[1][2]->gene_IDs->N*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*GR[1][2]->gene_IDs->N*cos(M_PI/3.+M_PI/20.) );
    EpsSetRgb(f,1./dim,(147/255.)/dim,(137/255.)/dim);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*get_NR_up(GR[1][2])*sin(M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[1][2])*cos(M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*get_NR_up(GR[1][2])*sin(M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[1][2])*cos(M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[1],
            GR[1][2]->gene_IDs->N,get_NR_down(GR[1][2]),get_NR_up(GR[1][2]));
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*GR[1][2]->gene_IDs->N*sin(M_PI/3.+M_PI/20.)+5,
                             YMAX/4.+ scalef*GR[1][2]->gene_IDs->N*cos(M_PI/3.+M_PI/20.)+5,bla1);
    
    // Kidney
    EpsSetRgb(f,1,138/255.,1);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*GR[1][3]->gene_IDs->N*sin(M_PI-M_PI/20.),
                    YMAX/4.+ scalef*GR[1][3]->gene_IDs->N*cos(M_PI-M_PI/20.) ,
                    3*XMAX/4.+ scalef*GR[1][3]->gene_IDs->N*sin(M_PI+M_PI/20.),
                    YMAX/4.+ scalef*GR[1][3]->gene_IDs->N*cos(M_PI+M_PI/20.) );
    EpsSetRgb(f,1/dim,(138/255.)/dim,1/dim);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*get_NR_up(GR[1][3])*sin(M_PI-M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[1][3])*cos(M_PI-M_PI/20.) ,
                    3*XMAX/4.+ scalef*get_NR_up(GR[1][3])*sin(M_PI+M_PI/20.),
                    YMAX/4.+ scalef*get_NR_up(GR[1][3])*cos(M_PI+M_PI/20.) );
    sprintf(bla1, "%s: %ld genes (%u U; %u  D )",ED.Tissues[2],
            GR[1][3]->gene_IDs->N,get_NR_down(GR[1][3]),get_NR_up(GR[1][3]));
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*GR[1][3]->gene_IDs->N*sin(M_PI+M_PI/20.),
                             YMAX/4.+ scalef*GR[1][3]->gene_IDs->N*cos(M_PI+M_PI/20.)-10,bla1);
    
    // comparative
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[1][1]->gene_IDs,GR[1][2]->gene_IDs,GR[1][4]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Condition --  LPS, Tie2 vs WT:   \n");
    printf("Li_o=%ld   Lung_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Lu=%ld   Li_and_Diff=%ld   Lu_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef*2);
    EpsSetRgb(f,168/255.,1,3/255.);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/20.) ,
                    3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(+M_PI/20.) );
    EpsDrawEllipseArc(f, 3*XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      30+180/40.,150-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*abc->N*sin(-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-M_PI/20.) ,
                    3*XMAX/4.+ scalef*abc->N*sin(+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(+M_PI/20.) );
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(+M_PI/20.) -25,
                             YMAX/4.+ scalef*(ab->N+abc->N)*cos(+M_PI/20.) + 5,bla1);
    
    EpsSetLinewidth(f, scalef*(double)GR[1][0]->gene_IDs->N*2);
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/3.+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.));
    EpsDrawLine(f, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(M_PI/3.-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.));
    EpsSetRgb(f,168/255.,1,3/255.);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[1],ac->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.)+5,
                             YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.) + 5,bla1);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[0],bc->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.)+5,
                             YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.) + 5,bla1);
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f, 3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(-M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(-M_PI/3.+M_PI/20.));
    
    EpsDrawLine(f, 3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(M_PI/3.-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(M_PI/3.-M_PI/20.));
    
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc;delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[1][1]->gene_IDs,GR[1][3]->gene_IDs,GR[1][5]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Condition --  LPS, Tie2 vs WT:   \n");
    printf("Li_o=%ld   Ki_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Ki=%ld   Li_and_Diff=%ld   Ki_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef*2);
    EpsSetRgb(f,1,0,0);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.+M_PI/20.) );
    EpsDrawEllipseArc(f, 3*XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      150+180/40.,270-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*abc->N*sin(-2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-2*M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*abc->N*sin(-2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(-2*M_PI/3.+M_PI/20.) );
    
    EpsSetLinewidth(f, scalef*(double)GR[1][0]->gene_IDs->N*2);
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-2*M_PI/3.+M_PI/20.) - 50,
                             YMAX/4.+ scalef*(ab->N+abc->N)*cos(-2*M_PI/3.+M_PI/20.)+5,bla1);
    
    
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI/3.-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.));
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.));
    EpsSetRgb(f,1,0,0);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[2],ac->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.)-5,
                             YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.) - 5,bla1);
    sprintf(bla1, "%ld not in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[0],bc->N);
    if(textok) EpsDrawString(f, 0,  3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.) - 15,
                             YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.) - 25,bla1);
    
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(-M_PI/3.-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(-M_PI/3.-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(-M_PI/3.-M_PI/20.));
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(-M_PI+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(-M_PI+M_PI/20.));
    
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc;delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[1][2]->gene_IDs,GR[1][3]->gene_IDs,GR[1][6]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Condition --  LPS, Tie2 vs WT:   \n");
    printf("Lu_o=%ld   Ki_o=%ld   Diff_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Lu_and_Ki=%ld   Lu_and_Diff=%ld   Ki_and_Diff=%ld   \n",ab->N,ac->N,bc->N);
    printf("All=%ld\n",abc->N);
    
    EpsSetLinewidth(f, scalef*2);
    EpsSetRgb(f,0,102/255.,1);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.+M_PI/20.) );
    EpsDrawEllipseArc(f, 3*XMAX/4.,YMAX/4.,scalef*(ab->N+abc->N),scalef*(ab->N+abc->N),
                      270+180/40.,30-180/40.);
    EpsSetRgb(f,0,0,0);
    EpsFillTriangle(f, 3*XMAX/4.,YMAX/4., 3*XMAX/4.+ scalef*abc->N*sin(2*M_PI/3.-M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(2*M_PI/3.-M_PI/20.) ,
                    3*XMAX/4.+ scalef*abc->N*sin(2*M_PI/3.+M_PI/20.),
                    YMAX/4.+ scalef*abc->N*cos(2*M_PI/3.+M_PI/20.) );
    sprintf(bla1, "%ld shared (%ld opp.)",ab->N+abc->N,abc->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(2*M_PI/3.+M_PI/20.)+10,
                             YMAX/4.+ scalef*(ab->N+abc->N)*cos(2*M_PI/3.+M_PI/20.)+5,bla1);
    
    EpsSetLinewidth(f, scalef*(double)GR[1][0]->gene_IDs->N*2);
    EpsSetRgb(f,0.9,0.9,0.9);
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(M_PI/3.+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.));
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N)*cos(-M_PI-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.));
    EpsSetRgb(f,0,102/255.,1);
    sprintf(bla1, "%ld NOT in %s (%ld stat. sig. diff.)",ao->N+ac->N,ED.Tissues[2],ac->N);
    if(textok) EpsDrawString(f, 0, 3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.) - 5,
                             YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.) - 5,bla1);
    sprintf(bla1, "%ld NOT in %s (%ld stat. sig. diff.)",bo->N+bc->N,ED.Tissues[1],bc->N);
    if(textok) EpsDrawString(f, 0,  3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.)+15,
                             YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.) + 5,bla1);
    
    EpsSetRgb(f,0.5,0.5,0.5);
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N)*cos(M_PI/3.+M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*sin(M_PI/3.+M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+ao->N+ac->N)*cos(M_PI/3.+M_PI/20.));
    EpsDrawLine(f,  3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N)*cos(-M_PI-M_PI/20.),
                3*XMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*sin(-M_PI-M_PI/20.),
                YMAX/4.+ scalef*(ab->N+abc->N+bo->N+bc->N)*cos(-M_PI-M_PI/20.));
    
    delete ao;delete bo;delete co1;delete ab;delete ac;delete bc; delete abc;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); co1=new longint_ARRAY();
    ab=new longint_ARRAY(); ac=new longint_ARRAY(); bc=new longint_ARRAY();
    abc=new longint_ARRAY();
    VENN_DIAGRAM_in_three(GR[1][1]->gene_IDs,GR[1][2]->gene_IDs,GR[1][3]->gene_IDs,ao,bo,co1,ab,ac,bc,abc);
    printf("Condition --  LPS, Tie2 vs WT:   \n");
    printf("Li_o=%ld   Ku_o=%ld   Ki_o=%ld   \n",ao->N,bo->N,co1->N);
    printf("Li_and_Lu=%ld   Li_and_Ki=%ld   Lu_and_Ki=%ld   \n",ab->N,ac->N,bc->N);
    printf("All_3_Organs=%ld\n",abc->N);//getchar();
    EpsSetRgb(f, 0,0,0);
    EpsDrawCircle(f, 3*XMAX/4., YMAX/4., scalef*abc->N);
    //getchar();
    
    if(textok) sprintf(fn,"%s/Condition_summary_text", path);
    else sprintf(fn,"%s/Condition_summary_raw", path);
    close_EPSDRAW_file(f,fn,2);
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); ab=new longint_ARRAY();
    VENN_DIAGRAM(GR[0][1]->gene_IDs,GR[1][1]->gene_IDs,ao,bo,ab);
    printf("Tie2 vs WT in:   \n");
    printf("Liver Control_only=%ld   LPS_only=%ld   Both=%ld   \n",ao->N,bo->N,ab->N);getchar();
    delete ao;delete bo;delete ab;
    
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); ab=new longint_ARRAY();
    VENN_DIAGRAM(GR[0][2]->gene_IDs,GR[1][2]->gene_IDs,ao,bo,ab);
    printf("Lung Control_only=%ld   LPS_only=%ld   Both=%ld   \n",ao->N,bo->N,ab->N);getchar();
    delete ao;delete bo;delete ab;
    ao=new longint_ARRAY(); bo=new longint_ARRAY(); ab=new longint_ARRAY();
    VENN_DIAGRAM(GR[0][3]->gene_IDs,GR[1][3]->gene_IDs,ao,bo,ab);
    printf("Kidney Control_only=%ld   LPS_only=%ld   Both=%ld   \n",ao->N,bo->N,ab->N);getchar();
    delete ao;delete bo;delete ab;
    
    return(GR);
}

longint_ARRAY *get_collection_of_all_relevant_genes(Gene_Results_Matrix *GR_Gen,Gene_Results_Matrix *GR_Cond,Exp_Design ED){
    longint_ARRAY *gl;
    
    gl=new longint_ARRAY();
    for(int i=0;i<ED.NC;i++)
        for(int j=1;j<=6;j++){
            gl->add_longint_ARRAY(GR_Cond[i][j]->gene_IDs);
        }
    for(int i=0;i<ED.NG;i++)
        for(int j=1;j<=6;j++){
            gl->add_longint_ARRAY(GR_Gen[i][j]->gene_IDs);
        }
    return (gl);
}

void Ingrid_read_gene_lists(const char *path,double FC, double padj){
    Exp_Design ED;
    ED.NT=3;
    ED.NC=2;
    ED.NG=2;
    ED.Tissues = new charpointer[ED.NT];
    ED.Conditions = new charpointer[ED.NC];
    ED.Genotypes = new charpointer[ED.NG];
    ED.Tissues[0] = new char[100]; sprintf(ED.Tissues[0],"Liver");
    ED.Tissues[1] = new char[100]; sprintf(ED.Tissues[1],"Lung");
    ED.Tissues[2] = new char[100]; sprintf(ED.Tissues[2],"Kidney");
    ED.Conditions[0] = new char[100]; sprintf(ED.Conditions[0],"Control");
    ED.Conditions[1] = new char[100]; sprintf(ED.Conditions[1],"LPS");
    ED.Genotypes[0] = new char[100]; sprintf(ED.Genotypes[0],"WT");
    ED.Genotypes[1] = new char[100]; sprintf(ED.Genotypes[1],"Tie2");
    
    FILE *f;
    char fn[300];
    char attr_list[100][5000];
    unsigned long int Nf;
    
    sprintf(fn,"%s/Gene_list.txt",path);
    f=fopen(fn,"r");
    printf("%s\n",fn);//getchar();
    fscanf(f,"%ld",&N_GENES);
    printf("%s\nN=%ld\n",fn,N_GENES);//getchar();
    fgets(GSE_line,MAX_LINE,f);
    fgets(GSE_line,MAX_LINE,f);

    G_names  =  new GeneList [N_GENES+1];
    for (unsigned long int i=1; i<=N_GENES; i++) {
        fgets(GSE_line,MAX_LINE,f);
        Nf=separate_line(GSE_line,attr_list,'\t',5000);
        G_names[i]= new Gene;
        G_names[i]->ID=atol(attr_list[1]);
        strcpy(G_names[i]->ENSMUSG,attr_list[2]);
        strcpy(G_names[i]->Symbol,attr_list[3]);
        strcpy(G_names[i]->GeneName,attr_list[4]);
        G_names[i]->ENTREZ=atol(attr_list[5]);
        if(i/5000==i/5000.){
            printf("i=%ld\t%s\t%s\t%s\t%ld\n",i,G_names[i]->ENSMUSG,G_names[i]->Symbol,G_names[i]->GeneName,G_names[i]->ENTREZ);//getchar();
        }
        //fscanf(f,"%ld\t%s\t%s",&(G_names[i]->ID),G_names[i]->ENSMUSG,G_names[i]->Symbol);
    }
    fclose(f);
    for (unsigned int i=1; i<=10; i++) {
        printf("%lu\t%s\t%s\n",G_names[i]->ID,G_names[i]->ENSMUSG,G_names[i]->Symbol);
    }
   // getchar();
    
    Gene_Results_Matrix *GR_Gen,*GR_Cond;
    
    GR_Gen=process_Genotype_cut(ED,path,FC,padj,0,0);
    GR_Cond=process_Condition_cut(ED,path,FC,padj,0,0);
  
    longint_ARRAY *gl_all;
    gl_all=get_collection_of_all_relevant_genes(GR_Gen,GR_Cond,ED);
   
    for(int i=0;i<ED.NG;i++)
        for(int j=1;j<=6;j++){
            gl_all->add_longint_ARRAY(GR_Gen[i][j]->gene_IDs);
        }
    
    for(int i=0;i<ED.NC;i++)
        for(int j=1;j<=6;j++){
            gl_all->add_longint_ARRAY(GR_Gen[i][j]->gene_IDs);
        }
 
    sprintf(fn,"%s/__ALL_Diff_Expressed_Genes_Compendium.txt",path);
    f=fopen(fn,"w");
    fprintf(f, "NR\tENSMUSG\tSYMBOL\tGENE_NAME\tENTREZ_ID");
    for(int i=0;i<ED.NG;i++){
        fprintf(f,"\t%s__%s_to_%s_in_%s",ED.Genotypes[i],ED.Conditions[0],ED.Conditions[1],ED.Tissues[0]);
        fprintf(f,"\t%s__%s_to_%s_in_%s",ED.Genotypes[i],ED.Conditions[0],ED.Conditions[1],ED.Tissues[1]);
        fprintf(f,"\t%s__%s_to_%s_in_%s",ED.Genotypes[i],ED.Conditions[0],ED.Conditions[1],ED.Tissues[2]);
        
        fprintf(f,"\t%s__%s_to_%s_in_%s_over_%s",
                ED.Genotypes[i],ED.Conditions[0],ED.Conditions[1],ED.Tissues[0],ED.Tissues[1]);
        fprintf(f,"\t%s__%s_to_%s_in_%s_over_%s",
                ED.Genotypes[i],ED.Conditions[0],ED.Conditions[1],ED.Tissues[0],ED.Tissues[2]);
        fprintf(f,"\t%s__%s_to_%s_in_%s_over_%s",
                ED.Genotypes[i],ED.Conditions[0],ED.Conditions[1],ED.Tissues[1],ED.Tissues[2]);
    }
    for(int i=0;i<ED.NC;i++){
        fprintf(f,"\t%s__%s_to_%s_in_%s",ED.Conditions[i],ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[0]);
        fprintf(f,"\t%s__%s_to_%s_in_%s",ED.Conditions[i],ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[1]);
        fprintf(f,"\t%s__%s_to_%s_in_%s",ED.Conditions[i],ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[2]);
        fprintf(f,"\t%s__%s_to_%s_in_%s_over_%s",
                ED.Conditions[i],ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[0],ED.Tissues[1]);
        fprintf(f,"\t%s__%s_to_%s_in_%s_over_%s",
                ED.Conditions[i],ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[0],ED.Tissues[2]);
        fprintf(f,"\t%s__%s_to_%s_in_%s_over_%s",
                ED.Conditions[i],ED.Genotypes[0],ED.Genotypes[1],ED.Tissues[1],ED.Tissues[2]);
    }
    fprintf(f, "\n");
    
    unsigned long int id1;
    
    for(unsigned long int k=1;k<=gl_all->N;k++){
        fprintf(f, "%ld\t%s\t%s\t%s\t%ld",k,G_names[gl_all->A[k]]->ENSMUSG,G_names[gl_all->A[k]]->Symbol,G_names[gl_all->A[k]]->GeneName,G_names[gl_all->A[k]]->ENTREZ);
        for(int i=0;i<ED.NG;i++)
            for(int j=1;j<=6;j++){
                id1=GR_Gen[i][j]->gene_IDs->check_for_element(gl_all->A[k]);
                fprintf(f, "\t%lg",GR_Gen[i][j]->FC[id1]);
            }
        for(int i=0;i<ED.NC;i++)
            for(int j=1;j<=6;j++){
                id1=GR_Cond[i][j]->gene_IDs->check_for_element(gl_all->A[k]);
                fprintf(f, "\t%lg",GR_Cond[i][j]->FC[id1]);
            }
        fprintf(f, "\n");
    }
    
    fclose(f);
}

#endif /* Ingrid_NGS_Li_Lu_Ki_LPS_Tie2_h */
