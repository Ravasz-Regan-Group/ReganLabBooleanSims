#include<math.h>

void write_Undir_Dot_to_dir(sgraph *g,char pn[]){
	 FILE *net;
	 long long i;
	 char col[10];
	 char dotname[100];
	 snode::slinklist::slink *sc2;

	 if(g->N>500) {printf("Graph too large for GraphViz picture\n");return;}
	 sprintf(dotname,"%s/%s.dot",pn,g->gname);
	 net=fopen(dotname,"w");fprintf(net,"graph G {\n	graph [overlap=false];\n"); 
	 for(i=1;i<=g->N;i++){
		if(g->node[i].node_props==NULL)
			fprintf(net,"%lld [fillcolor=\"Red\",style=filled,shape=\"circle\"];\n",i);
		else 
			if(g->node[i].node_draw_props==NULL) {
				set_RGB_hex(g->node[i].node_props->nw*255,0,0,col);
				fprintf(net,"%lld [fillcolor=\"%s\",style=filled];\n",
					i,col);
			}
			else {fprintf(net,"%lld [label=\"%s\"];\n",i,g->node[i].node_draw_props->node_label);	
			  //printf("node %d has props: %s\n",sc1->ID,sc1->node_props->node_rgb.RGB_hex); getchar();
			 }
	 }
	 for(i=1;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;//printf("%d ",sc2->sneighbor->ID);getchar();
		//printf("\nlinks of %d:\t",sc1->ID);
		//sc1->links_of_node.list_node_linklist();	
		while(sc2!=NULL){
			//printf("1: %d ",sc2->sneighbor->ID);getchar();
			if(sc2->link_props==NULL) fprintf(net,"\t%lld -- %lld ;\n",i,sc2->sneighbor);
//				else fprintf(net,"\t%d -- %d [color = \"%s\",label=\"%lg\"];\n",i,sc2->sneighbor,sc2->link_props->color,
//								sc2->link_props->lw);
			else fprintf(net,"\t%lld -- %lld [label=\"%lg\"];\n",i,sc2->sneighbor, sc2->link_props->lw);
			//printf("2:  %d ",sc2->sneighbor->ID);getchar();
			sc2=sc2->next(sc2);//getchar();
		}
	}	
	fprintf(net,"}"); fflush(net);
	fclose(net);
//	printf("End.\n");//getchar();
}

void write_simple_Pajek_file(sgraph *g,const char pn[]){
	FILE *f;
	long long i;
	char dotname[100];
	snode::slinklist::slink *sc2;
	double m;
	
	sprintf(dotname,"%s/%s.net",pn,g->gname);
	f=fopen(dotname,"w");
	
	fprintf(f,"*Vertices %lld\r\n",g->N);
	for(i=1;i<=g->N;i++){
		//fprintf(f,"  %d  %d  %lf  %lf  ic Gray20 bc Gray20\n",i,i,x[i],y[i]);
		if(g->node[i].node_draw_props!=NULL)
			fprintf(f,"  %lld  \"%s\" \r\n",i,g->node[i].node_draw_props->node_label);
		else 	fprintf(f,"  %lld  %lld \r\n",i,g->node[i].ID);
	}
	
	fprintf(f,"*Edges\r\n");
	
	for(i=1;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;//printf("%d ",sc2->sneighbor->ID);getchar();
		//printf("\nlinks of %d:\t",sc1->ID);
		//sc1->links_of_node.list_node_linklist();	
		while(sc2!=NULL){
			//printf("1: %d ",sc2->sneighbor->ID);getchar();
			if(sc2->link_props==NULL) 
				 fprintf(f,"  %lld  %lld c Black\r\n",i,sc2->sneighbor);
			else {
				if(g->node[i].node_draw_props->basic->nw<g->node[sc2->sneighbor].node_draw_props->basic->nw)
					m=g->node[i].node_draw_props->basic->nw;
				else m=g->node[sc2->sneighbor].node_draw_props->basic->nw;
				if((sc2->link_props->lw/m)>0.25)
					fprintf(f,"  %lld  %lld  %lg c Black\r\n",i,sc2->sneighbor, sc2->link_props->lw/m);
			}
			//printf("2:  %d ",sc2->sneighbor->ID);getchar();
			sc2=sc2->next(sc2);//getchar();
		}
	}	
	fclose(f);
}		


void write_Undir_Dot(sgraph *g,char dotname[100]){
	 FILE *net;
	 long long i;
	 char col[10];
	 snode::slinklist::slink *sc2;

	 if(g->N>500) {printf("Graph too large for GraphViz picture\n");return;}
	 net=fopen(dotname,"w");fprintf(net,"graph G {\n	graph [overlap=false];\n"); 
	 for(i=1;i<=g->N;i++){
		if(g->node[i].node_props==NULL)
			fprintf(net,"%lld [fillcolor=\"Red\",style=filled,shape=\"circle\"];\n",i);
		else 
			if(g->node[i].node_draw_props==NULL){ 
				set_RGB_hex(g->node[i].node_props->nw*255,0,0,col);
				fprintf(net,"%lld [fillcolor=\"%s\",style=filled,label=\"%llu: %.3lf\"];\n",
					i,col,g->node[i].node_props->ID,g->node[i].node_props->nw);
				}
			else {fprintf(net,"%lld [fillcolor=\"%s\",style=filled,shape=\"%s\",label=\"%s\"];\n",
					 i,g->node[i].node_draw_props->node_rgb.RGB_hex,
					 g->node[i].node_draw_props->node_shape,g->node[i].node_draw_props->node_label);	
			  //printf("node %d has props: %s\n",sc1->ID,sc1->node_props->node_rgb.RGB_hex); getchar();
			 }
	 }
	 for(i=1;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;//printf("%d ",sc2->sneighbor->ID);getchar();
		//printf("\nlinks of %d:\t",sc1->ID);
		//sc1->links_of_node.list_node_linklist();	
		while(sc2!=NULL){
			//printf("1: %d ",sc2->sneighbor->ID);getchar();
			if(sc2->link_props==NULL) fprintf(net,"\t%lld -- %llu ;\n",i,sc2->sneighbor);
				else fprintf(net,"\t%lld -- %llu [color = \"%s\"];\n",i,sc2->sneighbor,sc2->link_props->color);
			//printf("2:  %d ",sc2->sneighbor->ID);getchar();
			sc2=sc2->next(sc2);//getchar();
		}
	}	
	fprintf(net,"}"); fflush(net);
	fclose(net);
//	printf("End.\n");//getchar();
}

void write_Dir_Dot(sgraph *g,char dotname[100]){
	 FILE *net;
	 int i;
	 snode::slinklist::slink *sc2;
	 if(g->N>500) {printf("Graph too large for GraphViz picture\n");return;}
	 
	 net=fopen(dotname,"w");fprintf(net,"digraph G {\n"); 
	 for(i=1;i<=g->N;i++){
		if(g->node[i].node_draw_props==NULL)
			{fprintf(net,"%d [fillcolor=\"Red\",style=filled,shape=\"circle\"];\n",i);
			 //printf("node %d has no props!\n",sc1->ID); getchar();
			}
		else {fprintf(net,"%d [fillcolor=\"%s\",style=filled,label=\"%s\"];\n",
					 i,g->node[i].node_draw_props->node_rgb.RGB_hex,
					 g->node[i].node_draw_props->node_label);	
			  //printf("node %d has props: %s\n",sc1->ID,sc1->node_props->node_rgb.RGB_hex); getchar();
			 }
	 }
	 for(i=1;i<=g->N;i++){
		//if(sc1==NULL) {printf("How on earth?!?\n");getchar();}
		sc2=g->node[i].Li.first_link_pt;//printf("%d ",sc2->sneighbor->ID);getchar();
		//printf("\nlinks of %d:\t",sc1->ID);
		//sc1->links_of_node.list_node_linklist();	
		while(sc2!=NULL){
			//printf("1: %d ",sc2->sneighbor->ID);getchar();
			if(sc2->link_props==NULL) fprintf(net,"\t%d -> %llu ;\n",i,sc2->sneighbor);
				else fprintf(net,"\t%d -> %llu [color = \"%s\"];\n",
					i,sc2->sneighbor,sc2->link_props->color);
			//printf("2:  %d ",sc2->sneighbor->ID);getchar();
			sc2=sc2->next(sc2);//getchar();
		}
		//if(sc1!=NULL) printf("\tgoing to %d!",sc1->ID);
			//else {printf("Should stop here!\n");sc1=NULL;getchar();}
	 }
	fprintf(net,"}"); fflush(net);
	fclose(net);
}

void write_linklist(sgraph *g,char dotname[100]){
	 FILE *net;
	 unsigned long long i;
	 snode::slinklist::slink *sc2;
	 
	 net=fopen(dotname,"w");
	 for(i=1;i<=g->N;i++){
	 	//if(sc1==NULL) {printf("How on earth?!?\n");getchar();}
		sc2=g->node[i].Li.first_link_pt;//printf("%d ",sc2->sneighbor->ID);getchar();
		//printf("\nlinks of %d:\t",sc1->ID);
		//sc1->links_of_node.list_node_linklist();	
		while(sc2!=NULL){
			//printf("1: %d ",sc2->sneighbor->ID);getchar();
			//if(sc2->link_props==NULL) 
				fprintf(net,"%lld\t%lld\n",i,sc2->sneighbor);
			//	else fprintf(net,"%lld\t%lld\t%lf\n",i,sc2->sneighbor,sc2->link_props->lw);
			//printf("2:  %d ",sc2->sneighbor->ID);getchar();
			sc2=sc2->next(sc2);//getchar();
		}
	 }
	fclose(net);
}

void write_linklist_weights(sgraph *g,const char path[]){
	 FILE *net;
	 unsigned long long i;
	 snode::slinklist::slink *sc2;
	 char dotname[400];
	 
	 sprintf(dotname,"%s/%s.net",path,g->gname);
	 net=fopen(dotname,"w");
	 for(i=1;i<=g->N;i++){
	 //	if(g->node[i].node_props==NULL){
		//	printf("Node with no props! Wrong routine to write this network out.");
		//	getchar();return;
		//}
		sc2=g->node[i].Li.first_link_pt;	
		while(sc2!=NULL){
			if(sc2->link_props==NULL) 
				fprintf(net,"%lld %lf\t%lld %lf\n",
						i,g->node[i].node_props->nw,
						sc2->sneighbor,g->node[sc2->sneighbor].node_props->nw);
			//	else fprintf(net,"%lld %lf\t%lld %lf\t%lf\n",
			//			i,g->node[i].node_props->nw,
			//			sc2->sneighbor,g->node[sc2->sneighbor].node_props->nw,
			//			sc2->link_props->lw);
            else fprintf(net,"%lld\t%lld\t%lf\n",
                         i,sc2->sneighbor,sc2->link_props->lw);
            sc2=sc2->next(sc2);
		}
	 }
	fclose(net);
}

void write_linklist_weight_ONLY(sgraph *g,char path[]){
	FILE *net;
	unsigned long long i;
	snode::slinklist::slink *sc2;
	char dotname[400];
	
	sprintf(dotname,"%s/%s.net",path,g->gname);
	net=fopen(dotname,"w");
	for(i=1;i<=g->N;i++){
		//	if(g->node[i].node_props==NULL){
		//	printf("Node with no props! Wrong routine to write this network out.");
		//	getchar();return;
		//}
		sc2=g->node[i].Li.first_link_pt;	
		while(sc2!=NULL){
			if(sc2->link_props==NULL) 
				fprintf(net,"%lld\t%lld\n",
						i,sc2->sneighbor);
			else fprintf(net,"%lld %llu\t%g\n",
						 i,sc2->sneighbor,sc2->link_props->lw);
			sc2=sc2->next(sc2);
		}
	}
	fclose(net);
}

struct sgraph *read_linklist_weights(char di[],char ne[]){
	unsigned long long i = 0,j,t,id1,id2;
	unsigned long long NR=0;
	char nev[300];
	sgraph *gall;
	FILE *be;
	double lw,b1,b2;
	
		
		sprintf(nev,"%s/%s.net",di,ne);
		be=fopen(nev,"r");
		j=1;
		while(j!=EOF){
			j=fscanf(be,"%lld%lf %lld%lf %lf",&id1,&b1,&id2,&b2,&lw);
			printf("%lld %lf %lld %lf %lf\n",id1,b1,id2,b2,lw);
			if(id1>NR) NR=id1;
			if(id2>NR) NR=id2;
		}
		fclose(be);
		
		gall=new class sgraph(NR,ne);	
		printf("NR=%lld\n",gall->N);//getchar();
		printf("Network in: %s\n",nev);
		be=fopen(nev,"r");
		t=1;
		while(t!=EOF){
			t=fscanf(be,"%lld%lf %lld%lf %lf",&id1,&b1,&id2,&b2,&lw);
			//printf("%d\n",id1);getchar();
			if(t!=EOF){
				gall->add_link(id1,id2);
				if(gall->node[id1].node_props==NULL){
					gall->node[id1].node_props=new node_attr(i,1+log(b1));
					gall->node[id1].node_props->cluster_ID=1;
				}
				if(gall->node[id2].node_props==NULL){
					gall->node[id2].node_props=new node_attr(j,1+log(b2));
					gall->node[id2].node_props->cluster_ID=1;
				}
				gall->node[id1].set_link_w(id2,1+log(lw));	
			}
		}
		fclose(be);
		return(gall);
}

struct sgraph *read_linklist_linkweights_only(char di[],char ne[]){
	unsigned long long j,t,id1,id2;
	unsigned long long NR=0;
	char nev[300];
	sgraph *gall;
	FILE *be;
	double lw;
	
	sprintf(nev,"cd %s\nls\n",di);system(nev);
	
	sprintf(nev,"%s/%s.net",di,ne);
	
	be=fopen(nev,"r");
	j=1;
	while(j!=EOF){
		j=fscanf(be,"%lld %lld %lg",&id1,&id2,&lw);
		printf("%lld  %lld  %lg\n",id1,id2,lw);
		if(id1>NR) NR=id1;
		if(id2>NR) NR=id2;
	}
	fclose(be);
	
	gall=new class sgraph(NR,ne);	
	printf("NR=%lld\n",gall->N);//getchar();
	printf("Network in: %s\n",nev);
	be=fopen(nev,"r");
	t=1;
	while(t!=EOF){
		t=fscanf(be,"%lld %lld %lg",&id1,&id2,&lw);
		//printf("%d\n",id1);getchar();
		if(t!=EOF){
			gall->add_link(id1,id2);
			gall->node[id1].set_link_w(id2,lw);	
		}
	}
	fclose(be);
	return(gall);
}


void write_nD(sgraph *g,char path[300]){
	FILE *net;
	unsigned long long i;
	char nev[400];
	snode::slinklist::slink *sc2;
	int k;
	
	if(g->D<1){
		printf("Not a spatially rendered graph!\n");getchar();return;
	}
	
	sprintf(nev,"%s/%s.network_%dD",path,g->gname,g->D);
	printf("writing %s\n",nev);//getchar();
	net=fopen(nev,"w");
	for(i=1;i<=g->N;i++){
		fprintf(net,"%lld\t",i);
		for(k=0;k<g->D;k++)
			fprintf(net,"%lf ",g->node[i].node_dD_props->poz[k]);
		fprintf(net,"\n");
	}
			
	fprintf(net,"0\t");for(k=0;k<g->D;k++) fprintf(net,"0.0 ");fprintf(net,"\n");

	// indicates break between nodes and links for easy reading
	for(i=1;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;	
		while(sc2!=NULL){
			if(sc2->link_props==NULL)
				fprintf(net,"%lld\t%lld\t%lf\n",i,sc2->sneighbor,1.);
			else 
				fprintf(net,"%lld\t%lld\t%lf\n",i,sc2->sneighbor,sc2->link_props->lw);
			sc2=sc2->next(sc2);
		}
	}
	fclose(net);
}

void write_network_draw_attr(sgraph *g,char path[300]){
	FILE *net;
	unsigned long long i;
	char nev[400];
	snode::slinklist::slink *sc2;
		
	sprintf(nev,"%s/%s.labels",path,g->gname);
	printf("writing %s\n",nev);//getchar();
	net=fopen(nev,"w");
	fprintf(net,"%lld\n\n",g->N);
	for(i=1;i<=g->N;i++){
		if(g->node[i].node_draw_props!=NULL)
			fprintf(net,"%lld\t%lld\t%lld\t%lld\t%lld\t%lg\t%s\t%s\n",i,
				g->node[i].ID,
				g->node[i].node_draw_props->basic->ID,
				g->node[i].node_draw_props->basic->cluster_ID,
				g->node[i].node_draw_props->basic->basin,
				g->node[i].node_draw_props->basic->nw,
				export_string(g->node[i].node_draw_props->node_label),
				g->node[i].node_draw_props->color);
        else {if(g->node[i].node_props!=NULL)
                fprintf(net,"%lld\t%lld\t%lld\t%lld\n",i,
					 g->node[i].ID,
					 g->node[i].node_props->ID,
					 g->node[i].node_props->cluster_ID);
                else  fprintf(net,"%lld\t%lld\n",i,g->node[i].ID);
            }
	}

	// indicates break between nodes and links for easy reading
	for(i=1;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;	
		while(sc2!=NULL){
			if(sc2->link_props==NULL)
				fprintf(net,"%lld\t%lld\t%lg\n",i,sc2->sneighbor,1.);
			else 
				fprintf(net,"%lld\t%lld\t%lg\n",i,sc2->sneighbor,sc2->link_props->lw);
			sc2=sc2->next(sc2);
		}
	}
	fclose(net);
}

void write_weighted_network(sgraph *g,char path[300]){
	FILE *net;
	unsigned long long i;
	char nev[400];
	snode::slinklist::slink *sc2;
	
	sprintf(nev,"%s/%s.labels",path,g->gname);
	printf("writing %s\n",nev);//getchar();
	net=fopen(nev,"w");
	fprintf(net,"%lld\n\n",g->N);
	for(i=1;i<=g->N;i++){
		if(g->node[i].node_props!=NULL)
			fprintf(net,"%lld\t%lld\t%lld\t%lld\t%lld\t%lg\n",i,
					g->node[i].ID,
					g->node[i].node_props->ID,
					g->node[i].node_props->cluster_ID,
					g->node[i].node_props->basin,
					g->node[i].node_props->nw);
		else fprintf(net,"%lld\t%lld\n",i,g->node[i].ID);
	}
	
	// indicates break between nodes and links for easy reading
	for(i=1;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;	
		while(sc2!=NULL){
			if(sc2->link_props==NULL)
				fprintf(net,"%lld\t%lld\t%lg\n",i,sc2->sneighbor,1.);
			else 
				fprintf(net,"%lld\t%lld\t%lg\n",i,sc2->sneighbor,sc2->link_props->lw);
			sc2=sc2->next(sc2);
		}
	}
	fclose(net);
}

void CFinder_label_output(sgraph *g,char path[300],int lw_by_node){
	FILE *net;
	unsigned long long i,m;
	char nev[400];
	snode::slinklist::slink *sc2;
	
	sprintf(nev,"%s/%s.Cfind",path,g->gname);
	net=fopen(nev,"w");
	
	for(i=1;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;	
		while(sc2!=NULL){
			if(sc2->link_props==NULL)
				fprintf(net,"%s\t%s\n",export_string(g->node[i].node_draw_props->node_label),
						export_string(g->node[sc2->sneighbor].node_draw_props->node_label));
			else {
				if(lw_by_node==0)
					fprintf(net,"%s\t%s\t%lf\n",export_string(g->node[i].node_draw_props->node_label),
							export_string(g->node[sc2->sneighbor].node_draw_props->node_label),sc2->link_props->lw);
				else{
					if(g->node[i].node_draw_props->basic->nw<g->node[sc2->sneighbor].node_draw_props->basic->nw)
						m=g->node[i].node_draw_props->basic->nw;
					else m=g->node[sc2->sneighbor].node_draw_props->basic->nw;
					fprintf(net,"%s\t%s\t%lf\n",export_string(g->node[i].node_draw_props->node_label),
							export_string(g->node[sc2->sneighbor].node_draw_props->node_label),sc2->link_props->lw/(double)m);	
				}	
			}
			sc2=sc2->next(sc2);
		}
	}
	fclose(net);
}

struct sgraph *read_nD(char di[],char ne[],unsigned int d){
	unsigned long long i = 0,id1,id2;
	unsigned long long NR=0;
	int j,t,k;
	char nev[300];
	sgraph *gall;
	FILE *be;
	double lw,*x;
	
	sprintf(nev,"%s/%s.network_%dD",di,ne,d);
	printf("opening %s\n",nev);
	be=fopen(nev,"r");
	j=1;
	x=new double[d];
	while(j!=EOF){
		j=fscanf(be,"%lld",&id1);
		for(k=0;k<d;k++)
			fscanf(be,"%lf",&x[0]);
		if(id1>NR) NR=id1;
		if(id1==0) j=EOF;
	}
	fclose(be);
	
	gall=new class sgraph(NR,ne);
	gall->D=d;	
	printf("NR=%lld\n",gall->N);//getchar();
	printf("Network in: %s\n",nev);
	be=fopen(nev,"r");
	id1=1;
	while(id1!=0){
		t=fscanf(be,"%lld",&id1);
		for(k=0;k<d;k++)
			fscanf(be,"%lf",&x[k]);
		//printf("%d\n",id1);//getchar();
		if(id1!=0)
			if(gall->node[id1].node_dD_props==NULL){
				gall->node[id1].node_dD_props=new node_dD_attr(i,1,1);
				gall->node[id1].node_dD_props->set_node_poz(d,x);
				gall->node[id1].node_props=gall->node[id1].node_dD_props->basic;
			}
				//printf("%d\t%lf\n",id1,gall->node[id1].node_dD_props->poz[0]);getchar();	
	}
	j=1;
	while(j!=EOF){
		j=fscanf(be,"%lld %lld %lf",&id1,&id2,&lw);
		//printf("%d %d %lf\n",id1,id2,lw);getchar();
		if(j!=EOF){
			gall->add_link(id1,id2);
			gall->node[id1].set_link_w(id2,1+log(lw));	
		}
	}
	fclose(be);
	//printf("%d %d %lf\n",id1,id2,lw);getchar();
	delete[] x; x=NULL;		
	return(gall);
}

struct sgraph *read_network_draw_attr(char di[],const char ne[],int draw){
	unsigned long long i,id1,id2,ID,clust,basin,nid;
	unsigned int NR=0;
	int j;
	char nev[400];
	sgraph *gall;
	FILE *be;
	double lw,nw = 0.0;
	char label[501],color[100];
	
	sprintf(nev,"%s/%s.labels",di,ne);
	printf("opening %s\n",nev);
	be=fopen(nev,"r");
	if(be==NULL) return(NULL);
	fscanf(be,"%d",&NR);
		
	gall=new class sgraph(NR,ne);
	printf("NR=%lld\n",gall->N);//getchar();
//	printf("Network in: %s\n",nev);
	for(i=1;i<=gall->N;i++){
		if(draw==1){
			fscanf(be,"%lld\t%lld\t%lld\t%lld\t%lld\t%lg%s\t%s\n",&id1,&nid,&ID,&clust,&basin,&nw,label,color);
		//	printf("%lld\t%lld\t%lld\t%lld\t%lld\t%lg%s\t%s\n",id1,nid,ID,clust,basin,nw,label,color);
			gall->node[i].node_draw_props=new node_draw_attr(ID,nw,clust);
			gall->node[i].node_draw_props->basic->basin=basin;
			gall->node[i].node_draw_props->basic->nw=nw;
			sprintf(gall->node[i].node_draw_props->node_label,"%s",nice_string(label));
			sprintf(gall->node[i].node_draw_props->color,"%s",color);
			gall->node[i].node_props=gall->node[i].node_draw_props->basic;
			gall->node[i].ID=nid;
		}
		else{
			fscanf(be,"%lld\t%lld\t%lld\t%lld\n",&id1,&nid,&ID,&clust);
			gall->node[i].node_draw_props=new node_draw_attr(ID,1,clust);
			gall->node[i].node_draw_props->basic->nw=nw;
			gall->node[i].node_props=gall->node[i].node_draw_props->basic;
			gall->node[i].ID=nid;
		}
		/*if((i/100==i/100.)||(i==gall->N)||(ID==1)) {
			printf("%lld\t%lld\t%lld\t%lld\t%lg\t%s\t%s\n",i,
				gall->node[i].node_draw_props->basic->ID,
				gall->node[i].node_draw_props->basic->cluster_ID,
				gall->node[i].node_draw_props->basic->basin,
				gall->node[i].node_draw_props->basic->nw,
				gall->node[i].node_draw_props->node_label,
				gall->node[i].node_draw_props->color);
			//getchar();
		}*/
	}
	j=1;
	while(j!=EOF){
		j=fscanf(be,"%lld %lld %lf",&id1,&id2,&lw);
		//printf("%lld %lld %lf\n",id1,id2,lw);getchar();
		if(j!=EOF){
			gall->add_link(id1,id2);
			gall->node[id1].set_link_w(id2,lw);	
		}
	}
	fclose(be);
	//printf("%d %d %lf\n",id1,id2,lw);getchar();
	//delete[] x; x=NULL;		
	return(gall);
}

void write_largest_cluster_linklist_weights(sgraph *g,char path[100],graph_clusters *erd){
	 FILE *net;
	 int i;
	 char nev[400];
	 snode::slinklist::slink *sc2;
	 unsigned long long *sorszam,hany;
	 
	 sorszam=new unsigned long long[g->N+1];
	 hany=0;
	 for(i=1;i<=g->N;i++){
		if(g->node[i].node_props->cluster_ID==erd->max_ID){
			hany++;
			sorszam[i]=hany;
		}
		else sorszam[i]=0;
	 }
	sprintf(nev,"%s/%s_maxCL.network",path,g->gname);
	net=fopen(nev,"w");
	for(i=1;i<=g->N;i++)
	   if(sorszam[i]>0){
	 	if(g->node[i].node_props==NULL){
			printf("Node with no props! Wrong routine to write this network out.");
			getchar();return;
		}
		sc2=g->node[i].Li.first_link_pt;	
		while(sc2!=NULL){
			fprintf(net,"%llu %d %lf\t%llu %llu %lf\t%lf\n",
							sorszam[i],i,g->node[i].node_props->nw,
							sorszam[sc2->sneighbor],sc2->sneighbor,g->node[sc2->sneighbor].node_props->nw,
							sc2->link_props->lw);
			sc2=sc2->next(sc2);
		}
	 }
	fclose(net);
}

struct sgraph *read_largest_cluster_linklist_weights(char di[],char ne[]){
	unsigned long long i,j,t,id1,id2;
	unsigned long long NR=0;
	char nev[300];
	sgraph *gall;
	FILE *be;
	double lw,b1,b2;
	
	printf("Network in: %s\n",nev);//getchar();
	
	sprintf(nev,"%s/%s.network",di,ne);
	be=fopen(nev,"r");
	j=1;
	while(j!=EOF){
		j=fscanf(be,"%lld%lld%lf %lld%lld%lf %lf",&id1,&i,&b1,&id2,&t,&b2,&lw);
		if(id1>NR) NR=id1;
		if(id2>NR) NR=id2;
	}
	fclose(be);
		
	gall=new class sgraph(NR,ne);	
	//printf("NR=%ld\n",gall->N);getchar();
	//printf("Network in: %s\n",nev);
	be=fopen(nev,"r");
	t=1;
	while(t!=EOF){
		t=fscanf(be,"%lld%lld%lf %lld%lld%lf %lf",&id1,&i,&b1,&id2,&j,&b2,&lw);
		//printf("%d\n",id1);getchar();
		if(t!=EOF){
			gall->add_link(id1,id2);
			if(gall->node[id1].node_props==NULL){
				gall->node[id1].node_props=new node_attr(i,1+log(b1));
				gall->node[id1].node_props->cluster_ID=1;
			}
			if(gall->node[id2].node_props==NULL){
				gall->node[id2].node_props=new node_attr(j,1+log(b2));
				gall->node[id2].node_props->cluster_ID=1;
			}
			gall->node[id1].set_link_w(id2,1+log(lw));	
		}
	}
	fclose(be);
	return(gall);
}

struct sgraph *read_onecluster_linklist_NO_weights(char di[],char ne[]){
	unsigned long long j,t,id1,id2;
	unsigned long long NR=0;
	char nev[300];
	sgraph *gall;
	FILE *be;
	
		
		sprintf(nev,"%s/%s",di,ne);
		//printf("Network in: %s\n",nev);getchar();
		be=fopen(nev,"r");
		j=1;
		while(j!=EOF){
			//j=fscanf(be,"%ld%ld%d",&id1,&id2,&bla);
			j=fscanf(be,"%lld%lld",&id1,&id2);
			//printf("%ld %ld\n",id1,id2);getchar();
			if(id1>NR) NR=id1;
			if(id2>NR) NR=id2;
		}
		fclose(be);
		//printf("NR=%ld\n",NR);
		gall=new class sgraph(NR,ne);	
		//printf("NR=%ld\n",gall->N);getchar();
		//printf("Network in: %s\n",nev);
		be=fopen(nev,"r");
		t=1;
		while(t!=EOF){
			//t=fscanf(be,"%ld%ld%d",&id1,&id2,&bla);
			t=fscanf(be,"%lld%lld",&id1,&id2);
			//printf("%d\n",id1);getchar();
			if(t!=EOF){
				gall->add_link(id1,id2);
			}
		}
		fclose(be);
		//printf("NR=%ld\n",gall->N);
		return(gall);
}

struct sgraph *read_any_linklist_NO_weights(char di[],char ne[]){
	unsigned long long j,t,id1,id2,NR=0,NR_ID=0;
	typedef unsigned long long * plist;
	plist *nonzero;
	char nev[300];
	sgraph *gall;
	FILE *be;
	
	
	sprintf(nev,"%s/%s",di,ne);
	//printf("Network in: %s\n",nev);getchar();
	be=fopen(nev,"r");
	j=1;
	while(j!=EOF){
		//j=fscanf(be,"%ld%ld%d",&id1,&id2,&bla);
		j=fscanf(be,"%lld%lld",&id1,&id2);
			id1++;id2++;
		//printf("%lld %lld\n",id1,id2);getchar();
		if(id1>NR_ID) NR_ID=id1;
		if(id2>NR_ID) NR_ID=id2;
	}
	fclose(be);

	printf("NR_ID=%lld\n",NR_ID);

	nonzero=new plist[NR_ID+1];
	for(j=0;j<=NR_ID;j++) nonzero[j]=NULL;
	be=fopen(nev,"r");

	j=1;NR=0;
	while(j!=EOF){
		//j=fscanf(be,"%ld%ld%d",&id1,&id2,&bla);
		j=fscanf(be,"%lld%lld",&id1,&id2);
			//id1++;id2++;
		//printf("%ld %ld\n",id1,id2);getchar();
		if(nonzero[id1]==NULL){
			NR++;
			nonzero[id1]=new unsigned long long;
			*nonzero[id1]=NR;
		}
		if(nonzero[id2]==NULL){
			NR++;
			nonzero[id2]=new unsigned long long;
			*nonzero[id2]=NR;
		}
	}
	fclose(be);
	
	printf("NR=%llu\n",NR);
	gall=new class sgraph(NR,ne);	
	//printf("NR=%ld\n",gall->N);getchar();
	//printf("Network in: %s\n",nev);
	be=fopen(nev,"r");
	t=1;
	while(t!=EOF){
		//t=fscanf(be,"%ld%ld%d",&id1,&id2,&bla);
		t=fscanf(be,"%lld%lld",&id1,&id2);
			//id1++;id2++;
		//printf("%d\n",id1);getchar();
		if(t!=EOF){
			gall->add_link(*nonzero[id1],*nonzero[id2]);
			gall->node[*nonzero[id1]].ID=id1;
			gall->node[*nonzero[id2]].ID=id2;
		}
	//	printf("adding link %lld-> %lld (read in %lld -> %lld)\n",*nonzero[id1],*nonzero[id2],id1,id2);getchar();
	}
	fclose(be);
	//printf("NR=%ld\n",gall->N);
	for(j=0;j<=NR_ID;j++)
		if(nonzero[j]!=NULL) {
			delete[] nonzero[j];nonzero[j]=NULL;
		}
	delete[] nonzero; nonzero=NULL;
	
	return(gall);
}

void write_largest_cluster_linklist_weights_3D(sgraph *g,char path[300],unsigned int max_cl){
	 FILE *net;
	 char nev[400];
	 snode::slinklist::slink *sc2;
	 unsigned long long i,*sorszam,hany;
	 
	 if(g->D!=3){
		printf("Not a 3D rendered graph!\n");getchar();return;
	 }

	 sorszam=new unsigned long long[g->N+1];
	 hany=0;
	 for(i=1;i<=g->N;i++){
		if(g->node[i].node_props->cluster_ID==max_cl){
			hany++;
			sorszam[i]=hany;
		}
		else sorszam[i]=0;
	 }
	 sprintf(nev,"%s/%s.network_3D",path,g->gname);
	 printf("writing %s\n",nev);//getchar();
	 net=fopen(nev,"w");
	 for(i=1;i<=g->N;i++)
	   if(sorszam[i]>0){
			if(g->node[i].node_props==NULL){
				printf("Node with no props! Wrong routine to write this network out.");
				getchar();return;
			}
			//if(g->node[i].node_props->nw>0)
				fprintf(net,"%lld %lld %lld %lg\t%lg  %lg  %lg\n", sorszam[i],i,g->node[i].node_props->ID,g->node[i].node_props->nw,
						g->node[i].node_dD_props->poz[0],g->node[i].node_dD_props->poz[1],g->node[i].node_dD_props->poz[2]);
			//else 			
			//	fprintf(net,"%d %ld %lf\t%lf  %lf  %lf\n", sorszam[i],i,1.,
			//			g->node[i].node_dD_props->poz[0],g->node[i].node_dD_props->poz[1],g->node[i].node_dD_props->poz[2]);
	   }

	fprintf(net,"0 0 0.0\t0.0  0.0  0.0\n");
	// indicates break between nodes and links for easy reading
	 for(i=1;i<=g->N;i++)
	   if(sorszam[i]>0){
		sc2=g->node[i].Li.first_link_pt;	
		while(sc2!=NULL){
			if(sc2->link_props==NULL)
				fprintf(net,"%lld\t%lld\t%lg\n",sorszam[i],sorszam[sc2->sneighbor],1.);
			else 
				fprintf(net,"%lld\t%lld\t%lg\n",sorszam[i],sorszam[sc2->sneighbor],sc2->link_props->lw);
			sc2=sc2->next(sc2);
		}
	 }
	fclose(net);
}

void write_largest_cluster_linklist_weights_3D_nodelabel(sgraph *g,char path[300],unsigned int max_cl){
	FILE *net;
	char nev[400];
	snode::slinklist::slink *sc2;
	unsigned long long i,*sorszam,hany;
	
	if(g->D!=3){
		printf("Not a 3D rendered graph!\n");getchar();return;
	}
	
	sorszam=new unsigned long long[g->N+1];
	hany=0;
	for(i=1;i<=g->N;i++){
		if(g->node[i].node_props->cluster_ID==max_cl){
			hany++;
			sorszam[i]=hany;
		}
		else sorszam[i]=0;
	}
	sprintf(nev,"%s/%s.network_3D",path,g->gname);
	printf("writing %s\n",nev);//getchar();
	net=fopen(nev,"w");
	for(i=1;i<=g->N;i++)
		if(sorszam[i]>0){
			if(g->node[i].node_props==NULL){
				printf("Node with no props! Wrong routine to write this network out.");
				getchar();return;
			}
			//if(g->node[i].node_props->nw>0)
			fprintf(net,"%lld %lld %lld %lg  %s\t%lg  %lg  %lg\n", sorszam[i],i,
					g->node[i].node_props->ID,g->node[i].node_props->nw,
					g->node[i].node_draw_props->node_label,
					g->node[i].node_dD_props->poz[0],g->node[i].node_dD_props->poz[1],g->node[i].node_dD_props->poz[2]);
			//else 			
			//	fprintf(net,"%d %ld %lf\t%lf  %lf  %lf\n", sorszam[i],i,1.,
			//			g->node[i].node_dD_props->poz[0],g->node[i].node_dD_props->poz[1],g->node[i].node_dD_props->poz[2]);
		}
	
	fprintf(net,"0 0 0.0 __\t0.0  0.0  0.0\n");
	// indicates break between nodes and links for easy reading
	for(i=1;i<=g->N;i++)
		if(sorszam[i]>0){
			sc2=g->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				if(sc2->link_props==NULL)
					fprintf(net,"%lld\t%lld\t%lg\n",sorszam[i],sorszam[sc2->sneighbor],1.);
				else 
					fprintf(net,"%lld\t%lld\t%lg\n",sorszam[i],sorszam[sc2->sneighbor],sc2->link_props->lw);
				sc2=sc2->next(sc2);
			}
		}
	fclose(net);
}

void write_largest_cluster_3D(sgraph *g,char path[300],unsigned int max_cl){
	FILE *net;
	unsigned long long i;
	char nev[400];
	snode::slinklist::slink *sc2;
	unsigned long long *sorszam,hany;
	
	if(g->D!=3){
		printf("Not a 3D rendered graph!\n");getchar();return;
	}
	
	sorszam=new unsigned long long[g->N+1];
	hany=0;
	for(i=1;i<=g->N;i++){
		if((g->node[i].node_props->cluster_ID==max_cl)||(max_cl==0)){
			hany++;
			sorszam[i]=hany;
		}
		else sorszam[i]=0;
	}
	sprintf(nev,"%s/%s.network_3D",path,g->gname);
	printf("writing %s\n",nev);//getchar();
		net=fopen(nev,"w");
		for(i=1;i<=g->N;i++)
			if(sorszam[i]>0){
				fprintf(net,"%lld %lld\t%lf  %lf  %lf\n", sorszam[i],i,
						g->node[i].node_dD_props->poz[0],g->node[i].node_dD_props->poz[1],g->node[i].node_dD_props->poz[2]);
				//else 			
				//	fprintf(net,"%d %ld %lf\t%lf  %lf  %lf\n", sorszam[i],i,1.,
				//			g->node[i].node_dD_props->poz[0],g->node[i].node_dD_props->poz[1],g->node[i].node_dD_props->poz[2]);
			}
				
				fprintf(net,"0 0\t0.0  0.0  0.0\n");
		// indicates break between nodes and links for easy reading
		for(i=1;i<=g->N;i++)
			if(sorszam[i]>0){
				sc2=g->node[i].Li.first_link_pt;	
				while(sc2!=NULL){
					if(sc2->link_props==NULL)
						fprintf(net,"%lld\t%lld\t%lf\n",sorszam[i],sorszam[sc2->sneighbor],1.);
					else 
						fprintf(net,"%lld\t%lld\t%lf\n",sorszam[i],sorszam[sc2->sneighbor],sc2->link_props->lw);
					sc2=sc2->next(sc2);
				}
			}
		fclose(net);
}

struct sgraph *read_largest_cluster_linklist_weights_3D(char di[],char ne[]){
	unsigned long long i,id1,id2,nodeid,NR=0,t;
	int j;
	char nev[300];
	sgraph *gall;
	FILE *be;
	double lw,b1,x[3];
	
	sprintf(nev,"%s/%s.network_3D",di,ne);
	printf("Network in: %s\n",nev);//getchar();
	be=fopen(nev,"r");
	j=1;
	while(j!=EOF){
		j=fscanf(be,"%lld%lld%lld%lg %lg%lg%lg",&id1,&i,&nodeid,&b1,&x[0],&x[1],&x[2]);
		if(id1>NR) NR=id1;
		if(id1==0) j=EOF;
	}
	fclose(be);
		
	gall=new class sgraph(NR,ne);
	gall->D=3;	
	printf("NR=%lld\n",gall->N);//getchar();
	printf("Network in: %s\n",nev);
	be=fopen(nev,"r");
	id1=1;
	while(id1!=0){
		t=fscanf(be,"%lld%lld%lld%lg %lg%lg%lg",&id1,&i,&nodeid,&b1,&x[0],&x[1],&x[2]);
		//printf("%d\n",id1);//getchar();
		if(id1!=0)
			if(gall->node[id1].node_dD_props==NULL){
			   gall->node[id1].node_dD_props=new node_dD_attr(i,b1,1);
			   gall->node[id1].node_dD_props->set_node_poz(3,x);
			   gall->node[id1].node_props=gall->node[id1].node_dD_props->basic;
			   gall->node[id1].node_props->ID=nodeid;
			}
		//printf("%d\t%lf\n",id1,gall->node[id1].node_dD_props->poz[0]);getchar();	
	}
	j=1;
	while(j!=EOF){
		j=fscanf(be,"%lld %lld %lg",&id1,&id2,&lw);
		//printf("%d %d %lf\n",id1,id2,lw);getchar();
		if(j!=EOF){
			gall->add_link(id1,id2);
			//gall->node[id1].set_link_w(id2,1+log(lw));
			gall->node[id1].set_link_w(id2,lw);
		}
	}
	fclose(be);
	//printf("%d %d %lf\n",id1,id2,lw);getchar();
	return(gall);
}

struct sgraph *read_largest_cluster_linklist_weights_3D_nodelabel(char di[],char ne[]){
	unsigned long long i,id1,id2,nodeid,NR=0,t;
	int j;
	char nev[300],lab[100];
	sgraph *gall;
	FILE *be;
	double lw,b1,x[3];
	
	sprintf(nev,"%s/%s.network_3D",di,ne);
	printf("Network in: %s\n",nev);//getchar();
	be=fopen(nev,"r");
	j=1;
	while(j!=EOF){
		j=fscanf(be,"%lld %lld %lld %lg  %s\t%lg  %lg  %lg\n",&id1,&i,&nodeid,&b1,lab,&x[0],&x[1],&x[2]);
		if(id1>NR) NR=id1;
		if(id1==0) j=EOF;
	}
	fclose(be);
	
	gall=new class sgraph(NR,ne);
	gall->D=3;	
	printf("NR=%lld\n",gall->N);//getchar();
	printf("Network in: %s\n",nev);
	be=fopen(nev,"r");
	id1=1;
	while(id1!=0){
		t=fscanf(be,"%lld %lld %lld %lg  %s\t%lg  %lg  %lg\n",&id1,&i,&nodeid,&b1,lab,&x[0],&x[1],&x[2]);
		//printf("%d\n",id1);//getchar();
		if(id1!=0)
			if(gall->node[id1].node_dD_props==NULL){
				gall->node[id1].node_dD_props=new node_dD_attr(i,b1,1);
				gall->node[id1].node_dD_props->set_node_poz(3,x);
				gall->node[id1].node_props=gall->node[id1].node_dD_props->basic;
				gall->node[id1].node_props->ID=nodeid;
				gall->node[id1].node_draw_props=new node_draw_attr(nodeid,b1,1);
				sprintf(gall->node[id1].node_draw_props->node_label,"%s",lab);
				//printf("%lld\t%s\n",id1,gall->node[id1].node_draw_props->node_label);//getchar();	
			}
	}
	j=1;
	while(j!=EOF){
		j=fscanf(be,"%lld %lld %lg",&id1,&id2,&lw);
		//printf("%lld %lld %lg\n",id1,id2,lw);getchar();
		if(j!=EOF){
			gall->add_link(id1,id2);
			//gall->node[id1].set_link_w(id2,1+log(lw));
			gall->node[id1].set_link_w(id2,lw);
		}
	}
	fclose(be);
	//printf("%d %d %lf\n",id1,id2,lw);getchar();
	return(gall);
}

struct sgraph *read_largest_cluster_3D(char di[],char ne[]){
	unsigned long long i,id1,id2;
	unsigned long long NR=0;
	int j,t;
	char nev[300];
	sgraph *gall;
	FILE *be;
	double lw,x[3];
	
	sprintf(nev,"%s/%s.network_3D",di,ne);
//	printf("Network in: %s\n",nev);//getchar();
		be=fopen(nev,"r");
		j=1;
		while(j!=EOF){
			j=fscanf(be,"%lld%lld %lf%lf%lf",&id1,&i,&x[0],&x[1],&x[2]);
			if(id1>NR) NR=id1;
			if(id1==0) j=EOF;
		}
		fclose(be);
		
		gall=new class sgraph(NR,ne);
		gall->D=3;	
		printf("NR=%lld\n",gall->N);//getchar();
			printf("Network in: %s\n",nev);
			be=fopen(nev,"r");
			id1=1;
			while(id1!=0){
				t=fscanf(be,"%lld%lld %lf%lf%lf",&id1,&i,&x[0],&x[1],&x[2]);
				//printf("%d\n",id1);//getchar();
				if(id1!=0)
					if(gall->node[id1].node_dD_props==NULL){
						gall->node[id1].node_dD_props=new node_dD_attr(i,1,1);
						gall->node[id1].node_dD_props->set_node_poz(3,x);
						gall->node[id1].node_props=gall->node[id1].node_dD_props->basic;
					}
						//printf("%d\t%lf\n",id1,gall->node[id1].node_dD_props->poz[0]);getchar();	
			}
				j=1;
				while(j!=EOF){
					j=fscanf(be,"%lld %lld %lf",&id1,&id2,&lw);
					//printf("%d %d %lf\n",id1,id2,lw);getchar();
					if(j!=EOF){
						gall->add_link(id1,id2);
						gall->node[id1].set_link_w(id2,1+log(lw));	
					}
				}
				fclose(be);
				//printf("%d %d %lf\n",id1,id2,lw);getchar();
				return(gall);
}

void write_Grad_over_ER_Dot(sgraph *er,sgraph *gr,char dotname[100],double rain){
	 FILE *net;
	 int i;
	 snode::slinklist::slink *sc2;
	 if(er->N>500) {printf("Graph too large for GraphViz picture\n");return;}
	 
	 net=fopen(dotname,"w");fprintf(net,"graph G {\n	graph [overlap=false];\n"); 
	 for(i=1;i<=er->N;i++){
		if(gr->node[i].Li.get_dirlink(i))
			{if(er->node[i].node_props->nw<=rain)
			   {er->node[i].node_draw_props->set_node_rgb(0,(unsigned int)((1.-er->node[i].node_props->nw)*255),0);
			   fprintf(net,"%d [fillcolor=\"%s\",style=filled,label=\"M-%s\"];\n",
					 i,er->node[i].node_draw_props->node_rgb.RGB_hex,er->node[i].node_draw_props->node_label);	
			  //printf("node %d has props: %s\n",sc1->ID,sc1->node_props->node_rgb.RGB_hex); getchar();
			   }
			 else fprintf(net,"%d [fillcolor=\"Grey\",style=filled,label=\"M-%s\"];\n",
					 i,er->node[i].node_draw_props->node_label);  
			}
		else {if(er->node[i].node_props->nw<=rain){
				er->node[i].node_draw_props->set_node_rgb((unsigned int)((1.-er->node[i].node_props->nw)*255),0,0);
			    fprintf(net,"%d [fillcolor=\"%s\",style=filled,label=\"%s\"];\n",
					i,er->node[i].node_draw_props->node_rgb.RGB_hex,er->node[i].node_draw_props->node_label);	
				}
			  else fprintf(net,"%d [fillcolor=\"White\",style=filled,label=\"%s\"];\n",
					 i,er->node[i].node_draw_props->node_label);  
			  //printf("node %d has props: %s\n",sc1->ID,sc1->node_props->node_rgb.RGB_hex); getchar();
			 }
	 }
	 for(i=1;i<=er->N;i++){
		sc2=er->node[i].Li.first_link_pt;//printf("%d ",sc2->sneighbor->ID);getchar();
		//printf("\nlinks of %d:\t",sc1->ID);
		//sc1->links_of_node.list_node_linklist();	
		while(sc2!=NULL){
			if(gr->linked(i,sc2->sneighbor)) 
				fprintf(net,"\t%d -- %llu [color=\"Green\"];\n",i,sc2->sneighbor);
			else {if( (er->node[i].node_props->nw<=rain)
				       &&(er->node[sc2->sneighbor].node_props->nw<=rain))
						{fprintf(net,"\t%d -- %llu [color=\"Black\"];\n",i,sc2->sneighbor);}
				     else fprintf(net,"\t%d -- %llu [color=\"White\"];\n",i,sc2->sneighbor);
			}
		sc2=sc2->next(sc2);//getchar();
		}
	}	
	fprintf(net,"}"); fflush(net);
	fclose(net);
}

void visualize_node_neighborhood(sgraph *g,longint_ARRAY *gl,const char cutname[],const char community_color[],const char path[],charpointer *genename){
    FILE *net;
    unsigned long long i;
    char nev[500];
    snode::slinklist::slink *sc2;
    double B=0.3;
    longint_ARRAY *gl_firstneighbors;
    
    sprintf(nev,"%s/%s_%s.eps",path,g->gname,cutname);
    net=fopen(nev,"w");
    setcolors_PAJEK_grow();
    int W=1000;
    EpsInit(net,0,0,W*(1+2*B),W*(1+2*B));
    
    double x_min=g->node[gl->A[1]].node_draw_props->poz[0];
    double x_max=g->node[gl->A[1]].node_draw_props->poz[0];
    double y_min=g->node[gl->A[1]].node_draw_props->poz[1];
    double y_max=g->node[gl->A[1]].node_draw_props->poz[1];
    
    
    int C_max=0;
    for (unsigned int ii=1; ii<=gl->N; ii++) {
        if (g->node[gl->A[ii]].node_draw_props->poz[0]<x_min) x_min=g->node[gl->A[ii]].node_draw_props->poz[0];
        if (g->node[gl->A[ii]].node_draw_props->poz[0]>x_max) x_max=g->node[gl->A[ii]].node_draw_props->poz[0];
        if (g->node[gl->A[ii]].node_draw_props->poz[1]<y_min) y_min=g->node[gl->A[ii]].node_draw_props->poz[1];
        if (g->node[gl->A[ii]].node_draw_props->poz[1]>y_max) y_max=g->node[gl->A[ii]].node_draw_props->poz[1];
    }
    printf("X: %lg to %lg\nY: %lg to %lg\n",x_min,x_max,y_min,y_max);
    EpsSetLinewidth(net,0.5);
    EpsSetFont(net,"Times-Bold",20);
    EpsSetRgb(net,0.7,0.7,0.7);
    double sx=x_max-x_min;double sy=y_max-y_min;
    
    
    int ok1,ok2;
    unsigned long int a,c;
    gl_firstneighbors=new longint_ARRAY();
    
    for(i=1;i<=g->N;i++){
        if (  (g->node[i].node_draw_props->poz[0]>=x_min-B*W)&&(g->node[i].node_draw_props->poz[0]<=x_max+B*W)
            &&(g->node[i].node_draw_props->poz[1]>=y_min-B*W)&&(g->node[i].node_draw_props->poz[1]<=y_max+B*W)){
            sc2=g->node[i].Li.first_link_pt;
            //printf("%lld ",sc2->sneighbor);getchar();
            //printf("\nlinks of %d:\t",sc1->ID);
            //sc1->links_of_node.list_node_linklist();
            if (g->node[i].node_props->cluster_ID>C_max)  C_max=(int)g->node[i].node_props->cluster_ID;
            a=gl->check_for_element(i);
            if(a>0)  ok1=1; else ok1=0;
            while(sc2!=NULL){
                //printf("1: %d ",sc2->sneighbor->ID);getchar();
                if (  (g->node[sc2->sneighbor].node_draw_props->poz[0]>=x_min-B*W)
                    &&(g->node[sc2->sneighbor].node_draw_props->poz[0]<=x_max+B*W)
                    &&(g->node[sc2->sneighbor].node_draw_props->poz[1]>=y_min-B*W)
                    &&(g->node[sc2->sneighbor].node_draw_props->poz[1]<=y_max+B*W)){
                   
                    if (g->node[sc2->sneighbor].node_props->cluster_ID>C_max)
                        C_max=(int)g->node[sc2->sneighbor].node_props->cluster_ID;
                    if (ok1) gl_firstneighbors->add_element(sc2->sneighbor);
                    
                    c=gl->check_for_element(sc2->sneighbor);
                    if(c>0) {ok2=1; }
                        else ok2=0;
                   
                    if ((ok1)||(ok2)) {EpsDrawLine(net,B*W+W*(g->node[i].node_draw_props->poz[0]-x_min)/sx,
                                    B*W+W*(g->node[i].node_draw_props->poz[1]-y_min)/sy,
                                    B*W+W*(g->node[sc2->sneighbor].node_draw_props->poz[0]-x_min)/sx,
                                    B*W+W*(g->node[sc2->sneighbor].node_draw_props->poz[1]-y_min)/sy);
                        
                        // printf("Link %ld  %ld\n",a,c);//getchar();
                        }
                }
                sc2=sc2->next(sc2);//getchar();
            }
        }
    }
    
    double *r,*gr,*b;
    
    r=new double[C_max+1];
    gr=new double[C_max+1];
    b=new double[C_max+1];
    for (int i=0; i<=C_max; i++) {r[i]=0;gr[i]=0;b[i]=0;}
    FILE *f=fopen(community_color,"r");
    int j;
    
    for (int i=0; i<=C_max; i++) fscanf(f,"%d%lg%lg%lg",&j,&(r[i]),&(gr[i]),&(b[i]));
    fclose(f);
    
    double R;
    
    EpsSetRgb(net,0,0,0);
    R=2;
    
    for(i=1;i<=g->N;i++){
        if (  (g->node[i].node_draw_props->poz[0]>=x_min-B*W)&&(g->node[i].node_draw_props->poz[0]<=x_max+B*W)
            &&(g->node[i].node_draw_props->poz[1]>=y_min-B*W)&&(g->node[i].node_draw_props->poz[1]<=y_max+B*W)){
       
            EpsSetRgb(net,0.5+r[g->node[i].node_props->cluster_ID]/2.,0.5+gr[g->node[i].node_props->cluster_ID]/2.,0.5+b[g->node[i].node_props->cluster_ID]/2.);
            EpsFillCircle(net,B*W+ W*(g->node[i].node_draw_props->poz[0]-x_min)/sx,
                              B*W+ W*(g->node[i].node_draw_props->poz[1]-y_min)/sy,
                              R);
        }
    }
    printf("N_firstneighbors =%ld\n",gl_firstneighbors->N);//getchar();
    R=4;
    for(i=1;i<=g->N;i++){
        if (  (g->node[i].node_draw_props->poz[0]>=x_min-B*W)&&(g->node[i].node_draw_props->poz[0]<=x_max+B*W)
            &&(g->node[i].node_draw_props->poz[1]>=y_min-B*W)&&(g->node[i].node_draw_props->poz[1]<=y_max+B*W)){
            
            a=gl->check_for_element(i);
            c=gl_firstneighbors->check_for_element(i);
            if ((a)||(c)) {
                EpsSetRgb(net,r[g->node[i].node_props->cluster_ID],gr[g->node[i].node_props->cluster_ID],b[g->node[i].node_props->cluster_ID]);
                EpsFillCircle(net,B*W+ W*(g->node[i].node_draw_props->poz[0]-x_min)/sx,
                          B*W+ W*(g->node[i].node_draw_props->poz[1]-y_min)/sy,
                          R);
            }
        }
    }
    EpsSetRgb(net,0,0,0);
   
    for(i=1;i<=gl->N;i++){
      //  EpsSetRgb(net,r[g->node[gl->A[i]].node_props->cluster_ID],gr[g->node[gl->A[i]].node_props->cluster_ID],b[g->node[gl->A[i]].node_props->cluster_ID]);
        
        EpsDrawCircle(net,  B*W+ W*(g->node[gl->A[i]].node_draw_props->poz[0]-x_min)/sx,
                            B*W+ W*(g->node[gl->A[i]].node_draw_props->poz[1]-y_min)/sy, g->node[gl->A[i]].k);
        EpsDrawString(net,0,B*W+ W*(g->node[gl->A[i]].node_draw_props->poz[0]-x_min)/sx,
                            B*W+ W*(g->node[gl->A[i]].node_draw_props->poz[1]-y_min)/sy,genename[i]);
    }

    sprintf(nev,"%s/%s_%s",path,g->gname,cutname);
    close_EPSDRAW_file(net,nev,2);

}


void visualize_dense_network_2D(sgraph *g,char community_color[],char path[]){
    FILE *net;
    unsigned long long i;
    char nev[500];
    snode::slinklist::slink *sc2;
    double B=0.4;
    
    sprintf(nev,"%s/%s.eps",path,g->gname);
    net=fopen(nev,"w");
    setcolors_PAJEK_grow();
    int W=1000;
    EpsInit(net,0,0,W*(1+2*B),W*(1+2*B));
    
    double x_min=g->node[1].node_draw_props->poz[0];
    double x_max=g->node[1].node_draw_props->poz[0];
    double y_min=g->node[1].node_draw_props->poz[1];
    double y_max=g->node[1].node_draw_props->poz[1];
    
    
    int C_max=0;
    for (unsigned int i=1; i<=g->N; i++) {
        if (g->node[i].node_draw_props->poz[0]<x_min) x_min=g->node[i].node_draw_props->poz[0];
        if (g->node[i].node_draw_props->poz[0]>x_max) x_max=g->node[i].node_draw_props->poz[0];
        if (g->node[i].node_draw_props->poz[1]<y_min) y_min=g->node[i].node_draw_props->poz[1];
        if (g->node[i].node_draw_props->poz[1]>y_max) y_max=g->node[i].node_draw_props->poz[1];
    }
    printf("X: %lg to %lg\nY: %lg to %lg\n",x_min,x_max,y_min,y_max);
    EpsSetLinewidth(net,0.5);
    EpsSetFont(net,"Times-Bold",20);
    EpsSetRgb(net,0.7,0.7,0.7);
    double sx=x_max-x_min;double sy=y_max-y_min;
    
    
    for(i=1;i<=g->N;i++){
        if (  (g->node[i].node_draw_props->poz[0]>=x_min-B*W)&&(g->node[i].node_draw_props->poz[0]<=x_max+B*W)
            &&(g->node[i].node_draw_props->poz[1]>=y_min-B*W)&&(g->node[i].node_draw_props->poz[1]<=y_max+B*W)){
            sc2=g->node[i].Li.first_link_pt;
            //printf("%lld ",sc2->sneighbor);getchar();
            //printf("\nlinks of %d:\t",sc1->ID);
            //sc1->links_of_node.list_node_linklist();
            if (g->node[i].node_props->cluster_ID>C_max)  C_max=(int)g->node[i].node_props->cluster_ID;
            while(sc2!=NULL){
                //printf("1: %d ",sc2->sneighbor->ID);getchar();
                if (  (g->node[sc2->sneighbor].node_draw_props->poz[0]>=x_min-B*W)
                    &&(g->node[sc2->sneighbor].node_draw_props->poz[0]<=x_max+B*W)
                    &&(g->node[sc2->sneighbor].node_draw_props->poz[1]>=y_min-B*W)
                    &&(g->node[sc2->sneighbor].node_draw_props->poz[1]<=y_max+B*W)){
                    
                    if (g->node[sc2->sneighbor].node_props->cluster_ID>C_max)
                        C_max=(int)g->node[sc2->sneighbor].node_props->cluster_ID;
                    
                    /* EpsDrawLine(net,B*W+W*(g->node[i].node_draw_props->poz[0]-x_min)/sx,
                     B*W+W*(g->node[i].node_draw_props->poz[1]-y_min)/sy,
                     B*W+W*(g->node[sc2->sneighbor].node_draw_props->poz[0]-x_min)/sx,
                     B*W+W*(g->node[sc2->sneighbor].node_draw_props->poz[1]-y_min)/sy);
                     */
                }
                sc2=sc2->next(sc2);//getchar();
            }
        }
    }
    
    double *r,*gr,*b;
    
    r=new double[C_max+1];
    gr=new double[C_max+1];
    b=new double[C_max+1];
    for (int i=0; i<=C_max; i++) {r[i]=0;gr[i]=0;b[i]=0;}
    FILE *f=fopen(community_color,"r");
    int j;
    
    for (int i=0; i<=C_max; i++) fscanf(f,"%d%lg%lg%lg",&j,&(r[i]),&(gr[i]),&(b[i]));
    fclose(f);
    
    double R;
    
    EpsSetRgb(net,0,0,0);
    R=2;
    
    for(i=1;i<=g->N;i++){
        if (  (g->node[i].node_draw_props->poz[0]>=x_min-B*W)&&(g->node[i].node_draw_props->poz[0]<=x_max+B*W)
            &&(g->node[i].node_draw_props->poz[1]>=y_min-B*W)&&(g->node[i].node_draw_props->poz[1]<=y_max+B*W)){
            
            EpsSetRgb(net,r[g->node[i].node_props->cluster_ID],gr[g->node[i].node_props->cluster_ID],b[g->node[i].node_props->cluster_ID]);
            EpsFillCircle(net,B*W+ W*(g->node[i].node_draw_props->poz[0]-x_min)/sx,
                          B*W+ W*(g->node[i].node_draw_props->poz[1]-y_min)/sy,
                          R);
        }
    }
    EpsSetRgb(net,0,0,0);
    

    sprintf(nev,"%s/%s",path,g->gname);
    close_EPSDRAW_file(net,nev,2);
 }



/*
void write_Undir_Circle_EPS(sgraph *g,char dotname[100]){
	 FILE *net;
	 long long i;
	 snode::slinklist::slink *sc2;

	 if(g->N>500) {printf("Graph too large for GraphViz picture\n");return;}
	 net=fopen(dotname,"w");
	 setcolors_PAJEK();
	 EpsInit(net,-10,-10,4*g->N+10,4*g->N+10);
	 EpsSetLinewidth(net,1);
	 EpsSetRgb(net,0,0,0);
	 for(i=1;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;//printf("%d ",sc2->sneighbor->ID);getchar();
		//printf("\nlinks of %d:\t",sc1->ID);
		//sc1->links_of_node.list_node_linklist();	
		while(sc2!=NULL){
			//printf("1: %d ",sc2->sneighbor->ID);getchar();
			if(sc2->link_props==NULL) 
				EpsDrawLine(net,2*g->N+2*g->N*cos(i*2*M_PI/(double)g->N),2*g->N+2*g->N*sin(i*2.*M_PI/(double)g->N),
							    2*g->N+2*g->N*cos(sc2->sneighbor*2.*M_PI/(double)g->N),2*g->N+2*g->N*sin(sc2->sneighbor*2.*M_PI/(double)g->N));
				else {EpsSetLinewidth(net,sc2->link_props->lw);
					  EpsSetColor_PAJEK(net,sc2->link_props->color);
					 }
			sc2=sc2->next(sc2);//getchar();
		}
	}	

	for(i=1;i<=g->N;i++){
		if(g->node[i].node_draw_props==NULL){EpsSetRgb(net,1,0,0);
			EpsFillCircle(net,2*g->N+2*g->N*cos(i*2.*M_PI/(double)g->N),2*g->N+2*g->N*sin(i*2.*M_PI/(double)g->N),5);
		}
		else { EpsSetRgb(net,g->node[i].node_draw_props->node_rgb.r/255.,g->node[i].node_draw_props->node_rgb.g/255.,g->node[i].node_draw_props->node_rgb.b/255.);
			   EpsFillCircle(net,2*g->N+2*g->N*cos(i*2.*M_PI/(double)g->N),2*g->N+2*g->N*sin(i*2.*M_PI/(double)g->N),5);
		}
		EpsSetRgb(net,0,0,0);
		EpsDrawCircle(net,2*g->N+2*g->N*cos(i*2.*M_PI/(double)g->N),2*g->N+2*g->N*sin(i*2.*M_PI/(double)g->N),5);
	 }
	EpsClose(net);
	//	printf("End.\n");//getchar();
}

void draw_RG(sgraph *g,double R,char path[]){
	 FILE *net;
	 unsigned long long i,k_min,k_max;
	 char nev[500];
	 snode::slinklist::slink *sc2;
	 sgraph *e;

	 sprintf(nev,"%s/%s.eps",path,g->gname);	
	 net=fopen(nev,"w");
	 setcolors_PAJEK_grow();
	 EpsInit(net,0,0,500,500);
	 
	 EpsSetLinewidth(net,0.5);
	 EpsSetRgb(net,0.7,0.7,0.7);
	 e=gradiens_network_undir(g);
		 
	 for(i=1;i<=e->N;i++){
		sc2=e->node[i].Li.first_link_pt;
		 //printf("%lld ",sc2->sneighbor);getchar();
		//printf("\nlinks of %d:\t",sc1->ID);
		//sc1->links_of_node.list_node_linklist();	
		while(sc2!=NULL){
			//printf("1: %d ",sc2->sneighbor->ID);getchar();
			if(   (fabs(g->node[i].node_dD_props->poz[0]-g->node[sc2->sneighbor].node_dD_props->poz[0])<=0.5)
				&&(fabs(g->node[i].node_dD_props->poz[1]-g->node[sc2->sneighbor].node_dD_props->poz[1])<=0.5)
				//&&(labs(g->node[i].k-g->node[sc2->sneighbor].k)==1)
				  )
					{EpsDrawLine(net,500*g->node[i].node_dD_props->poz[0],500*g->node[i].node_dD_props->poz[1],
								500*g->node[sc2->sneighbor].node_dD_props->poz[0],500*g->node[sc2->sneighbor].node_dD_props->poz[1]);
			}
			//else	{EpsSetRgb(net,0.8,0,0);EpsSetLinewidth(net,0.25);
			//		 EpsDrawLine(net,500*g->node[i].node_dD_props->poz[0],500*g->node[i].node_dD_props->poz[1],
			//					500*g->node[sc2->sneighbor].node_dD_props->poz[0],500*g->node[sc2->sneighbor].node_dD_props->poz[1]);
			//}
			sc2=sc2->next(sc2);//getchar();
		}
	}	
	delete e; e=NULL;
	
	k_min=100000;
	k_max=0;
	for(i=1;i<=g->N;i++){
		if(k_min>g->node[i].k) k_min=g->node[i].k;
		if(k_max<g->node[i].k) k_max=g->node[i].k;
	}
	//printf("k_min=%lld\tk_max=%lld\n",k_min,k_max);getchar();
	for(i=1;i<=g->N;i++){
			//printf("k=%lld (%lld-%lld)\tcolor nr %d\n",g->node[i].k,k_min,k_max,(unsigned int)(12*(g->node[i].k-k_min)/(k_max-k_min)));
			EpsSetColor_PAJEK(net,szinnev[(int)(13*(g->node[i].k-k_min)/(k_max-k_min))]);
			EpsFillCircle(net,500*g->node[i].node_dD_props->poz[0],
							  500*g->node[i].node_dD_props->poz[1],
							  500*R);
							  
	}
	EpsClose(net);
	net=fopen("last_color_bar.eps","w");
	EpsInit(net,0,0,120,270);
	EpsSetLinewidth(net,1);
	EpsSetFont(net,"Times-Bold",12);
	for(i=0;i<13;i++){
		EpsSetColor_PAJEK(net,szinnev[i]);
		EpsFillRectangle(net,0,(13-i)*20,50,(13-i)*20-15);
		sprintf(nev,"%d - %d",(int)(k_min+(i)*(k_max-k_min)/13.),(int)(k_min+(i+1)*(k_max-k_min)/13.));
		EpsSetRgb(net,0,0,0);
		EpsDrawString(net,0,60,(13-i)*20-12,nev);
	}			
	EpsClose(net);

	printf("End.\n");//getchar();
}

void draw_RG_basins(sgraph *g,double R,char path[]){
	FILE *net;
	char nev[500];
	snode::slinklist::slink *sc2;
	sgraph *e;
	graph_clusters *erdo;
	unsigned long long i;
	
	sprintf(nev,"%s/%s_basins.eps",path,g->gname);	
	net=fopen(nev,"w");
	setcolors_PAJEK_grow();
	EpsInit(net,0,0,500,500);
	
	EpsSetLinewidth(net,0.5);
	EpsSetRgb(net,0.8,0.8,0.8);
	e=gradiens_network_undir(g);
	
	for(i=1;i<=e->N;i++){
		sc2=e->node[i].Li.first_link_pt;
		//printf("%lld ",sc2->sneighbor);getchar();
		//printf("\nlinks of %d:\t",sc1->ID);
		//sc1->links_of_node.list_node_linklist();	
		while(sc2!=NULL){
			//printf("1: %d ",sc2->sneighbor->ID);getchar();
			if(   (fabs(g->node[i].node_dD_props->poz[0]-g->node[sc2->sneighbor].node_dD_props->poz[0])<=0.5)
				  &&(fabs(g->node[i].node_dD_props->poz[1]-g->node[sc2->sneighbor].node_dD_props->poz[1])<=0.5)
				  //&&(labs(g->node[i].k-g->node[sc2->sneighbor].k)==1)
				  )
			{EpsDrawLine(net,500*g->node[i].node_dD_props->poz[0],500*g->node[i].node_dD_props->poz[1],
						 500*g->node[sc2->sneighbor].node_dD_props->poz[0],500*g->node[sc2->sneighbor].node_dD_props->poz[1]);
			}
			//else	{EpsSetRgb(net,0.8,0,0);EpsSetLinewidth(net,0.25);
			//		 EpsDrawLine(net,500*g->node[i].node_dD_props->poz[0],500*g->node[i].node_dD_props->poz[1],
			//					500*g->node[sc2->sneighbor].node_dD_props->poz[0],500*g->node[sc2->sneighbor].node_dD_props->poz[1]);
			//}
			sc2=sc2->next(sc2);//getchar();
		}
	}

	erdo=measure_gradnet_clusters(e);
	transfer_clusternumber(e,g);//printf("clust transfer Megvan\n");getchar();
	delete e; e=NULL;

	for(i=1;i<=g->N;i++){
	//printf("k=%lld (%lld-%lld)\tcolor nr %d\n",g->node[i].k,k_min,k_max,(unsigned int)(12*(g->node[i].k-k_min)/(k_max-k_min)));
		if(g->node[i].node_props->cluster_ID<=COL_ORDER_MAX)
			EpsSetColor_PAJEK(net,szinnev[g->node[i].node_props->cluster_ID]);
		EpsFillCircle(net,500*g->node[i].node_dD_props->poz[0],
				  500*g->node[i].node_dD_props->poz[1],
				  500*R);
	
	}
	EpsClose(net);
//	printf("End.\n");//getchar();
}
*/

void export_pieces_of_network(sgraph *g,char pa[]){
	sgraph *g_part;
	char fn[100];
	graph_clusters *erdo;
	FILE *f;
	unsigned long long i,j;
	
	sprintf(fn,"%s/Membership_%s.inf",pa,g->gname);
	f=fopen(fn,"w");
	erdo=measure_clusters(g,fn);
	for(i=1;i<=g->N;i++)
		fprintf(f,"%lld\t%lld\n",i,g->node[i].node_props->cluster_ID);
	fclose(f);
	for(j=1;j<=erdo->cluster_NR;j++){
		if(erdo->cl_size[j]>1){
			printf("j=%lld\tcl_size=%lld\n",j,erdo->cl_size[j]);//getchar();
			g_part=get_cluster_component(g,erdo,j);
			write_network_draw_attr(g_part,pa);
			//if(j!=erdo->max_ID) Static_3D_net_writeout(g_part,1,pa);
			delete g_part;g_part=NULL;
		}
	}
	fclose(f);
	delete g;g=NULL;
	delete erdo;erdo=NULL;
}
