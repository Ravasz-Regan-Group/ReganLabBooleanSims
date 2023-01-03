void set_RGB_hex(double r,double g,double b,char rgbn[]){
		strcpy(rgbn,"");
		if(r<16) sprintf(rgbn,"#0%x",(int)r); else sprintf(rgbn,"#%2x",(int)r);
		if(g<16) sprintf(rgbn,"%s0%x",rgbn,(int)g); else sprintf(rgbn,"%s%2x",rgbn,(int)g);
		if(b<16) sprintf(rgbn,"%s0%x",rgbn,(int)b);else sprintf(rgbn,"%s%2x",rgbn,(int)b);
	}

class RGB{
	public:
		double r,g,b;
		char RGB_hex[12];
	RGB(){r=0,g=0,b=0;sprintf(RGB_hex,"#000000");}
	RGB(double rr,double gg,double bb){
		if(rr>255) rr=255;if(gg>255) gg=255;if(bb>255) bb=255;
		r=rr;g=gg;b=bb;strcpy(RGB_hex,"");
		if(r<16) sprintf(RGB_hex,"#0%x",(int)r); else sprintf(RGB_hex,"#%2x",(int)r);
		if(g<16) sprintf(RGB_hex,"%s0%x",RGB_hex,(int)g); else sprintf(RGB_hex,"%s%2x",RGB_hex,(int)g);
		if(b<16) sprintf(RGB_hex,"%s0%x",RGB_hex,(int)b);else sprintf(RGB_hex,"%s%2x",RGB_hex,(int)b);
	}
	void setRGB(double rr,double gg,double bb){
		if(rr>255) rr=255;if(gg>255) gg=255;if(bb>255) bb=255;
		r=rr;g=gg;b=bb;strcpy(RGB_hex,"");
		if(r<16) sprintf(RGB_hex,"#0%x",(int)r); else sprintf(RGB_hex,"#%2x",(int)r);
		if(g<16) sprintf(RGB_hex,"%s0%x",RGB_hex,(int)g); else sprintf(RGB_hex,"%s%2x",RGB_hex,(int)g);
		if(b<16) sprintf(RGB_hex,"%s0%x",RGB_hex,(int)b);else sprintf(RGB_hex,"%s%2x",RGB_hex,(int)b);
	}
};

class link_attr {
	public:
		double lw,l_prop;
		RGB link_rgb;		
		char color[20];
		char *link_label;
        int ID, distinc_link_labels;
		
        link_attr(double weight){
			lw=weight;l_prop=1;
			link_label=new char[2];
            strcpy(link_label,"");link_label[1]=0;
			sprintf(color,"Black");
            distinc_link_labels=0;
            ID=0;
		}
		link_attr(double weight,const char lab[]){
			lw=weight;l_prop=1;
			link_label=new char[strlen(lab)+2];
            sprintf(link_label,"%s",lab);
			link_label[strlen(lab)+1]=0;
            sprintf(color,"Black");
            distinc_link_labels=0;
            ID=0;
		}
		void set_link_weigth(double w){lw=w;}
    
        void add_to_link_weight(double w){lw+=w;}
		void set_link_color(const char co[]){sprintf(color,"%s",co);}
        void set_link_label(const char la[]){
            if (link_label!=NULL) delete[] link_label;
            link_label=new char[strlen(la)+2];
            sprintf(link_label,"%s",la);distinc_link_labels=1;
            link_label[strlen(la)+1]=0;
        }
        void append_link_label(char la[]){
            if(!strstr(link_label,la)){
                char *oldlabel;
                oldlabel=link_label;
                link_label=new char[strlen(oldlabel)+strlen(la)+4];
                sprintf(link_label,"%s; %s",oldlabel,la);
                distinc_link_labels++;
                link_label[strlen(oldlabel)+strlen(la)+3]=0;
                delete[] oldlabel;oldlabel=NULL;
            }
        }
        void set_link_rgb(unsigned short int r,unsigned short int g,unsigned short int b)
			{link_rgb.setRGB(r,g,b);}
};

class node_attr {
	public:
		unsigned long long ID;
		unsigned long long cluster_ID,basin;
		double nw;
		node_attr(double w){
			ID=0;cluster_ID=0;nw=w;basin=0;
		}
		node_attr(unsigned long long i){
			ID=i;cluster_ID=0;nw=0;basin=0;
		}
		node_attr(unsigned long long i,double w){
			ID=i;cluster_ID=0;nw=w;basin=0;
		}
		node_attr(unsigned long long i,double w,unsigned long long cl){
			ID=i;cluster_ID=cl;nw=w;basin=0;
		}
		void set_node_w(double w){nw=w;}
		void set_node_basin(unsigned long long w){basin=w;}
		void set_node_cl(unsigned long long cl){cluster_ID=cl;}
};

class node_draw_attr {
	public:
		node_attr *basic;
		RGB node_rgb;
		char color[20];
		char node_label[501];
		char node_shape[31];
		double poz[3];
		node_draw_attr(unsigned long long i){
			basic=new node_attr(i);
			sprintf(node_label,"%lld",i);
			sprintf(color,"Red");
			poz[0]=0;poz[1]=0;poz[2]=0;
		}
        node_draw_attr(unsigned long long i,char nname[]){
            basic=new node_attr(i);
            sprintf(node_label,"%s",nname);
            sprintf(color,"Red");
            poz[0]=0;poz[1]=0;poz[2]=0;
        }
        node_draw_attr(unsigned long long i,double w){
			basic=new node_attr(i,w);
			sprintf(node_label,"%lld",i);
			sprintf(color,"Red");
			poz[0]=0;poz[1]=0;poz[2]=0;
		}
		node_draw_attr(unsigned long long i,double w,unsigned long long cl){
			basic=new node_attr(i,w,cl);
			sprintf(node_label,"%lld",i);
			sprintf(color,"Red");
			poz[0]=0;poz[1]=0;poz[2]=0;
		}
		void set_node_weigth(double w)
			{basic->set_node_w(w);}
		void set_node_cluster(unsigned long long cl){basic->set_node_cl(cl);}
		void set_node_color(char co[]){sprintf(color,"%s",co);}
		void set_node_shape(char sh[]){strcpy(node_shape,sh);}
		void set_node_label(char la[]){strcpy(node_label,la);}
		void set_node_rgb(unsigned short int r,unsigned short int g,unsigned short int b)
			{node_rgb.setRGB(r,g,b);}
		void set_node_poz(double x, double y,double z){poz[0]=x,poz[1]=y,poz[2]=z;}
		//~node_draw_attr(void){delete basic;basic=NULL;}
};

class node_dD_attr {
	public:
		node_attr *basic;
		double *poz;
		node_dD_attr(unsigned long long i){
			basic=new node_attr(i);
			poz=NULL;
		}
		node_dD_attr(unsigned long long i,unsigned int d,double co[]){
			unsigned int j;
			basic=new node_attr(i);
			poz=new double[d];
			for(j=0;j<d;j++) poz[j]=co[j];	
		}
		node_dD_attr(node_attr *np,unsigned int d,double co[]){
			unsigned int j;
			basic=np;
			poz=new double[d];
			for(j=0;j<d;j++) poz[j]=co[j];	
		}
		node_dD_attr(unsigned long long i,double w){
			basic=new node_attr(i,w);
			poz=NULL;
		}
		node_dD_attr(unsigned long long i,double w,unsigned long long cl){
			basic=new node_attr(i,w,cl);
			poz=NULL;
		}
		void set_node_weigth(double w)
			{basic->set_node_w(w);}
		void set_node_cluster(unsigned long long cl){basic->set_node_cl(cl);}
		void set_node_poz(unsigned int d, double co[]){
			unsigned int i;
			poz=new double[d];
			for(i=0;i<d;i++) poz[i]=co[i];	
		}
		~node_dD_attr(void){
			delete basic;basic=NULL;
			delete poz; poz=NULL;
		}
};

class snode {
	public:
		unsigned long long k,k_in,k_out;
		node_attr *node_props;
		node_draw_attr *node_draw_props;
		node_dD_attr *node_dD_props;
		unsigned long long ID;

		class slinklist {
			public:
				class slink {
					public:
							unsigned long long sneighbor;
							link_attr *link_props;
					private:
						slink *next_link;
					public:
						slink() {sneighbor=0;link_props=NULL;}
						slink *next(slink *pt){return(pt->next_link);}
						~slink(void){delete link_props;link_props=NULL;}
					friend class slinklist;			
					friend class snode;			
				};
			public:
				slink *first_link_pt;
				
				slinklist(void) {first_link_pt=NULL;} //initialize slinklist;
				~slinklist(void) {
					slink *todel;
					todel=first_link_pt;
					while(todel!=NULL){
						first_link_pt=first_link_pt->next_link;
						delete todel;
						todel=first_link_pt;
					}
				} //initialize slinklist;
				unsigned short int get_dirlink(unsigned long long nodeID){
					slink *scroll;
					scroll=first_link_pt;
					while(scroll!=NULL){
					   if(scroll->sneighbor==nodeID) {return(1);}
					   scroll=scroll->next_link;
					}
					return(0);
					}
				
				double get_dirlink_weight(unsigned long long nodeID){
					slink *scroll;
					scroll=first_link_pt;
					while(scroll!=NULL){
						if(scroll->sneighbor==nodeID) {
							if(scroll->link_props!=NULL)  return(scroll->link_props->lw);
							else return(1.);
						}
						scroll=scroll->next_link;
					}
					return(0);
				}
				unsigned int delete_dirlink(unsigned long long nodeID){
					slink *scroll,*before;
					
					if(first_link_pt==NULL) return(0);
					scroll=first_link_pt;
					if(first_link_pt->sneighbor==nodeID){
						//printf("first\t");//getchar();
						first_link_pt=scroll->next_link;
						if(scroll->link_props!=NULL) 
								{delete scroll->link_props;scroll->link_props=NULL;}
						delete scroll;scroll=NULL;
						return(1);
					}
					scroll=first_link_pt->next_link;
					before=first_link_pt;
					while(scroll!=NULL){
						if(scroll->sneighbor==nodeID){
							before->next_link=scroll->next_link;
							if(scroll->link_props!=NULL) 
								{delete scroll->link_props;scroll->link_props=NULL;}
							delete scroll;scroll=NULL;
							return(1);
						}
					   before=scroll;	
					   scroll=scroll->next_link;
					}
					return(0);
				}

				link_attr *dirlink_prop(unsigned long long nodeID){
					slink *scroll;
					scroll=first_link_pt;
					while(scroll!=NULL){
					   if(scroll->sneighbor==nodeID) {
						//if(scroll->link_props!=NULL) printf("\t\t--->Props pointer exists: lw=%lf\n",scroll->link_props->lw);
						//else printf("\t\t---> Link but no pointer\n");	
						return(scroll->link_props);
						}
					   scroll=scroll->next_link;
					}
					//printf("\t\t---> Did not find dirlink!\n");
					return(NULL);
				}

				void set_dirlink_prop(unsigned long long nodeID,link_attr *pro){
					slink *scroll;
					scroll=first_link_pt;
					while(scroll!=NULL){
					   if(scroll->sneighbor==nodeID) {
						scroll->link_props=pro;
						return;
						}
					   scroll=scroll->next_link;
					}
				}
				void delete_all_linkprops(){
					slink *scroll;
					scroll=first_link_pt;
					while(scroll!=NULL){
					   if(scroll->link_props!=NULL){delete scroll->link_props;scroll->link_props=NULL;}
					   scroll=scroll->next_link;
					}
				}							
				void list_node_linklist(){
					slink *scroll;
					scroll=first_link_pt;
					printf(" links:{");
					while(scroll!=NULL){
					   if(scroll->next_link==NULL) printf("%lld",scroll->sneighbor);
					   else printf("%lld,",scroll->sneighbor);
					   scroll=scroll->next_link;
					}
					printf("}\n");
				}
//			friend class snode;
		};
	public: 
		struct slinklist Li;
		
		snode(void){ID=0;k=0;k_in=0;k_out=0;node_props=NULL;node_dD_props=NULL;node_draw_props=NULL;}
		
		unsigned short int add_node_dirlink(unsigned long long id2){
			slinklist::slink *new_l;
			if(Li.get_dirlink(id2)) return(0);
			new_l=new slinklist::slink;
			new_l->sneighbor=id2;
			new_l->next_link=Li.first_link_pt;
			Li.first_link_pt=new_l;	
			return(1);
		}
		void set_link_w(unsigned long long id2,double w){
			slinklist::slink *scroll;
			
			scroll=Li.first_link_pt;
			while((scroll!=NULL)&&(scroll->sneighbor!=id2)){
				scroll=scroll->next_link;
				}
			if(scroll==NULL) {printf("%lld is not a neighbor!\n",id2);getchar();
							  return;
			}	
			if(scroll->link_props==NULL)
				scroll->link_props=new link_attr(w);
			else scroll->link_props->set_link_weigth(w);
		}
	
	void add_to_link_weight(unsigned long long id2,double w){
		slinklist::slink *scroll;
		
		scroll=Li.first_link_pt;
		while((scroll!=NULL)&&(scroll->sneighbor!=id2)){
			scroll=scroll->next_link;
		}
		if(scroll==NULL) {printf("%lld is not a neighbor!\n",id2);getchar();
			return;
		}	
		if(scroll->link_props==NULL)
			scroll->link_props=new link_attr(w);
		else scroll->link_props->add_to_link_weight(w);
	}
	
	void set_link_w_ID_disc_prop(unsigned long long id2,double w,int idnow,int lp){
		slinklist::slink *scroll;
		
		scroll=Li.first_link_pt;
		while((scroll!=NULL)&&(scroll->sneighbor!=id2)){
			scroll=scroll->next_link;
		}
		if(scroll==NULL) {printf("%lld is not a neighbor!\n",id2);getchar();
			return;
		}	
		if(scroll->link_props==NULL)
			scroll->link_props=new link_attr(w);
		else scroll->link_props->set_link_weigth(w);
		scroll->link_props->distinc_link_labels=lp;
        scroll->link_props->ID=idnow;
	}
    
    void set_link_label(unsigned long long id2,double w,char la[]){
		slinklist::slink *scroll;
		
		scroll=Li.first_link_pt;
		while((scroll!=NULL)&&(scroll->sneighbor!=id2)){
			scroll=scroll->next_link;
		}
		if(scroll==NULL) {printf("%lld is not a neighbor!\n",id2);getchar();
			return;
		}
		if(scroll->link_props==NULL){
			scroll->link_props=new link_attr(w);
            scroll->link_props->set_link_label(la);
        }
		else {
            scroll->link_props->set_link_weigth(w);
            scroll->link_props->set_link_label(la);
        }
	}
    
    void append_link_label(unsigned long long id2,double w,char la[]){
		slinklist::slink *scroll;
		
		scroll=Li.first_link_pt;
		while((scroll!=NULL)&&(scroll->sneighbor!=id2)){
			scroll=scroll->next_link;
		}
		if(scroll==NULL) {printf("%lld is not a neighbor!\n",id2);getchar();
			return;
		}
		if(scroll->link_props==NULL){
			scroll->link_props=new link_attr(w);
            scroll->link_props->set_link_label(la);
        }
		else {
            scroll->link_props->append_link_label(la);
        }
	}

    
	
		double get_link_w(unsigned long long id2){
			slinklist::slink *scroll;

			scroll=Li.first_link_pt;
			while((scroll!=NULL)&&(scroll->sneighbor!=id2)){
				scroll=scroll->next_link;
				}
			if(scroll==NULL) {printf("%lld is not a neighbor!\n",id2);getchar();
							  return(-1);
			}	
			if(scroll->link_props==NULL) return(0);
			else return(scroll->link_props->lw);
		}
    
    double get_link_prop(unsigned long long id2){
        slinklist::slink *scroll;
        
        scroll=Li.first_link_pt;
        while((scroll!=NULL)&&(scroll->sneighbor!=id2)){
            scroll=scroll->next_link;
        }
        if(scroll==NULL) {printf("%lld is not a neighbor!\n",id2);getchar();
            return(-1);
        }
        if(scroll->link_props==NULL) return(0);
        else return(scroll->link_props->l_prop);
    }

		
		slinklist::slink *get_directed_link(unsigned long long id2){
			slinklist::slink *scroll;
			
			scroll=Li.first_link_pt;
			while((scroll!=NULL)&&(scroll->sneighbor!=id2)){
				scroll=scroll->next_link;
			}
			return(scroll);
		}
		
		~snode(void){
			delete node_props;node_props=NULL;
			delete node_draw_props;node_draw_props=NULL;
		}
	friend class sgraph;
};


class sgraph {
	friend class snode;
	friend class slinklist;
	public:
		long long unsigned int N;
		char gname[300];
		unsigned int D;
	public:
		snode *node;
		sgraph(unsigned long long n) {
			N=n;
			node=new snode[N+1];
			D=0;
		} // initialize graph
		sgraph(unsigned long long n,const char nam[]) {
			N=n;
			node=new snode[N+1];
			sprintf(gname,"%s",nam);
			D=0;
		} // initialize graph
        sgraph(const char nam[]) {
            N=0;
            node=NULL;
            sprintf(gname,"%s",nam);
            D=0;
        } // initialize graph
        
        ~sgraph(void){
			delete[] node;node=NULL;
		}
    
        unsigned long int add_node(char nname[]){
            snode *newp;
            unsigned long int i,j;
            
            j=check_for_node_with_name(nname);
            if(j==0){
                if(node==NULL){
                    N=1;
                    node=new snode[2];
                    node[1].node_draw_props=new node_draw_attr(1,nname);
                    node[1].node_props=new node_attr((double)1);
                    return(1);
                }
                else{
                    newp=node;node=NULL;
                    N++;
                    node=new snode[N+1];
                    for(i=1;i<N;i++)	node[i]=newp[i];
                    node[N].node_draw_props=new node_draw_attr(N,nname);
                    node[N].node_props=new node_attr((double)1);
                    //printf("about to delete element list. N=%ld\n",N);getchar();
                   // delete[] newp;
                    newp=NULL;
                    //printf("\t\tdone\n");getchar();
                    return(N);
                }
            }
            else node[j].node_props->nw++;
            return(j);
        }
        unsigned long int check_for_node_with_name(char nname[]){
            for (unsigned long int i=1; i<=N; i++) {
                if (node[i].node_draw_props!=NULL)
                    if (!strcmp(nname,node[i].node_draw_props->node_label)) return(i);
            }
            return(0);
        }
    
        void set_name(char nam[]){sprintf(gname,"%s",nam);}
		void add_node_props(unsigned long long i,node_attr *prop) {
			 if((i<1)||(i>N)||prop==NULL) return;
			 node[i].node_props=prop;
			}
	
		void add_node_draw_props(unsigned long long i,node_draw_attr *prop) {
			 if((i<1)||(i>N)||prop==NULL) return;
			 node[i].node_draw_props=prop;
			 node[i].node_props=node[i].node_draw_props->basic;
			}
    
		unsigned long long find_ID(unsigned long long id1){
			unsigned long long i;
			for(i=1;i<=N;i++){
				if(node[i].node_props!=NULL){
					if(node[i].node_props->ID==id1) return(i);
				}
				else return(0);
			}
            return 0;
		}
    
        unsigned long long index_by_nodeID(unsigned long long id1){
            unsigned long long i;
            for(i=1;i<=N;i++){
                if(node[i].ID==id1) return(i);
            }
            return 0;
        }
		
	    unsigned short int add_link(unsigned long long id1,unsigned long long id2){
			if(!node[id1].add_node_dirlink(id2)) return(0);
			node[id1].k_out++;
		    node[id2].k_in++;
		    if(!node[id2].Li.get_dirlink(id1)) 
				{node[id1].k++;node[id2].k++;}
		    if(id1==id2) node[id1].k++;	
			//printf("%lld: ",id1);
			//node[id1].Li.list_node_linklist(); 
			//getchar();
			return(1);
		}
		unsigned short int add_link_undir_check(unsigned long long id1,unsigned long long id2,link_attr *prop){
			if(linked(id1,id2)) return(0);
			if(!node[id1].add_node_dirlink(id2)) return(0);
			node[id1].k_out++;
		    node[id2].k_in++;
		    if(!node[id2].Li.get_dirlink(id1)) 
				{node[id1].k++;node[id2].k++;}
			if(prop!=NULL) node[id1].Li.set_dirlink_prop(id2,prop);	
		    if(id1==id2) node[id1].k++;	
			//printf("%lld: ",id1);
			//node[id1].Li.list_node_linklist(); 
			//getchar();
			return(1);
		}
		unsigned short int add_link_undir_check(unsigned long long id1,unsigned long long id2){
			if(linked(id1,id2)) return(0);
			if(!node[id1].add_node_dirlink(id2)) return(0);
			node[id1].k_out++;
		    node[id2].k_in++;
		    if(!node[id2].Li.get_dirlink(id1)) 
				{node[id1].k++;node[id2].k++;}
		    if(id1==id2) node[id1].k++;	
			//printf("%lld: ",id1);
			//node[id1].Li.list_node_linklist(); 
			//getchar();
			return(1);
		}

		unsigned short int add_link(unsigned long long id1,unsigned long long id2,link_attr *prop){
			if(!node[id1].add_node_dirlink(id2)) return(0);
			node[id1].k_out++;
		    node[id2].k_in++;
		    if(!node[id2].Li.get_dirlink(id1)) 
				{node[id1].k++;node[id2].k++;}
			if(prop!=NULL) node[id1].Li.set_dirlink_prop(id2,prop);	
		    if(id1==id2) node[id1].k++;	
			//printf("%lld: ",id1);
			//node[id1].Li.list_node_linklist(); 
			//getchar();
			return(1);
		}
		
		unsigned short int add_link(unsigned long long id1,unsigned long long id2,double lw){
			if(!node[id1].add_node_dirlink(id2)) return(0);
			node[id1].k_out++;
			node[id2].k_in++;
			if(!node[id2].Li.get_dirlink(id1)) 
			{node[id1].k++;node[id2].k++;}
			if(id1==id2) node[id1].k++;	
			//printf("%lld: ",id1);
			//node[id1].Li.list_node_linklist(); 
			//getchar();
			node[id1].set_link_w(id2,lw);
			return(1);
		}
		unsigned short int add_link_and_count(unsigned long long id1,unsigned long long id2){
			
			if(!node[id1].add_node_dirlink(id2)) {
				node[id1].add_to_link_weight(id2,1.);
				return(0);
			}
			node[id1].k_out++;
			node[id2].k_in++;
			if(!node[id2].Li.get_dirlink(id1)) 
			{node[id1].k++;node[id2].k++;}
			if(id1==id2) node[id1].k++;	
			//printf("%lld: ",id1);
			//node[id1].Li.list_node_linklist(); 
			//getchar();
			node[id1].set_link_w(id2,1.);
			return(1);
		}
	/*
	unsigned short int add_link(unsigned long long id1,unsigned long long id2,double lw,double lp){
		if(!node[id1].add_node_dirlink(id2)) return(0);
		node[id1].k_out++;
		node[id2].k_in++;
		if(!node[id2].Li.get_dirlink(id1)) 
		{node[id1].k++;node[id2].k++;}
		if(id1==id2) node[id1].k++;	
		//printf("%lld: ",id1);
		//node[id1].Li.list_node_linklist(); 
		//getchar();
		node[id1].set_link_w(id2,lw);
		return(1);
	}
    */
    
    unsigned short int add_link(unsigned long long id1,unsigned long long id2,double lw,int idnow,int lpnow){
        if(!node[id1].add_node_dirlink(id2)) return(0);
        node[id1].k_out++;
        node[id2].k_in++;
        if(!node[id2].Li.get_dirlink(id1))
        {node[id1].k++;node[id2].k++;}
        if(id1==id2) node[id1].k++;
        //printf("%lld: ",id1);
        //node[id1].Li.list_node_linklist();
        //getchar();
        node[id1].set_link_w_ID_disc_prop(id2,lw,idnow,lpnow);
        return(1);
    }
    
	unsigned short int add_link_with_label(unsigned long long id1,unsigned long long id2,double w,char la[]){
		if(!node[id1].add_node_dirlink(id2)) return(0);
		node[id1].k_out++;
		node[id2].k_in++;
		if(!node[id2].Li.get_dirlink(id1)) {node[id1].k++;node[id2].k++;}
		if(id1==id2) node[id1].k++;
		//printf("%lld: ",id1);
		//node[id1].Li.list_node_linklist();
		//getchar();
		node[id1].set_link_label(id2,w,la);
        
		return(1);
	}
    unsigned short int add_link_with_label_append(unsigned long long id1,unsigned long long id2,double w,char la[]){
		if(!node[id1].add_node_dirlink(id2)) {node[id1].append_link_label(id2,w,la); return(0);}
		node[id1].k_out++;
		node[id2].k_in++;
		if(!node[id2].Li.get_dirlink(id1)) {node[id1].k++;node[id2].k++;}
		if(id1==id2) node[id1].k++;
		//printf("%lld: ",id1);
		//node[id1].Li.list_node_linklist();
		//getchar();
		node[id1].append_link_label(id2,w,la);
        
		return(1);
	}
    

  
    
    
		void add_link_prop(unsigned long long id1,unsigned long long id2,link_attr *prop){
			if(prop!=NULL) node[id1].Li.set_dirlink_prop(id2,prop);	
		}
		
		link_attr *get_undir_link(unsigned long long id1,unsigned long long id2){
			link_attr *l;
			l=node[id1].Li.dirlink_prop(id2);
			if(l!=NULL) return(l);
			else return(node[id2].Li.dirlink_prop(id1));
		}
    
        link_attr *get_dir_link(unsigned long long id1,unsigned long long id2){
            link_attr *l;
            l=node[id1].Li.dirlink_prop(id2);
            if(l!=NULL) return(l);
            else return(NULL);
        }

		void delete_link(unsigned long long id1,unsigned long long id2){
			//node[id1].Li.list_node_linklist(); 
			unsigned int ok;
			ok=node[id1].Li.delete_dirlink(id2);
			if(ok==1) {node[id2].k_in--;node[id1].k_out--;}
			//node[id1].Li.list_node_linklist(); 
			//getchar();
		}
		void delete_link_undir(unsigned long long id1,unsigned long long id2){
			unsigned int ok;
			//printf("%lld: ",id1);node[id1].Li.list_node_linklist(); 
			//printf("%lld: ",id2);node[id2].Li.list_node_linklist(); 
			ok=node[id1].Li.delete_dirlink(id2);
			if(ok==1) {node[id2].k_in--;node[id1].k_out--;}
			ok=node[id2].Li.delete_dirlink(id1);
			if(ok==1) {node[id1].k_in--;node[id2].k_out--;}
			//printf("%lld: ",id1);node[id1].Li.list_node_linklist(); 
			//printf("%lld: ",id2);node[id2].Li.list_node_linklist(); 
			//getchar();
		}
		unsigned short int linked(unsigned long long id1, unsigned long long id2){
			if(node[id1].Li.get_dirlink(id2)!=0) return(1);
			if(node[id2].Li.get_dirlink(id1)!=0) return(1);
			return(0);		
		}
		link_attr *link_prop(unsigned long long id1, unsigned long long id2){
			link_attr *p;
			p=node[id1].Li.dirlink_prop(id2);
			if(p!=NULL) return(p);
			p=node[id2].Li.dirlink_prop(id1);
			if(p!=NULL) return(p);
			return(NULL);		
		}
		void get_neighbors(unsigned long long nd,unsigned long long neigh[]){
			unsigned long long i,j;
			
			if(node[nd].k==0) return;
			j=0;
			for(i=1;i<=N;i++){
				if((node[nd].Li.get_dirlink(i))||(node[i].Li.get_dirlink(nd))) 
					{neigh[j]=i;j++;}
			}
		}
	
		longlongint_ARRAY *get_neighbors_bothways(unsigned long long nd){
			unsigned long long i,j;
			longlongint_ARRAY *nbrs;
			
			nbrs=new longlongint_ARRAY();
			if(node[nd].k==0) return(nbrs);
			j=0;
			for(i=1;i<=N;i++){
				if((node[nd].Li.get_dirlink(i))||(node[i].Li.get_dirlink(nd))) {
					nbrs->add_element(i);
				}
			}
			return(nbrs);
		}
		void get_directed_neighbors_from_nodes_list(unsigned long long nd,unsigned long long neigh[]){
			unsigned long long j;
			snode::slinklist::slink *sc2;
			
			j=0;
			sc2=node[nd].Li.first_link_pt;
			while(sc2!=NULL){
				j++;
				neigh[j]=sc2->sneighbor;
				sc2=sc2->next(sc2);
			}
			neigh[0]=j;
			return;
		}
	
	
		double get_knn_for(unsigned long long nd){
			unsigned long long i,j;
			double knnav=0;
			
			if(node[nd].k==0) return(0);
			j=0;knnav=0;
			for(i=1;i<=N;i++){
				if((node[nd].Li.get_dirlink(i))||(node[i].Li.get_dirlink(nd))) 
					{knnav+=node[i].k;j++;}
			}
			if(j!=node[nd].k) {printf("error in get_knn_for: j!=node[nd].k!\n");getchar();}
			return(knnav/(double)j);
		}

		double clust_coeff(unsigned long long nd){
			unsigned long long n=0,i,j,*neigh;
			
			if(node[nd].k<2) return(0.);
			neigh=new unsigned long long[node[nd].k];
			j=0;
			for(i=1;i<=N;i++){	if(linked(nd,i)) {neigh[j]=i;j++;}}
			
			//printf("%lld registered =%lld, found is %lld\n",nd->ID,nd->k,i); 
			for(i=1;i<node[nd].k;i++)
				for(j=0;j<i;j++) {if(linked(neigh[i],neigh[j])) n++;}
			delete[] neigh;
			neigh=NULL;
			//printf("\t\t n=%lld\tC=%lf\n",n,(2.*n)/(double)(nd->k*(nd->k-1.)));getchar();
			return((2*n)/(double)(node[nd].k*(node[nd].k-1.)));
		}
		unsigned long long get_stored_link_number(){
			unsigned long long i,Ln=0;
			//printf("started...N=%lld",N);getchar();
			for(i=1;i<=N;i++) Ln+=node[i].k_out;
			return(Ln);
		}
    
    
		double get_DE_sum(){
			unsigned long long i;
			snode::slinklist::slink *lp;
			double de=0;
			
			for(i=1;i<=N;i++) {
				lp=node[i].Li.first_link_pt;
			    while(lp!=NULL){
				     de+=fabs(node[i].node_props->nw - node[lp->sneighbor].node_props->nw);
					 // printf("%lld->%lf\t%lld->%lf\tde=%lf\n",
					  //  i,node[i].node_props->nw,lp->sneighbor,node[lp->sneighbor].node_props->nw,de);
					 //getchar();
					 lp=lp->next(lp);
				}
			}
			return(de);
		}
		double overlap(unsigned long long i,unsigned long long j){
			unsigned long long k,n=0,mk;
			if(i==j) return(0);
			if((node[i].k==0)||(node[j].k==0)){
				printf("node with degree 0!!! largest cluster not checked!\n");
				return(0);
			}
			for(k=1;k<=N;k++) 
				if((k!=i)&&(k!=j)&&(linked(i,k))&&(linked(j,k))) n++;
			mk=node[i].k;
			if(node[j].k<mk) mk=node[j].k;
			if(linked(i,j))  return((n+1)/(double)(mk));
			else return((n+1)/(double)(mk+1));
		}

		void kmin_kmax(unsigned long long *kmin,unsigned long long *kmax){
			unsigned long long i;
			*kmin=*kmax=node[1].k;
			for(i=2;i<=N;i++){
				if(node[i].k<*kmin) *kmin=node[i].k;
				if(node[i].k>*kmax) *kmax=node[i].k;
			}
		}

		void IN_kmin_kmax(unsigned long long *kmin,unsigned long long *kmax){
			unsigned long long i;
			*kmin=*kmax=node[1].k_in;
			for(i=2;i<=N;i++){
				if(node[i].k_in<*kmin) *kmin=node[i].k_in;
				if(node[i].k_in>*kmax) *kmax=node[i].k_in;
			}
		}
        
        
};

typedef sgraph *sgraph_array;

void PK_in(sgraph *g,unsigned long long pk_in[]){
			unsigned long long i;
			for(i=1;i<=g->N;i++)
				pk_in[g->node[i].k_in]++;
}

void PK_in_out_undir(sgraph *g,unsigned long long pk_in[],unsigned long long pk_out[],unsigned long long pk_undir[]){
			unsigned long long i;
			for(i=1;i<=g->N;i++){
				pk_out[g->node[i].k_out]++;
				pk_in[g->node[i].k_in]++;
				pk_undir[g->node[i].k]++;
			}
}

void measure_PK_in(sgraph *g,char path[]){
	unsigned long long *pk,i;
	FILE *ki;
	char nev[300];
	
	pk=new unsigned long long[g->N+2];
	for(i=0;i<=g->N+1;i++) pk[i]=0;
//	printf("N=%lld, starting pk\n",N);getchar();
	for(i=1;i<=g->N;i++) pk[g->node[i].k_in]++;
	
	sprintf(nev,"%s/PK_in__%s.dat",path,g->gname);
	ki=fopen(nev,"w");
	for(i=0;i<=g->N+1;i++) 
		if(pk[i]) 
			fprintf(ki,"%lld\t%lg\t%lld\n",i,pk[i]/(double)g->N,pk[i]);
	fclose(ki);
	delete[] pk;
	pk=NULL;
}

void measure_PK_in_out(sgraph *g,char path[]){
	unsigned long long *pk_in,*pk_out,i;
	FILE *ki;
	char nev[300];
	
	pk_in=new unsigned long long[g->N+2];
	pk_out=new unsigned long long[g->N+2];
	for(i=0;i<=g->N+1;i++) {
		pk_in[i]=0;pk_out[i]=0;
	}
	for(i=1;i<=g->N;i++){
		pk_out[g->node[i].k_out]++;
		pk_in[g->node[i].k_in]++;
	}
	
	sprintf(nev,"%s/PK_in__%s.dat",path,g->gname);
	ki=fopen(nev,"w");
	for(i=0;i<=g->N+1;i++) 
		if(pk_in[i]) 
			fprintf(ki,"%lld\t%lg\t%lld\n",i,pk_in[i]/(double)g->N,pk_in[i]);
	fclose(ki);
	sprintf(nev,"%s/PK_out__%s.dat",path,g->gname);
	ki=fopen(nev,"w");
	for(i=0;i<=g->N+1;i++) 
		if(pk_out[i]) 
			fprintf(ki,"%lld\t%lg\t%lld\n",i,pk_out[i]/(double)g->N,pk_out[i]);
	fclose(ki);
	delete[] pk_in;	pk_in=NULL;
	delete[] pk_out;pk_out=NULL;
}

unsigned long long PK(sgraph *g,unsigned long long pk_undir[]){
	unsigned long long i,maxk;
	double kav=0;
	maxk=0;if(g==NULL) return(0);
	for(i=1;i<=g->N;i++) pk_undir[i]=0;
	for(i=1;i<=g->N;i++){
		  if(maxk<g->node[i].k) maxk=g->node[i].k;
		  pk_undir[g->node[i].k]++;
		  kav+=g->node[i].k;
		  //if(i/1000==i/1000.) printf("%lld\tkmax=%lld\n",i,maxk);//getchar();
	}
	printf("k_av=%lg\tK_max=%lld\n",kav/(double)g->N,maxk);
	return(maxk);
}

void PK_CK(sgraph *g,unsigned long long pk_undir[],double ck[]){
	unsigned long long i,maxk;
	
	maxk=0;if(g==NULL) return;
	for(i=1;i<=g->N;i++){
		  if(maxk<g->node[i].k) maxk=g->node[i].k;
		  pk_undir[g->node[i].k]++;
		  ck[g->node[i].k]+=g->clust_coeff(i);
		  if(i/5000==i/5000.) printf("%lld\tkmax=%lld\n",i,maxk);//getchar();
	}
	//printf("I'm through, kmax=%lld\n",maxk);getchar();
	for(i=0;i<=maxk;i++) if(pk_undir[i]) {ck[i]/=(double)pk_undir[i];
			//printf("%lld ck=%lf\n",i,ck[i]);
			}//getchar();
}

void PK_CK_large(sgraph *g,unsigned long long pk_undir[],double ck[],char path[]){
	unsigned long long i,maxk,N_ok;
	char ciname[400];
	FILE *f;
	int vege;
	double C;
	
	maxk=0;if(g==NULL) return;
	for(i=0;i<=g->N;i++){
		pk_undir[i]=0;
		ck[i]=0;
	}
	sprintf(ciname,"%s/%s_Cnode.dat",path,g->gname);
	f=fopen(ciname,"r");
	i=0;
	if(f!=NULL){
		vege=1;
		while(vege!=EOF){
			vege=fscanf(f,"%lld%lg",&i,&C);
			if(vege!=EOF){
				ck[g->node[i].k]+=C;
				if(maxk<g->node[i].k) 
					maxk=g->node[i].k;
				pk_undir[g->node[i].k]++;
			}
		}
		fclose(f);
	}
	N_ok=i;
	printf("For %lld nodes out of %lld,  C is already measured!\n",N_ok,g->N);

	if(N_ok<g->N){
		f=fopen(ciname,"a");
		for(i=N_ok+1;i<=g->N;i++){
			if(maxk<g->node[i].k) maxk=g->node[i].k;
			pk_undir[g->node[i].k]++;
			C=g->clust_coeff(i);
			fprintf(f,"%lld\t%lg\n",i,C);
			fflush(f);
			ck[g->node[i].k]+=C;
			if(i/5000==i/5000.) printf("%lld\tkmax=%lld\n",i,maxk);//getchar();
		}
		fclose(f);
	}
	//printf("I'm through, kmax=%lld\n",maxk);getchar();
	for(i=0;i<=maxk;i++) if(pk_undir[i]) {ck[i]/=(double)pk_undir[i];
		//printf("%lld ck=%lf\n",i,ck[i]);
	}//getchar();
}

double measure_PK_undir(sgraph *g,char path[101]){
	unsigned long long N,*pk,i;
	double k_av=0.;
	char pkname[310];
	FILE *ki;
	if(g==NULL){printf("Grapg does not exist!\n");return(0);}
	N=g->N;
	pk=new unsigned long long[N+2];
	for(i=0;i<=N+1;i++) pk[i]=0;
//	printf("N=%lld, starting pk\n",N);getchar();
	PK(g,pk);
	sprintf(pkname,"%s/%s_PK.dat",path,g->gname);
	ki=fopen(pkname,"w");
		for(i=0;i<=N+1;i++) if(pk[i]) {
								fprintf(ki,"%lld\t%lf\t%lld\n",i,pk[i]/(double)N,pk[i]);
								k_av+=pk[i]*i;
							}
	k_av/=(double)N;
	fclose(ki);
	delete[] pk;
	pk=NULL;
	return(k_av);
}

void measure_PK_CK(sgraph *g,double *k_av,double *C_av,char path[]){
	unsigned long long N,*pk,i;
	double *ck;
	char pkname[310],ckname[310];
	FILE *ki;
	if(g==NULL){printf("Grapg does not exist!\n");return;}
	N=g->N;
	pk=new unsigned long long[N+2];
	ck=new double[N+2];
	for(i=0;i<=N+1;i++) {ck[i]=0;pk[i]=0;}
	PK_CK(g,pk,ck);
	*k_av=0;*C_av=0;
	sprintf(pkname,"%s/PK_%s.dat",path,g->gname);
	sprintf(ckname,"%s/CK_%s.dat",path,g->gname);
	ki=fopen(pkname,"w");
		for(i=0;i<=N+1;i++) if(pk[i]) {
					fprintf(ki,"%lld\t%lf\t%lld\n",i,pk[i]/(double)N,pk[i]);
					*k_av+=pk[i]*i;
					}
	fclose(ki);
	ki=fopen(ckname,"w");
		for(i=0;i<=N+1;i++) if(pk[i]) {
			fprintf(ki,"%lld\t%lf\n",i,ck[i]);
			*C_av+=ck[i]*pk[i];		
		}
	*k_av/=(double)N;
	*C_av/=(double)N;
	fclose(ki);
	delete[] pk;
	pk=NULL;
	delete[] ck;
	ck=NULL;
}

void measure_PK_CK(sgraph *g,char pkname[101],char ckname[101],double *k_av,double *C_av){
	unsigned long long N,*pk,i;
	double *ck;
	FILE *ki;
	if(g==NULL){printf("Grapg does not exist!\n");return;}
	N=g->N;
	printf("N=%lld\n",N);
	
	pk=new unsigned long long[N+2];
	ck=new double[N+2];
	for(i=0;i<=N+1;i++) {ck[i]=0;pk[i]=0;}
	PK_CK(g,pk,ck);
	*k_av=0;*C_av=0;
	ki=fopen(pkname,"w");
		for(i=0;i<=N+1;i++) if(pk[i]) {
					fprintf(ki,"%lld\t%lf\t%lld\n",i,pk[i]/(double)N,pk[i]);
					*k_av+=pk[i]*i;
					}
	fclose(ki);
	ki=fopen(ckname,"w");
		for(i=0;i<=N+1;i++) if(pk[i]) {
			fprintf(ki,"%lld\t%lf\n",i,ck[i]);
			*C_av+=ck[i]*pk[i];		
		}
	*k_av/=(double)N;
	*C_av/=(double)N;
	fclose(ki);
	delete[] pk;
	pk=NULL;
	delete[] ck;
	ck=NULL;
}

void measure_CK_large(sgraph *g,double *k_av,double *C_av,char path[]){
	unsigned long long N,*pk,i;
	double *ck;
	char ckname[310];
	FILE *ki;
	if(g==NULL){printf("Grapg does not exist!\n");return;}
	N=g->N;
	pk=new unsigned long long[N+2];
	ck=new double[N+2];
	for(i=0;i<=N+1;i++) {ck[i]=0;pk[i]=0;}
	PK_CK_large(g,pk,ck,path);
	*k_av=0;*C_av=0;
	sprintf(ckname,"%s/CK_%s.dat",path,g->gname);
	ki=fopen(ckname,"w");
	for(i=0;i<=N+1;i++) if(pk[i]) {
		fprintf(ki,"%lld\t%lf\n",i,ck[i]);
		*C_av+=ck[i]*pk[i];		
	}
		*k_av/=(double)N;
	*C_av/=(double)N;
	fclose(ki);
	delete[] pk;
	pk=NULL;
	delete[] ck;
	ck=NULL;
}


void measure_kav_Cav(sgraph *g,double *k_av,double *C_av){
	unsigned long long N,*pk,i;
	double *ck;
	
	if(g==NULL){printf("Grapg does not exist!\n");return;}
	N=g->N;
	pk=new unsigned long long[N+2];
	ck=new double[N+2];
	for(i=0;i<=N+1;i++) {ck[i]=0;pk[i]=0;}
	PK_CK(g,pk,ck);
	*k_av=0;*C_av=0;
	for(i=0;i<=N+1;i++) if(pk[i]) *k_av+=pk[i]*i;
	for(i=0;i<=N+1;i++) if(pk[i]) *C_av+=ck[i]*pk[i];		
	*k_av/=(double)N;
	*C_av/=(double)N;
	delete[] pk;
	pk=NULL;
	delete[] ck;
	ck=NULL;
}

double Newman_assortativity(sgraph *g){ 
	unsigned long int i;
	double M=0,jk=0,j_k=0,j2_k2=0;
	snode::slinklist::slink *sc2;
	
	for(i=2;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;
		while(sc2!=NULL){
			M++;
			jk+=g->node[i].k * g->node[sc2->sneighbor].k;
			j_k+=g->node[i].k + g->node[sc2->sneighbor].k;
			j2_k2+=g->node[i].k * g->node[i].k+g->node[sc2->sneighbor].k * g->node[sc2->sneighbor].k;
			sc2=sc2->next(sc2);
		}
	}	
	return((M*jk-0.25*j_k*j_k)/(M*0.5*j2_k2-0.25*j_k*j_k));
}

void measure_Knn(sgraph *g,char path[]){
	unsigned long long N,*pk,i;
	double *k_nn,*sz_knn,knn0;
	char kname[310];
	FILE *ki;

	if(g==NULL){printf("Grapg does not exist!\n");return;}
	N=g->N;
	k_nn=new double[N+2];
	pk=new unsigned long long[N+2];
	sz_knn=new double[N+2];
	
	for(i=0;i<=N+1;i++) 
		{k_nn[i]=0;sz_knn[i]=0;pk[i]=0;}
	PK(g,pk);
	
	for(i=1;i<=N;i++){
		knn0=g->get_knn_for(i);
		k_nn[g->node[i].k]+=knn0;
		sz_knn[g->node[i].k]+=knn0*knn0;
		if(i/100==i/100.) printf("\t%lld\n",i);
	}
	for(i=1;i<=N;i++)
		if(pk[i]) {
			k_nn[i]/=(double)pk[i];
			sz_knn[i]/=(double)pk[i];
			sz_knn[i]=sqrt(sz_knn[i]-k_nn[i]*k_nn[i]);
		}
	sprintf(kname,"%s/Knn_av_%s.dat",path,g->gname);
	
	ki=fopen(kname,"w");
	for(i=0;i<=N+1;i++) 
		if(k_nn[i]) {
		  fprintf(ki,"%lld\t%lg\t%lg\n",i,k_nn[i],sz_knn[i]);
	}
	fclose(ki);
	delete[] pk;pk=NULL;
	delete[] k_nn;k_nn=NULL;
	delete[] sz_knn;sz_knn=NULL;
}


void measure_degree_correlations_3D(sgraph *g,char path[]){
	FILE *ki;
	unsigned long long *nk,kmax,kmin,i,j,sum;
	longnaturalpointer *nk1k2;
	char name[300];
	snode::slinklist::slink *sc2;

	nk=new unsigned long long[g->N+1];
	kmax=PK(g,nk);sum=0;
	kmin=g->N;
	for(i=1;i<=g->N;i++){
		if(g->node[i].k<kmin) 
			kmin=g->node[i].k;
	}
	printf("kmax=%lld\t\tkmin=%lld\n",kmax,kmin);//getchar();
		
	nk1k2=new longnaturalpointer[kmax+1];
	for(i=0;i<=kmax;i++){
		nk1k2[i]=new unsigned long long[kmax+1];
		for(j=0;j<=kmax;j++)
			nk1k2[i][j]=0;
	}
	for(i=1;i<=g->N;i++){
		sc2=g->node[i].Li.first_link_pt;
		while(sc2!=NULL){
			nk1k2[g->node[i].k][g->node[sc2->sneighbor].k]++;
			nk1k2[g->node[sc2->sneighbor].k][g->node[i].k]++;
			sc2=sc2->next(sc2);
		}
		//if(i/1000==i/1000.) printf("%lld\n",i);
	}
	//getchar();
	sprintf(name,"%s/KKcorel__%s.gnu",path,g->gname);
	ki=fopen(name,"w");
	for(i=kmin;i<=kmax;i++){
		fprintf(ki,"\n");	
		for(j=kmin;j<=kmax;j++)
			{if(nk[i]*nk[j]!=0)
				fprintf(ki,"%lld\t%lld\t%lg\n",i,j,(g->N*nk1k2[i][j])/(double)((i+j)*nk[i]*nk[j])-1);
			 else fprintf(ki,"%lld\t%lld\t%lg\n",i,j,0.0);
			}
	}
	fclose(ki);		
	delete[] nk;nk=NULL;
	for(i=0;i<=kmax;i++) {delete[] nk1k2[i];nk1k2[i]=NULL;}
	delete[] nk1k2;nk1k2=NULL;
}

class graph_clusters {
	public:
		unsigned long long cluster_NR,cl_leq_2;
		unsigned long long max_size,max_ID;
		double fragments_av,frag_above_1_av;
		unsigned long long *cl_size,*size_dist;
		
		graph_clusters(sgraph *g,unsigned long long n){
			unsigned long long i;
			
			cluster_NR=n;max_size=0;max_ID=0;fragments_av=0.;frag_above_1_av=0.;cl_leq_2=0;
			//printf("Cl Nr=%lld\n",cluster_NR);getchar();	
			cl_size=new unsigned long long[cluster_NR+1];
			for(i=0;i<=cluster_NR;i++) cl_size[i]=0;
			for(i=1;i<=g->N;i++)
				cl_size[g->node[i].node_props->cluster_ID]++;
			for(i=1;i<=cluster_NR;i++) if(cl_size[i]>max_size) {max_size=cl_size[i];max_ID=i;}
			size_dist=new unsigned long long[max_size+1];
			for(i=0;i<=max_size;i++) size_dist[i]=0;
			for(i=1;i<=cluster_NR;i++) {
				size_dist[cl_size[i]]++;
				if(i!=max_ID) {
					fragments_av+=cl_size[i];
					if(cl_size[i]>1){ frag_above_1_av+=cl_size[i];cl_leq_2++;}
				}
			}
			if(cluster_NR>1) fragments_av/=(double)(cluster_NR-1);
			if(cl_leq_2>1) frag_above_1_av/=(double)(cl_leq_2);
		}
		graph_clusters(sgraph *g,unsigned long long n,const char outname[200]){
			unsigned long long i;
			FILE *ki;
			
			cluster_NR=n;max_size=0;max_ID=0;fragments_av=0.;
			//printf("Cl Nr=%lld\n",cluster_NR);getchar();	
			cl_size=new unsigned long long[cluster_NR+1];
			for(i=0;i<=cluster_NR;i++) cl_size[i]=0;
			for(i=1;i<=g->N;i++)
				cl_size[g->node[i].node_props->cluster_ID]++;
			for(i=1;i<=cluster_NR;i++) if(cl_size[i]>max_size) {max_size=cl_size[i];max_ID=i;}
			size_dist=new unsigned long long[max_size+1];
			for(i=0;i<=max_size;i++) size_dist[i]=0;
			for(i=1;i<=cluster_NR;i++) {
				size_dist[cl_size[i]]++;
				if(i!=max_ID) fragments_av+=cl_size[i];
			}
			if(cluster_NR>1) fragments_av/=(double)(cluster_NR-1);
			//printf("Starting cluster fragment file now:\n",cluster_NR);getchar();	
			ki=fopen(outname,"w");
			for(i=1;i<=max_size;i++) 
				if(size_dist[i]) 
					fprintf(ki,"%lld\t%lf\t%lld\n",i,size_dist[i]/(double)cluster_NR,size_dist[i]);	
			fclose(ki);
		}
		~graph_clusters(){
			delete[] cl_size;cl_size=NULL;
			delete[] size_dist;size_dist=NULL;
		}
	// define it again with distribution file input
};

void burn(sgraph *gg,unsigned long long st,unsigned long long IDa){
	unsigned long long *neigh,i;
	gg->node[st].node_props->cluster_ID=IDa;
	//printf("\nCaught %lld into %lld, starting on %lld neighbors\t",st,IDa,gg->node[st].k);
	//getchar();
	neigh=new unsigned long long[gg->node[st].k+1];
	gg->get_neighbors(st,neigh);
	for(i=0;i<gg->node[st].k;i++){
		// printf("\t\tst=%lld\ti=%lld\t %lld\n",st,i,neigh[i]);
 		 if((neigh[i]!=st)&&(gg->node[neigh[i]].node_props->cluster_ID==0)) 
			{burn(gg,neigh[i],IDa);}
	}
	//printf("\t%lld over\n",st);
	delete[] neigh;neigh=NULL;
}


graph_clusters *measure_clusters(sgraph *g,const char clname[101]){
		unsigned long long j,clnr=0;
		graph_clusters *c;
		for(j=1;j<=g->N;j++){
			if(g->node[j].node_props==NULL){
				g->node[j].node_props=new node_attr(j,0.,0);
			}
			else g->node[j].node_props->cluster_ID=0;		
		}
		//printf("Starting to burn\n");getchar();
		for(j=1;j<=g->N;j++){
			if(g->node[j].node_props->cluster_ID==0){
				clnr++;
				burn(g,j,clnr);
				//printf("cl %lld\n",clnr);
			}
		//	if(j/100==j/100.) printf("j=%lld\tcl %lld\n",j,clnr);
		}
		printf("Burning done, %lld clusters\n",clnr);//getchar();
		
		if(clname==NULL) c=new graph_clusters(g,clnr);
			else c=new graph_clusters(g,clnr,clname);
		return(c);
}

sgraph *get_giant_component(sgraph *g){
	sgraph *net;
	int i;
	char nev[400];
	snode::slinklist::slink *sc2;
	unsigned long long *sorszam,hany;
	graph_clusters *c;
	
	c=measure_clusters(g,NULL);
	sorszam=new unsigned long long[g->N+1];
	hany=0;
	for(i=1;i<=g->N;i++){
		if(g->node[i].node_props->cluster_ID==c->max_ID){
			hany++;
			sorszam[i]=hany;
		}
		else sorszam[i]=0;
	}
	sprintf(nev,"%s_giant",g->gname);
	net=new sgraph(hany,nev);
	for(i=1;i<=g->N;i++)
		if(sorszam[i]>0){
			if(g->node[i].node_draw_props!=NULL){
				net->node[sorszam[i]].node_draw_props=new node_draw_attr(g->node[i].node_draw_props->basic->ID);
				net->node[sorszam[i]].node_draw_props->basic->ID=g->node[i].node_draw_props->basic->ID;
				strcpy(net->node[sorszam[i]].node_draw_props->node_label,"");
				sprintf(net->node[sorszam[i]].node_draw_props->node_label,"%s",g->node[i].node_draw_props->node_label);
				sprintf(net->node[sorszam[i]].node_draw_props->color,"%s",g->node[i].node_draw_props->color);
				net->node[sorszam[i]].node_props=new node_attr(g->node[i].node_draw_props->basic->ID,1,1);
			}
			else net->node[sorszam[i]].node_props=new node_attr(i,1,1);
		}
	for(i=1;i<=g->N;i++)
		if(sorszam[i]>0){
			sc2=g->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				if(sc2->link_props!=NULL) 
					net->add_link(sorszam[i],sorszam[sc2->sneighbor],sc2->link_props->lw);
				else net->add_link(sorszam[i],sorszam[sc2->sneighbor]);
				sc2=sc2->next(sc2);
			}
		}
	delete c;		
	return(net);	
}

sgraph *get_Strongly_Connected_Components_linked(sgraph *g){
    unsignedintpointer *Reachable;
    unsigned int i,j,k,ok;
    longint_ARRAY *SCC;
    sgraph *SCC_g;
    
    SCC=new longint_ARRAY();
    Reachable=new unsignedintpointer[g->N+1];
    for(i=0;i<=g->N;i++){
        Reachable[i]=new unsigned int[g->N+1];
        for(j=0;j<=g->N;j++) Reachable[i][j]=0;
        Reachable[i][i]=1;
    }
    ok=0;
    while(ok==0) {
        ok=1;
        for(i=1;i<=g->N;i++){
            for(j=1;j<=g->N;j++)
                if((i!=j) && (Reachable[i][j]==0)){
                    for(k=1;k<=g->N;k++)
                        if((Reachable[i][k]==1)&&(g->node[k].Li.get_dirlink(j)))
                        {Reachable[i][j]=1; ok=0;}
                }
        }
    }
    for(i=2;i<=g->N;i++)
        for(j=1;j<i;j++)
            if ((Reachable[i][j]==1) && (Reachable[j][i]==1)) {
                SCC->add_element(i);
                SCC->add_element(j);
            }
    printf("cycling has %lld nodes\n",g->N);
     char nev[400];
    sprintf(nev,"%s_StronglyConnectedComp",g->gname);
    SCC_g=new sgraph(0,nev);
    unsigned long long int id1;
    
    for(int i=1;i<=SCC->N;i++){
        if(g->node[SCC->A[i]].node_draw_props!=NULL)
             id1 = SCC_g->add_node(g->node[SCC->A[i]].node_draw_props->node_label);
        else id1 = SCC_g->add_node("test");
        SCC_g->add_node_props(id1,g->node[SCC->A[i]].node_props);
        SCC_g->node[id1].ID = g->node[SCC->A[i]].ID;
     }
    
    for(int i=1;i<=SCC->N;i++)
        for(int k=1;k<=SCC->N;k++)
            if(g->node[SCC->A[i]].Li.get_dirlink(SCC->A[k])>0)
                SCC_g->add_link(i,k,g->get_dir_link(SCC->A[i],SCC->A[k]));
  
    printf("SCC has %lld nodes\n",SCC_g->N);
    return(SCC_g);
}

sgraph *get_Strongly_Connected_Components_Disconneceted(sgraph *g){
    unsignedintpointer *Reachable;
    unsigned int i,j,k,ok;
    longint_ARRAY *SCC;
    sgraph *SCC_g;
    
    SCC=new longint_ARRAY();
    Reachable=new unsignedintpointer[g->N+1];
    for(i=0;i<=g->N;i++){
        Reachable[i]=new unsigned int[g->N+1];
        for(j=0;j<=g->N;j++) Reachable[i][j]=0;
        Reachable[i][i]=1;
    }
    ok=0;
    while(ok==0) {
        ok=1;
        for(i=1;i<=g->N;i++){
            for(j=1;j<=g->N;j++)
                if((i!=j) && (Reachable[i][j]==0)){
                    for(k=1;k<=g->N;k++)
                        if((Reachable[i][k]==1)&&(g->node[k].Li.get_dirlink(j)))
                        {Reachable[i][j]=1; ok=0;}
                }
        }
    }
    for(i=2;i<=g->N;i++)
        for(j=1;j<i;j++)
            if ((Reachable[i][j]==1) && (Reachable[j][i]==1)) {
                SCC->add_element(i);
                SCC->add_element(j);
            }
    printf("cycling has %lld nodes\n",g->N);
     char nev[400];
    sprintf(nev,"%s_SCC_Disconnected",g->gname);
    SCC_g=new sgraph(0,nev);
    unsigned long long int id1;
    
    for(int i=1;i<=SCC->N;i++){
        if(g->node[SCC->A[i]].node_draw_props!=NULL)
             id1 = SCC_g->add_node(g->node[SCC->A[i]].node_draw_props->node_label);
        else id1 = SCC_g->add_node("test");
        SCC_g->add_node_props(id1,g->node[SCC->A[i]].node_props);
        SCC_g->node[id1].ID = g->node[SCC->A[i]].ID;
     }
    
    for(int i=1;i<=SCC->N;i++)
        for(int k=1;k<=SCC->N;k++)
            if((g->node[SCC->A[i]].Li.get_dirlink(SCC->A[k])>0) && (Reachable[SCC->A[i]][SCC->A[k]]==1) && (Reachable[SCC->A[k]][SCC->A[i]]==1))
                SCC_g->add_link(i,k,g->get_dir_link(SCC->A[i],SCC->A[k]));
  
    printf("SCC has %lld nodes\n",SCC_g->N);
    return(SCC_g);
}

sgraph *get_cluster_component(sgraph *g,graph_clusters *c,unsigned long long cl_id){
	sgraph *net;
	int i;
	char nev[400];
	snode::slinklist::slink *sc2;
	unsigned long long *sorszam,hany;
	
	
	sorszam=new unsigned long long[g->N+1];
	hany=0;
	for(i=1;i<=g->N;i++){
		if(g->node[i].node_props->cluster_ID==cl_id){
			hany++;
			sorszam[i]=hany;
		}
		else sorszam[i]=0;
	}
	sprintf(nev,"%s_cluster_%lld__size-%lld",g->gname,cl_id,c->cl_size[cl_id]);
	net=new sgraph(hany,nev);
	for(i=1;i<=g->N;i++)
		if(sorszam[i]>0){
			if(g->node[i].node_draw_props!=NULL){
				net->node[sorszam[i]].node_draw_props=new node_draw_attr(g->node[i].node_draw_props->basic->ID);
				net->node[sorszam[i]].node_draw_props->basic->ID=g->node[i].node_draw_props->basic->ID;
				strcpy(net->node[sorszam[i]].node_draw_props->node_label,"");
				sprintf(net->node[sorszam[i]].node_draw_props->node_label,"%s",g->node[i].node_draw_props->node_label);
				sprintf(net->node[sorszam[i]].node_draw_props->color,"%s",g->node[i].node_draw_props->color);
				net->node[sorszam[i]].node_props=new node_attr(g->node[i].node_draw_props->basic->ID,1,1);
				net->node[sorszam[i]].ID=i;
			}
			else {net->node[sorszam[i]].node_props=new node_attr(i,1,1);
				  net->node[sorszam[i]].ID=i;
			}
		}
	for(i=1;i<=g->N;i++)
		if(sorszam[i]>0){
			sc2=g->node[i].Li.first_link_pt;	
			while(sc2!=NULL){
				if(sc2->link_props!=NULL) 
					net->add_link(sorszam[i],sorszam[sc2->sneighbor],sc2->link_props->lw);
				else net->add_link(sorszam[i],sorszam[sc2->sneighbor]);
				sc2=sc2->next(sc2);
			}
		}
	return(net);	
}

unsigned long long follow(sgraph *gg,unsigned long long st){
	if(gg->node[st].node_props->cluster_ID!=0) {
		return(gg->node[st].node_props->cluster_ID);
	}
	else {//printf("following %lld\n",gg->node[st].Li.first_link_pt->sneighbor);getchar();
		  gg->node[st].node_props->cluster_ID=follow(gg,gg->node[st].Li.first_link_pt->sneighbor);
		 //printf("Got %lld, it's in %lld\n",st,gg->node[st].node_props->cluster_ID);getchar();
		 return(gg->node[st].node_props->cluster_ID);	
		 }
}

graph_clusters *measure_gradnet_clusters(sgraph *g,char clname[101]){
		unsigned long long j,clnr=0;;
		graph_clusters *c;
		for(j=1;j<=g->N;j++){
			if(g->node[j].Li.get_dirlink(j)){
				clnr++;
				//printf("New min %lld, cl nr %lld\n",j,clnr);getchar();
				if(g->node[j].node_props==NULL){
					g->node[j].node_props=new node_attr(j,0.,clnr);
				}
				else g->node[j].node_props->cluster_ID=clnr;
			}
			else if(g->node[j].node_props==NULL){
					g->node[j].node_props=new node_attr(j,0.,0);
				}
		}
		for(j=1;j<=g->N;j++){
			if(g->node[j].node_props->cluster_ID==0){
				//printf("Starting from %lld\n",j);getchar();
				g->node[j].node_props->cluster_ID=follow(g,j);
			}
		}
		//printf("Burning done, %lld clusters\n",clnr);getchar();
		c=new graph_clusters(g,clnr,clname);
		return(c);
}

graph_clusters *measure_gradnet_clusters(sgraph *g){
		unsigned long long j,clnr=0;;
		graph_clusters *c;
		for(j=1;j<=g->N;j++){
			if(g->node[j].Li.get_dirlink(j)){
				clnr++;
				//printf("New min %lld, cl nr %lld\n",j,clnr);getchar();
				if(g->node[j].node_props==NULL){
					g->node[j].node_props=new node_attr(j,0.,clnr);
				}
				else g->node[j].node_props->cluster_ID=clnr;
			}
			else if(g->node[j].node_props==NULL){
					g->node[j].node_props=new node_attr(j,0.,0);
				}
		}
		for(j=1;j<=g->N;j++){
			if(g->node[j].node_props->cluster_ID==0){
				//printf("Starting from %lld\n",j);getchar();
				g->node[j].node_props->cluster_ID=follow(g,j);
			}
		}
		//printf("Burning done, %lld clusters\n",clnr);getchar();
		c=new graph_clusters(g,clnr);
		return(c);
}

void transfer_clusternumber(sgraph *grad,sgraph *er){
	unsigned long long i;
	if(grad->N!=er->N) return;
	for(i=1;i<=grad->N;i++){
		if(er->node[i].node_props!=NULL)
			er->node[i].node_props->cluster_ID=grad->node[i].node_props->cluster_ID;
		else {er->node[i].node_props=new node_attr(i,0.0,grad->node[i].node_props->cluster_ID);
			 }
	}
}

void screen_cluster_props(graph_clusters *p){
	printf("Clustering properties:\tcluster_NR=%lld\n\tlargesr cluster is %lld\tsize=%lld\n",
		p->cluster_NR,p->max_ID,p->max_size);
	printf("\taverage fragment size %lf\n",p->fragments_av);	
}

void burn_directed(sgraph *gg,unsigned long long st,unsigned long long IDa){
	unsigned long long *neigh,i;
	gg->node[st].node_props->cluster_ID=IDa;
	//printf("\nCaught %lld into %lld, starting on %lld neighbors\t",st,IDa,gg->node[st].k);
	//getchar();
	neigh=new unsigned long long[gg->node[st].k+1];
	gg->get_directed_neighbors_from_nodes_list(st,neigh);
	for(i=1;i<=neigh[0];i++){
		//printf("%lld,",neigh[i]);
		if(gg->node[neigh[i]].node_props->cluster_ID==0)	
		{burn(gg,neigh[i],IDa);}
		else printf("Closed loop size %llu ending in node %lld\n",
					gg->node[neigh[i]].node_props->cluster_ID,gg->node[st].node_props->cluster_ID);
	}
	delete[] neigh;neigh=NULL;
}

sgraph  *export_levels_acyclic(sgraph *g,char pa[]){
	// linklist of node i contains the INCOMING lins to i!
	unsigned long long i;
	unsigned int *level,uress,ok,added,maxl,c;
	snode::slinklist::slink *sc2;
	sgraph *net;
	unsigned long long *sorszam,hany;
	char nm[200];
	
	uress=1;
	level=new unsigned int[g->N+1];
	for(i=1;i<=g->N;i++) level[i]=0;
	c=0;
	do{
		uress=1;added=0;
		for(i=1;i<=g->N;i++)
			if(level[i]==0){
				if((g->node[i].Li.first_link_pt==NULL)||((g->node[i].Li.first_link_pt->sneighbor==i)&&(g->node[i].Li.first_link_pt->next(g->node[i].Li.first_link_pt)==NULL))) {
					level[i]=1;added=1;c++;
					printf("%d\t Node %lld\t%s is on level 1\n",c,i,g->node[i].node_draw_props->node_label);
					//getchar();
				}
				else {
					ok=1;maxl=0;
					sc2=g->node[i].Li.first_link_pt;
					while(sc2!=NULL){
						if((level[sc2->sneighbor]==0)&&(sc2->sneighbor!=i)){
							ok=0;
							sc2=NULL;
						}
						else {
							maxl=level[sc2->sneighbor]+1;
							sc2=sc2->next(sc2);
						}
					}
					if(ok==1) {
						level[i]=maxl;added=1;c++;
						printf("%d\tNode %lld\t%s is on level %d\n",c,i,g->node[i].node_draw_props->node_label,maxl);
						//getchar();
					}
					else uress=0;
				}
			}
		if((uress==0)&&(added==0)){
			printf("Network %s is not acyclic\n",g->gname);
	
			sorszam=new unsigned long long[g->N+1];
			hany=0;
			for(i=1;i<=g->N;i++){
				if(level[i]==0){
					hany++;
					sorszam[i]=hany;
				}
				else sorszam[i]=0;
			}
			printf("core network has %lld nodes\n",hany);
			sprintf(nm,"%s_cyclic",g->gname);
			net=new sgraph(hany,nm);
			for(i=1;i<=g->N;i++)
				if(sorszam[i]>0){
					if(g->node[i].node_draw_props!=NULL){
						net->node[sorszam[i]].node_draw_props=new node_draw_attr(g->node[i].node_draw_props->basic->ID);
						net->node[sorszam[i]].node_draw_props->basic->ID=g->node[i].node_draw_props->basic->ID;
						strcpy(net->node[sorszam[i]].node_draw_props->node_label,"");
						sprintf(net->node[sorszam[i]].node_draw_props->node_label,"%s",g->node[i].node_draw_props->node_label);
						sprintf(net->node[sorszam[i]].node_draw_props->color,"%s",g->node[i].node_draw_props->color);
						net->node[sorszam[i]].node_props=new node_attr(g->node[i].node_draw_props->basic->ID,1,1);
					}
					else net->node[sorszam[i]].node_props=new node_attr(i,1,1);
				}
			for(i=1;i<=g->N;i++)
				if(sorszam[i]>0){
					sc2=g->node[i].Li.first_link_pt;	
					while(sc2!=NULL){
						if(sorszam[sc2->sneighbor]>0){
							if(sc2->link_props!=NULL) 
								net->add_link(sorszam[i],sorszam[sc2->sneighbor],sc2->link_props->lw);
							else net->add_link(sorszam[i],sorszam[sc2->sneighbor]);
						}
						sc2=sc2->next(sc2);
					}
				}
			return(net);
		}
		//getchar();
	}
	while(uress==0);
		
	//acyclic_network_draw(g,level,pa);
	return(NULL);
}
