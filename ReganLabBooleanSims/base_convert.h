void decimal_to_any_base(unsigned long long int x,unsigned short int base,unsigned int char_nr, unsigned int anybase[]){
	unsigned long long int i;
	//printf("%ld to base %d, char nr=%d\n",x,base,char_nr);getchar();
	for(i=0;i<char_nr;i++) anybase[i]=0;

	for(i=0;i<char_nr;i++){
		anybase[i]=(int)fmod(x,base);
		x-=anybase[i];
		x=x/base;
	}
}	
void decimal_to_any_base(unsigned long long int x,unsigned short int base,unsigned int char_nr, unsigned short int *anybase){
	unsigned long long int i;
	//printf("%ld to base %d, char nr=%d\n",x,base,char_nr);getchar();
	for(i=0;i<char_nr;i++) anybase[i]=0;
    
	for(i=0;i<char_nr;i++){
		anybase[i]=(int)fmod(x,base);
		x-=anybase[i];
		x=x/base;
	}
}
char *decimal_to_string(unsigned long long int x,unsigned short int base,unsigned int char_nr){
	unsigned long long int i;
	char *anybase;
	//int ch_l;
	
//	ch_l=(int)(log(x)/log(base)+1);
	anybase=new char[char_nr+1];
	
	//printf("%ld to base %d, char nr=%d\n",x,base,char_nr);getchar();
	for(i=0;i<=char_nr;i++) anybase[i]=0;
	
	for(i=0;i<char_nr;i++){
		anybase[char_nr-i-1]=(int)fmod(x,base);
		x-=anybase[char_nr-i-1];
		x=x/base;
	}
	for(i=0;i<char_nr;i++)
		anybase[char_nr-i-1]+=48;
	return(anybase);
}	

/*
void decimal_to_nbase(unsigned long long int x,unsigned short int base,unsigned int leng,char res[]){
	//char *res;
	int i,po;
	unsigned int r;
	sprintf(res,"");
	if(x==0) {for(i=1;i<=leng;i++) sprintf(res,"%s0",res);return;}
	//char_nr=(int)ceil(log(x)/log(base));
	//if(char_nr>leng) {printf("given base character nr too small!\n");getchar();return;}
	//printf("%ld to base %d, char nr=%d\n",x,base,char_nr);//getchar();
	//res=new char[(int)(log(x)/log(base))+3];
	if(base<2) {printf("Only base 2 to 10!!!\n");getchar();return;}
	//for(i=1;i<=leng-char_nr;i++) sprintf(res,"%s0",res);
	if(base==10) {sprintf(res,"%ld",x);	return;}
	else {
		for(i=0;i<leng;i++){
			sprintf(res,"%c%s",(int)fmod(x,base)+48,res);
			x-=(int)fmod(x,base);
			x=x/base;
		}
		//for(i=leng-1;i>=0;i--){
		//	po=(int)pow(base,i);
		//	r=(int)(floor(x/pow(base,i)));
		//	//printf("%d div %d is %d\n",x,(int)pow(base,i),(int)(floor(x/pow(base,i))));
		//	sprintf(res,"%s%c",res,r+48);
		//	x-=r*po;
		//	//printf("i=%d res=%s\n",i,res);getchar();
		//}
		//printf("res=%s\tleng=%d\n",res,leng);getchar();
	    return;	
	}	
}
*/

unsigned long long int nbase_to_decimal(char *code,unsigned short int base){
	unsigned long long int numb;
	int i,char_nr;
	//printf("code=%s\n",code);
	if(base<2) {return(0);}
	char_nr=0;i=0;numb=0;
	while(code[i]!=0){
		if(code[i]-48>=base) {printf("Only base 2 to 10!!!\n\tLarger digit found: i=%d code-48=%d char=%c\n",
										i,code[i],code[i]);getchar();return(0);}
		else {char_nr++;
			  //printf("%d %c\n",i,code[i]);getchar();
			 }
		i++;	 
	}
	//printf("char_nr=%d\n",char_nr);
	for(i=1;i<=char_nr;i++)
		{numb+=(code[char_nr-i]-48)*(int)pow(base,i-1);
			//printf("adder %d times %d\t%d\n",code[char_nr-i]-48,(int)pow(base,i-1),numb);
			//getchar();
		}
	return(numb);	
}

unsigned long long int anybase_to_decimal(unsigned int *code,unsigned short int n,unsigned short int base){
	unsigned long long int numb,p;
	int i,j;
	//printf("code=%s\n",code);
	if(base<2) {return(0);}
	i=0;numb=0;
	for(i=0;i<n;i++){
		if(code[i]>=base) {printf("Only base %d!!!\n\tLarger digit found: i=%d code=%d\n",
										base,i,code[i]);getchar();return(0);}
	}
	//printf("char_nr=%d\n",char_nr);
	for(i=1;i<=n;i++)
		{p=1;
		 for(j=1;j<=i-1;j++) p*=base;
		  numb+=p*(code[i-1]);
		 //printf("adder %d times %d\t%d\n",code[n-i],(int)pow(base,i-1),numb);
			//getchar();
		}
	return(numb);	
}
unsigned long long int anybase_to_decimal(unsigned short int *code,unsigned short int n,unsigned short int base){
	unsigned long long int numb,p;
	int i,j;
	//printf("code=%s\n",code);
	if(base<2) {return(0);}
	i=0;numb=0;
	for(i=0;i<n;i++){
		if(code[i]>=base) {printf("Only base %d!!!\n\tLarger digit found: i=%d code=%d\n",
                                  base,i,code[i]);getchar();return(0);}
	}
	//printf("char_nr=%d\n",char_nr);
	for(i=1;i<=n;i++)
    {p=1;
        for(j=1;j<=i-1;j++) p*=base;
        numb+=p*(code[i-1]);
        //printf("adder %d times %d\t%d\n",code[n-i],(int)pow(base,i-1),numb);
        //getchar();
    }
	return(numb);	
}

