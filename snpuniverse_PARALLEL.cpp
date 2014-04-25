/*
 *Last Update: April 25, 2014
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
 *SNP_universe with multithreading by chromosome
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sys/resource.h>
#include <unordered_map>
#include <vector>
#include <gzstream.h>
#include <time.h>
#include <pthread.h>
#include <sched.h> //for CPU affinity
#include <unistd.h> //for sysconf

using namespace std;

struct threadData{
    string dirpath;
    vector<string> vcf_list;
    vector<string> chrom;
    unordered_map<string,int> chr;
    unordered_map<string,int> samples;
    char* snp;
    char* ref;
    int* count;
    double* depth;
    double* qs;
    int tid;
};

static const char *options="p:P:i:I:x:X:cm";
static string outfile;
static string path;
static string input;
bool compressed=false;
bool multioutput=false;
int SIZE=45000000;
int SAMPLECOUNT=1001;
int NUMTHREADS = 2;
int DEPTH = 1;

void parseArgs(int argc, char**argv){
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
	case 'p':
	case 'P': path = optarg; break; //complete path to the input
        case 'i':
        case 'I': input = optarg; break; //directory
	case 'x':
	case 'X': outfile = optarg; break; //output file
        case 'c': compressed =true; break; //.gz files
        case 'm': multioutput=true; break; //multiple outputs files by chrom
	default: break;
        } // switch
    } // while
} // parseArgs

int checkAlt(char ref,string alt,void *thread_data,int snppos){
    threadData *t = (threadData*) thread_data;
    int temp1=0,temp2=1;
    if(!alt.compare(".")){ return 0;} //monomorphic
    else if(alt.size()==1 && alt[0]!=ref){ 
	if(t->snp[snppos]!='B' && t->snp[snppos]!='C'){ //disregard all positions with indels or structural variants
	    t->snp[snppos]='A';
	    t->ref[snppos]=ref;
  	    return 1;
	}else if(t->snp[snppos]=='B'){
	    t->snp[snppos]='C'; //SNP+Indels/strucVariants
	    t->ref[snppos]=ref;
  	    return 1;
	}
    }

    if(alt.size()>1){ 
	while(temp1<temp2){
	    temp2=alt.find_first_of(",",temp1+1);
	    if(temp2-temp1>2 || temp2==alt.npos){ //indels,structural variants
		if(t->snp[snppos]=='A') t->snp[snppos]='C'; 
		else if(t->snp[snppos]!='C') t->snp[snppos]='B';
		return 0; 
	    }
	    temp1=temp2+1;
	}
    }
    if(t->snp[snppos]!='B' && t->snp[snppos]!='C'){
	t->snp[snppos]='A';
	t->ref[snppos]=ref;
    }else if(t->snp[snppos]=='B'){
	    t->snp[snppos]='C'; //SNP+Indels/strucVariants
	    t->ref[snppos]=ref;
    }
    return 1; //multiple ALTs
}

int locateSNP_1(string filepath,void *thread_data){
    threadData *t = (threadData*) thread_data;
    string linestream,temp,alt;
    int idx1=0,idx2=0,snppos=0,chrpos=0;
    char ref;
    
    ifstream fp(filepath.c_str());
    if (!fp.is_open()) { //check input file
	printf("ERROR: Failed to open the input file %s", filepath.c_str());
	return 1;
    }
    
    for(int x=0;getline(fp,linestream);x++){
	if(linestream[0]!='#'){
	    idx1 = linestream.find_first_of("\t",idx1+1); //first column
            temp=linestream.substr(0,idx1); //get the chrom number
            idx2 = linestream.find_first_of("\t",++idx1); //second column
            snppos = atoi(linestream.substr(idx1,idx2-idx1).c_str()); //get the SNP pos
            chrpos = t->chr.find(temp)->second; //printf("%s %d %d\n",temp.c_str(),chrpos,snppos);
	    idx1 = linestream.find_first_of("\t",idx2+1); //skip id column
    	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
    
	    if(idx2-idx1==2){ //skip indels and structural variants
		ref=linestream[idx2-1]; 
		alt=linestream.substr(idx2+1,linestream.find_first_of("\t",idx2+1)-idx2-1); 
	        //printf("%d %c %s %d\t",snppos,ref,alt.c_str(),snp[chrpos][snppos]);
		checkAlt(ref,alt,t,snppos); //check if SNP occur in a position
                //printf("%d\n",snp[chrpos][snppos]);
                
	    }else{
		t->snp[snppos]='B'; //IGNORE INDEL(deletion)
	    }
	}
	idx1=0;
    }
    fp.close();
    return 0;
}

int locateSNP_2(string filepath,void *thread_data){
    threadData *t = (threadData*) thread_data;
    string linestream,temp,alt;
    int idx1=0,idx2=0,snppos=0,chrpos=0;
    char ref;
    
    igzstream fp(filepath.c_str());
    if (!fp.good()) { //check input file
	printf("ERROR: Failed to open the input file %s", filepath.c_str());
	return 1;
    }
    
    for(int x=0;getline(fp,linestream);x++){
	if(linestream[0]!='#'){
	    idx1 = linestream.find_first_of("\t",idx1+1); //first column
            temp=linestream.substr(0,idx1); //get the chrom number
            idx2 = linestream.find_first_of("\t",++idx1); //second column
            snppos = atoi(linestream.substr(idx1,idx2-idx1).c_str()); //get the SNP pos
            chrpos = t->chr.find(temp)->second; //printf("%s %d %d\n",temp.c_str(),chrpos,snppos);
	    idx1 = linestream.find_first_of("\t",idx2+1); //skip id column
    	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
    
	    if(idx2-idx1==2){ //skip indels and structural variants
		ref=linestream[idx2-1]; 
		alt=linestream.substr(idx2+1,linestream.find_first_of("\t",idx2+1)-idx2-1); 
	        //printf("%d %c %s %d\t",snppos,ref,alt.c_str(),snp[chrpos][snppos]);
		checkAlt(ref,alt,t,snppos); //check if SNP occur in a position
                //printf("%d\n",snp[chrpos][snppos]);
                
	    }else{
		t->snp[snppos]='B'; //IGNORE INDEL(deletion)
	    }
	}
	idx1=0;
    }
    fp.close();
    return 0;
}

void *readFolder(void *thread_data){
    DIR *dp;
    threadData *t = (threadData*) thread_data;
    string temp,exten,filename,curpath,tempchrom;
    struct dirent *ep;
    time_t start,end;
    int dotpos;

    int numcor = sysconf(_SC_NPROCESSORS_ONLN);
    if(t->tid>=numcor){
	printf("ERROR:Insufficient number of cores.\n");
	exit(EXIT_FAILURE);
    }
	
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(t->tid+1,&cpuset);

    sched_setaffinity(0,sizeof(cpuset),&cpuset);
    
    dp = opendir(t->dirpath.c_str());
    if(t->dirpath[t->dirpath.size()-1]!='/') curpath = t->dirpath + '/';
    struct stat filestat;
    if(dp!=NULL){
	while(ep = readdir(dp)){
            time(&start);
	    if(!strcmp(ep->d_name,".") || !strcmp(ep->d_name,"..")) continue;
	    temp = curpath + ep->d_name;
            stat(temp.c_str(), &filestat);
	    filename=ep->d_name;
	    if(S_ISDIR(filestat.st_mode)){ //recursively read directories
                t->dirpath=temp;
		readFolder(thread_data);
                t->dirpath=curpath;
	    }else if(S_ISREG(filestat.st_mode)){
		//check if vcf
		dotpos=filename.find_last_of(".");
		if(dotpos!=filename.npos){
		    exten=filename.substr(dotpos);
		    tempchrom = t->chrom[t->tid] + ".";
		    if(!exten.compare(".vcf") && !compressed && filename.find(tempchrom)!=filename.npos){
		        printf("thread%d -  %s\t",t->tid,filename.c_str());
			t->vcf_list.push_back(temp);
			locateSNP_1(temp,thread_data);
			time(&end);
            		printf("Time: %.f sec\n",difftime(end,start));
		    }else if(!exten.compare(".gz") && compressed && filename.find(tempchrom)!=filename.npos){
			exten=filename.substr(filename.find_last_of(".",dotpos-1)); 
		        if(!exten.compare(".vcf.gz")){
			    printf("thread%d -  %s\t",t->tid,filename.c_str());
			    t->vcf_list.push_back(temp);
			    locateSNP_2(temp,thread_data);
			    time(&end);
            		    printf("Time: %.f sec\n",difftime(end,start));
			}
	            }
		}
		
	    }
	    temp.clear();
	    filename.clear();
	}
	closedir (dp);
    }else{
	printf("ERROR: Can't find the directory.");
    }
}

int printSNPlist_1(int tid,void *thread_data,unordered_map<string,int> samples,FILE *output){
    threadData *t = (threadData*) thread_data;
    string linestream,alt,temp1,formatfield,formatval,temp2,samname;
    char ref,*token=NULL,tok_ar[40];
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0,chrpos=0,setsize=samples.size();
    int AD[128];
    if(tid==0){
    	for(int i=0;i<setsize;i++){ 
  	    for(int j=0;j<DEPTH;j++)idx2=t->vcf_list[i].find_last_of("/");
            idx1=t->vcf_list[i].find_last_of("/",idx2-1)+1;  
            samname = t->vcf_list[i].substr(idx1,idx2-idx1); 
            //sample list
	    fprintf(output,">%d %s\n",samples.find(samname)->second,samname.c_str()); 
    	}
    	fprintf(output,"SampleID\tChrom\tPos\tA\tT\tC\tG\n");//header
    }

    //count SNPs
    int snpcount=0;
    for(int x=0;x<SIZE;x++){
	if(t->snp[x]=='A') snpcount++;
    }

    printf("%s: %d SNPs.\n",t->chrom[tid].c_str(),snpcount);

    for(int i=0;i<setsize;i++){
    	ifstream fp(t->vcf_list[i].c_str());
    	if (!fp.is_open()) { //check input file
	    printf("ERROR: Failed to open the input file %s", t->vcf_list[i].c_str());
	    return 1;
    	}
        
  	//get sample name
       	for(int j=0;j<DEPTH;j++)idx2=t->vcf_list[i].find_last_of("/");
        idx1=t->vcf_list[i].find_last_of("/",idx2-1)+1;  
        samname = t->vcf_list[i].substr(idx1,idx2-idx1); 
    
        idx1=0;
    	for(int x=0;getline(fp,linestream);x++){
	    if(linestream[0]!='#'){
                idx1 = linestream.find_first_of("\t",idx1+1); //first column
             	temp1 = linestream.substr(0,idx1); //get the chrom number
            	idx2 = linestream.find_first_of("\t",++idx1); //second column
            	snppos = atoi(linestream.substr(idx1,idx2-idx1).c_str()); //get the SNP pos
                chrpos = t->chr.find(temp1)->second;
              
            	if(t->snp[snppos]=='A'){
		    fprintf(output,"%d\t%s\t%d\t",samples.find(samname)->second,temp1.c_str(),snppos);
                    idx1 = linestream.find_first_of("\t",idx2+1); //skip ID
            	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
		    if(idx2-idx1==2) ref=linestream[idx2-1];
		    else ref='.';
		    idx1=linestream.find_first_of("\t",++idx2);
	 	    alt=linestream.substr(idx2,idx1-idx2); //alt
	    	    for(int y=0;y<3;y++)idx1 = linestream.find_first_of("\t",idx1+1); //3columns
		    idx2 = linestream.find_first_of("\t",++idx1);
		    //format field IDs
    	    	    formatfield=linestream.substr(idx1,idx2-idx1); 
		    //format value 
		    formatval=linestream.substr(idx2+1); 
		    idx1=idx3=0;
	            while(idx1!=formatfield.npos){   
			idx2 = formatfield.find_first_of(":\0",++idx1);
			temp2=formatfield.substr(idx1,idx2-idx1); //AD field
			idx4 = formatval.find_first_of(":\0",++idx3);
			if(!temp2.compare("AD")){
			    if((!checkAlt(ref,alt,thread_data,snppos) && alt[0]!='.') || ref=='.'){
				fprintf(output,"0\t0\t0\t0 (Indel/Structural Variant)\n"); 
				printf("Indel/SV:Sample %d %s %c\n",i,temp1.c_str(),t->snp[snppos]);
			    }else{ //printf("%s\n",formatval.c_str());
                                AD['A']=AD['T']=AD['C']=AD['G']=0;
				strcpy(tok_ar,formatval.substr(idx3,idx4-idx3).c_str());//fprintf(output,"ref:%c %s %s\t",ref,alt.c_str(),tok_ar);
				token=strtok(tok_ar,","); //fprintf(output,"%s\n",token);
				AD[ref]=atoi(token);  //count ref AD
				token = strtok(NULL,","); //fprintf(output,"%s\t",token);
				for(int j=0;token!=NULL;j+=2){
				    AD[alt[j]]=atoi(token);
				    token = strtok(NULL,",");
			        }
			        fprintf(output,"%d\t%d\t%d\t%d\n",AD['A'],AD['T'],AD['C'],AD['G']);   
				
			    }
			    break;          
			}
			idx1=idx2;
			idx3=idx4;
            	    }
		    if(temp2.compare("AD")){ //Not a SNP in other sample
			fprintf(output,"0\t0\t0\t0\n");
		    }
	       	}
                idx1=0;
	    }
    	}
	fp.close();
    }
    
    return 0;
}

int printSNPlist_2(int tid,void *thread_data,unordered_map<string,int> samples,FILE *output){
    threadData *t = (threadData*) thread_data;
    string linestream,alt,temp1,formatfield,formatval,temp2,samname;
    char ref,*token=NULL,tok_ar[40];
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0,chrpos=0,setsize=samples.size();
    int AD[128],stack[4];
    if(tid==0){
    	for(int i=0;i<setsize;i++){ 
  	    for(int j=0;j<DEPTH;j++)idx2=t->vcf_list[i].find_last_of("/");
            idx1=t->vcf_list[i].find_last_of("/",idx2-1)+1;  
            samname = t->vcf_list[i].substr(idx1,idx2-idx1); 
            //sample list
	    fprintf(output,">%d %s\n",samples.find(samname)->second,samname.c_str()); 
    	}
    	fprintf(output,"SampleID\tChrom\tPos\tA\tT\tC\tG\n");//header
    }

    //count SNPs
    int snpcount=0;
    for(int x=0;x<SIZE;x++){
	if(t->snp[x]=='A') snpcount++;
    }

    printf("%s: %d SNPs.\n",t->chrom[tid].c_str(),snpcount);

    for(int i=0;i<setsize;i++){
    	igzstream fp(t->vcf_list[i].c_str());
    	if (!fp.good()) { //check input file
	    printf("ERROR: Failed to open the input file %s", t->vcf_list[i].c_str());
	    return 1;
    	}
        
  	//get sample name
       	for(int j=0;j<DEPTH;j++)idx2=t->vcf_list[i].find_last_of("/");
        idx1=t->vcf_list[i].find_last_of("/",idx2-1)+1;  
        samname = t->vcf_list[i].substr(idx1,idx2-idx1); 
    
        idx1=0;
    	for(int x=0;getline(fp,linestream);x++){
	    if(linestream[0]!='#'){
                idx1 = linestream.find_first_of("\t",idx1+1); //first column
             	temp1 = linestream.substr(0,idx1); //get the chrom number
            	idx2 = linestream.find_first_of("\t",++idx1); //second column
            	snppos = atoi(linestream.substr(idx1,idx2-idx1).c_str()); //get the SNP pos
                chrpos = t->chr.find(temp1)->second;
              
            	if(t->snp[snppos]=='A'){
		    fprintf(output,"%d\t%s\t%d\t",samples.find(samname)->second,temp1.c_str(),snppos);
                    idx1 = linestream.find_first_of("\t",idx2+1); //skip ID
            	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
		    if(idx2-idx1==2) ref=linestream[idx2-1];
		    else ref='.';
		    idx1=linestream.find_first_of("\t",++idx2);
	 	    alt=linestream.substr(idx2,idx1-idx2); //alt
	    	    for(int y=0;y<3;y++)idx1 = linestream.find_first_of("\t",idx1+1); //3columns
		    idx2 = linestream.find_first_of("\t",++idx1);
		    //format field IDs
    	    	    formatfield=linestream.substr(idx1,idx2-idx1); 
		    //format value 
		    formatval=linestream.substr(idx2+1); 
		    idx1=idx3=0;
	            while(idx1!=formatfield.npos){   
			idx2 = formatfield.find_first_of(":\0",++idx1);
			temp2=formatfield.substr(idx1,idx2-idx1); //AD field
			idx4 = formatval.find_first_of(":\0",++idx3);
			if(!temp2.compare("AD")){
			    if((!checkAlt(ref,alt,thread_data,snppos) && alt[0]!='.') || ref=='.'){
				fprintf(output,"0\t0\t0\t0 (Indel/Structural Variant)\n"); 
				printf("Indel/SV:Sample %d %s %c\n",i,temp1.c_str(),t->snp[snppos]);
			    }else{ //printf("%s\n",formatval.c_str());
                                AD['A']=AD['T']=AD['C']=AD['G']=0;
				strcpy(tok_ar,formatval.substr(idx3,idx4-idx3).c_str());//fprintf(output,"ref:%c %s %s\t",ref,alt.c_str(),tok_ar);
				token=strtok(tok_ar,","); //fprintf(output,"%s\n",token);
				AD[ref]=atoi(token);  //count ref AD
				token = strtok(NULL,","); //fprintf(output,"%s\t",token);
				for(int j=0;token!=NULL;j+=2){
				    AD[alt[j]]=atoi(token);
				    token = strtok(NULL,",");
			        }
			        fprintf(output,"%d\t%d\t%d\t%d\n",AD['A'],AD['T'],AD['C'],AD['G']);   
				
			    }
			    break;          
			}
			idx1=idx2;
			idx3=idx4;
            	    }
		    if(temp2.compare("AD")){ //Not a SNP in other sample
			fprintf(output,"0\t0\t0\t0\n");
		    }
	       	}
                idx1=0;
	    }
    	}
	fp.close();
    }
    
    return 0;
}

void *multiprint_2(void *thread_data){
    threadData *t = (threadData*) thread_data;
    string linestream,alt,temp1,formatfield,formatval,temp2,samname,snpname; 
    char ref,*token=NULL,tok_ar[40];
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0,chrpos=0,DP=0,setsize=t->samples.size();
    int AD[128];
    float qual;
    string pad("00000000");

    FILE *output1,*output2,*output3,*output4; 
    temp1=outfile+"_SNP_AD_"+ t->chrom[t->tid] +".txt"; 
    output1=fopen(temp1.c_str(),"w"); 
    temp1=outfile+"_SNP_POS_"+ t->chrom[t->tid] +".txt"; 
    output2=fopen(temp1.c_str(),"w"); 
    temp1=outfile+"_INDELorSTRUCvar_AD_"+ t->chrom[t->tid] +".txt"; 
    output3=fopen(temp1.c_str(),"w"); 
    temp1=outfile+"_INDELorSTRUCvar_POS_"+ t->chrom[t->tid] +".txt"; 
    output4=fopen(temp1.c_str(),"w"); 

    int numcor = sysconf(_SC_NPROCESSORS_ONLN);
    if(t->tid>=numcor){
	printf("ERROR:Insufficient number of cores.\n");
	exit(EXIT_FAILURE);
    }
	
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(t->tid+1,&cpuset);

    sched_setaffinity(0,sizeof(cpuset),&cpuset);

    	for(int i=0;i<setsize;i++){ 
  	    for(int j=0;j<DEPTH;j++)idx2=t->vcf_list[i].find_last_of("/");
            idx1=t->vcf_list[i].find_last_of("/",idx2-1)+1;  
            samname = t->vcf_list[i].substr(idx1,idx2-idx1); 
            //sample list
	    fprintf(output1,">%d %s\n",t->samples.find(samname)->second,samname.c_str()); 
    	}
    	fprintf(output1,"SampleID\tSNP_ID\tA\tC\tG\tT\tQS\n");//header

    //count and print SNP position
    int snpcount=0,padcount=0;
    bool aSNP = false;
    string chrpad;
    if(t->tid<=8) chrpad = "0"+to_string(static_cast<long long>(t->tid+1));
    else chrpad = to_string(static_cast<long long>(t->tid+1));
    fprintf(output2,"SNP_ID\tChr\tPos\tSource\tRef\n");
    fprintf(output3,"SNP_ID\tChr_Pos\tRef\tAlt\tQual\tAD\n");
    fprintf(output4,"SNP_ID\tChr\tPos\n");
    //SNP_ID:[ref](1)[chr1](2)[pos](8)
    for(int x=0;x<SIZE;x++){
	if(t->snp[x]=='A'){ 
	    snpcount++;
	    snpname = to_string(static_cast<long long>(x));
	    padcount = snpname.length();
	    if(padcount<8) snpname = "1"+chrpad+pad.substr(0,8-padcount)+snpname;
	    else snpname = "1"+chrpad+snpname;
	    fprintf(output2,"%s\t%s\t%d\t\t%c\n",snpname.c_str(),t->chrom[t->tid].c_str(),x,t->ref[x]);
	}else if(t->snp[x]=='C'){
            snpname = to_string(static_cast<long long>(x));
	    padcount = snpname.length();
	    if(padcount<8) snpname = "1"+chrpad+pad.substr(0,8-padcount)+snpname;
	    else snpname = "1"+chrpad+snpname;
	    fprintf(output4,"%s\t%s\t%d\n",snpname.c_str(),t->chrom[t->tid].c_str(),x);
	}
    }

    printf("%s: %d SNPs.\n",t->chrom[t->tid].c_str(),snpcount);

    for(int i=0;i<setsize;i++){
    	igzstream fp(t->vcf_list[i].c_str());
    	if (!fp.good()) { //check input file
	    printf("ERROR: Failed to open the input file %s", t->vcf_list[i].c_str());
	    exit(EXIT_FAILURE);
    	}
        
  	//get sample name
       	for(int j=0;j<DEPTH;j++)idx2=t->vcf_list[i].find_last_of("/");
        idx1=t->vcf_list[i].find_last_of("/",idx2-1)+1;  
        samname = t->vcf_list[i].substr(idx1,idx2-idx1); 
    
        idx1=0;
    	for(int x=0;getline(fp,linestream);x++){
	    if(linestream[0]!='#'){
                idx1 = linestream.find_first_of("\t",idx1+1); //first column
             	temp1 = linestream.substr(0,idx1); //get the chrom number
            	idx2 = linestream.find_first_of("\t",++idx1); //second column
            	snppos = atoi(linestream.substr(idx1,idx2-idx1).c_str()); //get the SNP pos
                chrpos = t->chr.find(temp1)->second+1;
              
            	if(t->snp[snppos]=='A'){
                    snpname = to_string(static_cast<long long>(snppos));
		    padcount = snpname.length(); //pad 0's
                    //concatenation of Chr1 and Pos
		    if(padcount<8) snpname = "1"+chrpad+pad.substr(0,8-padcount)+snpname;
	            else snpname = "1"+chrpad+snpname;
		    fprintf(output1,"%d\t%s\t",t->samples.find(samname)->second,snpname.c_str());
                    idx1 = linestream.find_first_of("\t",idx2+1); //skip ID
            	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
		    if(idx2-idx1==2) ref=linestream[idx2-1];
		    else ref='.';
		    idx1=linestream.find_first_of("\t",++idx2);
	 	    alt=linestream.substr(idx2,idx1-idx2); //alt
		    idx2=linestream.find_first_of("\t",++idx1);
		    temp1 = linestream.substr(idx1,idx2-idx1);
		    if(idx2!=idx1 && temp1.find_first_not_of("0123456789.") == string::npos){ 
			qual=atof(temp1.c_str());//qual
		    }else{
			qual=0;
			printf("WARNING: Invalid QUAL for \"%s\" at position: %s:\"%s\"\n",samname.c_str(),snpname.c_str(),temp1.c_str()); 
		    }
		    idx1=idx2;
	    	    for(int y=0;y<2;y++)idx1 = linestream.find_first_of("\t",idx1+1); //3columns
		    idx2 = linestream.find_first_of("\t",++idx1);
		    //format field IDs
    	    	    formatfield=linestream.substr(idx1,idx2-idx1); 
		    //format value 
		    formatval=linestream.substr(idx2+1); 

 		    //print allele depth
		    idx1=idx3=0;
	            while(idx1!=formatfield.npos){   
			idx2 = formatfield.find_first_of(":\0",++idx1); 
			temp2=formatfield.substr(idx1,idx2-idx1); //get field ID
			idx4 = formatval.find_first_of(":\0",++idx3); 
			if(!temp2.compare("AD")){
			    /*(if((!checkAlt(ref,alt,thread_data,snppos) && alt[0]!='.') || ref=='.'){
				fprintf(output1,"0\t0\t0\t0 (Indel/Structural Variant)\n"); 
				printf("Indel/SV:Sample %d %s %c\n",i,.c_str(),t->snp[snppos]);
			    }else{ //printf("%s\n",formatval.c_str());*/
                                AD['A']=AD['T']=AD['C']=AD['G']=0;
				strcpy(tok_ar,formatval.substr(idx3,idx4-idx3).c_str());//get field values
				token=strtok(tok_ar,","); //fprintf(output1,"%s\n",token);
				AD[ref]=atoi(token);  //count ref AD
				token = strtok(NULL,","); //fprintf(output1,"%s\t",token);
				for(int j=0;token!=NULL;j+=2){
				    AD[alt[j]]=atoi(token);
				    token = strtok(NULL,",");
			        }
			        fprintf(output1,"%d\t%d\t%d\t%d\t%.2f\n",AD['A'],AD['C'],AD['G'],AD['T'],qual);   
				t->count[snppos]++;	
				aSNP = true;		
				
			    //}       
			}else if(!temp2.compare("DP")){
			    temp1 = formatval.substr(idx3,idx4-idx3);
			    if(idx4!=idx3 && temp1.find_first_not_of("0123456789") == string::npos) 
				DP=atoi(temp1.c_str());
			    else{
			        DP=0;
				printf("WARNING: Invalid depth(DP) for \"%s\" at position: %s:\"%s\"\n",samname.c_str(),snpname.c_str(),temp1.c_str()); 
			    }
			}
			idx1=idx2;
			idx3=idx4;
            	    }
		    if(!aSNP){ //Not a SNP in other sample
			AD['A']=AD['T']=AD['C']=AD['G']=0;
			AD[ref]=DP; 
			fprintf(output1,"%d\t%d\t%d\t%d\t%.2f\n",AD['A'],AD['C'],AD['G'],AD['T'],qual);
		    }
		    if(DP>0) t->depth[DP]++;
                    else t->depth[0]++;
			
		    if(qual>0) t->qs[(int)qual/5]++;
                    else t->qs[0]++;
	       	}else if(t->snp[snppos]=='C'){
		    snpname = to_string(static_cast<long long>(snppos));
		    padcount = snpname.length(); //pad 0's
                    //concatenation of Chr1 and Pos
		    if(padcount<8) snpname = "1"+chrpad+pad.substr(0,8-padcount)+snpname;
	            else snpname = "1"+chrpad+snpname;
		    fprintf(output3,"%d\t%s\t",t->samples.find(samname)->second,snpname.c_str());
                    idx1 = linestream.find_first_of("\t",idx2+1); //skip ID
            	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
		    fprintf(output3,"%s\t",linestream.substr(idx1,idx2-idx1).c_str());
		    idx1=linestream.find_first_of("\t",++idx2);
	 	    alt=linestream.substr(idx2,idx1-idx2); //alt
		    idx2=linestream.find_first_of("\t",++idx1);
		    temp1 = linestream.substr(idx1,idx2-idx1);
		    if(idx2!=idx1 && temp1.find_first_not_of("0123456789.") == string::npos){ 
			qual=atof(temp1.c_str());//qual
		    }else{
			qual=0;
			printf("WARNING: Invalid QUAL for \"%s\" at position: %s:\"%s\"\n",samname.c_str(),snpname.c_str(),temp1.c_str()); 
		    }
                    fprintf(output3,"%s\t%.2f\t",alt.c_str(),qual);
                    idx1=idx2;
                    for(int y=0;y<2;y++)idx1 = linestream.find_first_of("\t",idx1+1); //3columns
		    idx2 = linestream.find_first_of("\t",++idx1);
		    //format field IDs
    	    	    formatfield=linestream.substr(idx1,idx2-idx1); 
		    //format value 
		    formatval=linestream.substr(idx2+1); 
		  
 		    //print allele depth
		    idx1=idx3=0;
	            while(idx1!=formatfield.npos){   
			idx2 = formatfield.find_first_of(":\0",++idx1); 
			temp2=formatfield.substr(idx1,idx2-idx1); //get field ID
			idx4 = formatval.find_first_of(":\0",++idx3); 
			if(!temp2.compare("AD")){
			     fprintf(output3,"%s\n",formatval.substr(idx3,idx4-idx3).c_str());	
			     aSNP=true;	   
			     break;
			}
			idx1=idx2;
			idx3=idx4;
            	    }
		    if(aSNP==false)fprintf(output3,"\n");//monomorphic, no AD
		}
                idx1=DP=0;
		aSNP=false;
	    }
    	}
	fp.close();
    }
    fclose(output1);
    fclose(output2);
    fclose(output3);
    fclose(output4);
}

int main(int argc, char **argv)
{
    const rlim_t STACK_SIZE = 1000*1024*1024*2; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);
    double *tally1,*tally2,*tally3;
    time_t start, end;
    pthread_t  threads[NUMTHREADS];
    threadData thread_data[NUMTHREADS];
    unordered_map<string,int> chrmap,samples,gid_map;
    vector<string> chrnames; 

    time(&start);
    parseArgs(argc,argv);
    if (path.empty() || input.empty() || outfile.empty()) {
        cerr << "Usage:\n" << *argv
                  << " -p path to directory\n"
		  << " -i directory name\n"
		  << " -x output file\n"
		  << " -c compressed VCFs (vcf.gz)\n"
		  << " -m multiple output files\n"
                  << endl;
        return 1;
    }

    //load list of GIDs
    int pos1,pos2,gid;
    string linestream;
    ifstream list("list_3k.csv");
    if(!list.is_open()){
	printf("ERROR: Cannot open the GID list.");
	return 1;
    }
    for(int x=0;getline(list,linestream)!=NULL;x++){
	gid=0;
	pos1 = linestream.find_first_of(",") +1;
        pos2 = linestream.find_first_of(",", pos1);
	if(pos1==pos2) continue; // skip null BOX_POSITION_CODE
	gid = atoi(linestream.substr(pos2+1,linestream.find_first_of(",",pos2+1)-pos2-1).c_str());
	gid_map.insert(pair<string,int>(linestream.substr(pos1,pos2-pos1),gid));
    }
    
    tally1=(double*)calloc(SAMPLECOUNT,sizeof(double));
    tally2=(double*)calloc(SAMPLECOUNT,sizeof(double));
    tally3=(double*)calloc(SAMPLECOUNT,sizeof(double));
    chrmap["Chr1"]=0;
    chrmap["Chr2"]=1;
    chrmap["Chr3"]=2;
    chrmap["Chr4"]=3;
    chrmap["Chr5"]=4;
    chrmap["Chr6"]=5;
    chrmap["Chr7"]=6;
    chrmap["Chr8"]=7;
    chrmap["Chr9"]=8;
    chrmap["Chr10"]=9;
    chrmap["Chr11"]=10;
    chrmap["Chr12"]=11;
    chrmap["ChrSy"]=12;
    chrmap["ChrUn"]=13;

    chrnames.push_back("Chr1");
    chrnames.push_back("Chr2");
    chrnames.push_back("Chr3");
    chrnames.push_back("Chr4");
    chrnames.push_back("Chr5");
    chrnames.push_back("Chr6");
    chrnames.push_back("Chr7");
    chrnames.push_back("Chr8");
    chrnames.push_back("Chr9");
    chrnames.push_back("Chr10");
    chrnames.push_back("Chr11");
    chrnames.push_back("Chr12");
    chrnames.push_back("ChrSy");
    chrnames.push_back("ChrUn");

    /*THREADING of the comparison*/	
    for(int i=0;i<NUMTHREADS; i++){
        thread_data[i].dirpath=path+input;
    	thread_data[i].chr = chrmap;
        thread_data[i].chrom = chrnames;
    	thread_data[i].snp = (char*)calloc(SIZE,sizeof(char));
	thread_data[i].ref = (char*)calloc(SIZE,sizeof(char));
	thread_data[i].count = (int*)calloc(SIZE,sizeof(int));
	thread_data[i].depth = (double*)calloc(SAMPLECOUNT,sizeof(double));
        thread_data[i].qs = (double*)calloc(SAMPLECOUNT,sizeof(double));
	if(thread_data[i].snp==NULL){ printf("Not enough space for snp array."); return 1;}
        if(thread_data[i].ref==NULL){ printf("Not enough space for ref array."); return 1;}
	if(thread_data[i].count==NULL){ printf("Not enough space for count array."); return 1;}
	if(thread_data[i].depth==NULL){ printf("Not enough space for depth array."); return 1;}
        if(thread_data[i].qs==NULL){ printf("Not enough space for qs array."); return 1;}
    	thread_data[i].tid = i;
	pthread_create(&threads[i],NULL,readFolder, (void*) &thread_data[i]);
    }

    for(int i=0;i<NUMTHREADS;i++){
	pthread_join(threads[i],NULL);
    }
    
    int idx1=0,idx2=0,setsize=thread_data[0].vcf_list.size();
    string samname;
    for(int i=0;i<setsize;i++){ 
  	for(int j=0;j<DEPTH;j++)idx2=thread_data[0].vcf_list[i].find_last_of("/");
        idx1=thread_data[0].vcf_list[i].find_last_of("/",idx2-1)+1;  
        samname = thread_data[0].vcf_list[i].substr(idx1,idx2-idx1); 
        samples.insert(pair<string,int>(samname,gid_map.find(samname)->second));
    }

    time(&end);
    gid_map.clear();
    printf("Union Time: %.f sec\n",difftime(end,start));
    
    //PRINTING
    if(multioutput){
	for(int i=0;i<NUMTHREADS; i++){
    	    thread_data[i].samples = samples;
	    pthread_create(&threads[i],NULL,multiprint_2, (void*) &thread_data[i]);
    	}
        for(int i=0;i<NUMTHREADS;i++){
	    pthread_join(threads[i],NULL);
    	}
        FILE *output1,*output2,*output3;
	string temp=outfile+"_SNPPOS_vs_VARcount.txt"; 
	output1=fopen(temp.c_str(),"w"); 
	temp=outfile+"_POS_vs_DEPTH.txt"; 
	output2=fopen(temp.c_str(),"w"); 
        temp=outfile+"_POS_vs_QS.txt"; 
	output3=fopen(temp.c_str(),"w"); 
	for(int i=0;i<NUMTHREADS;i++){
	    //SNPpos vs VarietyCount
	    for(int j=0;j<SIZE;j++)
		if(thread_data[i].count[j]!=0) tally1[thread_data[i].count[j]]++;
	    //POS vs Depth
	    for(int j=0;j<SAMPLECOUNT;j++){ 
		tally2[j]+=thread_data[i].depth[j];
		tally3[j]+=thread_data[i].qs[j];
	    }
	}
        fprintf(output1,"#ofVar\t#ofSNPpos\n");
        fprintf(output2,"Depth\t#ofCalls\n");
	fprintf(output3,"QS\t#ofCalls\n");
	for(int i=0;i<SAMPLECOUNT;i++){ 
	    if(tally1[i]) fprintf(output1,"%d\t%.0f\n",i,tally1[i]);
	    if(tally2[i]) fprintf(output2,"%d\t%.0f\n",i,tally2[i]);
            if(tally3[i]) fprintf(output3,"%d\t%.0f\n",i*5,tally3[i]); //i=QS/5
	}
	free(tally1);
        free(tally2);
        free(tally3);
	fclose(output1);
	fclose(output2);
        fclose(output3);
    }else{
	FILE *output; 
        outfile+="_SNP_AD.txt"; 
        output=fopen(outfile.c_str(),"w"); 
        for(int i=0;i<NUMTHREADS;i++){
       	    if(!compressed) printSNPlist_1(i,&thread_data[i],samples,output);
    	    else printSNPlist_2(i,&thread_data[i],samples,output);
    	}
    	fclose(output); 
    }

    time(&end);
    chrmap.clear();
    samples.clear();
    for(int i=0;i<NUMTHREADS;i++){
    	if(thread_data[i].snp!=NULL) free(thread_data[i].snp);
	if(thread_data[i].ref!=NULL) free(thread_data[i].ref);
	if(thread_data[i].count!=NULL) free(thread_data[i].count);
	if(thread_data[i].depth!=NULL) free(thread_data[i].depth);
    }
    

    printf("Printing Time: %.f sec\n",difftime(end,start));
    
    return 0;
}
