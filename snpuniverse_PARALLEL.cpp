/*
 *Last Update: Feb. 25, 2014
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

using namespace std;

struct threadData{
    string dirpath;
    vector<string> vcf_list;
    vector<string> chrom;
    unordered_map<string,int> chr;
    char* snp;
    int tid;
};

static const char *options="p:P:i:I:x:X:c";
static string outfile;
static string path;
static string input;
bool compressed=false;
int SIZE=45000000;
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
	default: break;
        } // switch
    } // while
} // parseArgs

int checkAlt(char ref,string alt,void *thread_data,int snppos){
    threadData *t = (threadData*) thread_data;
    int temp1=0,temp2=1;
    if(!alt.compare(".")){ return 0;} //homozygous
    else if(alt.size()==1 && alt[0]!=ref){ 
	if(t->snp[snppos]!='B') t->snp[snppos]='A'; return 1;
    }

    if(alt.size()>1){ 
	while(temp1<temp2){
	    temp2=alt.find_first_of(",",temp1+1);
	    if(temp2-temp1>2 || temp2==alt.npos){ //indels,structural variants
		t->snp[snppos]='B'; 
		return 0; 
	    }
	    temp1=temp2+1;
	}
    }
    if(t->snp[snppos]!='B') t->snp[snppos]='A';
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
    string temp,exten,filename,curpath;
    struct dirent *ep;
    time_t start,end;
    int dotpos;
    
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
		    if(!exten.compare(".vcf") && !compressed && filename.find(t->chrom[t->tid])!=filename.npos){
		        printf("thread%d -  %s\t",t->tid,filename.c_str());
			t->vcf_list.push_back(temp);
			locateSNP_1(temp,thread_data);
			time(&end);
            		printf("Time: %.f sec\n",difftime(end,start));
		    }else if(!exten.compare(".gz") && compressed && filename.find(t->chrom[t->tid])!=filename.npos){
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

int printSNPlist_1(int tid,void *thread_data,unordered_map<string,int> &samples,FILE *output){
    threadData *t = (threadData*) thread_data;
    string linestream,alt,temp1,formatfield,formatval,temp2,samname;
    char ref,*token=NULL,tok_ar[40];
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0,chrpos=0,setsize=t->vcf_list.size();
    int AD[128];
    if(tid==0){
    	for(int i=0;i<setsize;i++){ 
  	    for(int j=0;j<DEPTH;j++)idx2=t->vcf_list[i].find_last_of("/");
            idx1=t->vcf_list[i].find_last_of("/",idx2-1)+1;  
            samname = t->vcf_list[i].substr(idx1,idx2-idx1); 
            samples.insert(pair<string,int>(samname,samples.size()));
	    fprintf(output,">%lu %s\n",samples.size()-1,samname.c_str()); //sample list
    	}
   
    	fprintf(output,"SampleID\tChrom\tPos\tA\tT\tC\tG\n");//header
    }
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

        if(i==0){ //count SNPs
	    int snpcount=0;
	    for(int x=0;x<SIZE;x++){
		if(t->snp[x]=='A') snpcount++;
	    }
	    printf("Total of %d SNPs found.\n",snpcount);
	}
    
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

int printSNPlist_2(int tid,void *thread_data,unordered_map<string,int> &samples,FILE *output){
    threadData *t = (threadData*) thread_data;
    string linestream,alt,temp1,formatfield,formatval,temp2,samname;
    char ref,*token=NULL,tok_ar[40];
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0,chrpos=0,setsize=t->vcf_list.size();
    int AD[128];
    if(tid==0){
    	for(int i=0;i<setsize;i++){ 
  	    for(int j=0;j<DEPTH;j++)idx2=t->vcf_list[i].find_last_of("/");
            idx1=t->vcf_list[i].find_last_of("/",idx2-1)+1;  
            samname = t->vcf_list[i].substr(idx1,idx2-idx1); 
            samples.insert(pair<string,int>(samname,samples.size()));
	    fprintf(output,">%lu %s\n",samples.size()-1,samname.c_str()); //sample list
    	}
   
    	fprintf(output,"SampleID\tChrom\tPos\tA\tT\tC\tG\n");//header
    }
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

        if(i==0){ //count SNPs
	    int snpcount=0;
	    for(int x=0;x<SIZE;x++){
		if(t->snp[x]=='A') snpcount++;
	    }
	    printf("Total of %d SNPs found.\n",snpcount);
	}
    
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

int main(int argc, char **argv)
{
    const rlim_t STACK_SIZE = 1000*1024*1024*2; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);
    time_t start, end;
    pthread_t  threads[NUMTHREADS];
    threadData thread_data[NUMTHREADS];
    unordered_map<string,int> chrmap,samples;
    vector<string> chrnames; 

    time(&start);
    parseArgs(argc,argv);
    if (path.empty() || input.empty() || outfile.empty()) {
        cerr << "Usage:\n" << *argv
                  << " -p path to directory\n"
		  << " -i directory name\n"
		  << " -x output file\n"
		  << " -c compressed VCFs (vcf.gz)\n"
                  << endl;
        return 1;
    }

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
    	thread_data[i].tid = i;
	pthread_create(&threads[i],NULL,readFolder, (void*) &thread_data[i]);
    }

    for(int i=0;i<NUMTHREADS;i++){
	pthread_join(threads[i],NULL);
    }

    FILE *output; 
    outfile+="_SNP_AD.txt"; 
    output=fopen(outfile.c_str(),"w"); 
    
    for(int i=0;i<NUMTHREADS;i++){
       	if(!compressed) printSNPlist_1(i,&thread_data[i],samples,output);
    	else printSNPlist_2(i,&thread_data[i],samples,output);
    }
    fclose(output);
    time(&end);
    chrmap.clear();
    for(int i=0;i<NUMTHREADS;i++)
    	if(thread_data[i].snp!=NULL) free(thread_data[i].snp);
    

    printf("Time: %.f sec\n",difftime(end,start));
    return 0;
}
