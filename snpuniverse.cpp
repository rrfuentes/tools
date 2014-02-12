#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sys/resource.h>
#include <map>
#include <vector>
#include <gzstream.h>
#include <pthread.h>
using namespace std;

/*
http://pubs.opengroup.org/onlinepubs/009695299/basedefs/sys/stat.h.html
S_ISBLK(m)  Test for a block special file.
S_ISCHR(m)  Test for a character special file.
S_ISDIR(m)  Test for a directory.
S_ISFIFO(m) Test for a pipe or FIFO special file.
S_ISREG(m)  Test for a regular file.
S_ISLNK(m)  Test for a symbolic link.
S_ISSOCK(m) Test for a socket.
*/

static const char *options="p:P:i:I:x:X:cm:";
static string outfile;
static string path;
static string input;
bool compressed=false;
int multicore=0;
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
        case 'm': multicore=atoi(optarg); break; //utilizes pthreads
	default: break;
        } // switch
    } // while
} // parseArgs

int checkAlt(char ref,string alt,short *&count,int snppos){
    int temp1=0,temp2=1;
    if(!alt.compare(".")){ return 0;} //homozygous
    else if(alt.size()==1 && alt[0]!=ref) return 1;

    if(alt.size()>1){ 
	while(temp1<temp2){
	    temp2=alt.find_first_of(",",temp1+1);
	    if(temp2-temp1>2 || temp2==alt.npos){
		count[snppos]-=2; //IGNORE Indel(insertion/sv) (Remove line to report the pos if SNP is present in other sample(also del filter in locateSNP()); -2 in case the position is represented twice)
		return 0; //indels,structural variants
	    }
	    temp1=temp2+1;
	}
    }
    return 1; //multiple ALTs
}

int locateSNP_1(string filepath,map<string,int> &pos,bool *&snp,short *&count){
    string linestream,temp,alt;
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0;
    char ref;
    int inisize=pos.size();
    
    ifstream fp(filepath.c_str());
    pair<map<string,int>::iterator,bool> ret;
    map<string,int>::iterator it;
    if (!fp.is_open()) { //check input file
	printf("ERROR: Failed to open the input file %s", filepath.c_str());
	return 1;
    }
    
    for(int x=0;getline(fp,linestream);x++){
	if(linestream[0]!='#'){
	    idx1 = linestream.find_first_of("\t",idx1+1); //first column
            idx1 = linestream.find_first_of("\t",idx1+1); //second column
            temp=linestream.substr(0,idx1); //get the chr+pos concatenated string
	    if(inisize==0){ //create map when reading first vcf file
		snppos=pos.size();
	        ret = pos.insert(pair<string,int>(temp,snppos)); 
		if(ret.second!=false){ //avoids duplicate positions
		    count[snppos]++;
		}else{
 		    count[ret.first->second]-=2; //IGNORE INDEL(duplicated positions)
		    snppos=-1;
		}
	    }else{
	        it=pos.find(temp);
	        if(it!=pos.end()){
		    snppos=it->second;
		    count[snppos]++; //record if the position is map in all samples
		}else snppos=-1;
	        //printf("%d",snppos);
	    }
	    idx1 = linestream.find_first_of("\t",idx1+1); //skip id column
    	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
	    if(idx2-idx1==2 && snppos!=-1){ //skip indels and structural variants
		ref=linestream[idx2-1]; 
		alt=linestream.substr(idx2+1,linestream.find_first_of("\t",idx2+1)-idx2-1);
	        //printf("%s %d %c %s %d\t",temp.c_str(),snppos,ref,alt.c_str(),snp[snppos]);
		snp[snppos] = (checkAlt(ref,alt,count,snppos) || snp[snppos]); //check if SNP occur in a position
                //printf("%d\n",snp[snppos]);
	    }else if(snppos!=-1){
		count[snppos]-=2; //IGNORE INDEL(deletion)
	    }
	}
	idx1=idx2=0;
    }
    fp.close();
    return 0;
}

int locateSNP_2(string filepath,map<string,int> &pos,bool *&snp,short *&count){
    string linestream,temp,alt;
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0;
    char ref;
    int inisize=pos.size();
    
    igzstream fp(filepath.c_str());
    pair<map<string,int>::iterator,bool> ret;
    map<string,int>::iterator it;
    if (!fp.good()) { //check input file
	printf("ERROR: Failed to open the input file %s", filepath.c_str());
	return 1;
    }
    
    for(int x=0;getline(fp,linestream);x++){
	if(linestream[0]!='#'){
	    idx1 = linestream.find_first_of("\t",idx1+1); //first column
            idx1 = linestream.find_first_of("\t",idx1+1); //second column
            temp=linestream.substr(0,idx1); //get the chr+pos concatenated string
	    if(inisize==0){ //create map when reading first vcf file
		snppos=pos.size();
	        ret = pos.insert(pair<string,int>(temp,snppos)); 
		if(ret.second!=false){ //avoids duplicate positions
		    count[snppos]++;
		}else{
 		    count[ret.first->second]-=2; //IGNORE INDEL(duplicated positions)
		    snppos=-1;
		}
	    }else{
	        it=pos.find(temp);
	        if(it!=pos.end()){
		    snppos=it->second;
		    count[snppos]++; //record if the position is map in all samples
		}else snppos=-1;
	        //printf("%d",snppos);
	    }
	    idx1 = linestream.find_first_of("\t",idx1+1); //skip id column
    	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
	    if(idx2-idx1==2 && snppos!=-1){ //skip indels and structural variants
		ref=linestream[idx2-1]; 
		alt=linestream.substr(idx2+1,linestream.find_first_of("\t",idx2+1)-idx2-1);
	        //printf("%s %d %c %s %d\t",temp.c_str(),snppos,ref,alt.c_str(),snp[snppos]);
		snp[snppos] = (checkAlt(ref,alt,count,snppos) || snp[snppos]); //check if SNP occur in a position
                //printf("%d\n",snp[snppos]);
	    }else if(snppos!=-1){
		count[snppos]-=2; //IGNORE INDEL(deletion)
	    }
	}
	idx1=idx2=0;
    }
    fp.close();
    return 0;
}

int readFolder(string path,map<string,int> &pos,bool *&snp,short *&count,vector<string> &vcf_list){
    DIR *dp;
    string temp,exten,filename;
    struct dirent *ep;
    int dotpos;
    dp = opendir(path.c_str());
    if(path[path.size()-1]!='/') path += '/';
    struct stat filestat;
    if(dp!=NULL){
	while(ep = readdir(dp)){
	    if(!strcmp(ep->d_name,".") || !strcmp(ep->d_name,"..")) continue;
	    temp = path + ep->d_name;
            stat(temp.c_str(), &filestat);
	    filename=ep->d_name;
	    if(S_ISDIR(filestat.st_mode)){ //recursively read directories
		readFolder(temp,pos,snp,count,vcf_list);
	    }else if(S_ISREG(filestat.st_mode)){
		//check if vcf
		dotpos=filename.find_last_of(".");
		if(dotpos!=filename.npos){
		    exten=filename.substr(dotpos);
		    if(!exten.compare(".vcf")){
		        printf("%s\n",filename.c_str());
			vcf_list.push_back(temp);
			locateSNP_1(temp,pos,snp,count);
		    }else if(!exten.compare(".gz") && compressed){
			exten=filename.substr(filename.find_last_of(".",dotpos-1)); 
		        if(!exten.compare(".vcf.gz")){
			    printf("%s\n",filename.c_str());
			    vcf_list.push_back(temp);
			    locateSNP_2(temp,pos,snp,count);
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
	return 1;
    }
    return 0;
}

int printSNPlist_1(map<string,int> &pos,bool *&snp,short *&count,vector<string> vcf_list,string outfile){
    FILE *output; 
    outfile+="_SNP_AD.txt"; 
    output=fopen(outfile.c_str(),"w"); 
    string linestream,alt,temp1,formatfield,formatval,temp2;
    char ref,*token=NULL,tok_ar[40];
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0,setsize=vcf_list.size();
    int AD[128];
    map<string,int>::iterator it;
    for(int i=0;i<setsize;i++) fprintf(output,">%d %s\n",i,vcf_list[i].c_str()); //sample list
    fprintf(output,"SampleID\tChrom\tPos\tA\tT\tC\tG\n");//header
    for(int i=0;i<setsize;i++){
    	ifstream fp(vcf_list[i].c_str());
    	if (!fp.is_open()) { //check input file
	    printf("ERROR: Failed to open the input file %s", vcf_list[i].c_str());
	    return 1;
    	}
    
    	for(int x=0;getline(fp,linestream);x++){
	    if(linestream[0]!='#'){
            	idx1 = linestream.find_first_of("\t",idx1+1); //CHROM
            	idx1 = linestream.find_first_of("\t",idx1+1); //POS
            	temp1=linestream.substr(0,idx1); //get the chr+pos concatenated string
		it=pos.find(temp1);
		if(it!=pos.end()) snppos=it->second;
		else snppos=-1;
               
            	if(snppos!=-1 && (count[snppos]==setsize && snp[snppos])){
		    fprintf(output,"%d\t%s\t",i,it->first.c_str());
                    idx1 = linestream.find_first_of("\t",idx1+1); //skip ID
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
		    idx1=0;
	            while(idx1!=formatfield.npos){   
			idx2 = formatfield.find_first_of(":\0",++idx1);
			temp2=formatfield.substr(idx1,idx2-idx1); 
			idx4 = formatval.find_first_of(":\0",++idx3);
			if(!temp2.compare("AD")){
			    if(!checkAlt(ref,alt,count,snppos) || ref=='.'){
				fprintf(output,"0\t0\t0\t0 (Indel/Structural Variant)\n"); printf("Indel/SV:Sample %d %s %d\n",i,temp1.c_str(),count[snppos]);
			    }else{
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
	    }
	    idx1=idx3=0;
    	}
	fp.close();
    }
    fclose(output);
    return 0;
}

int printSNPlist_2(map<string,int> &pos,bool *&snp,short *&count,vector<string> vcf_list,string outfile){
    FILE *output; 
    outfile+="_SNP_AD.txt"; 
    output=fopen(outfile.c_str(),"w"); 
    string linestream,alt,temp1,formatfield,formatval,temp2;
    char ref,*token=NULL,tok_ar[40];
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0,setsize=vcf_list.size();
    int AD[128];
    map<string,int>::iterator it;
    for(int i=0;i<setsize;i++) fprintf(output,">%d %s\n",i,vcf_list[i].c_str()); //sample list
    fprintf(output,"SampleID\tChrom\tPos\tA\tT\tC\tG\n");//header
    for(int i=0;i<setsize;i++){
    	igzstream fp(vcf_list[i].c_str());
    	if (!fp.good()) { //check input file
	    printf("ERROR: Failed to open the input file %s", vcf_list[i].c_str());
	    return 1;
    	}
    
    	for(int x=0;getline(fp,linestream);x++){
	    if(linestream[0]!='#'){
            	idx1 = linestream.find_first_of("\t",idx1+1); //CHROM
            	idx1 = linestream.find_first_of("\t",idx1+1); //POS
            	temp1=linestream.substr(0,idx1); //get the chr+pos concatenated string
		it=pos.find(temp1);
		if(it!=pos.end()) snppos=it->second;
		else snppos=-1;
               
            	if(snppos!=-1 && (count[snppos]==setsize && snp[snppos])){
		    fprintf(output,"%d\t%s\t",i,it->first.c_str());
                    idx1 = linestream.find_first_of("\t",idx1+1); //skip ID
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
		    idx1=0;
	            while(idx1!=formatfield.npos){   
			idx2 = formatfield.find_first_of(":\0",++idx1);
			temp2=formatfield.substr(idx1,idx2-idx1); 
			idx4 = formatval.find_first_of(":\0",++idx3);
			if(!temp2.compare("AD")){
			    if(!checkAlt(ref,alt,count,snppos) || ref=='.'){
				fprintf(output,"0\t0\t0\t0 (Indel/Structural Variant)\n"); printf("Indel/SV:Sample %d %s %d\n",i,temp1.c_str(),count[snppos]);
			    }else{
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
	    }
	    idx1=idx3=0;
    	}
	fp.close();
    }
    fclose(output);
    return 0;
}

int main(int argc, char **argv)
{
    const rlim_t STACK_SIZE = 1000*1024*1024*2; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);
    parseArgs(argc,argv);
    if (path.empty() || input.empty() || outfile.empty()) {
        cerr << "Usage:\n" << *argv
                  << " -p path to directory\n"
		  << " -i directory name\n"
		  << " -x output file\n"
		  << " -c compressed VCFs (vcf.gz)\n"
		  << " -m multicore=#of_alloted_cores\n"
                  << endl;
        return 1;
    }
   
    DIR *dp;
    int SIZE=30000000,setsize=0;
    string curpath,temp;
    vector<string> vcf_list;
    bool *snp;
    map<string,int> pos;
    short *count;
    struct dirent *ep;
    if(path[path.size()-1]=='/'){
	curpath = path + input;
    }else{
	curpath = path + "/" + input;
    }
    snp=(bool*)calloc(SIZE,sizeof(bool));
    count=(short*)calloc(SIZE,sizeof(short));
    dp = opendir(curpath.c_str());
    if(curpath[curpath.size()-1]!='/') curpath += '/';
    
    //pthread_t threads[2];
    //threadData *thread_data = (threadData*)malloc(2*sizeof()threadData);
    struct stat filestat;
    if(dp!=NULL){
	while(ep = readdir(dp)){
	    if(!strcmp(ep->d_name,".") || !strcmp(ep->d_name,"..")) continue;
	    temp = curpath + ep->d_name;
	    stat(temp.c_str(), &filestat);
	    if(S_ISDIR(filestat.st_mode)){
		readFolder(temp,pos,snp,count,vcf_list);
	    }
	    //printf("%s\n",ep->d_name);
	    temp.clear();
	}
	
        closedir (dp); 
	if(vcf_list.size() && !compressed) printSNPlist_1(pos,snp,count,vcf_list,outfile);
	if(vcf_list.size() && compressed) printSNPlist_2(pos,snp,count,vcf_list,outfile);
    }else{
	printf("ERROR: Can't find the directory.");
    }
    pos.clear();
    free(snp);
    free(count);
    return 0;
}
