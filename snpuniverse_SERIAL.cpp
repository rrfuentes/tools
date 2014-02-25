/*
 *Last Update: Feb. 25, 2014
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
 *SNP_universe with serial process 
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

using namespace std;

static const char *options="p:P:i:I:x:X:c";
static string outfile;
static string path;
static string input;
bool compressed=false;
int SIZE1=14,SIZE2=30000000;

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

int checkAlt(char ref,string alt,char **&snp,int chrpos,int snppos){
    int temp1=0,temp2=1;
    if(!alt.compare(".")){ return 0;} //homozygous
    else if(alt.size()==1 && alt[0]!=ref){ 
	if(snp[chrpos][snppos]!='B') snp[chrpos][snppos]='A'; return 1;
    }

    if(alt.size()>1){ 
	while(temp1<temp2){
	    temp2=alt.find_first_of(",",temp1+1);
	    if(temp2-temp1>2 || temp2==alt.npos){ //indels,structural variants
		snp[chrpos][snppos]='B'; 
		return 0; 
	    }
	    temp1=temp2+1;
	}
    }
    if(snp[chrpos][snppos]!='B') snp[chrpos][snppos]='A';
    return 1; //multiple ALTs
}

int locateSNP_1(string filepath,unordered_map<string,int> chr,char **&snp){
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
            chrpos = chr.find(temp)->second; //printf("%s %d %d\n",temp.c_str(),chrpos,snppos);
	    idx1 = linestream.find_first_of("\t",idx2+1); //skip id column
    	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
    
	    if(idx2-idx1==2){ //skip indels and structural variants
		ref=linestream[idx2-1]; 
		alt=linestream.substr(idx2+1,linestream.find_first_of("\t",idx2+1)-idx2-1); 
	        //printf("%d %c %s %d\t",snppos,ref,alt.c_str(),snp[chrpos][snppos]);
		checkAlt(ref,alt,snp,chrpos,snppos); //check if SNP occur in a position
                //printf("%d\n",snp[chrpos][snppos]);
                
	    }else{
		snp[chrpos][snppos]='B'; //IGNORE INDEL(deletion)
	    }
	}
	idx1=0;
    }
    fp.close();
    return 0;
}

int locateSNP_2(string filepath,unordered_map<string,int> chr,char **&snp){
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
            chrpos = chr.find(temp)->second; //printf("%s %d %d\n",temp.c_str(),chrpos,snppos);
	    idx1 = linestream.find_first_of("\t",idx2+1); //skip id column
    	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
    
	    if(idx2-idx1==2){ //skip indels and structural variants
		ref=linestream[idx2-1]; 
		alt=linestream.substr(idx2+1,linestream.find_first_of("\t",idx2+1)-idx2-1); 
	        //printf("%d %c %s %d\t",snppos,ref,alt.c_str(),snp[chrpos][snppos]);
		checkAlt(ref,alt,snp,chrpos,snppos); //check if SNP occur in a position
                //printf("%d\n",snp[chrpos][snppos]);
                
	    }else{
		snp[chrpos][snppos]='B'; //IGNORE INDEL(deletion)
	    }
	}
	idx1=0;
    }
    fp.close();
    return 0;
}

int readFolder(string path,unordered_map<string,int> chr,char **&snp,vector<string> &vcf_list){
    DIR *dp;
    string temp,exten,filename;
    struct dirent *ep;
    time_t start,end;
    int dotpos;
    
    dp = opendir(path.c_str());
    if(path[path.size()-1]!='/') path += '/';
    struct stat filestat;
    if(dp!=NULL){
	while(ep = readdir(dp)){
            time(&start);
	    if(!strcmp(ep->d_name,".") || !strcmp(ep->d_name,"..")) continue;
	    temp = path + ep->d_name;
            stat(temp.c_str(), &filestat);
	    filename=ep->d_name;
	    if(S_ISDIR(filestat.st_mode)){ //recursively read directories
		readFolder(temp,chr,snp,vcf_list);
	    }else if(S_ISREG(filestat.st_mode)){
		//check if vcf
		dotpos=filename.find_last_of(".");
		if(dotpos!=filename.npos){
		    exten=filename.substr(dotpos);
		    if(!exten.compare(".vcf")){
		        printf("%s\t",filename.c_str());
			vcf_list.push_back(temp);
			locateSNP_1(temp,chr,snp);
			time(&end);
            		printf("Time: %.f sec\n",difftime(end,start));
		    }else if(!exten.compare(".gz") && compressed){
			exten=filename.substr(filename.find_last_of(".",dotpos-1)); 
		        if(!exten.compare(".vcf.gz")){
			    printf("%s\t",filename.c_str());
			    vcf_list.push_back(temp);
			    locateSNP_2(temp,chr,snp);
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
	return 1;
    }
    return 0;
}

int printSNPlist_1(unordered_map<string,int> chr,char **snp,vector<string> vcf_list,string outfile){
    FILE *output; 
    outfile+="_SNP_AD.txt"; 
    output=fopen(outfile.c_str(),"w"); 
    string linestream,alt,temp1,formatfield,formatval,temp2;
    char ref,*token=NULL,tok_ar[40];
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0,chrpos=0,setsize=vcf_list.size();
    int AD[128];
    for(int i=0;i<setsize;i++) fprintf(output,">%d %s\n",i,vcf_list[i].c_str()); //sample list
    fprintf(output,"SampleID\tChrom\tPos\tA\tT\tC\tG\n");//header
    for(int i=0;i<setsize;i++){
    	ifstream fp(vcf_list[i].c_str());
    	if (!fp.is_open()) { //check input file
	    printf("ERROR: Failed to open the input file %s", vcf_list[i].c_str());
	    return 1;
    	}
        
        if(i==0){ //count SNPs
	    int snpcount=0;
	    for(int x=0;x<SIZE1;x++){
		for(int y=0;y<SIZE2;y++){
		    if(snp[x][y]=='A') snpcount++;
		}
	    }
	    printf("Total of %d SNPs found.\n",snpcount);
	}
    
    	for(int x=0;getline(fp,linestream);x++){
	    if(linestream[0]!='#'){
                idx1 = linestream.find_first_of("\t",idx1+1); //first column
             	temp1 = linestream.substr(0,idx1); //get the chrom number
            	idx2 = linestream.find_first_of("\t",++idx1); //second column
            	snppos = atoi(linestream.substr(idx1,idx2-idx1).c_str()); //get the SNP pos
                chrpos = chr.find(temp1)->second;
              
            	if(snp[chrpos][snppos]=='A'){
		    fprintf(output,"%d\t%s\t%d\t",i,temp1.c_str(),snppos);
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
			    if((!checkAlt(ref,alt,snp,chrpos,snppos) && alt[0]!='.') || ref=='.'){
				fprintf(output,"0\t0\t0\t0 (Indel/Structural Variant)\n"); printf("Indel/SV:Sample %d %s %c\n",i,temp1.c_str(),snp[chrpos][snppos]);
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
    fclose(output);
    return 0;
}


int printSNPlist_2(unordered_map<string,int> chr,char **snp,vector<string> vcf_list,string outfile){
    FILE *output; 
    outfile+="_SNP_AD.txt"; 
    output=fopen(outfile.c_str(),"w"); 
    string linestream,alt,temp1,formatfield,formatval,temp2;
    char ref,*token=NULL,tok_ar[40];
    int idx1=0,idx2=0,idx3=0,idx4=0,snppos=0,chrpos=0,setsize=vcf_list.size();
    int AD[128];
    for(int i=0;i<setsize;i++) fprintf(output,">%d %s\n",i,vcf_list[i].c_str()); //sample list
    fprintf(output,"SampleID\tChrom\tPos\tA\tT\tC\tG\n");//header
    for(int i=0;i<setsize;i++){
    	igzstream fp(vcf_list[i].c_str());
    	if (!fp.good()) { //check input file
	    printf("ERROR: Failed to open the input file %s", vcf_list[i].c_str());
	    return 1;
    	}

        if(i==0){ //count SNPs
	    int snpcount=0;
	    for(int x=0;x<SIZE1;x++){
		for(int y=0;y<SIZE2;y++){
		    if(snp[x][y]=='A') snpcount++;
		}
	    }
	    printf("Total of %d SNPs found.\n",snpcount);
	}
    
    	for(int x=0;getline(fp,linestream);x++){
	    if(linestream[0]!='#'){
                idx1 = linestream.find_first_of("\t",idx1+1); //first column
             	temp1 = linestream.substr(0,idx1); //get the chrom number
            	idx2 = linestream.find_first_of("\t",++idx1); //second column
            	snppos = atoi(linestream.substr(idx1,idx2-idx1).c_str()); //get the SNP pos
                chrpos = chr.find(temp1)->second;
              
            	if(snp[chrpos][snppos]=='A'){
		    fprintf(output,"%d\t%s\t%d\t",i,temp1.c_str(),snppos);
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
			temp2=formatfield.substr(idx1,idx2-idx1); 
			idx4 = formatval.find_first_of(":\0",++idx3);
			if(!temp2.compare("AD")){
			    if((!checkAlt(ref,alt,snp,chrpos,snppos) && alt[0]!='.') || ref=='.'){
				fprintf(output,"0\t0\t0\t0 (Indel/Structural Variant)\n"); printf("Indel/SV:Sample %d %s %c\n",i,temp1.c_str(),snp[chrpos][snppos]);
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
	    }
	    idx1=0;
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
    time_t start, end;
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
   
    DIR *dp;
    string curpath,temp;
    vector<string> vcf_list;
    char **snp;
    unordered_map<string,int> chr;
    struct dirent *ep;
    if(path[path.size()-1]=='/'){
	curpath = path + input;
    }else{
	curpath = path + "/" + input;
    }
    snp=(char**)calloc(SIZE1,sizeof(char*));
    for(int i=0;i<SIZE1;i++) snp[i]=(char*)calloc(SIZE2,sizeof(char));
    
    chr["Chr1"]=0;
    chr["Chr2"]=1;
    chr["Chr3"]=2;
    chr["Chr4"]=3;
    chr["Chr5"]=4;
    chr["Chr6"]=5;
    chr["Chr7"]=6;
    chr["Chr8"]=7;
    chr["Chr9"]=8;
    chr["Chr10"]=9;
    chr["Chr11"]=10;
    chr["Chr12"]=11;
    chr["ChrSy"]=12;
    chr["ChrUn"]=13;

    dp = opendir(curpath.c_str());
    if(curpath[curpath.size()-1]!='/') curpath += '/';
    
    struct stat filestat;
    if(dp!=NULL){
	while(ep = readdir(dp)){
	    if(!strcmp(ep->d_name,".") || !strcmp(ep->d_name,"..")) continue;
	    temp = curpath + ep->d_name;
	    stat(temp.c_str(), &filestat);
	    if(S_ISDIR(filestat.st_mode)){
		readFolder(temp,chr,snp,vcf_list);
	    }
	    //printf("%s\n",ep->d_name);
	    temp.clear();
	}
	printf("Printing result.\n");
        closedir (dp); 
	if(vcf_list.size() && !compressed) printSNPlist_1(chr,snp,vcf_list,outfile);
	if(vcf_list.size() && compressed) printSNPlist_2(chr,snp,vcf_list,outfile);
    }else{
	printf("ERROR: Can't find the directory.");
    }
    chr.clear();
    for(int i=0;i<SIZE1;i++)
    	if(snp[i]!=NULL) free(snp[i]);
    free(snp);
    time(&end);
    printf("Time: %.f sec\n",difftime(end,start));
    return 0;
}
