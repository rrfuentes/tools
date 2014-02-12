#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <map>
#include <vector>

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

static const char *options="p:P:i:I:d:D:x:X:y:Y:s:S:a:b:";
static string directable;
static string filestable;
static string path;
static string input;
static string diskpath;
static int set=-1;
static int offset1 = 0;
static int offset2 =0;
void parseArgs(int argc, char**argv){
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
	case 'p':
	case 'P': path = optarg; break; //complete path to the input
        case 'i':
        case 'I': input = optarg; break; //directory
	case 'd':
	case 'D': diskpath = optarg; break; //location of storage disk
	case 'x':
	case 'X': directable = optarg; break; //output FOLDER table
        case 'y':
        case 'Y': filestable = optarg; break; //output FILES table
	case 's':
        case 'S': set = atoi(optarg); break; //set number from list 
	case 'a': offset1 = atoi(optarg); //offset for adding new sample directories for existing list; next ID after the last tuple in the DB
	case 'b': offset2 = atoi(optarg); //offset for files(any files within a sample directory) contained in each sample directory
	default: break;
        } // switch
    } // while
} // parseArgs

void readFolder(FILE *&output2,string indv_f, int folder_c,int &file_c,map<string,int> gid_map,map<string, string> hash_map){
    DIR *dp;
    FILE *fp;
    string file,temp;
    char buf[150],hash[33];
   
    struct dirent *ep; 
    struct stat filestat;
    dp = opendir(indv_f.c_str()); 
    //dir_name  = indv_f.substr(idx, len-idx); cout << dir_name;
    if (dp != NULL)
    {
     	while (ep = readdir (dp)){ 
	    if(!strcmp(ep->d_name,".") || !strcmp(ep->d_name,"..")) continue;
	    temp = indv_f + ep->d_name;
	    stat(temp.c_str(), &filestat);
            if(S_ISREG(filestat.st_mode)){ //if(ep->d_type==DT_REG){
		/*file = "md5sum " + indv_f + ep->d_name;
                fp = popen(file.c_str(),"r");
		if(fp==NULL){
		    printf("Can't open a file(%s) to generate checksum.",file.c_str());
		}
	        fgets(buf,sizeof(buf),fp);  
		//setvbuf(fp,buf,_IOLBF,1024);
		//temp += buf;
		//strcpy(hash,temp.substr(0,32).c_str());
                strncpy(hash,buf,32); hash[32]='\0';*/
		fprintf(output2,"%d,%d,%s,%s\n",++file_c,folder_c,ep->d_name,hash_map.find(ep->d_name)->second.c_str());
                //pclose(fp);
	    }
	}

    	closedir (dp);
    }
}

int main(int argc, char **argv)
{
    parseArgs(argc,argv);
    if (path.empty() || input.empty() || diskpath.empty() || directable.empty() || filestable.empty() || set==-1) {
        cerr << "Usage:\n" << *argv
                  << " -p path to directory\n"
		  << " -i directory name\n"
		  << " -d storage path\n"
		  << " -x output for FOLDER table\n"
                  << " -y output for FILES table\n"
		  << " -s set number\n"
			//(0) 3K_IRRI_sample_info
		        //(1) CAAS-set_ZLi_PinYin,
                        //(2) IRRI_BGI_2nd_shipment_1300entries,
			//(3) P01-125-genomes-project)
		  << " -a directory offset for additional sample set\n"
                  << " -a file offset for additional sample set\n"
                  << endl;
        return 1;
    }

    if((offset1>0 && offset2==0) || (offset1==0 && offset2>0)) 
	printf("ERROR:Both offset for directory and file keys must be provided.");

    DIR *dp;
    FILE *output1, *output2;
    string curpath,filepath1,filepath2,linestream;
    string indv_f; //folder for individual
    map<string,int> gid_map;
    map<string,string> hash_map;
    int folder_c=0,file_c=0,pos1,pos2,gid; //counter
    pair<map<string,int>::iterator,bool> ret1; //map return 
    pair<map<string,string>::iterator,bool> ret2; //map return 
    if(path[path.size()-1]=='/'){
	curpath = path + input;
    }else{
	curpath = path + "/" + input;
    }
    filepath1 = directable + ".csv";
    filepath2 = filestable + ".csv";
    struct dirent *ep;     
    dp = opendir (curpath.c_str());
    if(curpath[curpath.size()-1]!='/') curpath += '/'; 
    //Open files
    output1 = fopen(filepath1.c_str(),"w");
    output2 = fopen(filepath2.c_str(),"w");

    //load list of GIDs
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
	gid = atoi(linestream.substr(pos2+1,linestream.find_first_of(",",pos2+1)-pos2).c_str());
	ret1 = gid_map.insert(pair<string,int>(linestream.substr(pos1,pos2-pos1),gid));
	//if(ret.second==false) printf("NOTE: row %d has duplicate BOX_POSITION_CODE%s.\n",x,linestream.substr(pos1,pos2-pos1).c_str());
	//cout << gid_map.find(linestream.substr(pos1,pos2-pos1))->second << " ";
    }
	
    //load checksum list
    if(set==0){
    	string checkad = curpath + "md5.txt"; 
    	ifstream check_s(checkad.c_str());
	if(!list.is_open()){
	    printf("ERROR: Cannot open the md5.txt.");
	    return 1;
    	}
	for(int x=0;getline(check_s,linestream)!=NULL;x++){
	    pos1 = linestream.find_first_of("/",36)+1; //NOTE: direc name is 3-4 char long
	    ret2 = hash_map.insert(pair<string,string>(linestream.substr(pos1,linestream.length()-pos1),linestream.substr(0,32)));
	    //cout << hash_map.find(linestream.substr(41,linestream.length()-41))->second <<"\n";
        }
    }
    if(offset1>0){
	folder_c=offset1-1;    
	file_c=offset2-1;
    }

    struct stat filestat;
    if (dp != NULL)
    {   
     	while (ep = readdir (dp)){ 
	    if(!strcmp(ep->d_name,".") || !strcmp(ep->d_name,"..")) continue;
	    indv_f = curpath + ep->d_name;
	    stat(indv_f.c_str(), &filestat);
            if(S_ISDIR(filestat.st_mode)){//if(ep->d_type==DT_DIR){ 
                folder_c++;
                if(indv_f[indv_f.size()-1]!='/') indv_f+='/';  
		fprintf(output1,"%d,%d,%s,%s,.\n",folder_c,gid_map.find(ep->d_name)->second,diskpath.c_str(),indv_f.c_str());
	       	if(set==0) readFolder(output2,indv_f,folder_c,file_c,gid_map,hash_map);
		//cout << indv_f<<"\n";
	    }else{
                //printf("ERROR: All files must be located in a folder of each respective individual.");
	    }
	    indv_f.clear();
	}

    	closedir (dp);
    }
    else
    	perror ("Couldn't open the directory");
    fclose(output1);
    fclose(output2);
    gid_map.clear();
    hash_map.clear();
    return 0;
}
