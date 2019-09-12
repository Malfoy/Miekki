#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <math.h>
#include <unistd.h>
#include <bitset>
#include "Miekki.h"
#include "utils.h"



using namespace std;
using namespace chrono;







inline uint number_miss(const string str1,const string str2){
	uint res(0);
	for(uint i(0);i<str1.size();++i){
		if(str1[i]!=str2[i]){
			res++;
		}
	}
	return res;
}




inline char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



inline  string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



inline string getCanonical(const string& str){
	return (min(str,revComp(str)));
}



uint collisions_unordered(vector<uint64_t>& V1, const vector<uint64_t>& V2){
	uint res(0);
	sort(V1.begin(),V1.end());
	for(uint i(0);i<V2.size();++i){
		if(binary_search (V1.begin(), V1.end(), V2[i])){
			++res;
		}
	}
	return res;
}


uint collisions_sort(vector<uint64_t>& V1,  vector<uint64_t>& V2){
	uint res(0);
	sort(V1.begin(),V1.end());
	V1.erase( unique( V1.begin(), V1.end() ), V1.end() );
	sort(V2.begin(),V2.end());
	V2.erase( unique( V2.begin(), V2.end() ), V2.end() );
	for(uint i1(0),i2(0);i1<V1.size() and i2<V2.size();){
		if(V1[i1]==V2[i2]){
			++res;++i2;++i1;
		}else if(V1[i1]<V2[i2]){
			++i1;
		}else{
			++i2;
		}
	}
	return res;
}



void help(){
	cout<<"This is a help message"<<endl;
	cout<<"Input "<<endl;
	cout<<"-i load a constructed index from disk"<<endl;
	cout<<"-l construct an index from a list of file"<<endl;
	cout<<"-a query a fasta file"<<endl;
	cout<<"\nOutput "<<endl;
	cout<<"-o output file name (out.txt)"<<endl;
	cout<<"-d dump the index on disk"<<endl;
	cout<<"\nPerformances "<<endl;
	cout<<"-h use 2^h minimizers per sequence (20  for 1048576 minimizers)"<<endl;
	cout<<"-k kmer size (31)"<<endl;
	cout<<"-s minimal estimated intersection to be reported (100)"<<endl;
	cout<<"-t thread number (1)"<<endl;
	cout<<"\nAdvanced usage "<<endl;
	cout<<"-f fingerprint size "<<endl;
	cout<<"-b 2^b bits used for the bloom filter  "<<endl;
	cout<<"-e exact mode, real intersection will be computed on hits"<<endl;
}



int main(int argc, char ** argv){
	if(argc<2){
		help();
		exit(0);
	}
	string index_file,index_file_of_file,query_fastq,query_lines_of_file,query_files_of_list,outputFile("out.txt"),index_dump("");
	uint64_t H(17),core_number(8),kmer_size(31),bloom_size(33),fingerprint_size(3);
	double threshold(200);
	bool exact_mode(false);
	srand (time(NULL));
	Miekki *INDEX;
	char c;
	while ((c = getopt (argc, argv, "i:l:a:h:t:f:k:s:b:o:ed:A:")) != -1){
		switch(c){
			case 'i':
				index_file=optarg;
				break;
			case 'l':
				index_file_of_file=optarg;
				break;
			case 'a':
				query_lines_of_file=optarg;
				break;
			case 'A':
				query_files_of_list=optarg;
				break;
			//~ case 'q':
				//~ query_fastq=optarg;
				//~ break;
			case 'o':
				outputFile=optarg;
				break;
			case 'h':
				H=stoi(optarg);
				break;
			case 't':
				core_number=stoi(optarg);
				break;
			case 'k':
				kmer_size=stoi(optarg);
				break;
			case 's':
				threshold=stof(optarg);
				break;
			case 'f':
				fingerprint_size=stoi(optarg);
				break;
			case 'b':
				bloom_size=stoi(optarg);
				break;
			case 'e':
				exact_mode=true;
				break;
			case 'd':
				index_dump=optarg;
				break;
		}
	}

	uint32_t bit_per_min=(5+fingerprint_size);

	cout<<"Using "<<bit_per_min<<" bits per minimizer, "<<intToString(1<<H)<<" minimizers so "<<intToString(bit_per_min*(1<<H))<<" bits per sequences"<<endl;
	auto start = system_clock::now();
	if(index_file!=""){
		INDEX =new Miekki(index_file);
		INDEX->out= new ofstream(outputFile.c_str());
		INDEX->core_number=core_number;
		cout<<"I output results in "<<outputFile<<endl;
		cout<<"Load sucessful"<<endl;
		//TODO CHANGE WHAT CAN BE CHANGED
	}else if(index_file_of_file!=""){
		INDEX =new Miekki(kmer_size,H,bit_per_min,5,0,outputFile,bloom_size,threshold,core_number);
		INDEX->index_file_of_file(index_file_of_file);
		INDEX->compress_index(1);
	}else{
		cout<<"What am I supposed to index ? use either -i or -l options please"<<endl;
		help();
		exit(0);
	}
	if(index_dump!=""){
		cout<<"I write this index on the disk for later"<<endl;
		INDEX->dump_disk(index_dump);
	}

	auto endIndex = system_clock::now();
	duration<double> elapsed_seconds = endIndex - start;
	cout<< "elapsed time: " << elapsed_seconds.count() << "s\n";
	if(query_lines_of_file!=""){
		if(exact_mode){
			cout<<"running in exact mode, actual intersection will be computed on hits found by the index"<<endl;
			INDEX->query_file_exact(query_lines_of_file);
			//TODO HANDLE FASTQ FILES
		}else{
			cout<<"running in approx mode, intersection is estimated by the index"<<endl;
			INDEX->query_file(query_lines_of_file);
		}
	}else if(query_files_of_list!=""){
		if(exact_mode){
			cout<<"running in exact mode, actual intersection will be computed on hits found by the index"<<endl;
			INDEX->query_file_of_file_exact(query_files_of_list);
		}else{
			INDEX->query_file_of_file(query_files_of_list);
		}
	}else{
		cout<<"No query file, No queries"<<endl;
	}
	auto endQuery = system_clock::now();
	elapsed_seconds = endQuery -endIndex;
	cout<< "elapsed time: " << elapsed_seconds.count() << "s\n";
	cout<<"The end"<<endl;

	return 0;
}
