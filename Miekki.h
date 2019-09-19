#ifndef MIEKKI
#define MIEKKI



#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <unordered_set>
#include "utils.h"



using namespace std;



struct similarity_score{
	uint32_t sequence_identifier;
	uint32_t score;
	//~ uint32_t active_minimizer;
};



class Miekki{
public:
	uint32_t kmer_size;
	uint64_t offsetUpdatekmer;
	uint32_t number_minimizer;
	uint32_t number_minimizer_log2;
	uint32_t number_bit_minimizer;
	uint32_t number_bit_mantis;
	minimizer maximal_minimizer;
	uint32_t index_size;
	uint32_t bloom_size_log2;
	uint32_t number_hash;
	uint64_t bloom_size;
	uint32_t threshold;
	uint32_t core_number;

	bool Bloom_Filter_Protection;
	bool compressed;
	bool jaccard_estimation;
	bool containment_estimation;

	vector<string> index;
	vector<uint64_t> index_disk;
	vector<uint8_t> Bloom_Filter;
	//~ unordered_set<uint64_t> Bloom_Filter_Forever;
	vector<string> file_names;
	vector<uint32_t> sketch_size;
	vector<uint64_t> genome_size;
	omp_lock_t lock[10000];
	ofstream* out;
	vector<similarity_score> query_output;



	Miekki(uint32_t kmer_size_val, uint32_t number_minimizer_log2_val, uint32_t number_bit_minimizer_val, uint32_t number_bit_mantis_val,uint32_t output_mode,const string& output_File,uint32_t bloom_size_log2_val,uint32_t threshold_val, uint32_t core_number_val)

		: kmer_size(kmer_size_val),number_minimizer((uint32_t)1<<number_minimizer_log2_val), number_minimizer_log2(number_minimizer_log2_val), number_bit_minimizer(number_bit_minimizer_val), number_bit_mantis(number_bit_mantis_val), maximal_minimizer(-1),bloom_size_log2(bloom_size_log2_val),threshold(threshold_val),core_number(core_number_val),Bloom_Filter_Protection(true)
		{

		for(uint i(0);i<10000;++i){
			 omp_init_lock(&lock[i]);
		}
		out= new ofstream(output_File.c_str());
		cout<<"I output results in "<<output_File<<endl;
		offsetUpdatekmer=1;
		offsetUpdatekmer<<=2*kmer_size;
		containment_estimation=false;
		number_hash=5;
		index.assign(number_minimizer,{});
		index_size=0;
		bloom_size=0;
		compressed=false;

		if(Bloom_Filter_Protection){
			bloom_size=1;
			bloom_size<<=bloom_size_log2;
			Bloom_Filter.assign(bloom_size/8,0);
		}
	}



	Miekki(const string& input);



	//CORE FUNCTIONS
	void insert_sequence(const string& str, const string& title);
	void insert_sequences(const vector<pair<string,string>>& Vstr);
	vector<similarity_score> query_sequence(const string& str,uint32_t& lol);
	vector<vector<similarity_score>> query_sequences( vector<pair<string,uint32_t>>& batch);
	vector<vector<similarity_score>> query_sequences1(vector<pair<string,uint32_t>>& batch);
	vector<vector<similarity_score>> query_sequences2(vector<pair<string,uint32_t>>& batch);
	void query_file(const string& str);
	void query_whole_file(const string& str);
	void query_file_of_file(const string& str);
	void query_file_of_file_exact(const string& str);
	void query_file_exact(const string& str);
	void query_whole_file_exact(const string& str,unordered_map<string,vector<pair<pair<string,string>,pair<double,double>>>>& batch);
	void index_file(const string& str);
	void index_file_of_file(const string& str);
	void dump_disk(const string& output_file);
	void compress_index(uint32_t compression_level);
	void decompress_index();



	//USAGE FUNCTIONS
	minimizer mantis(uint64_t n);
	pair<vector<minimizer>,vector<uint64_t>> minhash_sketch_partition(const string& reference,uint32_t& active_minimizer);
	vector<minimizer>  minhash_sketch_partition_solid_kmers(const string& reference,uint32_t& active_minimizer);
	void ground_truth(const string& sequence,const string& file, double jax);
	void ground_truth_batch(vector<pair<pair<string,string>,pair<double,double>>>& V,string file);
	void insert_bloom(uint64_t);
	bool check_bloom(uint64_t num);
	void get_minimizers(const string& str,vector<minimizer>&);
	void add_index(minimizer m,uint32_t i);
	void update_kmer(kmer& min, char nuc);
	void update_kmer_RC(kmer& min, char nuc);
	kmer rcb(kmer min);
	void merge_indexes(Miekki* other_index);

};




#endif
