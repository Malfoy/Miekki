#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <bitset>
#include <omp.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include "zstr.hpp"
#include "Miekki.h"
#include "utils.h"



using namespace std;



uint64_t TP(0),FP(0),FN(0);
minimizer maximal_minimizer(-1);
uint64_t maximal_hash(-1);



uint64_t hash_str(const string& str){
	uint64_t res(1);
	for(uint i(0);i<str.size();++i){
		res*=str[i];
	}
	return res;
}



inline bool exists_test (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}



minimizer Miekki::mantis(uint64_t n){
	//~ return n%(maximal_minimizer);//REGULAR MINHASH
	uint64_t prefix(floor(log2(n)));
	int offset=prefix-(number_bit_minimizer-number_bit_mantis);
	offset=max(offset,0);
	uint64_t suffix (n-((uint64_t)1<<prefix));
	suffix>>=offset;
	uint64_t res=suffix;
	res+=((prefix)<<(number_bit_minimizer-number_bit_mantis));
	return res;
}



static uint8_t mask[8]={0b00000001,0b00000010,0b00000100,0b00001000,0b00010000,0b00100000,0b01000000,0b10000000};



void Miekki::insert_bloom(uint64_t num){
	for(uint32_t i(0);i<number_hash;++i){
		uint64_t hash(universal_hash(num,i)%bloom_size);
		uint8_t hit(hash%8);
		uint8_t cell(Bloom_Filter[hash/8]);
		if((cell && mask[hit])==0){
			Bloom_Filter[hash/8]+=(1<<hit);
		}
	}
}



bool Miekki::check_bloom(uint64_t num){
	for(uint32_t i(0);i<number_hash;++i){
		uint64_t hash(universal_hash(num,i)%bloom_size);
		uint8_t cell(Bloom_Filter[hash/8]);
		uint8_t hit(hash%8);
		if((cell && mask[hit])==0){
			return false;
		}
	}
	return true;
}



vector<minimizer>  Miekki::minhash_sketch_partition(const string& reference){
	uint64_t mask((uint64_t)1<<(64-number_minimizer_log2));
	vector<minimizer> res(number_minimizer,0);
	uint fill(0);
	for(uint i(0);i+kmer_size-1<reference.size();++i){
		uint64_t anc(str2num(reference.substr(i,kmer_size)));
		//~ anc=unrevhash64(anc);
		uint64_t bucket(anc/mask);
		uint64_t value(anc%(mask));
		value=mantis(value);
		if(res[bucket]==0){
			fill++;
		}
		res[bucket]=max((minimizer)value,res[bucket]);
	}
	return res;
}



pair<vector<minimizer>,vector<uint64_t>>  Miekki::minhash_sketch_partition_keep_kmer(const string& reference,uint32_t& active_minimizer){
	uint64_t mask((uint64_t)1<<(64-number_minimizer_log2));
	vector<minimizer> res(number_minimizer,maximal_minimizer);
	vector<uint64_t> res2(number_minimizer,maximal_hash);
	vector<bool> duplicate(number_minimizer,false);
	active_minimizer=(0);
	for(uint i(0);i+kmer_size-1<reference.size();++i){
		uint64_t anc(str2num(reference.substr(i,kmer_size)));
		anc=unrevhash64(anc);
		uint64_t bucket(anc/mask);
		uint64_t value(anc%(mask));
		value=mantis(value);
		if((minimizer)value<res[bucket]){
			duplicate[bucket]=false;
			if(res[bucket]==maximal_minimizer){
				++active_minimizer;
			}
			res[bucket]=(minimizer)value;
			res2[bucket]=anc;
		}else if((minimizer)value==res[bucket] and res2[bucket]!=anc){
			duplicate[bucket]=true;
		}
	}
	for(uint i(0);i<number_minimizer;++i){
		if(duplicate[i]){
			res[i]=res2[i]=maximal_minimizer;
			--active_minimizer;
		}
	}
	return {res,res2};
}



void Miekki::add_index(minimizer m,uint32_t i){
	if(number_bit_minimizer==16){
		index[i].push_back(m/256);
		index[i].push_back(m%256);
		//~ exit(0);
	}else if(number_bit_minimizer==8){
		index[i].push_back(m);
	}else{
		cout<<"not implemented"<<endl;
		exit(0);
	}
}



void Miekki::insert_sequence(const string& str, const string& title){
	uint32_t actived_minimizer;
	double approx_cardinality(0);
	double active_minimizer=0;
	if(not Bloom_Filter_Protection){
		vector<minimizer> sketch(minhash_sketch_partition(str));
		for(uint i(0);i<sketch.size();++i){
			omp_set_lock(&lock[i%10000]);
			add_index(sketch[i],i);
			omp_unset_lock(&lock[i%10000]);
		}
	}else{
		auto double_sketch(minhash_sketch_partition_keep_kmer(str,actived_minimizer));
		#pragma omp critical(update_index_structure)
		{
			for(uint i(0);i<double_sketch.first.size();++i){
				add_index(double_sketch.first[i],i);
				if(double_sketch.first[i]!=maximal_minimizer){
					approx_cardinality+=(1/(pow(2,(double_sketch.first[i]>>(number_bit_minimizer-number_bit_mantis)))));
					active_minimizer++;
					#pragma omp critical(bloom)
					{
						//~ Bloom_Filter_Forever.insert(double_sketch.second[i]);
						insert_bloom(double_sketch.second[i]);
					}
				}
			}
			file_names.push_back(title);
			query_output.push_back({((uint32_t)index[0].size()-(uint32_t)1),0});//TODO OPTIM
			index_size++;
			sketch_size.push_back(active_minimizer);
			approx_cardinality=((0.72134*(active_minimizer*active_minimizer)/approx_cardinality));
			if(approx_cardinality>str.size()){
				genome_size.push_back(str.size());
			}else{
				genome_size.push_back(approx_cardinality);
			}
		}
	}
}



bool compare_Similarity_score(const similarity_score &a, const similarity_score &b){
    return a.score > b.score;
}



vector<similarity_score> Miekki::query_sequence(const string& str,uint32_t& active_minimizer){
	vector<similarity_score> result=(query_output);
	if(not Bloom_Filter_Protection){
		vector<minimizer> sketch(minhash_sketch_partition(str));
		vector<minimizer> slice;
		for(uint32_t i(0);i<sketch.size();++i){
			minimizer query_minimizer(sketch[i]);
			if(query_minimizer!=maximal_minimizer){
				get_minimizers(index[i],slice);
				for(uint32_t j(0);j<slice.size();++j){
					if(query_minimizer==slice[j]){
						result[j].score++;
					}
				}
			}
		}
	}else{
		auto double_sketch(minhash_sketch_partition_keep_kmer(str,active_minimizer));
		active_minimizer=(0);
		vector<minimizer> slice;
		for(uint32_t i(0);i<double_sketch.first.size();++i){
			minimizer query_minimizer(double_sketch.first[i]);
			if(query_minimizer!=maximal_minimizer){
				//~ if(Bloom_Filter_Forever.count(double_sketch.second[i])==0){
				if(not check_bloom(double_sketch.second[i])){
					continue;
				}
				active_minimizer++;
				get_minimizers(index[i],slice);
				for(uint32_t j(0);j<slice.size();++j){
					if(query_minimizer==slice[j]){
						result[j].score++;
					}else{
					}
				}
			}
		}
	}
	sort(result.begin(),result.end(),compare_Similarity_score);
	return result;
}



void Miekki::query_file(const string& str){
	if(not exists_test(str)){cout<<"File problem"<<endl;return;}
	auto in=new zstr::ifstream(str);

	#pragma omp parallel num_threads(core_number)
	{
		while(not in->eof()){
			string ref,head;
			#pragma omp critical(readfile)
			{
				getline(*in,head);
				getline(*in,ref);
			}
			if(ref.size()<kmer_size){
				ref="";
				continue;
			}
			uint active_minimizer(0);
			auto result(query_sequence(ref,active_minimizer));
			string toWrite;
			for(uint32_t i(0);i<min((uint32_t)5,(uint32_t)result.size());++i){//RETURN THE 10 BEST HITS
				//FILTER ON HITS
				if(result[i].score<10){
					break;
				}
				double jaccard_value(((double)(result[i].score)/sketch_size[result[i].sequence_identifier]));
				double intersection_value(((double)(result[i].score)*genome_size[result[i].sequence_identifier])/sketch_size[result[i].sequence_identifier]);
				if(intersection_value>0.5*threshold){
					toWrite+=(file_names[i]+"	"+to_string(result[i].score)+"	"+to_string((uint)intersection_value)+"	"+ to_string(jaccard_value)+";");
				}else{
					break;
				}
			}
			if(toWrite.size()>0){
				#pragma omp critical(outfile)
				{
					*out<<head<<":"<<toWrite<<"\n";
				}
			}
		}
	}
	*out<<flush;
	delete in;
}



void Miekki::index_file(const string& str){
	if(not exists_test(str)){
		cout<<"Missed file: "<<str<<endl;;
		return;
	}
	auto in=new zstr::ifstream(str);
	string ref,line;
	while(not in->eof()){
		getline(*in,line);
		if(line[0]=='>'){
		}else{
			ref+=line;
		}
	}
	if(ref.size()>=kmer_size){
		insert_sequence(ref,str);
	}
	delete in;
}



void Miekki::index_file_of_file(const string& str){
	if(not exists_test(str)){
		cout<<"Missed file of file: "<<str<<endl;
		return;
	}
	auto in=new zstr::ifstream(str);

	//~ #pragma omp parallel
	while(not in->eof()){
		string ref;
		#pragma omp critical(fof)
		{
			getline(*in,ref);
		}
		if(ref.size()>3){
			index_file(ref);
			cout<<"-"<<flush;
		}
	}
	cout<<endl;
	cout<<"Reference indexed: "<<index_size<<endl;
	if(Bloom_Filter_Protection){
		cout<<"BF size:"<<intToString(bloom_size)<<endl;
	}
	delete in;
}



void Miekki::dump_disk(const string& output_file){
	auto out=new zstr::ofstream(output_file);
	out->write(reinterpret_cast<char*>(&kmer_size),sizeof(kmer_size));
	out->write(reinterpret_cast<char*>(&number_minimizer_log2),sizeof(number_minimizer_log2));
	out->write(reinterpret_cast<char*>(&number_bit_minimizer),sizeof(number_bit_minimizer));
	out->write(reinterpret_cast<char*>(&number_bit_mantis),sizeof(number_bit_mantis));
	out->write(reinterpret_cast<char*>(&index_size),sizeof(index_size));
	out->write(reinterpret_cast<char*>(&bloom_size),sizeof(bloom_size));
	out->write(reinterpret_cast<char*>(&jaccard_estimation),sizeof(jaccard_estimation));
	out->write(reinterpret_cast<char*>(&containment_estimation),sizeof(containment_estimation));
	out->write(reinterpret_cast<char*>(&Bloom_Filter_Protection),sizeof(Bloom_Filter_Protection));
	if(compressed){
		decompress_index();
	}
	for(uint32_t i(0);i<number_minimizer;++i){
		auto point =&(index[i][0]);
		out->write((char*)point,index_size*sizeof(minimizer));
	}
	if(bloom_size!=0){
		auto point2 =&Bloom_Filter[0];
		out->write((char*)point2,Bloom_Filter.size());
	}
	auto point =&(sketch_size[0]);
	out->write((char*)point,index_size*sizeof(uint32_t));
	delete out;
}



Miekki::Miekki(const string& input){
	if(not exists_test(input)){
		cout<<"File problem"<<endl;
		return;
	}
	auto in=new zstr::ifstream(input);
	in-> read(reinterpret_cast<char *>(&kmer_size), sizeof(kmer_size));
	in-> read(reinterpret_cast<char *>(&number_minimizer_log2), sizeof(number_minimizer_log2));
	number_minimizer=1;
	number_minimizer<<=number_minimizer_log2;
	in-> read(reinterpret_cast<char *>(&number_bit_minimizer), sizeof(number_bit_minimizer));
	in-> read(reinterpret_cast<char *>(&number_bit_mantis), sizeof(number_bit_mantis));
	in-> read(reinterpret_cast<char *>(&index_size), sizeof(index_size));
	in-> read(reinterpret_cast<char *>(&bloom_size), sizeof(bloom_size));
	in-> read(reinterpret_cast<char *>(&jaccard_estimation), sizeof(jaccard_estimation));
	in-> read(reinterpret_cast<char *>(&containment_estimation), sizeof(containment_estimation));
	in-> read(reinterpret_cast<char *>(&Bloom_Filter_Protection), sizeof(Bloom_Filter_Protection));
	string str(index_size*sizeof(minimizer),-1);
	index.assign(number_minimizer,str);
	number_hash=(5);
	maximal_minimizer=(-1);
	compressed=false;
	for(uint32_t i(0);i<number_minimizer;++i){
		in->read((char*)(index[i].data()),index_size*sizeof(minimizer));
	}
	if(bloom_size!=0){
		Bloom_Filter.assign(bloom_size/8,0);
		in->read((char*)(Bloom_Filter.data()),bloom_size/8);
	}
	sketch_size.assign(index_size,0);
	in->read((char*)sketch_size.data(),index_size*sizeof(uint32_t));
	query_output.resize(index_size,{((uint32_t)index[0].size()-(uint32_t)1),0});
}



void Miekki::query_file_exact(const string& str){
	if(not exists_test(str)){cout<<"File problem"<<endl;return;}
	auto in=new zstr::ifstream(str);
	#pragma omp parallel num_threads(core_number)
	{
		unordered_map<string,vector<pair<pair<string,string>,pair<double,double>>>> batch;
		while(not in->eof()){
			string ref,head;
			#pragma omp critical(readfile)
			{
				getline(*in,head);
				getline(*in,ref);
			}
			if(ref.size()<kmer_size){
				ref="";
				continue;
			}
			uint active_minimizer(0);
			auto result(query_sequence(ref,active_minimizer));
			for(uint32_t i(0);i<min((uint32_t)5,(uint32_t)result.size());++i){
				//FILTER ON HITS
				if(result[i].score<10){break;}
				//RESULT ESTIMATION
				double jax,inter;
				jax=((double)(result[i].score)/sketch_size[result[i].sequence_identifier]);
				inter=((double)(result[i].score)*genome_size[result[i].sequence_identifier])/sketch_size[result[i].sequence_identifier];
				//FILTER ON ration
				//~ if(jax<threshold){break;}
				string file_name(file_names[result[i].sequence_identifier]);
				batch[file_name].push_back({{ref,head},{jax,inter}});
				if(batch[file_name].size()>=100){
					ground_truth_batch(batch[file_name],file_name);
					batch[file_name].clear();
				}else{
				}
			}
		}
		for (auto itr = batch.begin(); itr != batch.end(); ++itr) {
			ground_truth_batch(itr->second,itr->first);
		}
	}

	*out<<flush;
	delete in;
}



void Miekki::ground_truth_batch(vector<pair<pair<string,string>,pair<double,double>>>& V,string file){
	unordered_set<uint64_t> set_B;
	//~ double ln(parse_line_number(file));
	if(not exists_test(file)){
		cout<<"File problem: "<<file<<endl;
		return;
	}
	auto in=new zstr::ifstream(file);
	string ref,line;
	uint actual_line(1);
	while(not in->eof()){
		getline(*in,line);
		if(line[0]=='>'){
			if(ref.size()>=kmer_size){
				//~ if(actual_line==ln){
					for(uint i(0);i+kmer_size-1<ref.size();++i){
						set_B.insert(str2num(ref.substr(i,kmer_size)));
					}
					ref="";
					//~ break;
				//~ }
				ref="";
				actual_line++;
			}
		}else{
			ref+=line;
		}
	}
	if(ref.size()>=kmer_size){
		for(uint i(0);i+kmer_size-1<ref.size();++i){
			set_B.insert(str2num(ref.substr(i,kmer_size)));
		}
		ref="";
	}
	delete in;
	uint i(0);
	#pragma omp parallel for num_threads(core_number)
	for(i=0;i<V.size();++i){
		unordered_set<uint64_t> set_A;
		double nb_inter(0),nb_union(set_B.size());
		string sequence(V[i].first.first);
		double jaccard_estimation(V[i].second.first);
		double intersection_estimation(V[i].second.second);
		for(uint i(0);i+kmer_size-1<sequence.size();++i){
			uint64_t num(str2num(sequence.substr(i,kmer_size)));
			if(set_A.count(num)==0){
				if(set_B.count(num)==1){
					nb_inter++;
				}else{
					nb_union++;
				}
				set_A.insert(num);
			}
		}
		if(nb_inter>0){
			//REAL VALUE COMPUTATION
			double real_jax((double)nb_inter);
			real_jax/=nb_union;
			//FILTER REAL DATA
			if(nb_inter>1*threshold){
				//OUTPUT RESULT
				#pragma omp critical(outfile)
				{
					*out<<real_jax<<"	"<<jaccard_estimation<<"	"<<nb_inter<<"	"<<intersection_estimation<<"	"<<V[i].first.second<<"	"<<file<<endl;
				}
			}
		}
		set_A.clear();
	}
}



void Miekki::compress_index(uint32_t compression_level){
	for(uint32_t i(0);i<number_minimizer;++i){
		index[i]=compress_string(index[i],compression_level);
	}
	compressed=true;
}



void Miekki::decompress_index(){
	for(uint32_t i(0);i<number_minimizer;++i){
		index[i]=decompress_string(index[i]);
	}
	compressed=false;
}



void Miekki::get_minimizers(const string& cstr,vector<minimizer>& V){
	string str=compressed?decompress_string(cstr):cstr;
	if(number_bit_minimizer==16){
		V.resize(str.size()/2);
	}else{
		V.resize(str.size());
	}
	for(uint32_t i(0);i<str.size();i+=2){
		if(number_bit_minimizer==16){
			V[i/2]=(uint8_t)str[i]*256+(uint8_t)str[i+1];
		}else if(number_bit_minimizer==8){
			V[i]=(uint8_t)str[i];
		}else{
			cout<<"not implemented"<<endl;
			exit(0);
		}
	}
}

