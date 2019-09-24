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
static constexpr minimizer maximal_minimizer(-1);
static constexpr uint64_t maximal_hash(-1);



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



void Miekki::update_kmer(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdatekmer;
}



void Miekki::update_kmer_RC(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*kmer_size-2));
}



kmer Miekki::rcb(kmer min){
	kmer res(0);
	kmer offset(1);
	offset<<=(2*kmer_size-2);
	for(uint i(0); i<kmer_size;++i){
		res+=(3-(min%4))*offset;
		min>>=2;
		offset>>=2;
	}
	return res;
}


void print_bin(uint64_t n){
	uint64_t mask=1;
	mask<<=63;
	for(uint i(0);i<64;++i){
		cout<<n/mask;
		if(n/mask==1){n-=mask;}
		mask>>=1;
	}
	cout<<"\n";
}


minimizer Miekki::mantis(uint64_t n){
	//~ return n%(maximal_minimizer);//REGULAR MINHASH
	if(n==0){return (maximal_minimizer);}
	int64_t prefix((asm_log2(n)));
	int64_t exp=max(prefix-32+(int64_t)number_minimizer_log2,(int64_t)0);
	//~ cout<<exp<<endl;
	//~ print_bin(n);
	//~ cout<<zero_chain<<endl;
	//~ if(zero_chain>31){return (maximal_minimizer);}
	//~ uint64_t prefix2((mylog2(n)));
	//~ uint64_t prefix3(floor(log2(n)));
	//~ if(prefix!=prefix2 or prefix2!=prefix3){
		//~ cout<<n<<" "<<prefix<<" "<<prefix2<<" "<<prefix3<<endl;
	//~ }
	int offset=prefix-(number_bit_minimizer-number_bit_mantis);
	offset=max(offset,0);
	uint64_t suffix (n-((uint64_t)1<<prefix));
	suffix>>=offset;
	uint64_t res=suffix;
	res+=((exp)<<(number_bit_minimizer-number_bit_mantis));
	//~ print_bin(res);cin.get();
	return res;
}



static uint8_t mask[8]={0b00000001,0b00000010,0b00000100,0b00001000,0b00010000,0b00100000,0b01000000,0b10000000};



void Miekki::insert_bloom(uint64_t num){
	for(uint32_t i(0);i<number_hash;++i){
		//~ uint64_t hash(universal_hash(num,i)%bloom_size);
		uint64_t hash(universal_hash(num,i)>>bloom_size_log2);
		uint8_t hit(hash%8);
		uint8_t cell(Bloom_Filter[hash/8]);
		if((cell && mask[hit])==0){
			Bloom_Filter[hash/8]+=(1<<hit);
		}
	}
}



bool Miekki::check_bloom(uint64_t num){
	for(uint32_t i(0);i<number_hash;++i){
		//~ uint64_t hash(universal_hash(num,i)%bloom_size);
		uint64_t hash(universal_hash(num,i)>>bloom_size_log2);
		uint8_t cell(Bloom_Filter[hash/8]);
		uint8_t hit(hash%8);
		if((cell && mask[hit])==0){
			return false;
		}
	}
	return true;
}



pair<vector<minimizer>,vector<uint64_t>>  Miekki::minhash_sketch_partition(const string& reference,uint32_t& active_minimizer){
	uint64_t mask((uint64_t)1<<(64-number_minimizer_log2));
	vector<minimizer> res(number_minimizer,maximal_minimizer);
	vector<uint64_t> res2(number_minimizer,maximal_hash);
	vector<bool> duplicate(number_minimizer,false);
	//~ cout<<reference.substr(0,100)<<endl;
	//~ cout<<revComp(reference.substr(0,100))<<endl;
	active_minimizer=(0);
	uint64_t S_kmer(str2numstrand(reference.substr(0,kmer_size-1)));
	//~ cout<<kmer2str(S_kmer,kmer_size)<<endl;;
	uint64_t RC_kmer(rcb(S_kmer));
	//~ cout<<kmer2str(RC_kmer,kmer_size)<<endl;;
	for(uint i(0);i+kmer_size<reference.size();++i){
		update_kmer(S_kmer,reference[i+kmer_size-1]);
		update_kmer_RC(RC_kmer,reference[i+kmer_size-1]);
		//~ cout<<kmer2str(S_kmer,kmer_size)<<endl;
		//~ cout<<kmer2str(RC_kmer,kmer_size)<<endl;cin.get();
		uint64_t anc(min(S_kmer,RC_kmer));
		anc=revhash64(anc);
		uint64_t bucket(anc/mask);
		uint64_t value(anc%(mask));
		value=mantis(value);
		if((minimizer)value<res[bucket]){
			duplicate[bucket]=false;
			if(res[bucket]==maximal_minimizer){
				++active_minimizer;
				//~ cout<<bucket<<" ";
			}
			res[bucket]=(minimizer)value;
			res2[bucket]=anc;
		}else if((minimizer)value==res[bucket] and res2[bucket]!=anc){
			duplicate[bucket]=true;
		}
	}
	//~ for(uint i(0);i<number_minimizer;++i){
		//~ if(duplicate[i]){
			//~ res[i]=res2[i]=maximal_minimizer;
			//~ --active_minimizer;
		//~ }
	//~ }
	//~ if(active_minimizer*5<reference.size() and active_minimizer*2<number_minimizer){
		//~ cout<<"WEIRDD"<<endl;
		//~ cout<<reference.size()<<" "<<active_minimizer<<endl;
		//~ cin.get();
		//~ cout<<reference<<endl;cin.get();
	//~ }
	return {res,res2};
}



vector<minimizer>  Miekki::minhash_sketch_partition_solid_kmers(const string& reference,uint32_t& active_minimizer){
	pair<vector<minimizer>,vector<uint64_t>> sketch( minhash_sketch_partition(reference,active_minimizer));
	for(uint i(0);i<sketch.first.size();++i){
		if(not check_bloom(sketch.second[i])){
			sketch.first[i]=maximal_minimizer;
			active_minimizer--;
		}
	}
	return sketch.first;
}



template<typename T> void Miekki::minhash_sketch_partition_solid_kmers(const string& reference,uint32_t& active_minimizer, T out){
    pair<vector<minimizer>,vector<uint64_t>> sketch( minhash_sketch_partition(reference,active_minimizer));
    for(uint i(0);i<sketch.first.size();++i){
        if(not check_bloom(sketch.second[i])){
            sketch.first[i]=maximal_minimizer;
            active_minimizer--;
        }
        out[i] = sketch.first[i];
    }

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
	auto double_sketch(minhash_sketch_partition(str,actived_minimizer));
	#pragma omp critical(update_index_structure)
	{
		for(uint i(0);i<double_sketch.first.size();++i){
			add_index(double_sketch.first[i],i);
			if(double_sketch.first[i]!=maximal_minimizer){
				approx_cardinality+=(1/(pow(2,(double_sketch.first[i]>>(number_bit_minimizer-number_bit_mantis)))));
				active_minimizer++;
				if(not check_bloom(double_sketch.second[i])){
					#pragma omp critical(bloom)
					{
						insert_bloom(double_sketch.second[i]);
					}
				}
			}
		}
		file_names.push_back(title);
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



void Miekki::insert_sequences(const vector<pair<string,string>>& Vstr){
	uint32_t actived_minimizer;
	vector<pair<vector<minimizer>,vector<uint64_t>>> Vsketch;
	for(uint g(0);g<Vstr.size();g++){
		auto double_sketch(minhash_sketch_partition(Vstr[g].first,actived_minimizer));
		Vsketch.push_back(double_sketch);
	}

	#pragma omp critical(update_index_structure)
	{
		for(uint g(0);g<Vstr.size();g++){
            double approx_cardinality(0);
            uint32_t active_minimizer=0;
			for(uint i(0);i<Vsketch[g].first.size();++i){
				add_index(Vsketch[g].first[i],i);
				if(Vsketch[g].first[i]!=maximal_minimizer){
					approx_cardinality+=(1/(pow(2,(Vsketch[g].first[i]>>(number_bit_minimizer-number_bit_mantis)))));
					active_minimizer++;
					if(not check_bloom(Vsketch[g].second[i])){
						#pragma omp critical(bloom)
						{
							insert_bloom(Vsketch[g].second[i]);
						}
					}
				}
			}
			file_names.push_back(Vstr[g].second);
			index_size++;
			sketch_size.push_back(active_minimizer);
			approx_cardinality=((0.72134*(active_minimizer*active_minimizer)/approx_cardinality));
			if(approx_cardinality>Vstr[g].first.size()){
				genome_size.push_back(Vstr[g].first.size());
			}else{
                genome_size.push_back(uint64_t(approx_cardinality));
			}
		}
	}
}



matrix::vector<score_t> Miekki::query_sequence(const string& str,uint32_t& active_minimizer){
    matrix::vector<score_t> result(index_size, 0);
	auto double_sketch(minhash_sketch_partition(str,active_minimizer));
	active_minimizer=0;
	vector<minimizer> slice;
	for(uint32_t i(0);i<double_sketch.first.size();++i){
		minimizer query_minimizer(double_sketch.first[i]);
		if(query_minimizer!=maximal_minimizer){
			if(not check_bloom(double_sketch.second[i])){
				continue;
			}
			active_minimizer++;
			get_minimizers(index[i],slice);
			for(uint32_t j(0);j<slice.size();++j){
				if(query_minimizer==slice[j]){
                    result[j]++;
				}else{
				}
			}
		}
	}
	return result;
}



matrix::matrix<score_t> Miekki::query_sequences(vector<pair<string,uint32_t>>& batch){
	cout<<"-"<<flush;

    matrix::matrix<minimizer> sketch_batch(batch.size(), number_minimizer);
	for(uint i(0);i<batch.size();++i){
        minhash_sketch_partition_solid_kmers(batch[i].first,batch[i].second, sketch_batch.row(i));
	}

    matrix::matrix<score_t> result(batch.size(), index_size, 0);
	vector<minimizer> slice;
	//FOREACH MINIMIZER
    for(idx_t i_mini = 0 ; i_mini < idx_t(number_minimizer) ; ++i_mini){
		slice.clear();
        for(idx_t i_batch = 0 ; i_batch < idx_t(batch.size()) ; ++i_batch){
            minimizer mini = sketch_batch(i_batch, i_mini);
            if(mini == maximal_minimizer){
                continue;
            }
            if(slice.empty()){ get_minimizers(index[i_mini],slice); }
            for(idx_t genome_id = 0 ; genome_id < idx_t(index_size) ; ++genome_id){
                if(mini == slice[genome_id]){
                    ++result(i_batch, genome_id);
                }
			}
		}
	}

    return result;
}



void Miekki::filter_results(matrix::refvec<score_t> results, size_t nresults, score_t min_score, double min_intersection, std::vector<similarity_score>& min_heap) {
    const auto compare = [](const similarity_score& a, const similarity_score& b) { return a.intersection > b.intersection; };

    for(idx_t genome_id = 0; genome_id < results.size(); genome_id++) {
        score_t score = results[genome_id];
        if(score < min_score) continue;
        double jaccard = double(score) / sketch_size[genome_id];
        double intersection = jaccard * genome_size[genome_id];
        if(intersection < min_intersection) continue;

        if(min_heap.size() >= nresults) {
            if(min_heap.front().intersection > intersection) continue; // New element is less than the minimum
            std::pop_heap(min_heap.begin(), min_heap.end(), compare);
            min_heap.pop_back();
        }

        min_heap.emplace_back(similarity_score{seqid_t(genome_id), score, jaccard, intersection});
        std::push_heap(min_heap.begin(), min_heap.end(), compare);
    }

    std::sort_heap(min_heap.begin(), min_heap.end(), compare);
}



std::vector<similarity_score> Miekki::filter_results(matrix::refvec<score_t> results, size_t nresults, score_t min_score, double min_intersection) {
    std::vector<similarity_score> min_heap;
    filter_results(results, nresults, min_score, min_intersection, min_heap);
    return min_heap;
}



std::vector<matrix::vector<similarity_score>> Miekki::filter_results(matrix::matrix<score_t> results, size_t nresults, score_t min_score, double min_intersection) {
    std::vector<similarity_score> min_heap;
    std::vector<matrix::vector<similarity_score>> filtered_results;
    filtered_results.reserve(results.rows());

    for(idx_t i=0 ; i < results.rows() ; i++) {
        filter_results(results.row(i), nresults, min_score, min_intersection, min_heap);
        filtered_results.emplace_back(min_heap.size());
        std::copy(min_heap.begin(), min_heap.end(), filtered_results.back().begin());
        min_heap.clear();
    }

    return filtered_results;
}



void Miekki::query_file(const string& str){
	if(not exists_test(str)){cout<<"File problem"<<endl;return;}
	auto in=new zstr::ifstream(str);

	#pragma omp parallel num_threads(core_number)
	{
		vector<pair<string,uint32_t>> batch;
		vector<string> names;
        string toWrite;

        auto do_batch = [&]() {
            auto results = filter_results(query_sequences(batch), 10, 10, 0.5*threshold);
            for(uint i_batch(0);i_batch<batch.size();i_batch++) {
                auto& result = results[i_batch];
                toWrite+=names[i_batch]+":";
                for(auto& sim: result) {
                    toWrite+=(to_string(sim.genome)+"\t"+to_string(sim.matches)+"\t"+to_string(uint(sim.intersection))+"\t"+ to_string(sim.jaccard)+";");
                }
                toWrite+="\n";
            }
            if(toWrite.size()>0){
                #pragma omp critical(outfile)
                {
                    *out<<toWrite;
                }
                toWrite.clear();
            }
            batch.clear();
            names.clear();
        };


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
			}else{
				batch.push_back({ref,0});
				names.push_back(head);
				if(batch.size()>200){
					do_batch();
				}
			}
		}

		if(not batch.empty()) {
			do_batch();
		}
	}
	*out<<flush;
	delete in;
}



void Miekki::query_whole_file(const string& str){
	if(not exists_test(str)){cout<<"File problem"<<endl;return;}
	auto in=new zstr::ifstream(str);
	string ref,line;
	while(not in->eof()){
		getline(*in,line);
		if(line[0]=='>'){
		}else{
			ref+=line;
		}
	}
	if(ref.size()<kmer_size){return;}
	uint active_minimizer(0);
    auto result = filter_results(query_sequence(ref,active_minimizer), 10, 10, 0.5*threshold);
	string toWrite;

    for(auto& sim: result) {
        toWrite+=(to_string(sim.genome)+"\t"+to_string(sim.matches)+"\t"+to_string(uint(sim.intersection))+"\t"+ to_string(sim.jaccard)+";");
    }
	if(toWrite.size()>0){
		#pragma omp critical(outfile)
		{
			*out<<str<<":"<<toWrite<<"\n";
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
	#pragma omp parallel num_threads(core_number)
	{
		vector<pair<string,string>> seq2index;
		while(not in->eof()){
			string file_name;
			#pragma omp critical(fof)
			{
				getline(*in,file_name);
			}
			if(file_name.size()>3){
				if(not exists_test(file_name)){
					cout<<"Missed file: "<<file_name<<endl;
				}else{
					auto in2=new zstr::ifstream(file_name);
					string ref2,line;
					while(not in2->eof()){
						getline(*in2,line);
						if(line[0]=='>'){
						}else{
							ref2+=line;
						}
					}
					delete in2;
					if(ref2.size()>=kmer_size){
						seq2index.push_back({ref2,file_name});
						if(seq2index.size()>10){
							insert_sequences(seq2index);
							seq2index.clear();
						}
						cout<<"-"<<flush;
					}
				}
			}
		}
		insert_sequences(seq2index);
	}
	cout<<endl;
	cout<<"Reference indexed: "<<index_size<<endl;
	if(Bloom_Filter_Protection){
		cout<<"BF size:"<<intToString(bloom_size)<<endl;
	}
	delete in;
}



void Miekki::query_file_of_file(const string& str){
	if(not exists_test(str)){
		cout<<"Missed file of file: "<<str<<endl;
		return;
	}
	auto in=new zstr::ifstream(str);

	#pragma omp parallel num_threads(core_number)
	while(not in->eof()){
		string ref;
		#pragma omp critical(fof)
		{
			getline(*in,ref);
		}
		if(ref.size()>3){
			query_whole_file(ref);
			cout<<"-"<<flush;
		}
	}
	delete in;
}



void Miekki::query_file_of_file_exact(const string& str){
	if(not exists_test(str)){
		cout<<"Missed file of file: "<<str<<endl;
		return;
	}
	auto in=new zstr::ifstream(str);

	#pragma omp parallel num_threads(core_number)
	{

		unordered_map<string,vector<pair<pair<string,string>,pair<double,double>>>> batch;
		while(not in->eof()){
			string ref;
			#pragma omp critical(fof)
			{
				getline(*in,ref);
			}
			if(ref.size()>3){
				query_whole_file_exact(ref,batch);
				cout<<"-"<<flush;
			}
		}
		for (auto itr = batch.begin(); itr != batch.end(); ++itr) {
			ground_truth_batch(itr->second,itr->first);
		}
	}

	*out<<flush;
	delete in;
}



void Miekki::dump_disk(const string& output_file){
	auto out=new zstr::ofstream(output_file);
	out->write(reinterpret_cast<char*>(&kmer_size),sizeof(kmer_size));
	out->write(reinterpret_cast<char*>(&number_minimizer_log2),sizeof(number_minimizer_log2));
	out->write(reinterpret_cast<char*>(&number_bit_minimizer),sizeof(number_bit_minimizer));
	out->write(reinterpret_cast<char*>(&number_bit_mantis),sizeof(number_bit_mantis));
	out->write(reinterpret_cast<char*>(&index_size),sizeof(index_size));
	out->write(reinterpret_cast<char*>(&bloom_size_log2),sizeof(bloom_size_log2));
	out->write(reinterpret_cast<char*>(&bloom_size),sizeof(bloom_size));
	out->write(reinterpret_cast<char*>(&jaccard_estimation),sizeof(jaccard_estimation));
	out->write(reinterpret_cast<char*>(&containment_estimation),sizeof(containment_estimation));
	out->write(reinterpret_cast<char*>(&threshold),sizeof(threshold));
	out->write(reinterpret_cast<char*>(&compressed),sizeof(compressed));
	if(compressed){
		decompress_index();
	}
	for(uint32_t i(0);i<number_minimizer;++i){
		auto point =&(index[i][0]);
		out->write((char*)point,index_size*sizeof(minimizer));
	}
	auto point =&(genome_size[0]);
	out->write((char*)point,index_size*sizeof(uint64_t));
	if(bloom_size!=0){
		auto point2 =&Bloom_Filter[0];
		out->write((char*)point2,Bloom_Filter.size());
	}
	auto point2 =&(sketch_size[0]);
	out->write((char*)point2,index_size*sizeof(uint32_t));
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
	in-> read(reinterpret_cast<char *>(&bloom_size_log2), sizeof(bloom_size_log2));
	in-> read(reinterpret_cast<char *>(&bloom_size), sizeof(bloom_size));
	in-> read(reinterpret_cast<char *>(&jaccard_estimation), sizeof(jaccard_estimation));
	in-> read(reinterpret_cast<char *>(&containment_estimation), sizeof(containment_estimation));
	in-> read(reinterpret_cast<char *>(&threshold), sizeof(threshold));
	in-> read(reinterpret_cast<char *>(&compressed), sizeof(compressed));
	string str(index_size*sizeof(minimizer),-1);
	index.assign(number_minimizer,str);
	number_hash=(5);
	maximal_minimizer=(-1);
	compressed=false;
	offsetUpdatekmer=1;
	offsetUpdatekmer<<=2*kmer_size;
	for(uint32_t i(0);i<number_minimizer;++i){
		in->read((char*)(index[i].data()),index_size*sizeof(minimizer));
	}
	genome_size.assign(index_size,0);
	in->read((char*)genome_size.data(),index_size*sizeof(uint64_t));
	if(bloom_size!=0){
		Bloom_Filter.assign(bloom_size/8,0);
		in->read((char*)(Bloom_Filter.data()),bloom_size/8);
	}
	sketch_size.assign(index_size,0);
	in->read((char*)sketch_size.data(),index_size*sizeof(uint32_t));
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
			if(ref.size()<kmer_size or (ref[0]!='A' and ref[0]!='C' and ref[0]!='G' and ref[0]!='T' and ref[0]!='N')){
				ref="";
				continue;
			}
			uint active_minimizer(0);
            auto result = filter_results(query_sequence(ref,active_minimizer), 5, 10, threshold);
            for(auto& sim: result) {
                string file_name(file_names[sim.genome]);
                batch[file_name].push_back({{ref,head},{sim.jaccard,sim.intersection}});
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



void Miekki::query_whole_file_exact(const string& str,unordered_map<string,vector<pair<pair<string,string>,pair<double,double>>>>& batch){
	if(not exists_test(str)){cout<<"File problem"<<endl;return;}
	auto in=new zstr::ifstream(str);
	//~ unordered_map<string,vector<pair<pair<string,string>,pair<double,double>>>> batch;
	string ref,head,line;
	getline(*in,head);
	while(not in->eof()){
		getline(*in,line);
		if(line.size()<kmer_size or (line[0]!='A' and line[0]!='C' and line[0]!='G' and line[0]!='T' and line[0]!='N')){
			continue;
		}else{
			ref+=line;
		}
	}
	uint active_minimizer(0);
    auto result = filter_results(query_sequence(ref,active_minimizer), 5, 5, threshold);
    for(auto& sim: result) {
        string file_name(file_names[sim.genome]);
        batch[file_name].push_back({{ref,head},{sim.jaccard,sim.intersection}});
		if(batch[file_name].size()>=100){
			ground_truth_batch(batch[file_name],file_name);
			batch[file_name].clear();
		}
	}
	delete in;
}



void Miekki::ground_truth_batch(vector<pair<pair<string,string>,pair<double,double>>>& V,string file){
	unordered_set<uint64_t> set_B;
	//~ cout<<"go"<<endl;
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
				for(uint i(0);i+kmer_size-1<ref.size();++i){
					set_B.insert(str2num(ref.substr(i,kmer_size)));
				}
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
			if(nb_inter>0){
			//~ if(nb_inter>1*threshold){
				//OUTPUT RESULT
				#pragma omp critical(outfile)
				{
					*out<<real_jax<<"	"<<jaccard_estimation<<"	"<<nb_inter<<"	"<<intersection_estimation<<"	"<<V[i].first.second<<"	"<<file<<"\n";
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
    for(uint32_t i(0);i<str.size();i+=(number_bit_minimizer==16 ? 2 : 1)){
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


void Miekki::merge_indexes(Miekki* other_index){
	uint32_t i;
	#pragma omp parallel for num_threads(core_number)
	for(i=(0);i<number_minimizer;++i){
		index[i].insert(index[i].end(),other_index->index[i].begin(),other_index->index[i].end());
	}
	//TODO MERGE BF
	delete other_index;

}

