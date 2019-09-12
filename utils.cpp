#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <zlib.h>
#include "utils.h"
#include <algorithm>





using namespace std;
using namespace chrono;



  uint64_t nuc2int(char c){
	switch(c){
		/*
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		*/
		//~ case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		//~ case 'N': return 0;
		default: return 0;
	}
	//~ cout<<"Unknow nucleotide: "<<c<<"!"<<endl;
	exit(0);
	return 0;
}



string kmer2str(uint64_t num,uint k){
	string res;
	uint64_t anc(1);
	anc<<=(2*(k-1));
	for(uint i(0);i<k;++i){
		uint nuc=num/anc;
		num=num%anc;
		if(nuc==3){
			res+="T";
		}
		if(nuc==2){
			res+="G";
		}
		if(nuc==1){
			res+="C";
		}
		if(nuc==0){
			res+="A";
		}
		if (nuc>=4){
			cout<<nuc<<endl;
			cout<<"WTF"<<endl;
		}
		anc>>=2;
	}
	return res;
}




uint64_t asm_log2(const uint64_t x) {
  uint64_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}


uint64_t mylog2 (uint64_t val) {
    if (val == 0) return 0;
    if (val == 1) return 0;
    uint64_t ret = 0;
    while (val > 1) {
        val >>= 1;
        ret++;
    }
    return ret;
}



  uint64_t nuc2intrc(char c){
	switch(c){
		/*
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		*/
		case 'A': return 3;
		case 'C': return 2;
		case 'G': return 1;
		//~ case 'T': return 0;
		//~ case 'N': return 0;
		default: return 0;
	}
	//~ cout<<"Unknow nucleotide: "<<c<<"!"<<endl;
	exit(0);
	return 0;
}



double parse_line_number(string& str){
	size_t pos(str.find("ln"));
	if(pos==string::npos){
		return 0;
	}
	uint i(1);
	while(i+pos+2<str.size()){
		++i;
	}
	double res(stof(str.substr(pos+2,i)));
	str=str.substr(0,pos);
	return res;
}



string intToString(uint64_t n){
	if(n<1000){
		return to_string(n);
	}
	string end(to_string(n%1000));
	if(end.size()==3){
		return intToString(n/1000)+","+end;
	}
	if(end.size()==2){
		return intToString(n/1000)+",0"+end;
	}
	return intToString(n/1000)+",00"+end;
}



 uint32_t revhash ( uint32_t x ) {
	x = ( ( x >> 16 ) ^ x ) * 0x2c1b3c6d;
	x = ( ( x >> 16 ) ^ x ) * 0x297a2d39;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



 uint32_t unrevhash ( uint32_t x ) {
	x = ( ( x >> 16 ) ^ x ) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ( ( x >> 16 ) ^ x ) * 0x64ea2d65;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



 uint64_t revhash64 ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



 uint64_t unrevhash64 ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



uint64_t universal_hash ( uint64_t x, uint32_t k) {
	return unrevhash64(x)+(k*69*revhash64(x))%1024;
}



char revCompChar(char c) {
	switch (c) {
		//~ case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		case 't': return 'A';
		case 'c': return 'G';
		case 'g': return 'C';
	}
	return 'T';
}



string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



void clean(string& str){
	for(uint i(0); i< str.size(); ++i){
		switch(str[i]){
			case 'a':break;
			case 'A':break;
			case 'c':break;
			case 'C':break;
			case 'g':break;
			case 'G':break;
			case 't':break;
			case 'T':break;
			//~ case 'N':break;
			//~ case 'n':break;
			default:
			//~ str[i]='A';
			str="";
			break;
		}
	}
	transform(str.begin(), str.end(), str.begin(), ::toupper);
}



uint64_t str2numstrand(const string& str){
	uint64_t res(0);
	for(uint i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			case 'T':res+=3;break;

			case 'a':res+=0;break;
			case 'c':res+=1;break;
			case 'g':res+=2;break;
			case 't':res+=3;break;
			default: return 0 ;break;
			//~ default:cout<<"bug:"<<str[i]<<"!"<<endl;exit(0);
			//~ default:res+=0;break;
		}
	}
	return (res);
}



uint64_t str2num(const string& str){
	return (min(str2numstrand(str),str2numstrand(revComp(str))));
}



char randNucle(char c='N'){
	switch (rand()%4){
	//~ switch (xs(seed)%4){
		case 0:
			if(c!='A'){
				return 'A';
			}
			return randNucle(c);
		case 1:
			if(c!='C'){
				return 'C';
			}
			return randNucle(c);
		case 2:
			if(c!='G'){
				return 'G';
			}
			return randNucle(c);
		case 3:
			if(c!='T'){
				return 'T';
			}
			return randNucle(c);
	}
	return randNucle(c);
}



string rand_sequence(uint64_t n){
	string res;
	for(uint32_t i(0);i<n;++i){
		res+=randNucle();
	}
	return res;
}



string compress_string(const std::string& str,int compressionlevel = Z_BEST_COMPRESSION)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, compressionlevel) != Z_OK)
        throw(std::runtime_error("deflateInit failed while compressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();           // set the z_stream's input

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // retrieve the compressed bytes blockwise
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, Z_FINISH);

        if (outstring.size() < zs.total_out) {
            // append the block to the output string
            outstring.append(outbuffer,zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
        throw(std::runtime_error(oss.str()));
    }

    return outstring;
}



/** Decompress an STL string using zlib and return the original data. */
std::string decompress_string(const std::string& str)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (inflateInit(&zs) != Z_OK)
        throw(std::runtime_error("inflateInit failed while decompressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // get the decompressed bytes blockwise using repeated calls to inflate
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = inflate(&zs, 0);

        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }

    } while (ret == Z_OK);

    inflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib decompression: (" << ret << ") "
            << zs.msg;
        throw(std::runtime_error(oss.str()));
    }

    return outstring;
}





