#ifndef UTILS
#define UTILS



#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <unordered_map>



using namespace std;



#define kmer uint64_t
#define minimizer uint8_t



uint32_t revhash ( uint32_t x );
uint32_t unrevhash ( uint32_t x );
uint64_t revhash64 ( uint64_t x ) ;
uint64_t unrevhash64 ( uint64_t x ) ;
uint64_t universal_hash ( uint64_t x, uint32_t k);
uint64_t str2num(const string& str);
char randNucle(char c);
string rand_sequence(uint64_t n);
string intToString(uint64_t n);
string compress_string(const string& str,int compressionlevel);
string decompress_string(const string& str);
double parse_line_number(string& str);
void clean(string& str);



#endif
