#include <iostream>
#include <bitset>
#include <vector>
#include <unordered_map>
#include <string>
#include <seqan/seq_io.h>

struct bitread{
	static const size_t maxlen = 640;
	static const std::bitset<maxlen> * init_masks () {
		static std::bitset<maxlen> mask[maxlen];
		mask[0].set();
		for (size_t i=1; i<maxlen; i++){
			mask[i] = mask[i-1] >> 1;
		}
		return mask;
	}

	static const std::bitset<maxlen>& mask(size_t i) {
		static auto m = init_masks();
		return m[i];
	}

	public:
	size_t len;
	std::bitset<maxlen> evenbit;
	std::bitset<maxlen> oddbit;
	
	template<typename T>
	bitread(const T &s) : len{seqan::length(s)} {
		for (size_t i = 0; i < maxlen; i++) {
			std::bitset<2> tmp = seqan::ordValue (s[i]);
			evenbit[i] = tmp[0]; //TODO Rewrite??
			oddbit[i] = tmp[1];
		}
	}
	
	bitread() {}

	int dist_mask(const bitread &s){
                size_t mlen = std::min<size_t>(s.len, len);
                std::bitset<maxlen> res = ((s.evenbit ^ evenbit) | (s.oddbit ^ oddbit)) & bitread::mask(maxlen-mlen);
                return res.count();
        }
	
	bitread operator << (int i) {
		bitread tmp;
		tmp.evenbit = evenbit << i;
		tmp.oddbit = oddbit << i;
		tmp.len = len - i;
		return tmp;
	}
	
	bitread operator >> (int i) {
		bitread tmp;
                tmp.evenbit = evenbit >> i;
                tmp.oddbit = oddbit >> i;
                tmp.len = len + i;
                return tmp;
         }

};

struct AllShifts{
        std::unordered_map<int, bitread > shifts;
        static const size_t min_overlap_len = 300;
        AllShifts ( bitread b){
		int max_shift = b.len - min_overlap_len;  
		shifts[0] = b;
                for (size_t i = 1; i <= max_shift; i++){
                        shifts[i] = b >> i;
			shifts[-i] = b << i;
		}       
        }
};

using std::vector;
using seqan::Dna5String;
using seqan::CharString;

int main ()
{
	
	seqan::SeqFileIn seqFileIn_reads("merged_reads.fastq");
	
	vector< vector < int > > graph; // 

	vector<CharString> read_ids;
	vector<Dna5String> reads;
	readRecords(read_ids, reads, seqFileIn_reads);
	int tau = 10;
	int tmp_len = reads.size();// !!
	for (size_t j = 0; j < 10; j++){
		for (size_t i = 0; i < tmp_len; i++){
			reads.push_back(reads[i]);
		}
	}	
	std::vector <bitread> bits;
	
	for (const auto &_ : reads) {
		bits.push_back(bitread(_));
	}
	
	/*AllShifts as (bits[0]); просмотр сдвигов
	int t1 =  (as.shifts.size()-1)/2;
	int t2 = - t1;
	std::cout << as.shifts.size() << ' ' << t1 << ' ' << t2 << std::endl;	
	for (int k = t2; k <= t1; k++){
		std::cout << k << ' ' << as.shifts[k].len << ' ' << as.shifts[k].evenbit << std::endl ;
	}*/

	int num_of_reads = bits.size();
	std::cout << num_of_reads << std::endl;
	for (size_t i = 0; i < num_of_reads/*bits.size()*/; i++){
		AllShifts as (bits[i]);
		//int t1 =  (as.shifts.size()-1)/2;
                //int t2 = - t1; // без этого не работает???
                std::cout << i << ' ' << as.shifts.size() << ' ' <<  std::endl;
		//for (size_t j = 0; j < num_of_reads/*bits.size()*/; j++){
			//for (int k = t2; k <= t1; k++){
				//if (bits[j].dist_mask (as.shifts[k]) <= tau){
					//in graph: 
					//graph[j].push_back(i);
					//std::cout << j << ' ' << i << ' ' << k << std::endl;
					//break;
				//}
			//}
		//}
		as.shifts.clear();	
	}
	std::cout << "-----------------"  << std::endl;	
} 
