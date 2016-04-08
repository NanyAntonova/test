#include <iostream>
#include <bitset>
#include <vector>
#include <string>
#include <seqan/seq_io.h>
#include <chrono>

struct bitread{
	static const size_t maxlen = 640;
	static const std::bitset<maxlen> * init_masks () {
		//int const maxlen = 3;
		//std::bitset<maxlen> mask[maxlen];
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

	int half_hamming_shift(const bitread &s){
        	size_t mlen = std::min<size_t>(s.len, len);
        	std::bitset<maxlen> res = ((s.evenbit ^ evenbit) | (s.oddbit ^ oddbit)) << (maxlen-mlen);
        	return res.count();
	}

	int half_hamming_mask(const bitread &s){
                size_t mlen = std::min<size_t>(s.len, len);
                std::bitset<maxlen> res = ((s.evenbit ^ evenbit) | (s.oddbit ^ oddbit)) & bitread::mask(maxlen-mlen);
                return res.count();
        }
	

};

template<typename Ts1, typename Ts2>
int old_half_hamming(const Ts1 &s1, const Ts2 &s2) {
    int len1 = length(s1), len2 = length(s2); 
    size_t len = std::min<size_t>(len1, len2);

    int res = 0;
    for (size_t i = 0; i < len; ++i) {
        //res += (s1[i] == s2[i]) ? 0 : 1;
    	if (s1[i] != s2 [i])
		res++;
	}
    return res;
}

using std::vector;
using seqan::Dna5String;
using seqan::CharString;

int main ()
{
	
	seqan::SeqFileIn seqFileIn_reads("merged_reads.fastq");

	vector<CharString> read_ids;
	vector<Dna5String> reads;
	readRecords(read_ids, reads, seqFileIn_reads);
	
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
	int dist1 = 0, dist2 = 0, dist3 = 0;
	auto start_time = std::chrono::steady_clock::now(); // new
	for (size_t i = 0; i < bits.size(); i++) {
		for (size_t j = 0/*i*/; j < bits.size(); j++){
			//std::cout << bits[i].half_hamming (bits[j]) << ' ';
			dist1+=(bits[i].half_hamming_shift (bits[j]));
		}
		
	}
	auto end_time = std::chrono::steady_clock::now();
	auto t = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
	std::cout << t.count() << std::endl; 
	
	auto start_time1 = std::chrono::steady_clock::now(); // new mask
        for (size_t i = 0; i < bits.size(); i++) {
                for (size_t j = 0/*i*/; j < bits.size(); j++){
                        //std::cout << bits[i].half_hamming (bits[j]) << ' ';
                        dist2+=(bits[i].half_hamming_mask (bits[j]));
                }

        }
        auto end_time1 = std::chrono::steady_clock::now();
        auto t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time1 - start_time1);
        std::cout << t1.count() << std::endl;

	auto start_time2 = std::chrono::steady_clock::now(); // old
        for (size_t i = 0; i < bits.size(); i++) {
                for (size_t j = 0/*i*/; j < bits.size(); j++){
                        //std::cout << old_half_hamming (reads[i],reads[j]) << ' ';
        		dist3+=(old_half_hamming (reads[i],reads[j])); 
	       }
                
        }
        auto end_time2 = std::chrono::steady_clock::now();
        auto t2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time2 - start_time2);
	std::cout << t2.count() << std::endl;
	if (dist1 == dist3)
		 std::cout << "ok" << reads.size() << std::endl;
	return 0;

}
