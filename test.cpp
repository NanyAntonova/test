#include <iostream>
#include <bitset>
#include <vector>
#include <string>


#include <seqan/seq_io.h>

int myordValue (char s)
{
	int x;
	switch (s){
		case 'a':
			x=1;
			break;
		case 'g':
			x=2;
			break;
		case 'c':
			x=3;
			break;
		case 't':
			x=4;
			break;
		default:
			x=0;
	}
	return x;
}
//int const maxlen = 3;





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
};


int half_hamming(const bitread &s1, const bitread &s2){
	size_t len = std::min<size_t>(s1.len, s2.len);	
	//std::bitset<s1.maxlen> res = ((s1.evenbit ^ s2.evenbit) | (s1.oddbit ^ s2.oddbit)) << (s1.maxlen-len);
	std::bitset<s1.maxlen> res = ((s1.evenbit ^ s2.evenbit) | (s1.oddbit ^ s2.oddbit)) & bitread::mask(s1.maxlen-len);
	return res.count();
}


int main ()
{
	std::vector <std::string> s;
	/*std::ifstream file("input.txt");
	  while(file) {
	  std::string str;
	  std::getline(file, str);
	  if(str!="")
	  s.push_back (str);
	  }*/	
	s.push_back  ("aag");
	s.push_back  ("aaa");
	s.push_back  ("gg");
	std::vector <bitread> bits;
	for (const auto &_ : s) {
		bits.push_back(bitread(_));
	}
	// masks ();
	int x = half_hamming (bits[0], bits[1]);
	std::cout << x << std::endl;

	return 0;

}
