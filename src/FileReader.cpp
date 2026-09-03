#include "FileReader.hpp"
#include "linearraccess/seq_utils.hpp"

#include <cctype>
#include <fstream>
#include <iostream>

// trim \r, \n in end of the line
inline void trim_end(string &s){
	while(!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
}

// read sequences from file_name
bool FileReader::read(const string file_name, vector<string> &seq, vector<string> &seq_name){
	seq.clear();
	seq_name.clear();
	ifstream ifs(file_name);
	if(!ifs){
		cout << "Error: cannot open input file: " << file_name << endl;
		return false;
	}

	string line;
	while(getline(ifs, line)){
		trim_end(line);
		if(line.empty()) continue;
		if(line[0] == '>'){
			if(!seq.empty() && seq.back().empty()){
				cout << "Error: empty FASTA sequence: " << seq_name.back() << endl;
				return false;
			}
			// new sequence
			seq_name.push_back(line.substr(1));
			seq.push_back("");
		}else{
			if(seq.empty()){
				cout << "Error: FASTA sequence data appears before a header" << endl;
				return false;
			}
			// append sequence
			seq.back() += line;
		}
	}

	ifs.close();
	if(seq.empty()){
		cout << "Error: no FASTA records found in: " << file_name << endl;
		return false;
	}
	if(seq.back().empty()){
		cout << "Error: empty FASTA sequence: " << seq_name.back() << endl;
		return false;
	}
	for(size_t idx = 0; idx < seq.size(); idx++){
		try{
			seq[idx] = lcr::seq::normalize_sequence(seq[idx]);
		}catch(const std::invalid_argument& error){
			cout << "Error: invalid FASTA sequence " << seq_name[idx]
			     << ": " << error.what() << endl;
			return false;
		}
	}
	return true;
}
