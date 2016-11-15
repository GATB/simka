/*
 * SimkaCommon.hpp
 *
 *  Created on: 14 nov. 2016
 *      Author: gbenoit
 */

#ifndef GATB_SIMKA_SRC_CORE_SIMKACOMMON_HPP_
#define GATB_SIMKA_SRC_CORE_SIMKACOMMON_HPP_


static void simka2_writeString(const string& s, ofstream& file){
	u_int16_t size = s.size();
	file.write((char const*)(&size), sizeof(size));
	file.write(s.c_str(), size);
}

static void simka2_readString(string& s, ifstream& file){
	u_int16_t size;
	file.read((char*)(&size), sizeof(size));
	std::vector<char> buffer(size);
	file.read(&buffer[0], buffer.size());
	s.assign(buffer.begin(), buffer.end());
	//return string linkedDatasetID( buffer.begin(), buffer.end() );
}



#endif /* GATB_SIMKA_SRC_CORE_SIMKACOMMON_HPP_ */
