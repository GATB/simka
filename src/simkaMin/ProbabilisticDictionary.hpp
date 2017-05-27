/*
 * ProbabilisticDictionary.h
 *
 *  Created on: 25 mai 2017
 *      Author: gbenoit
 */

#ifndef SIMKA2_SRC_SIMKAMIN_PROBABILISTICDICTIONARY_HPP_
#define SIMKA2_SRC_SIMKAMIN_PROBABILISTICDICTIONARY_HPP_


typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf< u_int64_t, hasher_t > boophf_t;









/*
template <class T>
class bfile_iterator_first : public std::iterator<std::forward_iterator_tag, T>{

	public:


	bfile_iterator_first() : _is(nullptr), _pos(0) ,_inbuff (0), _cptread(0){
		_buffsize = 10000;
		_buffer = (T *) malloc(_buffsize*sizeof(T));
	}


	bfile_iterator_first(const bfile_iterator_first& cr){
		_buffsize = cr._buffsize;
		_pos = cr._pos;
		_is = cr._is;
		_buffer = (T *) malloc(_buffsize*sizeof(T));
		memcpy(_buffer,cr._buffer,_buffsize*sizeof(T) );
		_inbuff = cr._inbuff;
		_cptread = cr._cptread;
		_elem = cr._elem;
	}


	bfile_iterator_first(FILE* is): _is(is) , _pos(0) ,_inbuff (0), _cptread(0){
		_buffsize = 10000;
		_buffer = (T *) malloc(_buffsize*sizeof(T));
		int reso = fseek(_is,0,SEEK_SET);
		advance();
	}


	~bfile_iterator_first(){
		if(_buffer!=NULL)
			free(_buffer);
	}


	T const& operator*()  {  return _elem;  }


	bfile_iterator_first& operator++(){
		advance();
		return *this;
	}


	friend bool operator==(bfile_iterator_first const& lhs, bfile_iterator_first const& rhs){
		if (!lhs._is || !rhs._is)  {  if (!lhs._is && !rhs._is) {  return true; } else {  return false;  } }
		assert(lhs._is == rhs._is);
		return rhs._pos == lhs._pos;
	}


	friend bool operator!=(bfile_iterator_first const& lhs, bfile_iterator_first const& rhs)  {  return !(lhs == rhs);  }


private:
	void advance()
	{
		++_pos;

		if(_cptread >= _inbuff){
			int res = fread(_buffer,sizeof(T),_buffsize,_is);
			_inbuff = res; _cptread = 0;

			if(res == 0){
				_is = nullptr;
				_pos = 0;
				return;
			}
		}

		_elem = _buffer[_cptread];
		++_cptread;
		++_cptread;
	}
	T _elem;
	FILE * _is;
	unsigned long _pos;

	T * _buffer; // for buffered read
	int _inbuff, _cptread;
	int _buffsize;
};

template <class T>
class file_binary_first{

	public:


	file_binary_first(const char* filename){
		_is = fopen(filename, "rb");
		if (!_is) {
			throw std::invalid_argument("Error opening " + std::string(filename));
		}
	}


	~file_binary_first(){
		fclose(_is);
	}


	bfile_iterator_first<T> begin() const{
		return bfile_iterator_first<T>(_is);
	}

	bfile_iterator_first<T> end() const {return bfile_iterator_first<T>(); }

	size_t size () const{
		return 0;
	}//todo ?

private:
	FILE * _is;
};
*/


// iterator from disk file of T with buffered read
template <class T>
class bfile_iterator : public std::iterator<std::forward_iterator_tag, T>{
public:

	bfile_iterator()
: _is(nullptr)
, _pos(0) ,_inbuff (0), _cptread(0)
{
		_buffsize = 10000;
		_buffer = (T *) malloc(_buffsize*sizeof(T));
}

	bfile_iterator(const bfile_iterator& cr)
	{
		_buffsize = cr._buffsize;
		_pos = cr._pos;
		_is = cr._is;
		_buffer = (T *) malloc(_buffsize*sizeof(T));
		memcpy(_buffer,cr._buffer,_buffsize*sizeof(T) );
		_inbuff = cr._inbuff;
		_cptread = cr._cptread;
		_elem = cr._elem;
	}

	bfile_iterator(FILE* is): _is(is) , _pos(0) ,_inbuff (0), _cptread(0)
	{
		_buffsize = 10000;
		_buffer = (T *) malloc(_buffsize*sizeof(T));
		int reso = fseek(_is,0,SEEK_SET);
		advance();
	}

	~bfile_iterator()
	{
		if(_buffer!=NULL)
			free(_buffer);
	}


	T const& operator*()  {  return _elem;  }

	bfile_iterator& operator++()
    				{
		advance();
		return *this;
    				}

	friend bool operator==(bfile_iterator const& lhs, bfile_iterator const& rhs)
    				{
		if (!lhs._is || !rhs._is)  {  if (!lhs._is && !rhs._is) {  return true; } else {  return false;  } }
		assert(lhs._is == rhs._is);
		return rhs._pos == lhs._pos;
    				}

	friend bool operator!=(bfile_iterator const& lhs, bfile_iterator const& rhs)  {  return !(lhs == rhs);  }
private:
	void advance()
	{
		++_pos;

		if(_cptread >= _inbuff)
		{
			int res = fread(_buffer,sizeof(T),_buffsize,_is);
			_inbuff = res; _cptread = 0;

			if(res == 0)
			{
				_is = nullptr;
				_pos = 0;
				return;
			}
		}

		_elem = _buffer[_cptread];
		++_cptread;
	}
	T _elem;
	FILE * _is;
	unsigned long _pos;

	T * _buffer; // for buffered read
	int _inbuff, _cptread;
	int _buffsize;
};




template <class T>
class file_binary{
public:

	file_binary(const char* filename)
{
		_is = fopen(filename, "rb");
		if (!_is) {
			throw std::invalid_argument("Error opening " + std::string(filename));
		}
}

	~file_binary()
	{
		fclose(_is);
	}

	bfile_iterator<T> begin() const
    				{
		return bfile_iterator<T>(_is);
    				}

	bfile_iterator<T> end() const {return bfile_iterator<T>(); }

	size_t        size () const  {  return 0;  }//todo ?

private:
	FILE * _is;
};










//class ProbabilisticDictionary {
//public:
//	ProbabilisticDictionary();
//	virtual ~ProbabilisticDictionary();
//};




































template <class FingerprintType>
class ProbabilisticDictionary
{

	//typedef u_int16_t FingerprintType;
private:
	vector<FingerprintType>* _fingerprints;
	boophf_t * _bphf;
	int _fingerprint_size;
	uint64_t _fingerprint_range;


public:
	u_int64_t _nbElements;


	ProbabilisticDictionary(u_int64_t nbElements, file_binary<u_int64_t>& itKeys, int nbCores){
		_nbElements = nbElements;
		_fingerprint_size = sizeof(FingerprintType)*8;
		_fingerprint_range = (uint64_t)1<<_fingerprint_size;

		_fingerprints = new vector<FingerprintType>(nbElements);

		createMPHF(nbCores, itKeys);
		createFingerprints(itKeys);
	}

	//bool contains(u_int64_t key){
	//	u_int64_t index = _bphf->lookup(key);
	//	if(index == ULLONG_MAX) return false;
	//	if(_fingerprint_size>0)
	//	return _prob_set->exists(index, key);
	//	return true;
	//}

	u_int64_t getIndex(u_int64_t key, bool &exists){
		u_int64_t index = this->_bphf->lookup(key);
		if(index == ULLONG_MAX){
			exists = false;
			return 0;
		}

        FingerprintType fingerprint = korenXor(key)%_fingerprint_range;
        if(fingerprint != _fingerprints->at(index)){
			exists = false;
        	return 0;
        }

		exists = true;
		return index;
		//return this->_values.get_i(index);
	}



private:

	void createMPHF(size_t nbCores, file_binary<u_int64_t>& itKeys){
		double gammaFactor = 2;
		this->_bphf = new boomphf::mphf<u_int64_t,hasher_t>(_nbElements, itKeys, nbCores, gammaFactor);
	}

	void createFingerprints(file_binary<u_int64_t>& itKeys){
		for (bfile_iterator<u_int64_t> it = itKeys.begin(); it != itKeys.end(); ++it){
			u_int64_t key = *it;
			u_int64_t index = this->_bphf->lookup(key);
			FingerprintType fingerprint = korenXor(key)%_fingerprint_range;
			//cout << index << ": " << key << endl;
			_fingerprints->at(index) = fingerprint;
			//if (this->_fingerprint_size>0){
			//	this->_prob_set->add(index, std::get<0>(key_value));
			//}
			//this->_values.set_i(index, std::get<1>(key_value));
		}
		//for(u_int64_t& key: itKeys){
		//}

		//void add(const uint64_t i, const uint64_t key){
		//        uint64_t fingerprint = korenXor(key)%_fingerprint_range;
		//       _bas.set_i(i,fingerprint);
		//}
	}

	inline uint64_t korenXor(uint64_t x)const{
		x ^= (x << 21);
		x ^= (x >> 35);
		x ^= (x << 4);
		return x;
	}
};



#endif /* SIMKA2_SRC_SIMKAMIN_PROBABILISTICDICTIONARY_H_ */
