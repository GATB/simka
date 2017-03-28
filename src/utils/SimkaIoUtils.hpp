/*
 * SimkaIoUtils.h
 *
 *  Created on: 18 nov. 2016
 *      Author: gbenoit
 */

#ifndef GATB_SIMKA_SRC_UTILS_SIMKAIOUTILS_H_
#define GATB_SIMKA_SRC_UTILS_SIMKAIOUTILS_H_


#include "../../thirdparty/zlib.h"



class SimkaIoUtils {
public:
	//SimkaIoUtils();
	//virtual ~SimkaIoUtils();

	static string getDatasetID(const string& kmerSpectrumDir){
		string datasetID = System::file().getBaseName(kmerSpectrumDir);
		//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());
		return datasetID;
	}



	static u_int64_t simka2_getFileSize(const string& filename){
		std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
		u_int64_t size = in.tellg();
		in.close();
		return size;
	}
	/*
	u_int64_t simka2_getDatasetSize(const string& kmerSpectrumDir){
		u_int64_t datasetSize = 0;
		for(size_t partitionID=0; partitionID<SIMKA2_NB_PARTITIONS; partitionID++){
			string partFilename = kmerSpectrumDir + "/" + Stringify::format("%i", partitionID) + ".gz";
			datasetSize += simka2_getFileSize(partFilename);
		}
		return datasetSize;
	}
	*/

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

	static void simka2_writeDatasetInfo(ofstream& file, const string& datasetID, u_int64_t nbReads, u_int64_t nbDistinctKmers, u_int64_t nbKmers, u_int64_t chord_N2){
	    simka2_writeString(datasetID, file);
	    file.write((char const*)(&nbReads), sizeof(nbReads));
	    file.write((char const*)(&nbDistinctKmers), sizeof(nbDistinctKmers));
	    file.write((char const*)(&nbKmers), sizeof(nbKmers));
	    file.write((char const*)(&chord_N2), sizeof(chord_N2));
	}

	static void simka2_readDatasetInfo(ifstream& file, string& datasetID, u_int64_t& nbReads, u_int64_t& nbDistinctKmers, u_int64_t& nbKmers, u_int64_t& chord_N2){
		simka2_readString(datasetID, file);
		file.read((char *)(&nbReads), sizeof(nbReads));
		file.read((char *)(&nbDistinctKmers), sizeof(nbDistinctKmers));
		file.read((char *)(&nbKmers), sizeof(nbKmers));
		file.read((char *)(&chord_N2), sizeof(chord_N2));
	}

	static void simka2_transferDatasetInfo(ifstream& source, ofstream& dest){
		string datasetID;
		u_int64_t nbReads;
		u_int64_t nbDistinctKmers;
		u_int64_t nbKmers;
		u_int64_t chord_N2;

		simka2_readString(datasetID, source);
		source.read((char *)(&nbReads), sizeof(nbReads));
		source.read((char *)(&nbDistinctKmers), sizeof(nbDistinctKmers));
		source.read((char *)(&nbKmers), sizeof(nbKmers));
		source.read((char *)(&chord_N2), sizeof(chord_N2));

	    simka2_writeString(datasetID, dest);
	    dest.write((char const*)(&nbReads), sizeof(nbReads));
	    dest.write((char const*)(&nbDistinctKmers), sizeof(nbDistinctKmers));
	    dest.write((char const*)(&nbKmers), sizeof(nbKmers));
	    dest.write((char const*)(&chord_N2), sizeof(chord_N2));
	}

};








/* \brief Bag implementation for file
 */
template <typename Item> class BagGzFileSimka : public Bag<Item>, public SmartPointer
{
public:

    /** Constructor. */
	BagGzFileSimka (const std::string& filename) : _filename(filename), _gzfile(0)
    {
        /** We first erase the file. */
        System::file().remove (filename);

        /** We get a handle on the file. */
        _gzfile =   gzopen(filename.c_str(),"wb1"); // (gzFile *)
    }

    /** Destructor. */
    ~BagGzFileSimka ()
    {
        if (_gzfile)  {  gzclose(_gzfile); }
    }

    /** */
    const std::string& getName () const { return _filename; }

    /**  \copydoc Bag::insert */
    void insert (const Item& item)  { gzwrite(_gzfile,&item,sizeof(Item)); }

    void insert (const std::vector<Item>& items, size_t length)
    {
        if (length == 0)  { length = items.size(); }
        //_file->fwrite (items.data(), sizeof(Item), length);
        gzwrite(_gzfile,items.data(),sizeof(Item)*length);
    }

    void write(u_int64_t item, size_t length){
        gzwrite(_gzfile, (char const*)(&item), length);
    }
    /** */
    void insert (const Item* items, size_t length)
    {
        gzwrite(_gzfile,items,sizeof(Item)*length);
    }

    /**  \copydoc Bag::flush */
    void flush ()  {
        gzflush(_gzfile,Z_FINISH); // not sure if necessary
    }

private:
    std::string _filename;
    IFile* _file;
    gzFile  _gzfile;
};




template <class Item> class IteratorGzFileSimka : public Iterator<Item>
{
public:

    /** Constructor. */
	IteratorGzFileSimka () : _gzfile(0), _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(0), _isDone(true) {}

	IteratorGzFileSimka (const IteratorGzFileSimka& it):
    _filename(it._filename), _gzfile(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(it._cacheItemsNb), _isDone(true)
    {
        _gzfile =   gzopen(_filename.c_str(),"rb");
        _buffer  = (Item*) MALLOC (sizeof(Item) * _cacheItemsNb);
    }

    /** Constructor. */
	IteratorGzFileSimka (const std::string& filename, size_t cacheItemsNb=10000) :
    _filename(filename), _gzfile(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(cacheItemsNb), _isDone(true)

    {
        _gzfile =   gzopen(_filename.c_str(),"rb");
        _buffer  = (Item*) MALLOC (sizeof(Item) * _cacheItemsNb);
    }

    /** Destructor. */
    ~IteratorGzFileSimka ()
    {
        if (_gzfile)  { gzclose(_gzfile);   }
        if (_buffer) { FREE (_buffer); }
    }

    /** Affectation. */
    IteratorGzFileSimka& operator= (const IteratorGzFileSimka& it)
    {
        if (this != &it)
        {
            if (_gzfile)    {  gzclose(_gzfile);  }
            if (_buffer)  { FREE(_buffer); }

            _filename     = it._filename;
            _cpt_buffer   = it._cpt_buffer;
            _idx          = it._idx;
            _cacheItemsNb = it._cacheItemsNb;
            _isDone       = it._isDone;

            _gzfile =   gzopen(_filename.c_str(),"r");
            _buffer  = (Item*) MALLOC (sizeof(Item) * it._cacheItemsNb);
        }
        return *this;
    }

    /** \copydoc dp::Iterator::first */
    void first()
    {
        gzseek(_gzfile,0,SEEK_SET);
        _cpt_buffer = 0;
        _idx        = 0;
        _isDone     = false;
        next ();
    }

    /** \copydoc dp::Iterator::next */
    void next()
    {
        if (_cpt_buffer==0)
        {
            _idx = 0;
            _cpt_buffer = (gzread(_gzfile, _buffer, sizeof(Item)*_cacheItemsNb)  /  sizeof(Item))  ;  // gzread returns number of bytes uncompressed returned in _buffer
            //printf("refreshing buffer from gzread(), cptbuffer now: %d\n",_cpt_buffer);
            if (_cpt_buffer < 0)
            { // On other errors, gzread() shall return a value less than 0 and and applications may examine the cause using gzerror().
                // FIXME: more tests are needed but on my system (R: linux64), gzread returns the fixed number of items, then a lower number of items (as it reaches the end), then it returns -1 instead of 0 and prints a "data error" but at this point it's fine to just return
                int err;
                fprintf(stderr, "gzread error: %s\n", gzerror(_gzfile, &err));
            }
            if (_cpt_buffer<=0)  { _isDone = true;  return; }
        }
        *(this->_item) =  _buffer[_idx];
        _cpt_buffer --;
        _idx ++;
    }

    /** \copydoc dp::Iterator::isDone */
    bool isDone()  { return _isDone; }

    /** \copydoc dp::Iterator::item */
    Item& item ()  { return *(this->_item); }

    /** */
    size_t fill (std::vector<Item>& vec, size_t len=0)
    {
        if (len==0)  { len = vec.size(); }

        return (gzread(_gzfile,vec.data(), sizeof(Item)*len)  /  sizeof(Item));

    }

private:
    std::string     _filename;
    gzFile  _gzfile;
    Item*           _buffer;
    int             _cpt_buffer;
    int             _idx;
    size_t          _cacheItemsNb;
    bool            _isDone;
};




template <class Item> class IterableGzFileSimka : public Iterable<Item>, public virtual SmartPointer
{
public:

    /** */
	IterableGzFileSimka (const std::string& filename, size_t cacheItemsNb=10000)
    :   _filename(filename), _cacheItemsNb (cacheItemsNb)  {}

    /** */
    ~IterableGzFileSimka () {}

    /** */
    Iterator<Item>* iterator ()  { return new IteratorGzFileSimka<Item> (_filename, _cacheItemsNb); }

    /** */
    int64_t getNbItems ()   {  return -1; } // does not know value

    int64_t estimateNbItems ()   {  return 3* (System::file().getSize(_filename) / sizeof(Item)); }

private:
    std::string     _filename;
    size_t          _cacheItemsNb;
};



#endif /* GATB_SIMKA_SRC_UTILS_SIMKAIOUTILS_H_ */
