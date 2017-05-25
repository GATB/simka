/*
 * SimkaMinUtils.h
 *
 *  Created on: 25 mai 2017
 *      Author: gbenoit
 */

#ifndef SIMKA2_SRC_SIMKAMIN_SIMKAMINUTILS_H_
#define SIMKA2_SRC_SIMKAMIN_SIMKAMINUTILS_H_

























//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
template <class Item> class SimkaInputIterator : public Iterator<Item>
{
public:

	/** Constructor.
	* \param[in] ref : the referred iterator
	* \param[in] initRef : will call 'first' on the reference if true
	*/
	SimkaInputIterator(Iterator<Item>* refs, size_t nbBanks, u_int64_t maxReads)
	:  _mainref(0) {

		setMainref(refs);
		_ref = _mainref->getComposition()[0];
		_isDone = false;
		_nbDatasets = nbBanks;
		_nbBanks = _mainref->getComposition().size() / _nbDatasets;
		_maxReads = maxReads;
		_nbReadProcessed = 0;
		_currentBank = 0;
		_currentInternalBank = 0;
		_currentDataset = 0;

	}


    bool isFinished(){
        if(_currentDataset == _nbDatasets){
                _isDone = true;
                return true;
        }
        return false;
    }

	void nextDataset(){
		_currentDataset += 1;

		if(isFinished()) return;

		_currentBank = _currentDataset * _nbBanks;

		_currentInternalBank = 0;
		_nbReadProcessed = 0;

		if(isFinished()) return;

		_ref = _mainref->getComposition()[_currentBank];
		_isDone = false;
		first();
		//nextBank();
	}

	void nextBank(){
		//cout << "next bank" << endl;
		//cout << "next bank "<< endl;
		_currentInternalBank += 1;
		if(_currentInternalBank == _nbBanks){
			nextDataset();
		}
		else{
			_isDone = false;
			_currentBank += 1;
			_ref = _mainref->getComposition()[_currentBank];
			first();
		}
	}

    void first()
    {

        _ref->first();

       // while (!_ref->isDone() && _filter(_ref->item())==false)
         //       _ref->next();

        _isDone = _ref->isDone();

        if(!_isDone) *(this->_item) = _ref->item();

    }

	void next(){


		if(isFinished()){
			_isDone = true;
			return;
		}

		//cout << "haha" << endl;

		_ref->next();
		//while (!_ref->isDone() && _filter(_ref->item())==false) _ref->next();

		_isDone = _ref->isDone();

		//cout << "haha" << endl;
		//if(!_isDone){
			//cout << _currentBank << "  " << _isDone << endl;

		//}

		//cout << _nbReadProcessed << "  " << _currentBank << "    " << _nbBanks << "   " << _maxReads << endl;


		if(_isDone){
			if(isFinished()){
				//cout << _nbReadProcessed << endl;
				return;
			}
			else{
				//cout << _nbReadProcessed << endl;
				nextBank();
				if(isFinished()){
					//cout << _nbReadProcessed << endl;
					return;
				}
			}
		}
		else{
			*(this->_item) = _ref->item();
			_nbReadProcessed += 1;
		}

		if(_maxReads && _nbReadProcessed >= _maxReads){
			if(isFinished())
				return;
			else
				nextDataset();
		}

	}

    /** \copydoc  Iterator::isDone */
    bool isDone()  {  return _isDone;  }

    /** \copydoc  Iterator::item */
    Item& item ()  {  return *(this->_item);  }


private:

    bool            _isDone;
    size_t _currentBank;
    //vector<Iterator<Item>* > _refs;
    Iterator<Item>* _ref;
    size_t _nbBanks;
    u_int64_t _maxReads;
    //Filter _filter;
    u_int64_t _nbReadProcessed;
    size_t _currentInternalBank;
	size_t _currentDataset;
	size_t _nbDatasets;

    Iterator<Item>* _mainref;
    void setMainref (Iterator<Item>* mainref)  { SP_SETATTR(mainref); }
};




#endif /* SIMKA2_SRC_SIMKAMIN_SIMKAMINUTILS_H_ */
