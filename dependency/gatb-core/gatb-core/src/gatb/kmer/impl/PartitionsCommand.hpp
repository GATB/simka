/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _PARTITIONSCOMMAND__HPP_
#define _PARTITIONSCOMMAND__HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>

#include <gatb/kmer/impl/PartiInfo.hpp>
#include <gatb/kmer/api/ICountProcessor.hpp>

#include <gatb/tools/collections/api/Iterable.hpp>

#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/misc/api/Enums.hpp>
#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Pool.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>

#include <queue>
#include <limits>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** \brief Class that counts the number of occurrences of a kmer in several banks.
 *
 * This class manages a vector of occurrences for a kmer, each index of the vector
 * holding the occurrences number for a bank.
 *
 * It is used during the process of kmers counting and eases the work for counting
 * kmers per bank.
 */
class CounterBuilder
{
public:

    /** Constructor.
     * \param[in] nbBanks : number of banks parsed during kmer counting.
     */
    CounterBuilder (size_t nbBanks=1)  :  _abundancePerBank(nbBanks)  {}

    /** Get the number of banks.
     * \return the number of banks. */
    size_t size() const  { return _abundancePerBank.size(); }

    /** Initialization of the counting for the current kmer. This method should be called
     * when a kmer is seen for the first time.
     * \param[in] idxBank : bank index where the new current kmer has been found. */
    void init (size_t idxBank=0)
    {
        for (size_t k=0; k<_abundancePerBank.size(); k++)  { _abundancePerBank[k]=0; }
        _abundancePerBank [idxBank]= 1;
    }

    /** Increase the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void increase (size_t idxBank=0)  {  _abundancePerBank [idxBank] ++;  }

    /** Set the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void set (CountNumber val, size_t idxBank=0)  {  _abundancePerBank [idxBank] = val;  }

    /** Get the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank
     * \return the abundance of the current kmer for the given bank. */
    CountNumber operator[] (size_t idxBank) const  { return _abundancePerBank[idxBank]; }

    /** */
    const CountVector& get () const { return _abundancePerBank; }

private:
    CountVector _abundancePerBank;
};

/********************************************************************************/
template<size_t span>
class PartitionsCommand : public gatb::core::tools::dp::ICommand, public system::SmartPointer
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type           Type;
    typedef typename Kmer<span>::Count          Count;
    typedef ICountProcessor<span> CountProcessor;

    /** Constructor. */
    PartitionsCommand (
        gatb::core::tools::collections::Iterable<Type>& partition,
        CountProcessor*                                 processor,
        size_t                                          cacheSize,
        gatb::core::tools::dp::IteratorListener*        progress,
        tools::misc::impl::TimeInfo&                    timeInfo,
        PartiInfo<5>&                                   pInfo,
        int                                             passi,
		int                                             parti,
		size_t                                          nbCores,
		size_t                                          kmerSize,
		gatb::core::tools::misc::impl::MemAllocator&    pool
    );

    /** Destructor. */
    ~PartitionsCommand();

    /** Get the class name (for statistics). */
    virtual const char* getName() const = 0;

protected:
    gatb::core::tools::collections::Iterable<Type>&         _partition;
    gatb::core::tools::dp::IteratorListener*                _progress;
	PartiInfo<5>&                                           _pInfo;
    int                                                     _pass_num;
	int                                                     _parti_num;
    size_t                                                  _nbCores;
	size_t                                                  _kmerSize;
    size_t                                                  _cacheSize;
	gatb::core::tools::misc::impl::MemAllocator&            _pool;
	
    void insert (const Type& kmer, const CounterBuilder& count);

    tools::misc::impl::TimeInfo& _globalTimeInfo;
    tools::misc::impl::TimeInfo  _timeInfo;

    CountProcessor* _processor;
    void setProcessor (CountProcessor* processor)  { SP_SETATTR(processor); }
};

/********************************************************************************/
/** */
template<size_t span>
class PartitionsByHashCommand : public PartitionsCommand<span>
{
public:

    /** Shortcut. */ /* R: don't know how to avoid this code duplication => R1: I'm afraid it's not possible. */
    typedef typename Kmer<span>::Type           Type;
    typedef typename Kmer<span>::Count          Count;
    typedef ICountProcessor<span> CountProcessor;

    /** Constructor. */
    PartitionsByHashCommand (
        gatb::core::tools::collections::Iterable<Type>& partition,
        CountProcessor*                                 processor,
        size_t                                          cacheSize,
        gatb::core::tools::dp::IteratorListener*        progress,
        tools::misc::impl::TimeInfo&                    timeInfo,
        PartiInfo<5>&                                   pInfo,
        int                                             passi,
        int                                             parti,
        size_t                                          nbCores,
        size_t                                          kmerSize,
        gatb::core::tools::misc::impl::MemAllocator&    pool,
        u_int64_t                                       hashMemory
    );

    /** Get the class name (for statistics). */
    const char* getName() const { return "hash"; }

    /** */
    void execute ();

private:
    u_int64_t _hashMemory;
};
		
/********************************************************************************/
/** */
template<size_t span>
class PartitionsByVectorCommand : public PartitionsCommand<span>
{
public:

    /** Shortcut. */ /* R: don't know how to avoid this code duplication */
    typedef typename Kmer<span>::Type           Type;
    typedef typename Kmer<span>::Count          Count;
    typedef ICountProcessor<span> CountProcessor;

    static const size_t KX = 4 ;

private:
    //used for the priority queue
    typedef std::pair<int, Type> kxp; //id pointer in vec_pointer , value
    struct kxpcomp { bool operator() (kxp l,kxp r) { return ((r.second) < (l.second)); } } ;

public:
    /** Constructor. */
    PartitionsByVectorCommand (
        gatb::core::tools::collections::Iterable<Type>& partition,
        CountProcessor*                                 processor,
        size_t                                          cacheSize,
        gatb::core::tools::dp::IteratorListener*        progress,
        tools::misc::impl::TimeInfo&                    timeInfo,
        PartiInfo<5>&                                   pInfo,
        int                                             passi,
        int                                             parti,
        size_t                                          nbCores,
        size_t                                          kmerSize,
        gatb::core::tools::misc::impl::MemAllocator&    pool,
        std::vector<size_t>&                            offsets
    );

    /** Destructor. */
    ~PartitionsByVectorCommand ();

    /** Get the class name (for statistics). */
    const char* getName() const { return "vector"; }

    /** */
    void execute ();

private:

    Type**     	       _radix_kmers;
    bank::BankIdType** _bankIdMatrix;
	uint64_t*          _radix_sizes;
	uint64_t*          _r_idx;

    tools::dp::IDispatcher* _dispatcher;

	void executeRead   ();
    void executeSort   ();
    void executeDump   ();

    std::vector<size_t> _nbItemsPerBankPerPart;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif
