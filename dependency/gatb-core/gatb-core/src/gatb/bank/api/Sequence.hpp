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

/** \file Sequence.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Definition of what a genomic sequence is.
 */

#ifndef _GATB_CORE_BANK_SEQUENCE_HPP_
#define _GATB_CORE_BANK_SEQUENCE_HPP_

/********************************************************************************/

#include <gatb/tools/misc/api/Data.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
/********************************************************************************/

/** \brief Structure holding genomic information.
 *
 * A sequence holds several data :
 *  - comment (as a text)
 *  - genomic data
 *  - quality information (for fastq format, empty in other cases).
 *
 *  The genomic data is hold in a tools::misc::Data attribute and is supposed to
 *  hold nucleotides.
 *
 *  Actually, the inner format may be of different kind (ASCII, INTEGER, BINARY) and
 *  depends on the type of the bank that provides Sequence objects. For instance:
 *      - a FASTA bank will provide Sequence instances whose data is in ASCII
 *      - a BINARY bank will provide Sequence instances whose data is in BINARY
 *
 *  The buffer holding the nucleotides is located in the tools::misc::Data attribute, so have a look
 *  there to have further details on where the buffer can be allocated. Note just here
 *  that the buffer could be stored in the Data object itself, or may be a reference to a
 *  buffer allocated in another place.
 *
 *  The class Sequence is closely related to the IBank interface.
 *
 *  Note that this class should not be instantiated directly by end users; it is more likely that end users will
 *  receive such objects through an iteration from a bank.
 *
 * Example of use:
 * \snippet bank16.cpp  snippet16_seq
 *
 *  \see IBank
 */
struct Sequence
{
    /** Constructor.
     * \param[in] encoding : encoding scheme of the genomic data of the sequence */
    Sequence (tools::misc::Data::Encoding_e encoding = tools::misc::Data::ASCII) : _data(encoding), _index(0)  {}

    /** Constructor. For testing mainly : allows to set the genomic data through an ascii representation.
     * For instance, one can provide "ACTTACGCAGAT" as argument of this constructor.
     * \param[in] seq : the genomic data as an ascii string */
    Sequence (char* seq) : _data(seq), _index(0)  {}

    /** Destructor. */
    virtual ~Sequence ()  { }

    /** \return description of the sequence */
    virtual const std::string& getComment ()  const  { return _comment; }

    /** \return description of the sequence until first space */
    virtual const std::string getCommentShort ()  const  { return _comment.substr(0, _comment.find(' ')); }

    /** \return quality of the sequence (set if the underlying bank is a fastq file). */
    virtual const std::string& getQuality ()  const  { return _quality; }
    
    /** \return the data as a Data structure. */
    virtual tools::misc::Data& getData () { return _data; }

    /** Return the raw buffer holding the genomic data. IMPORTANT : getting genomic data this way
     * implies that the user knows what is the underlying encoding scheme in order to decode it
     * (may be ASCII, INTEGER or BINARY)
     * \return buffer holding the genomic data as a raw buffer. */
    virtual char* getDataBuffer ()  const { return _data.getBuffer(); }

    /** \return number of nucleotides in the sequence. */
    virtual size_t getDataSize () const  { return _data.size(); }

    /** \return encoding scheme of the data. */
    virtual tools::misc::Data::Encoding_e getDataEncoding () const  { return _data.getEncoding(); }

    /** Return the index of the sequence. It may be the index of the sequence in the database that
     * holds the sequence.
     * \return index of the sequence. */
    virtual size_t getIndex () const  { return _index; }

    /** Set the genomic data as a reference on a Data object (more precisely on a range in this data).
     * This method may be used when one wants that the genomic data of the sequence points to an
     * already existing buffer of nucleotides, which means that the sequence doesn't allocate any
     * memory for storing the genomic data, it only relies on data stored somewhere else.
     * This is mainly a shortcut to the gatb::core::tools::misc::Data::setRef method.
     * \param[in] ref : the referred Data instance holding the genomic data
     * \param[in] offset : starting index in the referred data
     * \param[in] length : length of the genomic data of the current sequence. */
    void setDataRef (tools::misc::Data* ref, int offset, int length)  {  _data.setRef (ref, offset, length);  }

    /** Set the index of the sequence. Typically, it should be called by a IBank iterator that knows what is
     * the index of the currently iterated sequence.
     * \param[in] index : index of the sequence */
    void setIndex (size_t index)  { _index = index; }

    /** Get an ascii representation of the sequence. IMPORTANT ! this implementation supposes that the
     * format of the Data attribute is ASCII. No conversion is done in case of other formats.
     * \return the ascii representation of the sequence. */
    std::string toString () const { return std::string (this->getDataBuffer(), this->getDataSize()); }

    /** Set the comment of the sequence (likely to be called by a IBank iterator).
     * \param[in] cmt : comment of the sequence */
    void setComment (const std::string& cmt)  { _comment = cmt; }

    /** Set the quality string of the sequence (likely to be called by a fastq iterator).
     * \param[in] qual : quality string of the sequence. */
    void setQuality (const std::string& qual)  { _quality = qual; }

    /** Comment attribute (note: should be private with a setter and getter). */
    std::string _comment;

    /** Quality attribute (note: should be private with a setter and getter). */
    std::string _quality;

private:

    /** Object holding the genomic data of the sequence (ie a succession of nucleotides). */
    tools::misc::Data _data;

    /** Index of the sequence (likely to be set by a IBank iterator). */
    size_t _index;
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_SEQUENCE_HPP_ */
