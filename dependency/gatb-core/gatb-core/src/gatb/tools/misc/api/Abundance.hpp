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

/** \file Abundance.hpp
 *  \brief Abundance definition
 *  \date 01/03/2013
 *  \author edrezen
 */

/********************************************************************************/

#ifndef _GATB_CORE_TOOLS_MISC_ABUNDANCE_HPP_
#define _GATB_CORE_TOOLS_MISC_ABUNDANCE_HPP_

/********************************************************************************/

#include <sys/types.h>
#include <hdf5/hdf5.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** \brief Define a specific kind of array.
 *
 * This structure provides an array definition as an array attribute 'value' of type 'Type',
 * holding 'precision' items.
 *
 * This structure is used as a basis for the LargeInt implementation.
 *
 * \see math::LargeInt
 */
template<typename Type, int precision>
struct ArrayData
{
    Type value[precision];
};

/********************************************************************************/

/** \brief Define a type that associates a value and an abundance.
 *
 * We have often to count kmers, so we define a specific structure for this.
 *
 * The structure has two templates types:
 *  - Type : the type of items (likely kmer type)
 *  - Number : abundance associated to the item
 */
template<typename Type, typename Number=u_int16_t> struct Abundance
{
    /** Constructor.
     * \param[in] val : value of the item
     * \param[in] abund : abundance of the item.
     */
    Abundance (const Type& val, const Number& abund) : value(val), abundance(abund) {}
    Abundance (const Type& val) : value(val), abundance(0) {}
    Abundance () :  abundance(0) { Type zero; zero.setVal(0); value = zero;}

    /** Affectation operator.
     * \param[in] a : object to be copied.
     * \return the copied instance.
     */
    Abundance& operator=(const Abundance& a)
    {
        if (&a != this)  {  value = a.value;  abundance=a.abundance;  }
        return *this;
    }

    /** Get the abundance of the object
     * \return the abundance.
     */
    const Number& getAbundance() const { return abundance; }

    /** Get the value of the item
     * \return the value.
     */
    const Type&   getValue()     const { return value;     }

    /** Equality operator. alues and abundances must be equals
     * \param[in] other : object to be compared to
     * \return true if values and abundances are the same
     */
    bool operator== (const Abundance& other) const  {  return value == other.value && abundance == other.abundance;  }

    /** Creates a HDF5 type identifier for the [kmer,abundance] structure. This type will be used
     * for dumping Count instances in a HDF5 file (like SortingCount algorithm does).
     * \param[in] isCompound : tells whether the structure is compound (SHOULD BE OBSOLETE IN THE FUTURE)
     * \return the HDF5 identifier for the type. */
    static hid_t hdf5 (bool& isCompound)
    {
        hid_t abundanceType = H5T_NATIVE_UINT16;

             if (sizeof(Number)==1) { abundanceType = H5T_NATIVE_UINT8;  }
        else if (sizeof(Number)==2) { abundanceType = H5T_NATIVE_UINT16; }
        else if (sizeof(Number)==4) { abundanceType = H5T_NATIVE_UINT32; }
        else if (sizeof(Number)==8) { abundanceType = H5T_NATIVE_UINT64; }
        else { throw "Bad type size for Abundance HDF5 serialization";   }

        hid_t result = H5Tcreate (H5T_COMPOUND, sizeof(Abundance));
        H5Tinsert (result, "value",      HOFFSET(Abundance, value),     Type::hdf5(isCompound));
        H5Tinsert (result, "abundance",  HOFFSET(Abundance, abundance), abundanceType);

        isCompound = true;

        return result;
    }

    Type    value;
    Number  abundance;
};

/********************************************************************************/
}}}}
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_ABUNDANCE_HPP_ */
