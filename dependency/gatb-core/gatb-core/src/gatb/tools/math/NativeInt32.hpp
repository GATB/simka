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

/** \file NativeInt32.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int32_t type
 */

#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_32_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_32_HPP_

/********************************************************************************/

#include <iostream>
#include <gatb/system/api/types.hpp>
#include <gatb/tools/misc/api/Abundance.hpp>
#include <hdf5/hdf5.h>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
/** \brief Math package */
namespace math  {
/********************************************************************************/

class NativeInt32 : private misc::ArrayData<u_int16_t, 1>
{
public:

    typedef ArrayData<u_int16_t, 1> POD;

    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    NativeInt32 (const u_int8_t& c=0)  {  value[0] = c;  }

    static const char* getName ()  { return "NativeInt32"; }

    static const size_t getSize ()  { return 8*sizeof(u_int16_t); }

    NativeInt32 operator+  (const NativeInt32& other)   const   {  return value[0] + other.value[0];  }
    NativeInt32 operator-  (const NativeInt32& other)   const   {  return value[0] - other.value[0];  }
    NativeInt32 operator|  (const NativeInt32& other)   const   {  return value[0] | other.value[0];  }
    NativeInt32 operator^  (const NativeInt32& other)   const   {  return value[0] ^ other.value[0];  }
    NativeInt32 operator&  (const NativeInt32& other)   const   {  return value[0] & other.value[0];  }
    NativeInt32 operator&  (const char& other)          const   {  return value[0] & other;        }
    NativeInt32 operator~  ()                           const   {  return ~value[0];               }
    NativeInt32 operator<< (const int& coeff)           const   {  return value[0] << coeff;       }
    NativeInt32 operator>> (const int& coeff)           const   {  return value[0] >> coeff;       }
    bool        operator!= (const NativeInt32& c)       const   {  return value[0] != c.value[0];     }
    bool        operator== (const NativeInt32& c)       const   {  return value[0] == c.value[0];     }
    bool        operator<  (const NativeInt32& c)       const   {  return value[0] < c.value[0];      }
    bool        operator<= (const NativeInt32& c)       const   {  return value[0] <= c.value[0];     }
    bool        operator>= (const NativeInt32& c)       const   {  return value[0] >= c.value[0];     }
    NativeInt32& operator+=  (const NativeInt32& other)    {  value[0] += other.value[0]; return *this; }
    NativeInt32& operator^=  (const NativeInt32& other)    {  value[0] ^= other.value[0]; return *this; }

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const NativeInt32 & l)
    {
        s << std::hex << l.value[0] << std::dec;  return s;
    }

    /********************************************************************************/
    inline static hid_t hdf5 (bool& isCompound)
    {
        return H5Tcopy (H5T_NATIVE_UINT32);
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_NATIVE_16_HPP_ */
