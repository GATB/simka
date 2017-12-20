#include <algorithm>
#include <iostream>
#include <cstring>
#include <cmath>
#include <memory>
#include <stdint.h>
#include <type_traits>
#include <utility>
#include <cassert>
#include "bias_corrected_estimate.hpp"
#include "linear_counting.hpp"
#include "../MurmurHash3.h"


#ifndef _HLL_RAW_H_
#define _HLL_RAW_H_

enum class Format {NORMAL, COMPACT_6BITS, COMPACT_5BITS, COMPACT_4BITS};

/**
 * T is the hashable type. The class can be used to count various types.
 * H is a class deriving from Hash<T>
 */
template<typename T>//, typename H = MurMurHash<T> >
class HllRaw {
/**
 * The line below is the fanciest thing I've ever done.
 * This enforces in compile time that the H parameter is a subclass of Hash<T>.
 *
 * Weirdly, static_assert is not a member of std::
 */
  //static_assert(std::is_base_of<Hash<T>, H>::value,
   // "Hll's H parameter has to be a subclass of Hash<T>");
private:
  /**
   * The below member variables' order aims at minimizing the memory fooprint
   * Please don't change it if it's not 100% necessary.
   */
  uint64_t numberOfBuckets;
  uint64_t bucketSize;
  uint64_t bucketMask;
  uint64_t valueMask;

  uint8_t bucketBits;
  uint8_t valueBits;
  uint8_t* synopsis;

  uint8_t leftMostSetBit(uint64_t hash) const {
    // We set the bucket bits to 0
    if (hash == 0)
      return 0;
    else
      /**
       * clz returns number of leading zero bits couting from the MSB
       * we have to add 1 to count the set bit
       * then we subtract bucketBits, since they are zeroed by valueMask
       */
      return (uint8_t)__builtin_clzll(hash & valueMask) + 1 - bucketBits;
  }



  void init(uint8_t bucketBits) {
    this -> bucketBits = bucketBits;
    this -> valueBits = 64 - bucketBits;
    this -> numberOfBuckets = 1UL << bucketBits;


    // Ex: with a 4 bits bucket on a 2 bytes size_t we want 1111 0000 0000 0000
    // So that's 10000 minus 1 shifted with 12 zeroes 1111 0000 0000 0000
    this -> bucketMask = (( 1UL << bucketBits ) - 1UL) << valueBits;

    // Ex: with a 4 bits bucket on a 16 bit size_t  we want 0000 1111 1111 1111
    // So that's 1 with 12 ( 16 - 4 ) zeroes 0001 0000 0000 0000, minus one = 0000 1111 1111 1111
    this -> valueMask = (1UL << valueBits ) - 1UL;

    /**
     * No matter what the serialization method is, we store the buckets as an array of bytes.
     */
    this -> synopsis = new uint8_t[numberOfBuckets];
  }
public:



  uint64_t bucket(uint64_t hash) {
    // Get the most significant bits and shift that
    // For example on a 16-bit platform, with a 4 bit bucketMask (valueBits = 12)
    // The value 0110 0111 0101 0001 is first bit masked to 0110 0000 0000 0000
    // And right shifted by 12 bits, which gives 0110, bucket 6.
    return (hash & bucketMask) >> valueBits;
  }

  /**
   * Below we implement the C++11's rule of five.
   * Since this class uses manual memory allocation, we have to to define:
   * - copy constructor,
   * - move constructor,
   * - copy assignment,
   * - move assignmnet,
   * - destructor.
   *
   * See: https://en.wikipedia.org/wiki/Rule_of_three_(C%2B%2B_programming)
   */
  HllRaw(uint8_t bucketBits) { 
    init(bucketBits);
    memset( synopsis, 0, numberOfBuckets );
  }


  /** Copy constructor */
  HllRaw(const HllRaw& other) : numberOfBuckets(other.numberOfBuckets),
    bucketMask(other.bucketMask),
    valueMask(other.valueMask),
    bucketBits(other.bucketBits),
    valueBits(other.valueBits),
    synopsis(new uint8_t[other.numberOfBuckets]) {
    memcpy(synopsis, other.synopsis, numberOfBuckets);
  }

  /** Move constructor */
  HllRaw(HllRaw&& other) noexcept : numberOfBuckets(other.numberOfBuckets),
    bucketMask(other.bucketMask),
    valueMask(other.valueMask),
    bucketBits(other.bucketBits),
    valueBits(other.valueBits) {
      // we swap the values, so that we get the synopsis from the temporary
      // object, and the current synopsis will be destroyed together with
      // the temporary object.
      std::swap(synopsis, other.synopsis);
    }

  /** Copy assignment operator */
  HllRaw& operator=(const HllRaw& other) {
    HllRaw tmp(other);
    *this = std::move(tmp);
    return *this;
  }

  /** Move assignment operator */
  HllRaw& operator=(HllRaw&& other) noexcept {
    delete[] synopsis;
    synopsis = other.synopsis;
    other.synopsis = nullptr;
    return *this;
  }

  ~HllRaw() {
    delete[] synopsis;
  }

  void add(T hashValue) {
    //H hashFunction;
    //uint64_t hashValue = hashFunction(value);

    const uint32_t dstBucket = bucket(hashValue);
    // We store in the synopsis for the the biggest leftmost one
    synopsis[dstBucket] = std::max(synopsis[dstBucket], leftMostSetBit(hashValue));
  }

  /*
  void add(const uint8_t __restrict__ otherSynopsis[]) {
    const uint32_t numberOfBucketsConst = numberOfBuckets;
    uint8_t* __restrict__ synopsis_ = synopsis;

    // note: LOOP VECTORIZED
    for (uint64_t i = 0; i < numberOfBucketsConst; i++) {
      synopsis_[i] = std::max(synopsis_[i], otherSynopsis[i]);
    }
  }

  void add(const HllRaw<T>& other) {
    assert(this->numberOfBuckets == other.numberOfBuckets);
    this->add(other.synopsis); 
  }
*/
  uint8_t getBucketBits() const {
    return bucketBits;
  }

  uint8_t* getCurrentSynopsis() {
    return this -> synopsis;
  }

  uint64_t getNumberOfBuckets() const {
    return this -> numberOfBuckets;
  }

  uint32_t getSynopsisSize(Format format) {
    uint32_t ret = 0;
    if(format == Format::NORMAL) {
      ret = numberOfBuckets;

    } else if(format == Format::COMPACT_6BITS) {
      uint8_t outputBucketSizeBits = static_cast<uint8_t>(std::log2(valueBits) + 0.5); //6
      uint32_t outputArraySizeBits = numberOfBuckets * outputBucketSizeBits;
      uint32_t outputArraySize = outputArraySizeBits >> 3; //12288
      ret = outputArraySize;

    } else if(format == Format::COMPACT_5BITS) {
      // we need 5 bits for every bucket
      // also, these buckets are stored in bytes, so we divide by 8
      ret = (5 * numberOfBuckets)/8;
    } else if(format == Format::COMPACT_4BITS) {
      ret = (4 * numberOfBuckets) / 8;
    } else {
      assert(0);
    }
    return ret;
  }

  uint32_t emptyBucketsCount() const {
    uint32_t emptyBuckets = 0;
    // note: LOOP VECTORIZED
    for (uint64_t i = 0; i < numberOfBuckets; i++) {
      emptyBuckets += static_cast<uint32_t>(synopsis[i] == 0);
    }
    return emptyBuckets;
  }

  /**
   * This function is a direct implementation of the formula presented in
   * Flajolet's paper. According to it the cardinality estimation is calculated
   * as follows:
   *
   * E = \alpha_m * m^2 * sum_{j=1}^m{2^{-M[j]}}
   *
   * where:
   *    E is the estimation
   *    \alpha_m is a constant derived from another formula
   *    m is number of buckets
   *    M[j] is count from the j-th bucket
   *
   */
  uint64_t estimate() const {
    double harmonicMean = 0.0;
    double alpha;
    // Alpha computation (see paper)
    switch (bucketBits) {
      case 4:
        alpha = 0.673;
        break;
      case 5:
        alpha = 0.697;
        break;
      case 6:
        alpha = 0.709;
      // since 11 and 14 are most likely cases, we make it computable in compile time
      case 11:
        alpha = (0.7213 / (1.0 + (1.079 / (1<<11)))); //0.72092
      case 14:
        alpha = (0.7213 / (1.0 + (1.079 / (1<<14)))); //0.721253
      default:
        alpha = (0.7213 / (1.0 + (1.079 / static_cast<double>(numberOfBuckets))));
    }
    for (uint64_t i = 0; i < numberOfBuckets; i++)
    {
        harmonicMean += 1.0 / (1 << (synopsis[i]));
    }
    harmonicMean = numberOfBuckets / harmonicMean;
    // std::llround returns a long long
    // note: other rounding functions return a floating point or shorter types
    uint64_t hllEstimate = std::llround(0.5 + alpha * harmonicMean * numberOfBuckets);
    return hllEstimate;
  }

  void deserialize8Bits(const char* byteArray1) {
    const unsigned char* __restrict__ byteArray = reinterpret_cast<const unsigned char*>(byteArray1);
    uint8_t* __restrict__ synopsis_ = this->synopsis;
    const uint32_t numberOfBucketsConst = getNumberOfBuckets();

    // note: LOOP VECTORIZED
    for(uint32_t i=0; i<numberOfBucketsConst; ++i) {
      synopsis_[i] = byteArray[i];
    }
  }

  void serialize8Bits(char* byteArray1) const {
    unsigned char* byteArray = reinterpret_cast<unsigned char*>(byteArray1);
    const uint32_t numberOfBucketsConst = getNumberOfBuckets();
    uint8_t* __restrict__ synopsis_ = this->synopsis;

    // note: LOOP VECTORIZED
    for(uint32_t i=0; i<numberOfBucketsConst; ++i) {
      byteArray[i] = synopsis_[i];
    }
  }
/**
 * We use dense representation for storing the buckets.
 * In the Hll class we store the synopsis as an array of 2^14 bytes. In fact,
 * since the counters use only 6 bits, when storing the synopsis on a permanent
 * storage we can compress it to 2^14 * 6 bits =  12288 bytes.
 * The idea is to split the buckets into groups of 4 and to store them in 3 bytes
 * (since 4 * 6 bits = 24 bits = 3 bytes).
 *
 * Hence, we use the following representation:
 *
 * +--------+--------+--------+---//
 * |00000011|11112222|22333333|4444
 * +--------+--------+--------+---//
 *
 * In order to deserialize an array of bytes into buckets, we iterate over
 * groups of three bytes which get stored in four buckets.
 * Likewise, to serialize the buckets, we iterate over groups of 4 buckets
 * which get stored in three bytes.
 *
 * In the functions below it's essential to set value of a single byte in only
 * one assignment. Otherwise we risk encountering a write after write hazard,
 * which could make the operation significantly slower.
 */

  void deserialize6Bits(const char* byteArray1) {
    //bgidx stands for bucket group index
    //
    const unsigned char* __restrict__ byteArray = reinterpret_cast<const unsigned char*>(byteArray1);
    uint8_t* __restrict__ synopsis_ = this->synopsis;
    const uint32_t numberOfBucketsConst = getNumberOfBuckets();

    //note: LOOP VECTORIZED
    for(uint32_t bgidx = 0; bgidx < numberOfBucketsConst/4; ++bgidx) {
      synopsis_[bgidx*4] = byteArray[bgidx*3] >> 2;
      synopsis_[bgidx*4+1] = ((byteArray[bgidx*3] & 0x3) << 4) | (byteArray[bgidx*3+1] >> 4);
      synopsis_[bgidx*4+2] = ((byteArray[bgidx*3+1] & 0xF) << 2) | (byteArray[bgidx*3+2] >> 6);
      synopsis_[bgidx*4+3] = (byteArray[bgidx*3+2] & 0x3F);
    }
  }

  void serialize6Bits(char* byteArray1) const {
    //bgidx stands for bucket group index
    unsigned char* __restrict__ byteArray = reinterpret_cast<unsigned char*>(byteArray1);
    uint8_t* __restrict__ synopsis_ = this->synopsis;
    const uint32_t numberOfBucketsConst = getNumberOfBuckets();

    //note: LOOP VECTORIZED
    for(uint32_t bgidx = 0; bgidx < numberOfBucketsConst/4; ++bgidx) {
      byteArray[bgidx*3]   = (synopsis_[bgidx*4] << 2)    | (synopsis_[bgidx*4+1] >> 4);
      byteArray[bgidx*3+1] = (synopsis_[bgidx*4+1] << 4)  | (synopsis_[bgidx*4+2] >> 2);
      byteArray[bgidx*3+2] = (synopsis_[bgidx*4+2] << 6 ) | synopsis_[bgidx*4+3];
    }
  }

/**
 * Following functions serialize and deserialize the buckets using
 * 5 bits per bucket. To fit the buckets nicely in the array of bytes,
 * we have to split the buckets into groups of 8 and put them in 5
 * bytes as follows:
 *
 *   Byte 0   Byte 1   Byte 2   Byte 3   Byte 4
 * +--------+--------+--------+--------+--------+---//
 * +00000111|11222223|33334444|45555566|66677777|...
 * +--------+--------+--------+--------+--------+---//
 */
  void deserialize5BitsWithBase(const char* byteArray1, uint8_t base) {
    // a cast to be OK with types
    const unsigned char* __restrict__ byteArray = reinterpret_cast<const unsigned char*>(byteArray1);
    uint8_t* __restrict__ synopsis_ = this->synopsis;
    // we iterate over groups of 5 bytes and put them in 8 buckets
    const uint32_t numberOfBucketsConst = getNumberOfBuckets();
    for(uint32_t bgidx = 0; bgidx < numberOfBucketsConst/8; ++bgidx) {
      synopsis_[bgidx*8]   = base +   (byteArray[bgidx*5] >> 3);
      synopsis_[bgidx*8+1] = base + (((byteArray[bgidx*5]   & 0x07) << 2) | (byteArray[bgidx*5+1] >> 6));
      synopsis_[bgidx*8+2] = base +  ((byteArray[bgidx*5+1] & 0x3E) >> 1);
      synopsis_[bgidx*8+3] = base + (((byteArray[bgidx*5+1] & 0x01) << 4) | (byteArray[bgidx*5+2] >> 4));
      synopsis_[bgidx*8+4] = base + (((byteArray[bgidx*5+2] & 0x0F) << 1) | (byteArray[bgidx*5+3] >> 7));
      synopsis_[bgidx*8+5] = base +  ((byteArray[bgidx*5+3] & 0x7C) >> 2);
      synopsis_[bgidx*8+6] = base + (((byteArray[bgidx*5+3] & 0x03) << 3) | (byteArray[bgidx*5+4] >> 5));
      synopsis_[bgidx*8+7] = base +   (byteArray[bgidx*5+4] & 0x1F);
    }   
  }

  uint8_t serialize5BitsWithBase(char* byteArray1) const {
    uint8_t base = *std::min_element(synopsis, synopsis+numberOfBuckets);
    unsigned char*  byteArray = reinterpret_cast<unsigned char*>(byteArray1);

    // we iterate over groups of 8 buckets
    const uint32_t numberOfBucketsConst = getNumberOfBuckets();
    uint8_t  buckets[8];
    uint8_t*  synopsis_ = this->synopsis;
    for(uint32_t bgidx = 0; bgidx < numberOfBucketsConst/8; ++bgidx) {

      // normalize the buckets, i.e. subtract the base (min. value) and
      // make sure that the remainder fits into 5 bits (and cut off if not)
      for(uint32_t bidx = 0; bidx < 8; ++bidx) {
        const uint32_t bgidx_ = bgidx;
        uint8_t normBucket = synopsis_[bgidx_*8 + bidx]-base; // normalized bucket
        const uint8_t maxValIn5Bits = ((1<<5)-1); // max value fitting 5 bits
        buckets[bidx] = normBucket > maxValIn5Bits ? maxValIn5Bits : normBucket;
      }
      byteArray[bgidx*5]   = (buckets[0] << 3) | (buckets[1] >> 2);
      byteArray[bgidx*5+1] = (buckets[1] << 6) | (buckets[2] << 1) | (buckets[3] >> 4); 
      byteArray[bgidx*5+2] = (buckets[3] << 4) | (buckets[4] >> 1);
      byteArray[bgidx*5+3] = (buckets[4] << 7) | (buckets[5] << 2) | (buckets[6] >> 3);
      byteArray[bgidx*5+4] = (buckets[6] << 5) | (buckets[7]);
    }
    return base;
  }

/*
 * +--------+--------+---//
 * +00001111|22223333|...
 * +--------+--------+---//
 */
  void deserialize4BitsWithBase(const char* byteArray1, uint8_t base) {
    // a cast to be OK with the types
    const unsigned char* __restrict__ byteArray = reinterpret_cast<const unsigned char*>(byteArray1);

    uint8_t* __restrict__ synopsis_ = synopsis;
    const uint32_t numberOfBucketsConst = getNumberOfBuckets();

    // note: LOOP VECTORIZED
    for(uint32_t bgidx = 0; bgidx < numberOfBucketsConst/2; ++bgidx) {
      synopsis_[bgidx*2]   = base + (byteArray[bgidx] >> 4);
      synopsis_[bgidx*2+1] = base + (byteArray[bgidx] & 0x0f);
    }
  }

  uint8_t serialize4BitsWithBase(char* byteArray1) const {
    uint8_t base = *std::min_element(synopsis, synopsis+numberOfBuckets);
    unsigned char* __restrict__ byteArray = reinterpret_cast<unsigned char*>(byteArray1);

    const uint8_t maxValIn4Bits = ((1<<4)-1); // max value fitting 4 bits
    // we iterate over pairs of buckets
    
    uint8_t* __restrict__ synopsis_ = synopsis;
    const uint32_t numberOfBucketsConst = getNumberOfBuckets();

    // note: LOOP VECTORIZED
    for(uint32_t bgidx = 0; bgidx < numberOfBucketsConst/2; ++bgidx) {
      uint8_t normBucket1 = synopsis_[2*bgidx] - base;
      normBucket1 = (normBucket1 > maxValIn4Bits ? maxValIn4Bits : normBucket1);

      uint8_t normBucket2 = synopsis_[2*bgidx+1] - base;
      normBucket2 = (normBucket2 > maxValIn4Bits ? maxValIn4Bits : normBucket2);

      byteArray[bgidx] = (normBucket1 << 4) | normBucket2;
    }
    return base;
  }
};

#endif
