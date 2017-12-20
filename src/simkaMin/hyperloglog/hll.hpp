#include "hll_raw.hpp"
#include <exception>

#ifndef _HLL_H_
#define _HLL_H_



struct HLLHdr {
  char magic[2] = {'H','L'};
  char format;
  uint8_t bucketBase;
  char padding[4] = {'\0','\0','\0','\0'}; // padding to reach 8 bytes in length
};


struct SerializationError : public virtual std::runtime_error {
  SerializationError(const char* message) : std::runtime_error(message) {}
};


template<typename T>//, typename H = MurMurHash<T> >
class Hll {

  HllRaw<T> hll;
  LinearCounting linearCounting;
  uint32_t lcThreshold;
  uint32_t biasCorrectedThreshold;

  static uint8_t formatToCode(Format format) {
    uint8_t ret;
    if(format == Format::NORMAL) ret = 0x01;
    else if(format == Format::COMPACT_6BITS) ret = 0x02;
    else if(format == Format::COMPACT_5BITS) ret = 0x04;
    else if(format == Format::COMPACT_4BITS) ret = 0x08;
    else assert(0);
    return ret;
  }

public:
  Hll(uint8_t bucketBits, uint8_t lcBits) : 
      hll(bucketBits),
      linearCounting(lcBits) {
    // Google's paper suggests to set the threshold to this value 
   this -> biasCorrectedThreshold = hll.getNumberOfBuckets()*5;
   this -> lcThreshold = linearCounting.getLinearCountingThreshold(bucketBits);
  }

  Hll(uint8_t bucketBits) : Hll(bucketBits, bucketBits-4) {}

  void deserialize(const char* byteArray, Format format) {
    HLLHdr hdr = *(reinterpret_cast<const HLLHdr*>(byteArray));
    if(hdr.format != formatToCode(format))
      throw SerializationError("Requested serialization format is not the same as format"
        " in the synopsis' header.");

    const char* byteArrayHll = byteArray + sizeof(HLLHdr);
    if(format == Format::NORMAL) {
      hll.deserialize8Bits(byteArrayHll);

    } else if (format == Format::COMPACT_6BITS) {
      hll.deserialize6Bits(byteArrayHll);

    } else if (format == Format::COMPACT_5BITS) {
      uint8_t base = hdr.bucketBase;
      hll.deserialize5BitsWithBase(byteArrayHll, base);

    } else if (format == Format::COMPACT_4BITS) {
      uint8_t base = hdr.bucketBase;
      hll.deserialize4BitsWithBase(byteArrayHll, base);

    } else {
      throw SerializationError("Unknown format parameter in deserialize().");
    }
  }

  void serialize(char* byteArray, Format format) const {
    // for the time being we skip the header and serialize it once
    // the buckets are written down
    HLLHdr hdr;

    char* byteArrayHll = byteArray + sizeof(HLLHdr);
    uint8_t base = 0;
    if(format == Format::NORMAL) {
      hll.serialize8Bits(byteArrayHll);
    } else if (format == Format::COMPACT_6BITS) {
      hll.serialize6Bits(byteArrayHll);
    } else if (format == Format::COMPACT_5BITS) {
      base = hll.serialize5BitsWithBase(byteArrayHll);
    } else if (format == Format::COMPACT_4BITS) {
      base = hll.serialize4BitsWithBase(byteArrayHll);
    } else {
      throw SerializationError("Unknown format parameter in serialize().");
    }
    // serialize the header as well
    hdr.bucketBase = base;
    hdr.format = formatToCode(format);
    memcpy(byteArray, reinterpret_cast<char*>(&hdr), sizeof(HLLHdr));
  }

  void add(const Hll& other) {
    this->hll.add(other.hll); 
  }

  void add(T value) {
    hll.add(value);
    linearCounting.add(value);
  }

  uint32_t getSynopsisSize(Format format) {
    return hll.getSynopsisSize(format) + sizeof(HLLHdr);
  }

/**
 * Hll's error becomes significant for small cardinalities. For instance, when
 * the cardinality is 0, HLL(p=14) estimates it to ~11k.
 * To circumvent 
 */
  uint64_t approximateCountDistinct() {
    uint64_t e = this->hll.estimate();
    uint64_t ee;

    if(e <= biasCorrectedThreshold) {
      ee = BiasCorrectedEstimate::estimate(e, hll.getBucketBits());
    } else {
      ee = e;
    }

    uint64_t h;
    if(this->hll.emptyBucketsCount() != 0) {
      double v = static_cast<double>(hll.getNumberOfBuckets())/static_cast<double>(hll.emptyBucketsCount());
      h = hll.getNumberOfBuckets() * log(v);
    } else {
      h = ee;
    }

    if(h <= lcThreshold) {
      return h;
    } else {
      return ee;
    }
  }

};

#endif
