#include <algorithm>
#include <cassert>
#include <cstring>
#include <iterator>
#include <iostream>
#include <utility>
#include "linear_counting.hpp"
#include <cmath>


LinearCounting::LinearCounting(uint8_t prec) : precision(prec) {
  // LC bitmap has to be at least 8 bytes long. To this end,
  // precision ought to be greater or equal to 6.
  assert(prec >= 6);
  bitmap = new uint8_t[bitmapSize()];
  memset(bitmap, 0, bitmapSize());
}

LinearCounting::~LinearCounting() { delete[] bitmap; }

LinearCounting::LinearCounting(const LinearCounting& other) : precision(other.precision) {
  bitmap = new uint8_t[bitmapSize()];
  memcpy(bitmap, other.bitmap, 1 << (precision-3));
}

LinearCounting::LinearCounting(LinearCounting&& other) noexcept {
  std::swap(bitmap, other.bitmap);
}

LinearCounting& LinearCounting::operator=(const LinearCounting& other) {
  LinearCounting tmp(other); //reuse copy constructor
  *this = std::move(tmp);
  return *this;
}

LinearCounting& LinearCounting::operator=(LinearCounting&& other) {
  delete[] bitmap;
  bitmap = other.bitmap;
  other.bitmap = nullptr;
  return *this;
}


uint32_t LinearCounting::bitmapSize() const {
  return 1 << (precision-3);
}

uint32_t LinearCounting::getBits(uint64_t hashValue) const {
  return static_cast<uint32_t>(hashValue >> (64 - precision)); 
} 

/**
 * Count unset bits in the bitmap. To make it faster we convert the array
 * of bytes to an array of 64 bit integers and use __builtin_popcount.
 */
uint32_t LinearCounting::countUnsetBits() const {
  uint32_t count = 0;
  const uint64_t* begin = reinterpret_cast<const uint64_t*>(bitmap);
  const uint64_t* end = reinterpret_cast<const uint64_t*>(bitmap+bitmapSize());

  for (const uint64_t* ptr = begin; ptr != end; ++ptr) {
    // count number of unset bits in in the array
    // __builtin_popcountll counts number of set bits so we have to invert the table bitwise
    // NOTE: there is no set-bit counterpart of popcount
    count += __builtin_popcountll(~(*ptr));
  } 
  return count;
}

void LinearCounting::setBit(uint32_t bitIdx) {
  uint32_t byte_idx = bitIdx / 8;
  uint8_t bit_in_byte = bitIdx % 8;
  bitmap[byte_idx] |= (1 << bit_in_byte);
}

/**
* NOTE: This function expects a hash value and not a real value.
*/
void LinearCounting::add(uint64_t hashValue) {
  uint32_t relevantBits = static_cast<uint32_t>(hashValue >> (64 - precision)); 
  setBit(relevantBits);
}

uint64_t LinearCounting::getLinearCountingThreshold(uint8_t hllPrecision) {
  assert(hllPrecision >= 4 && hllPrecision <= 18);
  // we subtract 4 from precision since it's the minimal value
  return linearCountingThreshold[hllPrecision-4];
}

/**
 * LinearCounting gives its estimate according to the following formula:
 * E = m * ln(m/S)
 *
 * where:
 *   m is size of the bitmap
 *   S is number of unset bits (i.e. whose value =0)
 *   
 */
uint64_t LinearCounting::estimate() const {
  uint32_t unsetBits = countUnsetBits();
  uint32_t bitsInBitmap = 1 << precision; 
  return bitsInBitmap * log(static_cast<double>(bitsInBitmap) / static_cast<double>(unsetBits));
}

/** 
 * These values come from the Google paper appendix. Based on this array one can tell when
 * LinearCounting ought to be applied. For instance, if precision is 10 (i.e. 10 bits are used 
 * to choose the bucket), then array[10-4] = array[6] = 900. This means, that if the Bias-corrected
 * cardinality is below 900, LinearCounting's estimate should be used instead of the HLL's one.
 *
 * See: https://docs.google.com/document/d/1gyjfMHy43U9OWBXxfaeG-3MjGzejW1dlpyMwEYAAWEI/
 */

const uint64_t LinearCounting::linearCountingThreshold[] =
  {10, 20, 40, 80, 220, 400, 900, 1800, 3100, 6500, 11500, 20000, 50000, 120000, 350000};

