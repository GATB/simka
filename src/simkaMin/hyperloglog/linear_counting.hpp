#ifndef _LINEAR_COUNTING_H_
#define _LINEAR_COUNTING_H_


/**
 * This algorithm was proposed in Flajolet's paper as a remedy to a tremendous
 * error one can observe when the cardinality being estimated is low.
 *
 * In LinearCounting we use a bimap of 2^LINEAR_COUNTING_BITS bits. Whenever
 * a new value is added, we take the first LINEAR_COUNTING_BITS and use them
 * to set a bit in the bitmap.
 *
 * Eventually, the ratio of set to all bits is used to estimate cardinality
 * of the measured set. 
 */
class LinearCounting {
private:
  uint8_t precision;
  uint8_t* bitmap;

  static const uint64_t linearCountingThreshold[];
  uint32_t getBits(uint64_t hashValue) const;
  uint32_t countUnsetBits() const;

  void setBit(uint32_t bitIdx);
  uint32_t bitmapSize() const;

public:
  LinearCounting(uint8_t precision);
  virtual ~LinearCounting();
  LinearCounting(const LinearCounting& other);
  LinearCounting(LinearCounting&& other) noexcept;
  LinearCounting& operator=(const LinearCounting& other);
  LinearCounting& operator=(LinearCounting&& other);

  void add(uint64_t hashValue);
  static uint64_t getLinearCountingThreshold(uint8_t hllPrecision);
  uint64_t estimate() const;

};

#endif
