#ifndef _BIAS_CORRECTED_ESTIMATE_H_
#define _BIAS_CORRECTED_ESTIMATE_H_

#include <vector>

class BiasCorrectedEstimate {
  static double kNeighborsInterpolationBias(uint64_t rawEstimate, uint8_t precision);

  static const std::vector<std::vector<double>> rawEstimateData;
  static const std::vector<std::vector<double>> biasData;

public:
  static uint64_t estimate(uint64_t rawEstimate, uint8_t precision);
};


#endif