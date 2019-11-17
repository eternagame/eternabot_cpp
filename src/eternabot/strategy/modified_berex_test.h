//
//  berex_test.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake_mod_berex_test__
#define __RNAMake_mod_berex_test__

#include <stdio.h>

#include "eternabot/strategy.h"

namespace eternabot {

class ModifiedBerexTest: public Strategy {
public:
    ModifiedBerexTest() {
        params_ = std::vector<float>(12);
        params_[0] = 0.233331320365;
        params_[1] = 67.1928215814;
        params_[2] = 0.0846125555161;
        params_[3] = 97.8926617326;
        params_[4] = 0.124731583984;
        params_[5] = 76.3152913746;
        params_[6] = -64.8000139266;
        params_[7] = -25.1005947191;
        params_[8] = 1.2097925221;
        params_[9] = 33.3729652592;
        params_[10] = 171.705337159;
        params_[11] = 1.28554296532;
        mean_ = 84.0125821249;
        stdev_ = 8.91633847502;
        name_ = "ModifiedBerexTest";
    }
    
    ~ModifiedBerexTest() {}
    
public:
    float
    score(FeaturesOP const & features) {
      
        float score = 100;

        if(features->length > 30) {
            score -= abs(float(features->g_count) / float(features->length) - params_[0]) * params_[1];
            score -= abs(float(features->u_count) / float(features->length) - params_[2]) * params_[3];
            score -= abs(float(features->c_count) / float(features->length) - params_[4]) * params_[5];
        }

        float weight = exp(-((features->length - 100)*(features->length - 100))/5000)/(sqrt(2*3.14*1))*2.5;
        if     (features->fe < params_[6]) {
           score -= abs(features->fe - params_[6]) * params_[8] * weight;
        }
        else if(features->fe > params_[7]) {
           score -= abs(features->fe - params_[7]) * params_[8] * weight;
        }

        if     (features->meltpoint < params_[9]) {
            score -= abs(features->meltpoint - params_[9]) * params_[11] * weight;
        }
        else if(features->meltpoint > params_[10]) {
            score -= abs(features->meltpoint - params_[10]) * params_[11] * weight;
        }

        /*if     (features->fe < params_[6]) {
            score -= abs(features->fe - params_[6]) * params_[8];
        }
        else if(features->fe > params_[7]) {
            score -= abs(features->fe - params_[7]) * params_[8];
        }*/
        
        /*if     (features->meltpoint < params_[9]) {
            score -= abs(features->meltpoint - params_[9]) * params_[11];
        }
        else if(features->meltpoint > params_[10]) {
            score -= abs(features->meltpoint - params_[10]) * params_[11];
        }*/
        
        return score;
    }

};

}


#endif /* defined(__RNAMake__berex_test__) */
