//
//  a_basic_test.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__a_basic_test__
#define __RNAMake__a_basic_test__

#include <stdio.h>

#include "eternabot/strategy.h"

namespace eternabot {

class ABasicTest : public Strategy {
public:
    ABasicTest() {
        params_ = Floats(7);
        params_[0] = 0.423966515526;
        params_[1] = 93.0846762218;
        params_[2] = -1.87306343871;
        params_[3] = 1.15267359738;
        params_[4] = 0.9458;
        params_[5] = 63.60;
        params_[6] = 102.0;
        mean_ = 83.5007560083;
        stdev_ = 10.5290224709;
        name_ = "ABasicTest";
    }

    ~ABasicTest() {}

    inline
    float
    score(FeaturesOP const & features) {
        float total_pairs = features->gc + features->gu + features->ua;
        float score = 100;
        if(total_pairs > 0) {
            score -= fabs(features->ua / total_pairs - params_[0]) * params_[1];
        }
        float target_fe = params_[2] * total_pairs;
        score -= fabs(target_fe - features->fe) * params_[3];

        if(features->meltpoint < params_[5]) {
            score -= fabs(features->meltpoint - params_[5]) * params_[4];
        }
        else if(features->meltpoint > params_[6]) {
            score -= fabs(features->meltpoint - params_[6]) * params_[4];
        }
        return score;
    }

};

}

#endif /* defined(__RNAMake__a_basic_test__) */