//
//  strategy.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__strategy__
#define __RNAMake__strategy__

#include <stdio.h>

#include "eternabot/feature_generator.h"

namespace eternabot {
    
class Strategy {
public:
    
    Strategy() {}
    
    virtual
    ~Strategy() {}
    
public:
    
    virtual
    float
    score(FeaturesOP const &) = 0;
    
public:
    inline
    float
    mean() const { return mean_; }
    
    inline
    float
    stdev() const { return stdev_; }

    inline
    String const &
    name() { return name_; };

protected:
    float mean_, stdev_;
    Floats params_;
    String name_;
    
};

typedef std::shared_ptr<Strategy> StrategyOP;
typedef std::vector<StrategyOP>   StrategyOPs;
    
}

#endif /* defined(__RNAMake__strategy__) */
