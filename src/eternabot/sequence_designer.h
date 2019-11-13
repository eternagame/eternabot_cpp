//
//  sequence_designer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_designer__
#define __RNAMake__sequence_designer__

#include <stdio.h>

#include <base/option.h>
#include <util/random_number_generator.h>
#include <util/monte_carlo.h>
#include <secondary_structure/sequence_constraint.h>
#include <eternabot/scorer.h>

namespace eternabot {

struct SequenceDesignerResult {
    inline
    SequenceDesignerResult(
            String const & n_sequence,
            float n_score,
            float n_bp_diff_score):
            sequence(n_sequence),
            score(n_score),
            bp_diff_score(n_bp_diff_score){}

    String sequence;
    float score;
    float bp_diff_score;
};
    
typedef std::shared_ptr<SequenceDesignerResult> SequenceDesignerResultOP;
typedef std::vector<SequenceDesignerResultOP> SequenceDesignerResultOPs;

struct sequence_designer_result_less_than_key {
    inline bool operator() (
        SequenceDesignerResultOP const & r1,
        SequenceDesignerResultOP const & r2) {
    
        return (r1->score > r2->score);
    }
};

class MonteCarloMove {
public:
    MonteCarloMove() = default;

    virtual
    ~MonteCarloMove() = default;

    virtual
    MonteCarloMove *
    clone() const = 0;

public:
    virtual
    int
    move(
            secondary_structure::PoseOP) = 0;

    virtual
    void
    undo(
            secondary_structure::PoseOP) = 0;
};

class MutateBPMove : public MonteCarloMove {
public:
    MutateBPMove(
            secondary_structure::BasepairOPs designable_bps,
            std::vector<secondary_structure::ResTypes> possible_rt_types):
            designable_bps_(designable_bps),
            possible_rt_types_(possible_rt_types),
            rng_(util::RandomNumberGenerator()) {
        last_res_types_ = secondary_structure::ResTypes(2);
    }

    ~MutateBPMove() = default;

    MonteCarloMove *
    clone() const { return new MutateBPMove(*this); }

public:
    int
    move(
            secondary_structure::PoseOP p) {

        if(designable_bps_.size() == 0) { return 0; }

        current_ = designable_bps_[rng_.randrange((int)designable_bps_.size())];
        last_res_types_[0] = current_->res1()->res_type();
        last_res_types_[1] = current_->res2()->res_type();

        while(true) {
            current_res_types_ = possible_rt_types_[rng_.randrange((int) possible_rt_types_.size())];
            current_->res1()->res_type(current_res_types_[0]);
            current_->res2()->res_type(current_res_types_[1]);
            if(current_res_types_[0] != last_res_types_[0]) { break; }

        }

        return 1;
    }

    void
    undo(
            secondary_structure::PoseOP p) {
        current_->res1()->res_type(last_res_types_[0]);
        current_->res2()->res_type(last_res_types_[1]);
    }

private:
    secondary_structure::BasepairOPs designable_bps_;
    secondary_structure::BasepairOP current_;
    std::vector<secondary_structure::ResTypes> possible_rt_types_;
    util::RandomNumberGenerator rng_;
    secondary_structure::ResTypes last_res_types_, current_res_types_;

};

class MutateUnpairedResMove : public MonteCarloMove {
public:
    MutateUnpairedResMove(
            secondary_structure::ResidueOPs designable_unpaired_res):
            designable_unpaired_res_(designable_unpaired_res) {
        possible_res_types_ = secondary_structure::ResTypes {
            secondary_structure::ResType::ADE,
            secondary_structure::ResType::CYT,
            secondary_structure::ResType::GUA,
            secondary_structure::ResType::URA
        };

    }

    ~MutateUnpairedResMove() = default;

    MonteCarloMove *
    clone() const { return new MutateUnpairedResMove(*this); }

public:
    int
    move(
            secondary_structure::PoseOP p) {

        if(designable_unpaired_res_.size() == 0) { return 0; }

        current_ = designable_unpaired_res_[rng_.randrange((int)designable_unpaired_res_.size())];
        last_res_type_ = current_->res_type();
        while(true) {
            current_res_type_ = possible_res_types_[rng_.randrange((int)possible_res_types_.size())];
            current_->res_type(current_res_type_);
            if(current_res_type_ != last_res_type_) { break; }
        }
        return 1;
    }

    void
    undo(
            secondary_structure::PoseOP p) {
        current_->res_type(last_res_type_);

    }

private:
    secondary_structure::ResidueOP current_;
    secondary_structure::ResidueOPs designable_unpaired_res_;
    secondary_structure::ResTypes possible_res_types_;
    secondary_structure::ResType last_res_type_, current_res_type_;
    util::RandomNumberGenerator rng_;

};


typedef std::shared_ptr<MonteCarloMove> MonteCarloMoveOP;
typedef std::vector<MonteCarloMoveOP> MonteCarloMoveOPs;

    
class SequenceDesigner {
public:
    SequenceDesigner();
    
    ~SequenceDesigner() = default;

public:
    
    void
    setup();
    
    SequenceDesignerResultOPs const &
    design(secondary_structure::PoseOP const &);
    
    
public: //option wrappers

    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }
    
    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }
    
    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }
    
    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }
    
    inline
    bool
    has_option(String const & name) { return options_.has_option(name); }
    
    template<typename T>
    void
    set_option_value(
        String const & name,
        T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }
    
protected:
    
    void
    setup_options();
    
    void
    update_var_options();

private:
    void
    _find_designable_bps(
            secondary_structure::PoseOP);

    void
    _generate_inital_sequence(
            secondary_structure::PoseOP);

    bool
    _new_sequence_violations(
            Ints const &,
            Ints const &);


private: // new and possibly badly made functions
    void
    _set_initial_helix_sequence(
            secondary_structure::MotifOP);

    void
    _get_random_res_type_pair_gc_cap(
            secondary_structure::ResTypes &);

    void
    _get_random_res_type_pair(
            secondary_structure::ResTypes &);

    float
    _optimize_substructure(
            secondary_structure::PoseOP,
            int);

    float
    _bp_list_diff(
            secondary_structure::PoseOP,
            std::vector<std::vector<int>> const &,
            size_t,
            FeaturesOP);

private:
    struct Parameters {
        bool biased_gc_caps;
    };


private:
    secondary_structure::ResidueOPs designable_res_;
    secondary_structure::ResidueOPs designable_unpaired_res_;
    std::vector<Strings> possible_bps_;
    std::vector<secondary_structure::ResTypes> possible_rt_bps_;
    base::Options options_;
    util::RandomNumberGenerator rng_;
    util::MonteCarlo mc_;
    secondary_structure::BasepairOPs designable_bps_;
    Scorer scorer_;
    SequenceDesignerResultOPs results_;

    // tracking sequence constraints
    Ints current_violations_, next_violations_;
    secondary_structure::SequenceConstraints seq_constraints_;

    // current solutions
    std::map<secondary_structure::MotifOP, secondary_structure::ResTypes> current_restypes_;

    // directly keep track of vienna pair map
    std::vector<std::vector<int>> pair_map_;
    size_t pair_map_entries_;
    vienna::Vienna v_;


    int designs_, steps_;
    float temperature_;

    Parameters parameters_;
    
    
};
    
}

#endif /* defined(__RNAMake__sequence_designer__) */
