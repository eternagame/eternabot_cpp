//
//  sequence_designer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include <base/log.h>
#include <secondary_structure/sequence_tools.h>
#include <eternabot/sequence_designer.h>

namespace eternabot {


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// setup functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SequenceDesigner::SequenceDesigner():
        scorer_(Scorer()),
        results_(SequenceDesignerResultOPs()),
        rng_(util::RandomNumberGenerator()) {

    parameters_.biased_gc_caps = true;

    possible_bps_ = std::vector<Strings>({{"A", "U"}, {"U", "A"}, {"G", "C"}, {"C", "G"}});
    possible_rt_bps_ = std::vector<secondary_structure::ResTypes>{
            {secondary_structure::ResType::ADE, secondary_structure::ResType::URA},
            {secondary_structure::ResType::URA, secondary_structure::ResType::ADE},
            {secondary_structure::ResType::CYT, secondary_structure::ResType::GUA},
            {secondary_structure::ResType::GUA, secondary_structure::ResType::CYT}};


    // generate constraints
    auto disallowed_sequences = Strings{"AAAA", "CCCC", "GGGG", "UUUU"};
    for(auto const & seq : disallowed_sequences) { seq_constraints_.add_disallowed_sequence(seq); }
    seq_constraints_.add_gc_helix_stretch_limit(3);
    current_violations_ = Ints(seq_constraints_.num_constraints());
    next_violations_    = Ints(seq_constraints_.num_constraints());

    setup_options();
    temperature_ = 4.0f;
    mc_ = util::MonteCarlo(temperature_);
}


void
SequenceDesigner::setup() {}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// option functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
SequenceDesigner::setup_options() {
    options_.add_option("designs", 1, base::OptionType::INT);
    options_.add_option("steps", 1000, base::OptionType::INT);
    options_.lock_option_adding();
    update_var_options();
    
}
    
void
SequenceDesigner::update_var_options() {
    designs_ = options_.get_int("designs");
    steps_   = options_.get_int("steps");
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
SequenceDesignerResultOPs const &
SequenceDesigner::design(
        secondary_structure::PoseOP const & p) {
    results_ = SequenceDesignerResultOPs();

    auto motifs = p->helices();

    _find_designable_bps(p);       // find all basepairs with N-N residues
    for(auto const & r : p->residues()) {
        if(r->res_type() == secondary_structure::ResType::NONE) {
            designable_res_.push_back(r);
            auto bps = p->get_basepair(r->uuid());
            if(bps.size() == 0) {
                designable_unpaired_res_.push_back(r);
            }

            r->res_type(secondary_structure::ResType::ADE);
        }
    }

    for(auto const & h : motifs) {
        _set_initial_helix_sequence(h);
    }

    // get all non bp-step motifs
    for(auto const & m : p->motifs()) {
        if(m->mtype() == util::MotifType::HELIX) { continue; }
        motifs.push_back(m);
    }
    scorer_.setup(p);


    auto current_move = MonteCarloMoveOP(nullptr);
    auto moves = MonteCarloMoveOPs(2);
    moves[0] = std::make_shared<MutateBPMove>(designable_bps_, possible_rt_bps_);
    moves[1] = std::make_shared<MutateUnpairedResMove>(designable_unpaired_res_);

    std::cout << p->sequence() << std::endl;

    auto current_score = scorer_.score_secondary_structure(p);
    auto current_sequence = p->sequence();
    auto next_score = 0.0f;
    auto best_score = current_score;
    std::cout << current_score << std::endl;
    auto best_sequence = String("");
    auto pos = 0;
    exit(0);

    if(designable_unpaired_res_.size() == 0 && designable_bps_.size() == 0) {
        results_.push_back(std::make_shared<SequenceDesignerResult>(p->sequence(), current_score));
        return results_;
    }

    for(int i = 0; i < 10000; i++) {
        pos = rng_.randrange(1000);
        if(pos < 500) {
            current_move = moves[0];
        }
        else {
            current_move = moves[1];
        }

        // move didn't do anything, try again
        if(current_move->move(p) == 0) {
            i--;
            continue;
        }


        next_score = scorer_.score_secondary_structure(p);

        if(next_score > best_score) {
            best_score = next_score;
            best_sequence = p->sequence();
        }
    }

    std::cout << best_sequence << " " << best_score << std::endl;

    exit(0);

    /*

    for(auto const & m : p->motifs()) {
        if(m->mtype() == util::MotifType::HELIX) { continue; }
        motifs.push_back(m);
    }
    auto levels = std::map<secondary_structure::MotifOP, int>();
    auto children = std::map<secondary_structure::MotifOP, secondary_structure::MotifOPs>();
    auto used_motifs = std::map<secondary_structure::MotifOP, int>();
    auto head_motifs = secondary_structure::MotifOPs();

    for(auto const & m : motifs) {
        if(m->mtype() == util::MotifType::HAIRPIN) {
            levels[m] = 0;
            used_motifs[m] = 1;
            head_motifs.push_back(m);
        }
    }


    while(motifs.size() != used_motifs.size()) {
        for(auto const & m : motifs) {
            if(used_motifs.find(m) != used_motifs.end()) { continue; }
            for(auto const kv : used_motifs) {
                for(auto const & end : m->ends()) {
                    if (kv.first->ends()[0] == end) {
                        levels[m] = levels[kv.first] + 1;
                        children[kv.first].push_back(m);
                        used_motifs[m] = 1;
                    }
                }
            }
        }
     }

    auto sub_structures = std::vector<secondary_structure::MotifOPs>();

    for(auto const & m : head_motifs) {
        auto seen = std::map<secondary_structure::MotifOP, int>();
        auto current_children = children[m];
        std::cout << current_children.size() << std::endl;
    }

    exit(0);

    _generate_inital_sequence(p);  // generate initial sequence, get rid as many violations as possible
    scorer_.setup(p);

    auto current_score = scorer_.score_secondary_structure(p);
    auto current_sequence = p->sequence();
    auto best_score = current_score;
    auto best_sequence = p->sequence();
    auto next_score = 0.0f;
    auto last_pair = Strings{"", ""};

    for(int i = 0; i < steps_; i++) {
        auto pos = rng_.randrange(designable_bps_.size());
        auto & pair = _get_random_pair();
        last_pair[0] = designable_bps_[pos]->res1()->name();
        last_pair[1] = designable_bps_[pos]->res2()->name();
        _set_bp_sequence(pair, designable_bps_[pos]);
        next_violations_ = seq_constraints_.violations(p);

        if(_new_sequence_violations(current_violations_, next_violations_)) {
            _set_bp_sequence(last_pair, designable_bps_[pos]);
            continue;
        }

        next_score = scorer_.score_secondary_structure(p);
        if(mc_.accept(current_score, next_score)) {
            current_score = next_score;
        }
        else {
            _set_bp_sequence(last_pair, designable_bps_[pos]);
            continue;
        }

        if(current_score > best_score) {
            best_score = current_score;
            best_sequence = p->sequence();
        }
    }

    */

    //results_.push_back(std::make_shared<SequenceDesignerResult>(best_sequence, best_score));
    return results_;
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Strings const &
SequenceDesigner::_get_random_pair() {
    return possible_bps_[rng_.randrange(possible_bps_.size())];
}

void
SequenceDesigner::_set_bp_sequence(
        Strings const & pair,
        secondary_structure::BasepairOP bp) {
    bp->res1()->name(pair[0]);
    bp->res2()->name(pair[1]);
}


void
SequenceDesigner::_find_designable_bps(
        secondary_structure::PoseOP p) {
    for(auto const & bp : p->basepairs()) {
        if(bp->res1()->res_type() == secondary_structure::ResType::NONE &&
           bp->res2()->res_type() == secondary_structure::ResType::NONE ) {
            designable_bps_.push_back(bp);
        }
    }

}


bool
SequenceDesigner::_new_sequence_violations(
        Ints const & current_violations,
        Ints const & next_violations) {
    for(int i = 0; i < current_violations.size(); i++) {
        if(current_violations[i] < next_violations[i]) { return true; }
    }
    return false;
}

void
SequenceDesigner::_generate_inital_sequence(
        secondary_structure::PoseOP p) {

    auto count = 0;
    current_violations_ = seq_constraints_.violations(p);

    // inital sequence filling
    for (auto & bp : designable_bps_) {
        auto & p = _get_random_pair();
        _set_bp_sequence(p, bp);
    }

    next_violations_ = seq_constraints_.violations(p);
    while(_new_sequence_violations(current_violations_, next_violations_)) {
        count += 1;
        auto pos = rng_.randrange(designable_bps_.size());
        auto & pair = _get_random_pair();
        _set_bp_sequence(pair, designable_bps_[pos]);
        if(count > 1000000) {
            LOG_WARNING << "cannot find initial sequence that does not have sequence violations! ";
            current_violations_ = next_violations_;
            break;
        }
        next_violations_ = seq_constraints_.violations(p);
    }

}


void
SequenceDesigner::_set_initial_helix_sequence(
        secondary_structure::MotifOP h) {
    int n_bps = h->residues().size()/2;
    auto designable_1 = Ints(n_bps);
    auto designable_2 = Ints(n_bps);
    int i = 0;
    for(auto const & r : h->chains()[0]->residues()) {
        if(std::find(designable_res_.begin(), designable_res_.end(), r) != designable_res_.end()) {
            designable_1[i] = 1;
        }
        else {
            designable_1[i] = 0;
        }
        i++;
    }

    for(auto const & r : h->chains()[1]->residues()) {
        i--;
        if(std::find(designable_res_.begin(), designable_res_.end(), r) != designable_res_.end()) {
            designable_2[i] = 1;
        }
        else {
            designable_2[i] = 0;
        }
    }

    auto res_types_1 = secondary_structure::ResTypes(n_bps);
    auto res_types_2 = secondary_structure::ResTypes(n_bps);
    auto res_type_bp = secondary_structure::ResTypes(2);

    auto gc_count = 0;
    auto c_violations = seq_constraints_.violations(h);
    auto n_violations = Ints(c_violations);
    auto iter = 0;
    while(true) {
        i = 0;
        for(auto & bp : h->basepairs()) {
            // can design both sides of the helix
            if(designable_1[i] && designable_2[i]) {
                if ((i == 0 || i == n_bps - 1) && parameters_.biased_gc_caps) {
                    _get_random_res_type_pair_gc_cap(res_type_bp);
                } else {
                    _get_random_res_type_pair(res_type_bp);
                }
            }
            // cannot design either side
            else if(designable_1[i] == 0 && designable_2[i] == 0) {
                res_type_bp[0] = bp->res1()->res_type();
                res_type_bp[1] = bp->res2()->res_type();
            }
            // cannot design first residue
            else if(designable_1[i] == 0) {
                res_type_bp[0] = bp->res1()->res_type();
                res_type_bp[1] = secondary_structure::get_complement_res_type(bp->res1()->res_type());
            }
            // cannot design second residue
            else if(designable_2[i] == 0) {
                res_type_bp[0] = secondary_structure::get_complement_res_type(bp->res2()->res_type());
                res_type_bp[1] = bp->res2()->res_type();
            }

            res_types_1[i] = res_type_bp[0];
            res_types_2[i] = res_type_bp[1];

            bp->res1()->res_type(res_type_bp[0]);
            bp->res2()->res_type(res_type_bp[1]);

            i += 1;
        }

        iter += 1;
        if(iter > 1000) {
            LOG_WARNING << "could not find suitable starting sequence for helix";
            break;
        }

        n_violations = seq_constraints_.violations(h);
        if(_new_sequence_violations(c_violations, n_violations)) {
            continue;
        }

        break;


    }

    //std::cout << h->sequence() << std::endl;


}

void
SequenceDesigner::_get_random_res_type_pair_gc_cap(
        secondary_structure::ResTypes & pair) {
    // biased base pair selection for caped ends, 80% chance to be GC/CG over AU/UA
    // will do GCs
    if(rng_.randrange(1000) > 200) {
        // selecting GC
        if(rng_.randrange(1000) > 500) {
            pair[0] = secondary_structure::ResType::GUA;
            pair[1] = secondary_structure::ResType::CYT;
        }
        // selecting CG
        else {
            pair[0] = secondary_structure::ResType::CYT;
            pair[1] = secondary_structure::ResType::GUA;
        }
    }

    else {
        // selecting AU
        if(rng_.randrange(1000) > 500) {
            pair[0] = secondary_structure::ResType::ADE;
            pair[1] = secondary_structure::ResType::URA;
        }
        // selecting UA
        else {
            pair[0] = secondary_structure::ResType::URA;
            pair[1] = secondary_structure::ResType::ADE;
        }
    }
}

void
SequenceDesigner::_get_random_res_type_pair(
        secondary_structure::ResTypes & pair) {
    pair = possible_rt_bps_[rng_.randrange(possible_rt_bps_.size())];
}


}























