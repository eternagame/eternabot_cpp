//
// Created by Joseph Yesselman on 2019-11-12.
//

#ifndef RNAMAKE_NEW_SCORE_RNAS_H
#define RNAMAKE_NEW_SCORE_RNAS_H

#include <base/application.hpp>
#include <eternabot/scorer.h>

class ScoreRNAsApp : base::Application {
public:

    ScoreRNAsApp() : base::Application() {}

    ~ScoreRNAsApp() {}

public:

    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

private:
    void
    _setup_pairmap(
            secondary_structure::PoseOP);

    float
    _bp_list_diff(
            secondary_structure::PoseOP,
            std::vector<std::vector<int>> const &,
            size_t,
            eternabot::FeaturesOP);

private:
    struct Parameters {
        String csv, out_file;
        bool test_load;
        int start;
        int max_size;
    };

private:
    Parameters parameters_;

    std::vector<std::vector<int>> pair_map_;
    size_t pair_map_entries_;

};


#endif //RNAMAKE_NEW_SCORE_RNAS_H
