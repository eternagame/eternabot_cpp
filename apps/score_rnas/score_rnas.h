//
// Created by Joseph Yesselman on 2019-11-12.
//

#ifndef RNAMAKE_NEW_SCORE_RNAS_H
#define RNAMAKE_NEW_SCORE_RNAS_H

#include <base/application.hpp>

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
    struct Parameters {
        String csv, out_file;
    };

private:
    Parameters parameters_;

};


#endif //RNAMAKE_NEW_SCORE_RNAS_H
