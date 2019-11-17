//
// Created by Joseph Yesselman on 2019-11-12.
//

#include "score_rnas.h"

#include <base/file_io.h>
#include <base/backtrace.hpp>
#include <base/log.h>
#include "base/cl_option.h"
#include "secondary_structure/secondary_structure_parser.h"
#include "eternabot/sequence_designer.h"

void
ScoreRNAsApp::setup_options() {
    add_option("csv", String(""), base::OptionType::STRING, false);
    add_option("out_file", "eternabot.csv", base::OptionType::STRING);

}

void
ScoreRNAsApp::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);
    parameters_.csv       = get_string_option("csv");
    parameters_.out_file  = get_string_option("out_file");

}

void
ScoreRNAsApp::run() {
    base::init_logging();
    auto scorer = eternabot::Scorer(
            Strings{"BerexTest"},
            Floats{1});

    auto out = std::ofstream();
    out.open(parameters_.out_file);
    out << "sequence,structure,score,";
    auto strategy_names = scorer.strategy_names();
    for(int i = 0; i < strategy_names.size(); i++) {
        out << strategy_names[i];
        if(i != strategy_names.size()-1) { out << ","; }
    }
    out << std::endl;

    auto parser = secondary_structure::Parser();
    auto lines = base::get_lines_from_file(parameters_.csv);
    auto scores = Floats();
    for(int i = 1; i < lines.size(); i++) {
        if(lines[i].length() < 5) { break; }
        auto spl = base::split_str_by_delimiter(lines[i], ",");
        try {
            auto p = parser.parse_to_pose(spl[0], spl[1]);
            scorer.setup(p);
            scores = scorer.get_scores(p);
        }
        catch(...) { continue; }
        out << spl[0] << "," << spl[1] << "," << spl[2] << ",";
        for(int j = 0; j < scores.size(); j++) {
            out << scores[j];
            if(j != scores.size()-1) { out << ",";}

        }
        out << std::endl;

    }
    out.close();

}


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    auto app = ScoreRNAsApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;


}
