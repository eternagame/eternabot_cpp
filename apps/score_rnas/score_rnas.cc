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
    add_option("out_file", "eternabot.scores", base::OptionType::STRING);

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

    auto scorer = eternabot::Scorer();
    auto parser = secondary_structure::Parser();
    auto lines = base::get_lines_from_file(parameters_.csv);
    for(int i = 1; i < lines.size(); i++) {
        if(lines[i].length() < 5) { break; }
        auto spl = base::split_str_by_delimiter(lines[i], ",");
        auto p = parser.parse_to_pose(spl[0], spl[1]);


    }

    exit(0);

    /*auto designer = eternabot::SequenceDesigner();
    designer.set_option_value("steps", parameters_.steps);
    designer.setup();

    auto out = std::ofstream();
    out.open(parameters_.out_file);
    out << "opt_num,opt_score,opt_sequence,longest_gc_stretch" << std::endl;

    auto parser = secondary_structure::Parser();
    for(int i = 0; i < parameters_.n; i++) {
        auto p = parser.parse_to_pose(parameters_.seq, parameters_.ss);
        auto results = designer.design(p);
        p->replace_sequence(results[0]->sequence);
        std::cout << results[0]->score << " " << results[0]->sequence << " " << p->dot_bracket() << std::endl;
        out << i << "," << results[0]->score << "," << results[0]->sequence << ",";
        out << secondary_structure::find_longest_gc_helix_stretch(p) << std::endl;

    }*/
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
