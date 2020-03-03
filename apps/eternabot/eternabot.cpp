//
//  eternabot.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "eternabot.h"

#include <base/backtrace.hpp>
#include <base/log.h>
#include "base/cl_option.h"
#include "secondary_structure/secondary_structure_parser.h"
#include "eternabot/sequence_designer.h"

void
EternabotApp::setup_options() {
    add_option("seq", String(""), base::OptionType::STRING, false);
    add_option("ss", String(""), base::OptionType::STRING, true);
    add_option("steps", 100, base::OptionType::INT);
    add_option("n", 1, base::OptionType::INT);
    add_option("out_file", "eternabot.csv", base::OptionType::STRING);
    add_option("not_unique", false, base::OptionType::BOOL);

}

void
EternabotApp::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);
    parameters_.seq       = get_string_option("seq");
    parameters_.ss        = get_string_option("ss");
    parameters_.steps     = get_int_option("steps");
    parameters_.n         = get_int_option("n");
    parameters_.out_file  = get_string_option("out_file");
    parameters_.not_unique  = get_bool_option("not_unique");
    if(parameters_.seq == String("")) {
       for(int i = 0; i < parameters_.ss.length(); i++) {
           parameters_.seq += "N";
       }
       //std::cout << parameters_.ss << std::endl;
       //std::cout << parameters_.seq << std::endl;
    }
}

void
EternabotApp::run() {
    base::init_logging();

    auto designer = eternabot::SequenceDesigner();
    auto v = vienna::Vienna();
    designer.set_option_value("steps", parameters_.steps);
    designer.setup();

    auto out = std::ofstream();
    out.open(parameters_.out_file);
    out << "opt_num,structure,opt_score,bp_diff_score,opt_sequence,longest_gc_stretch" << std::endl;

    auto previous_solutions = Strings();
    auto parser = secondary_structure::Parser();
    for(int i = 0; i < parameters_.n; i++) {
        auto p = parser.parse_to_pose(parameters_.seq, parameters_.ss);
        designer.set_previous_solutions(previous_solutions);
        auto results = designer.design(p);
        if(!parameters_.not_unique) {
            previous_solutions.push_back(results[0]->sequence);
        }
        p->replace_sequence(results[0]->sequence);
        v.fold(results[0]->sequence);
        std::cout << results[0]->score << " " << results[0]->bp_diff_score << " " <<  results[0]->sequence << " " << p->dot_bracket() << " " << v.get_structure() << std::endl;
        out << i << "," << p->dot_bracket() << "," << results[0]->score << "," << results[0]->bp_diff_score << ",";
        out << results[0]->sequence << ",";
        out << secondary_structure::find_longest_gc_helix_stretch(p) << std::endl;

    }
}


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(base::print_backtrace);

    auto app = EternabotApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;


}
