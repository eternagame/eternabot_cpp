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
    add_option("test_load", false, base::OptionType::BOOL);
    add_option("start", 0, base::OptionType::INT);
    add_option("max_size", 99999, base::OptionType::INT);

}

void
ScoreRNAsApp::parse_command_line(
        int argc,
        const char ** argv) {
    base::Application::parse_command_line(argc, argv);
    parameters_.csv       = get_string_option("csv");
    parameters_.out_file  = get_string_option("out_file");
    parameters_.test_load = get_bool_option("test_load");
    parameters_.start     = get_int_option("start");
    parameters_.max_size  = get_int_option("max_size");

}

void
ScoreRNAsApp::run() {
    base::init_logging();
    auto scorer = eternabot::Scorer();
    //auto scorer = eternabot::Scorer(
    //        Strings{"ModifiedBerexTest"},
    //        Floats{1});

    auto out = std::ofstream();
    out.open(parameters_.out_file);
    if(!parameters_.test_load) {
        out << "sequence,structure,score,";
        auto strategy_names = scorer.strategy_names();
        for (int i = 0; i < strategy_names.size(); i++) {
            out << strategy_names[i];
            if (i != strategy_names.size() - 1) { out << ","; }
        }
        out << ",total_score,bp_diff,predicted_structure";
        out << std::endl;
    }
    else {
        out << "sequence,structure,name,pass" << std::endl;
    }

    auto parser = secondary_structure::Parser();
    auto lines = base::get_lines_from_file(parameters_.csv);
    auto scores = Floats();
    for(int i = 1; i < lines.size(); i++) {
        if(i < parameters_.start) {
            continue;
        }
        //std::cout << i << std::endl;
        if(lines[i].length() < 5) { break; }
        scorer = eternabot::Scorer();
        auto score = 0.0;
        auto spl = base::split_str_by_delimiter(lines[i], ",");
        if(spl[0].length() > parameters_.max_size) {
            continue;
        }
        auto p = secondary_structure::PoseOP(nullptr);
        try {
            p = parser.parse_to_pose(spl[0], spl[1]);
            if(parameters_.test_load) {
                out << spl[0] << "," << spl[1] << "," << spl[2] << "," << 1 << std::endl;
                continue;

            }
            scorer.setup(p);
        }
        catch(...) {
            if(parameters_.test_load) {
                out << spl[0] << "," << spl[1] << "," << spl[2] << "," << 0 << std::endl;
            }

            continue;
        }

        try {
            score = scorer.score_secondary_structure(p);
            scores = scorer.get_scores(p);
            _setup_pairmap(p);
        }
        catch(secondary_structure::Exception const & e) {
            std::cout << spl[0] << " " << spl[1] << e.what() << std::endl;
            continue;
        }


        out << spl[0] << "," << spl[1] << "," << spl[2] << ",";
        for(int j = 0; j < scores.size(); j++) {
            out << scores[j];
            if(j != scores.size()-1) { out << ",";}

        }
        out << "," << score << "," << _bp_list_diff(p, pair_map_, pair_map_entries_, scorer.features());
        out << "," << scorer.features()->structure;
        out << std::endl;

    }
    out.close();

}

void
ScoreRNAsApp::_setup_pairmap(
        secondary_structure::PoseOP p) {
    pair_map_ = std::vector<std::vector<int>>(p->residues().size()+1);
    for(int i = 0; i < p->residues().size()+1; i++) {
        pair_map_[i] = std::vector<int>(p->residues().size()+1);
    }
    for(auto const & bp : p->basepairs()) {
        pair_map_[bp->res1()->num()][bp->res2()->num()] = 1;
        pair_map_[bp->res2()->num()][bp->res1()->num()] = 1;
    }
    pair_map_entries_ = p->residues().size()*p->residues().size();
}

float
ScoreRNAsApp::_bp_list_diff(
        secondary_structure::PoseOP p,
        std::vector<std::vector<int>> const & pair_map,
        size_t pair_list_size,
        eternabot::FeaturesOP features) {

    int pi = 0, pj = 0;
    auto score = 0.0f;
    auto & plist = features->dotplot;
    for(int i = 0 ; i < pair_list_size; i++) {
        if(plist[i].p < 0.001) { continue; }
        pi = plist[i].i;
        pj = plist[i].j;
        score += abs(pair_map[pi][pj] - plist[i].p);
    }
    return score;

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
