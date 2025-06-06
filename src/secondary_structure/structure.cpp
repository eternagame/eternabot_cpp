//
//  structure.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/structure.h"

namespace secondary_structure {

void
Structure::_setup_chains(
    String const & sequence,
    String const & dot_bracket) {
    
    if(sequence.length() != dot_bracket.length()) {
        throw Exception("cannot construct new SecondaryStructure object: new sequence and dot bracket are not the same length");
    }
    
    if(sequence.length() == 0) {
        throw Exception("cannot construct new SecondaryStructure object: sequence is of lenght zero!");
    }
    
    if(dot_bracket[0] != '(' && dot_bracket[0] != '.' &&
       dot_bracket[0] != '&' && dot_bracket[0] != '+') {
        throw Exception("cannot construct new SecondaryStructure object: dot bracket notation for secondary structure is not valid. perhaps you flipped seq and ss?");
    }
    
    ResidueOPs res;
    int count = 1;
    String chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXZ";
    auto valid_seq = String("AGUCTNWSMKRYBDHV&+-");
    auto valid_ss = String("[{(.)}]");
    int i = -1, ci = 0;
    for(auto & s : sequence) {
        i++;
        if(s != '&' && s != '+') {
            String name = "", db = "", chain_id = "";
            name += s; db += dot_bracket[i]; chain_id += chain_ids[ci];
            if(valid_seq.find(name) == std::string::npos) {
                throw Exception(
                    name + " is not a valid name for a residue valid names are: AGUCTNWSMKRYBDHV");
            }
            
            if(valid_ss.find(db) == std::string::npos) {
                throw Exception(
                    db + " is not a valid structure element for a residue valid elements are: [{(.)}]");
            }
            
            auto r = std::make_shared<Residue>(name, db, count, chain_id, util::Uuid());
            res.push_back(r);
            count++;
        }
        else {
            ci += 1;
            chains_.push_back(std::make_shared<Chain>(res));
            res = ResidueOPs();
            
            // unlikely but hit max chains
            if(ci == chain_ids.length()-1) { ci = 0; }
        }
    }
    
    if(res.size() > 0) {
        chains_.push_back(std::make_shared<Chain>(res));
    }
    
}
    

ResidueOP
Structure::get_residue(
    int const & num,
    String const & chain_id,
    String const & i_code) {
    for( auto & c : chains_) {
        for (auto & r : c->residues() ){
            if (num == r->num() &&
                chain_id == r->chain_id()  &&
                i_code == r->i_code()) {
                return r;
            }
        }
    }
    return ResidueOP(NULL);
}

ResidueOP 
Structure::get_residue(
    util::Uuid const & uuid) {
    for( auto & c : chains_) {
        for (auto & r : c->residues() ){
            if ( r->uuid() == uuid) { return r; }
        }
    }
    
    return ResidueOP(NULL);
}
 
}
