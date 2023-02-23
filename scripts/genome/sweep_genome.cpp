//
// Created by aleferna on 23-02-2023.
//
#include <iostream>
#include "misc.h"
#include "genome.h"

using namespace std;

const string genome_attrs = "innr,regnr,outnr,in_scale_list,reg_threshold_list,reg_w_innode_list,"
                            "reg_w_regnode_list,out_threshold_list,out_w_regnode_list";

vector<Genome> readGenomes(const string &filename) {
    ifstream file(filename);
    if (not file)
        throw runtime_error("Failed to open file");

    string line;
    getline(file, line);
    if (line != genome_attrs)
        throw runtime_error("Inadequate file headers for genome input. Should be:\n" + genome_attrs);

    vector<Genome> genomes{};
    while (getline(file, line)) {
        auto attrs = stringToVector<string>(line, ',');
        auto it = attrs.begin();
        Genome genome;
        Genome::readGenomeInfo(it, genome);
        genomes.push_back(genome);
    }
    return genomes;
}

int main(int argc, char *argv[]) {
    vector<string> args(argv + 1, argv + argc);
    if (args.size() != 2)
        throw runtime_error("Inadequate arguments");

    auto genomes = readGenomes(argv[1]);
    //TODO: Implement network input sweeping function
}