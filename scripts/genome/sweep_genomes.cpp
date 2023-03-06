//
// Created by aleferna on 23-02-2023.
//
#include <iostream>
#include <vector>
#include <algorithm>
#include "misc.h"
#include "genome.h"

struct GenomeTableEntry {
    double chem;
    double food;
    int tau;
    int jkey_dec;
    int jlock_dec;
    int id;
};

using namespace std;
using Input = pair<double, double>;
using Output = pair<int, pair<int, int>>;

const string genome_headers = "innr,regnr,outnr,in_scale_list,reg_threshold_list,reg_w_innode_list,"
                              "reg_w_regnode_list,out_threshold_list,out_w_regnode_list";
const string out_headers = "chem,food,tau,jkey_dec,jlock_dec,id";

vector<Genome> readGenomes(const string &inputfile) {
    ifstream file(inputfile);
    if (not file)
        throw runtime_error("Failed to open file");

    string line;
    getline(file, line);
    if (line != genome_headers)
        throw runtime_error("Inadequate file headers for genome input. Should be:\n" + genome_headers);

    vector<Genome> genomes {};
    while (getline(file, line)) {
        auto attrs = stringToVector<string>(line, ',');
        auto it = attrs.begin();
        Genome genome;
        Genome::readGenomeInfo(it, genome);
        genomes.push_back(std::move(genome));
    }

    return genomes;
}

void writeGenomeTable(vector<GenomeTableEntry> &genome_table, string &outputfile) {
    ofstream file(outputfile);
    if (not file)
        throw runtime_error("Failed to open file");

    file << out_headers << endl;
    for (auto &row : genome_table) {
        file << row.chem << ',' << row.food << ',' << row.tau << ',' << row.jkey_dec << ',' << row.jlock_dec << ','
             << row.id << endl;
    }
}

vector<Input> makeInputs(
    double min_chem,
    double max_chem,
    double step_chem,
    double min_foodparc,
    double max_foodparc,
    double step_foodparc
) {
    vector<Input> inputs {};
    // Don't use <= comparisons for the loops because of double arithmetic imprecision
    for (double chem = min_chem; chem < max_chem + step_chem; chem += step_chem)
        for (double food = min_foodparc; food < max_foodparc + step_foodparc; food += step_foodparc)
            inputs.emplace_back(chem, food);
    return inputs;
}

Output getOutput(Genome &genome, Input &input, int MCSs) {
    genome.ResetGenomeState();
    for (int i = 0; i < MCSs; ++i) {
        genome.UpdateGeneExpression(array<double, 2>{input.first, input.second}, true);
        genome.FinishUpdate();
    }
    // Tau here is slightly inaccurate because we are not accounting for the fact cells take time to enter dividing state
    return {genome.outputnodes[0].Boolstate + 1, genome.calculateJdecs()};
}

vector<GenomeTableEntry> makeGenomeTable(vector<Genome> &genomes, vector<Input> &inputs, int MCSs) {
    vector<GenomeTableEntry> genome_table {};
    for (int i = 0; i < genomes.size(); ++i) {
        auto genome = genomes[i];
        for (auto &input: inputs) {
            auto output = getOutput(genome, input, MCSs);
            GenomeTableEntry row{
                input.first,
                input.second,
                output.first,
                output.second.first,
                output.second.second,
                i
            };
            genome_table.push_back(row);
        }
    }
    return genome_table;
}

int main(int argc, char *argv[]) {
    vector<string> args(argv + 1, argv + argc);
    for (auto &arg : args) {
        if (arg == "-h" or arg == "--help") {
            cout << "Usage: sweep_genome <inputfile> <outputfile> <min_chem> "
                    "<max_chem> <step_chem> <min_foodparc> <max_foodparc> <step_foodparc> [MCSs]"
                    "\n\n"
                    "Where:\n"
                    "\t-'inputfile' must be a CSV file similar to the ones used to backup "
                    "cell data but only containing a single row and only genome attributes "
                    "(" << genome_headers << ")\n"
                    "\t-'MCSs' controls for how many time-steps the simulated genome will receive the same inputs "
                    "(default: 50, type: INT)\n"
                    "\t-Food parameters should be given in terms of division parcels: "
                    "food parcels = food / (scaling_cell_to_ca_time * (divtime + divdur) / metabperiod)\n"
                    "\t-All numeric arguments besides 'MCSs' can be doubles or integers\n"
                    "\t-The input range is inclusive in the interval [min_val, max_val]"
                    "\n\n"
                    "Description: Creates a CSV file with the output of a genome for all "
                    "combinations of inputs given by the input parameters of the program" << endl;
            return EXIT_SUCCESS;
        }
    }
    if (args.size() < 8) {
        cerr << "Inadequate arguments, try: sweep_genome -h" << endl;
        return EXIT_FAILURE;
    }
    int MCSs = 50;
    if (args.size() == 9) {
        MCSs = stoi(args[8]);
        if (MCSs < 1) {
            cerr << "'MCSs' argument must be 1 or higher" << endl;
            return EXIT_FAILURE;
        }
    }

    auto genomes = readGenomes(args[0]);
    auto inputs = makeInputs(
        stod(args[2]),
        stod(args[3]),
        stod(args[4]),
        stod(args[5]),
        stod(args[6]),
        stod(args[7])
    );
    auto genome_table = makeGenomeTable(genomes, inputs, MCSs);
    writeGenomeTable(genome_table, args[1]);
    return EXIT_SUCCESS;
}