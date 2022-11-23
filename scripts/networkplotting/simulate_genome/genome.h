#ifndef GenomeHeader
#define GenomeHeader

#include <cstdio>
#include <array>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include "gene.h"

using namespace std;

class Genome
{
  public:


  int innr, regnr, outnr;

  //transformation of inputs for normalisation
  vector <double> inputscale;

  //Storage of the the GRN state
  //vector < Gene > inputnodes;
  vector < Gene > regnodes;
  vector < Gene > outputnodes;

  //constructor and destructor
  Genome();
  Genome(int in, int reg, int out);
  void InitGenome(int in, int reg, int out); //constructor body for handy init in classes
  ~Genome();

  Genome(const Genome &Parent); //copy constructor

  //read and write genome
  void ResetGenomeState(void);
  void ReadFromFile(char* filename);
  void WriteToFile(char* filename);
  void PrintGenome(char * filename); //makes a dot plot
  void OutputGeneState(void);

  void GetOutput(array<int,2> &out);
  //genomic functions
  void UpdateGeneExpression(const array<double,2> &input, bool sync_cells);
  void FinishUpdate(void);

  void MutateGenome(double mu, double mustd);

  //output functions


};


#endif
