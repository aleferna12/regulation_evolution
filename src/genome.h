#ifndef GenomeHeader
#define GenomeHeader

#include "header.hh"
#include "gene.hh"

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
  void ReadFromFile(char* filename);
  void WriteToFile(char* filename);

  //genomic functions
  void UpdateGeneExpression(const vector<double> &input);

  void MutateGenome(double mu, double mustd);

  //output functions
  void OutputGenome(void);
  void OutputGeneState(void);

};


#endif
