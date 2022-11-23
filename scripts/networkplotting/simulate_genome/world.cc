#include "header.hh"
#include "genome.h"
#include "random.h"
#include "random.h"

// g++ -Wall random.cpp gene.cpp genome.cpp world.cc -o <executable>
int main(int argc, char **argv)
{
  int i,j;
  array <double,2> input; //,0.,0.
  Seed(atoi(argv[1]));

  Genome testgenome;

  testgenome.ReadFromFile(argv[2]);

  char filename[300];
  sprintf(filename, "written.dot");

  testgenome.PrintGenome(filename);
  input[0]=atof(argv[3]);
  input[1]=atof(argv[4]);

  //cout<<"input is "<<input[0]<<" "<<input[1]<<endl;

  //development
  for (j=0; j<100; j++){

    testgenome.UpdateGeneExpression(input, false);
    testgenome.OutputGeneState();
  }

  return 0;
}
