#include "genome.h"
#include "random.h"

//default constructor required for init in other class
Genome::Genome()
{
  innr=0;
  regnr=0;
  outnr=0;
}

Genome::Genome(int in, int reg, int out)
{
  InitGenome(in, reg, out);
}

void Genome::InitGenome(int in, int reg, int out)
{
  int i;

  innr=in;
  regnr=reg;
  outnr=out;

  //how to scale the input
  for(i=0; i<innr; i++){
    inputscale.push_back(0.1);
  }

  //for (i=0; i<in; i++)
  //{
  //  inputnodes.push_back(*(new Gene(0, i, innr, regnr)));
  //}

  for (i=0; i<reg; i++){
    //regnodes.push_back(*(new Gene(1, i+innr, innr, regnr)));
    regnodes.emplace_back(1, i+innr, innr, regnr);
  }

  for (i=0; i<out; i++){
    //outputnodes.push_back(*(new Gene(2, i+innr+regnr, innr, regnr)));
    outputnodes.emplace_back(2, i+innr+regnr, innr, regnr);
  }
}

Genome::~Genome()
{
 //cout<<"destructed"<<endl;
}
//copy constructor
Genome::Genome(const Genome &Parent)
{
  innr=Parent.innr;
  regnr=Parent.regnr;
  outnr=Parent.outnr;

  inputscale=Parent.inputscale;
  //inputnodes=Parent.inputnodes;
  regnodes=Parent.regnodes;
  outputnodes=Parent.outputnodes;
  for (auto &n: regnodes){
    n.Boolstate=0;
    n.NewBool=0;
  }
  for (auto &n: outputnodes){
    n.Boolstate=0;
    n.NewBool=0;
  }
}

void Genome::ResetGenomeState(void)
{
  for (auto &n: regnodes){
    n.Boolstate=0;
    n.NewBool=0;
  }
  for (auto &n: outputnodes){
    n.Boolstate=0;
    n.NewBool=0;
  }
}

void Genome::ReadFromFile(char *filename)
{
  std::ifstream ifs;
  string line;
  int i, nodetype, nodenr;
  double scale, thresh;

  ifs.open( filename , std::ifstream::in );

  if (ifs.is_open()){

    getline(ifs, line);
    stringstream strstr(line);

    //first read the nr of nodes
    strstr >> innr >> regnr >> outnr;

    for (i=0; i<regnr; i++){
      regnodes.push_back(*(new Gene(1, i+innr, innr, regnr)));
    }

    for (i=0; i<outnr; i++){
      outputnodes.push_back(*(new Gene(2, i+innr+regnr, innr, regnr)));
    }

    getline(ifs, line);
    strstr.clear();
    strstr.str(std::string());
    strstr<<line;
    for(i=0;i<innr; i++){
      strstr>>scale;
      inputscale.push_back(scale);
    }

    //now read all the interaction strengths
    getline(ifs, line);
    while (line.length()){
      strstr.clear();
      strstr.str(std::string());
      strstr<<line;
      //stringstream strstr(line);
      //read the straightforward cell variables from the line
      strstr>> nodetype >>nodenr>>thresh;

      //if node is regulatory
      if(nodetype==1){
        regnodes[nodenr-innr].threshold=thresh;
        for (auto &w: regnodes[nodenr-innr].w_innode) {
          strstr >> w;
        }
        for (auto &w: regnodes[nodenr-innr].w_regnode) {
          strstr >> w;
        }
       }
      //if node is output
      else{
        outputnodes[nodenr-innr-regnr].threshold=thresh;
        for (auto &w: outputnodes[nodenr-innr-regnr].w_regnode) {
          strstr >> w;
        }
      }
      //get next line
      getline(ifs, line);
    }

  }
  else {
    cerr << "Genome::ReadFromFile error: could not open file. exiting ..."<<endl;
    exit(1);
  }
}

void Genome::WriteToFile(char *filename)
{
  std::ofstream ofs;

  ofs.open( filename , std::ofstream::out);

  //first read the nr of nodes
  ofs << innr <<" "<< regnr <<" "<< outnr<<endl;
  for(auto i: inputscale){
    ofs<<i<<" ";
  }
  ofs<<endl;
  //now write all the interaction strengths

  for (auto n: regnodes){
    ofs << "1 "<< n.genenr<<" "<<n.threshold;
    for (const auto w: n.w_innode) {
      ofs << " "<< w;
    }
    for (const auto w: n.w_regnode) {
      ofs << " "<< w;
    }
    ofs<<endl;
  }
  for (auto n: outputnodes){
    ofs << "2 "<< n.genenr<<" "<<n.threshold;
    for (const auto w: n.w_regnode) {
      ofs << " "<< w;
    }
    ofs<<endl;
  }
  ofs.flush();
  ofs.close();
}

void Genome::UpdateGeneExpression(const array<double,2> &input, bool sync_cells)
{
  int i,j;
  double newval;
  vector<double> v_input;

  //transpose input data for normalisation
  for (i=0; i<input.size(); i++){
    v_input.push_back(inputscale[i]*(double)input[i]);
  }
  //cerr <<"updating..."<<endl;
  //update the regulatory genes
  for (i=0; i<regnodes.size(); i++){
    newval=0.;
    //regulation by input nodes
    for (j=0; j<regnodes[i].w_innode.size(); j++){
      newval+=v_input[j]*regnodes[i].w_innode[j];
    }
    //regulation by regulation nodes
    for (j=0; j<regnodes[i].w_regnode.size(); j++){
      newval+=(double)(regnodes[j].Boolstate)*regnodes[i].w_regnode[j];
    }
    //regulation node is on or off depending on regulation
    if(newval>regnodes[i].threshold){
      regnodes[i].NewBool=1;
    }
    else{
      regnodes[i].NewBool=0;
    }
  }

  //update output nodes
  for (i=0; i<outputnodes.size(); i++){
    newval=0.;
    for (j=0; j<outputnodes[i].w_regnode.size(); j++){
      newval+=(double)(regnodes[j].Boolstate)*outputnodes[i].w_regnode[j];
    }

    //output node is on or off depending on regulation
    if(newval>outputnodes[i].threshold){
      outputnodes[i].NewBool=1;
    }
    else{
      outputnodes[i].NewBool=0;
    }

  }

  //only update BoolState if we update cells asynchronously
  if(!sync_cells) FinishUpdate();

}

//this function puts the NewBool state of nodes into BoolState
void Genome::FinishUpdate(void)
{
  int i;
  //update the actual expression
  for (i=0; i<regnodes.size();i++){
    regnodes[i].EndCycle();
  }
  for (i=0; i<outputnodes.size();i++){
    outputnodes[i].EndCycle();
  }

}

void Genome::MutateGenome(double mu, double mustd)
{
  int i,j;

  for(auto &n: inputscale){
    if(RANDOM()<mu){
      n+=RANDNORMAL(0., mustd);
    }
  }

  for (i=0; i<regnodes.size();i++){
    regnodes[i].Mutate(mu, mustd);
    regnodes[i].Boolstate=0;
    regnodes[i].NewBool=0;
  }

  for (i=0; i<outputnodes.size(); i++){
    outputnodes[i].Mutate(mu, mustd);
    outputnodes[i].Boolstate=0;
    outputnodes[i].NewBool=0;
  }

}

//prints a dot format network to standard output.
//madness, I know...
void Genome::PrintGenome(char *filename)
{

  std::ofstream ofs;

  ofs.open( filename , std::ofstream::out);

  ofs<< "digraph G { "<<endl;
  ofs<< "layout=\"dot\""<<endl;

  ofs<<"node [fontname=\"arial\",fontsize=18,style=filled];"<<endl;
  ofs<<"bgcolor=\"#FFFFFF\";"<<endl;

  int count;
  for (auto n: regnodes){
    ofs << n.genenr<<" [label=\""<<n.genenr<<", "<<n.threshold<<"\"];"<<endl;
    count=0;
    for (const auto w: n.w_innode) {
      ofs << count <<"-> "<<n.genenr<<"[label=\" "<<w<<"\"];"<<endl;
      count++;
    }
    count=0;
    for (const auto w: n.w_regnode) {
      ofs << count+innr <<"-> "<<n.genenr<<"[label=\" "<<w<<"\"];"<<endl;
      count++;
    }
  }
  ofs <<" "<<endl;
  for (auto n: outputnodes){
    ofs << n.genenr<<" [label=\""<<n.genenr<<", "<<n.threshold<<"\"];"<<endl;
    count=0;
    for (const auto w: n.w_regnode) {
      ofs << count+innr <<"-> "<<n.genenr<<"[label=\" "<<w<<"\"];"<<endl;
      count++;
    }
  }

  ofs<<"}"<<endl;
  ofs<<""<<endl;
  ofs<<""<<endl;

  ofs.flush();
  ofs.close();
}

void Genome::OutputGeneState(void)
{
int i;
  cout<<"gene states"<<endl;

  for ( const auto &n:regnodes){
    cout <<n.genenr<<"="<<n.Boolstate<<endl;
  }
  for (const auto &n:outputnodes){
    cout <<n.genenr<<"="<<n.Boolstate<<endl;
  }

}

void Genome::GetOutput(array<int,2> &out)
{
  for (auto &n: outputnodes){
    out[n.genenr-(innr+regnr)]=n.Boolstate;
  }
}
