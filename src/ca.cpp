/*

Copyright 1995-2006 Roeland Merks, Nick Savill

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/


/* CA.cpp: implementation of Glazier & Graner's Cellular Potts Model */

// This code derives from a Cellular Potts implementation written around 1995
// by Nick Savill

#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include "sticky.h"
#include "random.h"
#include "ca.h"
#include "parameter.h"
#include "dish.h"
#include "sqr.h"
#include "crash.h"
#include "hull.h"

#define ZYGFILE(Z) <Z.xpm>
#define XPM(Z) Z ## _xpm
#define ZYGXPM(Z) XPM(Z)

//Leonie's sign function
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/* define default zygote */
/* NOTE: ZYGOTE is normally defined in Makefile!!!!!! */
// #ifndef ZYGOTE
// #define ZYGOTE init
// //#include "init.xpm"
// #else
// #include ZYGFILE(ZYGOTE)
// #endif

/* STATIC DATA MEMBER INITIALISATION */
double copyprob[BOLTZMANN];

const int CellularPotts::nx[25] = {0, 0, 1, 0,-1, 1, 1,-1,-1, 0, 2, 0, -2, 1, 2, 2, 1,-1,-2,-2,-1, 0, 2, 0,-2 };
const int CellularPotts::ny[25] = {0,-1, 0, 1, 0,-1, 1, 1,-1,-2, 0, 2,  0,-2,-1, 1, 2, 2, 1,-1,-2,-2, 0, 2, 0 };

const int CellularPotts::nbh_level[5] = { 0, 4, 8, 20, 24 }; //0:self; 1:van Neumann; 2:Moore; etc...
int CellularPotts::shuffleindex[9]={0,1,2,3,4,5,6,7,8};

extern Parameter par;


/** PRIVATE **/

using namespace std;
void CellularPotts::BaseInitialisation(vector<Cell> *cells) {
  CopyProb(par.T);
  cell=cells;
  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
  else
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4]).";

}

CellularPotts::CellularPotts(vector<Cell> *cells,
			     const int sx, const int sy) {

  sigma=0;
  frozen=false;
  thetime=0;
  zygote_area=0;


  BaseInitialisation(cells);
  sizex=sx;
  sizey=sy;

  AllocateSigma(sx,sy);


  // fill borders with special border state
  for (int x=0;x<sizex;x++) {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey;y++) {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }

  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
  else
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";
}

CellularPotts::CellularPotts(void) {

  sigma=0;
  sizex=0; sizey=0;
  frozen=false;
  thetime=0;
  zygote_area=0;

  CopyProb(par.T);

  // fill borders with special border state
  for (int x=0;x<sizex;x++) {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey;y++) {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }
  if (par.neighbours>=1 && par.neighbours<=4)
    n_nb=nbh_level[par.neighbours];
  else
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";
}

// destructor (virtual)
CellularPotts::~CellularPotts(void) {
  if (sigma) {
    free(sigma[0]);
    free(sigma);
    sigma=0;
  }
}

void CellularPotts::AllocateSigma(int sx, int sy) {

  sizex=sx; sizey=sy;

  sigma=(int **)malloc(sizex*sizeof(int *));
  if (sigma==NULL)
    MemoryWarning();

  sigma[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (sigma[0]==NULL)
    MemoryWarning();


  {for (int i=1;i<sizex;i++)
    sigma[i]=sigma[i-1]+sizey;}

  /* Clear CA plane */
   {for (int i=0;i<sizex*sizey;i++)
     sigma[0][i]=0; }

}

//this function is used in ReadBackup in Dish to fill the ca plane
//AND to set the cell's moments (assumes cells have been initialised!)
int CellularPotts::SetNextSigma(int sig) {
  //the plane has a 1 px boundary on all size, therefore we place the pixels
  //within that
  static int xcount=1, ycount=1;

  if(xcount>=sizex-1 ||ycount>=sizey-1){
    return 1;
  }

  sigma[xcount][ycount]=sig;
  if (sig){
    (*cell)[sig].area++;
    (*cell)[sig].AddSiteToMoments(xcount, ycount);
  }

  ycount++;
  if(ycount==sizey-1){
    ycount=1;
    xcount++;
  }
  return 0;

}

void CellularPotts::IndexShuffle() {

  int i;
  int temp;
  int index1,index2;

  for (i=0;i<9;i++) {

    index1=RandomNumber(8);
    index2=RandomNumber(8);

    temp=shuffleindex[index1];
    shuffleindex[index1]=shuffleindex[index2];
    shuffleindex[index2]=temp;

  }
}


double sat(double x) {

  return x/(par.saturation*x+1.);
  //return x;

}

void CellularPotts::InitializeEdgeList(bool init_single_cell_center){
  // unordered_set<long long int> edgeSet;
	int neighbour;
	int x,y, xp, yp;
	int c, cp;
  int xstart =0;
  int xend=sizex-1;
  int ystart=0;
  int yend=sizey-1;

//  bool init_single_cell_center=true;

  //cout<<"edgevector is size "<<edgeSetVector.size_map()<<endl;;

  if (init_single_cell_center && sizex>200 && sizey>200){
    xstart=(sizex-1)/2-100;
    xend=(sizex-1)/2+100;
    ystart=(sizey-1)/2-100;
    yend=(sizey-1)/2+100;
  }


	for (x=xstart; x< xend; x++){
    for ( y=ystart; y<yend; y++){
      for (int nb=1; nb<=n_nb; nb++){
    		c = sigma[x][y]; //c = sigma here
    		xp = nx[nb]+x;  //check neighbour
    		yp = ny[nb]+y;

		      if (par.periodic_boundaries) {

          // since we are asynchronic, we cannot just copy the borders once
          // every MCS

    			if (xp<=0)
    				xp=sizex-2+xp;
        	if (yp<=0)
    				yp=sizey-2+yp;
       	  if (xp>=sizex-1)
    				xp=xp-sizex+2;
       	  if (yp>=sizey-1)
    				yp=yp-sizey+2;

          cp=sigma[xp][yp];
          }
    		else if (xp<=0 || yp<=0 || xp>=sizex-1 || yp>=sizey-1)
    			cp=-1;
        else
    			cp=sigma[xp][yp]; //if boundaries deal appropriately, else set the neighbour
    //if c != cp and cp is not boundary and sgn(x) is the sign function
    // (check for boundary state in the neigh.)
    if (cp != c && cp != -1 && c!=-1 && sgn(cp)+sgn(c)>=0){
      //edgeSetpair.insert({x,y,xp,yp});
      edgeSetVector.insert({x,y,xp,yp});
      // cout << edgeSetVector.subset_size() <<endl;
      // edgeSetVector.push_back({x,y,xp,yp});
      // edgeSetVector.insert_in_subset(edgeSetVector.subset_size());
      // cout << edgeSetVector.subset_size() << endl;
      //edgeSetpair.insert({xp,yp,x,y});
      edgeSetVector.insert({xp,yp,x,y});
      // edgeSetVector.push_back({xp,yp,x,y});
      // edgeSetVector.insert_in_subset(edgeSetVector.subset_size());
      // edgeVector.push_back({x,y,xp,yp});
		}
	}
}}
}

int CellularPotts::DeltaHWithMedium(int x,int y, PDE *PDEfield)
{
  double DH = 0.;
  int i, sxy, sxyp;
  int neighsite;

  /* Compute energydifference *IF* the copying were to occur */
  sxy = sigma[x][y];
  sxyp = MEDIUM; // sigma[xp][yp];     // ********** THIS IS MEDIUM

  // DH due to cell adhesion
  // for all neighbours n_i we do Delta( nrg(medium),nrg(n_i)) - Delta(nrg(focal),nrg(n_i) )
  // and sum everything, this amounts to calculate ( Nrg(future) - Nrg(now) )
  for (i=1;i<=n_nb;i++) {
    int xp2,yp2;
    xp2=x+nx[i]; yp2=y+ny[i];
    if (par.periodic_boundaries) {

      // since we are asynchronic, we cannot just copy the borders once every MCS
      //actually, haven't we dealt with this already in AmoebaeMove? No, this is about all the other neighbours

      if (xp2<=0) xp2=sizex-2+xp2;
      if (yp2<=0) yp2=sizey-2+yp2;
      if (xp2>=sizex-1) xp2=xp2-sizex+2;
      if (yp2>=sizey-1) yp2=yp2-sizey+2;

      neighsite=sigma[xp2][yp2];


    } else {
      if (xp2<=0 || yp2<=0 || xp2>=sizex-1 || yp2>=sizey-1)
        neighsite=-1;
      else
        neighsite=sigma[xp2][yp2];
    }

    if (neighsite==-1) { // border - for "fixed" boundary conditions
      DH += (sxyp==0?0:par.border_energy)-(sxy==0?0:par.border_energy);
    } else {
      //DH += (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) - (*cell)[sxy].EnergyDifference((*cell)[neighsite]);
      // notice that sxyp is medium, so there is no need of calling the function.
      DH += (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) - Adhesion_Energy(sxy , neighsite);
    }
  }



  // lambda is determined by chemical 0

  //cerr << "[" << lambda << "]";



  //THIS IS THE ONLY CASE HERE -->  if ( sxyp == MEDIUM ) {
  DH += (int)(par.lambda *  (1. - 2. * (double) ( (*cell)[sxy].Area() - (*cell)[sxy].TargetArea()) ));

    //}
  //else if ( sxy == MEDIUM ) {
  //  DH += (int)((par.lambda * (1. + 2. * (double) ( (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()) )));
  // }
  //else
   // DH += (int)((par.lambda * (2.+  2.  * (double) (  (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()
   // - (*cell)[sxy].Area() + (*cell)[sxy].TargetArea() )) ));


  // DH is determined at this point in cell_evolution, and the following does not apply - at least for now.

  //NO chemotaxis in cell_evolution... or at least not yet
  /* Chemotaxis */
//   if (PDEfield && (par.vecadherinknockout || (sxyp==0 || sxy==0))) {
//
//     // copying from (xp, yp) into (x,y)
//     // If par.extensiononly == true, apply CompuCell's method, i.e.
//     // only chemotactic extensions contribute to energy change
//     if (!( par.extensiononly && sxyp==0)) {
//       int DDH=(int)(par.chemotaxis*(sat(PDEfield->Sigma(0,x,y))-sat(PDEfield->Sigma(0,xp,yp))));
//       DH-=DDH;
//     }
//   }



  // WHEN YOU ARE DONE DEBUGGING THE ANISOTROPY OF THE CELLS, THIS CAN BE RESTORED
  // lambda part and cell migr.

//   const double lambda2=par.lambda2; // par.lambda2=0 in cell_evolution, so the following doesn't apply
//
//   /* Length constraint */
//   // sp is expanding cell, s is retracting cell
//
//
//   if ( sxyp == MEDIUM ) {
//     DH -= (int)(lambda2*( DSQR((*cell)[sxy].Length()-(*cell)[sxy].TargetLength())
//     - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x,y) -
//     (*cell)[sxy].TargetLength()) ));
//
//   }
//   else if ( sxy == MEDIUM ) {
//     DH -= (int)(lambda2*(DSQR((*cell)[sxyp].Length()-(*cell)[sxyp].TargetLength())
//     -DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x,y)-(*cell)[sxyp].TargetLength())));
//
//
//   }
//   else {
//     DH -= (int)(lambda2*( (DSQR((*cell)[sxyp].Length()-(*cell)[sxyp].TargetLength())
//     -DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x,y)-(*cell)[sxyp].TargetLength())) +
//     ( DSQR((*cell)[sxy].Length()-(*cell)[sxy].TargetLength())
//     - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x,y) -
//     (*cell)[sxy].TargetLength()) ) ));
//   }
//
//
   /*cell migration */
  //Joost's method
  double ax, ay;
   if((*cell)[sxy].getMu()>0.0001 || (*cell)[sxyp].getMu()>0.0001){
       double smeanx = (*cell)[sxy].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other side
       double smeany = (*cell)[sxy].getYpos();

       if(par.periodic_boundaries){
         if( (x-smeanx)>0 && (x-smeanx)>(smeanx-(x-(par.sizex-2))) ) {
           smeanx+=(par.sizex-2);
//          cerr<<"dhwm hello"<<endl;
           //cerr<<"passb"<<endl;
         }
         else if( (smeanx-x)>0 && (smeanx-x)>(x+(par.sizex-2)-smeanx) ){
           smeanx-=(par.sizex-2);
//            cerr<<"dhwm hello"<<endl;
           //cerr<<"passb"<<endl;
        }
         if( (y-smeany)>0 && (y-smeany)>(smeany-(y-(par.sizey-2))) ){
           smeany+=(par.sizey-2);
           //cerr<<"passb"<<endl;

        }
         else if( (smeany-y)>0 && (smeany-y)>(y+(par.sizey-2)-smeany) ){
           smeany-=(par.sizey-2);
           //cerr<<"passb"<<endl;
         }
       }

       ax=x-smeanx;
       ay=y-smeany;
       DH+=(*cell)[sxy].getMu()*(ax*(*cell)[sxy].getXvec() + ay*(*cell)[sxy].getYvec())/hypot(ax,ay);


     // cout << "Migrating1!"<<endl;
     //cout<< sxy<<" "<<(*cell)[sxy].getMu()<<" "<<sxyp<<" "<<(*cell)[sxyp].getMu()<<endl;

       // WAS THIS BELOW, BUGGY WITH WRAPPED BOUNDARIES
       //ax=x-(*cell)[sxy].getXpos();
	//ay=y-(*cell)[sxy].getYpos();
	//DH+=(*cell)[sxy].getMu()*(ax*(*cell)[sxy].getXvec() + ay*(*cell)[sxy].getYvec())/hypot(ax,ay);
   }

//Similarly to Joost's method, a bias due to chemokine gradient
    if((*cell)[sxy].getChemMu()>0.0001 || (*cell)[sxyp].getChemMu()>0.0001){
        double smeanx = (*cell)[sxy].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other side
        double smeany = (*cell)[sxy].getYpos();

        if(par.periodic_boundaries){
          if( (x-smeanx)>0 && (x-smeanx)>(smeanx-(x-(par.sizex-2))) ) {
            smeanx+=(par.sizex-2);
   //          cerr<<"dhwm hello"<<endl;
            //cerr<<"passb"<<endl;
          }
          else if( (smeanx-x)>0 && (smeanx-x)>(x+(par.sizex-2)-smeanx) ){
            smeanx-=(par.sizex-2);
   //            cerr<<"dhwm hello"<<endl;
            //cerr<<"passb"<<endl;
         }
          if( (y-smeany)>0 && (y-smeany)>(smeany-(y-(par.sizey-2))) ){
            smeany+=(par.sizey-2);
            //cerr<<"passb"<<endl;

         }
          else if( (smeany-y)>0 && (smeany-y)>(y+(par.sizey-2)-smeany) ){
            smeany-=(par.sizey-2);
            //cerr<<"passb"<<endl;
          }
        }

        ax=x-smeanx;
        ay=y-smeany;
        DH+=(*cell)[sxy].getChemMu()*(ax*(*cell)[sxy].getChemXvec() + ay*(*cell)[sxy].getChemYvec())/hypot(ax,ay);


      // cout << "Migrating1!"<<endl;
      //cout<< sxy<<" "<<(*cell)[sxy].getMu()<<" "<<sxyp<<" "<<(*cell)[sxyp].getMu()<<endl;

        // WAS THIS BELOW, BUGGY WITH WRAPPED BOUNDARIES
        //ax=x-(*cell)[sxy].getXpos();
   //ay=y-(*cell)[sxy].getYpos();
   //DH+=(*cell)[sxy].getMu()*(ax*(*cell)[sxy].getXvec() + ay*(*cell)[sxy].getYvec())/hypot(ax,ay);
    }
  return DH;
}



int CellularPotts::DeltaH(int x,int y, int xp, int yp, PDE *PDEfield)
{
  double DH = 0.;
  int i, sxy, sxyp;
  int neighsite;

  /* Compute energydifference *IF* the copying were to occur */
  sxy = sigma[x][y];
  sxyp = sigma[xp][yp];

  // cerr<< "x ,y : " << x<<" "<<y<<endl;
  // cerr<< "xp,yp: " << xp<<" "<<yp<<endl;
  //
  /* DH due to cell adhesion */
  for(i=1;i<=n_nb;i++) {
    int xp2,yp2;
    xp2=x+nx[i]; yp2=y+ny[i];
    if (par.periodic_boundaries) {
      // since we are asynchronic, we cannot just copy the borders once
      // every MCS
      //actually, haven't we dealt with this already in AmoebaeMove? No, this is about all the other neighbours
      if (xp2<=0) xp2=sizex-2+xp2;
      if (yp2<=0) yp2=sizey-2+yp2;
      if (xp2>=sizex-1) xp2=xp2-sizex+2;
      if (yp2>=sizey-1) yp2=yp2-sizey+2;
      neighsite=sigma[xp2][yp2];


    }else{
      if (xp2<=0 || yp2<=0 || xp2>=sizex-1 || yp2>=sizey-1)
        neighsite=-1;
      else
        neighsite=sigma[xp2][yp2];
    }

    if (neighsite==-1) { // border - for "fixed" boundary conditions
      DH += (sxyp==0?0:par.border_energy)-(sxy==0?0:par.border_energy);
    } else {
      //cerr<<"Hello DH 0.1"<<endl;

      // MAKE A FUNCTION that takes all values: J's and regulation and sigmas
      //get to the point where you write:
      //DH += Adhesion_Energy( s1 , s2, J value)
      //    - Adhesion_Energy( same but for the other pair)
      //DH += Adhesion_Energy( (*cell)[neighsite].GetExtProtExpress_Fraction() , (*cell)[sxyp].GetExtProtExpress_Fraction() , sxyp , neighsite, (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) )
      //    - Adhesion_Energy( (*cell)[neighsite].GetExtProtExpress_Fraction() , (*cell)[sxy].GetExtProtExpress_Fraction() , sxy , neighsite, (*cell)[sxy].EnergyDifference((*cell)[neighsite]) )
      // cerr<< "sigmas sxyp, neigh: "<<sxyp<<" "<<neighsite<<": "<< Adhesion_Energy(sxyp , neighsite)<<endl;
      // cerr<< "sigmas sxy, neigh: "<<sxy<<" "<<neighsite<<": "<< Adhesion_Energy(sxy , neighsite)<<endl;
      // cerr<< "DH += " << (int) (Adhesion_Energy(sxyp , neighsite) - Adhesion_Energy(sxy , neighsite))<<endl ;
      DH += Adhesion_Energy(sxyp , neighsite)
          - Adhesion_Energy(sxy , neighsite);
      // NOTICE THAT DH is an integer... dammit! This can create problems when multiplying with a factor
      //DH += (*cell)[sxyp].EnergyDifference((*cell)[neighsite])
      //      - (*cell)[sxy].EnergyDifference((*cell)[neighsite]);

    }
  }

  // cerr<<"End of Adh, DH = "<< DH<<endl;

  // lambda is determined by chemical 0

  //cerr << "[" << lambda << "]";
  if ( sxyp == MEDIUM ) {
    DH += (int)(par.lambda *  (1. - 2. * (double) ( (*cell)[sxy].Area() - (*cell)[sxy].TargetArea()) ));
  }
  else if ( sxy == MEDIUM ) {
    DH += (int)((par.lambda * (1. + 2. * (double) ( (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()) )));
  }
  else
    DH += (int)((par.lambda * (2.+  2.  * (double)(  (*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea() - (*cell)[sxy].Area() + (*cell)[sxy].TargetArea() )) ));


  // DH is determined at this point in cell_evolution, and the following does not apply - at least for now.

  //NO chemotaxis in cell_evolution... or at least not yet
  /* Chemotaxis */
//   if (PDEfield && (par.vecadherinknockout || (sxyp==0 || sxy==0))) {
//
//     // copying from (xp, yp) into (x,y)
//     // If par.extensiononly == true, apply CompuCell's method, i.e.
//     // only chemotactic extensions contribute to energy change
//     if (!( par.extensiononly && sxyp==0)) {
//       int DDH=(int)(par.chemotaxis*(sat(PDEfield->Sigma(0,x,y))-sat(PDEfield->Sigma(0,xp,yp))));
//       DH-=DDH;
//     }
//   }


  // WHENEVER YOU'LL NEED THIS, UNCOMMENT - BUT DO NOT USE IT TOGETHER WITH periodic_boundaries

  //const double lambda2=par.lambda2; // par.lambda2=0 in cell_evolution, so the following doesn't apply
  /* Length constraint */
  // sp is expanding cell, s is retracting cell

  /*
  if ( sxyp == MEDIUM ) {
    DH -= (int)(lambda2*( DSQR((*cell)[sxy].Length()-(*cell)[sxy].TargetLength())
		       - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x,y) -
			      (*cell)[sxy].TargetLength()) ));

  }
  else if ( sxy == MEDIUM ) {
    DH -= (int)(lambda2*(DSQR((*cell)[sxyp].Length()-(*cell)[sxyp].TargetLength())
			 -DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x,y)-(*cell)[sxyp].TargetLength())));


  }
  else {
    DH -= (int)(lambda2*( (DSQR((*cell)[sxyp].Length()-(*cell)[sxyp].TargetLength())
		     -DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x,y)-(*cell)[sxyp].TargetLength())) +
		    ( DSQR((*cell)[sxy].Length()-(*cell)[sxy].TargetLength())
		      - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x,y) -
			     (*cell)[sxy].TargetLength()) ) ));
  }*/

  // WAS LIKE THIS BEFORE, BUGGY WITH WRAPPED BOUNDARIES
//   if((*cell)[sxy].getMu()>0.0001 || (*cell)[sxyp].getMu()>0.0001){
//     //cout << "Migrating!"<<endl;
//     //cerr<<"sigma focal="<< sxy<<" "<<(*cell)[sxy].getXvec()<<" "<<(*cell)[sxy].getYvec()<<" ";
//     //cerr<<"sigma copy="<< sxyp<<" "<<(*cell)[sxyp].getXvec()<<" "<<(*cell)[sxyp].getYvec()<<endl;
//
//     if(sxy!=MEDIUM){
//       ax=x-(*cell)[sxy].getXpos(); //getXpos() returns meanx WHICH i HAVE
//       ay=y-(*cell)[sxy].getYpos(); //getYpos() returns meany
//       DH+=(*cell)[sxy].getMu()*(ax*(*cell)[sxy].getXvec() + ay*(*cell)[sxy].getYvec())/hypot(ax,ay);
//     }
//     if(sxyp!=MEDIUM){
//       ax=x-(*cell)[sxyp].getXpos(); //returns meanx
//       ay=y-(*cell)[sxyp].getYpos(); //returns meany
//       DH-=(*cell)[sxyp].getMu()*(ax*(*cell)[sxyp].getXvec() + ay*(*cell)[sxyp].getYvec())/hypot(ax,ay);
//     }
//   }
/*cell migration */
//Joost's method
double ax, ay;
  if((*cell)[sxy].getMu()>0.0001 || (*cell)[sxyp].getMu()>0.0001){
    if(sxy!=MEDIUM){
      //cerr<<"tvecx: "<<(*cell)[sxy].getXvec()<<", tvecy: "<< (*cell)[sxy].getYvec() <<endl;
      double smeanx = (*cell)[sxy].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other side
      double smeany = (*cell)[sxy].getYpos();

      if(par.periodic_boundaries){
        // if x is on the right and meanx is on the left
        // and if by moving meanx to the right we diminish this distance
        if( (x-smeanx)>0 && (x-smeanx)>(smeanx+(par.sizex-2)-x) ){
          smeanx+=(par.sizex-2);
          //cerr<<"dh s hello"<<endl;
          //cerr<<"passb"<<endl;
        }
        else if( (smeanx-x)>0 && (smeanx-x)>(x+(par.sizex-2)-smeanx) ) {
          smeanx-=(par.sizex-2);
//            cerr<<"dh s hello"<<endl;
          //cerr<<"passb"<<endl;
        }
        if( (y-smeany)>0 && (y-smeany)>(smeany-(y-(par.sizey-2))) ){
          smeany+=(par.sizey-2);
          //cerr<<"passb"<<endl;
        }
        else if( (smeany-y)>0 && (smeany-y)>(y+(par.sizey-2)-smeany) ){
          smeany-=(par.sizey-2);
          //cerr<<"passb"<<endl;
        }
      }

      ax=x-smeanx;
      ay=y-smeany;
      DH+=(*cell)[sxy].getMu()*(ax*(*cell)[sxy].getXvec() + ay*(*cell)[sxy].getYvec())/hypot(ax,ay);
    }
    if(sxyp!=MEDIUM){
       double spmeanx = (*cell)[sxyp].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other side
       double spmeany = (*cell)[sxyp].getYpos();

       if(par.periodic_boundaries){
         if( (x-spmeanx)>0 && (x-spmeanx)>(spmeanx-(x-(par.sizex-2))) ){
           spmeanx+=(par.sizex-2);
//            cerr<<"dh sp hello"<<endl;
           //cerr<<"passb"<<endl;
        }
         else if( (spmeanx-x)>0 && (spmeanx-x)>(x+(par.sizex-2)-spmeanx) ) {
           spmeanx-=(par.sizex-2);
//            cerr<<"dh sp hello"<<endl;
           //cerr<<"passb"<<endl;
         }
         if( (y-spmeany)>0 && (y-spmeany)>(spmeany-(y-(par.sizey-2))) ){
           spmeany+=(par.sizey-2);
           //cerr<<"passb"<<endl;
         }
         else if( (spmeany-y)>0 && (spmeany-y)>(y+(par.sizey-2)-spmeany) ){
           spmeany-=(par.sizey-2);
           //cerr<<"passb"<<endl;
         }
       }
       //ax=x-(*cell)[sxyp].getXpos(); //returns meanx
       //ay=y-(*cell)[sxyp].getYpos(); //returns meany
       ax=x-spmeanx;
       ay=y-spmeany;
       DH-=(*cell)[sxyp].getMu()*(ax*(*cell)[sxyp].getXvec() + ay*(*cell)[sxyp].getYvec())/hypot(ax,ay);
    }
  }


  //Similar to Joost's method, but for chemotaxis (no persistence!)
    if((*cell)[sxy].getChemMu()>0.0001 || (*cell)[sxyp].getChemMu()>0.0001){
      if(sxy!=MEDIUM){
        //cerr<<"tvecx: "<<(*cell)[sxy].getXvec()<<", tvecy: "<< (*cell)[sxy].getYvec() <<endl;
        double smeanx = (*cell)[sxy].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other side
        double smeany = (*cell)[sxy].getYpos();

        if(par.periodic_boundaries){
          // if x is on the right and meanx is on the left
          // and if by moving meanx to the right we diminish this distance
          if( (x-smeanx)>0 && (x-smeanx)>(smeanx+(par.sizex-2)-x) ){
            smeanx+=(par.sizex-2);
            //cerr<<"dh s hello"<<endl;
            //cerr<<"passb"<<endl;
          }
          else if( (smeanx-x)>0 && (smeanx-x)>(x+(par.sizex-2)-smeanx) ) {
            smeanx-=(par.sizex-2);
  //            cerr<<"dh s hello"<<endl;
            //cerr<<"passb"<<endl;
          }
          if( (y-smeany)>0 && (y-smeany)>(smeany-(y-(par.sizey-2))) ){
            smeany+=(par.sizey-2);
            //cerr<<"passb"<<endl;
          }
          else if( (smeany-y)>0 && (smeany-y)>(y+(par.sizey-2)-smeany) ){
            smeany-=(par.sizey-2);
            //cerr<<"passb"<<endl;
          }
        }

        ax=x-smeanx;
        ay=y-smeany;
        DH+=(*cell)[sxy].getChemMu()*(ax*(*cell)[sxy].getChemXvec() + ay*(*cell)[sxy].getChemYvec())/hypot(ax,ay);
      }
      if(sxyp!=MEDIUM){
         double spmeanx = (*cell)[sxyp].getXpos(); //getXpos() returns meanx - which I have to wrap if pixel's on the other side
         double spmeany = (*cell)[sxyp].getYpos();

         if(par.periodic_boundaries){
           if( (x-spmeanx)>0 && (x-spmeanx)>(spmeanx-(x-(par.sizex-2))) ){
             spmeanx+=(par.sizex-2);
  //            cerr<<"dh sp hello"<<endl;
             //cerr<<"passb"<<endl;
          }
           else if( (spmeanx-x)>0 && (spmeanx-x)>(x+(par.sizex-2)-spmeanx) ) {
             spmeanx-=(par.sizex-2);
  //            cerr<<"dh sp hello"<<endl;
             //cerr<<"passb"<<endl;
           }
           if( (y-spmeany)>0 && (y-spmeany)>(spmeany-(y-(par.sizey-2))) ){
             spmeany+=(par.sizey-2);
             //cerr<<"passb"<<endl;
           }
           else if( (spmeany-y)>0 && (spmeany-y)>(y+(par.sizey-2)-spmeany) ){
             spmeany-=(par.sizey-2);
             //cerr<<"passb"<<endl;
           }
         }
         //ax=x-(*cell)[sxyp].getXpos(); //returns meanx
         //ay=y-(*cell)[sxyp].getYpos(); //returns meany
         ax=x-spmeanx;
         ay=y-spmeany;
         DH-=(*cell)[sxyp].getChemMu()*(ax*(*cell)[sxyp].getChemXvec() + ay*(*cell)[sxyp].getChemYvec())/hypot(ax,ay);
      }
    }
//
//
  return DH;
}

double CellularPotts::Adhesion_Energy(int sigma1, int sigma2)
{
  double Jval = (*cell)[sigma1].EnergyDifference((*cell)[sigma2]);
  if(sigma1==sigma2) return 0;
  if(! (sigma1 && sigma2) )
    return Jval;

  //double fr1 = (*cell)[sigma1].GetExtProtExpress_Fraction();
  //double fr2 = (*cell)[sigma2].GetExtProtExpress_Fraction();
  //double minfr = (fr1<fr2)?fr1:fr2;
  // minfr is in [0,1],
  // I want to map it so that if minfr=0, I return 43, and if minfr=1 I return Jval
  // very simple function would be linear:
  // a line between (0,43) and (1,Jval) looks like:
  // y - 43 = (Jval-43)/(1-0) * (minfr1 - 0) ; i.e.
  // y = 43 + minfr*(Jval-43)

//  double  Jtoreturn = 43. - minfr*(43. - Jval);
  double Jtoreturn=Jval;

  return Jtoreturn;
}

bool CellularPotts::Probability(int DH)
{
  if ( DH > BOLTZMANN-1 )
    return false;
  else if ( DH < 0 || RANDOM() < copyprob[DH] )
    return true;
   return false;
}


//this the function that changes the sigma
// it also does book-keeping of everything
void CellularPotts::ConvertSpinToMedium(int x,int y)
{
  int tmpcell;
  if ( (tmpcell=sigma[x][y]) ) { // if tmpcell is not MEDIUM - this should be excluded outside, so it is good check
    (*cell)[tmpcell].DecrementArea(); // returns --area
    (*cell)[tmpcell].RemoveSiteFromMoments(x,y);

    //if cell's area is zero -> kill it
    if (!(*cell)[tmpcell].Area()) {
      (*cell)[tmpcell].Apoptose();
      //cerr << "Cell " << tmpcell << " apoptosed\n";
    }
  }else{
    cerr<<"ConvertSpinToMedium(): Error. Trying to convert spin of medium with medium"<<endl;
  }

//   if ( (tmpcell=sigma[xp][yp]) ) {// if tmpcell is not MEDIUM, it gains one pixel
//     (*cell)[tmpcell].IncrementArea();
//     (*cell)[tmpcell].AddSiteToMoments(x,y);
//
//   }

  // Book-keeping of contacts
  int check=0;
  int point;
  int idSelf=sigma[x][y];     // Surely not medium <- we check this in this special case
  int idNeigh=MEDIUM;  // MEDIUM IN MY CASE

  ///use this if you need info about contacts with neighbours
  // go through neighbourhood
  for(int k=1; k<=n_nb; k++){
    // Do you not update contacts at boundaries???
    int neix=x+nx[k];
    int neiy=y+ny[k];
    if(neix<=0 || neix>=sizex-1 || neiy<=0 || neiy>=sizey-1){
      if( par.periodic_boundaries ){
        if(neix<=0) neix=sizex-2+neix;
        if(neix>=sizex-1) neix=neix-sizex+2;
        if(neiy<=0) neiy=sizey-2+neiy;
        if(neiy>=sizey-1) neiy=neiy-sizey+2;
      }else{
        continue;
      }
    }


    //sigmaneigh = sigma[ neix ][ neiy ];



    //if(x+nx[k]>0 && x+nx[k]<sizex-1 && y+ny[k]>0 && y+ny[k]<sizey-1 ){
    //  point=sigma[x+nx[k]][y+ny[k]]; //take sigma of neighbour, this can be medium
    point=sigma[neix][neiy]; //take sigma of neighbour, this can be medium
    // if neigh is not the same as self (you don't have boundaries with yourself, are you thereby really free alone? :P )
      // if neigh is not same as self, there are a few cases:
      // self is not medium <- we update
      if(point!=idSelf){
        //if idSelf not medium, we will remove this point from its contact (because we are going to copy neigh in it)
        if(idSelf)
          check+=(*cell)[idSelf].updateNeighbourBoundary(point, -1);
        if(point)
          check+=(*cell)[point].updateNeighbourBoundary(idSelf, -1); //and if point is not medium, we are removing self from contacts of neigh
      }
      // if neigh is not the same as
      if(point!=idNeigh){
        if(idNeigh)
          check+=(*cell)[idNeigh].updateNeighbourBoundary(point, 1);
        if(point)
          check+=(*cell)[point].updateNeighbourBoundary(idNeigh, 1);
      }
      if(check){
        printf("error in ConvertSpinToMedium(): wrongly updating neighbours of copy event idSelf %d and idNeigh %d\n", idSelf, idNeigh);
        //printf("agent nr %d\n", agentid);
        exit(1);
      }
    //}
  }

  //cout<<"in ConvertSpin, after: "<<(*cell)[1].neighbours[0].first<<endl;


  sigma[x][y] = MEDIUM;    // THE SPIN COPYING
  //cerr<<"FYI: ConvertSpinToMedium() happened" << endl;
  //exit(1);

}



//this the function that changes the sigma
// it also does book-keeping of everything
void CellularPotts::ConvertSpin(int x,int y,int xp,int yp)
{
  int tmpcell;
  if( (tmpcell=sigma[x][y]) ) { // if tmpcell is not MEDIUM
    (*cell)[tmpcell].DecrementArea(); // returns --area
    (*cell)[tmpcell].RemoveSiteFromMoments(x,y);

    //if cell's area is zero -> kill it
    if ( !(*cell)[tmpcell].Area() ){
      (*cell)[tmpcell].Apoptose();
      //cerr << "Cell " << tmpcell << " apoptosed\n";
    }
  }

  if( (tmpcell=sigma[xp][yp]) ) {// if tmpcell is not MEDIUM, it gains one pixel
    (*cell)[tmpcell].IncrementArea();
    (*cell)[tmpcell].AddSiteToMoments(x,y);

  }

  // Book-keeping of contacts
  int check=0;
  int point;
  int idSelf=sigma[x][y];
  int idNeigh=sigma[xp][yp];
  ///use this if you need info about contacts with neighbours
  for(int k=1; k<=n_nb; k++)
  {
    // Do you not update contacts at boundaries???
    int neix=x+nx[k];
    int neiy=y+ny[k];
    if(neix<=0 || neix>=sizex-1 || neiy<=0 || neiy>=sizey-1){
      if( par.periodic_boundaries ){
        if(neix<=0) neix=sizex-2+neix;
        if(neix>=sizex-1) neix=neix-sizex+2;
        if(neiy<=0) neiy=sizey-2+neiy;
        if(neiy>=sizey-1) neiy=neiy-sizey+2;
      }else{
        continue;
      }
    }


    //sigmaneigh = sigma[ neix ][ neiy ];



    //if(x+nx[k]>0 && x+nx[k]<sizex-1 && y+ny[k]>0 && y+ny[k]<sizey-1 ){
    //  point=sigma[x+nx[k]][y+ny[k]]; //take sigma of neighbour, this can be medium
    point=sigma[neix][neiy]; //take sigma of neighbour, this can be medium

    //Error goes that point = 658 = idNeigh, idSelf=MEDIUM
    if(point!=idSelf)
    {
      if(idSelf){
        check+=(*cell)[idSelf].updateNeighbourBoundary(point, -1);
        if(check) cerr<<"Here1"<<endl;
      }
      if(point){
        check+=(*cell)[point].updateNeighbourBoundary(idSelf, -1);
        if(check) cerr<<"Here2"<<endl;
      }
    }
    if(point!=idNeigh)
    {
      if(idNeigh)
        check+=(*cell)[idNeigh].updateNeighbourBoundary(point, 1);
      if(point)
        check+=(*cell)[point].updateNeighbourBoundary(idNeigh, 1);
    }
    if(check){
      printf("error in ConvertSpin(): wrongly updating neighbours of copy event idSelf %d and idNeigh %d\n", idSelf, idNeigh);
      cerr<<"Extra info: idSelf=sigma[x][y], where x and y are "<<x<<","<<y<<endl;
      cerr<<"Extra info: idNeigh=sigma[xp][yp], where xp and yp are "<<xp<<","<<yp<<endl;
      cerr<<"Extra info: neix,neiy "<<neix<<","<<neiy<<endl;
      //printf("agent nr %d\n", agentid);

      exit(1);
    }
   //}
  }

  //cout<<"in ConvertSpin, after: "<<(*cell)[1].neighbours[0].first<<endl;

  sigma[x][y] = sigma[xp][yp]; // = idNeigh - THE SPIN COPYING


}


/** PUBLIC **/
// stiff = conn_diss i.e. the Delta H coming from break connectivity of cell, it is zero in debugging
int CellularPotts::CopyvProb(int DH,  double stiff) {

  double dd;
  int s;
  s=(int)stiff;
  if (DH<=-s) return 2; //accept copy if DH+s <= 0

  // if DH becomes extremely large, calculate probability on-the-fly
  // common values - from 0 to 1023 - are stored in copyprob - checked
  // BOLTZMANN = 1024 and is defined in sticky.h
  if(DH+s > BOLTZMANN-1)
    dd=exp( -( (double)(DH+s)/par.T ));
  else
    dd=copyprob[DH+s];

  //for(int i=0;i<1024;i++) cerr<<i<<" "<<copyprob[i]<<endl;
  //exit(1);

  if (RANDOM()<dd) return 1; else return 0;
}

void CellularPotts::CopyProb(double T) {
  int i;
  for ( i = 0; i < BOLTZMANN; i++ )
    copyprob[i] = exp( -( (double)(i)/T ) );
}

void CellularPotts::FreezeAmoebae(void)
{
  if (frozen)
    frozen=FALSE;
  else
    frozen=TRUE;
}

#include <fstream>
//! Monte Carlo Step. Returns summed energy change
// In this version there is a chance that you copy medium (from outer space)
int CellularPotts::AmoebaeMove2(PDE *PDEfield)
{
  double chanceofmedium=par.chancemediumcopied;  //this is the chance that we test for medium instead of actual pixel
                              // (i.e. that we raise mediumflag for that)
                              //this is NOT probability that medium is inserted,
                              // that follows the normal process and depends on nrg.
                              //BTW, we raise this after we are sure that k != kp: we don't want medium inside the cell :P

  bool mediumflag, converted, inout1, inout2;

  // cout<<"I'm in AmoebaeMove"<<endl;
  int loop,p;
  //int updated=0;
  thetime++;
  int SumDH=0;

  if (frozen)
    return 0;

  //variables used in loop
  int x,y, xp,yp, xn, yn;
  int k,kp;
  int edgesize, rindex;

  loop = edgeSetVector.size_map()/n_nb;
  //cout <<" loop size "<<loop<<endl;
  for (int i=0;i<loop;i++) {
    mediumflag=false;   // what happens if this is not re-set to false all the times?
                        // that the first time it goes true it will be true for the rest of the loop!
    converted=false;
    edgesize=edgeSetVector.size_map(); //are there still edges to update?
    if (edgesize==0){break;}
    else{
      rindex=RandomNumber(edgesize)-1; //pick random edge (Randomnumber goes from 1-N)
      auto it = edgeSetVector[rindex]; //this returns an iterator to an unordered map element out of the vector edgevector in
      x=(*it).first[0];
      y=(*it).first[1];
      xp=(*it).first[2];
      yp=(*it).first[3];
      k=sigma[x][y];

     if (par.periodic_boundaries) {
      // Boundary management:
      // since we are asynchronic, we cannot just copy the borders once
      // every MCS
      if (x<=0) x=sizex-2+x;
      if (y<=0) y=sizey-2+y;
      if (x>=sizex-1) x=x-sizex+2;
      if (y>=sizey-1) y=y-sizey+2;
      if(xp<=0) xp=sizex-2+xp;  //sizex-2 because both borders
      if(yp<=0) yp=sizey-2+yp;
      if(xp>=sizex-1) xp=xp-sizex+2;
      if(yp>=sizey-1) yp=yp-sizey+2;
      k=sigma[x][y];
      kp=sigma[xp][yp];

    } else {

      if (xp<=0 || yp<=0 || xp>=sizex-1 || yp>=sizey-1)
        kp=-1;
      else
        kp=sigma[xp][yp];
    }
    // test for border state (relevant only if we do not use periodic boundaries)
    // test always passed with periodic boundaries
    if (kp!=-1) {
      // Don't even think of copying the special border state into you!
      if ( k  != kp ) {
        // connectivity dissipation:
        int H_diss=0;

        if (!ConnectivityPreservedP(x,y))
          H_diss=par.conn_diss;

        //if k is different from medium and with small chance = chancemedium,
        // we propose to copy medium into k
        if(k!= MEDIUM && RANDOM() < chanceofmedium){
          kp=MEDIUM; //and so kp = MEDIUM
          mediumflag=true;
        }
        //IF we are not copying medium from outer space
        if(!mediumflag){
          int D_H=DeltaH(x,y,xp,yp,PDEfield);

          if( (p=CopyvProb(D_H,H_diss))>0 ){
            ConvertSpin( x,y,xp,yp );
            converted=true;
            SumDH+=D_H;
          }
        }else{
          int D_H=DeltaHWithMedium(x,y,PDEfield);

          if ((p=CopyvProb(D_H,H_diss))>0) {
            ConvertSpinToMedium( x,y );
            converted=true;
            SumDH+=D_H;
          }

        }
        if (converted){
          //go through neighbourhood
          for (int j = 1; j <= n_nb; j++){
    				xn = nx[j]+x;
    				yn = ny[j]+y;

    				if (par.periodic_boundaries) {
    					 	// since we are asynchronic, we cannot just copy the borders once
    					 	// every MCS
    					if (xn<=0)
    						xn=sizex-2+xn;
    					if(yn<=0)
    						yn=sizey-2+yn;
    					if (xn>=sizex-1)
    						xn=xn-sizex+2;
    				  if (yn>=sizey-1)
    						yn=yn-sizey+2;
    				}
    				if (xn>0 && yn>0 && xn<sizex-1 && yn<sizey-1){//if the neighbour site is within the lattice

    					// if we should add the edge to the edgelist, add it
              if (sigma[xn][yn] != sigma[x][y] && sgn(sigma[xn][yn])+sgn(sigma[x][y])>=0){ //if we should add the edge to the edgelist, add it
                inout1=edgeSetVector.insert({x,y,xn,yn});
                inout2=edgeSetVector.insert({xn,yn,x,y});
                if (inout1) loop += (double)1./(double)n_nb; // (double)2/n_nb IS A DOUBLE
                if (inout2) loop += (double)1./(double)n_nb; // (double)2/n_nb IS A DOUBLE
    					}
    					// if the sites have the same celltype and they have an edge, remove it. not sure what the sgn does
              if ((sigma[xn][yn] == sigma[x][y] || sgn(sigma[xn][yn])+sgn(sigma[x][y])<0)){
                inout1=edgeSetVector.erase({x,y,xn,yn});
                inout2=edgeSetVector.erase({xn,yn,x,y});
                if (inout1) loop -= (double)1./(double)n_nb; // (double)2/n_nb IS A DOUBLE
                if (inout2) loop -= (double)1./(double)n_nb;

    					}
    				} //end if neighbour in lattice
    			}//end neighbourhood loop
        }//end if converted

      } //end if neighbour is different pixel type
    }//end if neighbour is in the field
  } //if edgelist still full
  }//end for loop

  return SumDH;

}

// // In this version there is a chance that you copy medium (from outer space)
// int CellularPotts::AmoebaeMove2(PDE *PDEfield)
// {
//   double chanceofmedium=par.chancemediumcopied;  //this is the chance that we test for medium instead of actual pixel
//                               // (i.e. that we raise mediumflag for that)
//                               //this is NOT probability that medium is inserted,
//                               // that follows the normal process and depends on nrg.
//                               //BTW, we raise this after we are sure that k != kp: we don't want medium inside the cell :P
//
//   bool mediumflag;
//
//   // cout<<"I'm in AmoebaeMove"<<endl;
//   int loop,p;
//   //int updated=0;
//   thetime++;
//   int SumDH=0;
//
//   if (frozen)
//     return 0;
//
//   loop=(sizex-2)*(sizey-2);
//
//   for (int i=0;i<loop;i++) {
//     mediumflag=false;   // what happens if this is not re-set to false all the times?
//                         // that the first time it goes true it will be true for the rest of the loop!
//     // take a random site
//     int xy = (int)(RANDOM()*(sizex-2)*(sizey-2));
//     int x = xy%(sizex-2)+1;
//     int y = xy/(sizex-2)+1;
//
//     // take a random neighbour
//     // effect of copying could be that medium is introduced
//
//     int xyp=(int)(n_nb*RANDOM()+1);
//     int xp = nx[xyp]+x;
//     int yp = ny[xyp]+y;
//
//     int k=sigma[x][y];  // k is the sigma of focal grid point
//
//     int kp;
//     if (par.periodic_boundaries) {
//       // Boundary management:
//       // since we are asynchronic, we cannot just copy the borders once
//       // every MCS
//
//       if(xp<=0) xp=sizex-2+xp;  //sizex-2 because both borders
//       if(yp<=0) yp=sizey-2+yp;
//       if(xp>=sizex-1) xp=xp-sizex+2;
//       if(yp>=sizey-1) yp=yp-sizey+2;
//
//       kp=sigma[xp][yp]; // kp is the sigma of neighbour grid point
//
//     } else {
//
//       if (xp<=0 || yp<=0 || xp>=sizex-1 || yp>=sizey-1)
//         kp=-1;
//       else
//         kp=sigma[xp][yp];
//     }
//
//     // test for border state (relevant only if we do not use periodic boundaries)
//     // test always passed with periodic boundaries
//     if (kp!=-1) {
//       // Don't even think of copying the special border state into you!
//
//       // Try to copy if sites do not belong to the same cell
//       // if k == kp it would be pointless to copy
//       if ( k  != kp ) {
//         // connectivity dissipation:
//         int H_diss=0;
//
//         if (!ConnectivityPreservedP(x,y))
//           H_diss=par.conn_diss;
//
//         //Here we already know that k != kp,
//         //if k is different from medium and with small chance = chancemedium,
//         // we propose to copy medium into k
//         if(k!= MEDIUM && RANDOM() < chanceofmedium){
//           //cerr<<"AmoabeMove2, error during debugging we shouldn't be here"<<endl;
//           kp=MEDIUM; //and so kp = MEDIUM
//           mediumflag=true;
//         }
//         //IF we are not copying medium from outer space
//         if(!mediumflag){
//
// //           if(sigma[x][y]==sigma[xp][yp]){
// //             cerr<<"AmoebaeMove2: Error. Same sigmas"<<endl;
// //             exit(1);
// //           }
//
//           int D_H=DeltaH(x,y,xp,yp,PDEfield);
//
//           //cerr<<"Hello AmoebaeMove2: D_H = "<<D_H<<endl;
//           if( (p=CopyvProb(D_H,H_diss))>0 ){
//             //cerr<<"Hello before ConvertSpin"<<endl;
//             ConvertSpin( x,y,xp,yp );
//
//             //cerr<<"Hello after ConvertSpin"<<endl;
//             SumDH+=D_H;
//           }
//         }else{
// //           cerr<<"Hello AmoebaeMove2: medium flag activated"<<endl;
//           int D_H=DeltaHWithMedium(x,y,PDEfield);
//
// //           cerr<<"Hello AmoebaeMove2 0.4"<<endl;
//           if ((p=CopyvProb(D_H,H_diss))>0) {
//             //cerr<<"this shouldn't happen if chanceofmedium =0 "<<endl;
//             //exit(1);
//             //cerr<<"Hello before ConvertSpinToMedium"<<endl;
//             ConvertSpinToMedium( x,y );
//             //cerr<<"Hello after ConvertSpinToMedium"<<endl;
//             SumDH+=D_H;
//           }
//
//         }
//       }
//     }
//   }
//
//   // cout<<"Exiting AmoebaeMove"<<endl;
//   return SumDH;
//
// }


//! Monte Carlo Step. Returns summed energy change
int CellularPotts::AmoebaeMove(PDE *PDEfield)
{
 // cout<<"I'm in AmoebaeMove"<<endl;
  int loop,p;
  //int updated=0;
  thetime++;
  int SumDH=0;

  if (frozen)
    return 0;

  loop=(sizex-2)*(sizey-2);

  for (int i=0;i<loop;i++) {

    // take a random site
    int xy = (int)(RANDOM()*(sizex-2)*(sizey-2));
    int x = xy%(sizex-2)+1;
    int y = xy/(sizex-2)+1;

    // take a random neighbour


    int xyp=(int)(n_nb*RANDOM()+1);
    int xp = nx[xyp]+x;
    int yp = ny[xyp]+y;

    int k=sigma[x][y];  // k is the sigma of focal grid point

    int kp;
    if (par.periodic_boundaries) {
      // Boundary management:
      // since we are asynchronic, we cannot just copy the borders once
      // every MCS

      if(xp<=0) xp=sizex-2+xp;
      if(yp<=0) yp=sizey-2+yp;
      if(xp>=sizex-1) xp=xp-sizex+2;
      if(yp>=sizey-1) yp=yp-sizey+2;

      kp=sigma[xp][yp]; // kp is the sigma of neighbour grid point

    } else {

      if (xp<=0 || yp<=0 || xp>=sizex-1 || yp>=sizey-1)
        kp=-1;
      else
        kp=sigma[xp][yp];

    }


    // test for border state (relevant only if we do not use periodic boundaries)
    // test always passed with periodic boundaries
    if (kp!=-1) {
      // Don't even think of copying the special border state into you!

      // Try to copy if sites do not belong to the same cell
      // if k == kp it would be pointless to copy
      if ( k  != kp ) {
        // connectivity dissipation:
        int H_diss=0;
        if (!ConnectivityPreservedP(x,y))
          H_diss=par.conn_diss;
        int D_H=DeltaH(x,y,xp,yp,PDEfield);
        if ((p=CopyvProb(D_H,H_diss))>0) {
          ConvertSpin( x,y,xp,yp );
          SumDH+=D_H;
        }
      }
    }
  }

  // cout<<"Exiting AmoebaeMove"<<endl;
  return SumDH;

}

/** A simple method to plot all sigma's in window
    without the black lines */
void CellularPotts::PlotSigma(Graphics *g, int mag) {

  for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++) {
      for (int xm=0;xm<mag;xm++)
	for (int ym=0;ym<mag;ym++)
      g->Point( sigma[x][y], mag*x+xm, mag*y+ym);
  }

}

//function to plot the direction of the target vector as a cell colour
void CellularPotts::CellAngleColour(Graphics *g)
{

  for ( int i = 0; i < sizex; i++ )
    for ( int j = 0; j < sizey; j++ ) {
      int colour;

      if(i==0 || i== sizex-1 || j==0 || j == sizey){
        colour=0;
        g->Point( colour, 2*i, 2*j);
        g->Point( colour, 2*i+1, 2*j);
        g->Point( colour, 2*i, 2*j+1);
        g->Point( colour, 2*i+1, 2*j+1);
        continue;
      }

      if (sigma[i][j]<=0) {
        colour=0;
      }else{
        colour = (*cell)[sigma[i][j]].AngleColour();
        //colour = sigma[i][j];
      }

      //colour point if this is a cell
      if (sigma[i][j]>0){

        g->Point( colour, 2*i, 2*j); //draws 2i,2j
        //check if the other 3 pixels in the image should be coloured as boundary
        //if this cell different from what is on i+1,j
        //south
        if ( sigma[i][j] != sigma[i+1][j] && i+1 < sizex-1){
          g->Point( 1, 2*i+1, 2*j );
        }
        else{
          g->Point( colour, 2*i+1, 2*j );
        }
        //east
        if( sigma[i][j] != sigma[i][j+1] && j+1<sizey-1){
          g->Point( 1, 2*i, 2*j+1 );
        }
        else {
          g->Point( colour, 2*i, 2*j+1 );
        }
        //southeast
        if (i+1<sizex-1 && j+1<sizey-1 && (sigma[i][j]!=sigma[i+1][j+1] || sigma[i+1][j]!=sigma[i][j+1]) ){
          g->Point( 1, 2*i+1, 2*j+1 );
        }
        else {
          g->Point( colour, 2*i+1, 2*j+1 );
        }
      }//if this is a cell

    } //end of for loop




}

void CellularPotts::CellOrderColour(Graphics *g)
{
  vector <double> cellorder;
  double temp;

//determine order value for each cell (mapped to value between 0 and 1, 0 being fully aligned and 1 being
//fully counteraligned)
  for(auto c: (*cell)){
    temp=0.;
    for (auto n: c.neighbours){
      temp+=(c.tvecx*(*cell)[n.first].tvecx+c.tvecy*(*cell)[n.first].tvecy)+1.;
    }
    cellorder.push_back(temp/(2*c.neighbours.size()));
  }

  for ( int i = 0; i < sizex; i++ )
    for ( int j = 0; j < sizey; j++ ) {
      int colour;

      if(i==0 || i== sizex-1 || j==0 || j == sizey){
        colour=0;
        g->Point( colour, 2*i, 2*j);
        g->Point( colour, 2*i+1, 2*j);
        g->Point( colour, 2*i, 2*j+1);
        g->Point( colour, 2*i+1, 2*j+1);
        continue;
      }

      if (sigma[i][j]<=0) {
        colour=0;
      }else{
        colour = (int)(cellorder[sigma[i][j]]*127+2);
        //colour = sigma[i][j];
      }

      //colour point if this is a cell
      if (sigma[i][j]>0){

        g->Point( colour, 2*i, 2*j); //draws 2i,2j
        //check if the other 3 pixels in the image should be coloured as boundary
        //if this cell different from what is on i+1,j
        //south
        if ( sigma[i][j] != sigma[i+1][j] && i+1 < sizex-1){
          g->Point( 1, 2*i+1, 2*j );
        }
        else{
          g->Point( colour, 2*i+1, 2*j );
        }
        //east
        if( sigma[i][j] != sigma[i][j+1] && j+1<sizey-1){
          g->Point( 1, 2*i, 2*j+1 );
        }
        else {
          g->Point( colour, 2*i, 2*j+1 );
        }
        //southeast
        if (i+1<sizex-1 && j+1<sizey-1 && (sigma[i][j]!=sigma[i+1][j+1] || sigma[i+1][j]!=sigma[i][j+1]) ){
          g->Point( 1, 2*i+1, 2*j+1 );
        }
        else {
          g->Point( colour, 2*i+1, 2*j+1 );
        }
      }//if this is a cell

    } //end of for loop




}

int **CellularPotts::SearchNandPlot(Graphics *g, bool get_neighbours)
{
  int i, j,q;
  int **neighbours=0;


  // Allocate neighbour 2D matrix - of size ncells x ncells
  // Notice that this does not seem to be done because the function is called with
  // get_neighbours=false
  if (get_neighbours) {
    neighbours=(int **)malloc((cell->size()+1)*sizeof(int *));
    if (neighbours==NULL)
      MemoryWarning();

    neighbours[0]=(int *)malloc((cell->size()+1)*(cell->size()+1)*sizeof(int));
    if (neighbours[0]==NULL)
      MemoryWarning();

    for (i=1;i<(int)cell->size()+1;i++)
      neighbours[i]=neighbours[i-1]+(cell->size()+1);

    /* Clear this matrix */
    for (i=0;i<((int)cell->size()+1)*((int)cell->size()+1);i++)
      neighbours[0][i]=EMPTY;
  }

  //for ( i = 0; i < sizex-1; i++ )
  //  for ( j = 0; j < sizey-1; j++ ) {
  for ( i = 0; i < sizex; i++ )
    for ( j = 0; j < sizey; j++ ) {
      int colour;

      if(i==0 || i== sizex-1 || j==0 || j == sizey){
        colour=0;
        g->Point( colour, 2*i, 2*j);
        g->Point( colour, 2*i+1, 2*j);
        g->Point( colour, 2*i, 2*j+1);
        g->Point( colour, 2*i+1, 2*j+1);
        continue;
      }

      if (sigma[i][j]<=0) {
        colour=0;
      }else{
        if(!par.divisioncolour){
          colour = (*cell)[sigma[i][j]].Colour()+7*((*cell)[sigma[i][j]].getTau()-1);
          // if (colour !=2 && colour !=3) cerr<<" Wrong colour! "<< colour<<endl;
        }else{
          colour = (*cell)[sigma[i][j]].getTau()+1+( (*cell)[sigma[i][j]].getTau()-1 ) *4+(*cell)[sigma[i][j]].TimesDivided();
        }


        //colour = sigma[i][j];
      }

      if (g && sigma[i][j]>0)  /* if draw */
        g->Point( colour, 2*i, 2*j); //draws 2i,2j

      // WAS -> if ( sigma[i][j] != sigma[i+1][j] )  /* if cellborder */ /* etc. etc. */
      if ( sigma[i][j] != sigma[i+1][j] && i+1 < sizex-1) //if this cell different from what is on i+1,j
      {
        if(g) g->Point( 1, 2*i+1, 2*j ); //draws 2i+1,2j black because cells neighbours
        if (get_neighbours){
          if (sigma[i][j]>0) {
            for (q=0;q<(int)cell->size();q++)
              if (neighbours[sigma[i][j]][q]==EMPTY) {
                neighbours[sigma[i][j]][q]=sigma[i+1][j];
                break;
              }else if(neighbours[sigma[i][j]][q]==sigma[i+1][j])
                break;
	    }
	    if (sigma[i+1][j]>0) {
	      for (q=0;q<(int)cell->size();q++)
            if (neighbours[sigma[i+1][j]][q]==EMPTY) {
              neighbours[sigma[i+1][j]][q]=sigma[i][j];
              break;
            }else if(neighbours[sigma[i+1][j]][q]==sigma[i][j])
              break;
	    }
	  } //endif get_neighbours
	// else same sigma
	}else if (g && sigma[i][j]>0)
      g->Point( colour, 2*i+1, 2*j ); //we colour the point i+1,j
      //we check sigma at point i,j+1
      if( sigma[i][j] != sigma[i][j+1] && j+1<sizey-1){
        if(g) g->Point( 1, 2*i, 2*j+1 );
        //check neighbours
        if (get_neighbours){
          if (sigma[i][j]>0) {
            for (q=0;q<(int)cell->size();q++)
              if (neighbours[sigma[i][j]][q]==EMPTY) {
            neighbours[sigma[i][j]][q]=sigma[i][j+1];
            break;
              }
              else
            if (neighbours[sigma[i][j]][q]==sigma[i][j+1])
              break;
          }

          if (sigma[i][j+1]>0) {

            for (q=0;q<(int)cell->size();q++)
              if (neighbours[sigma[i][j+1]][q]==EMPTY) {
            neighbours[sigma[i][j+1]][q]=sigma[i][j];
            break;
              }
              else
            if (neighbours[sigma[i][j+1]][q]==sigma[i][j])
              break;
          }
        }//end get_neighbours
      }
      else if (g && sigma[i][j]>0)
        g->Point( colour, 2*i, 2*j+1 );

        /* Cells that touch eachother's corners are NO neighbours */
        if (i+1<sizex-1 && j+1<sizey-1 && (sigma[i][j]!=sigma[i+1][j+1] || sigma[i+1][j]!=sigma[i][j+1]) ){
          if(g) g->Point( 1, 2*i+1, 2*j+1 );
        }else if(g && sigma[i][j]>0)
          g->Point( colour, 2*i+1, 2*j+1 );
    } //end of for loop

  if (get_neighbours)
    return neighbours;
  else
    return 0;

}


// void CellularPotts::ReadZygotePicture(void) {
//
//
//
//   int pix,cells,i,j,c,p,checkx,checky;
//   char **pixelmap;
//   char pixel[3];
//
//   sscanf(ZYGXPM(ZYGOTE)[0],"%d %d %d %d",&checkx,&checky,&cells,&pix);
//
//   if ((checkx>sizex)||(checky>sizey)) {
//     std::cerr <<  "ReadZygote: The included xpm picture is smaller than the grid!\n";
//     std::cerr << "\n Please adjust either the grid size or the picture size.\n";
//     std::cerr << sizex << "," << sizey << "," << checkx << "," << checky << "\n";
//     exit(1);
//   }
//
//   pixelmap=(char **)malloc(cells*sizeof(char *));
//   if (pixelmap==NULL) MemoryWarning();
//
//   pixelmap[0]=(char *)malloc(cells*3*sizeof(char));
//   if (pixelmap[0]==NULL) MemoryWarning();
//
//   for(i=1;i<cells;i++)
//     pixelmap[i]=pixelmap[i-1]+3;
//
//   for (i=0;i<cells;i++) {
//     for (j=0;j<pix;j++)
//       pixelmap[i][j]=ZYGXPM(ZYGOTE)[i+1][j];
//     pixelmap[i][pix]='\0';
//   }
//
//   for (i=0;i<sizex*sizey;i++) sigma[0][i]=0;
//   fprintf(stderr,"[%d %d]\n",checkx,checky);
//
//   int offs_x, offs_y;
//   offs_x=(sizex-checkx)/2;
//   offs_y=(sizey-checky)/2;
//
//   for (i=0;i<checkx;i++)
//     for (j=0;j<checky;j++) {
//       for (p=0;p<pix;p++)
//         pixel[p]=ZYGXPM(ZYGOTE)[cells+1+j][i*pix+p];
//
//       pixel[pix]='\0';
//
//       for (c=0;c<cells;c++) {
// 	if (!(strcmp(pixelmap[c],pixel))) {
// 	  if ( (sigma[offs_x+i][offs_y+j]=c) ) {
//
// 	    // if c is _NOT_ medium (then c=0)
// 	    // assign pixel values from "sigmamax"
// 	    sigma[offs_x+i][offs_y+j]+=(Cell::MaxSigma()-1);
// 	  }
// 	}
//
//       }
//     }
//
//   free(pixelmap[0]);
//   free(pixelmap);
// }


void CellularPotts::ConstructInitCells(Dish &beast){

  // Get the maximum cell ID (mostly equal to the cell number)
  //int loop=sizex*sizey;
  int cells=0;
  for (int i=1;i<sizex-1;i++) for(int j=1;j<sizey-1;j++){
    if(cells<sigma[i][j]) cells=sigma[i][j];
  }
  cerr<<"nr cells placed "<<cells<<endl;

  cerr << "[ cells = " << cells << "]\n";

  // construct enough cells for the zygote.  "cells", contains the
  // number of colours (excluding background).
  {
    for (int i=0; i<cells; i++) {
      cell->push_back(Cell(beast));
    }
  }

  //initialises mean x and y to some values, does not matter what,
  // as long as meanx and y get a defined number BEFORE we calculate MeasureCellSize()
  // which, in turn, if periodic_boundaries = true, depends on meanx and y
  for(vector<Cell>::iterator c=cell->begin(); c!=cell->end();c++){
    c->InitMeanX(sizex/2.);
    c->InitMeanY(sizey/2.);
  }


  // Set the area and target area of the cell
  // makes use of the pointer to the Cell pointer of Dish
  // which is a member of CellularPotts
  MeasureCellSizes();

  //for (vector<Cell>::iterator c=cell->begin(); c!=cell->end();c++) {
  //  cerr<<"sigma: "<<c->Sigma()<<". New meanx and y calculated as: "<<c->getXpos()<<" "<<c->getYpos()<<endl;
  //}

  // set zygote_area to mean cell area.
  int mean_area=0;
  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    mean_area+=c->Area();
  }
  if (cells!=0)
    mean_area/=cells;

  zygote_area=mean_area;

  cerr << "mean_area = " << mean_area << "\n";
  // set all cell areas to the mean area
  {
    for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
      c->setMu(0.0);
      c->SetEatProb(par.eatprob);

      if (par.target_area) {
	c->SetTargetArea(par.target_area);
      } else {
	c->SetTargetArea(mean_area);

      }
    }
  }
  cerr << "ConstructInitCells is done"<<endl;
}

void CellularPotts::MeasureCellSizes(void) {

  // Clean areas of all cells, including medium
  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    c->SetTargetArea(0);
    c->area = 0;
    c->CleanMoments();
  }

  // calculate the area of the cells
  for (int x=1;x<sizex-1;x++) {
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	(*cell)[sigma[x][y]].IncrementTargetArea();
	(*cell)[sigma[x][y]].IncrementArea();
	(*cell)[sigma[x][y]].AddSiteToMoments(x,y);

      }
    }
  }

  // set the actual area to the target area
  {
  for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
    c->SetAreaToTarget();
  }
  }
}

void CellularPotts::MeasureCellSize(Cell &c) {

  c.CleanMoments();

  // calculate the area of the cell
  for (int x=1;x<sizex-1;x++) {
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y] == c.sigma) {
	(*cell)[sigma[x][y]].IncrementTargetArea();
	(*cell)[sigma[x][y]].IncrementArea();
	(*cell)[sigma[x][y]].AddSiteToMoments(x,y);

      }
    }
  }

//   // set the actual area to the target area
//   {
//   for (vector<Cell>::iterator c=cell->begin();c!=cell->end();c++) {
//     c->SetAreaToTarget();

//   }

}

// Rewritten version of the funciton below
// it uses the cell meanx and meany for the center of mass of the cell,
// IN THE FUTURE -> only calculates PCA for cells that need to be split...
Dir *CellularPotts::FindCellDirections3(void) const
{
  // double array - allocated (why c style?) the size population
  double *sumx=0,*sumy=0;
  double *sumxx=0,*sumxy=0,*sumyy=0;
  double *n=0;

  double xmean=0,ymean=0,sxx=0,sxy=0,syy=0;
  double D,lb1=0,lb2=0;

  Dir *celldir;


  /* Allocation of sufficient memory space */
  if( (sumx= (double *)malloc(cell->size()*sizeof(double)))==NULL)
    MemoryWarning();
  if( (sumy= (double *)malloc(cell->size()*sizeof(double)))==NULL)
    MemoryWarning();
  if ((sumxx=(double *)malloc(cell->size()*sizeof(double)))==NULL)
    MemoryWarning();
  if((sumxy=(double *)malloc(cell->size()*sizeof(double)))==NULL)
    MemoryWarning();
  if((sumyy=(double *)malloc(cell->size()*sizeof(double)))==NULL)
    MemoryWarning();

  if((n=(double *)malloc(cell->size()*sizeof(double)))==NULL)
    MemoryWarning();

  if ( !(celldir=new Dir[cell->size()]) )
    MemoryWarning();


  /* Initialization of the variables */
  for (int i=0;i<(int)cell->size();i++) {
    sumx[i]=0.;
    sumy[i]=0.;
    sumxx[i]=0.;
    sumxy[i]=0.;
    sumyy[i]=0.;
    n[i]=0L;
  }

  /* Find sumx, sumy, sumxx and sumxy for all cells */
  for (int x=1;x<sizex-1;x++){
    for (int y=1;y<sizey-1;y++){
      //for (int x=0;x<sizex;x++)
      //  for (int y=0;y<sizey;y++)
      if(sigma[x][y]>0){
        //cerr<<sigma[x][y]<<endl;

        sumx[0]+=(double)x;  //pos 0 contains global statistics
        sumy[0]+=(double)y;
        sumxx[0]+=(double)x*x;
        sumxy[0]+=(double)x*y;
        sumyy[0]+=(double)y*y;
        n[0]++;

        int isigma = sigma[x][y];

        double tmpx=x; //because we may change them if wrapped
        double tmpy=y;

        double meanx = (*cell)[isigma].meanx;
        double meany = (*cell)[isigma].meany;

        if(par.periodic_boundaries){
          //Now the wrapping - in this function we wrap around meanx:
          //if x is closer to running average when wrapped, we wrap it
          if( x-meanx>0 && x-meanx >(meanx+(sizex-2)-x ) ) {
            tmpx-=(sizex-2);
            //cerr<<"directions passb1"<<endl;
          }else if( meanx-x>0 && meanx-x > (x+(sizex-2)-meanx) ){
            tmpx+=(sizex-2);
            //cerr<<"directions passb1"<<endl;
          }
          //same for y
          if( y-meany >0 && y-meany > (meany+(sizey-2) -y) ){
            tmpy-=(sizey-2);
            //cerr<<"directions passb1"<<endl;
          }else if( meany - y>0 && meany-y>(y+(sizey-2)-meany) ){
            tmpy+=(sizey-2);
            //cerr<<"directions passb1"<<endl;
          }
        }

        sumx[isigma]+=tmpx; //pos indicised by sigma - local stats about x direction
        sumy[isigma]+=tmpy;

        //Different model here, we just do <(x-meanx)^2>
        sumxx[isigma]+=( pow( tmpx - meanx ,2.) );
        sumyy[isigma]+=( pow( tmpy - meany ,2.) );
        sumxy[isigma]+=( (tmpx - meanx)*(tmpy - meany) );

        // This is how it was
        //sumxx[isigma]+=(double)tmpx*tmpx;
        //sumxy[isigma]+=(double)tmpx*tmpy; // mixed coordinate
        //sumyy[isigma]+=(double)tmpy*tmpy;

        n[isigma]++;

      }
    }
  }
  //exit(1);
  // Compute the principal axes for all cells from eigenvalues of correlation matrix C
  //     ( sxx sxy )
  // C = <          >
  //     ( sxy syy )
  // We diagonalise this and find eigenvalues lb1 and lb2
  //recalculate the means while we're at it


  double small_number = 0.0000001;
  double large_enough_number = 1./small_number;

  for(int i=0;i<(int)cell->size();i++){
    celldir[i].meanx=(*cell)[i].meanx;
    celldir[i].meany=(*cell)[i].meany;
    //cerr<<"celldir[i].bb1 = "<<celldir[i].bb1<<endl;  // there are some problems here and at initilisation (see video)
    double xmean = celldir[i].meanx;
    double ymean = celldir[i].meany;
    sxx=0.;
    syy=0.;
    sxy=0.;
    // why this? of course moments are ill defined for small cells
    // maybe it's no problem because they don't divide

    if(n[i]>10){
      sxx=sumxx[i]/((double)n[i]); //or maybe n[i]-1 ? Would be strange because this is all the data
      syy=sumyy[i]/((double)n[i]);
      sxy=sumxy[i]/((double)n[i]);

      // Now we have the elements of the matrix, we diagonalise by solving eigenvalue problem (call x the eigenvalues)
      // x^2 - (sxx+syy)x + (sxx*syy-sxy*sxy) = 0
      D=sqrt( (sxx+syy)*(sxx+syy)-4.*(sxx*syy-sxy*sxy) ); // this is the discriminant
      lb1=(sxx+syy+D)/2.;    // these are the two solutions, i.e. the two eigenvalues
      lb2=(sxx+syy-D)/2.;
      // Now lb1 > lb2 because ... well, look at it, then lb1 is largest eigenvalue,
      // so its eigenvector is the principal axis of the cell
      celldir[i].lb1=lb1; celldir[i].lb2=lb2;
    }
    //Now we get the eigenvectors:
    // first eigenvector v1 is the solution of C - lb1*v1 =0
    // C is covar. matrix, lb1 is eigenv 1 just calculated.
    // v1 has an x and a y component, here ve express v1 as x/y = sxy/(lb1-sxx)

    //this case is when there is no covariance, so the cartesian axis are already the best basis
    if(sxy!=0.0){
//       cerr<<"sxy!=0, lb1: "<<lb1<<", syy: "<<syy<<", sxx: "<<sxx<<", sxy: "<<sxy<<endl;
      celldir[i].bb1=sxy/(lb1-syy); // APPARENTLY THIS IS FINE- This is y/x - shoudn't it be = syy/(sxy-lb1)?
      // I think this is correct : celldir[i].bb1=syy/(lb1-sxy);
      // and NOT WHAT IS THERE.. UNLESS I am wrong :P
      if (fabs(celldir[i].bb1)<small_number) {
        if (celldir[i].bb1>0.)
          celldir[i].bb1=small_number;
        else
          celldir[i].bb1=-small_number;
      }

      celldir[i].aa1=ymean-xmean*celldir[i].bb1; //this is the intercept to the y axis of the line with slope first eigenv.
      // which passes through xmean and ymean
      celldir[i].bb2= (-1.)/celldir[i].bb1; // bb2 is the direction perpendicular to the first eigenvector
      // (because the perpend. to a line y=mx+q has slope -1/m)

      celldir[i].aa2=ymean-celldir[i].bb2*xmean; // this is the intercept to y axis of the line of the second eigenvector
      // along this line we cut the cell !!!
    }else{
//       cerr<<"sxy=0, ";
      // USED TO BE
      //celldir[i].bb1=1.; WHICH IS DEFINITELY WRONG

      //Because later we are doing operations on the slope, we should choose a large enough number that does not overflow
      // a good idea could be to choose a slope so that the almost vertical line makes less than epsilon=0.1 error across the whole field
      // so that division is effectively vertical (or horizontal)
      // this is a line that has slope par.sizex/epsilon, let's add a little bit to be extra safe (times 2)
      // with a field size of 1000 and epsilon = 0.1 -> m= 10000 that's ok small for more calculations

      //double large_enough_number = (2.*(double)par.sizex)/0.1;
      //double large_enough_number = 1./0.0000001;
      double random_plus_or_minus_1 = -1+2*(int)(2.*RANDOM());
      if(sxx>syy){
        celldir[i].bb1 = 0.;
        celldir[i].aa1=ymean;
        celldir[i].bb2 = random_plus_or_minus_1*large_enough_number;
        celldir[i].aa2=ymean-celldir[i].bb2*xmean;
//         cerr<<"sxx>syy"<<endl;
      }
      else if(syy>sxx){
        celldir[i].bb1 = random_plus_or_minus_1*large_enough_number;
        celldir[i].aa1=ymean-xmean*celldir[i].bb1;
        //celldir[i].bb2 = 0.;
        celldir[i].bb2 = -1.*random_plus_or_minus_1*small_number;
        celldir[i].aa2=ymean;
//         cerr<<"syy>sxx"<<endl;
      }
      else{
//         cerr<<"syy=sxx"<<endl;
        celldir[i].bb1 = (RANDOM() <0.5)? 0. : (random_plus_or_minus_1*large_enough_number); //if sxx==syy we randomise vertical or horizontal
        if(celldir[i].bb1 > 1. || celldir[i].bb1 < 1.){
          celldir[i].aa1=ymean-xmean*celldir[i].bb1;
          //celldir[i].bb2 = 0.;
          celldir[i].bb2 = -1.*random_plus_or_minus_1*small_number;
          celldir[i].aa2=ymean;
        }else{
          celldir[i].aa1=ymean;
          celldir[i].bb2 = random_plus_or_minus_1*large_enough_number;
          celldir[i].aa2=ymean-celldir[i].bb2*xmean;
        }
      }
    }
//     cerr<<"Sigma: "<<i<<", bb2: "<<celldir[i].bb2<<", aa2: "<<celldir[i].aa2<<endl;
  }

  //}

  /* bevrijd gealloceerd geheugen */
  free(sumx);
  free(sumy);
  free(sumxx);
  free(sumxy);
  free(sumyy);
  free(n);

  return celldir;

}



//this version of the function patches up the one below for correct moment calculation with wrapped boundaries:
// we keep track of the running average x and y, if a pixel's distance to this mean is larger than its wrapped distance,
// we wrap it, and use that for calculating the moments:
// if at the end mean(x) <0 or > nrow-1, we wrap it -> this should give the correct meanx and y
// moreover, by wrapping we should ge the right variance, hope this works also for the sumxy :P
// this should work all the times - except maybe for some limit cases when a cell is as large as the field
Dir *CellularPotts::FindCellDirections2(void) const
{
  // double array - allocated (why c style?) the size population
  double *sumx=0,*sumy=0;
  double *sumxx=0,*sumxy=0,*sumyy=0;
  double *n=0;

  double xmean=0,ymean=0,sxx=0,sxy=0,syy=0;
  double D,lb1=0,lb2=0;

  Dir *celldir;

  /* Allocation of sufficient memory space */
  if( (sumx= (double *)malloc(cell->size()*sizeof(double)))==NULL)
    MemoryWarning();
  else
    if( (sumy= (double *)malloc(cell->size()*sizeof(double)))==NULL)
      MemoryWarning();
    else
      if ((sumxx=(double *)malloc(cell->size()*sizeof(double)))==NULL)
        MemoryWarning();
      else
        if((sumxy=(double *)malloc(cell->size()*sizeof(double)))==NULL)
          MemoryWarning();
        else
          if((sumyy=(double *)malloc(cell->size()*sizeof(double)))==NULL)
            MemoryWarning();
          else
            if((n=(double *)malloc(cell->size()*sizeof(double)))==NULL)
              MemoryWarning();


   if ( !(celldir=new Dir[cell->size()]) )
              MemoryWarning();


   /* Initialization of the variables */

   for (int i=0;i<(int)cell->size();i++) {

     sumx[i]=0.;
     sumy[i]=0.;
     sumxx[i]=0.;
     sumxy[i]=0.;
     sumyy[i]=0.;
     n[i]=0L;

   }


   /* Find sumx, sumy, sumxx and sumxy for all cells */
   for (int x=1;x<sizex-1;x++)
     for (int y=1;y<sizey-1;y++)
   //for (int x=0;x<sizex;x++)
   //  for (int y=0;y<sizey;y++)
       if (sigma[x][y]>0) {
         //cerr<<sigma[x][y]<<endl;


         sumx[0]+=(double)x;  //pos 0 contains global statistics
         sumy[0]+=(double)y;
         sumxx[0]+=(double)(x*x);
         sumxy[0]+=(double)(x*y);
         sumyy[0]+=(double)(y*y);
         n[0]++;

         int isigma = sigma[x][y];
         int tmpx=x; //because we may change them if wrapped
         int tmpy=y;

         if(par.periodic_boundaries){
           //Now the wrapping:
           if(n[isigma]!=0){
              double curr_avrg_x = sumx[isigma]/((double)(n[isigma]));
              double curr_avrg_y = sumy[isigma]/((double)(n[isigma]));

              //if x is closer to running average when wrapped, we wrap it
              if( x-curr_avrg_x>0 && x-curr_avrg_x >(curr_avrg_x+(sizex-2)-x ) ) {
                tmpx-=(sizex-2); //NO, sizex-2 -> WRONG: sizex -1 because only 1 border
                //cerr<<"directions passb1"<<endl;
              }
              else if( curr_avrg_x-x>0 && curr_avrg_x-x > (x+(sizex-2)-curr_avrg_x) ){
                tmpx+=(sizex-2);
                //cerr<<"directions passb1"<<endl;
              }
              //same for y
              if( y-curr_avrg_y >0 && y-curr_avrg_y > (curr_avrg_y+(sizey-2) -y) ){
                tmpy-=(sizey-2);
                //cerr<<"directions passb1"<<endl;
              }
              else if( curr_avrg_y - y>0 && curr_avrg_y-y>(y+(sizey-2)-curr_avrg_y) ){
                tmpy+=(sizey-2);
                //cerr<<"directions passb1"<<endl;
              }
           }
         }

         sumx[isigma]+=(double)tmpx; //pos indicised by sigma - local stats about x direction
         sumy[isigma]+=(double)tmpy;

         sumxx[isigma]+=(double)(tmpx*tmpx);
         sumxy[isigma]+=(double)(tmpx*tmpy); // mixed coordinate
         sumyy[isigma]+=(double)(tmpy*tmpy);

         n[isigma]++;

       }
       //exit(1);
       // Compute the principal axes for all cells from eigenvalues of correlation matrix C
       //     ( sxx sxy )
       // C = <          >
       //     ( sxy syy )
       // We diagonalise this and find eigenvalues lb1 and lb2

         for (int i=0;i<(int)cell->size();i++) {
           // why this? of course moments are ill defined for small cells
           // maybe it's no problem because they don't divide
           if (n[i]>10) {

             xmean=((double)sumx[i])/((double)n[i]); //mean x
             ymean=((double)sumy[i])/((double)n[i]); //mean y

             sxx=(double)(sumxx[i]) - ((double)(sumx[i]*sumx[i]))/(double)n[i]; //normalise <- This could be buggy across boundaries because sumxx might get weird across boundaries
             sxx=sxx/(double)(n[i]-1);  //standard deviation x (estimate = (E(X^2) - E(X)^2)/ (N*(N-1)) )
             //sxx=sxx/(double)(n[i]);

             sxy=(double)(sumxy[i]) - ((double)(sumx[i]*sumy[i]))/(double)n[i];
             sxy=sxy/(double)(n[i]-1);  //standard deviation xy
             //sxy=sxy/(double)(n[i]);

             syy=(double)(sumyy[i]) - ((double)(sumy[i]*sumy[i]))/(double)n[i];
             syy=syy/(double)(n[i]-1);  //standard deviation y
             //syy=syy/(double)(n[i]);

             //cerr<<"sumxx[i]: "<<sumxx[i]<<", sumyy[i]: "<<sumyy[i]<<", sumxy[i]: "<<sumxy[i]<<endl;
             //cerr<<"xmean: "<<xmean<<", ymean: "<<ymean<<", sxx: "<<sxx<<", syy: "<<syy<<", sxy: "<<sxy<<endl;

             // HANG ON HERE :
             // sumx contains as much info as meanx does, so we can wrap it no problem
             // notice that below here sumx is not used anymore, only meanx so we wrap only that
             if(par.periodic_boundaries){
               //it could be that with our wrapping meanx and y are now <0 or >nrow and ncol, we fix this here
               if(xmean<1) xmean+=sizex-2;
               else if(xmean>=sizex-1) xmean-=sizex-2;
               if(ymean<1) ymean+=sizey-2;
               else if(ymean>=sizey-1) ymean-=sizey-2;
             }
             celldir[i].meanx=xmean;
             celldir[i].meany=ymean;

             // Now we have the elements of the matrix, we diagonalise by solving eigenvalue problem (call x the eigenvalues)
             // x^2 - (sxx+syy)x + (sxx*syy-sxy*sxy) = 0
             D=sqrt( (sxx+syy)*(sxx+syy)-4.*(sxx*syy-sxy*sxy) ); // this is the discriminant
             lb1=(sxx+syy+D)/2.;    // these are the two solutions, i.e. the two eigenvalues
             lb2=(sxx+syy-D)/2.;
             // Now lb1 > lb2 because ... well, look at it, then lb1 is largest eigenvalue,
             // so its eigenvector is the principal axis of the cell
             celldir[i].lb1=lb1; celldir[i].lb2=lb2;
           }
           //Now we get the eigenvectors:
           // first eigenvector v1 is the solution of C - lb1*v1 =0
           // C is covar. matrix, lb1 is eigenv 1 just calculated.
           // v1 has an x and a y component, here ve express v1 as x/y = sxy/(lb1-sxx)

           //this case is when there is no covariance, so the cartesian axis are already the best basis
           if(sxy!=0.0){
             //cerr<<"sxy!=0, lb1: "<<lb1<<", syy: "<<syy<<", sxx: "<<sxx<<", sxy: "<<sxy<<endl;
             celldir[i].bb1=sxy/(lb1-syy); // APPARENTLY THIS IS FINE- This is y/x - shoudn't it be = syy/(sxy-lb1)?
                                           // I think this is correct : celldir[i].bb1=syy/(lb1-sxy);
                                           // and NOT WHAT IS THERE.. UNLESS I am wrong :P
             if (fabs(celldir[i].bb1)<.00001) {
               if (celldir[i].bb1>0.)
                 celldir[i].bb1=.00001;
               else
                 celldir[i].bb1=-.00001;
             }
             //cerr<<"celldir[i].bb1 = "<<celldir[i].bb1<<endl;  // there are some problems here and at initilisation (see video)

             celldir[i].aa1=ymean-xmean*celldir[i].bb1; //this is the intercept to the y axis of the line with slope first eigenv.
                                                      // which passes through xmean and ymean
             celldir[i].bb2= (-1.)/celldir[i].bb1; // bb2 is the direction perpendicular to the first eigenvector
                                                 // (because the perpend. to a line y=mx+q has slope -1/m)

             celldir[i].aa2=ymean-celldir[i].bb2*xmean; // this is the intercept to y axis of the line of the second eigenvector
                                                      // along this line we cut the cell !!!
           }else{
             // USED TO BE
             //celldir[i].bb1=1.; WHICH IS DEFINITELY WRONG

             //Because later we are doing operations on the slope, we should choose a large enough number that does not overflow
             // a good idea could be to choose a slope so that the almost vertical line makes less than epsilon=0.1 error across the whole field
             // so that division is effectively vertical (or horizontal)
             // this is a line that has slope par.sizex/epsilon, let's add a little bit to be extra safe (times 2)
             // with a field size of 1000 and epsilon = 0.1 -> m= 10000 that's ok small for more calculations
             double large_enough_number = (2.*(double)par.sizex)/0.1;
             double random_plus_or_minus_1 = -1+2*(int)(2.*RANDOM());
             if(sxx>syy){
               celldir[i].bb1 = 0.;
               celldir[i].aa1=ymean;
               celldir[i].bb2 = random_plus_or_minus_1*large_enough_number;
               celldir[i].aa2=ymean-celldir[i].bb2*xmean;
             }
             else if(syy>sxx){
               celldir[i].bb1 = random_plus_or_minus_1*large_enough_number;
               celldir[i].aa1=ymean-xmean*celldir[i].bb1;
               celldir[i].bb2 = random_plus_or_minus_1*0.00001;
               celldir[i].aa2=ymean;
             }
             else{
               celldir[i].bb1 = (RANDOM() <0.5)? 0. : (random_plus_or_minus_1*large_enough_number); //if sxx==syy we randomise vertical or horizontal
               if(celldir[i].bb1 > 1. || celldir[i].bb1 < 1.){
                 celldir[i].aa1=ymean-xmean*celldir[i].bb1;
                 celldir[i].bb2 = random_plus_or_minus_1*0.00001;
                 celldir[i].aa2=ymean;
              }else{
                celldir[i].aa1=ymean;
                celldir[i].bb2 = random_plus_or_minus_1*large_enough_number;
                celldir[i].aa2=ymean-celldir[i].bb2*xmean;
              }
             }
          }
         }

       //}

       /* bevrijd gealloceerd geheugen */
       free(sumx);
       free(sumy);
       free(sumxx);
       free(sumxy);
       free(sumyy);
       free(n);

       return celldir;

}

// See function above
//Buggy function, if cells are across border they don't get the right moments.
// (e.g. standard dev will be spread throughout the whole field).
// because the moments are calculated without taking cell continuity into account
// Needs a better strategy:
// How about we go through field and
// idea 1) fill up arrays with cell position
// idea two, when we find a new idea we spin around it in circle until we get the whole cell
// ... how much overhead is that? - far more than just one loop through the CA for sure.

//maybe a better heuristic could be:
// keep track of the previous x you added, if it is much larger or much smaller than
// the one we find now, we subtract (or add) one field size, this should help with continuity
// (unless we get giant cells (or very small fields))
Dir *CellularPotts::FindCellDirections(void) const
{
  // double array - allocated (why c style?) the size population
  double *sumx=0,*sumy=0;
  double *sumxx=0,*sumxy=0,*sumyy=0;
  double *n=0;

  double xmean=0,ymean=0,sxx=0,sxy=0,syy=0;
  double D,lb1=0,lb2=0;

  Dir *celldir;

  /* Allocation of sufficient memory space */
  if( (sumx= (double *)malloc(cell->size()*sizeof(double)))==NULL)
    MemoryWarning();
  else
    if( (sumy= (double *)malloc(cell->size()*sizeof(double)))==NULL)
      MemoryWarning();
    else
      if ((sumxx=(double *)malloc(cell->size()*sizeof(double)))==NULL)
	MemoryWarning();
      else
	if((sumxy=(double *)malloc(cell->size()*sizeof(double)))==NULL)
	  MemoryWarning();
	else
	  if((sumyy=(double *)malloc(cell->size()*sizeof(double)))==NULL)
	    MemoryWarning();
	  else
	    if((n=(double *)malloc(cell->size()*sizeof(double)))==NULL)
	      MemoryWarning();


  if ( !(celldir=new Dir[cell->size()]) )
    MemoryWarning();


  /* Initialization of the variables */

  for (int i=0;i<(int)cell->size();i++) {

    sumx[i]=0.;
    sumy[i]=0.;
    sumxx[i]=0.;
    sumxy[i]=0.;
    sumyy[i]=0.;
    n[i]=0L;

  }


  /* Find sumx, sumy, sumxx and sumxy for all cells */

  for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++)
      if (sigma[x][y]>0) {
	sumx[0]+=(double)x;  //pos 0 contains global statistics
	sumy[0]+=(double)y;
	sumxx[0]+=(double)x*x;
	sumxy[0]+=(double)x*y;
	sumyy[0]+=(double)y*y;

	n[0]++;

	sumx[sigma[x][y]]+=(double)x; //pos indicised by sigma - local stats about x direction
	sumy[sigma[x][y]]+=(double)y;

	sumxx[sigma[x][y]]+=(double)x*x;
	sumxy[sigma[x][y]]+=(double)x*y; // mixed coordinate
	sumyy[sigma[x][y]]+=(double)y*y;

	n[sigma[x][y]]++;

      }

  /* Compute the principal axes for all cells */

  {
    for (int i=0;i<(int)cell->size();i++) {
      // why this? of course moments are ill defined for small cells
      // maybe it's no problem because they don't divide
      if (n[i]>10) {

	xmean=((double)sumx[i])/((double)n[i]); //mean x
    ymean=((double)sumy[i])/((double)n[i]); //mean y

	sxx=(double)(sumxx[i])-((double)(sumx[i]*sumx[i]))/(double)n[i];
	sxx=sxx/(double)(n[i]-1);  //standard deviation x (estimate = (E(X^2) - E(X)^2)/ (N*(N-1)) )

	sxy=(double)(sumxy[i])-((double)(sumx[i]*sumy[i]))/(double)n[i];
    sxy=sxy/(double)(n[i]-1);  //standard deviation y

	syy=(double)(sumyy[i])-((double)(sumy[i]*sumy[i]))/(double)n[i];
	syy=syy/(double)(n[i]-1);

	D=sqrt( (sxx+syy)*(sxx+syy)-4.*(sxx*syy-sxy*sxy) ); // this is a discriminant
	lb1=(sxx+syy+D)/2.;lb2=(sxx+syy-D)/2.;              // this is a PCA ???
	celldir[i].lb1=lb1; celldir[i].lb2=lb2;
      }
      if (sxy==0.0)
	celldir[i].bb1=1.;
      else
	celldir[i].bb1=sxy/(lb1-syy);

      if (fabs(celldir[i].bb1)<.00001) {
	if (celldir[i].bb1>0.)
	  celldir[i].bb1=.00001;
	else
	  celldir[i].bb1=-.00001;
      }

      celldir[i].aa1=ymean-xmean*celldir[i].bb1;
      celldir[i].bb2= (-1.)/celldir[i].bb1;

      celldir[i].aa2=ymean-celldir[i].bb2*xmean;
    }

  }

  /* bevrijd gealloceerd geheugen */
  free(sumx);
  free(sumy);
  free(sumxx);
  free(sumxy);
  free(sumyy);
  free(n);

  return celldir;

}

void CellularPotts::ShowDirections(Graphics &g, const Dir *celldir) const
{
  int i;

  if (cell->size()>1)
    for (i=1;i<(int)cell->size();i++)
      g.Line(0,(int)(2*celldir[i].aa1),sizex*2,(int)((celldir[i].aa1+celldir[i].bb1*sizey)*2),2);

}

//a make-over for the function DivideCells
// DON'T USE IT, at best it's the same as the other, at worst it might be incomplete
// commented out to make sure you don't accidentally use it!
// vector<int> CellularPotts::DivideCells2(vector<bool> which_cells)
// {
//   //cerr<<"Hello begin DivideCells"<<endl;
//   int sigmaneigh;
//   // for the cell directions
//   Dir *celldir=0;
//
//   // Allocate space for divisionflags
//   vector<int> divflags( cell->size()*2 + 5 ); // automatically initialised to zero
//
//   //curious: here it also complains when which_cells contains as many cells as the whole cell vector
//   //but further down, it will just divide all cells if which_cells is empty...
//   // the comment above really is not true
//   if ( !(which_cells.size()==0 || which_cells.size() >= cell->size()) ) {
//     cerr<<"which_cells.size()="<<which_cells.size()<<endl;
//     throw "In CellularPotts::DivideCells, Too few elements in vector<int> which_cells.";
//   }
//   // division
//   {
//     celldir=FindCellDirections3(); //find cell directions here
//     //celldir=FindCellDirections2();  // notice that this is called once, the first time any cell is divided
//     //celldir=FindCellDirections(); BUGGY in two ways: doesn't handle division plane across boundaries at all, it misplaces the normal one as well.
//
//     for (int i=1;i<sizex-1;i++) for (int j=1;j<sizey-1;j++)
//       if (sigma[i][j]>0) // i.e. not medium and not border state (-1)
//       {
//         // Pointer to mother. Warning: Renew pointer after a new
//         // cell is added (push_back). Then, the array *cell is relocated and
//         // the pointer will be lost...
//         Cell *motherp=&((*cell)[sigma[i][j]]); //mother points to the cell holding this pixel
//         Cell *daughterp;
//
//         // Divide if NOT medium and if DIV bit set or divide_always is set
//         // if which_cells is given, di+nx>0 && i+nx[k]<sizex-1 divide only if the cell
//         // is marked in which_cells.
//         if( !which_cells.size() || which_cells[motherp->sigma] ){
//           //divflags is a vector that works like this: divflags[ mother_sigma ] = daughter_sigma
//           //if first time we get this mother then divflags at pos mother_sigma is 0
//           if( !(divflags[ motherp->Sigma() ]) ){
//             // add daughter cell, copying states of mother
//
//             //we first check if we can recycle some position already exisiting in the vector
//             //such position would come from a cell that has previously apoptosed
//             vector<Cell>::iterator c;
//             bool replaced=false;
//             for ( c=cell->begin(), c++ ; c!=cell->end(); c++){
//               if(c->AliveP() == false && c->TargetArea() <= 0 && c->Area() == 0 ){
//                 //we recycle this sigma for the new cell
//                 daughterp=new Cell(*(motherp->owner), motherp->getTau(), c->Sigma()); //set recycled sigma
//                 daughterp->CellBirth(*motherp);
//                 *c = *daughterp;    // notice that the operator = (equal) is overloaded, see cell.h
//                 replaced=true;
//                 break;
//               }
//             }
//             if( !replaced ){
//               daughterp=new Cell(*(motherp->owner));  //this calls  Cell(...){ owner=&who; ConstructorBody()}
//               // this is what prints Tomato2
//               daughterp->CellBirth(*motherp);
//               cell->push_back(*daughterp);  //this calls default copy constructor Cell(const Cell &src)
//               // prints "Tomato"
//               //this puts new cells at the end of array if there was no space to recycle
//               // renew pointer to mother (because after push_back memory might be relocated)
//               motherp=&((*cell)[sigma[i][j]]);
//             }
//
//             divflags[ motherp->Sigma() ]=daughterp->Sigma(); //daughtersigma is set to newest sigma in ConstructorBody of Cell
//             delete daughterp;
//             // array may be relocated after "push_back"
//             // renew daughter pointers
//             if(replaced) daughterp=&(*c);
//             else daughterp=&(cell->back());
//             //daughterp=&(cell->back());
//           }else{
//             daughterp=&( (*cell)[ divflags[motherp->Sigma()] ] );
//           }
//
//           // if site is below the minor axis of the cell: sigma of new cell
//           // to properly choose this we have to check where this pixel is
//           int checki=i;
//           int checkj=j;
//           if(par.periodic_boundaries){
//             //Check if this pixel is closer to mean pos when wrapped, we wrap it
//             double meanx=celldir[motherp->sigma].meanx;
//             double meany=celldir[motherp->sigma].meany;
//
//             if( (checki-meanx)>0 && (checki-meanx)>(meanx-(checki-(par.sizex-2))) ) {
//               checki-=(par.sizex-2);
//               //cerr<<"celldiv passb1"<<endl;
//             }
//             else if( (meanx-checki)>0 && (meanx-checki)>(checki+(par.sizex-2)-meanx) ){
//               checki+=(par.sizex-2);
//               //cerr<<"celldiv passb2"<<endl;
//             }
//             if(  checkj-meany>0  &&  checkj-meany >(meany-(checkj-(par.sizey-2))) ){
//               checkj-=(par.sizey-2);
//               //cerr<<"celldiv passb3"<<endl;
//             }
//             else if( meany-checkj>0   &&  meany-checkj >(checkj+(par.sizey-2)-meany) ){
//               checkj+=(par.sizey-2);
//               //cerr<<"celldiv passb4"<<endl;
//             }
//           }
//           if( checkj>( (int)( celldir[motherp->sigma].aa2 + celldir[motherp->sigma].bb2*(double)checki) ) ){
//             motherp->DecrementArea();
//             motherp->DecrementTargetArea();
//             motherp->RemoveSiteFromMoments(i,j);
//
//             sigma[i][j]=daughterp->Sigma();  // WHERE is daughterp->Sigma() defined?
//             daughterp->IncrementArea();
//             daughterp->IncrementTargetArea();
//             daughterp->AddSiteToMoments(i,j);
//
//             //go through neighbourhood to update contacts
//             // to new daughter contacts we now pass duration from mother
//             // sigma[i][j] is daughter, sigmaneigh can be daughter, mother, medium, someone else
//             for (int k=1; k<=n_nb; k++){
//               //if wrapped boundaries we wrap i+nx[k] and j+ny[k] around (if needed)
//               //if fiexed boundaries, we exclude them from the neigh counting
//               int neix=i+nx[k];
//               int neiy=j+ny[k];
//               if(neix<=0 || neix>=sizex-1 || neiy<=0 || neiy>=sizey-1){
//                 if( par.periodic_boundaries ){
//                   if(neix<=0) neix=sizex-2+neix;
//                   if(neix>=sizex-1) neix=neix-sizex+2;
//                   if(neiy<=0) neiy=sizey-2+neiy;
//                   if(neiy>=sizey-1) neiy=neiy-sizey+2;
//                 }else{
//                   continue;
//                 }
//               }
//
//
//               sigmaneigh = sigma[ neix ][ neiy ];
//               //if sigmaneigh is not sigma, we update the contact of daughter cell with it,
//               //and the contact of that cell with daughter (provided it is not medium)
//               if( sigmaneigh  != sigma[i][j] ){
//                 //cout<<sigmaneigh<<" "<<sigma[i][j]<<" +1"<<endl;
//                 //cerr<<"Hello 0.1"<<endl;
//                 (*cell)[sigma[i][j]].updateNeighbourBoundary(sigmaneigh,1);
//                 //take duration from mother iff sigmaneigh is not mother
//                 if(sigmaneigh!=motherp->Sigma() && sigmaneigh!=MEDIUM)
//                   (*cell)[sigma[i][j]].SetNeighbourDurationFromMother(sigmaneigh, motherp->returnDuration(sigmaneigh) );
//
//                 //cerr<<"Hello 0.2"<<endl;
//                 //also cell to which sigmaneigh belongs must be updated, if it is not medium
//                 if(sigmaneigh){
//                   //cerr<<"Hello 0.21, pos i,j: "<< i<<","<<j<<endl;
//                   //cerr<<"Sigmaneigh (method) is: " << (*cell)[sigmaneigh].Sigma()<< " and should be (from CA) "<< sigmaneigh;
//                   //cerr<< " from pos neix,neiy "<< neix<<","<<neiy<< endl;
//                   //cerr<< "Also before calculations nei x and y: "<<i+nx[k]<<","<<j+ny[k]<<endl;
//                   (*cell)[sigmaneigh].updateNeighbourBoundary(sigma[i][j],1);
//                   //cerr<<"Hello 0.22"<<endl;
//                   if( sigmaneigh!=motherp->Sigma() ){
//                     //cerr<<"Hello 0.23"<<endl;
//                     (*cell)[sigmaneigh].SetNeighbourDurationFromMother( sigma[i][j], motherp->returnDuration(sigmaneigh) );
//                   }
//                 }
//                 //cerr<<"Hello 0.3"<<endl;
//                 if (sigmaneigh!=motherp->Sigma()){
//                   //cout<<sigmaneigh<<" "<<motherp->Sigma()<<" -1"<<endl;
//                   motherp->updateNeighbourBoundary(sigmaneigh,-1);
//                   //cerr<<"Hello 0.4"<<endl;
//                   if(sigmaneigh)
//                     (*cell)[sigmaneigh].updateNeighbourBoundary(motherp->Sigma(),-1);
//                   //cerr<<"Hello 0.5"<<endl;
//                 }
//               }else//sigmaneigh==sigma[i][j] This pixel has already become a daughter pixel,
//                 //remove from contacts between mother and daughter
//               {
//                 motherp->updateNeighbourBoundary(sigmaneigh,-1);
//                 (*cell)[sigmaneigh].updateNeighbourBoundary(motherp->Sigma(),-1);
//               }
//             }
//           }
//         }
//       }
//   }
//   if (celldir)
//     delete[] (celldir);
//
//   return divflags;
// }


// TODO: This should probably be in cell.h
// TODO: This should probably be in cell.h
int CellularPotts::DivideCell(int cell_sigma, BoundingBox box) {

  int sigmaneigh;

  // for the cell directions
  Dir *celldir=0;
  vector<int> toprint;

  Cell *motherp = &((*cell)[cell_sigma]);
  Cell *daughterp;
  /* division */

  //we first check if we can recycle some position already exisiting in the vector
  //such position would come from a cell that has previously apoptosed
  vector<Cell>::iterator c;
  bool replaced=false;
  for ( c=cell->begin(), c++ ; c!=cell->end(); c++){
    if(!c->AliveP() && c->TargetArea() <= 0 && c->Area() == 0 ){
      //we recycle this sigma for the new cell
      //set recycled sigma
      daughterp=new Cell(*(motherp->owner), motherp->getTau(), c->Sigma());
      daughterp->CellBirth(*motherp);
      *c = *daughterp;    // notice that the operator = (equal) is overloaded, see cell.h
      replaced=true;
      break;
    }
  }
  if( !replaced ){
    //cout << "We do not recycle, calling function new"<< endl;
    // THIS USED TO BE ABOVE, WHERE THE SIGN *** IS !!!
    //MAKES NEW CELL AT THE END OF ARRAY
    daughterp=new Cell(*(motherp->owner));  //this calls  Cell(...){ owner=&who; ConstructorBody()}
    int momcol=motherp->colour;
    daughterp->CellBirth(*motherp);
    cell->push_back(*daughterp);  //this calls default copy constructor Cell(const Cell &src)
    // prints "Tomato"
    //this puts new cells at the end of array if there was no space to recycle
    // renew pointer to mother (because after push_back memory might be relocated)
    motherp=&((*cell)[cell_sigma]);
    if(motherp->colour!=momcol) cerr<<"this is the problem"<<endl;
  }
  // renew pointers
  delete daughterp;
  if(replaced){
    daughterp=&(*c);
    //cerr<<"mother sigma: "<<motherp->Sigma()<<", daughter sigma"<<daughterp->Sigma()<<endl;
  }
  else {
    //cerr<<"mother sigma: "<<motherp->Sigma()<<", daughter sigma"<<daughterp->Sigma()<<endl;
    daughterp=&(cell->back());
  }

  int pixel_count = 0;
  for (int i=box.minx; i <= box.maxx; i++) for (int j = box.miny; j <= box.maxy; j++) {
    if (sigma[i][j] == motherp->sigma) {
      ++pixel_count;
      /* Now the actual division takes place */

      /* If celldirections where not yet computed: do it now */
      if (!celldir)
        celldir = FindCellDirections3();
      // if site is below the minor axis of the cell: sigma of new cell
      // to properly choose this we have to check where this pixel is
      int checki = i;
      int checkj = j;

      if (checkj > ((int) (celldir[motherp->sigma].aa2 + celldir[motherp->sigma].bb2 * (double) checki))) {
        motherp->DecrementArea();
        motherp->RemoveSiteFromMoments(i, j);
        sigma[i][j] = daughterp->Sigma();
        daughterp->IncrementArea();
        daughterp->AddSiteToMoments(i, j);
        //go through neighbourhood to update contacts
        // to new daughter contacts we now pass duration from mother
        // sigma[i][j] is daughter, sigmaneigh can be daughter, mother, medium, someone else
        for (int k = 1; k <= n_nb; k++) {
          //if wrapped boundaries we wrap i+nx[k] and j+ny[k] around (if needed)
          //if fiexed boundaries, we exclude them from the neigh counting
          int neix = i + nx[k];
          int neiy = j + ny[k];
          sigmaneigh = sigma[neix][neiy];
          //if sigmaneigh is not sigma, we update the contact of daughter cell with it,
          //and the contact of that cell with daughter (provided it is not medium)
          if (sigmaneigh != sigma[i][j]) {

            //update the edgeSetVector
            edgeSetVector.insert({i, j, neix, neiy});
            edgeSetVector.insert({neix, neiy, i, j});

            (*cell)[sigma[i][j]].updateNeighbourBoundary(sigmaneigh, 1);
            //take duration from mother iff sigmaneigh is not mother
            if (sigmaneigh != motherp->Sigma() && sigmaneigh != MEDIUM)
              (*cell)[sigma[i][j]].SetNeighbourDurationFromMother(sigmaneigh, motherp->returnDuration(sigmaneigh));

            //also cell to which sigmaneigh belongs must be updated, if it is not medium
            if (sigmaneigh) {
              (*cell)[sigmaneigh].updateNeighbourBoundary(sigma[i][j], 1);
              if (sigmaneigh != motherp->Sigma()) {
                (*cell)[sigmaneigh].SetNeighbourDurationFromMother(sigma[i][j], motherp->returnDuration(sigmaneigh));
              }
            }

            if (sigmaneigh != motherp->Sigma()) {
              motherp->updateNeighbourBoundary(sigmaneigh, -1);

              if (sigmaneigh)
              (*cell)[sigmaneigh].updateNeighbourBoundary(motherp->Sigma(), -1);
            }
          } else//sigmaneigh==sigma[i][j] This pixel has already become a daughter pixel,
            //remove from contacts between mother and daughter
          {
            //update the edgeSetVector
            edgeSetVector.erase({i, j, neix, neiy});
            edgeSetVector.erase({neix, neiy, i, j});

            motherp->updateNeighbourBoundary(sigmaneigh, -1);
            (*cell)[sigmaneigh].updateNeighbourBoundary(motherp->Sigma(), -1);
          }
        } // end neighbour loop
      }
    }
  }
  int expected_area = motherp->Area() + daughterp->Area();
  if (pixel_count != expected_area) {
    cerr << "Something went wrong when dividing cell " << motherp->Sigma() << endl;
    cerr << "Expected cell area is " << expected_area << " but only " << pixel_count << " pixels were found inside bounding box";
    cerr << "Terminating the program";
    exit(1);
  }
  if (celldir)
    delete[] (celldir);

  return daughterp->sigma;
}


vector<int> CellularPotts::DivideCells(vector<bool> which_cells)
{

  int sigmaneigh;

  // for the cell directions
  Dir *celldir=0;

  // Allocate space for divisionflags
  vector<int> divflags( cell->size()*2 + 5 ); // automatically initialised to zero
  vector<int> toprint;

  //curious: here it also complains when which_cells contains as many cells as the whole cell vector
  //but further down, it will just divide all cells if which_cells is empty...
  // the comment above really is not true
  if ( !(which_cells.size()==0 || which_cells.size() >= cell->size()) ) {
    throw "In CellularPotts::DivideCells, Too few elements in vector<int> which_cells.";
  }

  /* division */
  {
    for (int i=1;i<sizex-1;i++) for (int j=1;j<sizey-1;j++)
      if (sigma[i][j]>0) // i.e. not medium and not border state (-1)
	    {
	      Cell *motherp=&((*cell)[sigma[i][j]]); //mother points to the cell holding this pixel
	      Cell *daughterp;
        //if(motherp->colour!=2 &&motherp->colour!=3) cerr<<"mother has wrong colour: "<<motherp->colour<<" "<<motherp->Sigma()<<endl;
	      // Divide if NOT medium and if DIV bit set or divide_always is set; if which_cells is given  divide only if the cell is listed
	      if( !which_cells.size() || which_cells[motherp->sigma] ){
          //divflags is a vector that works like this: divflags[ mother_sigma ] = daughter_sigma
          //if first time we get this mother then divflags at pos mother_sigma is 0
          if( !(divflags[ motherp->Sigma() ]) ){
            //we first check if we can recycle some position already exisiting in the vector
            //such position would come from a cell that has previously apoptosed
            vector<Cell>::iterator c;
            bool replaced=false;
            for ( c=cell->begin(), c++ ; c!=cell->end(); c++){
              if(c->AliveP() == false && c->TargetArea() <= 0 && c->Area() == 0 ){
                //we recycle this sigma for the new cell
                //set recycled sigma
                daughterp=new Cell(*(motherp->owner), motherp->getTau(), c->Sigma());
                daughterp->CellBirth(*motherp);
                *c = *daughterp;    // notice that the operator = (equal) is overloaded, see cell.h
                replaced=true;
                break;
              }
            }
            if( !replaced ){
              //cout << "We do not recycle, calling function new"<< endl;
              // THIS USED TO BE ABOVE, WHERE THE SIGN *** IS !!!
              //MAKES NEW CELL AT THE END OF ARRAY
              daughterp=new Cell(*(motherp->owner));  //this calls  Cell(...){ owner=&who; ConstructorBody()}
              int momcol=motherp->colour;
              daughterp->CellBirth(*motherp);
              cell->push_back(*daughterp);  //this calls default copy constructor Cell(const Cell &src)
                                            // prints "Tomato"
                                            //this puts new cells at the end of array if there was no space to recycle
              // renew pointer to mother (because after push_back memory might be relocated)
              motherp=&((*cell)[sigma[i][j]]);
              if(motherp->colour!=momcol) cerr<<"this is the problem"<<endl;
            }

        divflags[ motherp->Sigma() ]=daughterp->Sigma(); //daughtersigma is set to newest sigma in ConstructorBody of Cell
        delete daughterp;

        //cout << "Hello2" << endl;
		// array may be relocated after "push_back"

		// renew daughter pointers

        //trying stuff, line below is how it was
        if(replaced){
          daughterp=&(*c);
          //cerr<<"mother sigma: "<<motherp->Sigma()<<", daughter sigma"<<daughterp->Sigma()<<endl;
        }
        else {
          //cerr<<"mother sigma: "<<motherp->Sigma()<<", daughter sigma"<<daughterp->Sigma()<<endl;
          daughterp=&(cell->back());
        }
        //daughterp=&(cell->back());

		/* administration on the onset of mitosis */

		/* Ancestry is taken care of in copy constructor of Cell
		   see cell.hh: Cell(const Cell &src, bool newcellP=false) : Cytoplasm(src) {} */

		/* inherit  polarity of mother */
		// All that needs to be copied is copied in the copy constructor
		// of Cell and in the default copy constr. of its base class Cytoplasm
		// note: also the celltype is inherited
	      }else{
            daughterp=&( (*cell)[ divflags[motherp->Sigma()] ] );
        }
	      /* Now the actual division takes place */

	      /* If celldirections where not yet computed: do it now */
	      if (!celldir)
            celldir=FindCellDirections3();
            //celldir=FindCellDirections2();  // notice that this is called once, the first time any cell is divided
            //celldir=FindCellDirections(); BUGGY in two ways: doesn't handle division plane across boundaries at all, it misplaces the normal one as well.


	      // if site is below the minor axis of the cell: sigma of new cell
          // to properly choose this we have to check where this pixel is
          int checki=i;
          int checkj=j;
          if(par.periodic_boundaries){
            //I guess it is the usual check, if this pixel is closer to mean pos when wrapped, we wrap it
            // remember that i indices to sizex, and j to sizey
            //int meanx=motherp->meanx; //THIS MOVES THE MEAN while mother and daugther are created, which is buggy, add mean to directions2
            //int meany=motherp->meany;
            double meanx=celldir[motherp->Sigma()].meanx;
            double meany=celldir[motherp->Sigma()].meany;

            if( ((checki-meanx)>0) && ((checki-meanx)>(meanx-(checki-(par.sizex-2)))) ) {
              checki-=(par.sizex-2);
            }
            else if( ((meanx-checki)>0) && ((meanx-checki)>(checki+(par.sizex-2)-meanx)) ){
              checki+=(par.sizex-2);
            }
            if(  ((checkj-meany)>0)  &&  ((checkj-meany) >(meany-(checkj-(par.sizey-2)))) ){
              checkj-=(par.sizey-2);
            }
            else if( ((meany-checkj)>0)   &&  ((meany-checkj) >(checkj+(par.sizey-2)-meany)) ){
              checkj+=(par.sizey-2);
            }
          }

          if( checkj>( (int)( celldir[motherp->sigma].aa2 + celldir[motherp->sigma].bb2*(double)checki) ) ){
            motherp->DecrementArea();
            motherp->RemoveSiteFromMoments(i,j);
            sigma[i][j]=daughterp->Sigma();  // WHERE is daughterp->Sigma() defined?
            daughterp->IncrementArea();
            daughterp->AddSiteToMoments(i,j);
            //go through neighbourhood to update contacts
            // to new daughter contacts we now pass duration from mother
            // sigma[i][j] is daughter, sigmaneigh can be daughter, mother, medium, someone else
            for (int k=1; k<=n_nb; k++){
              //if wrapped boundaries we wrap i+nx[k] and j+ny[k] around (if needed)
              //if fiexed boundaries, we exclude them from the neigh counting
              int neix=i+nx[k];
              int neiy=j+ny[k];
              if(neix<=0 || neix>=sizex-1 || neiy<=0 || neiy>=sizey-1){
                if( par.periodic_boundaries ){
                  if(neix<=0) neix=sizex-2+neix;
                  if(neix>=sizex-1) neix=neix-sizex+2;
                  if(neiy<=0) neiy=sizey-2+neiy;
                  if(neiy>=sizey-1) neiy=neiy-sizey+2;
                }else{
                  continue;
                }
              }
              sigmaneigh = sigma[ neix ][ neiy ];
              //if sigmaneigh is not sigma, we update the contact of daughter cell with it,
              //and the contact of that cell with daughter (provided it is not medium)
              if( sigmaneigh  != sigma[i][j] ){

                //update the edgeSetVector
                edgeSetVector.insert({i,j,neix,neiy});
                edgeSetVector.insert({neix,neiy,i,j});

                (*cell)[sigma[i][j]].updateNeighbourBoundary(sigmaneigh,1);
                //take duration from mother iff sigmaneigh is not mother
                if(sigmaneigh!=motherp->Sigma() && sigmaneigh!=MEDIUM)
                  (*cell)[sigma[i][j]].SetNeighbourDurationFromMother(sigmaneigh, motherp->returnDuration(sigmaneigh) );

                //also cell to which sigmaneigh belongs must be updated, if it is not medium
                if(sigmaneigh){
                  (*cell)[sigmaneigh].updateNeighbourBoundary(sigma[i][j],1);
                  if( sigmaneigh!=motherp->Sigma() ){
                      (*cell)[sigmaneigh].SetNeighbourDurationFromMother( sigma[i][j], motherp->returnDuration(sigmaneigh) );
                    }
                }

                if (sigmaneigh!=motherp->Sigma()){
                  motherp->updateNeighbourBoundary(sigmaneigh,-1);

                  if(sigmaneigh)
                    (*cell)[sigmaneigh].updateNeighbourBoundary(motherp->Sigma(),-1);
                }
              }else//sigmaneigh==sigma[i][j] This pixel has already become a daughter pixel,
                      //remove from contacts between mother and daughter
                {
                  //update the edgeSetVector
                  edgeSetVector.erase({i,j,neix,neiy});
                  edgeSetVector.erase({neix,neiy,i,j});

                  motherp->updateNeighbourBoundary(sigmaneigh,-1);
                  (*cell)[sigmaneigh].updateNeighbourBoundary(motherp->Sigma(),-1);
                }
            } // end neighbour loop
	     }
	   }
	 }
  }
  if (celldir)
    delete[] (celldir);

  return divflags;

}

bool CellularPotts::PlaceOneCellsAtXY(int posx,int posy, int cellsize, int cellsigma)
{
  int radsq=(int)(((double)cellsize)/3.14);
  cerr<<"Placing cell at posx: "<<posx<<" posy:"<<posy<<", size: "<<par.size_init_cells<<endl;

  for (int x=posx-cellsize;x<posx+cellsize;x++){
    for (int y=posy-cellsize;y<posy+cellsize;y++){
      if( (x-posx)*(x-posx)+(y-posy)*(y-posy)<radsq && x>=1 && x<sizex-1 && y>=1 && y<sizey-1){  // circle
//       if( abs(x-posx)<sqrt(radsq) && abs(y-posy)<sqrt(radsq)+5 && x>=1 && x<sizex-1 && y>=1 && y<sizey-1){ //rectangle
//        if( ((x-posx)*(x-posx)/((double)radsq + 1500) +(y-posy)*(y-posy)/(-1500+(double)radsq))<1. && x>=1 && x<sizex-1 && y>=1 && y<sizey-1) { // ellipse
        if( sigma[x][y] ){
          cerr << "PlaceOneCellsAtXY(): Error. There is a cell here already?"<< endl;
          exit(1);
        }else{
          sigma[x][y]=cellsigma;
        }
      }
    }
  }
  return true;
}

int CellularPotts::PlaceOneCellRandomly(int sig, int cellsize)
{
  bool overlap=false;
  int radsq=(int)(((double)cellsize)/3.14);
  int check=0;
  int placedpix=0;
  int x0,y0;

  while (check<100)
  {
    overlap=false;
    placedpix=0;
    x0=radsq+RandomNumber(sizex-2*radsq); //should avoid putting them across boundaries, radsq
                                                  //(radius square is overdoing it,
                                                  //but it's ok given the bugs that happen
                                                  //when you put cells across boundaries :S )
    y0=radsq+RandomNumber(sizey-2*radsq);
    // check overlap
    for (int x=x0-cellsize;x<x0+cellsize;x++){
      for (int y=y0-cellsize;y<y0+cellsize;y++){
        if( (x-x0)*(x-x0)+(y-y0)*(y-y0)<radsq && x>=1 && x<sizex-1 && y>=1 && y<sizey-1)  // circle
        {
          if( sigma[x][y] ){
            overlap=true;
            check++;
            break;
          }
        }
      }
      if(overlap)
         break;
    }

    if(!overlap){
      for (int x=x0-cellsize;x<x0+cellsize;x++){
        for (int y=y0-cellsize;y<y0+cellsize;y++){
          if((x-x0)*(x-x0)+(y-y0)*(y-y0)<radsq && x>=1 && x<sizex-1 && y>=1 && y<sizey-1){ //circle
            sigma[x][y]=sig;
            placedpix++;
          }
        }
      }
      break;
    }
  }

  return placedpix;

}

int CellularPotts::PlaceCellsRandomly(int n, int cellsize)
{
    bool overlap=false;
    int radsq=(int)(((double)cellsize)/3.14);
    int count=0, check=0;

    while(count<n && check<100)
    {
        overlap=false;
        int x0=radsq+RandomNumber(sizex-2*radsq); //should avoid putting them across boundaries, radsq
                                                  //(radius square is overdoing it,
                                                  //but it's ok given the bugs that happen
                                                  //when you put cells across boundaries :S )
        int y0=radsq+RandomNumber(sizey-2*radsq);
        // check overlap
        for (int x=x0-cellsize;x<x0+cellsize;x++){
           for (int y=y0-cellsize;y<y0+cellsize;y++){
                if( (x-x0)*(x-x0)+(y-y0)*(y-y0)<radsq && x>=1 && x<sizex-1 && y>=1 && y<sizey-1)  // circle
//              if( abs(x-x0)<sqrt(radsq) && abs(y-y0)<sqrt(radsq)+5 && x>=1 && x<sizex-1 && y>=1 && y<sizey-1) //rectangle
//             if( ((x-x0)*(x-x0)/((double)radsq-20) +(y-y0)*(y-y0)/(100+(double)radsq))<1. && x>=1 && x<sizex-1 && y>=1 && y<sizey-1)  // ellipse

               {
                 if( sigma[x][y] ){
                       overlap=true;
                       check++;
                       //cerr << "Overlap. count: "<<n<<" check: "<<check<< endl;
                       break;
                 }
               }

           }
           if(overlap)
               break;
        }

        if(!overlap)
        {
            check=0;
            count++;
            //cerr << "No overlap. count: "<<n<<" check: "<<check<< endl;
            for (int x=x0-cellsize;x<x0+cellsize;x++)
            {
                for (int y=y0-cellsize;y<y0+cellsize;y++)
                {
//                   if( abs(x-x0)<sqrt(radsq) && abs(y-y0)<sqrt(radsq)+5 && x>=1 && x<sizex-1 && y>=1 && y<sizey-1)   //rectangle
                   if((x-x0)*(x-x0)+(y-y0)*(y-y0)<radsq && x>=1 && x<sizex-1 && y>=1 && y<sizey-1) //circle
//                  if( ((x-x0)*(x-x0)/((double)radsq-20) +(y-y0)*(y-y0)/(100+(double)radsq))<1.  && x>=1 && x<sizex-1 && y>=1 && y<sizey-1) //ellipses
                        sigma[x][y]=count;

                }

            }
        }
    }

    cerr<<"Placed "<<count<<" cells out of "<<n<<" requested."<<endl;
    //exit(1);
    return count;
}

// Places cells at regular distance from one another:
// For square cells of size s, the spatial occupation is sqrt(s).
// Since we have to place n of them, we use a square of space of size sqrt(n)*(sqrt(s) + a_little_bit),
// centered at the center of grid, which means that the upper left corner of the first cell is at
// x=(sizex-sqrt(n)*(sqrt(s) + a_little_bit))/2
int CellularPotts::PlaceCellsOrderly(int n_cells,int size_cells)
{
    int count=0;
    int a_little_bit=2;
    int beginx, endx, beginy, endy;
    int step, avrg_area;
    int smaller_dimension=( par.sizex < par.sizey)?par.sizex:par.sizey;
    int sqrt_n_cells = 1+sqrt(n_cells);
    //to avoid having 49 cells when you want 50, I'm rounding sqrt(n_cells) to +1
    if( (  sqrt_n_cells  *(sqrt(size_cells) + a_little_bit))>smaller_dimension ){
      std::cerr << "PlaceCellsOrderly(): Error. Too many cells or too large size?" << '\n';
      exit(1);
    }

    step = ( sqrt(size_cells) + a_little_bit );

    avrg_area=0;

    // each x,y point denotes the upper left corner of a cell
    // with i,j we run through the cell
    // for different initial conditions we will have oders of putting the cells
    //this is quite extendable
    vector<int> v_order_x;
    vector<int> v_order_y;

    if(par.cell_placement){
      std::cerr << "yay, placing cells where I want!"<<endl;
      switch( par.cell_placement ){  
        case 1:
          beginx = a_little_bit;
          endx =   a_little_bit + sqrt_n_cells*(sqrt(size_cells) + a_little_bit);
          beginy = (par.sizey -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
          endy =   (par.sizey +  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
          break;
        case 2:
          beginx = (par.sizex -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
          endx =   (par.sizex +  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
          beginy = a_little_bit;
          endy =   a_little_bit + sqrt_n_cells*(sqrt(size_cells) + a_little_bit);
          break;
        case 3:
          beginx = par.sizex -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit)- a_little_bit;
          endx =   par.sizex - a_little_bit;
          beginy = (par.sizey -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
          endy =   (par.sizey +  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
          break;
        case 4:
          beginx = (par.sizex -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
          endx =   (par.sizex +  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
          beginy = par.sizey -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit) - a_little_bit;
          endy =   par.sizey - a_little_bit;
          break;
      default:
          std::cerr << "PlaceCellsOrderly(): Error. Got an unusable value for par.cell_placement" << '\n';
          exit(1);  
      }
      for(int x = beginx ; x < endx ; x += step ) v_order_x.push_back( x );
      for(int y = beginy ; y < endy ; y += step ) v_order_y.push_back( y );
    }
    else{
      beginx = (par.sizex -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
      endx =   (par.sizex +  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
      beginy = (par.sizey -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
      endy =   (par.sizey +  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
    
      
      switch( par.init_cell_config ){
        case 0: {
          for(int x = beginx ; x < endx ; x += step ) v_order_x.push_back( x );
          for(int y = beginy ; y < endy ; y += step ) v_order_y.push_back( y );
          }
          break;
        case 1: {
          for(int x = endx ; x > beginx ; x -= step ) v_order_x.push_back( x );
          for(int y = beginy ; y < endy ; y += step ) v_order_y.push_back( y );
          }
          break;
        case 2: {
          for(int x = beginx ; x < endx ; x += step ) v_order_x.push_back( x );
          for(int y = endy ; y > beginy ; y -= step ) v_order_y.push_back( y );
          }
          break;
        case 3: {
          for(int x = endx ; x > beginx ; x -= step ) v_order_x.push_back( x );
          for(int y = endy ; y > beginy ; y -= step ) v_order_y.push_back( y );
          }
          break;
        default:
          std::cerr << "PlaceCellsOrderly(): Error. Got an unusable value for par.init_config" << '\n';
          exit(1);
      }
    }

    for(auto x: v_order_x){
      for(auto y: v_order_y){
    // for(int x = beginx ; x < endx ; x += step ){
    //   for(int y = beginy ; y < endy ; y += step ){
        std::cerr << "Cell will be placed at: "<< x<<","<<y << '\n';
      
        count++;
        int this_area=0;
        for(int i=0; i<sqrt(size_cells); i++){
          for(int j=0; j<sqrt(size_cells); j++){
            if(sigma[x+i][y+j]){
              std::cerr << "Sigma = "<< sigma[x+i][y+j] << '\n';
              std::cerr << "Grid point "<< x+i <<","<< y+j <<" is already occupied" << '\n';
              exit(1);
            }
            sigma[x+i][y+j]=count;
            this_area++;
            avrg_area++;
            if(this_area == size_cells) break;
          }
          if(this_area == size_cells) break;
        }
        if(count == n_cells) break;
      }
      if(count == n_cells) break;
    }


    cerr<<"Placed "<<count<<" cells out of "<<n_cells<<" requested; avrg area = "<< avrg_area/(double)count<<endl;
    //exit(1);
    return count;
}

int CellularPotts::Place2Groups(int placement,int size_cells, int groupsize)
{
    int count=0;
    int this_area=0;
    int a_little_bit=2;

    int smaller_dimension=( par.sizex < par.sizey)?par.sizex:par.sizey;
    int sqrt_n_cells = 1+sqrt(groupsize);
    //to avoid having 49 cells when you want 50, I'm rounding sqrt(n_cells) to +1
    if( (  sqrt_n_cells  *(sqrt(size_cells) + a_little_bit))>smaller_dimension ){
      std::cerr << "Place2Groups(): Error. Too many cells or too large size?" << '\n';
      exit(1);
    }

    // int begin = (smaller_dimension-  sqrt(n_cells)*(sqrt(size_cells) + a_little_bit))/2;
    // int end = (smaller_dimension +  sqrt(n_cells)*(sqrt(size_cells) + a_little_bit))/2;
    int beginx1 = a_little_bit;
    int endx1 =   sqrt_n_cells*(sqrt(size_cells) + a_little_bit);
    int beginx2 = (par.sizex -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit));
    int endx2 =   par.sizex-1;

    if (beginx2<endx1) {
      std::cerr<<"Place2Groups(): Error. Groups overlap, exiting."<<endl;
      exit(1);
    }
    int beginy, endy;
    if(placement==1){//the far end of the field
      beginy = (par.sizey -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit));
      endy =   par.sizey-1;
    }else if (placement==2){
      beginy = (par.sizey -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit));
      endy =   par.sizey-1;
      beginx1 = (par.sizex-2*sqrt_n_cells*(sqrt(size_cells)+a_little_bit))/2;
      endx1 =   par.sizex/2 -(sqrt(size_cells)-a_little_bit);
      beginx2 = par.sizex/2;
      endx2 =   (par.sizex+2*sqrt_n_cells*(sqrt(size_cells)+a_little_bit))/2;
    }else{
      beginy = (par.sizey -  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
      endy =   (par.sizey +  sqrt_n_cells*(sqrt(size_cells) + a_little_bit))/2;
    }

    int step = ( sqrt(size_cells) + a_little_bit );
    
    
    //place the cells
    for (int yy=beginy; yy<endy; yy+=step){
      for (int xx=beginx1; xx<endx1; xx+=step){

        this_area=0;
        count++;
        for(int i=0; i<sqrt(size_cells); i++){
          for(int j=0; j<sqrt(size_cells); j++){
            if(sigma[xx+i][yy+j]){
              std::cerr << "Sigma = "<< sigma[xx+i][yy+j] << '\n';
              std::cerr << "Grid point "<< xx+i <<","<< yy+j <<" is already occupied" << '\n';
              exit(1);
            }
            sigma[xx+i][yy+j]=count;
            this_area++;
            if(this_area == size_cells) break;
          }
          if(this_area == size_cells) break;
        }
      }

      for (int xx=beginx2; xx<endx2; xx+=step){

        this_area=0;
        count++;
        for(int i=0; i<sqrt(size_cells); i++){
          for(int j=0; j<sqrt(size_cells); j++){
            if(sigma[xx+i][yy+j]){
              std::cerr << "Sigma = "<< sigma[xx+i][yy+j] << '\n';
              std::cerr << "Grid point "<< xx+i <<","<< yy+j <<" is already occupied" << '\n';
              exit(1);
            }
            sigma[xx+i][yy+j]=count;
            this_area++;
            if(this_area == size_cells) break;
          }
          if(this_area == size_cells) break;
        }
      }
    }


    cerr<<"Placed "<<count<<" cells out of "<<groupsize*2<<" requested;"<<endl;
    //exit(1);
    return count;
}

/**! Fill the plane with initial cells
 \return actual amount of cells (some are not draw due to overlap) */
int CellularPotts::ThrowInCells(int n,int cellsize) {

  //  int gapx=(sizex-nx*cellsize)/(nx+1);
  //int gapy=(sizey-ny*cellsize)/(ny+1);

  int cellnum=1;

  for (int i=0;i<n;i++) {

    // draw a circle at x0, y0
    int x0=RandomNumber(sizex);
    int y0=RandomNumber(sizey);

    bool overlap=false;

    // check overlap
    for (int x=0;x<cellsize;x++)
      for (int y=0;y<cellsize;y++)
	if ( (
	      ( (x-cellsize/2)*(x-cellsize/2)+(y-cellsize/2)*(y-cellsize/2) )<
	      ( (cellsize/2)*(cellsize/2))) &&
	     ( x0+x<sizex && y0+y<sizey ) )
	  if (sigma[x0+x][y0+y]) {
	    overlap=true;
	    break;
	  }

    if (!overlap) {
      for (int x=0;x<cellsize;x++)
	for (int y=0;y<cellsize;y++)
	  if ( (
		( (x-cellsize/2)*(x-cellsize/2)+(y-cellsize/2)*(y-cellsize/2) )<
		( (cellsize/2)*(cellsize/2))) &&
	       ( x0+x<sizex && y0+y<sizey ) )
	    sigma[x0+x][y0+y]=cellnum;

      cellnum++;
    }
  }
  cerr << "[ cellnum = " << cellnum << "]";

  // repair borders
  // fill borders with special border state
  for (int x=0;x<sizex-1;x++) {
    sigma[x][0]=-1;
    sigma[x][sizey-1]=-1;
  }
  for (int y=0;y<sizey-1;y++) {
    sigma[0][y]=-1;
    sigma[sizex-1][y]=-1;
  }

  {for (int x=1;x<sizex-2;x++) {
      sigma[x][1]=0;
      sigma[x][sizey-2]=0;
    }}
  {for (int y=1;y<sizey-2;y++) {
      sigma[1][y]=0;
      sigma[sizex-2][y]=0;
    }}
  return cellnum;
}


int CellularPotts::GrowInCells(int n_cells, int cell_size, double subfield) {


  int sx = (int)((sizex-2)/subfield);
  int sy = (int)((sizey-2)/subfield);

  int offset_x = (sizex-2-sx)/2;
  int offset_y = (sizey-2-sy)/2;

  if (n_cells==1) {
    return GrowInCells(1, cell_size, sizex/2, sizey/2, 0, 0);
  } else {
    return GrowInCells(n_cells, cell_size, sx, sy, offset_x, offset_y);
  }
}

 //WARNING: Cells placed with this function will be much bigger than cell_size: algoritm doesn't work properly
int CellularPotts::GrowInCells(int n_cells, int cell_size, int sx, int sy, int offset_x, int offset_y) {

  // make initial cells using Eden Growth

  int **new_sigma=(int **)malloc(sizex*sizeof(int *));
  if (new_sigma==NULL)
    MemoryWarning();

  new_sigma[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (new_sigma[0]==NULL)
    MemoryWarning();

  for (int i=1;i<sizex;i++)
    new_sigma[i]=new_sigma[i-1]+sizey;

  /* Clear CA plane */
  { for (int i=0;i<sizex*sizey;i++)
     new_sigma[0][i]=0;
  }


  // scatter initial points, or place a cell in the middle
  // if only one cell is desired
  int cellnum=cell->size()-1;

  if (n_cells>1) {



    { for (int i=0;i<n_cells;i++) {

      sigma[RandomNumber(sx)+offset_x][RandomNumber(sy)+offset_y]=++cellnum;

    }}
  } else {
    sigma[sx][sy]=++cellnum;

  }

  // Do Eden growth for a number of time steps
  {for (int i=0;i<cell_size;i++) {
    for (int x=1;x<sizex-1;x++)
      for (int y=1;y<sizey-1;y++) {

	if (sigma[x][y]==0) {
	  // take a random neighbour
	  int xyp=(int)(8*RANDOM()+1);
	  int xp = nx[xyp]+x;
	  int yp = ny[xyp]+y;
	  int kp;
	  //  NB removing this border test yields interesting effects :-)
	  // You get a ragged border, which you may like!
	  if ((kp=sigma[xp][yp])!=-1)
	    if (kp>(cellnum-n_cells))
              new_sigma[x][y]=kp;
            else
	      new_sigma[x][y]=0;
	  else
	    new_sigma[x][y]=0;

	} else {
	  new_sigma[x][y]=sigma[x][y];
	}
      }

    // copy sigma to new_sigma, but do not touch the border!
	  {  for (int x=1;x<sizex-1;x++) {
      for (int y=1;y<sizey-1;y++) {
	sigma[x][y]=new_sigma[x][y];
      }
    }
  }}}
  free(new_sigma[0]);
  free(new_sigma);

  return cellnum;
}

// Hugely expensive function that removes cells forcefully from CA plane
// by looking at the whole plane... the whole frigging plane?
// for a cell that is -like- 5 pixels big when we remove it because it is dying
// 8O
// void CellularPotts::RemoveCell(Cell* thiscell)
// {
//   int sigmaneigh;
//   for (int x=1; x<sizex-1;x++)
//     for (int y=1; y<sizey-1;y++){
//       if (sigma[x][y]==thiscell->Sigma()){
//         sigma[x][y]=0;
//         thiscell->DecrementArea();
//         thiscell->RemoveSiteFromMoments(x,y);
//
//         //remove contact that neighbours have with this pixel
//         for (int k=1; k<=n_nb; k++){
//           int neix=x+nx[k];
//           int neiy=y+ny[k];
//           if(neix<=0 || neix>=sizex-1 || neiy<=0 || neiy>=sizey-1){
//             if( par.periodic_boundaries ){
//               if(neix<=0) neix=sizex-2+neix;
//               if(neix>=sizex-1) neix=neix-sizex+2;
//               if(neiy<=0) neiy=sizey-2+neiy;
//               if(neiy>=sizey-1) neiy=neiy-sizey+2;
//             }else{
//               continue;
//             }
//           }
//           sigmaneigh=sigma[neix][neiy];
//           if( sigmaneigh != thiscell->Sigma()){
//             edgeSetVector.insert({x,y,neix,neiy});
//             edgeSetVector.insert({neix,neiy,x,y});
//             if( sigmaneigh){
//               (*cell)[sigmaneigh].updateNeighbourBoundary(thiscell->Sigma(),-1);
//               (*cell)[sigmaneigh].updateNeighbourBoundary(0,1);
//             }
//           }
//           else{
//             edgeSetVector.erase({x,y,neix,neiy});
//             edgeSetVector.erase({neix,neiy,x,y});
//           }
//         }
//       }
//     }
//
// }

int CellularPotts::FancyloopX(int loopdepth, int meanx, int meany, int thissig, bool above){
  int removed=0;
  bool loop=true;
  int py, px;
  int sigmaneigh;
  py=meany-((int)above*2-1)*loopdepth;
  if(py<=0||py>=sizey-1){
    if( par.periodic_boundaries ){
      if(py<=0)
        py=sizey-2+py;
      else if (py>=sizey-1)
        py=py-sizey+2;
    }else{
      loop=false;
    }
  }

  if(loop)
    for(int x=meanx-loopdepth; x<=meanx+loopdepth;x++){
      px=x;
      if(x<=0 || x>=sizex-1){
        if( par.periodic_boundaries ){
          if(x<=0)
            px=sizex-2+x;
          else if (x>=sizex-1)
            px=x-sizex+2;
        }else{
          continue;
        }
      }
//         cerr<<"removing?"<<endl;
        if (sigma[px][py]==thissig){
//           cerr<<"going on with removing X"<<endl;
          sigma[px][py]=0;
          removed++;
          for (int k=1; k<=n_nb; k++){
            int neix=px+nx[k];
            int neiy=py+ny[k];
            if(neix<=0 || neix>=sizex-1 || neiy<=0 || neiy>=sizey-1){
              if( par.periodic_boundaries ){
                if(neix<=0) neix=sizex-2+neix;
                if(neix>=sizex-1) neix=neix-sizex+2;
                if(neiy<=0) neiy=sizey-2+neiy;
                if(neiy>=sizey-1) neiy=neiy-sizey+2;
              }else{
                continue;
              }
            }
            sigmaneigh=sigma[neix][neiy];
            if( sigmaneigh == 0){
              edgeSetVector.erase({px,py,neix,neiy});
              edgeSetVector.erase({neix,neiy,px,py});
            }
        }
      }
    }

  return removed;
}

int CellularPotts::FancyloopY(int loopdepth, int meanx, int meany, int thissig, bool left){
  int removed=0;
  bool loop=true;
  int py, px;
  int sigmaneigh;

  px=meanx-((int)left*2-1)*loopdepth;
  if(px<=0||px>=sizex-1){
    if( par.periodic_boundaries ){
      if(px<=0)
        px=sizex-2+px;
      else if (px>=sizex-1)
        px=px-sizex+2;
    }else{
      loop=false;
    }
  }

  if(loop)
    for(int y=meany-loopdepth+1; y<=meany+loopdepth-1;y++){
      py=y;
      if(y<=0 || y>=sizey-1){
        if( par.periodic_boundaries ){
          if(y<=0)
            py=sizey-2+y;
          else if (y>=sizey-1)
            py=y-sizey+2;
        }else{
          continue;
        }
      }

        if (sigma[px][py]==thissig){
          sigma[px][py]=0;
          removed++;
          for (int k=1; k<=n_nb; k++){
            int neix=px+nx[k];
            int neiy=py+ny[k];
            if(neix<=0 || neix>=sizex-1 || neiy<=0 || neiy>=sizey-1){
              if( par.periodic_boundaries ){
                if(neix<=0) neix=sizex-2+neix;
                if(neix>=sizex-1) neix=neix-sizex+2;
                if(neiy<=0) neiy=sizey-2+neiy;
                if(neiy>=sizey-1) neiy=neiy-sizey+2;
              }else{
                continue;
              }
            }
            sigmaneigh=sigma[neix][neiy];
            if( sigmaneigh == 0){
              edgeSetVector.erase({px,py,neix,neiy});
              edgeSetVector.erase({neix,neiy,px,py});
            }
        }
        }
      }

    return removed;
}


//Hopefully a little less expensive function than the one above
void CellularPotts::RemoveCell(Cell* thiscell,int min_area, int meanx, int meany)
{
  int sigmaneigh,thissig,thisarea;
  int nx, ny;
  int countpix=0; //to check if we removed all pixels
  bool loop=true;
  int loopdepth=1;
  thissig=thiscell->Sigma();
  thisarea=thiscell->Area();

//   cerr<<"meanx: "<<meanx<<" meany: "<<meany<<endl;
  if(!thisarea){
	cerr<<"CA-RemoveCell warning: attempting to remove dead cell. "<<endl;
  }

  if(sigma[meanx][meany]==thissig){
    sigma[meanx][meany]=0;
    countpix++;
//     cerr<<"Removed one, countpix ="<<countpix<<endl;
  }

  while(true){
    //top row
    //Here we go through progressively larger squares (only its boundary)
    //first through the top row, then the bottom row, then left and right sides
    //which are smaller.
    countpix+=FancyloopX(loopdepth,meanx,meany,thissig,true);
    if(countpix==thisarea) break;
    countpix+=FancyloopX(loopdepth,meanx,meany,thissig,false);
    if(countpix==thisarea) break;
    countpix+=FancyloopY(loopdepth,meanx,meany,thissig,true);
    if(countpix==thisarea) break;
    countpix+=FancyloopY(loopdepth,meanx,meany,thissig,false);
    if(countpix==thisarea) break;

//     cerr<<"Area: "<<thisarea<<" loopdepth: "<<loopdepth<<" countpix: "<<countpix<<endl;
//     if(loopdepth>10) exit(1);

    loopdepth++;

  }
   
   //take care that area is set to 0
    thiscell->DecrementAreaBy(thisarea);
    //removed the tracking of the moments. Check if this gives major problems

	
  //deal with the cells that have this cell as a neighbour
  for(auto &neigh:thiscell->neighbours){
    int signeigh=neigh.first;
    int blength=neigh.second.first;
  
    if(signeigh){
      (*cell)[signeigh].setNeighbour(thissig,0,0);
      (*cell)[signeigh].updateNeighbourBoundary(0,blength);
    }
    
  }
  
  thiscell->clearNeighbours();

//   exit(1);

}

// Predicate returns true when connectivity is locally preserved
// if the value of the central site would be changed
// note that at this point we already know that focal point is neighbouring a point with different sigma
bool CellularPotts::ConnectivityPreservedP(int x, int y) {

  // Use local nx and ny in a cyclic order (starts at upper left corner)
  // zeroth site is ignored and first site is repeated, for easier looping - that's why the array is 10 long
  // even though there are only 8 neighbours
  // (see below: int s_next_nb=sigma[x+cyc_nx[i+1]][y+cyc_ny[i+1]];)
  const int cyc_nx[10] = {-1, -1, 0, 1, 1, 1, 0, -1, -1, -1 };
  const int cyc_ny[10] = {0, -1,-1,-1, 0, 1, 1,  1,  0, -1 };

  int sxy=sigma[x][y]; // the central site
  if (sxy==0) return true;

  int n_borders=0; // to count the amount of sites in state sxy bordering a site !=sxy

  static int stack[8]; // stack to count number of different surrounding cells
  int stackp=-1;
  bool one_of_neighbors_medium=false;

  for (int i=1;i<=8;i++) {

    int s_nb=sigma[x+cyc_nx[i]][y+cyc_ny[i]]; //this is sigma of neighbour
    int s_next_nb=sigma[x+cyc_nx[i+1]][y+cyc_ny[i+1]]; //this is sigma of next neighbour on the list
    //if at least one of them == focal sigma, and they are different from each other
    // i.e. if one of them is same as sxy and other not (but whihc one does not matter)
    if ((s_nb==sxy || s_next_nb==sxy) && (s_nb!=s_next_nb)) {

      // check whether s_nb is adjacent to non-identical site,
      n_borders++; // count it
    }
    int j;
    bool on_stack_p=false;

    // we need the next heuristic to prevent stalling at
    // cell-cell borders
    // do not enforce constraint at two cell interface(no medium)
    if (s_nb) {
      for (j=stackp;j>=0;j--) {
        if (s_nb==stack[j]) {
          on_stack_p=true;
          break;
        }
      }
      if (!on_stack_p) {
        if(stackp>6) {
          cerr << "Stack overflow, stackp=" << stackp << "\n";
        }
        stack[++stackp]=s_nb;
      }
    } else {
      one_of_neighbors_medium=true;
    }
  }

  // number of different neighbours is stackp+1;
  if (n_borders>2 && ( (stackp+1)>2 || one_of_neighbors_medium) ) {
    return false;
  }
  else
    return true;

}


double CellularPotts::CellDensity(void) const {

  // return the density of cells
  int sum=0;
  for (int i=0;i<sizex*sizey;i++) {
    if (sigma[0][i]) {
      sum++;
    }
  }
  return (double)sum/(double)(sizex*sizey);

}

double CellularPotts::MeanCellArea(void) const {

  int sum_area=0, n=0;
  double sum_length=0.;
  vector<Cell>::iterator c=cell->begin(); ++c;

  for (;
	c!=cell->end();
	c++) {

    sum_area+=c->Area();
    sum_length+=c->Length();
    n++;
  }

  cerr << "Mean cell length is " << sum_length/((double)n) << endl;
  return (double)sum_area/(double)n;
}

void CellularPotts::ResetTargetLengths(void)  {
   vector<Cell>::iterator c=cell->begin(); ++c;

   for (;
        c!=cell->end();
        c++) {

     c->SetTargetLength(par.target_length);

}

}

void CellularPotts::SetRandomTypes(void) {

  // each cell gets a random type 1..maxtau

  vector<Cell>::iterator c=cell->begin(); ++c;

  for ( ; c!=cell->end(); c++) {

    int celltype = RandomNumber(par.maxtau);
    c->setTau(celltype);

    c->SetCellTypeProperties();  //depending on model, set different growth rates
  }

}

//cells are changed to a specific tau with some probability
void CellularPotts::SetRandomTypes(int whichtau, double probtau) {
  vector<Cell>::iterator c=cell->begin(); ++c;

  for ( ; c!=cell->end(); c++) {
    if( RANDOM()<probtau ){
      c->setTau(whichtau);
      //c->SetCellTypeProperties();   // set different growth rates
    }
  }

}

void CellularPotts::GrowAndDivideCells(int growth_rate) {

  vector<Cell>::iterator c=cell->begin(); ++c;
  vector<bool> which_cells(cell->size());

  for (;
       c!=cell->end();
       c++) {

    // only tumor cells grow and divide
    if (c->getTau()==2) {

      c->SetTargetArea(c->TargetArea()+growth_rate);

      if (c->Area()>par.target_area) {
	which_cells[c->Sigma()]=true;
      } else {
	which_cells[c->Sigma()]=false;
      }

      if (c->chem[1]<0.9) { //arbitrary oxygen threshold for the moment
	c->setTau(3);
      }
    } else {
      which_cells[c->Sigma()]=false;
    }

  }

  DivideCells(which_cells);

}

double CellularPotts::DrawConvexHull(Graphics *g, int color) {

  // Draw the convex hull of the cells
  // using Andrew's Monotone Chain Algorithm (see hull.cpp)

  // Step 1. Prepare data for 2D hull code

  // count number of points to determine size of array
  int np=0;
  for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	np++;
      }
    }

  Point *p=new Point[np];

  int pc=0;
  for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	p[pc++]=Point(x,y);
      }
    }

  // Step 2: call 2D Hull code
  Point *hull=new Point[np];
  int nph=chainHull_2D(p,np,hull);

  // Step 3: draw it
  for (int i=0;i<nph-1;i++) {
    g->Line(2*hull[i].x,2*hull[i].y,2*hull[i+1].x,2*hull[i+1].y, color);
  }


  // Step 4: calculate area of convex hull
  double hull_area=0.;
  for (int i=0;i<nph-1;i++) {
    hull_area+=hull[i].x*hull[i+1].y-hull[i+1].x*hull[i].y;
  }
  hull_area/=2.;

  //cerr << "Area = " << hull_area << "\n";

  delete[] p;
  delete[] hull;

  return hull_area;

}

double CellularPotts::Compactness(double *res_compactness, double *res_area, double *res_cell_area) {

  // Calculate compactness using the convex hull of the cells
  // We use Andrew's Monotone Chain Algorithm (see hull.cpp)

  // Step 1. Prepare data for 2D hull code

  // count number of points to determine size of array
  int np=0;
  for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	np++;
      }
    }

  Point *p=new Point[np];

  int pc=0;
  for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++) {
      if (sigma[x][y]) {
	p[pc++]=Point(x,y);
      }
    }

  // Step 2: call 2D Hull code
  Point *hull=new Point[np];
  int nph=chainHull_2D(p,np,hull);

  //// Step 3: draw it
  //for (int i=0;i<nph-1;i++) {
  //  g->Line(2*hull[i].x,2*hull[i].y,2*hull[i+1].x,2*hull[i+1].y, color);
  //}


  // Step 3: calculate area of convex hull
  double hull_area=0.;
  for (int i=0;i<nph-1;i++) {
    hull_area+=hull[i].x*hull[i+1].y-hull[i+1].x*hull[i].y;
  }
  hull_area/=2.;

  // Step 4: calculate total cell area
  double cell_area=0;

  vector<Cell>::const_iterator c;

  for ( (c=cell->begin(),c++);
       c!=cell->end();
       c++) {
    cell_area+=c->Area();
  }

  delete[] p;
  delete[] hull;


  // put intermediate results into optional pointers
  if (res_compactness) {
    *res_compactness = cell_area/hull_area;
  }
  if (res_area) {
    *res_area = hull_area;
  }
  if (res_cell_area) {
    *res_cell_area = cell_area;
  }

  // return compactness
  return cell_area/hull_area;

}
