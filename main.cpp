/**
 * ################################# *
 *                                   *
 *            DRAGONCELLO            *
 *                                   *
 *  -------------------------------- *
 *  Development hystory:             *
 *                                   *
 *  Daniele Gaggero (2015)           *
 *  Carmelo Evoli (2015)             *
 *  Silvio Sergio Cerri (2016)       *
 *                                   *
 *   - - - - - - - - - - - - - - -   *
 *  Reference paper:                 *
 *                                   *
 *  Cerri et al., JCAP 10:019 (2017) *
 *                                   *
 * ################################# *
 *
 */

#include <iostream>
#include <ctime>
#include <math.h>

#include "dragoncello.h"
#include "constants.h"

using namespace std;

int main(int argc, char** argv) {
  
  cout << "Welcome to DRAGONCELLO!" << endl;
  
  DRAGONCELLO* draghetto = new DRAGONCELLO(Zparticle, Aparticle);
  
  draghetto->initGrids();
   	  
  draghetto->initBfield();
  
  draghetto->initAnisoSpatialDiffusion();

  draghetto->initDriftLikeVelocities();

  draghetto->initSource();
  
  //draghetto->initEloss(); //IC-losses for leptons only (un-comment if needed)
  
  draghetto->initCN2Daniso();
  
  draghetto->Run2D();
    
  delete draghetto;
  
  return 0;
}
