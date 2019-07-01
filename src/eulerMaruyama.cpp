
#include <iostream>
#include <fstream>
#include "aaEnvironment.hpp"
#include "affineArithmeticClass.hpp"
#include <boost/random.hpp>
#include <ctime>
#include <cstdlib>
#include <cmath>





void eulerMaruyama() {
  // double tBegin, tEnd, dt, IC, theta, mu, sigma, sqrtdt, ysn;
  int N, i;
  
  AAEnvironment env(1);
  i_double dt(.0001); /* ten times as much as in the original Python program */
  
  

  i_double theta(1.0);
  i_double mu(1.2);
  i_double sigma(0.3);
  i_double(sqrtdt) = sqrt(dt);

  AffineArithmeticClass yn(env);
  AffineArithmeticClass ysn(env);
  i = 1;
  N = 1000;

  while (i<=N) {
    /* pb en fait rg a un support non borne, comment fait on? */
    /*!npk rg between -100000.0 and 100000.0 mean 0.0 */
    /*!npk cov rg rg 1.0 1.0 */
    AffineArithmeticClass rg = randomVariable(env, i_double(-100.0, 100.0), i_double(0.0), i_double(1.0));
    ysn = yn + dt*(theta*(mu-yn)) + sigma*sqrtdt*rg;
    yn = ysn;
    i = i+1;
   
  }

  cout << yn << endl;
  cout << "Range: " << yn.range() << endl;
  cout << "E: "<< yn.expectation() << endl;
  yn.performAnalysisForCMI(cout);
}

int main(){
  eulerMaruyama();
}

 
