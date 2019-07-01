#include <iostream>
#include <fstream>
#include "aaEnvironment.hpp"
#include "affineArithmeticClass.hpp"
#include <boost/random.hpp>
#include <ctime>
#include <cstdlib>
#include <cmath>



void filterExampleTest(){
  AAEnvironment env(1);
  AffineArithmeticClass  xn3(env),  yn(env), yn1(env), yn2(env);
  int nTrials = 100;

  // Modified to a less uniform distribution.
  AffineArithmeticClass xn1 = randomVariable(env, i_double(-1.0, 1.0), i_double(0.0), i_double(0.33333));
  AffineArithmeticClass xn2 = randomVariable(env, i_double(-1.0, 1.0), i_double(0.0), i_double(0.33333));
 

  for (int i =0; i < nTrials; ++i){
    /*!npk xn between -0.2  and 0.2  mean 0.0 */
    AffineArithmeticClass xn = randomVariable(env, i_double(-1.0, 1.0), i_double(0.0), i_double(0.3333));

      /* -- Cut and paster from Olivier's example: filter_roundingerrors.c --*/
    //yn = (0.7*xn - 1.3*xn1 + 1.1*xn2 + 1.4*yn1 - 0.45*yn2);
    yn = (i_double(0.7)*xn) - (i_double(1.3)*xn1) + (i_double(1.1)*xn2) + (i_double(1.4) *yn1 )- (i_double(0.7)* yn2) ;
    //yn = i_double(1.7196456)*yn - i_double(0.75435)*yn1 + i_double(0.8684)*xn - i_double(1.7368)*xn1 + i_double(0.8684)*xn2;
    yn.normalizeForm_assign();
    xn2=xn1;
    xn1=xn;
    yn2=yn1;
    yn1=yn;
    cout << "Iteration # " << i << " done ! " << endl;
  }

  cout << yn << endl;
  cout << "Range: " << yn.range() << endl;
  cout << "Expectation: " << yn.expectation() << endl;
  
  cout << "--- Analysis results for yn --- " << endl;
  yn.performAnalysisForCMI(cout);
  
  
}
void filterExampleTestWithRoundoff(){
  AAEnvironment env(1);
  AffineArithmeticClass  xn3(env),  yn(env), yn1(env), yn2(env);
  int nTrials = 500;

  AffineArithmeticClass xn1 = randomVariable(env, i_double(-1.0, 1.0), 0.0, 2/3);
  AffineArithmeticClass xn2 = randomVariable(env, i_double(-1.0, 1.0), 0.0, 2/3);
  

  for (int i =0; i < nTrials; ++i){
    /*!npk xn between -1.0  and 1.0  mean 0.2 */
    AffineArithmeticClass xn = randomVariable(env, i_double(-1.0, 1.0), 0.2, 2.0/3.0);
    AffineArithmeticClass en = randomVariable(env, i_double(-1E-06,1E-06), 0.0, 1E-10);
    AffineArithmeticClass en1 = randomVariable(env, i_double(-1E-06,1E-06), 0.0, 1E-10);
    AffineArithmeticClass en2 = randomVariable(env, i_double(-1E-06,1E-06), 0.0, 1E-10);
    AffineArithmeticClass en3 = randomVariable(env, i_double(-1E-06,1E-06), 0.0, 1E-10);
    AffineArithmeticClass en4 = randomVariable(env, i_double(-1E-06,1E-06), 0.0, 1E-10);
    AffineArithmeticClass en5 = randomVariable(env, i_double(-1E-06,1E-06), 0.0, 1E-10);
    AffineArithmeticClass en6 = randomVariable(env, i_double(-1E-06,1E-06), 0.0, 1E-10);
    AffineArithmeticClass en7 = randomVariable(env, i_double(-1E-06,1E-06), 0.0, 1E-10);
    AffineArithmeticClass en8 = randomVariable(env, i_double(-1E-06,1E-06), 0.0, 1E-10);
    /* -- Cut and paster from Olivier's example: filter_roundingerrors.c --*/
    // Do not forget to bracket this or C++ will give an error.
    yn = (0.7*xn - 1.3*xn1 + 1.1*xn2 + 1.4*yn1 - 0.45*yn2);
    // Do not forget to bracket the RHS or else C++ compiler will give an error.
    
    yn = ( yn + 0.7*xn*en + 1.3*xn1*en1 + (0.7*xn - 1.3*xn1)*en3 + 1.1*xn2*en2 + 
      ( 0.7*xn - 1.3*xn1 + 1.1*xn2)*en4 + 1.4*yn1*en5 + 
      (0.7*xn - 1.3*xn1 + 1.1*xn2 + 1.4*yn1)*en6 + 0.45*yn2*en7 +
	   (0.7*xn - 1.3*xn1 + 1.1*xn2 + 1.4*yn1 - 0.45*yn2)*en8  );
    xn2=xn1;
    xn1=xn;
    yn2=yn1;
    yn1=yn;
    cout << "Iteration # " << i << " done ! " << endl;
  }

  cout << yn << endl;
  cout << "Range: " << yn.range() << endl;
  cout << "Expectation: " << yn.expectation() << endl;
  
  cout << "--- Analysis results for yn --- " << endl;
  yn.performAnalysisForCMI(cout);
  
  
}

int main(){
  filterExampleTest();
}
