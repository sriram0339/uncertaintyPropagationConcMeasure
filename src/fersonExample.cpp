
#include <iostream>
#include <fstream>
#include "aaEnvironment.hpp"
#include "affineArithmeticClass.hpp"
#include <boost/random.hpp>
#include <ctime>
#include <cstdlib>
#include <cmath>



void fersonExample2(){
  AAEnvironment env(1);
  /*!npk t1 between 2.99 and 3.01  mean 3. */
  /*!npk cov t1 t1 0.002 0.0038 */

  // Let us go lower level to create these variables
  
  AffineArithmeticClass t1 = randomVariable(env, i_double(-0.01,.01), i_double(0.0), i_double(4E-06, 16E-06));

  t1 = i_double(3.0) + t1;
  //int t1ID = env.makeFreshVariable(i_double(2.99,3.01), i_double(3.0), i_double(9.002, 9.0038));
 
  
  
  
  /*!npk t2 between 0.99 and 1.01  mean 1.0 */
  /*!npk cov t2 t2 0.002 0.0038 */
  AffineArithmeticClass t2 = randomVariable(env, i_double(-0.01,0.01), i_double(0.0), i_double(4E-06, 16E-06));
  t2 = i_double(1.0) + t2;
  i_double I(1.1444, 1.1477);
  AffineArithmeticClass interv = randomVariable(env, I, I, square(I));
  AffineArithmeticClass t1Sq = square(t1);
  AffineArithmeticClass t2Sq = square(t2);
  AffineArithmeticClass t1_3 = t1Sq * t1;
  AffineArithmeticClass t2_3 = t2Sq * t2;
  AffineArithmeticClass t1_4 = square(t1Sq);
  AffineArithmeticClass t2_4 = square(t2Sq);
  AffineArithmeticClass t2_5 = t2_4*t2;
  AffineArithmeticClass t1_5 = t1_4*t1;

  // Mistake # 1 term: - i_double(31.9645099690168)*t1_3 should be t1_4?

  AffineArithmeticClass res = ( i_double(1.09467683743172) *t1_5 + //??? 
				i_double( 15.3223572746484) *t1_4*t2 // OK
				+i_double(86.0297202975092)*t1_3 *t2Sq // OK
				+i_double(390.260443722082)*t1*t2_4 // OK
				+i_double(218.60704232612)*t2_5 // OK
				- i_double(31.9645099690168)*t1_4 // ???
				- i_double(361.491077961158)*t1_3*t2 // OK
				- i_double(800.011819233185)*t1Sq*t2Sq
				- i_double(1334.22983617068)*t1*t2_3 //
				- i_double(2257.7540329627)*t2_4 + // OK
				i_double(375.296580066101)*t1_3 // OK
				+i_double(2465.09948122341)*t1Sq*t2 // OK
				+i_double(4103.19455646941)*t1*t2Sq // OK
				+i_double(6134.72163986329)*t2_3 //OK
				-i_double(1960.94158113753)*t1Sq // OK 
				-i_double(7497.63254682701)*t1*t2 // OK
				-i_double(9497.30109206786)*t2Sq // OK
				+i_double(4772.67824380499)*t1 // OK
				+i_double(9816.86840978377)*t2 // OK
				-i_double(4668.81238573154) // OK
				+interv ); // OK

  
  cout << "result = " << res << endl;
  cout << "range = " << res.range() << endl;
  cout << "expectation= " << res.expectation() << endl;
  
  res.performAnalysisForCMI(cout);
}


void fersonExample(){
  AAEnvironment env(1);
  
   AffineArithmeticClass t1 = randomVariable(env, i_double(-0.01,0.01), i_double(0.0), i_double(4E-06, 16E-06));

  t1 = i_double(3.0) + t1;
  //int t1ID = env.makeFreshVariable(i_double(2.99,3.01), i_double(3.0), i_double(9.002, 9.0038));
 
  
  
  
  /*!npk t2 between 0.99 and 1.01  mean 1.0 */
  /*!npk cov t2 t2 0.002 0.0038 */
  AffineArithmeticClass t2 = randomVariable(env, i_double(-0.01,0.01), i_double(0.0), i_double(4E-06, 16E-06));
  t2 = i_double(1.0) + t2;

  cout << t1 << endl;
  cout << t2 << endl;

  /*!npk t1 between 2.99 and 3.01  mean 3. */
  /*!npk cov t1 t1 0.002 0.0038 */

  // Let us go lower level to create these variables
  
  //AffineArithmeticClass t1 = randomVariable(env, i_double(2.99,3.01), i_double(3.0), i_double(9.002, 9.0038));

  //int t1ID = env.makeFreshVariable(i_double(2.99,3.01), i_double(3.0), i_double(9.002, 9.0038));
 
  
  
  
  /*!npk t2 between 0.99 and 1.01  mean 1.0 */
  /*!npk cov t2 t2 0.002 0.0038 */
  // AffineArithmeticClass t2 = randomVariable(env, i_double(0.99,1.01), i_double(1.0), i_double(1.002, 1.0038));

  
  
  /*!npk interv between 1.1444264300811058 and 1.1477539164945061 */
  i_double I(1.1444264, 1.1477539);
  AffineArithmeticClass interv = randomVariable(env, I, I, square(I));
  
  AffineArithmeticClass res =  ((((1.09467683743172*t1+(15.3223572746484*t2-31.9645099690168))*t1+((86.0297202975092*t2-361.491077961158)*t2+375.296580066101))*t1+((-800.011819233185*t2+2465.09948122341)*t2-1960.94158113753))*t1+((((390.260443722082*t2-1334.22983617068)*t2+4103.19455646941)*t2-7497.63254682701)*t2+4772.67824380499))*t1+((((218.60704232612*t2-2257.7540329627)*t2+6134.72163986329)*t2-9497.30109206786)*t2+9816.86840978377)*t2-4668.81238573154+interv;


  
  cout << "result = " << res << endl;
  cout << "range = " << res.range() << endl;
  cout << "expectation= " << res.expectation() << endl;
  
  res.performAnalysisForCMI(cout);
  

}

int main(){
  fersonExample();
}
