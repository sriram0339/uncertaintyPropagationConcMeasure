#include <iostream>
#include <fstream>
#include "aaEnvironment.hpp"
#include "affineArithmeticClass.hpp"
#include <boost/random.hpp>
#include <ctime>
#include <cstdlib>
#include <cmath>

/*-- Robotics running example #1 for paper --*/

void testRoboticsExample(){
  // This has been implemented elsewhere
  // Thanks!
  const double angles[10]= {10,60,110,160,140,100,60,20,10,0.0};
  int nAngles = 10;
  AAEnvironment aa(2);
  i_double pi(3.1415,3.1416);
  
  // x := TruncGaussian(0, 0.05, -0.2, 0.2 )
  i_double i0(-0.1,0.1), i1(0.0), i2(0.00248,0.00252);
  int y1 = aa.makeFreshVariable(i0,i1,i2);
  AffineArithmeticClass x (aa);
  x.setCoefficient(y1,i_double(1.0));

  // y := TruncGaussian(0, 0.1, -0.5, 0.5);
  i_double i3(-0.5,0.5), i4(0.0), i5(0.0098,0.012);
  int y2 = aa.makeFreshVariable(i3,i4,i5);
  AffineArithmeticClass y (aa);
  y.setCoefficient(y2, i_double(1.0));
  int nTrials=100;

  //for theta in angles
  for (int kk =0 ; kk < nTrials; ++kk){
    for (int i = 0; i < nAngles; ++i){
    
      // d := Uniform(0.98, 1.02)
      i_double i6(0.98,1.02), i7(1.0), i8(1.001,1.002);
      int yy3 = aa.makeFreshVariable(i6,i7,i8);
      AffineArithmeticClass d(aa);
      d.setCoefficient(yy3,i_double(1.0));
      // deg2rad(theta)
      i_double theta(angles[i]);
      theta = theta * pi/i_double(180.0);

      // 1 + TruncGaussian(0,0.01,-0.05,0.05)
      i_double i9(-0.05,0.05), i10(0), i11(0.0001,0.00011);
      int yy4 = aa.makeFreshVariable(i9,i10,i11);
      AffineArithmeticClass tmp(aa);
      tmp.setConstant(i_double(1.0));
      tmp.setCoefficient(yy4,i_double(1.0));

      // t := deg2rad(theta) * (1 + Trunc...)
      tmp.scale_assign(theta);

      // x := x + d * cos(t)

      AffineArithmeticClass tmp1(aa);
      tmp1.cosine_assign(tmp);
      tmp1.multiply_assign(d);
      x.add_assign(tmp1);

      // y := y + d * sin(t)
    
      AffineArithmeticClass tmp2(aa);
      tmp2.sine_assign(tmp);
      tmp2.multiply_assign(d);
      y.add_assign(tmp2);
    }
  }
  cout << "done " << endl;
  cout << aa << endl;
  cout << "x= " << x << endl;
  cout << "Range: " << x.range() << endl;
  cout << "Expectation: " << x.expectation() << endl;
  // AffineArithmeticClass tmp3(x);
  // tmp3.multiply_assign(x);
  // cout << "Second Moment: " << tmp3.expectation() << endl;
  cout << "Denominator Bound for CMI: "<< x.getDenominatorBoundForConcMeasure() << endl;
  
  cout << "y= " << y << endl;
  cout << "Range: " << y.range() << endl;
  cout << "Expectation: " << y.expectation() << endl;
  // AffineArithmeticClass tmp4(y);
  //tmp4.multiply_assign(y);
  //cout << "Second Moment: " << tmp4.expectation() << endl;
  
  cout << "Denominator Bound for CMI: "<< y.getDenominatorBoundForConcMeasure() << endl;

  // Now print the graph for x
  std::ofstream outGraphFile;
  outGraphFile.open("xGraph.dot");
  std::set<int> relVars;
  x.collectRelevantNoiseSymbols(relVars);
  Graph gRet;
  aa.computeDependencyGraph(gRet, relVars);
  outputGraphToDot(gRet,outGraphFile);
  int maxDeg = findMaximumDegree(gRet);
  cout << "Chromatic number upper bound for x = " << (1+maxDeg) << endl;
  outGraphFile.close();
  i_double var;
  double denom;
  double M;
  x.clusterNoiseSymbolsForChernoffHoeffding(gRet, denom, var, M);
  cout << " Denominator bound for Chernoff Hoeffding for x is " << denom << endl;
  cout << " Variance for x is " << var << endl;
  cout << " Max deviation " << M << endl;

  cout << "--- Analysis results for x --- " << endl;
  x.performAnalysisForCMI(cout);

  cout << "--- Analysis results for y --- " << endl;
  y.performAnalysisForCMI(cout);
  
}

/*-- End Robotics running exple --*/

void monteCarloSineCosineTest2(){
 
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<> dist1(-0.1,0.1);
  boost::random::uniform_real_distribution<> dist2(0.1,0.15);
  double totalS = 0.0;
  double varS = 0.0;
  double maxS = -1.0;
  double minS = 1.0;
  double totalC = 0.0;
  double varC = 0.0;
  double maxC = maxS;
  double minC = minS;
  double totalT = 0.0;
  double varT = 0.0;
  double maxT = -10000000.0;
  double minT = 1000000000.0;
  int N = 100000;
  double totalI = 0.0;
  double varI = 0.0;
  double minI = 100000000.0;
  double maxI = -1000000000.0;
  
  for (int i=0; i < N; ++i){
    double y = dist1(gen);
    double z = dist2(gen);
    double v3 = sin( z* (1.0 + y));
    totalS=totalS+v3;
    varS = varS + v3*v3;
    if (v3 < minS)
      minS = v3;
    if (v3 > maxS)
      maxS = v3;

    double v4 = cos(z * (1.0+y));
    totalC = totalC + v4;
    varC = varC + v4*v4;
     if (v4 < minC)
      minC = v4;
    if (v4 > maxC)
      maxC = v4;

    double v5 = tan(0.3 * z * (1.0 + y));
    totalT = totalT + v5;
    varT = varT + v5 * v5;
    if (v5 < minT)
      minT = v5;
    if (v5 > maxT)
      maxT = v5;

    double v6 = 1.0/(0.1* z * ( 1.0 + y));
    totalI = totalI + v6;
    varI = varI + v6 * v6;
    if (v6 < minI)
      minI = v6;
    if (v6 > maxI)
      maxI = v6;
    
  }

  varS = varS/N - (totalS/N)*(totalS/N);
  varC = varC/N - (totalC/N)*(totalC/N);
  varT = varT/N - (totalT/N) * (totalT/N);
  varI = varI/N - (totalI/N) * (totalI/N);
  cout << "MC100000 range (sine): " << minS << "," << maxS << endl;
  cout << "MC100000 expect (sine): " << totalS/N << endl;
  cout << "MC100000 variance (sine): " << varS << endl<<endl;
  cout << "MC100000 range (cosine): " << minC << "," << maxC << endl;
  cout << "MC100000 expect (cosine): " << totalC/N << endl;
  cout << "MC100000 variance (cosine): " << varC << endl<<endl;
  cout << "MC100000 range (tan): " << minT << "," << maxT << endl;
  cout << "MC100000 expect (tan): " << totalT/N << endl;
  cout << "MC100000 variance (tan): " << varT << endl<<endl;
  cout << "MC100000 range (inverse): " << minI << "," << maxI << endl;
  cout << "MC100000 expect (inverse): " << totalI/N << endl;
  cout << "MC100000 variance (inverse): " << varI << endl<<endl;
  
  
}

void test2(){
  // Test for sine and cosine.
  // Generate a noise symbol y modeling a 10% relative error.
  AAEnvironment aa(1);
  // range(y) = [0.1,0.15], expectation = [0.12,0.13], variance = [0.015, 0.016]
  i_double I1(-0.1, 0.1), I2(0), I3(0.003,0.004);
  int y = aa.makeFreshVariable(I1, I2, I3);
  AffineArithmeticClass f1(aa);
  f1.setConstant(i_double(1.0));
  f1.setCoefficient(y, i_double(1.0));

  // Make an angle that is between 0.1, 0.15 radians
  i_double I4(0.1,0.15), I5(0.124,0.126), I6(0.015,0.016);
  int z = aa.makeFreshVariable(I4,I5,I6);
  AffineArithmeticClass f2(aa);
  //Make an affine form of that angle
  f2.setCoefficient(z,i_double(1.0));

  AffineArithmeticClass f3(f2);
  f3.multiply_assign(f1);

  AffineArithmeticClass f4(aa);
  f4.sine_assign(f3);
  
  cout << "sin( z * (1+y) ) = "<< f4 << endl;
  cout << "Range: " << f4.range() << endl;
  cout << "Expectation: "<< f4.expectation() << endl;
  cout << "Variance : " << f4.variance_raw() << endl<<endl;
  
  AffineArithmeticClass f5(aa);
  f5.cosine_assign(f3);
  
  cout << "cos( z * (1+y) ) = "<< f5 << endl;
  cout << "Range: " << f5.range() << endl;
  cout << "Expectation: "<< f5.expectation() << endl;
  cout << "Variance : " << f5.variance_raw() << endl<<endl;

  

  AffineArithmeticClass f7(aa);
  AffineArithmeticClass f8(f3);
  f8.scale_assign(i_double(0.1));
  f7.inverse_assign(f8);
  
  cout << "1/( 0.1*z * (1+y) ) = "<< f7 << endl;
  cout << "Range: " << f7.range() << endl;
  cout << "Expectation: "<< f7.expectation() << endl;
  cout << "Variance : " << f7.variance_raw() << endl <<endl;
  

  
  AffineArithmeticClass f6(aa);
  f3.scale_assign(i_double(0.3));
  f6.tan_assign(f3);
  
  cout << "tan( 0.3* z * (1+y) ) = "<< f6 << endl;
  cout << "Range: " << f6.range() << endl;
  cout << "Expectation: "<< f6.expectation() << endl;
  cout << "Variance : " << f6.variance_raw() << endl << endl;

}




void test1() {
  // test1:
  //     define four independent noise symbols y1, y2, y3, y4
  //     
  //     define two affine forms
  //                   f1 = [-1, 2] + [2, 3] * y1 + [1,2] * y2 - [1,1.5] *y3
  //                   f2 = [-1, 1] + [-1, -0.5] * y1 + [-1, 0] * y4
  //
  //     compute the sum of the two forms as a third form f3
  //     compute the product of the two forms as a form f4

  AAEnvironment aa(0);
  int i;
  for (i=0; i <= 4; ++i) {
    double d = (double) (i+1);
    i_double I1(-d, d), I2(0), I3(0.2 *d*d, 0.25*d*d);
    aa.makeFreshVariable(I1, I2, I3);
  }

  i_double I1(0.1,0.15), I2(0), I3(0.014,0.016);
  int y = aa.makeFreshVariable(I1,I2, I3);
  
  
  // Print the enviroment.
  cout << "Environment -- > "<<endl;
  cout << aa << endl;
  // Make the two affine forms.

  AffineArithmeticClass f1(aa);
  f1.setConstant( i_double(-1.0, 2.0) );
  f1.setCoefficient(1, i_double(2.0,3.0));
  f1.setCoefficient(2, i_double(1.0,2.0));
  f1.setCoefficient(3, i_double(1.0,1.5));

  cout << "\t f1 = " << f1 << endl;

  AffineArithmeticClass f2(aa);
  f2.setConstant( i_double(-1.0,1.0));
  f2.setCoefficient(1, i_double(-1.0,-0.5));
  f2.setCoefficient(4, i_double(-1.0, 0.0));

  cout << "\t f2 = " << f2 << endl;

  AffineArithmeticClass f3(f1);
  f3.add_assign(f2);

  cout << "\t f1 + f2 = " << f3 << endl;

  AffineArithmeticClass f4(f1);
  f4.sub_assign(f2);

  cout << "\t f1 - f2 = " << f4 << endl;
  // Testing sum 

  AffineArithmeticClass f5(f2);
  f5.scale_assign(i_double(1,2));
  cout << "\t [1,2] * f2  = "<< f5 << endl;

  AffineArithmeticClass f6(f2);
  f6.multiply_assign(f1);
  cout << "\t f1 * f2 = " << f6 << endl;

  AffineArithmeticClass f7(aa);
  f6.clear();
  I1 = i_double(0.2, 0.22);
  f6.setConstant(I1);
  f6.setCoefficient(y, I1);
  cout << "\t f6 = " << f6 << endl;
  f7.sine_assign(f6);
  cout << "\t sin(f6) = " << f7 << endl;
  cout << aa << endl;
}

void test3(){
  AAEnvironment aa(3);
  AffineArithmeticClass x = randomVariable(aa, i_double(-0.5,0.5), i_double(0), i_double(0.08333333333333333));
  x = i_double(2.5) + x;
  cout << "x = " << x << endl;
  AffineArithmeticClass y = square(x) -x ;
  cout << "y = " << y << endl;
  cout << "Range of y = " << y.range() << endl;
  cout << "Expectation: "<< y.expectation() << endl;
  y.performAnalysisForCMI(cout);
  cout << "Environment = " << endl;
  cout << aa << endl;
}

void test4(){
  AAEnvironment aa(4);
  // range(y) = [0.1,0.15], expectation = [0.12,0.13], variance = [0.015, 0.016]
  i_double I1(-0.5, 0.5), I2(0.0), I3(0.2,0.2);
  std::cout << I1 << std::endl;
  int y = aa.makeFreshVariable(I1, I2, I3);
  AffineArithmeticClass x(aa);
  x.setConstant(i_double(1.5,1.5));
  x.setCoefficient(y, i_double(1.0,1.0));
  cout << "x = " << x << endl;
  cout << x.range() << endl; 
  AffineArithmeticClass z(aa);
  z.power_assign(x,7);
  cout << "(x)^7 = " << z << "range:" << z.range() << endl;
  cout << aa << std::endl;
  return;
}

int main( ) {
  //test1();
  test2();
  monteCarloSineCosineTest2();
  test4();
  // testRoboticsExample();
  //  test3();
}
