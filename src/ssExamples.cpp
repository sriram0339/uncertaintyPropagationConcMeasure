/* 
   ssExamples: Benchmarks from TACAS 2016 paper by Bouissou et al.
        Copyright (C) 2016-2019  Sriram Sankaranarayanan srirams@colorado.edu

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <https://www.gnu.org/licenses/>.

 */
#include <iostream>
#include <fstream>
#include "aaEnvironment.hpp"
#include "affineArithmeticClass.hpp"
#include <boost/random.hpp>
#include <boost/timer/timer.hpp>
#include <ctime>
#include <cstdlib>
#include <cmath>

using namespace boost::timer;

void simple2DArmMotionExample(){
  cpu_timer affineFormConstructionTime;
  const double angles[10]= {10,60,110,160,140,100,60,20,10,0.0};
  int nAngles = 10;
  AAEnvironment env(0);
  i_double pi(3.1415,3.1416);
  AffineArithmeticClass x = randomVariable(env, i_double(-0.1, 0.1), i_double(0.0), i_double(0.00248,0.00252));
  AffineArithmeticClass y = randomVariable(env, i_double(-0.5, 0.5), i_double(0.0), i_double(0.0098,0.012));

  int nReps = 100;

   for (int kk =0 ; kk < nReps; ++kk){
     for (int i = 0; i < nAngles; ++i){
       i_double i6(-0.02,0.02), i7(0.0), i8(0.000133);
       AffineArithmeticClass d = randomVariable(env, i6, i7, i8);
       d = i_double(1.0) + d;
       i_double theta = i_double(angles[i]) * pi/i_double(180.0);
       i_double i9(-0.05,0.05), i10(0), i11(0.0001,0.00011);
       AffineArithmeticClass tmp0 = randomVariable(env, i9, i10, i11);
       AffineArithmeticClass tmp = (i_double(1.0) + tmp0)*theta; 
       x = x + d  * cosine(tmp);
       y = y + d * sine(tmp);
       
     }
   }
   affineFormConstructionTime.stop();
   cpu_timer analysisTime;
   cout << "x = " << x << endl;
   cout << "Range: " << x.range() << endl;
   cout << "Expectation: " << x.expectation() << endl;
   x.performAnalysisForCMI(cout);
   analysisTime.stop();
   
    // cout << "y = " << y << endl;
    // cout << "Range: " << y.range() << endl;
    // cout << "Expectation: " << y.expectation() << endl;
   // y.performAnalysisForCMI(cout);
    cout << "Time taken to construct Affine form = " << affineFormConstructionTime.format()<< endl;
  cout << "Time taken to analyze Affine form = " << analysisTime.format()<< endl;
  
}

void carSteeringExample(){
  cpu_timer affineFormConstructionTime;
  AAEnvironment env(1);
  const double angles[25]={0, -3, 0, -3,0,-2,0,-1,-3,0,-3,0, 0, -3, 0, -3,0,-2.5,0,0,-2.5,0,-2.6,-2.7,0};
  const double times[25] = {15,15,10,15,20,30,10,4,30,20,10,30, 15,15,10,15,20,30,10,10,10,20,15,15,45};
  const double speeds[25] = {55,12, 60,10,60,12,50,10,25,60,10,45, 55,12, 40,10,60,12,20,10,25,60,10,45,50};
  int nMoves = 15;
  i_double delta(0.1);
  int nRounds = 3;
  i_double pi (3.1415);
  i_double i1 (-0.05,0.05);
  i_double i2(0.0);
  i_double i3(1e-05);

  AffineArithmeticClass x (env), y(env), theta(env);
  i_double L (25.0);
  for (int j = 0; j < nRounds ; ++j){
    cout << "Round # " << j << endl;
    for (int i = 0; i < nMoves; ++i){
      i_double a = i_double(angles[i])* (pi/i_double(180.0));
      int  T = times[i];
      i_double v = speeds[i];
      AffineArithmeticClass tanphi= randomVariable(env, i1, i2, i3);
      tanphi = i_double(1.0) + tanphi;
      tanphi.scale_assign(tan(a));
      AffineArithmeticClass vel= randomVariable(env, i1, i2, i3);
      vel  = i_double(1.0) + vel;
      vel.scale_assign(v);
      for (int k = 0; k < T; ++k){
	x = x + delta* vel * cosine(theta);
	y = y + delta * vel * sine(theta);
	theta = theta + (delta/L) * vel * tanphi;
	
      }

    }
  }

  affineFormConstructionTime.stop();
  cpu_timer analysisTime;
   cout << "x = " << x << endl;
   cout << "Range: " << x.range() << endl;
   cout << "Expectation: " << x.expectation() << endl;
   x.performAnalysisForCMI(cout);
   analysisTime.stop();
   cout << "Time taken to construct Affine form = " << affineFormConstructionTime.format()<< endl;
  cout << "Time taken to analyze Affine form = " << analysisTime.format()<< endl;
  
  
}


void modifiedTankFillingExample(){
  cpu_timer affineFormConstructionTime;
  AAEnvironment env(2);
  AffineArithmeticClass volMeasured(env);
  int N = 26;
  
  for (int i = 0; i < N; ++i){
    AffineArithmeticClass f = randomVariable(env, i_double(-0.03,0.03), i_double(0.0), i_double(0.0003));
    f = i_double(0.1) + f;
    AffineArithmeticClass e = randomVariable(env, i_double(-0.03, 0.03), i_double(0.0), i_double(0.0001));
    volMeasured = volMeasured + f + e;
  }
  affineFormConstructionTime.stop();
  cpu_timer analysisTime ;
  cout << " After " << N << " Iterations : " << volMeasured << endl;
  
  cout << "Range: " << volMeasured.range() << endl;
  cout << "Expectation: " << volMeasured.expectation() << endl;
  volMeasured.performAnalysisForCMI(cout);
  analysisTime.stop();
  cout <<"Time taken to construct affine form:"<<affineFormConstructionTime.format() << endl;
  cout << "Time taken to analyze Affine form = " << analysisTime.format()<< endl;
}

void anesthesiaInfusionExample(){
  cpu_timer affineFormConstructionTime;
  AAEnvironment env(3);
  double infusionTimings[7] = {20, 15, 15, 15, 15, 15, 45};
  double infusionRates[7] = { 3, 3.2, 3.3, 3.4, 3.2, 3.1, 3.0};
  AffineArithmeticClass x1(env), x2(env), x3(env), x4(env);
  i_double e0(-0.4, 0.4), e1(0.0), e2(0.006,0.0064);
  
  for (int i = 0; i < 7; ++i){
    double currentInfusion= 20.0*infusionRates[i];
    int curTime = infusionTimings[i];
    for (int j = 0; j < 40*curTime; ++j){
      AffineArithmeticClass e  = randomVariable(env, e0, e1, e2);
      e = e + i_double(1.0);
      AffineArithmeticClass u = e * i_double(currentInfusion);
      AffineArithmeticClass x1n = i_double(0.9012)* x1 + i_double(0.0304) * x2 + i_double(0.0031) * x3 + i_double(2.676e-1) * u;
      AffineArithmeticClass x2n = i_double(0.0139)* x1 + i_double(0.9857) * x2 + i_double(2e-3)*u;
      AffineArithmeticClass x3n = i_double(0.0015) * x1 + i_double(0.9985) * x3+ i_double(2e-4)*u ;
      AffineArithmeticClass x4n = i_double(0.0838) * x1 + i_double(0.0014) * x2 + i_double(0.0001) *x3 + i_double(0.9117) * x4 + i_double(12e-3) * u;
      x1 = x1n;
      x2 = x2n;
      x3 = x3n;
      x4 = x4n;
    }

  }
  
  affineFormConstructionTime.stop();
  cpu_timer analysisTimer;
  cout << "x4 = " << x4 << endl;
  cout << "Range: " << x4.range() << endl;
  cout << "Expectation: " << x4.expectation() << endl;
  x4.performAnalysisForCMI(cout);
  cout << "Time taken to construct Affine form = " << affineFormConstructionTime.format()<< endl;
  cout << "Time taken to analyze Affine form = " << analysisTimer.format()<< endl;
  
}

void tumorModelExample(){
  cpu_timer affineFormConstructionTime;
  i_double a(1.0), b0(1.0), beta(1.0);
  AAEnvironment env(4);
  AffineArithmeticClass X(env);
  X.setConstant(i_double(1.0));
  int nExecs = 100;
  for (int i = 0; i < nExecs; ++i){
    // AffineArithmeticClass omega = randomVariable(env,i_double(-0.1,0.1), i_double(0.0), i_double(0.0001));
    i_double I = X.range();
    set<int> vars;
    X.collectRelevantNoiseSymbols(vars);
    i_double rngXOmega = i_double(-0.1,0.1) * I;
    i_double expectXOmega = i_double(0.0); // Because of independence.
    i_double secondMoment = i_double(0.0001) * X.secondMoment_raw(); // Again, this is to avoid causing a huge blowup
    // omega = omega * X; // SS: This will blow up -- we need an alternative
    AffineArithmeticClass omega = depRandomVariable(env,vars,rngXOmega,expectXOmega,secondMoment);
    AffineArithmeticClass tmp1(env);
    tmp1.specialFunctionForTumorModel_assign(X);
    AffineArithmeticClass tmp2 = X - omega - tmp1;
    X = X + i_double(0.1) * tmp2;
  }
  affineFormConstructionTime.stop();

  cpu_timer analysisTimer;
   cout << "X = " << X << endl;
   cout << "Range: " << X.range() << endl;
   cout << "Expectation: " << X.expectation() << endl;
   X.performAnalysisForCMI(cout);
   cout << "Time taken to construct Affine form = " << affineFormConstructionTime.format()<< endl;
   cout << "Time taken to analyze Affine form = " << analysisTimer.format()<< endl;
}


void doubleWellExample(){
  cpu_timer affineFormConstructionTime;
  AAEnvironment env(5);
  i_double varepsilon(0.05);
  AffineArithmeticClass X(env);
  AffineArithmeticClass Y(env);
  X.setConstant(i_double(0.5));
  Y.setConstant(i_double(0.5));
  int nExecs = 50;
  for (int i = 0; i < nExecs; ++i){
    
    AffineArithmeticClass omegaX = randomVariable(env,i_double(-0.5,0.5), i_double(0.0), i_double(0.0001));
    AffineArithmeticClass omegaY = randomVariable(env,i_double(-0.5,0.5), i_double(0.0), i_double(0.0001));
    AffineArithmeticClass X_3 = power(X,3);
    AffineArithmeticClass Y_3 = power(Y,3);
    AffineArithmeticClass X1n = i_double(-0.1)*X + i_double( 0.1)*X_3 + varepsilon*omegaX;
    AffineArithmeticClass Y1n = i_double( 0.1)*Y + i_double(-0.1)*Y_3 + varepsilon*omegaY;
    X = X + X1n;
    Y = Y + Y1n;
  }
  affineFormConstructionTime.stop();

  cpu_timer analysisTimer;
   cout << "X = " << X << endl;
   cout << "Range: " << X.range() << endl;
   cout << "Expectation: " << X.expectation() << endl;
   X.performAnalysisForCMI(cout);

   cout << "Y = " << Y << endl;
   cout << "Range: " << Y.range() << endl;
   cout << "Expectation: " << Y.expectation() << endl;
   Y.performAnalysisForCMI(cout);
   
   cout << "Time taken to construct Affine form = " << affineFormConstructionTime.format()<< endl;
   cout << "Time taken to analyze Affine form = " << analysisTimer.format()<< endl;
}

void rimlessWheelExample(){
  cpu_timer affineFormConstructionTime;
  AAEnvironment env(6);
  i_double n_range (-8.0,8.0);
  i_double n_mean (0.0);
  i_double n_var (2.25);
  i_double offs(8.0);
  i_double pi(3.14159265359);
  AffineArithmeticClass x(env);
  int N=1000;
  i_double theta(30.0);
  i_double deg2rad = pi/i_double(180.0);
  i_double cos2theta = square(cos(theta*deg2rad));
  i_double g (10.0);
  i_double L(1.0);
  i_double twoGByL = (i_double(2.0)*g/L);

  x.setConstant(1.0);

  
  for(int i=0; i < N; ++i){
    AffineArithmeticClass gamma = randomVariable(env, n_range, n_mean, n_var);
    gamma = gamma + offs;
    AffineArithmeticClass beta1 = theta/i_double(2.0) + gamma;
    AffineArithmeticClass beta2 = theta/i_double(2.0) - gamma;
    beta1 = deg2rad*beta1;
    beta2 = deg2rad*beta2;
    x = (cos2theta * ( x  + twoGByL  *(i_double(1.0) - cosine(beta1)) )) - twoGByL * (i_double(1.0)- cosine(beta2));
  }

  affineFormConstructionTime.stop();
  cpu_timer analysisTimer;
  cout << "x = " << x << endl;
  cout << "Range: " << x.range() << endl;
  cout << "Expectation: " << x.expectation() << endl;
  x.performAnalysisForCMI(cout);
  
   cout << "Time taken to construct Affine form = " << affineFormConstructionTime.format()<< endl;
   cout << "Time taken to analyze Affine form = " << analysisTimer.format()<< endl;
  
}

void cartPoleBalanceModel(){
  bool cheapMul = true;
  cpu_timer affineFormConstructionTime;
  i_double delta(0.1);
  i_double nRange(-10.0,10.0);
  i_double nMean(0.0);
  i_double nVar(1.0);
  AAEnvironment env(8);
  AffineArithmeticClass x(env), dxdt(env), theta(env), dthetadt(env);
  int N = 10;
  i_double L (0.5);
  i_double g(9.8);
  i_double mp(1.0);
  i_double mc(10.0);

  for (int i = 0; i< N; ++i){
    cout << i << endl;
    AffineArithmeticClass w1 = randomVariable(env, nRange, nMean, nVar);
    AffineArithmeticClass w2 = randomVariable(env, nRange, nMean, nVar);
    AffineArithmeticClass w3 = randomVariable(env, nRange, nMean, nVar);
    AffineArithmeticClass w4 = randomVariable(env, nRange, nMean, nVar);
    AffineArithmeticClass feedbackTerm (env);

    feedbackTerm = i_double(10.0)*x + i_double(-289.83)*theta + i_double(19.53)*dxdt + i_double(-63.25)*dthetadt;

    
    w1 = w1 * i_double(0.01)* sqrt(delta);
    w2 = w2 * i_double(0.01) * sqrt(delta);
    w3 = w3 * i_double(0.03)* sqrt(delta);
    w4 = w4 * i_double(0.03) * sqrt(delta);
    
    AffineArithmeticClass tmp0(env);
    tmp0.power_assign(dthetadt,2);
    tmp0.scale_assign(L);

    // tmp1 = (mc + mp sin(theta)^2)
    AffineArithmeticClass tmp1(env);
    tmp1 = power(sine(theta),2);
    tmp1.scale_assign(mp);
    i_double c = tmp1.getConstant();
    c = c + mc;
    tmp1.setConstant(c);
    // denTerm = 1/(mc + mp sin(theta)^2)
    AffineArithmeticClass denTerm(env);
    denTerm.inverse_assign(tmp1);

    
    AffineArithmeticClass ddx(env);
    if (cheapMul){
      AffineArithmeticClass tmp4(env);
      tmp4 = mp * sine(theta);
      AffineArithmeticClass tmp5(env);
      tmp5 = tmp0 - g * cosine(theta);
      tmp4.linearized_multiply_assign(tmp5);
      ddx = feedbackTerm - tmp4 ;
      ddx.linearized_multiply_assign(denTerm);
    } else {
      ddx = (feedbackTerm - mp*sine(theta)*(tmp0 - g * cosine(theta)) ) * denTerm;
    }
    AffineArithmeticClass ddtheta(env);
    if (cheapMul){
      AffineArithmeticClass tmp3(feedbackTerm);
      tmp3.linearized_multiply_assign(cosine(theta));
      tmp0.linearized_multiply_assign(sine( i_double(2.0) * theta));
      tmp0.scale_assign (i_double(0.5)* mp);
      ddtheta = tmp3 - tmp0 + (mc+mp)*g*sine(theta);
      ddtheta.linearized_multiply_assign(denTerm);
      ddtheta.scale_assign(i_double(1.0)/L);
    } else {
      ddtheta = (feedbackTerm * cosine(theta) - mp * tmp0 * i_double(0.5)*sine(i_double(2.0)*theta)  + (mc+mp)*g * sine(theta)) * denTerm * (i_double(1.0)/L);
    }
    x = x + delta * dxdt + delta * w1;
    theta = theta + delta * dthetadt + delta * w2;
    dxdt = x + delta * ddx + delta * w3;
    dthetadt = dthetadt+ delta * ddtheta + delta * w4;
    
  }

  affineFormConstructionTime.stop();
  cpu_timer analysisTimer;
  cout << "theta = " << theta << endl;
  cout << "Range: " << theta.range() << endl;
  cout << "Expectation: " << theta.expectation() << endl;
  theta.performAnalysisForCMI(cout);
  
  cout << "Time taken to construct Affine form = " << affineFormConstructionTime.format()<< endl;
  cout << "Time taken to analyze Affine form = " << analysisTimer.format()<< endl;
   
  
  

}

void cartPoleBalanceModelPolynomial(){
  cpu_timer affineFormConstructionTime;
  i_double delta(0.1);
  i_double nRange(-10.0,10.0);
  i_double nMean(0.0);
  i_double nVar(1.0);
  AAEnvironment env(8);
  AffineArithmeticClass x(env), dxdt(env), theta(env), dthetadt(env);
  int N = 10;
 

  for (int i = 0; i< N; ++i){
    cout << i << endl;
    AffineArithmeticClass w1 = randomVariable(env, nRange, nMean, nVar);
    AffineArithmeticClass w2 = randomVariable(env, nRange, nMean, nVar);
    AffineArithmeticClass w3 = randomVariable(env, nRange, nMean, nVar);
    AffineArithmeticClass w4 = randomVariable(env, nRange, nMean, nVar);
    AffineArithmeticClass feedbackTerm (env);

    feedbackTerm = i_double(10.0)*x + i_double(-289.83)*theta + i_double(19.53)*dxdt + i_double(-63.25)*dthetadt;

    
    w1 = w1 * i_double(0.01)* sqrt(delta);
    w2 = w2 * i_double(0.01) * sqrt(delta);
    w3 = w3 * i_double(0.03)* sqrt(delta);
    w4 = w4 * i_double(0.03) * sqrt(delta);

    AffineArithmeticClass ddx(env);

    ddx = i_double(0.1)*feedbackTerm + i_double(0.98)*theta;
    AffineArithmeticClass tmp0(env);
    tmp0= i_double(-0.75)*theta + i_double(-0.1)*feedbackTerm;

    AffineArithmeticClass thetaSq(env);
    thetaSq.power_assign(theta,2);
    
    tmp0.linearized_multiply_assign( thetaSq);

    AffineArithmeticClass tmp1(env);
    tmp1.power_assign(dthetadt,2);
    tmp1.scale_assign(i_double(-0.05));
    tmp1.linearized_multiply_assign(theta);

    ddx = ddx + tmp0 + tmp1;

    AffineArithmeticClass ddtheta(env);
    ddtheta = i_double(0.2)* feedbackTerm + i_double(21.56) * theta;

    AffineArithmeticClass tmp2(env);
    tmp2 = i_double(-5.75)*theta + i_double(-0.12) * feedbackTerm;
    tmp2.linearized_multiply_assign(thetaSq);

    AffineArithmeticClass tmp3(env);
    tmp3.power_assign(dthetadt,2);
    tmp3.scale_assign(i_double(-0.1));

    ddtheta = ddtheta + tmp2 + tmp3;
    
    x = x + delta * dxdt + w1;
    theta = theta + delta * dthetadt + w2;
    dxdt = x + delta * ddx + w3;
    dthetadt = dthetadt+ delta * ddtheta + w4;
    

    
  }

    affineFormConstructionTime.stop();
  cpu_timer analysisTimer;
  cout << "theta = " << theta << endl;
  cout << "Range: " << theta.range() << endl;
  cout << "Expectation: " << theta.expectation() << endl;
  theta.performAnalysisForCMI(cout);
  
  cout << "Time taken to construct Affine form = " << affineFormConstructionTime.format()<< endl;
  cout << "Time taken to analyze Affine form = " << analysisTimer.format()<< endl;
  
}


int main(int argc, char * argv[]){
  int sel=-1;
  if (argc >= 2){
    sel =atoi(argv[1]);
  }
  switch(sel){
  case 1:
    cout << "Simple 2D arm motion example  " << endl;
    simple2DArmMotionExample();
    break;
  case 2: 
    cout << " Cart steering example " << endl;
    carSteeringExample();
    break;
  case 3:
    cout << "Tank filling " << endl;
    modifiedTankFillingExample();
    break;
  case 4:
    cout << "Anesthesia model simulation" << endl;
    anesthesiaInfusionExample();
    break;
  case 5:
    cout << "Tumor Model"  << endl;
    tumorModelExample();
    break;
  case 6:
    cout << "Double Well Model" << endl;
    doubleWellExample();
    break;
  case 7:
    cout << "Rimless wheel model" << endl;
    rimlessWheelExample();
    break;
  case 8:
    cout << "Cartpole balance model" << endl;
    cartPoleBalanceModel();
    break;

  case 9:
     cout << "Polynomial Cartpole balance model" << endl;
    cartPoleBalanceModelPolynomial();
    break;
  default:
    cout << "Unknown model" <<endl;
    
  }
}
