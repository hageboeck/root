// Tests for the RooRealBinding, and numeric integrators that use it.
// Author: Stephan Hageboeck, CERN  05/2020

#include "RooRealBinding.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooIntegrator1D.h"
#include "RooNumIntConfig.h"
#include "RooHelpers.h"
#include <Math/ProbFuncMathCore.h>
#include <Math/SpecFuncMathCore.h>

#include "OldRooIntegrator1D.cxx"

#include "gtest/gtest.h"

///
TEST(RooRealBinding, BatchEvalFeature) {
  RooRealVar a("a", "a", -100, 100);
  RooRealVar b("b", "b", -100, 100);
  RooFormulaVar formula("formula", "1.3*a + 1.4*b", RooArgList(a, b));

  std::vector<double> as;
  std::vector<double> bs;
  std::generate_n(std::back_inserter(as), 10, [](){ static double val = 0; return val += 0.3;});
  std::generate_n(std::back_inserter(bs), 10, [](){ static double val = 0; return val += 0.4;});

  std::vector<RooSpan<const double>> data;
  data.emplace_back(as);
  data.emplace_back(bs);

  RooRealBinding binding(formula, RooArgSet(a, b));
  auto result = binding.getValBatch(data);
  for (unsigned int i=0; i < 10; ++i) {
    EXPECT_DOUBLE_EQ(result[i], 1.3 * as[i] + 1.4 * bs[i]) << "result[" << i << "] a=" << as[i] << " b=" << bs[i];
  }
}


TEST(Roo1DIntegrator, RunFormulaVar_Trapezoid) {

  RooRealVar x("x", "x", -100, 100);
  RooRealVar a("a", "a", 0.2, -100, 100);
  RooRealVar b("b", "b", 0.3, -100, 100);
  RooFormulaVar formula("formula", "0.1 + x*(a + b*x)", RooArgList(x, a, b));
  auto solution = [&x](double a, double b){
    auto indefInt = [=](double y){
      return y*(0.1 + y*(1./2.*a + 1./3.*b * y));
    };
    return indefInt(x.getMax()) - indefInt(x.getMin());
  };
  RooRealBinding binding(formula, RooArgSet(x, a, b));

  // The integrators will warns, since we let them run until maxSteps
  RooHelpers::HijackMessageStream hijack(RooFit::WARNING, RooFit::Integration);

  // Test the recursion anchors of the Romberg integration
  {
    RooIntegrator1D oneStep(binding, RooIntegrator1D::Trapezoid, 1, 1.E-15);
    EXPECT_DOUBLE_EQ(oneStep.integral(), 0.5*200.*(2*0.1 + 2.*0.3*10000.));
    x = -100.;
    const double left = formula.getVal();
    x = 100.;
    const double right = formula.getVal();
    // Run integral again, also to make sure that setting x has no effect:
    EXPECT_DOUBLE_EQ(oneStep.integral(), 0.5*200.*(left + right));

    RooIntegrator1D twoStep(binding, RooIntegrator1D::Trapezoid, 2, 1.E-15);
    x = 0.;
    const double middle = formula.getVal();
    twoStep.applySeriesAcceleration(false);
    const double noAccel = twoStep.integral();
    EXPECT_DOUBLE_EQ(noAccel, 0.25*200.*(left + right) + 0.5*200.*middle);

    twoStep.applySeriesAcceleration(true);
    const double accel = twoStep.integral();
    EXPECT_LT(fabs(accel - solution(a.getVal(), b.getVal())), 0.8 * fabs(noAccel - solution(a.getVal(), b.getVal())))
        << "Expect with acceleration to be better than without.";
  }

  // Now run many steps
  {
    constexpr unsigned int nSteps = 25;
    constexpr double relEps = 1.E-50;
    RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, relEps);
    double inputs[] = {1., 3.123};
    double target = solution(1., 3.123);
    EXPECT_LT(fabs(integrator.integral(inputs) - target)/target, 1.E-14);

    target = solution(a.getVal(), b.getVal());
    EXPECT_LT(fabs(integrator.integral() - target)/target, 1.E-14);

    inputs[0] = 4.; inputs[1] = 5.;
    target = solution(4., 5.);
    EXPECT_LT(fabs(integrator.integral(inputs) - target)/target, 1.E-14);
  }
}


TEST(Roo1DIntegrator, RunQuarticFormulaVar) {
  constexpr unsigned int nSteps = 25;
  constexpr double relEps = 1.E-50;
  RooRealVar x("x", "x", -50, 50);
  RooRealVar a("a", "a", 0.2, -100, 100);
  RooRealVar b("b", "b", 0.3, -100, 100);
  RooRealVar c("c", "c", 0.4, -100, 100);
  RooRealVar d("d", "d", 0.5, -100, 100);
  RooFormulaVar formula("formula", "0.1 + x*(a + x*(b + x*(c + d * x)))", RooArgList(x, a, b, c, d));
  auto solution = [&x](double a, double b, double c, double d){
    auto indefInt = [=](double y){
      return y*(0.1 + y*(1./2.*a + y*(1./3.*b + y*(1./4.*c + 1./5.*d * y))));
    };
    return indefInt(x.getMax()) - indefInt(x.getMin());
  };
  RooRealBinding binding(formula, RooArgSet(x, a, b, c, d));
  RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, relEps);

  double target = solution(0.2, 0.3, 0.4, 0.5);
  EXPECT_LT(fabs(integrator.integral() - target)/target, 1.E-13);
}


TEST(Roo1DIntegrator, ConvergenceSettings) {
  constexpr unsigned int nSteps = 25;
  RooRealVar x("x", "x", 0.1, 50);
  RooRealVar a("a", "a", 0.2, -100, 100);
  RooFormulaVar formula("formula", "log(a*x)", RooArgList(x, a));
  auto solution = [&x](double a){
    auto indefInt = [=](double y){
      return 1./a * (y * log(y) - y);
    };
    return indefInt(a*x.getMax()) - indefInt(a*x.getMin());
  };
  RooRealBinding binding(formula, RooArgSet(x, a));

  for (double relEps : {0.1, 1.E-3, 1.E-6, 1.E-8}) {
    RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, relEps);

    double target = solution(0.2);
    double integral = integrator.integral();
    EXPECT_LT(fabs(integral - target)/target, relEps) << "With integral=" << integral << "\ttarget=" << target;
    EXPECT_GT(fabs(integral - target)/target, relEps/1000.) << "With integral=" << integral << "\ttarget=" << target;
  }
}


TEST(Roo1DIntegrator, RunVsOldIntegrator) {
  constexpr unsigned int nSteps = 25;
  constexpr double relEps = 1.E-50;
  RooRealVar x("x", "x", -100, 100);
  RooRealVar a("a", "a", 0.2, -100, 100);
  RooRealVar b("b", "b", 0.3, -100, 100);

  RooFormulaVar formula("formula", "0.1 + x*(a + b*x)", RooArgList(x, a, b));
  auto solution = [&x](double a, double b){
    auto indefInt = [=](double y){
      return y*(0.1 + y*(1./2.*a + 1./3.*b * y));
    };
    return indefInt(x.getMax()) - indefInt(x.getMin());
  };
  RooRealBinding binding(formula, RooArgSet(x, a, b));

  RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, relEps);
  OldRooIntegrator1D old1D(binding, OldRooIntegrator1D::Trapezoid, nSteps, relEps);

  double inputs[2];
  inputs[0] = 0.2; inputs[1] = 0.3;
  a = 0.2;
  b = 0.3;
  double target = solution(0.2, 0.3);
  EXPECT_LE(fabs(integrator.integral(inputs) - target), fabs(old1D.integral(inputs) - target));

  target = solution(4.4, 5.5);
  a = 4.4;
  b = 5.5;
  EXPECT_LE(fabs(integrator.integral() - target), fabs(old1D.integral() - target));
}


TEST(Roo1DIntegrator, RunErf) {
  const double min=0, max=1;
  RooRealVar theX("theX", "x", min, max);
  RooFormulaVar gaus("gaus", "ROOT::Math::gaussian_pdf(theX, 1, 0)", theX);
  RooRealBinding binding(gaus, theX);
  double targetError = 0.;

  RooHelpers::HijackMessageStream hijack(RooFit::WARNING, RooFit::Integration);

  for (unsigned int nSteps = 4; nSteps < 24; ++nSteps) {
    RooIntegrator1D integrator(binding, RooIntegrator1D::Trapezoid, nSteps, 1.E-20);
    const double integral = integrator.integral();
    const double error = fabs(integral - (ROOT::Math::gaussian_cdf(max, 1, 0) - ROOT::Math::gaussian_cdf(min, 1, 0)));
    if (nSteps == 4) {
      targetError = error;
    } else {
      // Error should go down faster than 2^nSteps because of series acceleration.
      targetError /= 3.;
      // But cannot be better than double precision
      EXPECT_LT(error, std::max(targetError, 1.E-16) )    << "For step " << nSteps << " with integral=" << integral;
    }
    if (nSteps > 10)
      EXPECT_LT(error / integral, 1.E-4) << "For step " << nSteps << " with integral=" << integral;
    if (nSteps > 15)
      EXPECT_LT(error / integral, 1.E-6) << "For step " << nSteps << " with integral=" << integral;
    if (nSteps > 21)
      EXPECT_LT(error / integral, 1.E-8) << "For step " << nSteps << " with integral=" << integral;
  }
}

TEST(Roo1DIntegrator, RunErf_Midpoint) {
  const double min=0, max=1;
  RooRealVar theX("theX", "x", min, max);
  RooFormulaVar gaus("gaus", "ROOT::Math::gaussian_pdf(theX, 1, 0)", theX);
  RooRealBinding binding(gaus, theX);
  double targetError = 0.;

  RooHelpers::HijackMessageStream hijack(RooFit::WARNING, RooFit::Integration);

  for (unsigned int nSteps = 4; nSteps < 24; ++nSteps) {
    RooIntegrator1D integrator(binding, RooIntegrator1D::Midpoint, nSteps, 1.E-20);
    const double integral = integrator.integral();
    const double error = fabs(integral - (ROOT::Math::gaussian_cdf(max, 1, 0) - ROOT::Math::gaussian_cdf(min, 1, 0)));
    if (nSteps == 4) {
      targetError = error;
    } else {
      // Error should go down faster than 2^nSteps because of series acceleration.
      targetError /= 3.;
      // But cannot be better than double precision
      EXPECT_LT(error, std::max(targetError, 1.E-16) )    << "For step " << nSteps << " with integral=" << integral;
    }
    if (nSteps > 10)
      EXPECT_LT(error / integral, 1.E-4) << "For step " << nSteps << " with integral=" << integral;
    if (nSteps > 15)
      EXPECT_LT(error / integral, 1.E-6) << "For step " << nSteps << " with integral=" << integral;
    if (nSteps > 21)
      EXPECT_LT(error / integral, 1.E-8) << "For step " << nSteps << " with integral=" << integral;
  }
}
