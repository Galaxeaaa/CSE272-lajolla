#pragma once

// We use double for most of our computation.
// Rendering is usually done in single precision Reals.
// However, lajolla is a educational renderer with does not
// put emphasis on the absolute performance. 
// We choose double so that we do not need to worry about
// numerical accuracy as much when we render.
// Switching to Realing point computation is easy --
// just set Real = Real.
using Real = double;

// Lots of PIs!
const Real c_PI = Real(3.14159265358979323846);
const Real c_INVPI = Real(1.0) / c_PI;
const Real c_TWOPI = Real(2.0) * c_PI;
const Real c_INVTWOPI = Real(1.0) / c_TWOPI;
const Real c_FOURPI = Real(4.0) * c_PI;
const Real c_INVFOURPI = Real(1.0) / c_FOURPI;
const Real c_PIOVERTWO = Real(0.5) * c_PI;
const Real c_PIOVERFOUR = Real(0.25) * c_PI;