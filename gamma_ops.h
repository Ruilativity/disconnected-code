#ifndef _GAMMA_OPS_H_
#define _GAMMA_OPS_H_

// Calcualte \Gamma, which constructs operators

SpinMatrix gamma_ops(int i) {
  // QDP uses Degrand-Rossi basis of gamma matrices,
  // which are transposes of the normal Euclidean gamma matrices.
  // This program calculate the Dirac matrices that constructs
  // the operators in EDM analysis in LANL

  // --- Gamma matrix ordering in Chroma ---
  // Here (-) sign denotes the correction go Gamma(i)
  // For example, Gamma(13) = g_0 g_2 g_3 = - g_1 g_0 g_1 g_2 g_3 = - g_1 g_5 = -axial
  //
  // Scalar       : i = 0 (gamma_0)
  // Pseudoscalar : i = 15 (gamma_5)
  // vector       : i = 1,2,4,8 (gamma_mu)
  // axial-vector : i = -7,11,-13,14 (gamma_mu gamma_5)
  // tensor       : i = 3,5,6,9,10,12 (gamma_mu gamma_nu)

  // Rearrange the gamma matrix order as follows:
  // Scalar       : i = 0 (gamma_0)
  // Pseudoscalar : i = 1 (gamma_5)
  // vector       : i = 2,3,4,5 (gamma_mu)
  // axial-vector : i = 6,7,8,9 (gamma_mu gamma_5)
  // tensor       : i = 10,11,12,13,14,15 (gamma_mu gamma_nu)

  // In the following switch-case statement, we compensate
  // overall sign due to axial-vector from Chroma Gamma(i)

  int g;
  SpinMatrix s = 1.0;
  switch(i) {
    case  0:  g =  0;         break;  //scalar; I
    case  1:  g = 15;         break;  //pseudoscalar; g_5
    case  2:  g =  1;         break;  //vector; g_0
    case  3:  g =  2;         break;  //vector; g_1
    case  4:  g =  4;         break;  //vector; g_2
    case  5:  g =  8;         break;  //vector; g_3
    case  6:  g = 14;         break;  //axial;  g_0 g_5
    case  7:  g = 13; s=-1.0; break;  //axial;  g_1 g_5
    case  8:  g = 11;         break;  //axial;  g_2 g_5
    case  9:  g =  7; s=-1.0; break;  //axial;  g_3 g_5
    case 10:  g =  9;         break;  //tensor; g_0 g_3
    case 11:  g = 10;         break;  //tensor; g_1 g_3
    case 12:  g = 12;         break;  //tensor; g_2 g_3
    case 13:  g =  3;         break;  //tensor; g_0 g_1
    case 14:  g =  6;         break;  //tensor; g_1 g_2
    case 15:  g =  5;         break;  //tensor; g_0 g_2
    default:
      std::cerr << "Erorr! (gamma_ops) out of range of i=" << i << std::endl;
      exit(1);
  } // switch(i)

  SpinMatrix ret = s*Gamma(g);

  return ret;
}

#endif
