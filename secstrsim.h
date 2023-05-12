// 0=Coil 1=Helix 2=strand 3=gap
char SS_code[]=" HE-";
float SS_sim_score[4][4]={
  // Ploop_C1_Mammoth.prop norm=933 PC_div=0.147
  //Coil    Helix    Strand  Gap
  { 0.930, -1.518, -1.502,  0.715}, // Coil
  {-1.518,  0.784, -7.512, -0.506}, // Helix
  {-1.502, -7.512,  1.342, -0.886}, // Strand
  { 0.715, -0.506, -0.886, -2.303}  // Gap
};
  /*  { 0.8524, -1.59,  -0.9157, 0.9273}, // Coil
  {-1.59,    0.593, -5.699, -1.117}, // Helix
  {-0.9157, -5.699,  1.868, -1.951}, // Strand
  { 0.9273, -1.117, -1.951, -2.303}  // Gap
  };*/
  /*
  {0.8257, -0.9366, -0.8824, 0.6194}, // Coil
  {-0.9366, 0.4889, -3.232, -0.08375},// Helix
  {-0.8824,-3.232,   1.696, -1.762},  // Strand
  {0.6194, -0.08375,-1.762, -2.303}   // Gap
  };*/
//# Bin 3: <PC_div>= 0.626 norm= 6.01e+03
//# Secondary structure: Frequency
//# Coil	Helix	Strand	Gap
//# 0.272	0.509	0.157	0.062
//# Secondary structure: Propensities
//#SS	Coil	Helix	Strand	Gap
//Coil	0.8257	-0.9366	-0.8824	0.6194
//Helix	-0.9366	0.4889	-3.232	-0.08375
//Strand	-0.8824	-3.232	1.696	-1.762
//Gap	0.6194	-0.08375	-1.762	-2.303
