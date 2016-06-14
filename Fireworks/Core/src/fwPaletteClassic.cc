
namespace fireworks {

static const float s_forWhite[][3] ={
{ 0.79, 0.79, 0.12 }, //yellow (made it a bit darker)
{ 0.47, 0.00, 0.64 }, //purple
{ 0.98, 0.70, 0.00 }, //yellowish-orange
{ 0.18, 0.00, 0.59 }, //purplish-blue
{ 0.98, 0.54, 0.00 }, //orange
{ 0.00, 0.11, 1.00 }, //blue
{ 0.99, 0.26, 0.01 }, //dark orange
{ 0.00, 0.80, 0.78 }, //cyan
{ 1.00, 0.06, 0.00 }, //red
{ 0.33, 0.64, 0.14 }, //green
{ 0.60, 0.06, 0.23 }, //burgundy
{ 0.65, 0.92, 0.17 }, //lime{ 0.99, 1.00, 0.39 },
{ 0.00, 0.46, 1.00 }, //azure+9
{ 1.00, 0.00, 0.40 }, //pink-3
{ 0.02, 1.00, 0.40 }, //teal+8
{ 0.40, 0.40, 0.40 }, //gray
{ 0.00, 0.00, 0.00 }, //black

{ 0.85, 0.85, 0.58 },
{ 0.87, 0.72, 0.92 },
{ 0.99, 0.88, 0.59 },
{ 0.79, 0.72, 0.90 },
{ 1.00, 0.82, 0.59 },
{ 0.71, 0.75, 0.99 },
{ 1.00, 0.80, 0.72 },
{ 0.71, 0.98, 0.95 },
{ 0.99, 0.74, 0.70 },
{ 0.77, 0.86, 0.65 },
{ 0.90, 0.74, 0.79 },
{ 0.67, 0.95, 0.52 },
{ 0.57, 0.78, 1.00 }, //azure+9
{ 1.00, 0.57, 0.74 }, //pink-5
{ 0.73, 1.00, 0.83 }, //teal+9
{ 0.80, 0.80, 0.80 }, //gray
{ 0.60, 0.60, 0.60 }  //blackish gray
};

static const float s_forBlack[][3] ={
{ 1.00, 1.00, 0.20 }, //yellow
{ 0.53, 0.00, 0.69 }, //purple
{ 0.98, 0.74, 0.01 }, //yellowish-orange
{ 0.24, 0.00, 0.64 }, //purplish-blue
{ 0.98, 0.60, 0.01 }, //orange
{ 0.01, 0.14, 1.00 }, //blue
{ 0.99, 0.33, 0.03 }, //dark orange
{ 0.01, 0.83, 0.81 }, //cyan
{ 1.00, 0.09, 0.00 }, //red
{ 0.40, 0.69, 0.20 }, //green
{ 0.65, 0.10, 0.29 }, //burgundy
{ 0.65, 0.92, 0.17 }, //lime
{ 0.00, 0.39, 0.79 }, //azure+9
{ 1.00, 0.00, 0.40 }, //pink-3
{ 0.02, 1.00, 0.40 }, //teal+8
{ 0.70, 0.70, 0.70 }, //gray
{ 1.00, 1.00, 1.00 }, //white

/*
{1.,0.,0.}, //red
{0.,0.,1.}, //blue
{0.,1.,1.}, //cyan
{0.,1.,0.}, //green
{1.,0.,1.}, //magenta
{1.,0.5,0.0},  //orange
{1.,1.,0.}, //yellow
{0.5,0.5,0.5}, //gray
*/
{ 0.27, 0.27, 0.04 },
{ 0.19, 0.00, 0.24 },
{ 0.19, 0.15, 0.00 },
{ 0.14, 0.00, 0.38 },
{ 0.19, 0.11, 0.00 },
{ 0.01, 0.05, 0.33 },
{ 0.17, 0.05, 0.02 },
{ 0.00, 0.33, 0.29 },
{ 0.34, 0.03, 0.01 },
{ 0.15, 0.24, 0.06 },
{ 0.24, 0.02, 0.11 },
{ 0.22, 0.30, 0.07 },
{ 0.00, 0.20, 0.26 }, //azure+8
{ 0.35, 0.00, 0.14 }, //pink-2
{ 0.00, 0.35, 0.12 }, //teal+9
{ 0.22, 0.22, 0.22 }, //gray
{ 0.36, 0.36, 0.36 }  //whitish gray
/*
{0.7,0.0,0.0},
{0.0,0.0,0.7},
{0.0,.7,0.7},
{0.0,.7,0.},
{.7,0.,.7},
{.7,0.4,0.0},
{.7,.7,0.0},
{0.3,0.3,0.3}
 */
};


const static unsigned int s_size = sizeof(s_forBlack)/sizeof(s_forBlack[0]);
}
