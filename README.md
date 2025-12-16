# ves (Virtual Element Spaces)

This is a library containing helper functions and structures for the implementation of virtual element methods.

The principle of the library is to provide the main projectors for different VEM formulations.

## Available formulations

### Two dimensions

For two-dimensional elements there are three formulations of the method available:

* [Modified VEM](https://www.sciencedirect.com/science/article/pii/S0898122113003179#s000015) (`V2D`)

* [Serendipity VEM](https://www.sciencedirect.com/science/article/pii/S0045793016300391) (`SV2D`)

* [Serendipity Enlarged Enhanced VEM](https://arxiv.org/abs/2210.02653) (`SE2V2D`)


The main workflow to work with these formulations is to provide them with a polygon (`std::vector<Eigen::Vector2d>`) and the desired element order, and have the projectors associated with this formulation computed and available.

The idea is for this to be used alongside the [Projector Assembly](https://www.scipedia.com/public/Moherdaui_et_al_2024a) procedure to instrumentalize into full blown implementations of VEM.

Example usage:
```cpp
#include "ves.h"

constexpr int nv = 5;
constexpr int order = 2;
std::vector<Eigen::Vector2d> poly = regularPolygon(nv);

ves::V2D VE(poly, k);
// H^1 semi-norm projector \Pi^\nabla_k
const Eigen::MatrixXd PiGrad = VE.PiGrad();
// L^2 projector \Pi^0_k
const Eigen::MatrixXd Pi0 = VE.Pi0();

ves::SV2D SVE(poly, k);
// Order of polynomial space associated with internal degrees of freedom
const int innerOrder = SVE.InnerOrder();
// L^2 projector \Pi^0_k
const Eigen::MatrixXd Pi0 = SVE.Pi0();
// L^2 projector of x-derivative \Pi^0_{k-1}\partial_x
const Eigen::MatrixXd Pi0Dx = SVE.Pi0Dx();
// L^2 projector of y-derivative \Pi^0_{k-1}\partial_y
const Eigen::MatrixXd Pi0Dy = SVE.Pi0Dy();

ves::SE2V2D SE2VE(poly, k);
// Order of polynomial space employed in derivative projectors.
const int gradOrder = SE2VE.SFGradOrder();
// Order of polynomial space associated with internal degrees of freedom
const int innerOrder = SE2VE.InnerOrder();
// L^2 projector \Pi^0_k
const Eigen::MatrixXd Pi0 = SVE.Pi0();
// L^2 projector of x-derivative \Pi^0_{gradOrder}\partial_x
const Eigen::MatrixXd Pi0Dx = SVE.Pi0Dx();
// L^2 projector of y-derivative \Pi^0_{gradOrder}\partial_y
const Eigen::MatrixXd Pi0Dy = SVE.Pi0Dy();

```

### Three dimensions

This part is still in development.

Ideally the same three formulations will be available.
