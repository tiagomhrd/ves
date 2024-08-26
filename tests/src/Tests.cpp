#include <iostream>
#include <catch_amalgamated.hpp>
#include <corecrt_math_defines.h>
#include "ves.h"
#include "mnl/include/mnl.hpp"


const double tol = 1e-8;
TEST_CASE("HHGTVEM Examples") {
    SECTION("Square"){
		std::vector<Eigen::Vector2d> poly;
		poly.reserve(4);
		poly.emplace_back(0., 0.);
		poly.emplace_back(1., 0.);
		poly.emplace_back(1., 1.);
		poly.emplace_back(0., 1.);

		const double d = sqrt(2.);

		SECTION("k=1"){
			const ves::V2D VE(poly, 1);

			SECTION("B"){
				Eigen::MatrixXd refB(3, 4);
				refB << 1, 1, 1, 1,
						-d, d, d,-d,
						-d,-d, d, d;
				refB *= .25;

				const Eigen::MatrixXd B = VE.BGrad();
				REQUIRE_THAT((B - refB).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}

			SECTION("D"){
				Eigen::MatrixXd refD(4, 3);
				refD << 4.,-d,-d,
						4., d,-d,
						4., d, d,
						4.,-d, d; 
				refD *= .25;

				const Eigen::MatrixXd D = VE.D();
				REQUIRE_THAT((D - refD).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}

			SECTION("G"){
				Eigen::MatrixXd refG(3, 3);
				refG << 2., 0., 0.,
						0., 1., 0.,
						0., 0., 1.;
				refG *= .5;

				const Eigen::MatrixXd G = VE.GGrad();
				REQUIRE_THAT((G - refG).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}
		}
		SECTION("k=2"){
			const ves::V2D VE(poly, 2);

			SECTION("B"){
				Eigen::MatrixXd refB(6, 9);
				refB <<  0, 0, 0, 0,   0,   0,   0,   0, 12,
						-d, d, d,-d,   0, 4*d,   0,-4*d,  0,
						-d,-d, d, d,-4*d,   0, 4*d,   0,  0,
						 1, 1, 1, 1,   0,   4,   0,   4,-12,
						 1,-1, 1,-1,   0,   0,   0,   0,  0,
						 1, 1, 1, 1,   4,   0,   4,   0,-12;
				refB /= 12.;

				const Eigen::MatrixXd B = VE.BGrad();
				REQUIRE_THAT((B - refB).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}

			SECTION("D"){
				Eigen::MatrixXd refD(9, 6);
				refD << 24,-6.*d,-6.*d, 3, 3, 3,
						24, 6.*d,-6.*d, 3,-3, 3,
						24, 6.*d, 6.*d, 3, 3, 3,
						24,-6.*d, 6.*d, 3,-3, 3,
						24,    0,-6.*d, 0, 0, 3,
						24, 6.*d,    0, 3, 0, 0,
						24,    0, 6.*d, 0, 0, 3,
						24,-6.*d,    0, 3, 0, 0,
						24,    0,    0, 1, 0, 1; 
				refD /= 24.;

				const Eigen::MatrixXd D = VE.D();
				REQUIRE_THAT((D - refD).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}

			SECTION("G"){
				Eigen::MatrixXd refG(6, 6);
				refG << 24,  0,  0, 1, 0, 1,
						 0, 12,  0, 0, 0, 0,
						 0,  0, 12, 0, 0, 0,
						 0,  0,  0, 2, 0, 0,
						 0,  0,  0, 0, 1, 0,
						 0,  0,  0, 0, 0, 2;
				refG /= 24.;

				const Eigen::MatrixXd G = VE.GGrad();
				REQUIRE_THAT((G - refG).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}
		}
	}
	SECTION("Pentagon"){
		std::vector<Eigen::Vector2d> poly;
		poly.reserve(5);
		poly.emplace_back(0., 0.);
		poly.emplace_back(3., 0.);
		poly.emplace_back(3., 2.);
		poly.emplace_back(1.5, 4.);
		poly.emplace_back(0., 4.);

		SECTION("k=1"){
			const ves::V2D VE(poly, 1);

			SECTION("B"){
				Eigen::MatrixXd refB(3, 5);
				refB << 4, 4, 4, 4, 4,
					   -8, 4, 8, 4,-8,
					   -6,-6, 3, 6, 3;
				refB /= 20.;

				const Eigen::MatrixXd B = VE.BGrad();
				REQUIRE_THAT((B - refB).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}

			SECTION("D"){
				Eigen::MatrixXd refD(5, 3);
				refD << 1470,-399,-532,
						1470, 483,-532,
						1470, 483,  56,
						1470,  42, 644,
						1470,-399, 644;
				refD /= 1470.;

				const Eigen::MatrixXd D = VE.D();
				REQUIRE_THAT((D - refD).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}

			SECTION("G"){
				Eigen::MatrixXd refG(3, 3);
				refG << 1050,  30,  40,
						   0, 441,   0,
						   0,   0, 441;
				refG /= 1050.;

				const Eigen::MatrixXd G = VE.GGrad();
				REQUIRE_THAT((G - refG).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}
		}
		SECTION("k=2"){
			const ves::V2D VE(poly, 2);

			SECTION("B"){
				Eigen::MatrixXd refB(6, 11);
				refB <<      0,     0,     0,     0,      0,     0,     0,     0,     0,      0,  88200,
						-11760,  5880, 11760,  5880, -11760,     0, 23520, 23520,     0, -47040,      0,
					     -8820, -8820,  4410,  8820,   4410,-35280,     0, 17640, 17640,      0,      0,
						  6384,  3864,  7728,   336,   6384,     0, 15456,  8400,     0,  25536, -74088,
						  6650, -5026,  1897,  2828,  -6349, -1008, -3808,  8750, -2142,  -1792,      0,
						  6384,  6384,   336,  7728,   3864, 25536,     0,  8400, 15456,      0, -74088;
				refB /= 88200.;

				const Eigen::MatrixXd B = VE.BGrad();
				REQUIRE_THAT((B - refB).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}

			SECTION("D"){
				Eigen::MatrixXd refD(11, 6);
				refD << 176400, -47880, -63840, 12996,  17328, 23104,
						176400,  57960, -63840, 19044, -20976, 23104,
						176400,  57960,   6720, 19044,   2208,   256,
						176400,   5040,  77280,   144,   2208, 33856,
						176400, -47880,  77280, 12996, -20976, 33856,
						176400,   5040, -63840,   144,  -1824, 23104,
						176400,  57960, -28560, 19044,  -9384,  4624,
						176400,  31500,  42000,  5625,   7500, 10000,
						176400, -21420,  77280,  2601,  -9384, 33856,
						176400, -47880,   6720, 12996,  -1824,   256,
						176400,      0,      0,  4770,  -1452,  8480;
				refD /= 176400.;

				const Eigen::MatrixXd D = VE.D();
				REQUIRE_THAT((D - refD).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}

			SECTION("G"){
				Eigen::MatrixXd refG(6, 6);
				refG << 4410000,       0,       0, 119250, -36300, 212000,
						      0, 1852200,       0,      0,      0,      0, 
						      0,       0, 1852200,      0,      0,      0,
							  0,       0,       0, 200340, -30492,      0,
							  0,       0,       0, -30492, 139125, -30492,
							  0,       0,       0,      0, -30492, 356160;
				refG /= 4410000.;

				const Eigen::MatrixXd G = VE.GGrad();
				const Eigen::MatrixXd difG = (G - refG);
				REQUIRE_THAT((G - refG).norm(), Catch::Matchers::WithinAbs(0.0, tol));
			}
		}
	}
}

/* Auxiliar function to generate regular polygons */
static std::vector<Eigen::Vector2d> regularPolygon(const int nv) {
	std::vector<Eigen::Vector2d> polygon;
	polygon.reserve(nv);
	const double halfTheta = M_PI / nv;
	const double r = pow(2. * sin(halfTheta), -1.);
	for (int i{}; i < nv; ++i)
		polygon.emplace_back(r * cos(2 * i * halfTheta), r * sin(2 * i * halfTheta));
	return polygon;
}

TEST_CASE("Polynomial Endomorphism") {
	std::stringstream ss;
	const int maxNumberSides = 10;
	const int maxOrder = 4;
	for (int nv = 3; nv <= maxNumberSides; ++nv) {
		ss << nv << "-gon"; 
		SECTION(ss.str()) {
			std::stringstream vss;
			auto poly = regularPolygon(nv);
			for (int k = 1; k <= maxOrder; ++k){
				vss << "V2D k=" << k;
				const int nk = mnl::PSpace2D::SpaceDim(k);
				SECTION(vss.str()){
					ves::V2D VE(poly, k);
					const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nk, nk);
					const Eigen::MatrixXd D = VE.D();

					const Eigen::MatrixXd PiGrad = VE.PiGrad();
					REQUIRE_THAT((PiGrad * D - I).norm(),Catch::Matchers::WithinAbs(0.0, tol)); 
					
					const Eigen::MatrixXd Pi0 = VE.Pi0();
					REQUIRE_THAT((Pi0 * D - I).norm(),Catch::Matchers::WithinAbs(0.0, tol)); 
				}
				vss.str("");
				vss << "SV2D k=" << k;
				SECTION(vss.str()){
					ves::SV2D VE(poly, k);
					const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nk, nk);
					const Eigen::MatrixXd D = VE.D();

					const Eigen::MatrixXd Pi0 = VE.Pi0();
					REQUIRE_THAT((Pi0 * D - I).norm(),Catch::Matchers::WithinAbs(0.0, tol)); 
				}
				vss.str("");
				vss << "SE2V2D k=" << k;
				SECTION(vss.str()) {
					ves::SE2V2D VE(poly, k);
					const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nk, nk);
					const Eigen::MatrixXd D = VE.D().leftCols(nk);

					const Eigen::MatrixXd Pi0 = VE.Pi0();
					REQUIRE_THAT((Pi0 * D - I).norm(), Catch::Matchers::WithinAbs(0.0, tol));
				}
			}
		}
		ss.str("");
	}
}


int main(int argc, char* argv[]) {
	int result = Catch::Session().run(argc, argv);
	return result;
}