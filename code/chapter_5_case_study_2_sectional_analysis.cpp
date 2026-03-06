// chapter_5_case_study_2_sectional_analysis.cpp
//
// Since for the second case study no laboratory data is available, a numerical simulation using
// Latin Hypercube Sampling (LHS) is performed. The output of this code has resulted in the data of
// Chapter 5 Sectional analysis.csv and is displayed in Figure 5.6 of the dissertation.
//
// Author: Rein de Vries <rein.devries@tno.nl>
// Date: 10 November 2025

#include "rdv/msvcr.hpp"
#include "rdv/interpolate.hpp"
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "rdv/random_variable.hpp"
#include "rdv/sampling.hpp"
#include "rdv/vector_ext.hpp"
#include <fstream>
#include <omp.h>

int main(int argc, char* argv[]) {
	// Initialisation
	rdv::tic();
	std::cout.setf(std::ios::left);
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
	
	// Geometry	
	double b = 9.4;                                  // width of cross-section
	auto rv_h = rdv::rv_lognormal::mv(0.470, 0.10);  // height of cross-section
	auto rv_c = rdv::rv_gamma::mv(0.030, 0.17);      // concrete cover
	
	// Reinforcement
	double A_s_bot1 = M_PI * rdv::sq(0.025 / 2.0) * b / 0.140;  // first layer, the very bottom
	double A_s_bot2 = M_PI * rdv::sq(0.025 / 2.0) * b / 0.280;  // second layer, closer to middle
	double A_s_top = M_PI * rdv::sq(0.025 / 2.0)  * b / 0.280;  // top layer
	
	// Yield stress of steel reinforcement
	std::vector<double> x_f_y = {279.1e6, 295.2e6, 297.3e6};
	auto rv_ln_f_y = rdv::rv_students_t(size(x_f_y) - 1, 
			mean(log(x_f_y)), stddev_s(log(x_f_y)) * sqrt(1.0 + 1.0 / size(x_f_y)));
	
	// Young's modulus of steel reinforcement
	double E_sm = 205e9;
	auto rv_E_s = rdv::rv_lognormal::mv(E_sm, 0.02);
	
	// Steel stress-strain curve parameters
	double f_ym = mean(x_f_y);
	double ε_ym = f_ym / E_sm;
	double ε_y_to_ε_y2 = 0.01 / ε_ym;
	double ε_y_to_ε_u = 0.2 / ε_ym;
	double f_y_to_f_u = 420.0e6 / f_ym;
	
	// Concrete compressive strength
	std::vector<double> x_f_c = {61.3e6, 68.2e6, 74.0e6, 75.6e6, 62.5e6, 79.4e6};
	x_f_c = 0.82 * x_f_c;  // cube to cilinder strength
	auto rv_ln_f_c = rdv::rv_students_t(size(x_f_c) - 1, 
			mean(log(x_f_c)), stddev_s(log(x_f_c)) * sqrt(1.0 + 1.0 / size(x_f_c)));
	
	// Non-linearity
	double c_1 = 0.75;
	auto rv_c_2 = rdv::rv_lognormal::mv(0.41, 0.1);
	double c_3 = 4.1;
	
	// Output calculated values
	std::cout << "A_s_bot1: " << A_s_bot1 * 1e6 << " mm²\n";
	std::cout << "A_s_bot2: " << A_s_bot2 * 1e6 << " mm²\n";
	std::cout << "A_s_top: " << A_s_top * 1e6 << " mm²\n";
	std::cout << "m_f_c: " << mean(x_f_c) / 1e6 << " MPa\n";
	std::cout << "V_f_c: " << cov_s(x_f_c) << "\n";
	std::cout << "m_f_y: " << mean(x_f_y) / 1e6 << " MPa\n";
	std::cout << "V_f_y: " << cov_s(x_f_y) << "\n";
	std::cout << "ε_ym = f_ym / E_sm: " << ε_ym * 1e6 << " · 10⁻⁶\n";
	
	// Create random latin hypercube
	std::cout << "Creating latin hypercube...\n";
	size_t n_samples = 100;
	std::vector<std::vector<double>> samples;
	rdv::latin_hyp_sampling(6, n_samples, samples);
	
	// Store yield moments
	std::vector<double> M_ys(n_samples);
	
	// Store moment-strain curves
	size_t curve_points = 101;
	double curve_ε_0 = 0.0;
	double curve_dε = 25e-6;
	std::vector<std::vector<double>> curves;
	resize(curves, n_samples, curve_points);
	
	// Loop over samples
	std::cout << "Calculating curves...\n";
	for (size_t i_sample = 0; i_sample < n_samples; i_sample++) {
		const std::vector<double>& sample = samples[i_sample];
		
		// Obtain realizations of random variables
		double h = rv_h.quantile(sample[0]);
		double c = rv_c.quantile(sample[1]);
		double f_y = exp(rv_ln_f_y.quantile(sample[2]));
		double E_s = rv_E_s.quantile(sample[3]);
		double f_c = exp(rv_ln_f_c.quantile(sample[4]));
		double c_2 = rv_c_2.quantile(sample[5]);
	
		// Derived properties
		double E_c = 21.5e9 * pow(f_c / 1e7, 1.0 / 3.0);
		double d_bot1 = h - c - 0.025 / 2.0;
		double d_bot2 = h - c - 0.025 - 0.014 - 0.025 / 2.0;
		double d_top = h - c - 0.025 / 2.0;
		
		// Concrete stress-strain curve (Thorenfeldt et al., 1987)
		double n_th = 0.8 + f_c / 17.24e6;
		double ε_0 = f_c / E_c * n_th / (n_th - 1.0);
		auto σ_c_th = [&](double ε_c) {
			if (ε_c < 0.0)
				return 0.0;
			double k_th = (ε_c < ε_0 ? 1.0 : 0.67 + f_c / 62.07e6);
			return f_c * n_th * ε_c / ε_0 / (n_th - 1.0 + pow(ε_c / ε_0, n_th * k_th));
		};
		
		// Steel stress-strain curve
		//
		//      |                      ____
		//      |    ______....----````    |
		//      |   /                      |
		//      |  /                       |
		//      | /                        |
		// _____ /_________________________|__
		//      /  ε_y    ε_y2            ε_u
		//     /|
		//    / |
		//   /  |
		//  /   |
		
		// Calculate points of curve
		double ε_y = f_y / E_s;
		double ε_y2 = ε_y_to_ε_y2 * ε_y;
		double ε_u = ε_y_to_ε_u * ε_y;
		double f_u = f_y_to_f_u * f_y;
		std::vector<double> steel_ε = {-1.0, 0.0, ε_y, ε_y2, ε_u, ε_u + 1e-6, 1.0};
		std::vector<double> steel_σ = {-E_s, 0.0, f_y, f_y,  f_u,        0.0, 0.0};
		
		// Function for integration of cross-section
		auto integrate = [&](double ε_c_top, double x, double& H, double& M) {
			// Calculate axial force and moment
			double κ = ε_c_top / x;
			H = 0.0;
			M = 0.0;
			
			// Concrete 'slices'
			int N = 100;
			double dz = h / N;
			
			for (int i = 0; i < N; i++) {
				double z = (i + 0.5) * dz;
				double ε_c = ε_c_top - κ * z;
				double F_c = -σ_c_th(ε_c) * dz * b;  // negative, because tensile force is positive
				
				H += F_c;
				M += F_c * z;
			}
			
			// Bottom steel
			{
				double z = d_bot1;
				double ε_c = ε_c_top - κ * z;
				double ε = -ε_c;
				
				double σ_s = rdv::interpolate(steel_ε, steel_σ, ε, rdv::extrapolation::linear);
				double F_s = σ_s * A_s_bot1;
				double F_c = -σ_c_th(ε_c) * A_s_bot1;  // negative, because tension is positive
				
				H += F_s;
				H -= F_c;  // subtract its earlier contribution (if any)
				
				M += F_s * z;
				M -= F_c * z;  // subtract its earlier contribution (if any)
			}
			{
				double z = d_bot2;
				double ε_c = ε_c_top - κ * z;
				double ε = -ε_c;
				
				double σ_s = rdv::interpolate(steel_ε, steel_σ, ε, rdv::extrapolation::linear);
				double F_s = σ_s * A_s_bot2;
				double F_c = -σ_c_th(ε_c) * A_s_bot2;  // negative, because tension is positive
				
				H += F_s;
				H -= F_c;  // subtract its earlier contribution (if any)
				
				M += F_s * z;
				M -= F_c * z;  // subtract its earlier contribution (if any)
			}
			
			// Top steel, account for its compressive contribution
			{
				double z = h - d_top;
				double ε_c = ε_c_top - κ * z;
				double ε = -ε_c;
				
				double σ_s = rdv::interpolate(steel_ε, steel_σ, ε, rdv::extrapolation::linear);
				double F_s = σ_s * A_s_top;
				double F_c = -σ_c_th(ε_c) * A_s_top;  // negative, because tension is positive
				
				H += F_s;
				H -= F_c;  // subtract its earlier contribution (if any)
				
				M += F_s * z;
				M -= F_c * z;  // subtract its earlier contribution (if any)
			}
		};
		
		// Cross-section analysis
		std::vector<double> εs;
		std::vector<double> Ms;
		εs.reserve(141);
		Ms.reserve(141);
		
		εs.push_back(0.0);
		Ms.push_back(0.0);
		
		for (double ε_c_top = 0.000025; ε_c_top < 0.0035001; ε_c_top += 0.000025) {
			// Bisection method; find x where horizontal forces are balanced (H = 0)
			double a = 0.0;
			double b = h;
					
			for (int i = 0; b - a > 1e-6 && i < 100; i++) {
				double x = (a + b) / 2.0;
				double H, M;
				integrate(ε_c_top, x, H, M);
				
				if (H > 0.0)
					a = x;
				else
					b = x;
			}
			
			double x = (a + b) / 2.0;
			
			// Calculate bottom strain
			double κ = ε_c_top / x;
			double ε_s_bot1 = -ε_c_top + κ * d_bot1;
			
			// Obtain corresponding moment
			double H, M;
			integrate(ε_c_top, x, H, M);
			
			// Include the strain modification and store point
			double ε = ε_s_bot1 * (c_1 + c_2 * pow(ε_s_bot1 / ε_y, c_3));
			εs.push_back(ε);
			Ms.push_back(M);
		}
		
		// Find yield moment
		M_ys[i_sample] = rdv::interpolate(εs, Ms, ε_y);
		
		// Interpolate for the final curves
		for (size_t i_point = 0; i_point < curve_points; i_point++) {
			double curve_ε = curve_ε_0 + curve_dε * i_point;
			curves[i_sample][i_point] = rdv::interpolate(εs, Ms, curve_ε);
		}
	}
	
	// Output the curves to file
	std::string file_name = "Moment-strain curves.csv";
	std::cout << "Writing result to: " << file_name << "\n";
	
	std::ofstream output_csv(file_name);
	output_csv << "M_y";
	for (size_t i_sample = 0; i_sample < size(samples); i_sample++)
		output_csv << "," << M_ys[i_sample] / 1e3;
	output_csv << "\n";
	
	for (size_t i_point = 0; i_point < curve_points; i_point++) {
		double curve_ε = curve_ε_0 + curve_dε * i_point;
		output_csv << curve_ε * 1e6;
		for (size_t i_sample = 0; i_sample < size(samples); i_sample++)
			output_csv << "," << curves[i_sample][i_point] / 1e3;
		output_csv << "\n";
	}
	
	// Print elapsed time
	rdv::toc();
	return 0;
}
