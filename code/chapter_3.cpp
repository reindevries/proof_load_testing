// chapter_3.cpp
//
// This code iteratively calculates the target proof load factors for CC2 and CC3 in bending and
// shear. These are the factors provided in data file Chapter 3 Factors CC2 and CC3.csv and are
// plotted in Figure 3.3 of the dissertation.
//
// Author: Rein de Vries <rein.devries@tno.nl>
// Date: 10 November 2025

#include "rdv/msvcr.hpp"
#include "rdv/golden_section.hpp"
#include "rdv/random_variable.hpp"
#include "rdv/reliability.hpp"
#include "influence_line.hpp"
#include "load_effect.hpp"

void calc_target_loads(const std::vector<double>& m_Q, const std::vector<double>& V_Q, 
		std::vector<double>& Q_PL_CC2, std::vector<double>& Q_PL_CC3) {
	// Prepare result vectors
	ASSERT(size(m_Q) == size(V_Q));
	Q_PL_CC2.resize(size(m_Q));
	Q_PL_CC3.resize(size(m_Q));
	
	// Random variables
	auto rv_θ_Q_PL = rdv::rv_lognormal::mv(1.0, 0.10);
	auto rv_θ_Q = rdv::rv_lognormal::mv(1.0, 0.11);  // no dynamic amplification
	auto rv_C_0Q = rdv::rv_lognormal::mv(1.0, 0.10);
	
	for (size_t i = 0; i < size(m_Q); i++) {
		// Variable traffic load effect
		auto rv_Q = rdv::rv_gumbel_max::mv(m_Q[i], V_Q[i]);
		
		// Solve for consequence class 2
		rdv::golden_section_minimize([&](double Q_PL) {
			double beta;
			rdv::FORM(4, [&](auto u) {
				double θ_Q_PL = rv_θ_Q_PL.from_std_norm(u[0]);
				double θ_Q = rv_θ_Q.from_std_norm(u[1]);
				double C_0Q = rv_C_0Q.from_std_norm(u[2]);
				double Q = rv_Q.from_std_norm(u[3]);
				return θ_Q_PL * Q_PL - θ_Q * C_0Q * Q;
			}, beta);
			return rdv::sq(beta - 3.4);
		}, m_Q[i], 3.0 * m_Q[i], Q_PL_CC2[i]);
		
		// Solve for consequence class 3
		rdv::golden_section_minimize([&](double Q_PL) {
			double beta;
			rdv::FORM(4, [&](auto u) {
				double θ_Q_PL = rv_θ_Q_PL.from_std_norm(u[0]);
				double θ_Q = rv_θ_Q.from_std_norm(u[1]);
				double C_0Q = rv_C_0Q.from_std_norm(u[2]);
				double Q = rv_Q.from_std_norm(u[3]);
				return θ_Q_PL * Q_PL - θ_Q * C_0Q * Q;
			}, beta);
			return rdv::sq(beta - 4.0);
		}, m_Q[i], 3.0 * m_Q[i], Q_PL_CC3[i]);
	}
}

int main(int argc, char *argv[]) {
	// Initialisation
	std::cout.setf(std::ios::left);
	
	// Span lengths
	std::vector<double> L(40);
	for (size_t i = 0; i < size(L); i++)
		L[i] = (i + 1) * 5.0;
	
	// Load effects following from Eurocode LM, AASTHO HL20 and AASTHO HL93
	std::vector<double> M_Q_LM1, V_Q_LM1;
	std::vector<double> M_Q_HS20, V_Q_HS20;
	std::vector<double> M_Q_HL93, V_Q_HL93;
	load_effect_LM1(L, M_Q_LM1, V_Q_LM1);
	load_effect_HS20(L, M_Q_HS20, V_Q_HS20);
	load_effect_HL93(L, M_Q_HL93, V_Q_HL93);
	
	// Bending moment statistical description (WIM data analysis)
	std::vector<double> m_M_Q = {0.41, 1.15, 2.22, 3.56, 5.06, 6.58, 8.11, 9.63, 11.1, 
			12.7, 14.2, 15.7, 17.2, 18.9, 20.6, 22.3, 24.3, 26.5, 28.8, 31.3, 33.8, 36.4, 
			38.9, 41.5, 44.3, 47.0, 49.8, 52.7, 55.6, 58.4, 61.4, 64.5, 67.6, 70.8, 73.9, 
			77.1, 80.3, 83.6, 86.9, 90.3};
	std::vector<double> V_M_Q = {0.071, 0.023, 0.025, 0.018, 0.022, 0.023, 0.025, 0.026, 
			0.025, 0.025, 0.024, 0.023, 0.023, 0.025, 0.028, 0.031, 0.036, 0.044, 0.052, 
			0.060, 0.067, 0.071, 0.073, 0.076, 0.078, 0.080, 0.081, 0.082, 0.083, 0.084, 
			0.084, 0.085, 0.086, 0.086, 0.086, 0.086, 0.086, 0.085, 0.085, 0.085};
	
	// Shear force statistical description (WIM data analysis)
	std::vector<double> m_V_Q = {0.26, 0.42, 0.57, 0.71, 0.81, 0.87, 0.91, 0.95, 0.98, 
			1.01, 1.04, 1.07, 1.11, 1.17, 1.22, 1.27, 1.31, 1.34, 1.38, 1.44, 1.45, 1.47, 
			1.52, 1.56, 1.59, 1.62, 1.63, 1.68, 1.70, 1.72, 1.76, 1.79, 1.81, 1.85, 1.87, 
			1.89, 1.93, 1.93, 1.98, 2.02};
	std::vector<double> V_V_Q = {0.045, 0.039, 0.030, 0.030, 0.032, 0.032, 0.031, 0.031, 
			0.030, 0.030, 0.032, 0.038, 0.045, 0.052, 0.060, 0.064, 0.068, 0.069, 0.072, 
			0.074, 0.074, 0.073, 0.074, 0.077, 0.079, 0.078, 0.079, 0.078, 0.078, 0.078, 
			0.079, 0.081, 0.083, 0.085, 0.085, 0.087, 0.087, 0.089, 0.091, 0.092};
	
	// Calculate target loads
	std::vector<double> M_Q_PL_CC2, M_Q_PL_CC3;
	std::vector<double> V_Q_PL_CC2, V_Q_PL_CC3;
	calc_target_loads(m_M_Q, V_M_Q, M_Q_PL_CC2, M_Q_PL_CC3);
	calc_target_loads(m_V_Q, V_V_Q, V_Q_PL_CC2, V_Q_PL_CC3);
	
	// Print CC2 factors
	std::cout << "CC2 factors:\n";
	std::cout << std::setw(15) << "L" <<
				 std::setw(15) << "M LM1" <<
				 std::setw(15) << "M HS20" <<
				 std::setw(15) << "M HL93" <<
				 std::setw(15) << "V LM1" <<
				 std::setw(15) << "V HS20" <<
				 std::setw(15) << "V HL93" << "\n";
	for (size_t i = 0; i < size(L); i++)
		std::cout << std::setw(15) << L[i] << 
					 std::setw(15) << M_Q_PL_CC2[i] * 1e6 / M_Q_LM1[i] << 
					 std::setw(15) << M_Q_PL_CC2[i] * 1e6 / M_Q_HS20[i] << 
					 std::setw(15) << M_Q_PL_CC2[i] * 1e6 / M_Q_HL93[i] << 
					 std::setw(15) << V_Q_PL_CC2[i] * 1e6 / V_Q_LM1[i] << 
					 std::setw(15) << V_Q_PL_CC2[i] * 1e6 / V_Q_HS20[i] << 
					 std::setw(15) << V_Q_PL_CC2[i] * 1e6 / V_Q_HL93[i] << "\n";
	
	// Print CC3 factors
	std::cout << "CC3 factors:\n";
	std::cout << std::setw(15) << "L" <<
				 std::setw(15) << "M LM1" <<
				 std::setw(15) << "M HS20" <<
				 std::setw(15) << "M HL93" <<
				 std::setw(15) << "V LM1" <<
				 std::setw(15) << "V HS20" <<
				 std::setw(15) << "V HL93" << "\n";
	for (size_t i = 0; i < size(L); i++)
		std::cout << std::setw(15) << L[i] << 
					 std::setw(15) << M_Q_PL_CC3[i] * 1e6 / M_Q_LM1[i] << 
					 std::setw(15) << M_Q_PL_CC3[i] * 1e6 / M_Q_HS20[i] << 
					 std::setw(15) << M_Q_PL_CC3[i] * 1e6 / M_Q_HL93[i] << 
					 std::setw(15) << V_Q_PL_CC3[i] * 1e6 / V_Q_LM1[i] << 
					 std::setw(15) << V_Q_PL_CC3[i] * 1e6 / V_Q_HS20[i] << 
					 std::setw(15) << V_Q_PL_CC3[i] * 1e6 / V_Q_HL93[i] << "\n";
	
	return 0;
}
