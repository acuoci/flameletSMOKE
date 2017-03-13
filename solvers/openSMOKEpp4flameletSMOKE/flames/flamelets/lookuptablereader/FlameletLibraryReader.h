/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2016 Alberto Cuoci                                       |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/
 
#ifndef OpenSMOKE_FlameletLibraryReader_H
#define OpenSMOKE_FlameletLibraryReader_H

#include "FlameletReader.h"
#include "flames/flamelets/lookuptablereader/logNormalDistribution.h"

namespace OpenSMOKE
{
	//!  A class to manage flamelets with same enthalpy defect, but different scalar dissipation rates
	/*!
	This class provides the tools to manage flamelets with same enthalpy defect, but different scalar dissipation rates
	*/

	class FlameletLibraryReader
	{

	public:

		enum pdfScalarDissipationRate { X_PDF_DIRAC_DELTA, X_PDF_LOG_NORMAL };

	public:

		/**
		*@brief Default constructor
		*/
		FlameletLibraryReader();

		/**
		*@brief Reads from the XML files the set of flamelets with same enthalpy defect, but different scalar dissipation rates
		*@param list_of_flamelets list of XML files corresponding to the different flamelets
		*/
		void Read(const std::vector<boost::filesystem::path> list_of_flamelets);

		/**
		*@brief Prints a summary on the screen
		*/
		void Summary();

		/**
		*@brief Returns the mean values
		*@param csi current mean mixture fraction
		*@param csiv2 current mean mixture fraction variance (normalized)
		*@param chi_st current stoichiometric scalar dissipation rate [1/s]
		*@param extracted returned mean values (returned)
		*/
		void GetMeanValues(const double csi, const double csiv2, const double chi_st, std::vector<double>& extracted);

		/**
		*@brief Returns the mean values of mass fractions of species
		*@param csi current mean mixture fraction
		*@param csiv2 current mean mixture fraction variance (normalized)
		*@param chi_st current stoichiometric scalar dissipation rate [1/s]
		*@param &omegaFavre returned mass fractions of species
		*/
		void ExtractMeanValues(const double csi, const double csiv2, const double chi_st, std::vector<double> &omegaFavre);

		/**
		*@brief Sets the log-normal distribution function (default is Delta Dirac)
		*@param sigma the scale factor of the distribution (default is 1.31)
		*/
		void SetLogNormalChiDistribution(const double sigma);

		/**
		*@brief Sets the list of species to be updated on the fly
		*@param names list of species to be updated on the fly
		*/
		void SetSpeciesToExtract(const std::vector<std::string> names);

		/**
		*@brief To enable a short summary on the screen
		*/
		void SetShowFlamelet();

		/**
		*@brief Returns the density of the fuel stream [kg/m3]
		*/
		inline double density_r_fuel() const;

		/**
		*@brief Returns the density of the oxidizer stream [kg/m3]
		*/
		inline double density_r_oxidizer() const;

		/**
		*@brief Returns the enthalpy of the fuel stream [J/kg]
		*/
		inline double enthalpy_f_fuel() const;

		/**
		*@brief Returns the enthalpy of the oxidizer stream [J/kg]
		*/
		inline double enthalpy_f_oxidizer() const;

		/**
		*@brief Returns the temperature of the fuel stream [K]
		*/
		inline double temperature_f_fuel() const;

		/**
		*@brief Returns the temperature of the oxidizer stream [K]
		*/
		inline double temperature_f_oxidizer() const;

	private:

		pdfScalarDissipationRate	 X_pdf_type_;					//!< type of distribution function for scalar dissipation rate

		unsigned int _n;											//!< number of flamelets

		std::vector<std::string>	_names_of_species_to_extract;	//!< list of species to be extracted on the fly


		std::vector<double> _chi_st;		//!< stoichiometric scalar dissipation rates (0-index based) [1/s]
		std::vector<double> _enthalpy_defect;	//!< enthalpy defects (1-index based) [J/kg]
		std::vector< FlameletReader >	_flame;	//!< vector of flamelets (1-index based)
		

		bool   _iMultiScalarDissipationRates;	//!< if true, more than one scalar dissipation rate is available
		bool 	_showFlamelet;			//!< if true, a short summary is printed on the screen

		double _temperature_f_max;		//!< maximum temperature [K]
		double _temperature_f_min;		//!< minimum temperature [K]
		double _density_r_max;			//!< maximum density [kg/m3]
		double _density_r_min;			//!< minimum density [kg/m3] 
			
		double _chi_st_min;			//!< minimum scalar dissipation rate available [1/s]
		double _chi_st_max;			//!< minimum scalar dissipation rate available [1/s]

		double _temperature_f_fuel;		//!< temperature of fuel stream [K]
		double _temperature_f_oxidizer;		//!< temperature of oxidizer stream [K]
		double _enthalpy_f_fuel;		//!< enthalpy of fuel stream [J/kg]
		double _enthalpy_f_oxidizer;		//!< enthalpy of oxidizer stream [J/kg]
		double _density_r_fuel;			//!< density of fuel stream [kg/m3]
		double _density_r_oxidizer;		//!< density of oxidizer stream [kg/m3]

	// Data reserved to log-normal distribution
	private:

		OpenSMOKE::logNormalDistribution	logNormalPDF_;		//!< log-normal distribution
		double 					_chi_log_normal_sigma;	//!< scalae factor for the log-normal distribution

		std::vector<double>			_t_favre;		//!< auxiliary vector to be used when the log-normal distribution is on
		std::vector<double>			_rho_reynolds;		//!< auxiliary vector to be used when the log-normal distribution is on
		std::vector< std::vector<double> >	_w_favre;		//!< auxiliary vector to be used when the log-normal distribution is on
		std::vector<double>			_cp_favre;		//!< auxiliary vector to be used when the log-normal distribution is on
		std::vector<double>			_mu_favre;		//!< auxiliary vector to be used when the log-normal distribution is on
		std::vector<double>			_as_favre;		//!< auxiliary vector to be used when the log-normal distribution is on
		std::vector<double>			_alpha_favre;		//!< auxiliary vector to be used when the log-normal distribution is on

	private:

		/**
		*@brief Returns the mean values when the delta-dirac distribution is adopted
		*@param csi current mean mixture fraction
		*@param csiv2 current mean mixture fraction variance (normalized)
		*@param chi_st current stoichiometric scalar dissipation rate [1/s]
		*@param extracted returned mean values (returned)
		*/
		void GetMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, std::vector<double> &extracted);

		/**
		*@brief Returns the mean values when the log-normal distribution is adopted
		*@param csi current mean mixture fraction
		*@param csiv2 current mean mixture fraction variance (normalized)
		*@param chi_st current stoichiometric scalar dissipation rate [1/s]
		*@param extracted returned mean values (returned)
		*/
		void GetMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, std::vector<double> &extracted);

		/**
		*@brief Returns the mean values of mass fractions of species when the delta-dirac distribution is adopted
		*@param csi current mean mixture fraction
		*@param csiv2 current mean mixture fraction variance (normalized)
		*@param chi_st current stoichiometric scalar dissipation rate [1/s]
		*@param &omegaFavre returned mass fractions of species
		*/
		void ExtractMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, std::vector<double> &omegaFavre);

		/**
		*@brief Returns the mean values of mass fractions of species when the log-normal distribution is adopted
		*@param csi current mean mixture fraction
		*@param csiv2 current mean mixture fraction variance (normalized)
		*@param chi_st current stoichiometric scalar dissipation rate [1/s]
		*@param &omegaFavre returned mass fractions of species
		*/
		void ExtractMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, std::vector<double> &omegaFavre);

		/**
		*@brief Returns a fatal error message on the screen and stops the simulation
		*@param error_message the returned error message
		*/
		void ErrorMessage(const std::string error_message);

		/**
		*@brief Returns a warning message on the screen
		*@param warning_message the returned warning message
		*/
		void WarningMessage(const std::string warning_message);
	};
}

#include "FlameletLibraryReader.hpp"

#endif // OpenSMOKE_FlameletLibraryReader_H
	
