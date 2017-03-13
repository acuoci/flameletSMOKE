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
|	License                                                           |
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
 
#ifndef OpenSMOKE_NonAdiabaticFlameletLibraryReader_H
#define OpenSMOKE_NonAdiabaticFlameletLibraryReader_H

#include "FlameletLibraryReader.h"
#include "LookUpTablesUtilities.h"

namespace OpenSMOKE
{

	//!  A class to manage the whole non-adiabatic flamelet library 
	/*!
	This class provides the tools to manage the whole non-adiabatic flamelet library
	*/

	class NonAdiabaticFlameletLibraryReader
	{

	public:

		/**
		*@brief Default constructor
		*/
		NonAdiabaticFlameletLibraryReader();

		/**
		*@brief Reads the whole lookup table from the location indicated by the folder_name_ path
		*/
		void Read();

		/**
		*@brief Reads the whole lookup table from the specified path
		*@param folder_name the path to the folder containing the whole library
		*/
		void Read(const boost::filesystem::path folder_name);

		/**
		*@brief Prints a short summary on the screen
		*/
		void Summary();

		/**
		*@brief Returns the mean values
		*@param csi current mean mixture fraction
		*@param csiv2 current mean mixture fraction variance (normalized)
		*@param chi_st current stoichiometric scalar dissipation rate [1/s]
		*@param phi the current enthalpy defect [J/kg]
		*@param extracted returned mean values (returned)
		*/
		void   GetMeanValues(const double csi, const double csiv2, const double chi_st, const double phi, std::vector<double>& extracted);

		/**
		*@brief Returns the mean values of mass fractions of species
		*@param csi current mean mixture fraction
		*@param csiv2 current mean mixture fraction variance (normalized)
		*@param chi_st current stoichiometric scalar dissipation rate [1/s]
		*@param phi the current enthalpy defect [J/kg]
		*@param &omegaFavre returned mass fractions of species
		*/
		void   ExtractMeanValues(const double csi, const double csiv2, const double chi_st, const double phi, std::vector<double> &omegaFavre);

		/**
		*@brief Reconstructs the enthalpy defect from the mixture fraction (i.e. composition) and temperature
		*@param csi current mean mixture fraction
		*@param csiv2 current mean mixture fraction variance (normalized)
		*@param chi_st current stoichiometric scalar dissipation rate [1/s]
		*@param T the current temperature [K]
		*/
		double GetEnthalpyDefectFromTemperature(const double csi, const double csiv2, const double chi_st, const double T);

		/**
		*@brief Sets the path to the folder containing the whole lookup table
		*@param folder_name the path to the folder containing the whole lookup table
		*/
		void SetLibraryFolder(const boost::filesystem::path folder_name);

		/**
		*@brief Sets the log-normal distribution fo the scalar dissipation rate (default is delta-dirac)
		*@param sigma the scale factor in the log-normal distribution (default 1.31)
		*/
		void SetLogNormalChiDistribution(const double sigma);

		/**
		*@brief Sets on the adiabatic mode (even if the lookup table is not adiabatic)
		*/
		void SetAdiabaticMode();

		/**
		*@brief If true, a short summary is printed on the screen for each flamelet
		*/
		void SetShowFlamelet();

		/**
		*@brief If true, a short summary is printed on the screen for each flamelet library
		*/
		void SetShowFlameletLibrary();

		/**
		*@brief Sets the list of species to be extracted on the fly
		*@param names list of species to be extracted on the fly
		*/
		void SetSpeciesToExtract(const std::vector<std::string> names);

		/**
		*@brief Sets the list of species to be extracted on the fly
		*@param names list of species to be extracted on the fly
		*/
		void SetSpeciesToExtract(const std::string list_of_names);

		/**
		*@brief Returns the number of species to be extracted on the fly
		*/
		int number_of_species();

		/**
		*@brief Returns the species to be extracted on the fly
		*/
		const std::vector<std::string> species();

		/**
		*@brief Returns the index of species corresponding to the given name
		*@param name name of species for which the index is requested
		*/
		int index_of_species(const std::string name);

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

		boost::filesystem::path 				folder_name_;	//!< path to the folder containing the whole lookup table
		FlameletLibraryReader::pdfScalarDissipationRate 	X_pdf_type_;	//!< type of distribution function for the scalar dissipation rate

		std::vector<FlameletLibraryReader> flamelet_libraries;		//!< vector of flamelet libraries

		unsigned int	   			_nphi;							//!< number of enthalpy defects
		std::vector<std::string>	names_of_species_to_extract_;	//!< list of species to be extracted on the fly
		std::vector<double>			_phi;							//!< list of enthalpy defects [J/kg]

	private:

		bool _adiabatic_mode;			//!< if true, calculations are carried out in adiabatic conditions
		bool _showFlamelet;				//!< if true, a short summary is printed on the screen for each flamelet
		bool _showFlameletLibrary;		//!< if true, a short summary is printed on the screen for each flamelet library

		double _temperature_f_fuel;		//!< temperature of fuel stream [K]
		double _temperature_f_oxidizer;		//!< temperature of oxidizer stream [K]
		double _enthalpy_f_fuel;		//!< enthalpy of fuel stream [J/kg]
		double _enthalpy_f_oxidizer;		//!< enthalpy of oxidizer stream [J/kg]
		double _density_r_fuel;			//!< density of fuel stream [kg/m3]
		double _density_r_oxidizer;		//!< density of oxidizer stream [kg/m3] 
		double _phi_max;			//!< maximum enthalpy defect [J/kg]
		double _phi_min;			//!< minimum enthalpy defect [J/kg] 
		int    _jPhiAdiabatic;			//!< index of adiabatic flamelet library

		double 	_chi_log_normal_sigma;		//!< scale factor for the log-normal distribution

	private:

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

#include "NonAdiabaticFlameletLibraryReader.hpp"

#endif // OpenSMOKE_NonAdiabaticFlameletLibraryReader_H
	
