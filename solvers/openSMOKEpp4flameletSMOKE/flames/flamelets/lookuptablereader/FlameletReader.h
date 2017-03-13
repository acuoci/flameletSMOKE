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
 
#ifndef OpenSMOKE_FlameletReader_H
#define OpenSMOKE_FlameletReader_H

#include "OpenSMOKEpp"
#include "math/OpenSMOKEUtilities.h"
#include <boost/filesystem.hpp>

namespace OpenSMOKE
{
	//!  A class to manage pre-processed flamelet solutions
	/*!
	This class provides the tools to manage pre-processed flamelet solutions
	*/

	class FlameletReader
	{

	public:

		/**
		*@brief Default constructor
		*/
		FlameletReader();

		/**
		*@brief Read the pre-processed flamelet solution from the corresponding XML file
		*@param file_name the XML file containing the flamelet solution
		*/
		void Read(const boost::filesystem::path file_name);

		/**
		*@brief Prints on the screen a summary about the flamelet solution
		*/
		void Summary();

		/**
		*@brief Extracts from the pre-processed flamelet solution the relevant averaged data
		*@param csi mixture fraction
		*@param csiv2 bormalized mixture fraction variance
		*@param extracted extracted data
		*/
		void GetMeanValues(const double csi, const double csiv2, std::vector<double>& extracted);

		/**
		*@brief Extracts from the pre-processed flamelet solution the mass fractions of selected species
		*@param csi mixture fraction
		*@param csiv2 bormalized mixture fraction variance
		*@param omegaFavre the extracted mass fraction of selected species
		*/
		void ExtractMeanValues(const double csi, const double csiv2, std::vector<double> &omegaFavre);

		/**
		*@brief Sets the names of species to be extracted
		*@param names list of species to be extracted
		*/
		void SetSpeciesToExtract(const std::vector<std::string> names);

		/**
		*@brief Returns the maximum temperature (favre) [K]
		*/
		double temperature_f_max() const;

		/**
		*@brief Returns the minimum temperature (favre) [K]
		*/
		double temperature_f_min() const;

		/**
		*@brief Returns the maximum density (reynolds) [kg/m3]
		*/
		double density_r_max() const;

		/**
		*@brief Returns the minimum density (reynolds) [kg/m3]
		*/
		double density_r_min() const;

		/**
		*@brief Returns true is the solution is cold (i.e. no flame)
		*/
		bool   cold() const;

		/**
		*@brief Returns the fuel density (reynolds) [kg/m3]
		*/
		double density_r_fuel() const;

		/**
		*@brief Returns the oxidizer density (reynolds) [kg/m3]
		*/
		double density_r_oxidizer() const;

		/**
		*@brief Returns the fuel enthalpy (favre) [J/kg]
		*/
		double enthalpy_f_fuel() const;

		/**
		*@brief Returns the oxidizer enthalpy (favre) [J/kg]
		*/
		double enthalpy_f_oxidizer() const;

		/**
		*@brief Returns the fuel temperature (favre) [K]
		*/
		double temperature_f_fuel() const;

		/**
		*@brief Returns the oxidizer temperature (favre) [K]
		*/
		double temperature_f_oxidizer() const;

		/**
		*@brief Returns the stoichiometric scalar dissipation rate [Hz]
		*/
		double stoichiometric_scalar_dissipation_rate() const;

		/**
		*@brief Returns the enthalpy defect (positive) [J/kg]
		*/
		double enthalpy_defect() const;

	private:

		std::string _name;			//!< XML file containing the flamelet solution

		double _chi_st;				//!< stoichiometric scalar dissipation rate [Hz]
		double _chi_max;			//!< maximum scalar dissipation rate [Hz]
		double _as;					//!< strain rate [Hz]
		double _enthalpy_defect;	//!< enthalpy defect [J/kg]

		unsigned int _n_csi;		//!< number of mixture fraction points
		unsigned int _n_variance;	//!< number of mixture fraction variance points;

		std::vector<double>	_csi;						//!< mixture fraction points (favre)
		std::vector<double>	_variance_normal;			//!< normal variance points (favre)

		std::vector< std::vector<double> > _mf_r;		//!< mixture fraction points (reynolds)
		std::vector< std::vector<double> > _mfv_r;		//!< normal variance points (reynolds)

		std::vector< std::vector<double> > _density_r;		//!< density (reynolds) [kg/m3]
		std::vector< std::vector<double> > _temperature_f;	//!< temperature (favre) [K]

		std::vector< std::vector<double> > _mw_f;			//!< molecular weight (favre) [kg/kmol]
		std::vector< std::vector<double> > _cp_f;			//!< constant pressure specific heat (favre) [J/kg]
		std::vector< std::vector<double> > _lambda_f;		//!< thermal conductivity (favre) [J/s/m/K]
		std::vector< std::vector<double> > _mu_f;			//!< dynamic viscosity (favre) [kg/m/s]
		std::vector< std::vector<double> > _as_f;			//!< absorption coefficient (favre) [1/m]
		std::vector< std::vector<double> > _alpha_f;		//!< thermal diffusivity (favre) [m2/s]

		std::vector< std::vector< std::vector<double> > > _w_f;		//!< Favre mass fractions (favre)

		unsigned int				_number_of_species;					//!< total number of species available
		std::vector<std::string>	_names_of_species_to_extract;		//!< list of species to be extracted

		double _temperature_f_max;	//!< maximum temperature (favre) [K]
		double _temperature_f_min;	//!< minimum temperature (favre) [K]
		double _density_r_max;		//!< maximum density (reynolds) [kg/m3]
		double _density_r_min;		//!< minimum density (reynolds) [kg/m3]
		bool   _cold;			//!< true if the solution is cold (i.e. no flame)
		double enthalpy_f_fuel_;	//!< enthalpy of fuel stream [J/kg]
		double enthalpy_f_oxidizer_;	//!< enthalpy of oxidizer stream [J/kg]


	private:

		void ErrorMessage(const std::string error_message);
		void WarningMessage(const std::string warning_message);
	};
}

#include "FlameletReader.hpp"

#endif // #ifndef OpenSMOKE_FlameletReader_H

	
