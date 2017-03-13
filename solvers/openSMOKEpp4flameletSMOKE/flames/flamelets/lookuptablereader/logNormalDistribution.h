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
 
#ifndef OpenSMOKE_logNormalDistribution_H
#define OpenSMOKE_logNormalDistribution_H

#include "OpenSMOKEpp"
#include "math/OpenSMOKEUtilities.h"
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>


namespace OpenSMOKE
{
	//!  A class to manage log-Normal distributions
	/*!
	This class provides the tools to manage log-Normal distributions
	*/

	class logNormalDistribution
	{

	public:

		/**
		*@brief Default constructor
		*/
		logNormalDistribution();

		/**
		*@brief Sets the scale factor for the log-normal distribution
		*@param sigma scale factor
		*/
		void SetSigma(const double sigma);

		/**
		*@brief Sets the minimum scalar dissipation rate which is allowed
		*@param X_min the minimum allowed scalar dissipation rate [1/s]
		*/
		void SetXMin(const double X_min);

		/**
		*@brief Sets the maximum scalar dissipation rate which is allowed
		*@param X_max the maximum allowed scalar dissipation rate [1/s]
		*/
		void SetXMax(const double X_max);

		/**
		*@brief Sets the list of available scalar dissipation rates (0-index based)
		*@param X list of available scalar dissipation rates [1/s]
		*/
		void SetX(const std::vector<double> X);

		/**
		*@brief Sets the mean stoichiometric scalar dissipation rate
		*@param Xmean the mean stoichiometric scalar dissipation rate [1/s]
		*/
		void SetXmean(const double Xmean);

		/**
		*@brief Returns the integrated mean value
		*@param values list of values corresponding to the assigned scalar dissipation rates (0-index based)
		*/
		double GetMeanValue(const std::vector<double> &values);

		/**
		*@brief Returns the integrated mean value
		*@param values list of values corresponding to the assigned scalar dissipation rates (0-index based)
		*/
		std::vector<double> GetMeanValue(const std::vector < std::vector<double> > &values);

	private:

		double sigma_;				//!< scale parameter for the distribution
		double sigma_squared_;			//!< squared scale parameter
		double X_min_;				//!< minimum allowed scalar dissipation rate [1/s]
		double X_max_;				//!< maximum allowed scalar dissipation rate [1/s]
		unsigned int n_;			//!< number of available scalar dissipation rates (inclding the additional one corresponding to X_max_)
		
		std::vector<double> X_;			//!< available scalar dissipation rates [1/s]
		std::vector<double> integral_0_;	//!< integrals corresponding to int(pdf,X(i),X(i+1))
		std::vector<double> integral_1_;	//!< integrals corresponding to int(X*pdf,X(i),X(i+1))

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

#include "logNormalDistribution.hpp"

#endif // #ifndef OpenSMOKE_logNormalDistribution_H

	
