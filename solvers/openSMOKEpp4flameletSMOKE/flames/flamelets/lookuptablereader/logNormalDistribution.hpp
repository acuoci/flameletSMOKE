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
 
namespace OpenSMOKE
{
	logNormalDistribution::logNormalDistribution()
	{
		X_min_  	= 1.e-16;
		X_max_  	= 1.e5;
		n_ 		= 0;
		sigma_ 		= 1.31;
		sigma_squared_ = sigma_*sigma_;
	}

	void logNormalDistribution::SetSigma(const double sigma)
	{
		if (sigma < 0.50 || sigma > 5)
			ErrorMessage("Sigma must be between 0.50 and 5 (default: 1.31)");
		sigma_ = sigma;
		sigma_squared_ = sigma_*sigma_;
	}

	void logNormalDistribution::SetXMin(const double X_min)
	{
		if (X_min <= 1.e-32 || X_min > 1.e-6)
			ErrorMessage("Lower limit of scalar dissipation grid must be between 1.e-32 and 1.e-6 (default: 1.e-16)");
		X_min_ = X_min;
	}

	void logNormalDistribution::SetXMax(const double X_max)
	{
		if (X_max < 1. || X_max > 1.e6)
			ErrorMessage("Lower limit of scalar dissipation grid must be between 1.e0 and 1.e6 (default: 1.e3)");
		X_max_ = X_max;
	}

	void logNormalDistribution::SetX(const std::vector<double> X)
	{
		if (X[0] != 0.)
			ErrorMessage("The first flamelet must be at thermodynamic equilibrium (i.e. X=0)");
		if (X[1] <= X_min_)
			ErrorMessage("The second flamelet scalar dissipation rate must be large thant the minimum allowed value");
		if (X[X.size()-1] >= X_max_)
			ErrorMessage("Please adjust Xmax in your input file, because the maximum X you provided is larger than Xmax");
		
		n_= X.size()+1;
		X_.resize(n_);
		for(unsigned int i=0;i<n_-1;i++)
			X_[i] = X[i];
		X_[n_-1] = X_max_;

		integral_0_.resize(n_-1);
		integral_1_.resize(n_-1);

	}

	void logNormalDistribution::SetXmean(const double Xmean)
	{
		const double XmeanCorrected = std::max(X_min_, Xmean);

		// Parameters of the distribution
		const double M = std::log(XmeanCorrected) - sigma_squared_/2.;
		const double b = 1./2./sigma_squared_;
		const double coeff = std::exp(1./4./b+M)/(2.*sigma_*std::sqrt(2.*b));

		// Create the distribution
 		boost::math::lognormal pdf(M,sigma_);

		// Calculation of integral #1
		for(unsigned int i=0;i<n_-1;i++)
			integral_0_[i] = boost::math::cdf(pdf, X_[i+1]) - boost::math::cdf(pdf, X_[i]);

		// Calculation of integral #2
		{
			double XMeanCurrent = 0.;
			for(unsigned int i=0;i<n_-1;i++)
			{
				const double XA = std::max(X_[i],X_min_);
				const double XB = std::max(X_[i+1], X_min_);
				const double IA = coeff * boost::math::erf( (2.*b*(std::log(XA)-M)-1.)/(2.*std::sqrt(b)));
				const double IB = coeff * boost::math::erf( (2.*b*(std::log(XB)-M)-1.)/(2.*std::sqrt(b)));

				integral_1_[i] = IB-IA;
				XMeanCurrent += integral_1_[i];
			}

			const double correction = XmeanCorrected/XMeanCurrent;

			for(unsigned int i=0;i<n_-1;i++)
				integral_1_[i] *= correction;

			if (Xmean>=500.) 
				std::cout << "Very high X: " << Xmean << " vs " << XMeanCurrent << std::endl;
		}
	}

	double logNormalDistribution::GetMeanValue(const std::vector<double> &values)
	{
		double sum = 0.;
		for(unsigned int i=0;i<n_-2;i++)
		{
			const double m = (values[i+1]-values[i])/(X_[i+1]-X_[i]);
			const double F = values[i]-m*X_[i];
			sum += F*integral_0_[i] + m*integral_1_[i];
		}

		// Last interval
		// Option 1: keep the last value
		{
			unsigned int i = n_-2;
			const double m = 0.;
			const double F = values[i]-m*X_[i];
			sum += F*integral_0_[i] + m*integral_1_[i];
		}

		return sum;
	}


	std::vector<double> logNormalDistribution::GetMeanValue(const std::vector < std::vector<double> > &values)
	{
		const unsigned s = values[0].size();
		std::vector<double> sum(s);

		for(unsigned j=1;j<s;j++)
		{
			sum[j] = 0.;
			for(unsigned int i=0;i<n_-2;i++)
			{
				const double m = (values[i+1][j]-values[i][j])/(X_[i+1]-X_[i]);
				const double F = values[i][j]-m*X_[i];
				sum[j] += F*integral_0_[i] + m*integral_1_[i];
			}

			// Last interval
			// Option 1: keep the last value
			{
				unsigned int i = n_-2;
				const double m = 0.;
				const double F = values[i][j]-m*X_[i];
				sum[j] += F*integral_0_[i] + m*integral_1_[i];
			}
		}

		return sum;
	}

	void logNormalDistribution::ErrorMessage(const std::string error_message)
	{
		std::cout << "ScalarDissipationRateDistribution Error" << std::endl;
		std::cout << "Error message: " << error_message << std::endl;
		exit(-1);
	}

	void logNormalDistribution::WarningMessage(const std::string warning_message)
	{
		std::cout << "ScalarDissipationRateDistribution Warning" << std::endl;
		std::cout << "Warning message: " << warning_message << std::endl;
		std::cout << "Press enter to continue..." << std::endl;
		getchar();
	}
}
