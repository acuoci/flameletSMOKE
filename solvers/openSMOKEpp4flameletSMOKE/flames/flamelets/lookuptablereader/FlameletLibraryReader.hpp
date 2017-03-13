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
	FlameletLibraryReader::FlameletLibraryReader()
	{
		_temperature_f_max = -1.e16;
		_temperature_f_min = 1.e16;
		_density_r_max = -1.e16;
		_density_r_min = 1.e16;
		_iMultiScalarDissipationRates = true;
		_showFlamelet = false;

		X_pdf_type_ = X_PDF_DIRAC_DELTA;
		_chi_log_normal_sigma = 1.31;
	}

	void FlameletLibraryReader::SetShowFlamelet()
	{
		_showFlamelet = true;
	}

	void FlameletLibraryReader::SetLogNormalChiDistribution(const double sigma)
	{
		X_pdf_type_ = X_PDF_LOG_NORMAL;
		_chi_log_normal_sigma = sigma;
	}

	void FlameletLibraryReader::SetSpeciesToExtract(const std::vector<std::string> names)
	{
		_names_of_species_to_extract = names;
	}

	void FlameletLibraryReader::Read(const std::vector<boost::filesystem::path> list_of_flamelets)
	{
		// Number of scalar dissipation rates
		_n = list_of_flamelets.size();

		// Read the flamelets
		_flame.resize(_n + 1);
		for (unsigned int i = 1; i <= _n; i++)
		{
			if (_names_of_species_to_extract.size() > 0)
				_flame[i].SetSpeciesToExtract(_names_of_species_to_extract);
			_flame[i].Read(list_of_flamelets[i-1]);
			if (_showFlamelet == true)
				_flame[i].Summary();
		}

		// Assign the stoichiometric scalar dissipation rates
		_chi_st.resize(_n);
		for (unsigned int i = 0; i < _n; i++)
			_chi_st[i] = _flame[i+1].stoichiometric_scalar_dissipation_rate();

		// Assign the enthalpy defect (must be the same for every flamelet)
		_enthalpy_defect.resize(_n + 1);
		for (unsigned int i = 1; i <= _n; i++)
			_enthalpy_defect[i] = _flame[i].enthalpy_defect();

		// Search for maxima and minima of relevant variables
		for (unsigned int i = 1; i <= _n; i++)
		{
			if (_temperature_f_max < _flame[i].temperature_f_max())
				_temperature_f_max = _flame[i].temperature_f_max();

			if (_temperature_f_min > _flame[i].temperature_f_min())
				_temperature_f_min = _flame[i].temperature_f_min();

			if (_density_r_max < _flame[i].density_r_max())
				_density_r_max = _flame[i].density_r_max();

			if (_density_r_min > _flame[i].density_r_min())
				_density_r_min = _flame[i].density_r_min();
		}

		// Check if the scalar dissipation rates have been provided in the right order
		for (unsigned int i = 0; i < _n - 1; i++)
			if (_chi_st[i] >= _chi_st[i + 1])
				ErrorMessage("The scalar dissipation rates are not in the right order!");

		// Check if the enthalpy defect is the same for all the flamelets
		for (unsigned int i = 1; i <= _n - 1; i++)
			if (_enthalpy_defect[i] != _enthalpy_defect[i + 1])
				ErrorMessage("The enthalpy defect is not the same for every flamelet!");

		// Minimum and maximum scalar dissipation rates available
		_chi_st_min = _chi_st[0];
		_chi_st_max = _chi_st[_n-1];

		// Check if only a single scalar dissipation rate is available
		if (_n == 1)
			_iMultiScalarDissipationRates = false;

		// Additional check on the number of scalar dissipation rates
		if (_n < 3 && X_pdf_type_ == X_PDF_LOG_NORMAL)
			ErrorMessage("Log-normal distribution function for scalar dissipation rate can be used only with at least three flamelets...");

		// Prepare the log-normal distribution function for the scalar dissipation rate
		if (X_pdf_type_ == X_PDF_LOG_NORMAL)
		{
			logNormalPDF_.SetSigma(_chi_log_normal_sigma);
			logNormalPDF_.SetX(_chi_st);

			_t_favre.resize(_n);
			_rho_reynolds.resize(_n);
			_cp_favre.resize(_n);
			_mu_favre.resize(_n);
			_as_favre.resize(_n);
			_alpha_favre.resize(_n);

			if (_names_of_species_to_extract.size() > 0)
			{
				_w_favre.resize(_n);
				for (unsigned int j = 0; j < _n; j++)
					_w_favre[j].resize(_names_of_species_to_extract.size() + 1);
			}
		}

		_temperature_f_fuel = _flame[1].temperature_f_fuel();
		_temperature_f_oxidizer = _flame[1].temperature_f_oxidizer();
		_enthalpy_f_fuel = _flame[1].enthalpy_f_fuel();
		_enthalpy_f_oxidizer = _flame[1].enthalpy_f_oxidizer();
		_density_r_fuel = _flame[1].density_r_fuel();
		_density_r_oxidizer = _flame[1].density_r_oxidizer();
	}

	void FlameletLibraryReader::GetMeanValues(const double csi, const double csiv2, const double chi_st, std::vector<double>& extracted)
	{
		if (_iMultiScalarDissipationRates == true)
		{
			if (X_pdf_type_ == X_PDF_LOG_NORMAL)
			{
				GetMeanValuesMultiScalarDissipationRatesLogNormal(csi, csiv2, chi_st, extracted);
			}
			else
			{
				GetMeanValuesMultiScalarDissipationRatesDirac(csi, csiv2, chi_st, extracted);
			}
		}
		else
		{
			_flame[1].GetMeanValues(csi, csiv2, extracted);
		}
	}

	void FlameletLibraryReader::GetMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, std::vector<double>& extracted)
	{
		if (chi_st <= _chi_st_min)
		{
			_flame[1].GetMeanValues(csi, csiv2, extracted);
			return;
		}
		else if (chi_st >= _chi_st_max)
		{
			_flame[_n].GetMeanValues(csi, csiv2, extracted);
			return;
		}
		else
		{
			int kChi = 0;
			double ratio;
			std::vector<double> extracted_1(extracted.size());
			std::vector<double> extracted_2(extracted.size());

			// Finding kChi
			for (unsigned int k = 2; k <= _n; k++)
				if (chi_st <= _chi_st[k-1])
				{
					kChi = k - 1;
					break;
				}

			// Finding flamelets mean values
			_flame[kChi].GetMeanValues(csi, csiv2, extracted_1);
			_flame[kChi + 1].GetMeanValues(csi, csiv2, extracted_2);

			// Interpolation ratio
			ratio = (chi_st - _chi_st[kChi-1]) / (_chi_st[kChi] - _chi_st[kChi-1]);

			// Mean values
			for (unsigned int k = 1; k <= extracted.size() - 1; k++)
				extracted[k] = extracted_1[k] + ratio*(extracted_2[k] - extracted_1[k]);
		}
	}

	void FlameletLibraryReader::GetMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, std::vector<double>& extracted)
	{
		logNormalPDF_.SetXmean(chi_st);

		for (unsigned int j=0;j<_n;j++)
		{
			_flame[j+1].GetMeanValues(csi, csiv2, extracted);

			_t_favre[j] = extracted[1];
			_rho_reynolds[j] = extracted[2];
			_as_favre[j] = extracted[3];
			_mu_favre[j] = extracted[4];
			_alpha_favre[j] = extracted[5];
			_cp_favre[j] = extracted[6];
		}

		extracted[1] = logNormalPDF_.GetMeanValue(_t_favre);
		extracted[2] = logNormalPDF_.GetMeanValue(_rho_reynolds);
		extracted[3] = logNormalPDF_.GetMeanValue(_as_favre);
		extracted[4] = logNormalPDF_.GetMeanValue(_mu_favre);
		extracted[5] = logNormalPDF_.GetMeanValue(_alpha_favre);
		extracted[6] = logNormalPDF_.GetMeanValue(_cp_favre);
	}


	void FlameletLibraryReader::ExtractMeanValues(const double csi, const double csiv2, const double chi_st, std::vector<double> &omegaFavre)
	{
		if (_iMultiScalarDissipationRates == true)
		{
			if (X_pdf_type_ == X_PDF_LOG_NORMAL)
				ExtractMeanValuesMultiScalarDissipationRatesLogNormal(csi, csiv2, chi_st, omegaFavre);
			else
				ExtractMeanValuesMultiScalarDissipationRatesDirac(csi, csiv2, chi_st, omegaFavre);
		}
		else
			_flame[1].ExtractMeanValues(csi, csiv2, omegaFavre);
	}

	void FlameletLibraryReader::ExtractMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, std::vector<double> &omegaFavre)
	{
		if (chi_st <= _chi_st_min)
		{
			_flame[1].ExtractMeanValues(csi, csiv2, omegaFavre);
			return;
		}
		else if (chi_st >= _chi_st_max)
		{
			_flame[_n].ExtractMeanValues(csi, csiv2, omegaFavre);
			return;
		}
		else
		{
			int kChi = 0;
			double ratio;
			std::vector<double> omegaFavre_1, omegaFavre_2;
			omegaFavre_1.resize(_names_of_species_to_extract.size() + 1);
			omegaFavre_2.resize(_names_of_species_to_extract.size() + 1);

			// Finding kChi
			for (unsigned int k = 2; k <= _n; k++)
				if (chi_st <= _chi_st[k-1])
				{
					kChi = k - 1;
					break;
				}

			// Finding flamelets mean values
			_flame[kChi].ExtractMeanValues(csi, csiv2, omegaFavre_1);
			_flame[kChi + 1].ExtractMeanValues(csi, csiv2, omegaFavre_2);

			// Interpolation ratio
			ratio = (chi_st - _chi_st[kChi-1]) / (_chi_st[kChi] - _chi_st[kChi-1]);

			// Mean values
			for (unsigned int j = 1; j <= _names_of_species_to_extract.size(); j++)
				omegaFavre[j] = omegaFavre_1[j] + ratio*(omegaFavre_2[j] - omegaFavre_1[j]);
		}
	}

	void FlameletLibraryReader::ExtractMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, std::vector<double> &omegaFavre)
	{
		// cuoci _chi_pdf.AssignMeanScalarDissipationRate(chi_st);
		logNormalPDF_.SetXmean(chi_st);

		for (unsigned int j=0;j<_n;j++)
			_flame[j+1].ExtractMeanValues(csi, csiv2, _w_favre[j]);

		//cuoci omegaFavre = _chi_pdf.ExtractMeanValue(_w_favre);
		omegaFavre = logNormalPDF_.GetMeanValue(_w_favre);
	}

	void FlameletLibraryReader::Summary()
	{
		std::cout << std::endl;
		std::cout << "Flamelet library properties            " << std::endl;
		std::cout << "  * Number of Flamelets:               " << _n << std::endl;
		std::cout << "  * Enthalpy defect:                   " << _enthalpy_defect[1]/1000.	<< " J/kg" << std::endl;
		std::cout << "  * Min stoic. scalar diss. rate:      " << _chi_st_min 			<< " 1/s" << std::endl;
		std::cout << "  * Max stoic. scalar diss. rate:      " << _chi_st_max 			<< " 1/s" << std::endl;
		std::cout << "  * Temperature min:                   " << _temperature_f_min 		<< " K" << std::endl;
		std::cout << "  * Temperature max:                   " << _temperature_f_max 		<< " K" << std::endl;
		std::cout << "  * Density min:                       " << _density_r_min 		<< " kg/m3" << std::endl;
		std::cout << "  * Density max:                       " << _density_r_max 		<< " kg/m3" << std::endl;
		std::cout << "  * Temperature fuel:                  " << _temperature_f_fuel 		<< " K" << std::endl;
		std::cout << "  * Temperature oxidizer:              " << _temperature_f_oxidizer 	<< " K" << std::endl;
		std::cout << "  * Density fuel:                      " << _density_r_fuel 		<< " kg/m3" << std::endl;
		std::cout << "  * Density oxidizer:                  " << _density_r_oxidizer 		<< " kg/m3" << std::endl;
		std::cout << "  * Enthalpy fuel:                     " << _enthalpy_f_fuel 		<< " J/kg" << std::endl;
		std::cout << "  * Enthalpy oxidizer:                 " << _enthalpy_f_oxidizer 		<< " J/kg" << std::endl << std::endl;
		std::cout << std::endl;
	}

	void FlameletLibraryReader::ErrorMessage(const std::string error_message)
	{
		std::cout << "FlameletLibraryReader" << std::endl;
		std::cout << "Error message: " << error_message << std::endl;
		exit(-1);
	}

	void FlameletLibraryReader::WarningMessage(const std::string warning_message)
	{
		std::cout << "FlameletLibraryReader" << std::endl;
		std::cout << "Warning message: " << warning_message << std::endl;
		std::cout << "Press enter to continue..." << std::endl;
		getchar();
	}

	double FlameletLibraryReader::density_r_fuel() const
	{
		return _density_r_fuel;
	}

	double FlameletLibraryReader::density_r_oxidizer() const
	{
		return _density_r_oxidizer;
	}

	double FlameletLibraryReader::enthalpy_f_fuel() const
	{
		return _enthalpy_f_fuel;
	}

	double FlameletLibraryReader::enthalpy_f_oxidizer() const
	{
		return _enthalpy_f_oxidizer;
	}

	double FlameletLibraryReader::temperature_f_fuel() const
	{
		return _temperature_f_fuel;
	}

	double FlameletLibraryReader::temperature_f_oxidizer() const
	{
		return _temperature_f_oxidizer;
	}
}


	
