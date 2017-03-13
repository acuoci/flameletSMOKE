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
 
#include <boost/math/constants/constants.hpp>

namespace OpenSMOKE
{
	FlameletReader::FlameletReader()
	{
		_name = "NULL";
		_temperature_f_max = -1.e16;
		_temperature_f_min = 1.e16;
		_density_r_max = -1.e16;
		_density_r_min = 1.e16;
		_cold = false;
	}

	void FlameletReader::SetSpeciesToExtract(const std::vector<std::string> names)
	{
		_names_of_species_to_extract = names;
	}

	void FlameletReader::Read(const boost::filesystem::path file_name)
	{
		std::string tag;

		_name = file_name.filename().string();

		rapidxml::xml_document<> xml_main_input;
		std::vector<char> xml_main_input_string;

		if (!boost::filesystem::exists(file_name))
			OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " file does not exist");

		OpenSMOKE::OpenInputFileXML(xml_main_input, xml_main_input_string, file_name);

		// Stoichiometric scalar dissipation rate
		{
			rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("stoichiometric-scalar-dissipation-rate");
			if (indices_node != 0)
			{
				std::stringstream values(indices_node->value());
				values >> _chi_st;
				_as = _chi_st*boost::math::constants::pi<double>();
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the stoichiometric-scalar-dissipation-rate leaf");
		}

		// Maximum scalar dissipation rate
		{
			rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("max-scalar-dissipation-rate");
			if (indices_node != 0)
			{
				std::stringstream values(indices_node->value());
				values >> _chi_max;
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the max-scalar-dissipation-rate leaf");
		}

		// Enthalpy defect
		{
			rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("enthalpy-defect");
			if (indices_node != 0)
			{
				std::stringstream values(indices_node->value());
				values >> _enthalpy_defect;
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the enthalpy-defect leaf");
		}

		// Enthalpy of fuel stream
		{
			rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("enthalpy-fuel");
			if (indices_node != 0)
			{
				std::stringstream values(indices_node->value());
				values >> enthalpy_f_fuel_;
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the enthalpy-fuel leaf");
		}

		// Enthalpy of oxidizer stream
		{
			rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("enthalpy-oxidizer");
			if (indices_node != 0)
			{
				std::stringstream values(indices_node->value());
				values >> enthalpy_f_oxidizer_;
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the enthalpy-oxidizer leaf");
		}

		// Number of mixture fraction points
		{
			rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("number-mixture-fraction-points");
			if (indices_node != 0)
			{
				std::stringstream values(indices_node->value());
				values >> _n_csi;
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the number-mixture-fraction-points leaf");
		}

		// Number of variance space points
		{
			rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("number-variances");
			if (indices_node != 0)
			{
				std::stringstream values(indices_node->value());
				values >> _n_variance;
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the number-variances leaf");
		}

		// Number of species available
		{
			rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("number-species");
			if (indices_node != 0)
			{
				std::stringstream values(indices_node->value());
				values >> _number_of_species;
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the number-species leaf");
		}

		// Memory allocation
		{
			_csi.resize(_n_csi + 1);
			_variance_normal.resize(_n_variance + 1);

			_mf_r.resize(_n_variance + 1);
			_mfv_r.resize(_n_variance + 1);
			_temperature_f.resize(_n_variance + 1);
			_density_r.resize(_n_variance + 1);
			_mw_f.resize(_n_variance + 1);
			_cp_f.resize(_n_variance + 1);
			_lambda_f.resize(_n_variance + 1);
			_mu_f.resize(_n_variance + 1);
			_as_f.resize(_n_variance + 1);
			_alpha_f.resize(_n_variance + 1);

			for (unsigned int j = 1; j <= _n_variance; j++)
			{
				_mf_r[j].resize(_n_csi + 1);
				_mfv_r[j].resize(_n_csi + 1);
				_temperature_f[j].resize(_n_csi + 1);
				_density_r[j].resize(_n_csi + 1);
				_mw_f[j].resize(_n_csi + 1);
				_cp_f[j].resize(_n_csi + 1);
				_lambda_f[j].resize(_n_csi + 1);
				_mu_f[j].resize(_n_csi + 1);
				_as_f[j].resize(_n_csi + 1);
				_alpha_f[j].resize(_n_csi + 1);
			}
		}

		// Read mixture fraction
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("mixture-fraction");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for(unsigned int i=0;i<_n_csi;i++)
					values >> _csi[_n_csi-i];
				for(unsigned int i=0;i<_n_csi;i++)
					_csi[_n_csi-i] = 1. - _csi[_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the mixture-fraction leaf");
		}

		// Read normalized variances
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("variances");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i = 0; i<_n_variance; i++)
					values >> _variance_normal[i + 1];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the variances leaf");
		}

		// Read mixture fraction (Reynolds)
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("mixture-fraction-reynolds");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i=0;i<_n_csi;i++)
					for (unsigned int j=0;j<_n_variance;j++)
						values >> _mf_r[j+1][_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the mixture-fraction-reynolds leaf");
		}

		// Read normalized variances (Reynolds)
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("variances-reynolds");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i = 0; i<_n_csi; i++)
					for (unsigned int j = 0; j<_n_variance; j++)
						values >> _mfv_r[j + 1][_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the variances-reynolds leaf");
		}

		// Read density (Reynolds)
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("density-reynolds");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i = 0; i<_n_csi; i++)
					for (unsigned int j = 0; j<_n_variance; j++)
						values >> _density_r[j + 1][_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the density-reynolds leaf");
		}
		
		// Read temperature (Favre)
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("temperature-favre");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i = 0; i<_n_csi; i++)
					for (unsigned int j = 0; j<_n_variance; j++)
						values >> _temperature_f[j + 1][_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the temperature-favre leaf");
		}

		// Read molecular weight (Favre)
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("mw-favre");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i = 0; i<_n_csi; i++)
					for (unsigned int j = 0; j<_n_variance; j++)
						values >> _mw_f[j + 1][_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the mw-favre leaf");
		}

		// Read constant pressure specific heat (Favre)
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("cp-favre");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i = 0; i<_n_csi; i++)
					for (unsigned int j = 0; j<_n_variance; j++)
						values >> _cp_f[j + 1][_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the cp-favre leaf");
		}

		// Read thermal conductivity (Favre)
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("lambda-favre");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i = 0; i<_n_csi; i++)
					for (unsigned int j = 0; j<_n_variance; j++)
						values >> _lambda_f[j + 1][_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the lambda-favre leaf");
		}

		// Read dynamic viscosity (Favre)
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("mu-favre");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i = 0; i<_n_csi; i++)
					for (unsigned int j = 0; j<_n_variance; j++)
						values >> _mu_f[j + 1][_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the mu-favre leaf");
		}

		// Read Planck mean absorption coefficient (Favre)
		{
			rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("k-planck-favre");

			if (profiles_size_node != 0)
			{
				std::stringstream values(profiles_size_node->value());
				for (unsigned int i = 0; i<_n_csi; i++)
					for (unsigned int j = 0; j<_n_variance; j++)
						values >> _as_f[j + 1][_n_csi-i];
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the k-planck-favre leaf");
		}

		// Reconstruct thermal diffusivity [m2/s]
		{
			for (unsigned int j = 0; j < _n_variance; j++)
				for (unsigned int i = 0; i < _n_csi; i++)
					_alpha_f[j + 1][i + 1] = _lambda_f[j + 1][i + 1] / _density_r[j + 1][i + 1] / _cp_f[j + 1][i + 1];
		}

		// Search for output species (if any)
		if (_names_of_species_to_extract.size() > 0)
		{
			std::vector<std::string> species_available(_number_of_species);

			// Read available species
			{
				rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node("species-names");

				if (profiles_size_node != 0)
				{
					std::stringstream values(profiles_size_node->value());
					for (unsigned int i = 0; i<_number_of_species; i++)
							values >> species_available[i];
				}
				else
					OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the species-names leaf");
			}

			// Allocate memory
			_w_f.resize(_names_of_species_to_extract.size() + 1);
			for (unsigned int k = 1; k <= _names_of_species_to_extract.size(); k++)
			{
				_w_f[k].resize(_n_variance + 1);
				for (unsigned int j = 1; j <= _n_variance; j++)
					_w_f[k][j].resize(_n_csi + 1);
			}

			// Read mass fractions of species (if available)
			for (unsigned int k = 0; k < _names_of_species_to_extract.size(); k++)
			{
				std::string name = "mass-fraction-" + _names_of_species_to_extract[k];

				// Read mass fractions
				{
					rapidxml::xml_node<>* profiles_size_node = xml_main_input.first_node("opensmoke")->first_node(name.c_str());

					if (profiles_size_node != 0)
					{
						std::stringstream values(profiles_size_node->value());
						for (unsigned int i = 0; i<_n_csi; i++)
							for (unsigned int j = 0; j<_n_variance; j++)
								values >> _w_f[k+1][j+1][_n_csi-i];
					}
					else
						OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the " + name + " leaf");
				}
			}
		}

		// Search for maxima and minima of relevant variables
		for (unsigned int j = 0; j < _n_variance; j++)
			for (unsigned int i = 0; i < _n_csi; i++)
			{
				if (_temperature_f_max < _temperature_f[j + 1][i + 1])
					_temperature_f_max = _temperature_f[j + 1][i + 1];

				if (_temperature_f_min > _temperature_f[j + 1][i + 1])
					_temperature_f_min = _temperature_f[j + 1][i + 1];

				if (_density_r_max < _density_r[j + 1][i + 1])
					_density_r_max = _density_r[j + 1][i + 1];

				if (_density_r_min > _density_r[j + 1][i + 1])
					_density_r_min = _density_r[j + 1][i + 1];
			}

		if (_temperature_f_max <= 1.05*std::max(_temperature_f[1][1], _temperature_f[1][_n_csi]))
			_cold = true;
	}

	void FlameletReader::GetMeanValues(const double csi, const double csiv2, std::vector<double>& extracted)
	{
		int iCsi = 0;
		int jCsiv2 = 0;

		for (unsigned int j = 2; j <= _n_variance; j++)
			if (csiv2 <= _variance_normal[j])
			{
				jCsiv2 = j - 1;
				break;
			}

		for (unsigned int i = 2; i <= _n_csi; i++)
			if (csi <= _csi[i])
			{
				iCsi = i - 1;
				break;
			}

		const double q = (_csi[iCsi + 1] - _csi[iCsi])*(_variance_normal[jCsiv2 + 1] - _variance_normal[jCsiv2]);
		const double q11 = (_csi[iCsi + 1] - csi) * (_variance_normal[jCsiv2 + 1] - csiv2);
		const double q21 = (csi - _csi[iCsi])   * (_variance_normal[jCsiv2 + 1] - csiv2);
		const double q12 = (_csi[iCsi + 1] - csi) * (csiv2 - _variance_normal[jCsiv2]);
		const double q22 = (csi - _csi[iCsi])   * (csiv2 - _variance_normal[jCsiv2]);

		// Temperature Favre [K]
		extracted[1] = (_temperature_f[jCsiv2][iCsi] * q11 + _temperature_f[jCsiv2][iCsi + 1] * q21 +
			_temperature_f[jCsiv2 + 1][iCsi] * q12 + _temperature_f[jCsiv2 + 1][iCsi + 1] * q22) / q;

		// Density Reynolds [kg/m3]
		extracted[2] = (_density_r[jCsiv2][iCsi] * q11 + _density_r[jCsiv2][iCsi + 1] * q21 +
			_density_r[jCsiv2 + 1][iCsi] * q12 + _density_r[jCsiv2 + 1][iCsi + 1] * q22) / q;

		// Absorption Coefficient Favre [1/m]
		extracted[3] = (_as_f[jCsiv2][iCsi] * q11 + _as_f[jCsiv2][iCsi + 1] * q21 +
			_as_f[jCsiv2 + 1][iCsi] * q12 + _as_f[jCsiv2 + 1][iCsi + 1] * q22) / q;

		// Dynamic Viscosity [kg/m/s]
		extracted[4] = (_mu_f[jCsiv2][iCsi] * q11 + _mu_f[jCsiv2][iCsi + 1] * q21 +
			_mu_f[jCsiv2 + 1][iCsi] * q12 + _mu_f[jCsiv2 + 1][iCsi + 1] * q22) / q;

		// Thermal diffusivity [m2/s]
		extracted[5] = (_alpha_f[jCsiv2][iCsi] * q11 + _alpha_f[jCsiv2][iCsi + 1] * q21 +
			_alpha_f[jCsiv2 + 1][iCsi] * q12 + _alpha_f[jCsiv2 + 1][iCsi + 1] * q22) / q;

		// Specific heat [J/kg/K]
		extracted[6] = (_cp_f[jCsiv2][iCsi] * q11 + _cp_f[jCsiv2][iCsi + 1] * q21 +
			_cp_f[jCsiv2 + 1][iCsi] * q12 + _cp_f[jCsiv2 + 1][iCsi + 1] * q22) / q;
	}

	void FlameletReader::ExtractMeanValues(const double csi, const double csiv2, std::vector<double> &omegaFavre)
	{
		int iCsi = 0;
		int jCsiv2 = 0;

		for (unsigned int j = 2; j <= _n_variance; j++)
			if (csiv2 <= _variance_normal[j])
			{
				jCsiv2 = j - 1;
				break;
			}

		for (unsigned int i = 2; i <= _n_csi; i++)
			if (csi <= _csi[i])
			{
				iCsi = i - 1;
				break;
			}

		double q = (_csi[iCsi + 1] - _csi[iCsi])*(_variance_normal[jCsiv2 + 1] - _variance_normal[jCsiv2]);
		double q11 = (_csi[iCsi + 1] - csi) * (_variance_normal[jCsiv2 + 1] - csiv2);
		double q21 = (csi - _csi[iCsi])   * (_variance_normal[jCsiv2 + 1] - csiv2);
		double q12 = (_csi[iCsi + 1] - csi) * (csiv2 - _variance_normal[jCsiv2]);
		double q22 = (csi - _csi[iCsi])   * (csiv2 - _variance_normal[jCsiv2]);

		for (unsigned int w = 1; w <= _names_of_species_to_extract.size(); w++)
			omegaFavre[w] = (_w_f[w][jCsiv2][iCsi] * q11 + _w_f[w][jCsiv2][iCsi + 1] * q21 +
				_w_f[w][jCsiv2 + 1][iCsi] * q12 + _w_f[w][jCsiv2 + 1][iCsi + 1] * q22) / q;
	}

	void FlameletReader::Summary()
	{

		std::cout << std::endl;
		std::cout << "Flamelet properties                    " << std::endl;
		std::cout << "  + File name:                         " << _name << std::endl;
		std::cout << "  + Mixture fraction points:           " << _n_csi << std::endl;
		std::cout << "  + Variance Mixture fraction points:  " << _n_variance << std::endl;
		std::cout << "  + Enthalpy defect:                   " << _enthalpy_defect/1.e3 << " kJ/kg" << std::endl;
		std::cout << "  + Strain rate:                       " << _as << " Hz" << std::endl;
		std::cout << "  + Stoic. scalar dissipation rate:    " << _chi_st << " Hz" << std::endl;
		std::cout << "  + Max scalar dissipation rate:       " << _chi_max << " Hz" << std::endl;
		std::cout << "  + Temperature min:                   " << _temperature_f_min << " K" << std::endl;
		std::cout << "  + Temperature max:                   " << _temperature_f_max << " K" << std::endl;
		std::cout << "  + Density min:                       " << _density_r_min << " kg/m3" << std::endl;
		std::cout << "  + Density max:                       " << _density_r_max << " kg/m3" << std::endl;
		std::cout << "  + Cold flame:                        " << _cold << std::endl;
		std::cout << std::endl;

	}

	void FlameletReader::ErrorMessage(const std::string error_message)
	{
		std::cout << "FlameletReader Error" << std::endl;
		std::cout << "File name:     " << _name << std::endl;
		std::cout << "Error message: " << error_message << std::endl;
		exit(-1);
	}

	void FlameletReader::WarningMessage(const std::string warning_message)
	{
		std::cout << "FlameletReader Warning" << std::endl;
		std::cout << "File name:     " << _name << std::endl;
		std::cout << "Warning message: " << warning_message << std::endl;
		std::cout << "Press enter to continue..." << std::endl;
		getchar();
	}

	double FlameletReader::temperature_f_max() const
	{
		return _temperature_f_max;
	}

	double FlameletReader::temperature_f_min() const
	{
		return _temperature_f_min;
	}

	double FlameletReader::density_r_max() const
	{
		return _density_r_max;
	}

	double FlameletReader::density_r_min() const
	{
		return _density_r_min;
	}

	bool FlameletReader::cold() const
	{
		return _cold;
	}

	double FlameletReader::density_r_fuel() const
	{
		return _density_r[1][_n_csi];
	}

	double FlameletReader::density_r_oxidizer() const
	{
		return _density_r[1][1];
	}

	double FlameletReader::enthalpy_f_fuel() const
	{
		return enthalpy_f_fuel_;
	}

	double FlameletReader::enthalpy_f_oxidizer() const
	{
		return enthalpy_f_oxidizer_;
	}

	double FlameletReader::temperature_f_fuel() const
	{
		return _temperature_f[1][_n_csi];
	}

	double FlameletReader::temperature_f_oxidizer() const
	{
		return  _temperature_f[1][1];
	}

	double FlameletReader::stoichiometric_scalar_dissipation_rate() const
	{
		return _chi_st;
	}
	
	double FlameletReader::enthalpy_defect() const
	{
		return _enthalpy_defect;
	}
}

	
