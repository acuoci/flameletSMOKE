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
	NonAdiabaticFlameletLibraryReader::NonAdiabaticFlameletLibraryReader()
	{
		X_pdf_type_ = FlameletLibraryReader::X_PDF_DIRAC_DELTA;
		_adiabatic_mode = false;
		_showFlamelet = false;
		_showFlameletLibrary = false;

		folder_name_ = "LibraryXML";
		_chi_log_normal_sigma = 1.31;
	}

	void NonAdiabaticFlameletLibraryReader::SetLibraryFolder(const boost::filesystem::path folder_name)
	{
		folder_name_ = folder_name;
	}

	void NonAdiabaticFlameletLibraryReader::SetLogNormalChiDistribution(const double sigma)
	{
		X_pdf_type_ = FlameletLibraryReader::X_PDF_LOG_NORMAL;
		_chi_log_normal_sigma = sigma;
	}

	void NonAdiabaticFlameletLibraryReader::SetSpeciesToExtract(const std::vector<std::string> names)
	{
		names_of_species_to_extract_ = names;
	}

	void NonAdiabaticFlameletLibraryReader::SetSpeciesToExtract(const std::string list_of_names)
	{
		std::vector<std::string> names;
		std::stringstream stream(list_of_names);
		while (!stream.eof())
		{
			std::string dummy;
			stream >> dummy;
			if (dummy == "" || dummy == "\n")
				continue;
			names.push_back(dummy);
		}
		SetSpeciesToExtract(names);
	}

	void NonAdiabaticFlameletLibraryReader::SetAdiabaticMode()
	{
		_adiabatic_mode = true;
	}

	void NonAdiabaticFlameletLibraryReader::SetShowFlamelet()
	{
		_showFlamelet = true;
	}
	void NonAdiabaticFlameletLibraryReader::SetShowFlameletLibrary()
	{
		_showFlameletLibrary = true;
	}

	int NonAdiabaticFlameletLibraryReader::number_of_species()
	{
		return names_of_species_to_extract_.size();
	}

	const std::vector<std::string> NonAdiabaticFlameletLibraryReader::species()
	{
		return names_of_species_to_extract_;
	}

	int NonAdiabaticFlameletLibraryReader::index_of_species(const std::string name)
	{
		for (unsigned int j = 1; j <= names_of_species_to_extract_.size(); j++)
			if (names_of_species_to_extract_[j-1] == name)
				return j;
		ErrorMessage("Species " + name + " cannot be extracted from the flamelet library");
		return 0;
	}

	void NonAdiabaticFlameletLibraryReader::Read()
	{
		Read(folder_name_);
	}

	void NonAdiabaticFlameletLibraryReader::Read(const boost::filesystem::path folder_name)
	{
		folder_name_ = folder_name;

		std::vector<GroupOfFlamelets> groups;
		SearchForXMLFiles(folder_name, groups);

		_nphi = groups.size();

		_phi.resize(_nphi + 1);
		for (unsigned int i = 1; i <= _nphi; i++)
			_phi[i] = -groups[i - 1].enthalpy_defect;

		for (unsigned int i = 1; i <= _nphi - 1; i++)
			if (_phi[i] <= _phi[i + 1])
				ErrorMessage("Enthalpy defects must be specified in decreasing order!");

		_jPhiAdiabatic = 0;
		for (unsigned int i = 1; i <= _nphi; i++)
			if (_phi[i] == 0.)
			{
				_jPhiAdiabatic = i;
				break;
			}

		if (_jPhiAdiabatic == 0)
			ErrorMessage("Missing adiabatic enthalppy defect!");

		// Maximum enthalpy defect
		_phi_max = _phi[_nphi];

		// Minimum enthalpy defect
		_phi_min = _phi[1];

		// Settings
		if (_nphi == 1)	SetAdiabaticMode();

		// Memory allocation
		flamelet_libraries.resize(_nphi + 1);

		// Read flamelet sub-libraries
		for (unsigned int i = 1; i <= _nphi; i++)
		{
			std::cout << " * Reading sub-library #" << i << " (Enthalpy defect: " << _phi[i]/1000. << " kJ/kg)" << std::endl;

			if (X_pdf_type_ == FlameletLibraryReader::X_PDF_LOG_NORMAL)
				flamelet_libraries[i].SetLogNormalChiDistribution(_chi_log_normal_sigma);

			if (names_of_species_to_extract_.size() > 0)
				flamelet_libraries[i].SetSpeciesToExtract(names_of_species_to_extract_);

			if (_showFlamelet == true)
				flamelet_libraries[i].SetShowFlamelet();

			flamelet_libraries[i].Read(groups[i-1].files);

			if (_showFlameletLibrary == true)
				flamelet_libraries[i].Summary();

			if (_adiabatic_mode == true)
				break;
		}

		// Adiabatic details
		_temperature_f_fuel = flamelet_libraries[_jPhiAdiabatic].temperature_f_fuel();
		_temperature_f_oxidizer = flamelet_libraries[_jPhiAdiabatic].temperature_f_oxidizer();
		_enthalpy_f_fuel = flamelet_libraries[_jPhiAdiabatic].enthalpy_f_fuel();
		_enthalpy_f_oxidizer = flamelet_libraries[_jPhiAdiabatic].enthalpy_f_oxidizer();
		_density_r_fuel = flamelet_libraries[_jPhiAdiabatic].density_r_fuel();
		_density_r_oxidizer = flamelet_libraries[_jPhiAdiabatic].density_r_oxidizer();
	}

	void NonAdiabaticFlameletLibraryReader::GetMeanValues(const double csi, const double csiv2, const double chi_st, const double phi, std::vector<double>& extracted)
	{
		if (_adiabatic_mode == true)
		{
			flamelet_libraries[1].GetMeanValues(csi, csiv2, chi_st, extracted);
		}
		else
		{
			if (phi >= _phi_min)
			{
				flamelet_libraries[1].GetMeanValues(csi, csiv2, chi_st, extracted);
				return;
			}
			else if (phi <= _phi_max)
			{
				flamelet_libraries[_nphi].GetMeanValues(csi, csiv2, chi_st, extracted);
				return;
			}
			else
			{
				int k = 0;
				double ratio;
				std::vector<double> extracted_1(extracted.size());
				std::vector<double> extracted_2(extracted.size());

				// Finding k
				for (unsigned int j = 2; j <= _nphi; j++)
					if (phi >= _phi[j])
					{
						k = j;
						break;
					}

				// Finding flamelets mean values
				flamelet_libraries[k].GetMeanValues(csi, csiv2, chi_st, extracted_1);
				flamelet_libraries[k - 1].GetMeanValues(csi, csiv2, chi_st, extracted_2);

				// Interpolation ratio
				ratio = (phi - _phi[k]) / (_phi[k - 1] - _phi[k]);

				// Mean values
				for (unsigned int i = 1; i <= extracted.size() - 1; i++)
					extracted[i] = extracted_1[i] + ratio*(extracted_2[i] - extracted_1[i]);
			}
		}
	}

	double NonAdiabaticFlameletLibraryReader::GetEnthalpyDefectFromTemperature(const double csi, const double csiv2, const double chi_st, const double T)
	{
		const int iTemperature = 1;
		std::vector<double> extracted(7);

		// Case: T>T0
		GetMeanValues(csi, csiv2, chi_st, _phi_min, extracted);
		double T0 = extracted[iTemperature];
		if (T > T0) return _phi_min;

		// Case: T<Tmin
		GetMeanValues(csi, csiv2, chi_st, _phi_max, extracted);
		double Tmin = extracted[iTemperature];
		if (T < Tmin) return _phi_max;

		// General case
		{
			double phiA = 0.;
			double phiB = 0.;
			double phiC, TC;

			// Finding k
			for (unsigned int j = 2; j <= _nphi; j++)
			{
				GetMeanValues(csi, csiv2, chi_st, _phi[j], extracted);
				double Tj = extracted[iTemperature];
				if (Tj <= T)
				{
					phiA = _phi[j];
					phiB = _phi[j - 1];
					break;
				}
			}

			int kmax = 15;
			for (int k = 1; k <= kmax; k++)
			{
				phiC = 0.50*(phiA + phiB);
				GetMeanValues(csi, csiv2, chi_st, phiC, extracted);
				TC = extracted[iTemperature];

				if (T <= TC)
				{
					phiB = phiC;
				}
				else
				{
					phiA = phiC;
				}

				if (fabs(TC - T) / T < 1.e-5)
					break;
			}

			return phiC;
		}
	}

	void NonAdiabaticFlameletLibraryReader::ExtractMeanValues(const double csi, const double csiv2, const double chi_st, const double phi, std::vector<double> &omegaFavre)
	{
		if (_adiabatic_mode == true)
		{
			flamelet_libraries[1].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre);
		}
		else
		{
			if (phi >= _phi_min)
			{
				flamelet_libraries[1].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre);
				return;
			}
			else if (phi <= _phi_max)
			{
				flamelet_libraries[_nphi].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre);
				return;
			}
			else
			{
				int k = 0;
				double ratio;
				std::vector<double> omegaFavre_1(omegaFavre.size());
				std::vector<double> omegaFavre_2(omegaFavre.size());

				// Finding k
				for (unsigned int j = 2; j <= _nphi; j++)
					if (phi >= _phi[j])
					{
						k = j;
						break;
					}

				// Finding flamelets mean values
				flamelet_libraries[k].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre_1);
				flamelet_libraries[k - 1].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre_2);

				// Interpolation ratio
				ratio = (phi - _phi[k]) / (_phi[k - 1] - _phi[k]);

				// Mean values
				for (unsigned int w = 1; w <= names_of_species_to_extract_.size(); w++)
					omegaFavre[w] = omegaFavre_1[w] + ratio*(omegaFavre_2[w] - omegaFavre_1[w]);
			}
		}
	}


	void NonAdiabaticFlameletLibraryReader::Summary()
	{
		std::cout << "Non adiabatic flamelet library summary" << std::endl;

		std::cout << "  + Folder:                      " << folder_name_.string() << std::endl;
		std::cout << "  + Number of enthalpy defects:  " << _nphi << std::endl;
		std::cout << "  + Minimum enthalpy defect:     " << _phi_min / 1000. << " kJ/kg" << std::endl;
		std::cout << "  + Maximum enthalpy defect:     " << _phi_max / 1000. << " kJ/kg" << std::endl;
		std::cout << "  + Adiabatic flamelets:         " << _jPhiAdiabatic << std::endl;
		std::cout << "  + Temperature fuel:            " << _temperature_f_fuel << " K" << std::endl;
		std::cout << "  + Temperature oxidizer:        " << _temperature_f_oxidizer << " K" << std::endl;
		std::cout << "  + Density fuel:                " << _density_r_fuel << " kg/m3" << std::endl;
		std::cout << "  + Density oxidizer:            " << _density_r_oxidizer << " kg/m3" << std::endl;
		std::cout << "  + Enthalpy fuel:               " << std::setprecision(6) << _enthalpy_f_fuel << " J/kg" << std::endl;
		std::cout << "  + Enthalpy oxidizer:           " << std::setprecision(6) << _enthalpy_f_oxidizer << " J/kg" << std::endl;
		std::cout << std::endl;
	}

	void NonAdiabaticFlameletLibraryReader::ErrorMessage(const std::string error_message)
	{
		std::cout << "NonAdiabaticFlameletLibraryReader Error" << std::endl;
		std::cout << "File name:     " << folder_name_.string() << std::endl;
		std::cout << "Error message: " << error_message << std::endl;
		exit(-1);
	}

	void NonAdiabaticFlameletLibraryReader::WarningMessage(const std::string warning_message)
	{
		std::cout << "NonAdiabaticFlameletLibraryReader Warning" << std::endl;
		std::cout << "File name:       " << folder_name_.string() << std::endl;
		std::cout << "Warning message: " << warning_message << std::endl;
		std::cout << "Press enter to continue..." << std::endl;
		getchar();
	}

	double NonAdiabaticFlameletLibraryReader::density_r_fuel() const
	{
		return _density_r_fuel;
	}

	double NonAdiabaticFlameletLibraryReader::density_r_oxidizer() const
	{
		return _density_r_oxidizer;
	}

	double NonAdiabaticFlameletLibraryReader::enthalpy_f_fuel() const
	{
		return _enthalpy_f_fuel;
	}

	double NonAdiabaticFlameletLibraryReader::enthalpy_f_oxidizer() const
	{
		return _enthalpy_f_oxidizer;
	}

	double NonAdiabaticFlameletLibraryReader::temperature_f_fuel() const
	{
		return _temperature_f_fuel;
	}

	double NonAdiabaticFlameletLibraryReader::temperature_f_oxidizer() const
	{
		return _temperature_f_oxidizer;
	}
}


	
