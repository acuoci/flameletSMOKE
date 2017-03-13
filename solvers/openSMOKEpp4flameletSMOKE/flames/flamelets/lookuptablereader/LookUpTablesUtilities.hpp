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
 
namespace OpenSMOKE
{
	void SearchForXMLFiles(const boost::filesystem::path folder_name, std::vector<GroupOfFlamelets>& groups)
	{
		std::vector<boost::filesystem::path> list_xml_files;
		{
			if (!boost::filesystem::exists(folder_name) || !boost::filesystem::is_directory(folder_name))
				OpenSMOKE::FatalErrorMessage("The specified " + folder_name.string() + " folder does not exist");

			boost::filesystem::recursive_directory_iterator it(folder_name);
			boost::filesystem::recursive_directory_iterator endit;

			while (it != endit)
			{
				if (boost::filesystem::is_regular_file(*it) && it->path().extension() == ".xml")
					list_xml_files.push_back(it->path());
				++it;
			}
		}

		std::vector<double> list_stoichiometric_scalar_dissipation_rates(list_xml_files.size());
		std::vector<double> list_max_scalar_dissipation_rates(list_xml_files.size());
		std::vector<double> list_enthalpy_defects(list_xml_files.size());
		std::vector<double> list_max_temperatures(list_xml_files.size());
		for (unsigned int i = 0; i < list_xml_files.size(); i++)
			ReadHeaderSingleFlameletFromXMLFile(list_xml_files[i], list_stoichiometric_scalar_dissipation_rates[i],
				list_max_scalar_dissipation_rates[i], list_enthalpy_defects[i], list_max_temperatures[i]);

		for (unsigned int i = 0; i < list_xml_files.size(); i++)
		{
			std::cout << std::left << std::setw(5) << i + 1;
			std::cout << std::left << std::setw(14) << std::scientific << list_stoichiometric_scalar_dissipation_rates[i];
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(2) << list_enthalpy_defects[i] / 1000.;
			std::cout << std::left << std::setw(14) << std::fixed << std::setprecision(2) << list_max_temperatures[i];
			std::cout << std::left << list_xml_files[i].filename();
			std::cout << std::endl;
		}

		// Search how many enthalpy defects are available
		std::vector<double> enthalpy_defects_groups;
		for (unsigned int i = 0; i < list_enthalpy_defects.size(); i++)
		{
			bool is_already_available = false;
			for (unsigned int j = 0; j < enthalpy_defects_groups.size(); j++)
				if (list_enthalpy_defects[i] == enthalpy_defects_groups[j])
				{
					is_already_available = true;
					break;
				}

			if (is_already_available == false)
				enthalpy_defects_groups.push_back(list_enthalpy_defects[i]);
		}
		std::sort(enthalpy_defects_groups.begin(), enthalpy_defects_groups.end());

		// Search how many scalar dissipation rates are available
		std::vector<double> scalar_dissipation_rates_groups;
		for (unsigned int i = 0; i < list_stoichiometric_scalar_dissipation_rates.size(); i++)
		{
			bool is_already_available = false;
			for (unsigned int j = 0; j < scalar_dissipation_rates_groups.size(); j++)
				if (list_stoichiometric_scalar_dissipation_rates[i] == scalar_dissipation_rates_groups[j])
				{
					is_already_available = true;
					break;
				}

			if (is_already_available == false)
				scalar_dissipation_rates_groups.push_back(list_stoichiometric_scalar_dissipation_rates[i]);
		}
		std::sort(scalar_dissipation_rates_groups.begin(), scalar_dissipation_rates_groups.end());

		// Organize groups
		const double threshold_temperature = 400.;
		groups.resize(enthalpy_defects_groups.size());
		for (unsigned int i = 0; i < groups.size(); i++)
		{
			bool is_cold_flame_already_in = false;

			groups[i].enthalpy_defect = enthalpy_defects_groups[i];

			for (unsigned int j = 0; j < scalar_dissipation_rates_groups.size(); j++)
			{
				for (unsigned int k = 0; k < list_xml_files.size(); k++)
				{
					if ((list_enthalpy_defects[k] == enthalpy_defects_groups[i]) &&
						(list_stoichiometric_scalar_dissipation_rates[k] == scalar_dissipation_rates_groups[j]))
					{
						if (list_max_temperatures[k] > threshold_temperature)
						{
							groups[i].files.push_back(list_xml_files[k]);
							groups[i].scalar_dissipation_rates.push_back(list_stoichiometric_scalar_dissipation_rates[k]);
							groups[i].max_temperatures.push_back(list_max_temperatures[k]);
						}
						else
						{
							if (is_cold_flame_already_in == false)
							{
								is_cold_flame_already_in = true;
								groups[i].files.push_back(list_xml_files[k]);
								groups[i].scalar_dissipation_rates.push_back(list_stoichiometric_scalar_dissipation_rates[k]);
								groups[i].max_temperatures.push_back(list_max_temperatures[k]);
							}
						}

					}
				}
			}
		}
		for (unsigned int i = 0; i < groups.size(); i++)
		{
			std::cout << "Group: " << std::scientific << groups[i].enthalpy_defect / 1000. << std::endl;
			for (unsigned int j = 0; j < groups[i].files.size(); j++)
				std::cout << " * "	<< std::scientific << groups[i].scalar_dissipation_rates[j] << " "
									<< std::scientific << groups[i].max_temperatures[j] << " "
									<< groups[i].files[j].filename().string() << std::endl;
		}
	}

	void ReadHeaderSingleFlameletFromXMLFile(	const boost::filesystem::path file_name,
												double& stoichiometric_scalar_dissipation_rate,
												double& max_scalar_dissipation_rate,
												double& enthalpy_defect,
												double& max_temperature)
	{
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
				values >> stoichiometric_scalar_dissipation_rate;
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
				values >> max_scalar_dissipation_rate;
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
				values >> enthalpy_defect;
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the enthalpy-defect leaf");
		}

		// Maximum temperature
		{
			rapidxml::xml_node<>* indices_node = xml_main_input.first_node("opensmoke")->first_node("max-temperature");
			if (indices_node != 0)
			{
				std::stringstream values(indices_node->value());
				values >> max_temperature;
			}
			else
				OpenSMOKE::FatalErrorMessage("The " + file_name.string() + " is corrupted. Impossible to access the max-temperature leaf");
		}
	}
}
