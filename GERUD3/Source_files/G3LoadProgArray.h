#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
using namespace std;

// This file is part of Gerud3.  This version is version 3.1, the first working version
// of Gerud3.  Any changes should be indicated here.  This version is ANSI compliant and
// should compile on any standard C++ compiler.  This file was last altered on August 6, 2015.

class PAloader
{
public:
	bool KnownMother;
	int NumberEmbryos;
	int NumberLoci;
	string locusname[50];
	int NumberRuns;			  // Maximum number of runs is 100
	int EmbryosAssayed[100];  // Number of embryos assayed in the given run
	int EmbryosFromMale[100][10];  // Number of embryos in the progeny array from each of the 10 possible males
	char buffer[256];
	string motherid;
	int motherallele1[50];
	int motherallele2[50];
	int **embryoallele1;
	int **embryoallele2;
	bool LoadErr;

	int LoadProgArray(char* PAfile)
	{
		int i, j, k;
		size_t parsepos;
		int substrbegin, substrend;
		string tempstring;
		string tempsubstr;

		embryoallele1 = new int*[50];
		for (i = 0; i < 50; i++)
			embryoallele1[i] = new int[1000];

		embryoallele2 = new int*[50];
		for (i = 0; i < 50; i++)
			embryoallele2[i] = new int[1000];

		ifstream p_file;
		p_file.open(PAfile);
		if (!p_file.good())
			return 0;
		string embryoid;
		LoadErr = false;

		KnownMother = true;
		getline(p_file, tempstring);
		if (tempstring[0] == 'U' || tempstring[0] == 'u')
			KnownMother = false;

		getline(p_file, tempstring);
		substrend = tempstring.find_first_not_of("1234567890");
		tempsubstr = tempstring.substr(0, substrend);
		for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
			buffer[k] = tempsubstr[k];
		buffer[k] = '\0';
		NumberEmbryos = atoi(buffer);

		getline(p_file, tempstring);
		substrend = tempstring.find_first_not_of("1234567890");
		tempsubstr = tempstring.substr(0, substrend);
		for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
			buffer[k] = tempsubstr[k];
		buffer[k] = '\0';
		NumberLoci = atoi(buffer);

		if (NumberLoci > 50)
			NumberLoci = 50;

		if (NumberEmbryos > 1000)
			NumberEmbryos = 1000;

		// Load Locus Names

		getline(p_file, tempstring);
		parsepos = 0;
		for (i = 0; i < NumberLoci; i++)
		{
			while (tempstring[parsepos] != '!')
			{
				parsepos++;
				if (parsepos > tempstring.length())
				{
					parsepos = tempstring.length() - 1;
					break;
				}
			}
			substrbegin = parsepos;
			while (tempstring[parsepos] != ' ' && tempstring[parsepos] != ',' && tempstring[parsepos] != '\t' 
				&& tempstring[parsepos] != '\n' && tempstring[parsepos] != '\r' && tempstring[parsepos] != '\0')
			{
				parsepos++;
				if (parsepos > tempstring.length())
				{
					parsepos = tempstring.length() - 1;
					break;
				}
			}
			substrend = parsepos;
			if (substrend > substrbegin && substrend <= tempstring.length() && substrbegin < tempstring.length())
			{
				locusname[i] = tempstring.substr(substrbegin, substrend - substrbegin);
			}
			else
			{
				locusname[i] = "locusnameerror";
				LoadErr = true;
			}
		}

		// If the mother is known, then the next line is the maternal genotype
		int numbernumbers = 0;
		vector<int> numbers;
		int buffercounter;
		char buff1[4];
		char buff2[4];
		if (KnownMother)
		{
			getline(p_file, tempstring);
			tempsubstr.clear();
			parsepos = 0;
			while (tempstring[parsepos] != ',' && tempstring[parsepos] != ' ' && tempstring[parsepos] != '\t')
				parsepos++;
			motherid = tempstring.substr(0, parsepos);

			while (parsepos < tempstring.length())
			{
				while (tempstring[parsepos] != '0' && tempstring[parsepos] != '1' && tempstring[parsepos] != '2' && tempstring[parsepos] != '3'
					&& tempstring[parsepos] != '4' && tempstring[parsepos] != '5' && tempstring[parsepos] != '6' && tempstring[parsepos] != '7'
					&& tempstring[parsepos] != '8' && tempstring[parsepos] != '9')
				{
					parsepos++;
					if (parsepos > tempstring.length())
					{
						parsepos = tempstring.length();
						break;
					}
				}

				buffercounter = 0;
				while (tempstring[parsepos] == '0' || tempstring[parsepos] == '1' || tempstring[parsepos] == '2' || tempstring[parsepos] == '3'
					|| tempstring[parsepos] == '4' || tempstring[parsepos] == '5' || tempstring[parsepos] == '6' || tempstring[parsepos] == '7'
					|| tempstring[parsepos] == '8' || tempstring[parsepos] == '9')
				{
					buffer[buffercounter] = tempstring[parsepos];
					parsepos++;
					buffercounter++;
				}
				buffer[buffercounter] = '\0';
				if (buffercounter > 0 && buffercounter < 6)
				{
					numbers.push_back(atoi(buffer));
					numbernumbers++;
				}

				if (buffercounter == 6)
				{
					for (i = 0; i < 3; i++)
						buff1[i] = buffer[i];
					for (i = 0; i < 3; i++)
						buff2[i] = buffer[i + 3];
					buff1[3] = '\0';
					buff2[3] = '\0';
					numbers.push_back(atoi(buff1));
					numbers.push_back(atoi(buff2));
					numbernumbers = numbernumbers + 2;
				}
			} // end of while loop

			if (numbernumbers >= NumberLoci * 2)
			{
				for (i = 0; i < NumberLoci; i++)
				{
					motherallele1[i] = numbers[i * 2];
					motherallele2[i] = numbers[i * 2 + 1];
				}
			}
			else
			{
				LoadErr = true;
			}

		} // end of if known mother


		// Now load the embryo genotypes
		int embryocounter = 0;
		tempstring.clear();
		for (j = 0; j < NumberEmbryos; j++)
		{
			if (p_file.eof())
			{
				LoadErr = true;
				break;
			}
			numbers.clear();
			numbernumbers = 0;
			if (!p_file.eof())
				getline(p_file, tempstring);
			tempsubstr.clear();
			parsepos = 0;
			while (tempstring[parsepos] != ',' && tempstring[parsepos] != ' ' && tempstring[parsepos] != '\t')
				parsepos++;
			embryoid = tempstring.substr(0, parsepos);

			while (parsepos < tempstring.length())
			{
				while (tempstring[parsepos] != '0' && tempstring[parsepos] != '1' && tempstring[parsepos] != '2' && tempstring[parsepos] != '3'
					&& tempstring[parsepos] != '4' && tempstring[parsepos] != '5' && tempstring[parsepos] != '6' && tempstring[parsepos] != '7'
					&& tempstring[parsepos] != '8' && tempstring[parsepos] != '9')
				{
					parsepos++;
					if (parsepos > tempstring.length())
					{
						parsepos = tempstring.length();
						break;
					}
				}

				buffercounter = 0;
				while (tempstring[parsepos] == '0' || tempstring[parsepos] == '1' || tempstring[parsepos] == '2' || tempstring[parsepos] == '3'
					|| tempstring[parsepos] == '4' || tempstring[parsepos] == '5' || tempstring[parsepos] == '6' || tempstring[parsepos] == '7'
					|| tempstring[parsepos] == '8' || tempstring[parsepos] == '9')
				{
					buffer[buffercounter] = tempstring[parsepos];
					parsepos++;
					buffercounter++;
				}
				buffer[buffercounter] = '\0';
				if (buffercounter > 0 && buffercounter < 6)
				{
					numbers.push_back(atoi(buffer));
					numbernumbers++;
				}

				if (buffercounter == 6)
				{
					for (i = 0; i < 3; i++)
						buff1[i] = buffer[i];
					for (i = 0; i < 3; i++)
						buff2[i] = buffer[i + 3];
					buff1[3] = '\0';
					buff2[3] = '\0';
					numbers.push_back(atoi(buff1));
					numbers.push_back(atoi(buff2));
					numbernumbers = numbernumbers + 2;
				}
			}

			if (numbernumbers >= NumberLoci * 2)
			{
				for (i = 0; i < NumberLoci; i++)
				{
					embryoallele1[i][j] = numbers[i * 2];
					embryoallele2[i][j] = numbers[i * 2 + 1];
				}
			}
			else
			{
				LoadErr = true;
			}


	
		} // end of j loop



		/*
		getline(p_file, tempstring);
		substrend = tempstring.find_first_not_of("0123456789");
		tempsubstr = tempstring.substr(0, substrend);
		for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
			buffer[k] = tempsubstr[k];
		buffer[k] = '\0';
		NumberIterations = atoi(buffer);

		if (NumberIterations > 10000)
			NumberIterations = 10000;

		cout << NumberIterations << " Iterations\n";

		getline(p_file, tempstring);
		substrend = tempstring.find_first_not_of("0123456789");
		tempsubstr = tempstring.substr(0, substrend);
		for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
			buffer[k] = tempsubstr[k];
		buffer[k] = '\0';
		NumberRuns = atoi(buffer);

		if (NumberRuns > 100)
			NumberRuns = 100;

		cout << NumberRuns << " Runs\n";

		getline(p_file, tempstring);

		for (i = 0; i < NumberRuns; i++)
		{
			getline(p_file, tempstring);
			substrbegin = tempstring.find_first_not_of("0123456789") + 1;

			substrend = tempstring.find_first_not_of("0123456789", substrbegin);
			tempsubstr = tempstring.substr(substrbegin, substrend);
			for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
				buffer[k] = tempsubstr[k];
			buffer[k] = '\0';
			EmbryosAssayed[i] = atoi(buffer);
			substrbegin = substrend + 1;

			for (j = 0; j < 10; j++)
			{
				substrend = tempstring.find_first_not_of("0123456789", substrbegin);
				tempsubstr = tempstring.substr(substrbegin, substrend);
				for (k = 0; k < static_cast<int>(tempsubstr.length()); k++)
					buffer[k] = tempsubstr[k];
				buffer[k] = '\0';
				EmbryosFromMale[i][j] = atoi(buffer);
				substrbegin = substrend + 1;
			} // j
		} // i

		for (i = 0; i < NumberRuns; i++)
		{
			cout << "Run " << i + 1 << "\t" << EmbryosAssayed[i] << " Embryos Assayed\t";
			for (j = 0; j < 10; j++)
			{
				cout << EmbryosFromMale[i][j] << " from Male" << j + 1 << "\t";
			}
			cout << "\n";
		}*/

		p_file.close();

		if (LoadErr)
			return 0;
		else
			return 1;
	} 

	void FreeMemory()
	{
		int i;
		for (i = 0; i < 50; i++)
			delete[] embryoallele1[i];
		delete[] embryoallele1;

		for (i = 0; i < 50; i++)
			delete[] embryoallele2[i];
		delete[] embryoallele2;

		embryoallele1 = NULL;
		embryoallele2 = NULL;
	}
};
