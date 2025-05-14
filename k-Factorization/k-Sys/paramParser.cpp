#include <fstream>
#include <algorithm>
#include <regex>
#include "TopGun.h"

using namespace std;

bool case_insensitive_match(char a, char b) {
    return std::tolower(static_cast<unsigned char>(a)) ==
        std::tolower(static_cast<unsigned char>(b));
}

size_t case_insensitive_find(const std::string& haystack, const std::string& needle) {
    const auto it = std::search(
        haystack.begin(), haystack.end(),
        needle.begin(), needle.end(),
        case_insensitive_match
    );
    
	if (it != haystack.end())
		return static_cast<size_t>(std::distance(haystack.begin(), it));
	else
		return string::npos;
}

bool iequals(const std::string& a, const std::string& b) {
	return a.size() == b.size() &&
		std::equal(a.begin(), a.end(), b.begin(), [](char c1, char c2) {
		return std::tolower(static_cast<unsigned char>(c1)) ==
			std::tolower(static_cast<unsigned char>(c2));
			});
}

static size_t getValue(string& tmp, size_t* pPos) {
	ltrim(tmp);
	if (!tmp.length()) {
		*pPos = string::npos;
		return 0;
	}

	if (tmp.front() != '"' || tmp.back() != '"') {
		size_t posE = tmp.find(' ');
		const size_t posE1 = tmp.find(';');
		if (posE1 != string::npos && posE1 < posE)
			posE = posE1;

		*pPos = posE != string::npos ? posE + 1 : string::npos;
		if (posE != string::npos)
			tmp = tmp.substr(0, posE);
	}
	else
		tmp = tmp.substr(1, tmp.length() - 2);

	return tmp.length();
}

static int getInteger(const string& str, size_t* pPos) {
	string tmp = str.substr(*pPos + 1);
	if (!getValue(tmp, pPos))
		return 1;

	const auto symb = tmp[0];
	if (isdigit(symb) || symb == '-' || symb == '+')
		return atoi(tmp.c_str());

	return tmp == "YES" ? 1 : 0;
}

static size_t find(const string& str, const char* strValue)
{
	const size_t pos = case_insensitive_find(str, strValue);
	return pos != string::npos ? pos + strlen(strValue) : pos;
}

template<typename T>
void setDefaultValue(int* pValue) {
	*pValue = 1;
}

template<typename T>
void setDefaultValue(string** pValue) {
	if (*pValue)
		**pValue = "";
}

template<typename T>
void setDefaultValue(tchar** pValue) {
	printfRed("No default parameter defined for tchar **");
}

template <typename T>
bool setParameter(int* pValue, const string& str, size_t* pos) {
	*pValue = getInteger(str, pos);
	return true;
}

template <typename T>
bool setParameter(string** pValue, const string& str, size_t* pos) {
	string tmp = str.substr(*pos + 1);
	if (getValue(tmp, pos)) {
		if (!*pValue)
			*pValue = new string();

		**pValue = tmp;
	}
	else
		setDefaultValue<string**>(pValue);

	return true;
}

template <typename T>
bool setParameter(tchar** pValue, const string& str, size_t* pos) {
	auto pStr = str.c_str();
	auto* chptr = strrchr(pStr, '}');
	if (!chptr)
		return false;

	string* groups[10];
	auto last = chptr - pStr;
	chptr = strchr(pStr, '{');
	if (!chptr)
		return false;

	auto first = chptr - pStr;
	string tmp1 = str.substr(first + 1, last - first - 1);
	string tmp = regex_replace(tmp1, regex("\\s"), "");
	std::regex patternA("([0-9]+,)+[0-9]+");
	std::regex patternB("[0-9]+");
	int ngrp = 0;
	while (tmp.length()) {
		pStr = tmp.c_str();
		chptr = strchr(pStr, '{');
		if (!chptr)
			break;

		first = chptr - pStr;
		chptr = strchr(pStr, '}');
		if (!chptr)
			return false;

		last = chptr - pStr;
		assert(ngrp < countof(groups));
		groups[ngrp] = new string(tmp.substr(first + 1, last - first - 1));
		if (!std::regex_match(*groups[ngrp], patternA) && !std::regex_match(*groups[ngrp], patternB))
			return false;

		ngrp++;
		tmp = tmp.substr(last + 1);
	}

	delete[] * pValue;
	const auto len = MAX_CYCLES_PER_SET * ngrp;
	*pValue = new tchar[len + ngrp + 1];
	**pValue = ngrp;
	auto pTo = *pValue + 1;
	memset(pTo, 0, len);
	bool retVal = true;
	int val;
	for (int j = 0; j < ngrp; j++) {
		auto& tmp = *groups[j];
		trim(tmp);
		int k = 0;
		while (tmp.length()) {
			if (isdigit(tmp[0])) {
				pStr = tmp.c_str();
				chptr = strchr(pStr, ',');
				if (chptr) {
					val = atoi(tmp.substr(0, chptr - pStr).c_str());
					tmp = tmp.substr(chptr - pStr + 1);
					ltrim(tmp);
				}
				else {
					val = atoi(tmp.c_str());
					tmp.clear();
				}

				if (k < MAX_CYCLES_PER_SET && val) {
					*(pTo + k++) = val;
				}
				else
					retVal = false;
			}
			else
				retVal = false;
		}

		pTo += MAX_CYCLES_PER_SET;
		delete groups[j];
	}

	return retVal;
}

template <typename T>
static int getParam(const string& str, const char* pKeyWord, T* pValue, size_t* pPos = nullptr) {
	size_t pos = find(str, pKeyWord);
	if (pos == string::npos)
		return 0;		// keyWord was not found

	const auto flg = pos == str.length();
	if (!flg && str[pos] != ' ' && (str[pos] != '=' && str[pos] != ':'))
		return 0;

	string tmp = str.substr(0, pos);
	if (!iequals(tmp, pKeyWord))
		return 0;

	tmp = str.substr(pos);
	ltrim(tmp);
	if (flg || (tmp[0] != '=' && tmp[0] != ':')) {
		// the value is determined by the presence of the keyword
		setDefaultValue<T>(pValue);
		if (!tmp.length())
			return -1;	// there is nothing after keyWords
	}
	else {
		if (!setParameter<T>(pValue, str, &pos)) {
			printfRed("Parameter is not set as expected: '%s'\n", str.c_str());
			exit(1);
		}
	}

	if (pPos)
		*pPos = pos;

	return 1;
}

int getParameter(string& line, const paramDescr* par, int nDescr, kSysParam& param) {
	const auto pos = line.find_first_of("=:");
	if (pos != string::npos) {
		auto beg = line.substr(0, pos);
		auto end = line.substr(pos + 1);
		rtrim(beg);
		ltrim(end);
		line = beg + "=" + end;
	}

	static const char* jobTask[] = { "RUN_JOB", "END_JOB" };
	for (int j = 0; j < 2; j++) {
		if (case_insensitive_find(line, jobTask[j]) != string::npos)
			return j + 1;
	}

	auto j = nDescr;
	while (j--) {
		auto paramNames = (par + j)->paramNames;
		bool rc = false;
		int i = (par + j)->numParams;
		while (!rc && i--) {
			if (j == 0)
				rc = getParam<int>(line, paramNames[i], param.val + i);
			else
				if (j == 1)
					rc = getParam<string*>(line, paramNames[i], &param.strVal[i]);
				else
					rc = getParam<tchar*>(line, paramNames[i], &param.u1fCycles[i]);
		}

		if (i >= 0)
			return 3;
	}

	if (j < 0)
		printfYellow("Sorry, the parameter for keyword '%s' was not found\n", line.c_str());

	return 0;
}

bool getParameters(ifstream& infile, const paramDescr* par, int nDescr, kSysParam& param, bool& firstSet, bool& endJob) {
	bool retVal = false;
	string line;
	while (getline(infile, line)) {		// For all the lines of the file
		trim(line);
		size_t pos = line.find("//");
		if (pos != string::npos)
			line = line.substr(0, pos);	// deleting a comment at the end of a line

		if (!line.size() || line[0] == ';')
			continue;					// Skip line if it is a comment OR empty

		if (line[0] == '/' && line[1] == '*') {
			// Commented out part of the input file
			// Let's find the closing 
			do {
				const size_t pos = line.find("*/");
				if (pos != string::npos) {
					line = line.substr(pos + 2);
					break;
				}
			} while (getline(infile, line));

			trim(line);
			if (!line.length())
				continue;

		}

		if (firstSet) {
			firstSet = !line.starts_with("PrintMatrices");
			if (firstSet) {
				firstSet = line != SECTION_PARAM_MAIN;
				continue;
			}
		}

		if (line == SECTION_PARAM_STRINGS || line == SECTION_PARAM_U1F_CONF)
			continue;

		//transform(line.begin(), line.end(), line.begin(), ::toupper);

		string subLine;
		while (!line.empty()) {
			const auto pos = line.find(";");
			if (pos != string::npos) {
				subLine = line.substr(0, pos);
				line = line.substr(pos + 1);
			}
			else {
				subLine = line;
				line = "";
			}

			trim(subLine);
			switch (getParameter(subLine, par, nDescr, param)) {
			case 2: endJob = true;
			case 1: return retVal;
			case 3: retVal = true;		// at least one parameter has been changed
			default: continue;
			}
		}
	}

	return retVal;
}

string getUF(const tchar* pU1F) {
	char group[32];
	tchar val;
	int k = 0;
	string str("{");
	auto nu = *pU1F++;
	while (nu--) {
		group[k] = '{';
		for (int i = 0; i < MAX_CYCLES_PER_SET && (val = pU1F[i]); i++) {
			int j = val < 10 ? 1 : val < 100 ? 10 : 100;
			do {
				group[++k] = '0' + val / j;
				val -= (val / j) * j;
			} while (j /= 10);

			group[++k] = ',';
		}
		pU1F += MAX_CYCLES_PER_SET;
		group[k] = '}';
		group[k + 1] = '\0';
		str.append(group);
		group[0] = ',';
		k = 1;
	}
	str.append("}");
	return str;
}
