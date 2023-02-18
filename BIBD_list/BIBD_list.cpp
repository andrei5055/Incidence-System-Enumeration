#include <iostream>
#include <algorithm>
#include <climits>
#include <sstream>


#ifndef WIN
#define sprintf_s(x, y, ...) sprintf(x, __VA_ARGS__)
#define strcpy_s(x, y, z)    strcpy(x, z)
#define memcpy_s(x, y, ...)  memcpy(x, __VA_ARGS__)
#endif

using namespace std;

#define SQUOT(x) "'" << x << "'" 

bool is_number(const std::string& s) {
    string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

typedef struct {
    const char* name;
    int limits[2];
    int* pOptValue;
} opt_descr;

int set_option(const string& inOpt, const char *arg, const opt_descr *opts, int num_opts, string &error) {
    if (inOpt[0] != '-')
        return 1;   // not an option

    std::stringstream buffer;
    for (int i = 0; i < num_opts; i++) {
        const auto& opt = opts[i];
        auto pos = inOpt.find(opt.name, 1);
        if (pos != string::npos) {
            if (inOpt[pos += strlen(opt.name)] != '=') {
                buffer << "Invalid option settings: " SQUOT(arg);
                buffer << "\nExpected '-" << opt.name << "= val', where 'val' is one of {";
                for (int i = opt.limits[0]; i <= opt.limits[1]; i++)
                    buffer << i;

                buffer << "}";
            }
            else {
                const string valStr = inOpt.substr(++pos);
                if (is_number(valStr)) {
                    const int val = atoi(valStr.c_str());
                    if (opt.limits[0] <= val && val <= opt.limits[1])
                        *opt.pOptValue = val;
                    else
                        buffer << "Invalid value " SQUOT(val) " for option " SQUOT(arg);
                }
                else
                    buffer << "Value for option -" << opt.name << " defined by " SQUOT(arg) " supposed to be a number, got " SQUOT(valStr);
            }
            break;
        }
    }

    error = buffer.str();
    return 0;
}

int main(int argc, char* argv[])
{
    const char param[] = { 'V', 'B', 'R', 'K', 'L' };
    const char delimiter[] = { '[', '(', ',', ':', ']', ')' };
    int minVal[5] = { 0, 0, 3, 3, 1 };
    int maxVal[5] = { INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX };
    int brc_condition = 1;  // the default value for checking and reporting 
                            // the fulfilment of Bruck-Ryser-Chowla conditions.
    opt_descr opts[] = { {"BRC", {-1, 2}, &brc_condition} };
    int cntr = argc > 1 ? 0 : -1;
    size_t prevPos, pos;
    string error;
    for (int i = 1; i < argc; i++) {
        string arg(argv[i]);
        transform(arg.begin(), arg.end(), arg.begin(), ::toupper);

        if (!set_option(arg, argv[i], opts, sizeof(opts) / sizeof(opts[1]), error)) {
            if (!error.empty()) {
                cout << error << "\n";
                exit(1);
            }
            continue;
        }
    
        // Parameter settings
        const auto symb = arg[0];
        int j = sizeof(param) / sizeof(param[0]);
        while (j-- && symb != param[j]);
        if (j < 0) {
            cout << "Can't parse parameter " SQUOT(argv[i]);
            exit(1);
        }

        prevPos = 0;
        int nPrev = 0;
        for (int m = 0; m < 6; m += 2) {    // loop for delimiters
            int n = 2;                      // loop for delimiters type
            while (n-- && (pos = arg.find(delimiter[m + n], prevPos)) == string::npos);
            if (n < 0) {
                cout << "Cannot find either of the two expected delimiters (" SQUOT(delimiter[m])
                    ", " SQUOT(delimiter[m]) ") in " SQUOT(argv[i]);
                exit(1);
            }

            if (m) {
                if (prevPos != pos) {
                    const string valStr = arg.substr(prevPos, pos - prevPos);
                    int val = atoi(valStr.c_str());
                    if (val > 0) {
                        if (m == 2) {
                            if (nPrev)
                                val++;

                            if (j < 2 || j > 3 || val >= minVal[j])
                                minVal[j] = val;
                            else {
                                cout << "The minimum value of the parameter " SQUOT(symb)
                                     " cannot be less than " << minVal[j];
                            }
                        } 
                        else {
                            maxVal[j] = val - (n ? 1 : 0);
                            if (maxVal[j] < minVal[j]) {
                                cout << "The maximum value (" << maxVal[j] << ") for parameter " 
                                     SQUOT(symb) " is less than its minimum value (" << minVal[j];
                            }
                        }
                    }
                    else {
                        cout << "Ignoring non-positive value " << val << " for parameter "
                            SQUOT(symb) " found in " SQUOT(argv[i]);
                    }
                }
                else {
                    // empty parameter value
                    // we don't need to change minVal/maxVal for parameter symb
                }
            }

            nPrev = n;
            prevPos = pos + 1;
        }
    }

    int total = 0;
    char buff[256], *pBuff = buff;
    const int rMax = 41;
    for (int r = 3; r <= rMax; r++) {
        for (int k = 3; k <= r; k++) {
            const auto lMax = r * (k - 1) / (2 * k - 1);
            for (int l = 1; l <= lMax; l++) {
                const auto vMinus1 = r * (k - 1) / l;
                if (r != l * vMinus1 / (k - 1))
                    continue;

                int v, b;
                const auto a = (v = vMinus1 + 1) * r;
                if ((a / k) * k != a)
                    continue;

                b = a / k;
                total++;
                const int val[] = { v, b, r, k, l };
                int j = sizeof(val) / sizeof(val[0]);
                while (j--) {
                    if (val[j] < minVal[j] || val[j] > maxVal[j])
                        break;
                }

                if (j >= 0)
                    continue; // j-th parameter is not in the predescribed range
 
                string comment;
                if (brc_condition == 2 && v != b)
                    continue;

                if (brc_condition >= 0 && v == b) {
                    // Check of Bruck-Ryser-Chowla conditions
                    comment += opts[0].name;
                    int num = 0;
                    if (!(v & 1)) {
                        // v is even, k - lambda should be a perfect square;
                        // Use the fact that every perfect square is a sum of consequtive 
                        // odd numbers starting from 1;
                        int odd = 1;
                        num = k - l;
                        while ((num -= odd) > 0)
                            odd += 2;
                    }
                    else {
                        // So far the check of condition for odd v is not imlemented
                        // we will mark it as "passed"
                    }

                    if (!brc_condition && num < 0 || brc_condition == 2 && !num)
                        continue;   // Not including this set of parameters in output list

                    comment += num ? "-" : "+";
                }

                if (cntr >= 0)
                    pBuff += sprintf_s(pBuff = buff, sizeof(buff), "%4d ", ++cntr);

                sprintf_s(pBuff, sizeof(buff) - (pBuff - buff), "#%4d: %4d %4d %2d %2d %2d  %s\n", total, v, b, r, k, l, comment.c_str());
                cout << buff;
            }
        }
    }
}

