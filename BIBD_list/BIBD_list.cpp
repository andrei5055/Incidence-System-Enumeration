﻿
#include <iostream>
#include <algorithm>
#include <climits>
#include <sstream>
#include <string.h>

#ifndef WIN
#define sprintf_s(x, y, ...) sprintf(x, __VA_ARGS__)
#define fopen_S(x, y, z)     x = fopen(y, z)
#define FPRINTF(x, y)        fprintf(x, "%s", y)
#define strcpy_s(x, y, z)    strcpy(x, z)
#define memcpy_s(x, y, ...)  memcpy(x, __VA_ARGS__)
#else
#define fopen_S(x, y, z)     fopen_s(&x, y, z)
#define FPRINTF(x, y)        fprintf(x, y)
#endif

using namespace std;

#define SQUOT(x) "'" << x << "'" 

bool is_number(const std::string& s) {
    string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

typedef struct {
    std::stringstream& buffer;
    FILE* file;
    bool coutFlag;
} out_struct;

template <typename T>
void output_func(const out_struct& out, T t) {
    out.buffer << t << std::endl;
    if (out.coutFlag)
        cout << out.buffer.str();

    if (out.file)
        FPRINTF(out.file, out.buffer.str().c_str());

    out.buffer.str("");
}

template<typename T, typename... Args>
void output_func(const out_struct& out, T t, Args... args) { // recursive variadic function
    out.buffer << t;
    output_func(out, args...);
}

class opt_descr;

typedef int (*check_fn)(const opt_descr& opt, int v, int b, int r, int k, int l, string& comment);
typedef void (*report_fn)(opt_descr& opt, const out_struct& out, int total);
typedef void (*constr_fn)(opt_descr& opt);

class opt_descr {
public:
    const string name;
    check_fn check_func;
    report_fn report_func;
    int limits[2];
    int intValue;
    string strValue;
    constr_fn constr_func[2];
    void* extra_param;
};

typedef enum {
    no_comment,        // Add to output without comment
    with_comment,      // Add to output with comment
    nothing            // Not including this set of parameters in output list
} output;

int set_option(const string& inOpt, const char *arg, opt_descr *opts, int num_opts, string &error) {
    if (inOpt[0] != '-')
        return 1;   // not an option

    std::stringstream buffer;
    for (int i = 0; i < num_opts; i++) {
        auto& opt = opts[i];
        auto pos = inOpt.find(opt.name, 1);
        if (pos == string::npos)
            continue;

        if (inOpt[pos += opt.name.length()] != '=') {
            buffer << "Invalid option settings: " SQUOT(arg);
            buffer << "\nExpected '-" << opt.name << "= val', where 'val' is one of {";
            for (int i = opt.limits[0]; i <= opt.limits[1]; i++)
                buffer << i;

            buffer << "}";
        }
        else {
            const string valStr = inOpt.substr(++pos);
            if (valStr.empty())
                continue;

            if (is_number(valStr)) {
                const int val = atoi(valStr.c_str());
                if (opt.limits[0] <= val && val <= opt.limits[1])
                    opt.intValue = val;
                else
                    buffer << "Invalid value " SQUOT(val) " for option " SQUOT(arg);
            }
            else {
                if (opt.intValue < opt.limits[0]) {
                    // string type option expected
                    opt.strValue = string(arg).substr(pos);
                } else
                    buffer << "Value for option -" << opt.name << " defined by " SQUOT(arg) " supposed to be a number, got " SQUOT(valStr);
            }
        }

        if (!(error = buffer.str()).empty())
            break;
    }

    return 0;
}

int gcd(int a, int b, int maxVal = INT_MAX) {
    int result = 1 + min(maxVal, min(a, b));
    while (--result > 1 && (a % result || b % result));
    return result;              // return gcd of and b
}

void init_Cntr(opt_descr& opt) {
    opt.extra_param = new int [3];
    memset(opt.extra_param, 0, 3 * sizeof(int));
}

int check_BRC(const opt_descr &opt, int v, int b, int r, int k, int λ, string& comment) {
    static int oddPrime[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };

    const auto brc_condition = opt.intValue;
    if (v != b)
       return brc_condition == 2? output::nothing : output::no_comment;

    // Check of Bruck-Ryser-Chowla conditions
    int num = 0;
    if (!(v % 2)) {
        // v is even, k - lambda should be a perfect square;
        // Use the fact that every perfect square is a sum of consequtive
        // odd numbers starting from 1;
        int odd = 1;
        num = k - λ;
        while ((num -= odd) > 0)
            odd += 2;
    }
    else {
        // We use the test from pp.25-26 of https://www.math.uwaterloo.ca/~cgodsil/pdfs/bigGeom.pdf
        // with an added check that the coefficients  (A, B, C) are square free.
        // Trying to find non-zero solution for the equation: A*x^2 + b*y^2 + C*z^2 = 0
        int m[] = { k - λ, (v >> 1) % 2 ? -λ : λ, -1 };  // { A, B, C }
        // Let us make these coefficients pairwise coprime.
        int cntr = 3;   // the number of pairs that still need to be checked for being coprime
        int j = 2;
        while (true) {
            const int j1 = j ? j - 1 : 2;
            const int j2 = 3 - j - j1;
            const int N = gcd(abs(m[j1]), abs(m[j2]));
            if (N > 1) {
                m[j1] /= N;
                m[j2] /= N;
                m[j] *= N;
                cntr = 3;
            }
            else
                if (!(--cntr))
                    break;

            if (!j--)
                j = 2;
        }

        for (int j = 0; j < 3; j++) {
            // make the number free of 2^2
            auto n = abs(m[j]);
            while (n >= 4 && ((n >> 2) << 2) == n) {
                n >>= 2;
            }

            // Make the number free of squaares of odd primes
            for (const auto p : oddPrime) {
                const auto p_quare = p * p;
                while (!(n % p_quare))
                    n /= p_quare;

                // Compare with the minimum of next possible value of p_quare
                // p_quare_next = p_next^2 <= (p + 2)^2 = p^2 + 2 * p + 2^2
                if (n < p_quare + (p<<1) + 4)
                    break;
            }

            m[j] = m[j] > 0 ? n : -n;
        }

        // For symmetrcal BIBD A, B, C cannot have the same sign,
        // Because C = -1 < 0 and A = k - lambda > 0, we don't need to check that mot all (A,B,C) have the same sign
        for (int j = 0; j < 3; j++) {
            const auto N = abs(m[j]);
            if (N < 3)
                continue;

            for (const auto p : oddPrime) {
                if (N % p)
                    continue;

                // val should be a square modulo p
                const int j1 = j ? j - 1 : 2;
                auto val = -m[j1] * m[3 - j - j1];
                val = val >= 0 ? val % p : p - (-val) % p;
                int q = -1;
                while (++q < p && (q * q) % p != val);
                if (q == p) {
                    // there are no integer solutions for A*x^2 + b*y^2 + C*z^2 = 0
                    num = 1;
                    j = 2;    // to exit the outer loop
                    break;
                }
            }
        }
    }

    if ((!brc_condition && num < 0) || (brc_condition == 2 && !num))
        return output::nothing;

    comment = opt.name + (num ? "-" : "+");
    auto* pntr = static_cast<int*>(opt.extra_param);
    ++*(pntr + (num? 1 : 0));
    return output::with_comment;
}

int check_Residual(const opt_descr& opt, int v, int b, int r, int k, int λ, string& comment) {
    static opt_descr opt_BRC = { "BRC", check_BRC, nullptr, { 1, 1 }, 1, ""};
    if (!opt_BRC.extra_param)
        opt_BRC.extra_param = opt.extra_param;

    if (λ > 2 || r != k + λ)
        return output::no_comment;

    string commentBRC;
    const auto retVal = check_BRC(opt_BRC, v + r, b + 1, r, k + λ, λ, commentBRC);
    if (retVal == output::with_comment)
        comment = commentBRC.back() + opt.name;
    else
        comment = "Something wrong with BRC";

    return output::with_comment;
}


int check_Simplicity(const opt_descr& opt, int v, int b, int r, int k, int λ, string& comment) {
    auto* pntr = static_cast<int*>(opt.extra_param);
    static const char* pComment[] = { "S", "S-", "!S"};
    int idx = 0;
    // Using conditions by Mann: the maximum multiplicity of a block m ≤ b/v
    if (b >= 2 * v && λ > 1) {
        if (k != 3 || λ <= v - 2) {
            // Calculate  gcd(B, r, λ)
            int val[] = { b, r, λ };
            auto a = λ;
            auto gcd_brλ = a;
            int cntr = 3;
            int j = 0;
            while (true) {
                int j1 = j < 2 ? j + 1 : 0;
                gcd_brλ = gcd(val[j1], val[3 - j++ - j1], a);
                if (gcd_brλ == 1)
                    break;

                if (a > gcd_brλ) {
                    a = gcd_brλ;
                    cntr = 3;
                }
                else
                    if (!--cntr)
                        break;
            }

            // By Lint and Ryser When b = 2 * v, multiplicity m >= 2 only when gcd(b, r, λ) is even
            if (b > 2 * v || !(gcd_brλ % 2)) {
                // Result by Lint and https://webspace.maths.qmul.ac.uk/l.h.soicher/repeated11.pdf
                // If 2 ≤ k ≤ v / 2, λ(v − 1) = r(k − 1), vr = bk, λ > 1, gcd(b, r, λ) = 1 and b > 2v, 
                // then there exist BIBD with repeated block
                if (gcd_brλ == 1 && b > 2 * v && k <= 4) {
                    idx = 1;  // Multiple blocks exist
                }
                else
                    return output::no_comment;
            }
        } else {
            // When we are here k = 3 and λ + 2 > v
            // Any such design will have a duplicated blocks
            // Proof: Concider the blocks having λ common point. Any other poins cannot below to two of these blocks, 
            // we would have two identical blocks. Thus this design will contain at least  λ + 2 >= v points
            idx = 2;
        }
    }

    ++*(pntr + idx);
    comment += pComment[idx];
    return output::with_comment;
}

void report_BRC(opt_descr& opt, const out_struct& out, int total) {
    auto* pntr = static_cast<int*>(opt.extra_param);
    output_func(out, "Bruck-Ryser-Chowla conditions results: out of ", *pntr + *(pntr + 1), " tests ",
        *pntr, " positive and ", *(pntr + 1), " negative.");

    delete[](int*)opt.extra_param;
}

void report_Residual(opt_descr& opt, const out_struct& out, int total) {
    auto* pntr = static_cast<int*>(opt.extra_param);
    output_func(out, "Residual condition was tested on ", *pntr + *(pntr + 1), " sets of parameters: ",
        *pntr, " positive and ", *(pntr + 1), " negative.");

    delete[](int*)opt.extra_param;
}

void report_Simplicity(opt_descr& opt, const out_struct& out, int total) {
    auto* pntr = static_cast<int*>(opt.extra_param);
    output_func(out, "Among all design sets: ", *pntr, " are simple, ",
        *(pntr + 1), " not simple, ", *(pntr + 2), " \"anti-simple\", ",
        total - (*pntr + *(pntr + 1) + *(pntr + 2)), " unknown.");

    delete[](int*)opt.extra_param;
}

int main(int argc, char* argv[])
{
    const wchar_t param[] = { 'V', 'B', 'R', 'K', 'L' };
    const char delimiter[] = { '[', '(', ',', ':', ']', ')' };
    int minVal[5] = { 0, 0, 3, 3, 1 };
    int maxVal[5] = { INT_MAX, INT_MAX, 41, INT_MAX, INT_MAX };
    opt_descr opts[] = { {"F", nullptr, nullptr, {}, -1, ""},
                         {"S", check_Simplicity, report_Simplicity, {-1, 2}, 1, "", init_Cntr},
                         {"BRC", check_BRC, report_BRC, {-1, 2}, 1, "", init_Cntr},
                         {"O", check_Residual, report_Residual, {-1,1}, 1, "", init_Cntr},
                       };

    const auto num_opts = sizeof(opts) / sizeof(opts[0]);

    int cntr = -1;  // by default this counter will not be used
    size_t prevPos, pos;
    string error;
    for (int i = 1; i < argc; i++) {
        string arg(argv[i]);
        transform(arg.begin(), arg.end(), arg.begin(), ::toupper);

        if (!set_option(arg, argv[i], opts, num_opts, error)) {
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

        if (j != 2 || minVal[2] != 3) {
            // We found a parameter with new limits that is not 'r'
            // OR it 'r', but its lower limit is different from its default value.
            cntr = 0;   // Activating an additional counter
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

    for (int j = 1; j < num_opts; j++) {
        const constr_fn constr_func = opts[j].constr_func[0];
        if (constr_func)
            (*constr_func)(opts[j]);
    }

    FILE* outFile = NULL;
    bool coutFlag = true;
    if (!opts[0].strValue.empty()) {
        const auto* pFileName = opts[0].strValue.c_str();
        coutFlag = (*pFileName == '+');
        if (coutFlag || *pFileName == '-')
            pFileName++;

        fopen_S(outFile, pFileName, "w");
    }

    int total = 0;
    char buff[256], *pBuff = buff;
    const int rMax = maxVal[2];
    for (int r = 3; r <= rMax; r++) {
        for (int k = 3; k <= r; k++) {
            const auto λ_max = r * (k - 1) / (2 * k - 1);
            for (int λ = 1; λ <= λ_max; λ++) {
                const auto vMinus1 = r * (k - 1) / λ;
                if (r != λ * vMinus1 / (k - 1))
                    continue;

                int v, b;
                const auto a = (v = vMinus1 + 1) * r;
                if ((a / k) * k != a)
                    continue;

                b = a / k;
                total++;
                const int val[] = { v, b, r, k, λ };
                int j = sizeof(val) / sizeof(val[0]);
                while (j--) {
                    if (val[j] < minVal[j] || val[j] > maxVal[j])
                        break;
                }

                if (j >= 0)
                    continue; // j-th parameter is not in the predescribed range
 
                string comment, tmp;
                for (int j = 1; j < num_opts; j++) {
                    const check_fn check_func = opts[j].check_func;
                    if (check_func) {
                        switch ((*check_func)(opts[j], v, b, r, k, λ, tmp)) {
                            case with_comment:  comment += " " + tmp;
                            case no_comment:    break;
                            case nothing:       continue;
                        }
                    }
                }

                pBuff = buff;
                if (cntr >= 0)
                    pBuff += sprintf_s(pBuff, sizeof(buff), "%4d ", ++cntr);

                sprintf_s(pBuff, sizeof(buff) - (pBuff - buff), "#%4d: %4d %4d %2d %2d %2d  %s\n", total, v, b, r, k, λ, comment.c_str());
                if (outFile)
                    FPRINTF(outFile, buff);

                if (coutFlag)
                    cout << buff;
            }
        }
    }

    std::stringstream buffer;
    buffer << std::endl;
    out_struct out = { buffer, outFile, coutFlag };

    for (int j = 1; j < num_opts; j++) {
        const report_fn report_func = opts[j].report_func;
        if (report_func)
            (*report_func)(opts[j], out, cntr < 0? total : cntr);
    }

    if (outFile)
        fclose(outFile);

    for (int j = 1; j < num_opts; j++) {
        const constr_fn destr_func = opts[j].constr_func[1];
        if (destr_func)
            (*destr_func)(opts[j]);
    }
}
