
#include <iostream>
#include <algorithm>
#include <climits>
#include <sstream>
#include <string.h>
#include <assert.h>

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

typedef enum {
    no_comment,        // Add to output without comment
    with_comment,      // Add to output with comment
    adjust_comments,   // Add output to previous comment
    replace_comments,  // Replace previous comments by a new one
    nothing            // Not including this set of parameters in output list
} output;

typedef output (*check_fn)(const opt_descr& opt, int v, int b, int r, int k, int l, string& comment);
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

int set_option(const string& inOpt, const char *arg, opt_descr *opts, int num_opts, string &error) {
    if (inOpt[0] != '-')
        return 1;   // not an option

    std::stringstream buffer;
    for (int i = 0; i < num_opts; i++) {
        auto& opt = opts[i];
        auto pos = inOpt.find(opt.name, 1);
        if (pos != 1)
            continue;

        if (inOpt[pos += opt.name.length()] != '=') {
            int k = i;
            while (++k < num_opts && inOpt.find(opts[k].name, 1) == string::npos);
            if (k < num_opts)   // There is another option with the same initial characters.
                continue;       // We will parce it later.

            int j = opt.limits[0];
            buffer << "Invalid option settings: " SQUOT(arg);
            buffer << "\nExpected '-" << opt.name << "= val', where 'val' is one of {" << j++;
            for (; j <= opt.limits[1]; j++)
                buffer << ", " << j;

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
    opt.extra_param = new int [4];
    memset(opt.extra_param, 0, 4 * sizeof(int));
}

void init_DPS(opt_descr& opt) {

    if (opt.intValue != 2) {
        const auto pos = opt.strValue.find_first_of(",|:");
        if (pos != string::npos) {
            init_Cntr(opt);
            auto valStr = opt.strValue.substr(0, pos);
            auto* pntr = static_cast<int*>(opt.extra_param);
            if (is_number(valStr))
                pntr[0] = atoi(valStr.c_str());
            else {
                transform(valStr.begin(), valStr.end(), valStr.begin(), ::toupper);
                if (valStr[pos] == '+' || valStr.find("ALL"))
                    pntr[0] = INT_MAX;
                else
                    cout << "Cannot parse parameter `" << opt.strValue
                    << "' set to option '" << opt.name << "'";
            }
        }
    } else
        opt.strValue = "solutions.txt";

    if (!opt.strValue.empty()) {
        FILE* f;
        fopen_S(f, opt.strValue.c_str(), "w");
        fclose(f);
    }
}


int check_DPS_inequality(int b, int r, int k, int λ, int d = 0, int limitM = 1)
{
// Check inequality by P.Dobcsanyi, D.A.Preece, L.H.Soicher:
//    m*(k-y)*(k-y-1) <= (y+1)*y*(b-d) - 2*y*k*r + k*(k-1)*λ
// where m is the multiplicity of a block and y is any integer.
    const auto yMax = 2 * (k + 1);
    int i = 0;    // we don't need to check 0, because for y = 0, m <= λ and it's already covered
    while (++i < yMax) {
        int y = -3 * i;
        while ((y += 2 * i) <= i) {
            if (y == k || (y == k - 1))
                continue;     // to avoid division by 0

            const auto dividend = (y + 1) * y * (b - d) - 2 * y * k * r + k * (k - 1) * λ;
            const auto divider = (k - y) * (k - y - 1);
            if (dividend / divider <= limitM)
                return 0; // there are no multiple blocks
        }
    }
    return 1;
}

int check_Simplicity(int v, int b, int r, int k, int λ) {
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
                else {
                    if (!check_DPS_inequality(b, r, k, λ))
                        return 0; // there are no multiple blocks

                    if (!check_DPS_inequality(b, r, k, λ, 1))
                        return 3; // got some informations regarding parameter d for non-simple designs. 

                    return -1;
                }
            }
        }
        else {
            // When we are here k = 3 and λ + 2 > v
            // Any such design will have a duplicated blocks
            // Proof: Concider the blocks having λ common point. Any other poins cannot below to two of these blocks, 
            // we would have two identical blocks. Thus this design will contain at least  λ + 2 >= v points
            idx = 2;
        }
    }

    return idx;
}

output check_Simplicity(const opt_descr& opt, int v, int b, int r, int k, int λ, string& comment) {
   static const char* pComment[] = { "S", "S-", "!S", "DPS: m=2, d=0"};
   int idx = check_Simplicity(v, b, r, k, λ);
   if (idx < 0 /* || idx == 3 */) {
#if CHECK_ON_COMPLEMENT
       idx = check_Simplicity(v, b, b - r, v - k, b - 2 * r + λ);
       if (idx < 0)
#endif
           return output::no_comment;
   }

    auto* pntr = static_cast<int*>(opt.extra_param);
    ++*(pntr + idx);
    comment += pComment[idx];
    return output::with_comment;
}

bool is_simple(const string& comment, bool& flag, int* pVal = nullptr) {
    static const char* tag = "DPS: m";
    static auto lenTag = strlen(tag);
    const auto comStr = comment.c_str();
    auto pos = comment.find(tag);
    flag = false;
    if (pos != string::npos) {
        pos += lenTag;
        if (comStr[pos++] != '=')
            pos++;

        const auto posEnd = comment.find(" ");
        string valStr = comment.substr(pos);
        if (posEnd > pos)
            valStr = valStr.substr(0, posEnd - pos);

        const int val = atoi(valStr.c_str());
        flag = val == 2;
        if (pVal)
            *pVal = val;

        return false;
    }

    pos = comment.find('S');
    if (pos != string::npos) {
        if (!(comStr[pos + 1] == '-' || (pos > 0 && comStr[pos - 1] == '!'))) {
            if (pVal)
                *pVal = 1;

            return true;
        }
    }

    return false;
}

output check_DPS_condition(const opt_descr& opt, int v, int b, int r, int k, int λ, string& comment) {
    bool flag;
    if (is_simple(comment, flag))
        return output::no_comment;

    string opt_comment, saved_comment;

    // Check inequality by P.Dobcsanyi, D.A.Preece, L.H.Soicher:
    //    m*(k-y)*(k-y-1) ≤ (y+1)*y*b - 2*y*k*r + k*(k-1)*λ
    // where m is the multiplicity of a block and y is any integer.
    const auto yMax = 2 * (k + 1);

    // Loop over all possible multiplicities
    // It is obvious that if the inequality with fixed (b,r,k,λ) fails for pair (m,y), it also fails for (m+1,y).
    // Therefore, we must exit the next loop when we first see that the inequality is not satisfied
    bool inequality_holds = true;
    const int mMax = flag ? 2 : λ;
    int mLast = -1;

    int m = 1;
    while (inequality_holds && (++m <= mMax)) {
        int i = 0;  // we don't need to check 0, because for y = 0, m <= λ and it's already covered
        int lastReported = -1;
        while (inequality_holds && (++i < yMax)) {
            int y = -3 * i;
            // Loop for y = -i, i
            while ((y += 2 * i) <= i) {
                if (y == k || y == k - 1) {
                    // The BPS inequality can be rewritten as
                    //  0 ≤ (k+1)*k*b - 2*(k^2)*r + k*(k-1)*λ    for y=k or 
                    //  0 ≤ k*(k-1)*b - 2*(k-1)*k*r + k*(k-1)*λ  for y=k-1
                    // 
                    // Since k > 1, we can rewrite them as
                    //  0 ≤ (k+1)*b - 2*k*r + (k-1)*λ    for y=k or 
                    //  0 ≤ b - 2*r + λ                  for y=k-1
                    //
                    // OR
                    //  0 ≤ k*(b - 2*r + λ) + b - λ    for y=k or 
                    //  0 ≤ b - 2*r + λ                for y=k-1
                    //
                    // They both hold, because when 2*k ≤ v, we do have 2*r ≤ b, and 1≤ λ ≤ b
                    // No need to check
                    continue;
                }

                const auto dividend = (y + 1) * y * b - 2 * y * k * r + k * (k - 1) * λ;
                const auto divider = (k - y) * (k - y - 1);
                inequality_holds = m * divider <= dividend;
                if (!inequality_holds) {
                    if (opt_comment != saved_comment)
                        opt_comment = saved_comment;

                    if (!flag && mLast != m - 1)
                        opt_comment += " m<=" + to_string(m - 1);

                    // Define the maximum for disjoint block number for previous value of m
                    break;
                }

                if (dividend == m * divider) {
                    if (!opt_comment.empty())
                        opt_comment += ",";

                    if (!flag)
                        opt_comment += " m=" + to_string(mLast = m);

                    int yOut = y;
                    if (λ > m) {
                        //    In that case there are some nonnegative integer solutions for
                        // x{0} +         x{y} +         x{y+1} = b - m
                        //              y*x{y} +   (y+1)*x{y+1} = k(r - m)
                        //        y*(y-1)*x{y} * (y+1)*y*x{y+1} = k(k − 1)(λ − m)
                        //
                        assert(m != λ);
                        const auto const1 = k * (k - 1) * (λ - m);
                        assert(const1 % y == 0);
                        const auto x1 = k * (r - m) - (const1 / y);
                        const auto const2 = k * (r - m) * (1 - y) + const1;
                        assert(const2 % (y + 1) == 0);
                        const auto x2 = const2 / (y + 1);
                        assert(x1 + x2 == b - m);
                        assert(x1 >= 0);
                        assert(x2 >= 0);
                        if (x1 && x2) {
                            opt_comment += " (" + to_string(y) + "," + to_string(y + 1) + ")";
                            continue;
                        }

                        if (x1) {
                            if (lastReported == yOut) {
                                // Removing last ',' from previous output
                                const auto pos = opt_comment.rfind(',');
                                if (pos != string::npos)
                                    opt_comment = opt_comment.substr(0, pos);

                                continue;
                            }
                        } else
                            yOut++;
                    } 
                    else {
                        assert(y == 1);
                        assert(b - λ == (r - λ) * k);
                    }

                    opt_comment += " (" + to_string(lastReported = yOut) + ")";
                }
            }
        }

        saved_comment = opt_comment;
    }

    auto* pntr = static_cast<int*>(opt.extra_param);
    if (opt_comment.empty()) {
        ++*(pntr + 1);
        return output::no_comment;
    }


    ++*(pntr);
    comment = flag? opt_comment : opt.name + ":" + opt_comment;
    return flag? output::adjust_comments : output::with_comment;
}

#define USE_CHECK   0

#if USE_CHECK
void check(const int *solution, const int* right_part, int last_dx = 0) {
    static int r_part[3], last;
    if (!solution) {
        last = last_dx;
        for (int i = 0; i < 3; i++)
            r_part[i] = right_part[i];

        return;
    }

    int sum[3] = { solution[0] + solution[1], solution[1], 0 };
    for (int i = 2; i <= last; i++) {
        sum[0] += solution[i];
        sum[1] += i * solution[i];
        sum[2] += i * (i - 1) * solution[i] / 2;
    }

    for (int i = 0; i < 3; i++) {
        if (sum[i] + right_part[i] != r_part[i])
            assert(false);
    }
}
#else
#define check(x, y, ...)
#endif

#define MAX_SOLUTION_NUMB   10

void outSolution(const int *solution, const int len, const int shift, string& info_s, string& info_m) {
    info_s.erase();
    for (int i = 0; i <= len; i++) {
        if (solution[i])
            info_s += info_m = " n" + to_string(i + shift) + "=" + to_string(solution[i]);
    }
}

output solve_DPS_system(const opt_descr& opt, int v, int b, int r, int k, int λ, string& comment) {
    if (v == b || λ == 1)
        return output::no_comment;

    bool flag;
    int mMax = λ;
    const auto simple = is_simple(comment, flag, &mMax);
//    if (simple)
//        return output::no_comment;

    // Investigating integer solutions for following system:
    // 
    // x{0} +         x{y} +         x{y+1} = b - m
    //              y*x{y} +   (y+1)*x{y+1} = k(r - m)
    //        y*(y-1)*x{y} * (y+1)*y*x{y+1} = k(k − 1)(λ − m)

    int nd = 1;
    const char *dps_comment = "DPS: m=2, d=0";
    auto pos = comment.find(dps_comment);
    const auto rc = pos != string::npos;
    if (rc) {
        mMax = 2;
        string valStr = comment.substr(pos + strlen(dps_comment));
        pos = valStr.rfind("(");
        if (pos != string::npos) {
            valStr = valStr.substr(pos + 1);
            pos = valStr.rfind(")");
            valStr = valStr.substr(0, pos);
            if (!is_number(valStr)) {
                pos = valStr.find(",");
                assert(pos != string::npos);
                valStr = valStr.substr(pos+1);
            }

            nd = atoi(valStr.c_str());
        }
    }

    string cmnt;
    if (λ == 2 || nd > 1) {
        const int xd = nd > 1
            ? k * ((k - 1) * (λ - 2) - (r - 2) * (nd - 2)) / nd
            : simple
            ? k * (k - 1) * (λ - 1) / 2
            : k * (r - 2);
        const int x1 = nd > 1
                       ? (k * (r - 2) - nd * xd) / (nd - 1)
                       : simple
                       ? k * (r - 1) - 2 * xd
                       : 0;
        const int x0 = b - xd - mMax - x1;
        assert(x0 >= 0);
        string cmnt1 = opt.name;
        if (rc) {
            // When we are here, x0 = 0, k elements which belong to these two
            // identical blocks with all remaining (b - 2) blocks form a BIBD with parameters
            // (k, b-2, r - 2, k(r-2)/(b-2), λ-2)
            const auto k1 = k * (r - 2) / (b - 2);
            assert(k1 - 1 == (λ - 2) * (k - 1) / (r - 2));
            cmnt1 += "+DPS";
        }

        cmnt = " m=" + to_string(simple ? 1 : 2) + " (n";
        if (x0)
            cmnt += "0=" + to_string(x0) + ", n";

        if (simple) {
            comment = cmnt1 + ":" + cmnt + "1=" + to_string(x1) + ", n2=" + to_string(xd) + ")";
            return rc ? output::replace_comments : output::with_comment;
        }


        if (x1)
           cmnt += to_string(nd - 1) + "=" + to_string(x1) + ", n";

        cmnt += to_string(nd) + "=" + to_string(xd) + ")";
    }

    int idx = 3;
    pos = comment.rfind("m<=");
    if (pos == string::npos) {
        idx = 2;
        pos = comment.rfind("m=");
    }

    if (pos != string::npos) {
        string valStr = comment.substr(pos + idx);
        if (!is_number(valStr)) {
            pos = valStr.find_last_of(", (");
            assert(pos != string::npos);
            mMax = atoi(valStr.substr(0, pos).c_str());
        } else
            mMax = atoi(valStr.c_str());
    }

    // limits for block intersection from Connor's theorem
    int cMin = k + λ - r;
    if (cMin < 0)
        cMin = 0;

    int cMax = (2 * k * λ + r * (r - k - λ)) / r;
    if (cMax > k)
        cMax = k; // the intersection in k elements will be considered only in terms of multiplicities

    int rightPart[3];
    int integers[256]; 
    auto len = k + 1;
    // We will keep left part coefficient only for 3d equation
    int* leftPart = 2 * len < sizeof(integers) / sizeof(integers[0])? integers : new int[2 * len];
    int* solution = leftPart + len;
    memset(leftPart, 0, len * sizeof(leftPart[0]));

    // idx0, n0 = n1 =...= n{idx0-1} = 0; n{idx0} is the lowest, which is not necessarily 0
    for (int i = 0; i <= cMax - cMin; i++) {
        leftPart[i] = (i - 1) * i / 2;
    }

    int idx0Adj = 0;
    bool useAdj = false;
    if (!cmnt.empty()) {
        // We already know all solutions for the maximum possible multiplicity
        mMax--;
    }
    else {
        useAdj = (!simple) && comment.find("d=0") != string::npos;
        if (useAdj)
            idx0Adj = 1;
    }


    FILE* f = NULL;
    if (!opt.strValue.empty()) {
        fopen_S(f, opt.strValue.c_str(), "a");
        if (f)
            fprintf(f, "Solutions for (v, b, r, k, lambda) = (%d, %d, %d, %d, %d)\n", v, b, r, k, λ);
    }

    int m = 0;
    string info, info_m, info_s;
    while (++m <= mMax) {
        int cntr = 0;
        int idx0 = cMin;
        if (useAdj && m == mMax && idx0 < idx0Adj)
            idx0 = idx0Adj;

        rightPart[0] = b - m;
        rightPart[1] = k * (r - m);
        rightPart[2] = k * (k - 1) * (λ - m) / 2;
        for (int i = 0; i < idx0; i++) {
            rightPart[1] -= rightPart[0];
            rightPart[2] -= rightPart[1];
            assert(rightPart[1] >= 0);
            assert(rightPart[2] >= 0);
        }

        len = cMax - cMin;
        if (cMax == k)
            len--;

        bool run_loop = true;;
        int numSolutions = 0;
        if (!rightPart[2]) {
            solution[0] = rightPart[0] - (solution[1] = rightPart[1]);
            if (solution[0] >= 0) {
                outSolution(solution, 1, idx0, info_s, info_m);
                numSolutions = 1;
            }

            run_loop = false;
        }

        memset(solution, 0, (len+1) * sizeof(solution[0]));

        const bool check_validity = m == λ - 1 && !idx0;
        const char* pNumSolPrefix = "#:";
        int i = len + 1;
        check(NULL, rightPart, len);
        int step = -1;
        int val, val1;
        while (run_loop) {
            while ((i += step) >= 0 && i <= len) {
                if (step > 0) {
                    if (solution[i]) {
                        --solution[i];
                        rightPart[2] += leftPart[i];
                        rightPart[1] += i;
                        rightPart[0]++;
                        check(solution, rightPart);
                        step = -1;
                    }
                    else {
                        if (i >= len) {       // the current highest coordinate has been set to 0
                            break;
                            if (--len == 2)   // decrease the number of the last modified coordinate
                                break;
                            step = -1;
                        }
                    }

                    continue;
                }

                switch (i) {
                case 2:     val = rightPart[2];
                            if (rightPart[0] >= val && rightPart[1] >= val<<1) {
                                rightPart[2] = 0;
                                rightPart[1] -= val<<1;
                                rightPart[0] -= solution[2] = val;
                                check(solution, rightPart);
                            }
                            else {
                                step = 1;
                                if (val = solution[i = 3]) {
                                    rightPart[2] += val * leftPart[i];
                                    rightPart[1] += val * i;
                                    rightPart[0] += val;
                                    solution[i] = 0;
                                    check(solution, rightPart);
                                }
                            }
                            break;
                case 1:     val = rightPart[1];
                            if (rightPart[0] >= val) {
                                rightPart[1] = 0;
                                rightPart[0] -= solution[1] = val;
                                check(solution, rightPart);
                            }
                            else {
                                step = 1;
                                if (val = solution[i = 2]) {
                                    rightPart[2] += val;
                                    rightPart[1] += val << 1;
                                    rightPart[0] += val;
                                    solution[i] = 0;
                                    check(solution, rightPart);
                                }
                            }
                            break;
                case 0:     solution[0] = rightPart[0];
                            rightPart[0] = 0;
                            break;
                default:
                            val = rightPart[2] / leftPart[i];
                            val1 = rightPart[1] / i;
                            if (val > val1)
                                val = val1;

                            if (val > rightPart[0])
                                val = rightPart[0];

                            // A necessary condition for the existence of non-negative solutions
                            // which will satisfy the second equation of our system
                            //
                            // Should we decide to keep the value for solution[i] on the next steps
                            // we will need to cover the deficit:
                            auto deficit1 = rightPart[1] - val * i;
                            const auto maxNonZeroCoordinates = rightPart[0] - val;
                            const auto mult = rightPart[2] ? i - 1 : 1;
                            if (maxNonZeroCoordinates * mult < deficit1) {
                                // We are using the maximal possible value for solution[i-1]
                                // with its coefficient leftPart[1][i - 1]
                                step = 1;
                                continue;
                            }

#if 0  // This works but the execution time is a bit worse
                            // If the previous criterium did not work, we could try a more sophisticated one, if rightPart[2] > 0;
                            // Define the maximal index of coordinate which could help us to cover the mentioned deficite.

                            int j = 2;
                            while (j < i && rightPart[2] / leftPart[j])
                                j++;

                            if (maxNonZeroCoordinates * (j - 1) < deficit) {
                                step = 1;
                                continue;
                            }
#endif
                            auto deficit2 = rightPart[2] - val * leftPart[i];
                            auto deficit0 = rightPart[0] - val;
                            while (val && deficit2 < deficit1 && deficit0 < deficit1 - deficit2) {
                                val--;
                                deficit0++;
                                deficit1 += i;
                                deficit2 += leftPart[i];
                            }

                            if ((rightPart[0] - val) * mult < deficit1) {
                                step = 1;
                                continue;
                            }

                            if (solution[i] = val) {
                                rightPart[2] -= val * leftPart[i];
                                rightPart[1] -= val * i;
                                rightPart[0] -= val;
                                check(solution, rightPart);
                            }
                }

                if (!rightPart[0]) {
                    if (rightPart[1]) {
                        // It is impossible to satisfy the second equation
                        // We need to reject the coordinate just assigned
                        val = solution[i];
                        rightPart[2] += val * leftPart[i];
                        rightPart[1] += val * i;
                        rightPart[0] += val;
                        solution[i] = 0;
                        check(solution, rightPart);

                        while (i < len && !solution[++i]);
                        if (i == len)
                            break;
                    }
                    else
                        if (!rightPart[2])
                            break;

                    step = 1;
                    if (!i--)
                        rightPart[i = 0] += solution[0];
                }
            }


            while (len >= 0 && !solution[len])
                len--;

            if (!rightPart[0] && !rightPart[1] && !rightPart[2]) {
                // Solution found
                bool valid_solution = true;
                if (check_validity) {
                    // the solution for m = λ - 1 is not valid if there is a pair
                    // of different indices (i, j): i != j, i + j > k and n{i} != 0 n{j} != 0
                    // OR n{i} > 1 for 2 * i > k
#if 1
                    // This check is stronger than the other one
                    int nBlock = 0;
                    int x = k;
                    for (int i = len + 1; --i > 2; ) {
                        for (int t = solution[i]; t-- && i > nBlock;)
                            x -= i - nBlock++;
                    }

                    if (x < 0)
                        valid_solution = false;
#else
                    const auto half_k = (k + 1) >> 1;
                    int i = len + 1;
                    while (--i > half_k) {
                        const auto ni = solution[i];
                        if (!ni)
                            continue;

                        if (ni > 1)
                            break;

                        int j;
                        const auto jMin = k - (j = i) + 1;
                        while (--j > jMin && !solution[j]);

                        if (j > jMin)
                            break;
                    }
                    if (i > half_k)
                        valid_solution = false;
#endif
                }

                if (valid_solution) {
                    if (f) {
                        char buffer[512], * pBuff;
                        const auto lenBuff = sizeof(buffer) / sizeof(buffer[0]);
                        if (!cntr++) {
                            pBuff = buffer;
                            fprintf(f, "m = %d", m);
                            for (int i = 0; i <= len; i++)
                                pBuff += sprintf_s(pBuff, lenBuff - (pBuff - buffer), "   n%d", i + idx0);

                            fprintf(f, "%s\n", buffer);
                        }

                        pBuff = buffer;
                        pBuff += sprintf_s(pBuff, lenBuff - (pBuff - buffer), "#%3d:", cntr);
                        for (int i = 0; i <= len; i++)
                            pBuff += sprintf_s(pBuff, lenBuff - (pBuff - buffer), "%5d", solution[i]);

                        fprintf(f, "%s\n", buffer);
                    }

                    if (!numSolutions++)
                        outSolution(solution, len, idx0, info_s, info_m);

                    if (numSolutions >= MAX_SOLUTION_NUMB) {
                        pNumSolPrefix = "#:>=";
                        if (!f)
                            break;
                    }
                }
            }

            step = 1;
            rightPart[2] += solution[2];
            rightPart[0] += solution[0] + solution[1] + solution[2];
            rightPart[1] += solution[1] + (solution[2] << 1);
            solution[0] = solution[1] = solution[2] = 0;
            check(solution, rightPart);

            if (len <= 2)   // no reasons to change x{2}, it is the oly one in the 3-d equation
                break;      // no more solutions for current m

            i = 2;
            while (!solution[++i]);
            i--;
        }

        if (!info.empty())
            info += ";";

        info += " m=" + to_string(m) + " ";
        if (numSolutions <= 1) {
            if (numSolutions)
                info += "(" + info_s.substr(1) + ")";
            else
                info += " NO_SOLUTIONS";
        } else
            info += pNumSolPrefix + to_string(numSolutions) + info_m;
    }

    if (leftPart != integers)
        delete[] leftPart;

    if (!info.empty() && !cmnt.empty())
        info += ";";

    if (f)
        fclose(f);

    comment = opt.name + ":" + info + cmnt;
    return rc ? output::replace_comments : output::with_comment;
}

output check_BRC(const opt_descr& opt, int v, int b, int r, int k, int λ, string& comment) {
    static int oddPrime[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };

    const auto brc_condition = opt.intValue;
    if (v != b)
        return brc_condition == 2 ? output::nothing : output::no_comment;

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
                if (n < p_quare + (p << 1) + 4)
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
    ++* (pntr + (num ? 1 : 0));
    return output::with_comment;
}

output check_Residual(const opt_descr& opt, int v, int b, int r, int k, int λ, string& comment) {
    static opt_descr opt_BRC = { "BRC", check_BRC, nullptr, { 1, 1 }, 1, "" };
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

void summary_title(opt_descr& opt, const out_struct& out, int total) {
    output_func(out, "Among all ", total, " BIBDs parameter sets:");
}

void report_Simplicity(opt_descr& opt, const out_struct& out, int total) {
    auto* pntr = static_cast<int*>(opt.extra_param);
    output_func(out, "   ", *pntr, " correspond to simple, ",
        *(pntr + 1), " not simple, ", *(pntr + 2), " \"anti-simple\" designs, ",
        total - (*pntr + *(pntr + 1) + *(pntr + 2)), " unknown.");

    delete[](int*)opt.extra_param;
}

void report_DPS(opt_descr& opt, const out_struct& out, int total) {
    auto* pntr = static_cast<int*>(opt.extra_param);
    output_func(out, "   ", *pntr, " (out of ", *pntr + *(pntr+1), " tested) have limitations for block intersections.");
    delete[](int*)opt.extra_param;
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

int parsingParameters(int argc, char* argv[], int num_opts, opt_descr* opts,  
    const wchar_t* param, int numParam, int *minVal, int* maxVal, int* pCntr) {
    const char delimiter[] = { '[', '(', ',', ':', ']', ')' };
    string error;
    size_t prevPos, pos;
    for (int i = 1; i < argc; i++) {
        string arg(argv[i]);
        transform(arg.begin(), arg.end(), arg.begin(), ::toupper);

        if (arg == "-H" || arg == "-HELP") {
            string prog_name(argv[0]);
            const auto pos = prog_name.rfind('\\');
            if (pos != string::npos)
                prog_name = prog_name.substr(pos+1);
            
            cout << "Description of PARAMETERS:";
            cout << "\nBy default, the program " << prog_name << " will generate all possible sets";
            cout << "\nof BIBD parameters (v,b,r,k,lambda) for 3<=r<= 41, 3<=k<=r lexicografically ordered by (r, k, lambda).";
            cout << "\n\nThese limits are subject to change. To do that, use one of Param=value, Param=[min,max] or Param[min,max]";
            cout << "\nwhere Param is one of {v,V,b,B,r,R,k,K,l,L}";
            cout << "\nand \"value\" or [min, max] define, respectively, the value OR range of values for that parameter.";
            cout << "\nWhen one of the {min, max} is omitted, the program will use the appropriate default value.";
            cout << "\n\nFor instance, for r[10,] v[,128], the program will construct all sets of parameters for 10<=r<=41 and v<=128.";
            cout << "\n\nEach set of parameters generated by the program is provided by its number";
            cout << "\nfrom \"Handbook of Combinatorial Designs\" by Charles J.Colbourn, Jeffrey H.Dinitz";
            cout << "\n\nIf for at least one parameter a non-default value for range is specified,"; 
            cout << "\nthe program will provide a separate numbering for the sets of parameters (v,b,r,k,lambda),";
            cout << "\nand for every set it will still indicate its number from \"Handbook\"";
            cout << "\n\n\n";
            cout << "Description of OPTIONS:   (any options can be specified in lowercase or uppercase letters, but without spaces)";
            cout << "\n  -f=[+/-]fileName : Outputting results to \"fileName\".  When one of +/- precedes \"fileName\"";
            cout << "\n                     the output will also (respectively, will not) be displayed on the screen.";
            cout << "\n                     By default, the program prints results to the screen, unless -f=\"fileName\" is used";
            cout << "\n";
            cout << "\n   -brc[=N]        : Bruck-Ryser-Chowla conditions";
            return -1;
        }

        if (!set_option(arg, argv[i], opts, num_opts, error)) {
            if (!error.empty()) {
                cout << error << "\n";
                return -1;
            }
            continue;
        }

        // Parameter settings
        const auto symb = arg[0];
        int j = numParam;
        while (j-- && symb != param[j]);
        if (j < 0) {
            cout << "Can't parse parameter " SQUOT(argv[i]);
            return -1;
        }

        if (j != 2 || minVal[2] != 3) {
            // We found a parameter with new limits that is not 'r'
            // OR it 'r', but its lower limit is different from its default value.
            *pCntr = 0;   // Activating an additional counter
        }

        pos = arg.find("=");
        if (pos != string::npos) {
            const string valStr = arg.substr(pos+1);
            if (is_number(valStr)) {
                minVal[j] = maxVal[j] = atoi(valStr.c_str());
                continue;
            }
        }

        prevPos = 0;
        int nPrev = 0;
        for (int m = 0; m < 6; m += 2) {    // loop for delimiters
            int n = 2;                      // loop for delimiters type
            while (n-- && (pos = arg.find(delimiter[m + n], prevPos)) == string::npos);
            if (n < 0) {
                cout << "Cannot find either of the two expected delimiters (" SQUOT(delimiter[m])
                    ", " SQUOT(delimiter[m+1]) ") in " SQUOT(argv[i]);
                return -1;
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
    return 0;
}

bool set_comments(const opt_descr* opts, int num_opts, int v, int b, int r, int k, int λ, string& comment) {
    string tmp;
    for (int j = 1; j < num_opts; j++) {
        const check_fn check_func = opts[j].check_func;
        if (!check_func)
            continue;

        tmp = comment;
        switch ((*check_func)(opts[j], v, b, r, k, λ, tmp)) {
            case with_comment:      comment += (comment.empty()? " " : "; ") + tmp;
                                    break;
            case adjust_comments:   comment += tmp;
            case no_comment:        break;
            case replace_comments:  comment = " " + tmp;
                                    break;
            case nothing:           return false;
        }
    }

    return true;
}

int main(int argc, char* argv[])
{
    const wchar_t param[] = { 'V', 'B', 'R', 'K', 'L' };
    const char delimiter[] = { '[', '(', ',', ':', ']', ')' };
    int minVal[5] = { 0, 0, 3, 3, 1 };
    int maxVal[5] = { INT_MAX, INT_MAX, 41, INT_MAX, INT_MAX };
    opt_descr opts[] = { {"F", nullptr, summary_title, {}, -1, ""},
                         {"S", check_Simplicity, report_Simplicity, {-1, 2}, 1, "", init_Cntr},
                         {"DPS", check_DPS_condition, report_DPS, {-1, 2}, 1, "", init_Cntr}, //  Dobcsanyi, Preecec, Soicherc condition for the equality holding in their inequality
                         {"SYS", solve_DPS_system, nullptr, {0, 2}, -1, "", init_DPS},
                         {"BRC", check_BRC, report_BRC, {-1, 2}, 1, "", init_Cntr},
                         {"O", check_Residual, report_Residual, {-1,1}, 1, "", init_Cntr},
                       };

    const auto num_opts = sizeof(opts) / sizeof(opts[0]);

    const auto start = clock();

    int cntr = -1;  // by default this counter will not be used
    const int numParam = sizeof(param) / sizeof(param[0]);
    if (parsingParameters(argc, argv, num_opts, opts, param, numParam, minVal, maxVal, &cntr) < 0)
        exit(1);

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

    std::stringstream buffer;
    buffer << std::endl;
    out_struct out = { buffer, outFile, coutFlag };

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
 
                string comment;
                if (!set_comments(opts, num_opts, v, b, r, k, λ, comment))
                    continue;

                pBuff = buff;
                if (cntr >= 0)
                    pBuff += sprintf_s(pBuff, sizeof(buff), "%4d ", ++cntr);

                sprintf_s(pBuff, sizeof(buff) - (pBuff - buff), "#%4d: %4d %4d %2d %2d %2d  ", total, v, b, r, k, λ);
                output_func(out, buff, comment);
            }
        }
    }


    for (int j = 0; j < num_opts; j++) {
        const report_fn report_func = opts[j].report_func;
        if (report_func)
            (*report_func)(opts[j], out, cntr < 0? total : cntr);
    }

    auto runTime = (unsigned long)(clock() - start) / CLOCKS_PER_SEC;
    const auto sec = runTime % 60;
    auto min = (runTime -= sec) / 60;
    if (min > 60) {
        const auto hours = (min - (min % 60)) / 60;
        output_func(out, "Run Time: ", hours, " hrs ", min % 60, " min ", sec, " sec");
    } else
        output_func(out, "Run Time: ", min, " min ", sec, " sec");

    if (outFile)
        fclose(outFile);

    for (int j = 0; j < num_opts; j++) {
        const auto destr_func = opts[j].constr_func[1];
        if (destr_func)
            (*destr_func)(opts[j]);
    }
}

