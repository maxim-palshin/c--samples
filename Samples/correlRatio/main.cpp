/*

    calc correlation Ratio for

    dataset

    age|brand|
    ___|_____|
    27 |  1  |
    33 |  2  |
    16 |  3  |
    29 |  3  |
    32 |  2  |
    23 |  1  |
    25 |  2  |
    28 |  1  |
    22 |  3  |
    18 |  3  |
    26 |  2  |
    26 |  1  |
    15 |  3  |
    29 |  2  |
    26 |  3  |




 */


#include <iostream>
#include <map>
#include <numeric>
#include <vector>
#include <cmath>

//calc mean
double mean(const std::vector<double> &x) {
    auto sum = std::accumulate(x.begin(), x.end(), 0.0);
    return sum / static_cast<double>(x.size());
}

//X
//a vector containing the quantitative variable.
//y
//a vector containing the qualitative variable (e.g. a factor).
double correlRatio(const std::vector<double> &x, const std::vector<std::size_t> &y) {

    if (x.size() != y.size()) {
        throw std::runtime_error("error size not eq");
    }

    std::map<std::size_t, std::vector<double>> groups;

    for(std::size_t i = 0; i != x.size(); ++i){
        groups[y[i]].push_back(x[i]);
    }


    double intraclass_variance = 0;
    double groups_mean = 0;
    double result=0;

    for (auto [_, value]: groups) {
        auto avg = mean(value);
        groups_mean += avg/static_cast<double>(groups.size());
        double variance = 0;
        for (auto i: value) {
            variance += std::pow(i - avg, 2);
        }
        intraclass_variance += variance;
    }

    for(auto [_, value]: groups){
        std::size_t  count = value.size();
        result += static_cast<double>(count) * std::pow(mean(value)-groups_mean,2);
    }


    return result/(result+intraclass_variance);
}






int main() {

    std::vector<double> age {27, 33, 16, 29, 32, 23, 25, 28, 22, 18, 26, 26, 15, 29, 26};
    std::vector<std::size_t> brand {1, 2, 3, 3, 2, 1, 2, 1, 3, 3, 2, 1, 3, 2, 3};

    std::cout << correlRatio(age, brand) << '\n';

    return 0;
}
