/// calc Remez for e^x [-1,0,1]

#include <iostream>
#include <vector>
#include <numbers>


#include <eigen3/Eigen/Dense>

template<typename T,typename F>
double derivative(T f, F x, T h) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

template<typename T,typename F>
double derivative(T f, F x) {
    double h = sqrt(std::numeric_limits<F>::epsilon());
    return (f(x + h) - f(x - h)) / (2.0 * h);
}


template<typename T>
double  get_max_fun(T f, double a, double b){

    std::size_t  i = 1000;

    double startx = 0.5;

    while(startx <= 1 and startx >= -1 and i != 0){
        auto d = derivative(f,  startx);
        startx += d * 0.1;
        --i;
    }
    return startx;
}

template<typename T = double >
std::vector<T> polif(T a, T b, std::size_t n) {
    auto tn = static_cast<T>(n);
    auto pi = std::numbers::pi;

    std::vector<T> points(n + 1);

    for (std::size_t i = 0; auto &point: points) {
        auto ti = static_cast<T>(i);
        point = (a + b) / 2. + (b - a) / 2. * std::cos(  (ti * pi) / (tn + 1.0));
        i++;
    }


    return points;
}

template<typename T>
using dmatrix = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;

template<typename F>
auto remez_step(F f, const dmatrix<double>& A, const Eigen::VectorXd& ys){

    Eigen::VectorXd x = A.colPivHouseholderQr().solve(ys);

    auto aprox_f = [coef = x](auto x){ return coef[0] + coef[1]*x;};

    auto err = [&aprox_f,&f](auto x){return  std::abs(f(x)-aprox_f(x));};

    auto e = get_max_fun(err, ys[0],ys[ys.size()-1]);

    Eigen::VectorXd r(x.size()-1);

    for(std::size_t i = 0; i != x.size()-1; ++i){
        r[i] = x[i];
    }

    return std::make_tuple(r,e);
}

template<typename F>
auto remez_init(F f, const std::vector<double> &points){

    const auto n = points.size();

    std::vector<double> ys(points.size());

    std::ranges::transform(points, ys.begin(), f);


    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(n,n);
    Eigen::VectorXd v(n);
    for(std::size_t i = 0; i != n; ++i){
        a(i,0) = 1;
        a(i,n-1) = std::pow(-1,i);
        a(i,1) = points[i];
        v[i] = f(points[i]);
    }





    return remez_step(f,a,v);
}





int main() {


    //auto xs = polif(0., 5., 2);
    std::vector<double>xs{-1,0,1};


    auto [a,b] = remez_init([](auto x){return std::exp(x);}, xs);

    std::cout << a << "\n\n" << b << '\n';




    return 0;
}
