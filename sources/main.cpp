#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <thread>

#include <boost/timer/timer.hpp>
#include <boost/system/system_error.hpp>
#include <boost/core/demangle.hpp>
#include <boost/core/typeinfo.hpp>

#include <include/Matrix_type.hpp>

static constexpr int NI = 3;
static constexpr int NJ = 3;
static constexpr int NK = 3;
static constexpr int NL = 3;


static const std::array<size_t, 1> shape1D = { {NI} };
static const std::array<size_t, 2> shape2D = { {NI, NJ} };
static const std::array<size_t, 3> shape3D = { {NI, NJ, NK} };
static const std::array<size_t, 4> shape4D = { {NI, NJ, NK, NL} };

struct print_type
{
    template <class T>
    void operator() (T) const
    {
        auto const& ti = BOOST_CORE_TYPEID(T);
        std::cout << boost::core::demangled_name(ti) << std::endl;
    }
};

int main(int argc, char** argv)
{

    try {

        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_2D ..." << std::endl;
            MATRIX<2, int> A2(shape2D);
            std::cout << "CTOR MATRIX_2D ..." << std::endl;
            A2.Random(1, 8);
            std::cout << "RANDOM MATRIX_2D ..." << std::endl;
            std::cout << A2 << std::endl;
            std::cout << "DETERMINAT 2D = " << A2.DET() << std::endl;
            std::cout << "END TEST MATRIX_2D " << tmr.format() << std::endl;
        }
        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_3D ..." << std::endl;
            MATRIX<3, int> A3(shape3D);
            A3.Random(1, 12);
            size_t count3D{};
            do {
                count3D++;
            } while (A3.DET_FULL() != 0);
            std::cout << "Детерминант матрицы 3D ориентации i = " << A3.DET_orient('i') << std::endl;
            std::cout << "Детерминант матрицы 3D ориентации j = " << A3.DET_orient('j') << std::endl;
            std::cout << "Детерминант матрицы 3D ориентации k = " << A3.DET_orient('k') << std::endl;
            std::cout << "Полный детерминант матрицы 3D = " << A3.DET_FULL() << std::endl;
            for(const auto& x : A3.transversal('j'))
            std::cout << "MATRIX_3D трансверсальные сечения = \n" << x << std::endl;
            std::cout << "END TEST MATRIX_3D, ITERATION = " << count3D << tmr.format() << std::endl;

        }
        {
            boost::timer::cpu_timer tmr;
            std::cout << "START TEST MATRIX_4D ..." << std::endl;
            MATRIX<4, int> A4(shape4D);
            A4.Random(1, 7);
            size_t count4D{};
            do {
                count4D++;
            } while (A4.DET_FULL() != 0);
            std::cout << "Детерминант матрицы 4D ориентации i = " << A4.DET_orient('i') << std::endl;
            std::cout << "Детерминант матрицы 4D ориентации j = " << A4.DET_orient('j') << std::endl;
            std::cout << "Детерминант матрицы 4D ориентации k = " << A4.DET_orient('k') << std::endl;
            std::cout << "Детерминант матрицы 4D ориентации l = " << A4.DET_orient('l') << std::endl;
            std::cout << "Полный детерминант матрицы 4D = " << A4.DET_FULL() << std::endl;
            for(const auto& x : A4.transversal('j'))
            std::cout << "MATRIX_4D трансверсальные сечения по j = \n" << x << std::endl;
            std::cout << "END TEST MATRIX_4D, ITERATION = " << count4D << tmr.format() << std::endl;
        }

    }
    catch (boost::system::system_error& ec) {
        std::cerr << "Error ocured! # " << ec.code() << " Messaga "
            << ec.what() << std::endl;
    }


    return 0;
}
